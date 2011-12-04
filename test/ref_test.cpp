/*
 * ===========================================================================
 *
 *       Filename:  ref_test.cpp
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  12/01/2011 11:30:52 PM
 *
 *    Description:  tester for ref_seq
 *
 *       Revision:  none
 *
 * ===========================================================================
 */

#include <gtest/gtest.h>
#include <ref_seq.h>

TEST(base_vote, basic) {
    base_vote vote('A');
    EXPECT_EQ('A', vote.winner());
    EXPECT_EQ(1, vote.max_vote());

    vote.add_char('C');
    vote.add_char('C');
    EXPECT_EQ('C', vote.winner());
    EXPECT_EQ(2, vote.max_vote());

    vote.add_code(2);
    vote.add_code(2);
    vote.add_code(2);
    EXPECT_EQ('G', vote.winner());
    EXPECT_EQ(3, vote.max_vote());

    vote.reset();
    vote.add_code(3);
    EXPECT_EQ('T', vote.winner());
    EXPECT_EQ(1, vote.max_vote());
}

TEST(vote_box, basic) {
    vote_box box('T');
    EXPECT_EQ(true, box.is_valid(0.5));
    EXPECT_EQ(false, box.has_supply(0.5));
    EXPECT_EQ('T', box.get_vote());

    box.ignore();
    EXPECT_EQ(false, box.is_valid(0.5));
    EXPECT_EQ(false, box.has_supply(0.5));

    box.select('C');
    EXPECT_EQ(false, box.is_valid(0.5));
    EXPECT_EQ(false, box.has_supply(0.5));

    box.select('C');
    box.select('C');
    EXPECT_EQ(true, box.is_valid(0.5));
    EXPECT_EQ(false, box.has_supply(0.5));
    EXPECT_EQ('C', box.get_vote());
    EXPECT_EQ(5, box.total);

    box.supply('T');
    box.supply('T');
    EXPECT_EQ(true, box.is_valid(0.5));
    EXPECT_EQ(false, box.has_supply(0.5));
    box.supply('T');
    EXPECT_EQ(true, box.is_valid(0.5));
    EXPECT_EQ(true, box.has_supply(0.5));
    EXPECT_EQ('T', box.get_supply());
}

char dna_txt[]  = "ACGTAACCGGTTAAACCCGGGTTTTGCAAAAAAAAAAAAAAAA";
char dna_txt1[] = "ACGTAACCGGTTAAACCCGGGTGTTGCAAAAAAAAAAAAAAAA";
char dna_txt2[] = "ACGTAACCGGTTAAACCCGGGTTGTTGCAAAAAAAAAAAAAAAA";
char dna_txt3[] = "ACGTAACCGGTTAAACCCGGGTTGGTTGCAAAAAAAAAAAAAAAA";
char dna_txt4[] = "ACGTAACCGGTTAAACCCGGGTTGTTGCAAAAAAAAAAAAAAAAGGCCTTAA";
char dna_txt5[] = "ACGTAACCGGTTAAACCCGGGTTGTTGCAAAAAAAAAAAAAAAAGGCCTTAAC";
char dna_txt6[] = "TTTTACGTAACCGGTTAAACCCGGGTTGTTGCAAAAAAAAAAAAAAAA";
char dna_txt7[] = "TTTTTACGTAACCGGTTAAACCCGGGTTGTTGCAAAAAAAAAAAAAAAA";
char dna_post[] = "CGT";
char dna_pre[] = "TGC";

class ref_test : public testing::Test {
protected:
    virtual void SetUp() {
        unsigned char bseg[24];
        sz = strlen(dna_txt);
        dna_seq::text2bin(dna_txt, bseg, 24);
        pref = new ref_seq(bseg);
        paligner = new seq_aligner();
    }
    virtual void TearDown() {
        if (pref) delete pref;
        if (paligner) delete paligner;
    }
    int sz;
    ref_seq *pref;
    seq_aligner *paligner;
};

TEST_F(ref_test, basic) {
    EXPECT_EQ(sz, pref->length());

    // test contained
    EXPECT_EQ(false, pref->contained(-1));
    EXPECT_EQ(true, pref->contained(0));
    EXPECT_EQ(true, pref->contained(sz-1));
    EXPECT_EQ(false, pref->contained(sz));

    // test get_accessor
    seq_accessor facsr = pref->get_accessor(0, true);
    for (int i = 0; i < sz; ++i) {
        EXPECT_EQ(dna_txt[i], facsr.next());
    }
    seq_accessor bacsr = pref->get_accessor(sz-1, false);
    for (int i = sz-1; i >= 0; --i) {
        EXPECT_EQ(dna_txt[i], bacsr.next());
    }

    // test get_seedmap
    hash_table seedmap;
    pref->get_seedmap(seedmap, 0xFFFFFFFF);
    // -1 for the tailing sequence with all A
    EXPECT_EQ(sz-15-1, seedmap.size());
    for (int i = 0; i < sz-16; ++i) {
        EXPECT_EQ(true, seedmap.find(dna_seq::encode(dna_txt+i)) 
                != seedmap.end());
    }
    EXPECT_EQ(true, seedmap.find(dna_seq::encode(dna_txt+sz-15)) == seedmap.end());
}

TEST_F(ref_test, grow) { 
    int sz_post = strlen(dna_post);
    pref->append(dna_post, sz_post);
    EXPECT_EQ(true, pref->contained(sz + sz_post - 1));
    EXPECT_EQ(false, pref->contained(sz + sz_post));

    int sz_pre = strlen(dna_pre);
    pref->prepend(dna_pre, sz_pre);
    EXPECT_EQ(true, pref->contained(-sz_pre));
    EXPECT_EQ(false, pref->contained(-sz_pre-1));
    EXPECT_EQ(sz, pref->length());
}

TEST_F(ref_test, change) { 
    seq_accessor ac_seg(dna_txt1, true, strlen(dna_txt1));
    EXPECT_EQ(true, pref->try_align(paligner, 0, &ac_seg));
    ac_seg.reset(0);
    EXPECT_EQ(true, pref->try_align(paligner, 0, &ac_seg));

    pref->evolve();
    seq_accessor ac_ref = pref->get_accessor(0, true);
    for (ac_seg.reset(0); ac_seg.has_more(); )
        EXPECT_EQ(ac_seg.next(), ac_ref.next());
}

TEST_F(ref_test, remove) {
    seq_accessor fac_seg(dna_txt+1, true, sz-1);
    EXPECT_EQ(true, pref->try_align(paligner, 0, &fac_seg));
    fac_seg.reset(0);
    EXPECT_EQ(true, pref->try_align(paligner, 0, &fac_seg));
    pref->evolve();
    EXPECT_EQ(sz-1, pref->length());
    fac_seg = pref->get_accessor(0, true);
    for (int i = 1; i <= sz-1; ++i) 
        EXPECT_EQ(dna_txt[i], fac_seg.next());
}

TEST_F(ref_test, insert) {
    int sz_seg = strlen(dna_txt2);
    seq_accessor fac_seg(dna_txt2, true, sz_seg);
    EXPECT_EQ(true, pref->try_align(paligner, 0, &fac_seg));
    EXPECT_EQ(sz_seg, paligner->nedit);
    fac_seg.reset(0);
    EXPECT_EQ(true, pref->try_align(paligner, 0, &fac_seg));
    pref->evolve();
    EXPECT_EQ(sz_seg, pref->length());
    fac_seg = pref->get_accessor(0, true);
    for (int i = 0; i < sz_seg; ++i) 
        EXPECT_EQ(dna_txt2[i], fac_seg.next());
}

TEST_F(ref_test, insert2) {
    int sz_seg = strlen(dna_txt3);
    seq_accessor fac_seg(dna_txt3, true, sz_seg);
    EXPECT_EQ(true, pref->try_align(paligner, 0, &fac_seg));
    pref->evolve();
    EXPECT_EQ(sz+1, pref->length());
    fac_seg = pref->get_accessor(0, true);
    for (int i = 0; i <= sz; ++i) 
        EXPECT_EQ(dna_txt2[i], fac_seg.next());
}

TEST_F(ref_test, back_insert) {
    int sz_seg = strlen(dna_txt2);
    seq_accessor bac_seg(dna_txt2 + sz_seg - 1, false, sz_seg);
    EXPECT_EQ(true, pref->try_align(paligner, sz-1, &bac_seg));
    EXPECT_EQ(sz_seg, paligner->nedit);
    bac_seg.reset(0);
    EXPECT_EQ(true, pref->try_align(paligner, sz-1, &bac_seg));
//    fac_seg.reset(0);
//    for (int i = 0; i < paligner->nedit; ++i) 
//        printf("%d %c, ", paligner->edits[i].op, paligner->edits[i].val);
//    printf("\n");
    pref->evolve();
    EXPECT_EQ(sz_seg, pref->length());
    bac_seg = pref->get_accessor(sz_seg-1, false);
    for (int i = sz_seg-1; i >= 0; --i) 
        EXPECT_EQ(dna_txt2[i], bac_seg.next());
}

TEST_F(ref_test, back_insert2) {
    int sz_seg = strlen(dna_txt3);
    seq_accessor bac_seg(dna_txt3 + sz_seg - 1, false, sz_seg);
    EXPECT_EQ(true, pref->try_align(paligner, sz-1, &bac_seg));
    EXPECT_EQ(sz_seg, paligner->nedit);
    pref->evolve();
    EXPECT_EQ(sz+1, pref->length());
    bac_seg = pref->get_accessor(sz, false);
    for (int i = sz; i >= 0; --i) 
        EXPECT_EQ(dna_txt2[i], bac_seg.next());
}

TEST_F(ref_test, append) {
    int sz4 = strlen(dna_txt4);
    seq_accessor fac_seg(dna_txt4, true, sz4);
    EXPECT_EQ(true, pref->try_align(paligner, 0, &fac_seg));
    EXPECT_EQ(true, pref->contained(sz+1));
    int sz5 = strlen(dna_txt5);
    seq_accessor fac_seg2(dna_txt5, true, sz5);
    EXPECT_EQ(true, pref->try_align(paligner, 0, &fac_seg2));
    pref->evolve();
    EXPECT_EQ(sz5, pref->length());
    fac_seg2 = pref->get_accessor(0, true);
    for (int i = 0; i < sz5; ++i)
        EXPECT_EQ(dna_txt5[i], fac_seg2.next());
}

TEST_F(ref_test, prepend) {
    int sz6 = strlen(dna_txt6);
    seq_accessor bac_seg6(dna_txt6 + sz6 - 1, false, sz6);
    EXPECT_EQ(true, pref->try_align(paligner, sz-1, &bac_seg6));
    EXPECT_EQ(true, pref->contained(-1));
    int sz7 = strlen(dna_txt7);
    seq_accessor bac_seg7(dna_txt7 + sz7 - 1, false, sz7);
    EXPECT_EQ(true, pref->try_align(paligner, sz-1, &bac_seg7));
//    for (int i = 0; i < paligner->nedit; ++i) 
//        printf("%d %c, ", paligner->edits[i].op, paligner->edits[i].val);
//    printf("\n");
    pref->evolve();
    EXPECT_EQ(sz7, pref->length());
    bac_seg7 = pref->get_accessor(sz7-1, false);
    for (int i = sz7-1; i >= 0; --i) 
        EXPECT_EQ(dna_txt7[i], bac_seg7.next());
}
