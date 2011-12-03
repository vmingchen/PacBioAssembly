/*
 * ===========================================================================
 *
 *       Filename:  aligner_test.cpp
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/26/2011 09:31:12 PM
 *
 *    Description:  tester for local_aligner
 *
 *       Revision:  none
 *
 * ===========================================================================
 */

#include <gtest/gtest.h>
#include <seq_aligner.h>
#include	<fstream>
#include	<string>

char dna_ref[] = "ACGTAACCGGTT";
char dna_seg1[] = "CGTAAGC";
char dna_seg2[] = "GTAACGGGTTAA";
char dna_seg3[] = "TCGTAAC";

seq_aligner *paligner = NULL;

class aligner_test : public testing::Test {
protected:
    void edit_tester(seq_accessor *pref, seq_aligner *paligner) { 
        int j = 0;      // index into dna_ref
        for (int i = 0; i < paligner->nedit; ++i) {
            char op = paligner->edits[i].op;
            if (op == MATCH) { 
                EXPECT_EQ(pref->at(j), paligner->edits[i].val);
                ++j;
            } else if (op == INSERT) {
                EXPECT_EQ(pref->at(j), paligner->edits[i].val);
                ++j; 
            } 
        }
    };
};

TEST_F(aligner_test, forward) {
    paligner = new seq_aligner();
    seq_accessor ref1((char*)dna_ref, true, 7);
    seq_accessor seg1((char*)dna_seg1, true, 6);
    EXPECT_LE(6, paligner->align(&seg1, &ref1));
    EXPECT_GE(7, paligner->align(&seg1, &ref1));
    EXPECT_EQ(2, paligner->final_cost());
    edit_tester(&ref1, paligner);

    seq_accessor ref2((char*)dna_ref, true, 8);
    seq_accessor seg2((char*)dna_seg1, true, 7);
    EXPECT_EQ(7, paligner->align(&seg2, &ref2));
    EXPECT_EQ(2, paligner->final_cost());
    edit_tester(&ref2, paligner);

    seq_accessor ref3((char*)dna_ref, true, 8);
    seq_accessor seg3((char*)dna_seg3, true, 7);
    EXPECT_EQ(7, paligner->align(&seg3, &ref3));
    EXPECT_EQ(1, paligner->final_cost());
    edit_tester(&ref3, paligner);
}

TEST_F(aligner_test, backward) {
    seq_accessor ref((char*)dna_ref + 7, false, 7);
    seq_accessor seg((char*)dna_seg1 + 6, false, 7);
    EXPECT_EQ(7, paligner->align(&seg, &ref));
    EXPECT_EQ(1, paligner->final_cost());
    edit_tester(&ref, paligner);
}

TEST_F(aligner_test, overlay) {
    seq_accessor ref((char*)dna_ref + 2, true, 10);
    seq_accessor seg((char*)dna_seg2, true, 12);
    EXPECT_EQ(10, paligner->align(&seg, &ref));
    EXPECT_EQ(1, paligner->final_cost());
    edit_tester(&ref, paligner);
}

TEST_F(aligner_test, remove) {
    seq_accessor ref(dna_ref, true, 10);
    seq_accessor seg(dna_ref+1, true, 9);
    EXPECT_EQ(10, paligner->align(&seg, &ref));
    EXPECT_EQ(10, paligner->nedit);
    EXPECT_EQ(INSERT, paligner->edits[0].op);
    EXPECT_EQ(1, paligner->final_cost());
    edit_tester(&ref, paligner);

    ref.reset(0);
    seg.reset(0);
    EXPECT_EQ(9, paligner->align(&ref, &seg));
    EXPECT_EQ(10, paligner->nedit);
    EXPECT_EQ(DELETE, paligner->edits[0].op);
    EXPECT_EQ(1, paligner->final_cost());
    edit_tester(&seg, paligner);
}

TEST_F(aligner_test, sample) {
    std::ifstream fin("test/real_align.txt");
    std::string ref_str, seg_str;
    fin >> ref_str >> seg_str;
    seq_accessor ref((char*)ref_str.c_str() + ref_str.length()-1, 
            false, ref_str.length());
    seq_accessor seg((char*)seg_str.c_str() + seg_str.length()-1,  
            false, seg_str.length());
    EXPECT_LT(0, paligner->align(&seg, &ref));
    edit_tester(&ref, paligner);

    fin >> ref_str >> seg_str;
    seq_accessor ref2((char*)ref_str.c_str(), 
            false, ref_str.length());
    seq_accessor seg2((char*)seg_str.c_str(),  
            false, seg_str.length());
    EXPECT_EQ(-1, paligner->align(&seg2, &ref2));
}

