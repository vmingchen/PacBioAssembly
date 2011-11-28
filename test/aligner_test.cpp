/*
 * ===========================================================================
 *
 *       Filename:  aligner_test.cpp
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/26/2011 09:31:12 PM
 *
 *    Description:  test local_aligner
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
char dna_seg2[] = "GTAACGGGTTAAAA";

class aligner_test : public testing::Test {
protected:
    virtual void SetUp() {};
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
    seq_aligner aligner;
    seq_accessor ref1((char*)dna_ref, true, 12);
    seq_accessor seg1((char*)dna_seg1, true, 6);
    EXPECT_LE(6, aligner.align(&seg1, &ref1));
    EXPECT_GE(7, aligner.align(&seg1, &ref1));
    EXPECT_EQ(2, aligner.final_cost());
    edit_tester(&ref1, &aligner);

    seq_accessor ref2((char*)dna_ref, true, 12);
    seq_accessor seg2((char*)dna_seg1, true, 7);
    EXPECT_EQ(7, aligner.align(&seg2, &ref2));
    EXPECT_EQ(2, aligner.final_cost());
    edit_tester(&ref2, &aligner);
}

TEST_F(aligner_test, backward) {
    seq_aligner aligner;
    seq_accessor ref((char*)dna_ref + 7, false, 7);
    seq_accessor seg((char*)dna_seg1 + 6, false, 7);
    EXPECT_EQ(7, aligner.align(&seg, &ref));
    EXPECT_EQ(1, aligner.final_cost());
    edit_tester(&ref, &aligner);
}

TEST_F(aligner_test, overlay) {
    seq_aligner aligner;
    seq_accessor ref((char*)dna_ref + 2, true, 10);
    seq_accessor seg((char*)dna_seg2, true, 14);
    EXPECT_EQ(10, aligner.align(&seg, &ref));
    EXPECT_EQ(1, aligner.final_cost());
    edit_tester(&ref, &aligner);
}

TEST_F(aligner_test, sample) {
    seq_aligner aligner;
    std::ifstream fin("test/real_align.txt");
    std::string ref_str, seg_str;
    fin >> ref_str >> seg_str;
//    std::cout << ref_str << "\n" << seg_str << "\n";
    seq_accessor ref((char*)ref_str.c_str() + ref_str.length()-1, 
            false, ref_str.length());
    seq_accessor seg((char*)seg_str.c_str() + seg_str.length()-1,  
            false, seg_str.length());
    EXPECT_LT(0, aligner.align(&seg, &ref));
    edit_tester(&ref, &aligner);
}

