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

char dna_ref[] = "ACGTAACCGGTT";
char dna_seg1[] = "CGTAAG";

TEST(seq_aligner, basic1) {
    seq_aligner aligner;
    seq_accessor ref((char*)dna_ref, true, 12);
    seq_accessor seg((char*)dna_seg1, true, 6);
    EXPECT_EQ(6, aligner.align(&seg, &ref));
    EXPECT_EQ(21, aligner.get_score(5, 6));

    int j = 0;      // index into dna_ref
    for (int i = 0; i < aligner.nedit; ++i) {
        char op = aligner.edits[i].op;
        if (op == MATCH) { 
            EXPECT_EQ(ref.at(j), aligner.edits[i].val);
            ++j;
        } else if (op == INSERT) {
            EXPECT_EQ(ref.at(j), aligner.edits[i].val);
            ++j; 
        } 
    }
}
