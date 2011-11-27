/*
 * ===========================================================================
 *
 *       Filename:  dna_test.cpp
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/26/2011 07:34:54 PM
 *
 *    Description:  
 *
 *       Revision:  none
 *
 * ===========================================================================
 */

#include <gtest/gtest.h>
#include <dna.h>

char dna_str[] = "ACGTGTCA";

TEST(dna_accessor, length) {
    dna_accessor da((char *)dna_str, true, 4);
    EXPECT_EQ(4, da.length());
}
