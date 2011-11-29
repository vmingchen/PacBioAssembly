/*
 * ===========================================================================
 *
 *       Filename:  dna_test.cpp
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/26/2011 07:34:54 PM
 *
 *    Description:  test dna_seq and seq_accessor
 *
 *       Revision:  none
 *
 * ===========================================================================
 */

#include <gtest/gtest.h>
#include <dna_seq.h>

char dna_str[] = "ACGTGTCATCGGATCAACCGGTT";

TEST(dna_seq, binary) {
    unsigned char bin_buf[10+4];
    char txt_buf[41];
    EXPECT_EQ(10, dna_seq::text2bin(dna_str, bin_buf, 14));
    EXPECT_EQ(23, dna_seq::bin2text(bin_buf, txt_buf, 41));
    EXPECT_STREQ(dna_str, txt_buf);
    printf("%08x\n", dna_seq::seed_at(bin_buf, 0));
    printf("%08x\n", *((unsigned*)(bin_buf+4)));
    printf("%08x\n", dna_seq::encode(dna_str));
    EXPECT_EQ(0x34DAB41B, dna_seq::seed_at(bin_buf, 0));
    EXPECT_EQ(0x41A34DBB, dna_seq::seed_at(bin_buf, 2));
    EXPECT_EQ(0xAF058D36, dna_seq::seed_at(bin_buf, 7));
}

TEST(seq_accessor, forward) {
    seq_accessor da((char *)dna_str, true, 4);
    EXPECT_EQ(4, da.length());
    EXPECT_EQ('A', da.next());
    EXPECT_EQ('C', da.next());
    EXPECT_EQ('G', da.next());
    EXPECT_EQ('T', da.next()); EXPECT_EQ('G', da.at(2));
    EXPECT_EQ(false, da.has_more());
    da.reset(2);
    EXPECT_EQ(true, da.has_more());
    EXPECT_EQ('G', da.next());
    EXPECT_EQ('T', da.next());
    EXPECT_EQ(false, da.has_more());
}

TEST(seq_accessor, backward) {
    seq_accessor da((char *)dna_str + 4, false, 3);
    EXPECT_EQ(3, da.length());
    EXPECT_EQ('G', da.next());
    EXPECT_EQ('T', da.next());
    EXPECT_EQ(true, da.has_more());
    EXPECT_EQ('G', da.next());
    EXPECT_EQ(false, da.has_more());
    EXPECT_EQ('T', da.at(1));
    da.reset(1);
    EXPECT_EQ('T', da.next());
    EXPECT_EQ('G', da.next());
    EXPECT_EQ(false, da.has_more());
}
