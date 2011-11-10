/*
 * ===========================================================================
 *
 *       Filename:  binary_parser.cpp
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/08/2011 11:59:36 AM
 *
 *    Description:  Parse DNA sequence in Binary format
 *
 *       Revision:  none
 *
 * ===========================================================================
 */
#include "binary_parser.h"

//unsigned text2bin(const char *ptext, unsigned char *pbin, unsigned buflen)
//{
//    int tlen = strlen(ptext);
//    int blen = (tlen+sizeof(unsigned)+4-1)/4; 
//    char *p = pbin;
//    
//    assert(buflen >= blen);
//    *((unsigned*)p) = tlen;
//    p += sizeof(unsigned);
//
//    *p = 0;
//    for (int i = 0; i < tlen; ++i) {
//        c = (ptext[i] == 'T') ? 3 : ((ptext[i] == 'G') ? 2 : (ptext[i] == 'C'));
//        int byte_offset = (i & 0x3) << 1;
//        *p |= (c << byte_offset);
//        if (byte_offset == 0x6) *++p = 0; 
//    }
//    return blen;
//}
//
//unsigned bin2text(const unsigned char *pbin, char *ptext, unsigned buflen)
//{
//    static const char codes[4] = {'A', 'C', 'G', 'T'};
//    int tlen = *((unsigned *)pbin);
//    unsigned char *p = pbin + sizeof(unsigned);
//
//    assert(buflen > tlen);
//    for (int i = 0; i < tlen; ++i) {
//        ptext[i] = codes[(p[i >> 2] >> ((i & 3)<<1)) & 3];
//    }
//
//    return tlen;
//}
//
//unsigned binary_parser::parse(const char *path) 
//{
//    return 0;
//}
//
//
//
