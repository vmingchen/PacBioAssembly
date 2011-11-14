/*
 * ===========================================================================
 *
 *       Filename:  binary_parser.h
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/08/2011 11:59:58 AM
 *
 *    Description:  
 *
 *       Revision:  none
 *
 * ===========================================================================
 */
#ifndef BINARY_PARSER_H
#define BINARY_PARSER_H
#include	<assert.h>

#define TR(x) ((x == 'A') ? 0 : ((x == 'C') ? 1 : (x == 'G' ? 2 : 3)))

// number of sequnce in a word
#define N_SEQ_WORD 16
// number of sequnce in a byte
#define N_SEQ_BYTE 4


class binary_parser {
public:
//    // create binary file from text file 
//    void create(const char *in_path, const char *out_path);
    // read binary file
    bool parse(const char *path) { return false; };
    // is empty or not
    bool empty(){ return true; };
    // read a record from the binary
    const unsigned* read(unsigned *plen) const { return NULL; };
    static unsigned encode(const char *ptext) {
        unsigned tseg = 0;
        unsigned char *t = (unsigned char*)&tseg;
        char *p = ptext;
        for (size_t j = 0; j < N_SEQ_WORD; ++j, ++p) {
            *t = (*t << 2) | TR(*p);
            if ((j & 0x3) == 0x3) ++t;
        }
        return tseg;
    }
    static void decode(unsigned code, char *ptext) {

    }
    // convert text file to binary file, return length of result
    static unsigned text2bin(const char *ptext, unsigned char *pbin, unsigned buflen) {
        size_t tlen = strlen(ptext);
        size_t blen = sizeof(unsigned) + (tlen+4-1)/4; 
        unsigned char *p = pbin;
        
        assert(buflen >= blen);
        *((unsigned*)p) = tlen;
        p += sizeof(unsigned);

        *p = 0;
        for (size_t i = 0; i < tlen; ++i) {
            int c = (ptext[i] == 'T') ? 3 : ((ptext[i] == 'G') ? 2 : (ptext[i] == 'C'));
            *p = (*p << 2) | c;
            if ((i & 0x3) == 0x3) *++p = 0;
        }
        if (tlen & 3) *p <<= ((4 - (tlen & 3))<<1);
        return blen;
    };
    static unsigned bin2text(const unsigned char *pbin, char *ptext, unsigned buflen) {
        static const char codes[4] = {'A', 'C', 'G', 'T'};
        size_t tlen = *((unsigned *)pbin);
        const unsigned char *p = pbin + sizeof(unsigned);

        assert(buflen > tlen);
        for (size_t i = 0; i < tlen; ++i) {
            ptext[i] = codes[(*p >> ((3 - (i & 3)) << 1)) & 3];
            if ((i & 0x3) == 0x3) ++p;
        }
        ptext[tlen] = '\0';

        return tlen;
    }
private:
    static unsigned char t2b(const char *pt, size_t tlen) {
        unsigned char b = 0;
        if (tlen == 1) {
            b = TR(pt[0]) << 6;
        } else if (tlen == 2) {
            b = TR(pt[0]) << 6 | TR(pt[1]) << 4;
        } else if (tlen == 3) {
            b = TR(pt[0]) << 6 | TR(pt[1]) << 4 | TR(pt[2]) << 2;
        } else {    // tlen >= 4
            b = TR(pt[0]) << 6 | TR(pt[1]) << 4 | TR(pt[2]) << 2 | TR(pt[3]);
        }
        return b;
    }
};

#endif
