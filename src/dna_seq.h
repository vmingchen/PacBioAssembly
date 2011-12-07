/*
 * ===========================================================================
 *
 *       Filename:  dna_seq.h
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/08/2011 11:59:58 AM
 *
 *    Description:  This file defines types related a DNA sequence
 *    including dna_seq, seq_accessor and a couple of macros. 
 *
 *       Revision:  none
 *
 * ===========================================================================
 */
#ifndef DNA_SEQ_H
#define DNA_SEQ_H
#include	<assert.h>
#include	"common.h"

//! convert DNA base to binary number
#define C2I(x) ((x == 'A') ? 0 : ((x == 'C') ? 1 : (x == 'G' ? 2 : 3)))
//! convert binary number to DNA base
#define I2C(x) ((x == 0) ? 'A' : ((x == 1) ? 'C' : (x == 2 ? 'G' : 'T')))

//! number of DNA bases in a word
#define N_SEQ_WORD 16
//! number of DNA bases in a byte
#define N_SEQ_BYTE 4

static const char codes[4] = {'A', 'C', 'G', 'T'};

struct bin_seq {
    int quality;
    int length;
    unsigned char *data;
};

struct txt_seq {
    int quality;
    char *data;
};

/**
 * Represents a DNA sequence. It is most used for its static functions
 * because in our progrom a DNA sequence is either a char pointer (text) or
 * a unsigned char pointer (binary). 
 **/
class dna_seq {
public:
//    // create binary file from text file 
//    void create(const char *in_path, const char *out_path);
    // read binary file
    bool parse(const char *path) { return false; };
    // is empty or not
    bool empty(){ return true; };
    // read a record from the binary
    const unsigned* read(unsigned *plen) const { return NULL; };

    /**
     * Return the seed at 'pos' of binary sequence pbin. 
     **/
    static t_seed seed_at(unsigned char *pbin, int pos) {
        pbin += sizeof(unsigned);
        if ((pos & 0x3) == 0) return *((unsigned*)(pbin+pos));

        unsigned char pseed[4];
        pbin += (pos >> 2);
        unsigned ls = (pos & 0x3) << 1;
        unsigned rs = 0x8 - ls;
        pseed[0] = (*pbin << ls) | (*++pbin >> rs);
        pseed[1] = (*pbin << ls) | (*++pbin >> rs);
        pseed[2] = (*pbin << ls) | (*++pbin >> rs);
        pseed[3] = (*pbin << ls) | (*++pbin >> rs);

        return *((unsigned*)pseed);
    }

    static char value_at(unsigned char bv, int idx) {
        return codes[(bv >> ((~idx & 0x3) << 1)) & 0x3];
    }

    /**
     * Return the seed pointed by ptext, which is supposed to be at least
     * 16-character-long. 
     */
    static unsigned encode(const char *ptext) {
        unsigned tseg = 0;
        unsigned char *t = (unsigned char*)&tseg;

        t[0] = t2b(ptext, 4);
        t[1] = t2b(ptext+4, 4);
        t[2] = t2b(ptext+8, 4);
        t[3] = t2b(ptext+12, 4);

        return tseg;
    };

    /**
     * Write back the DNA characters represents by code into ptext. 
     **/
    static void decode(unsigned code, char *ptext) {
        unsigned char *t = (unsigned char*)&code;
        b2t(*t++, ptext, 4);
        b2t(*t++, ptext+4, 4);
        b2t(*t++, ptext+8, 4);
        b2t(*t++, ptext+12, 4);
    };

    /**
     * Convert text sequence to binary sequence, return the length of the
     * resulting binary sequence.  
     */
    static unsigned text2bin(const char *ptext, unsigned char *pbin, unsigned buflen) {
        size_t tlen = strlen(ptext);
        size_t blen = sizeof(unsigned) + (tlen+4-1)/4; 
        unsigned char *pb = pbin;
        
        assert(buflen >= blen);
        *((unsigned*)pb) = tlen;
        pb += sizeof(unsigned);

        const char *pt = ptext;
        for (int i = tlen; i > 0; i-=4, pt+=4) {
            *pb++ = t2b(pt, i);
        }
        return blen;
    };

    /**
     * Convert binary sequence to text sequence. Return length of
     * sequence. 
     **/
    static unsigned bin2text(const unsigned char *pbin, char *ptext, unsigned buflen) {
        size_t tlen = *((unsigned *)pbin);
        const unsigned char *pb = pbin + sizeof(unsigned);
        char *pt = ptext; 

        assert(buflen > tlen);
        for (int i = tlen; i > 0; i-=4, pt+=4) {
            b2t(*pb++, pt, i);
        }
        ptext[tlen] = '\0';

        return tlen;
    }
private:
    static unsigned char t2b(const char *pt, size_t tlen) {
        unsigned char b = 0;
        if (tlen == 1) {
            b = C2I(pt[0]) << 6;
        } else if (tlen == 2) {
            b = C2I(pt[0]) << 6 | C2I(pt[1]) << 4;
        } else if (tlen == 3) {
            b = C2I(pt[0]) << 6 | C2I(pt[1]) << 4 | C2I(pt[2]) << 2;
        } else {    // tlen >= 4
            b = C2I(pt[0]) << 6 | C2I(pt[1]) << 4 | C2I(pt[2]) << 2 | C2I(pt[3]);
        }
        return b;
    }
    static void b2t(unsigned char bv, char *pt, size_t tlen) {
        if (tlen == 1) {
            pt[0] = codes[(bv >> 6) & 0x3];
        } else if (tlen == 2) {
            pt[0] = codes[(bv >> 6) & 0x3];
            pt[1] = codes[(bv >> 4) & 0x3];
        } else if (tlen == 3) {
            pt[0] = codes[(bv >> 6) & 0x3];
            pt[1] = codes[(bv >> 4) & 0x3];
            pt[2] = codes[(bv >> 2) & 0x3];
        } else {    // tlen >= 4
            pt[0] = codes[(bv >> 6) & 0x3];
            pt[1] = codes[(bv >> 4) & 0x3];
            pt[2] = codes[(bv >> 2) & 0x3];
            pt[3] = codes[bv & 0x3];
        }
    }
};

/**
 * An helper class facilitate the visit into a text sequence. It supports many
 * ways of access like backward/forward iteration, position access and
 * pointer access. 
 *
 **/
class seq_accessor {
public:
    /**
     * Construct a seq_accessor with *p be the underlying text sequence. The
     * direction and the length of the accessor are defined by f and l. 
     **/
    seq_accessor(char *p, bool f, int l) : pdna(p), pcur(p), forward(f), len(l), cnt(0) {};

    /**
     * Length of the underlying text sequence. 
     **/
    int length() { return len; };

    /**
     * Determine the direction of the accessor. 
     **/
    bool is_forward() { return forward; }

    /**
     * Check if there is more data when in iteration mode. 
     **/
    bool has_more() { return cnt < len; };

    /**
     * Return the next DNA base in iteration mode.
     **/
    char next() { ++cnt; return forward ? *pcur++ : *pcur--; };

    /**
     * Reset the number of read bases to 0. 
     **/
    void reset(int pos) { cnt = pos; pcur = forward ? pdna + pos: pdna - pos; };

    /**
     * Directly access the base at pos i. 
     **/
    char at(int i){ return forward ? *(pdna+i) : *(pdna-i); };

    /**
     * Return a pointer the base at pos i. 
     **/
    char* pt(int i) { return forward ? (pdna+i) : (pdna-i); };
private:
    char *pdna;
    char *pcur;
    int len;
    int cnt;
    bool forward;
};

#endif
