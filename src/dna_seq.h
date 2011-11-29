/*
 * ===========================================================================
 *
 *       Filename:  dna_seq.h
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/08/2011 11:59:58 AM
 *
 *    Description:  
 *
 *       Revision:  none
 *
 * ===========================================================================
 */
#ifndef DNA_SEQ_H
#define DNA_SEQ_H
#include	<assert.h>

#define C2I(x) ((x == 'A') ? 0 : ((x == 'C') ? 1 : (x == 'G' ? 2 : 3)))
#define I2C(x) ((x == 0) ? 'A' : ((x == 1) ? 'C' : (x == 2 ? 'G' : 'T')))

// number of sequnce in a word
#define N_SEQ_WORD 16
// number of sequnce in a byte
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
    static unsigned seed_at(unsigned char *pbin, int pos) {
        unsigned char pseed[4];

        pbin += sizeof(unsigned);
        unsigned *p = (unsigned*)(pbin + sizeof(unsigned));

        if ((pos & 0xf) == 0) return *((unsigned*)(pbin+pos));

        pbin += (pos >> 2);
        unsigned ls = pos & 0x3;
        unsigned rs = 0x4 - ls;
        ls <<= 1; rs <<= 1;
        pseed[0] = (*pbin << ls) | (*++pbin >> rs);
        pseed[1] = (*pbin << ls) | (*++pbin >> rs);
        pseed[2] = (*pbin << ls) | (*++pbin >> rs);
        pseed[3] = (*pbin << ls) | (*++pbin >> rs);

        return *((unsigned*)pseed);
    }
    static char value_at(unsigned char bv, int idx) {
        return codes[(bv >> ((~idx & 0x3) << 1)) & 0x3];
    }
    static unsigned encode(const char *ptext) {
        unsigned tseg = 0;
        unsigned char *t = (unsigned char*)&tseg;

        t[0] = t2b(ptext, 4);
        t[1] = t2b(ptext+4, 4);
        t[2] = t2b(ptext+8, 4);
        t[3] = t2b(ptext+12, 4);

        return tseg;
    };
    static void decode(unsigned code, char *ptext) {
        unsigned char *t = (unsigned char*)&code;
        b2t(*t++, ptext, 4);
        b2t(*t++, ptext+4, 4);
        b2t(*t++, ptext+8, 4);
        b2t(*t++, ptext+12, 4);
    };
    // convert text file to binary file, return length of result
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

class seq_accessor {
public:
    seq_accessor(char *p, bool f, int l) : pdna(p), pcur(p), forward(f), len(l), cnt(0) {};
    int length() { return len; };
    bool has_more() { return cnt < len; };
    char next() { ++cnt; return forward ? *pcur++ : *pcur--; };
    void reset(int pos) { cnt = pos; pcur = forward ? pdna + pos: pdna - pos; };
    char at(int i){ return forward ? *(pdna+i) : *(pdna-i); };
    char* pt(int i) { return forward ? (pdna+i) : (pdna-i); };
private:
    char *pdna;
    char *pcur;
    int len;
    int cnt;
    bool forward;
};

#endif
