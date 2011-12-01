/*
 * ===========================================================================
 *
 *       Filename:  ref_seq.h
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/29/2011 12:19:14 AM
 *
 *    Description:  
 *
 *       Revision:  none
 *
 *
 * ===========================================================================
 */

#ifndef REF_SEQ_H
#define REF_SEQ_H

#include	<deque>

#include	"dna_seq.h"
#include	"common.h"

class vote {
public:
    unsigned short acgt[4];
    unsigned total; 
    vote() { acgt[0] = acgt[1] = acgt[2] = acgt[3] = total = 0; };
    vote(char c) { add_char(c); }
    char get() {
        if (A > C && A > G && A > T)
            return 'A';
        else if (C > A && C > G && C > T)
            return 'C';
        else if (G > A && G > C && G > T)
            return 'G';
        else
            return 'T';
    }
    void add_char(char  c) { ++acgt[C2I(c)]; ++total; }
    void add_code(int c) { ++acgt[c]; ++total; }
}; 

struct ref_seq {
    ref_seq(const t_bseq *pseq) {};
    void append(char *pseg, int len) {
    };
    void prepend(char *pseg, int len) {
    };
    bool contained(int pos) {
        return pos >= pre && pos < post;
    };

    int beg;
    int end; 
    int pre;
    int post; 
    char seq_buf[3*MAX_SEQ_LEN];
    std::deque<Vote> consensus; 
};

#endif
