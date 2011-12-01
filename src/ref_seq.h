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

class base_vote {
public:
    unsigned short acgt[4];
    base_vote() { acgt[0] = acgt[1] = acgt[2] = acgt[3] = total = 0; };
    base_vote(char c) { add_char(c); }
    int max_vote() { 
        return std::max(acgt[0], std::max(acgt[1], std::max(acgt[2], acgt[3])));
    }
    char winner() {
        int mv = max_vote();
        return mv == acgt[0] ? 'A' 
            : (mv == acgt[1] ? 'C' : (mv == acgt[2] ? 'G' : 'T'))
    }
    void add_char(char  c) { ++acgt[C2I(c)]; }
    void add_code(int c) { ++acgt[c]; }
}; 

class vote_box {
public:
    vote_box() : total(0) {}
    base_vote selection;
    base_vote suppliment;
    int total; 
    bool select(int c) { selection.add_code(c); ++total}
    bool supply(int c) { suppliment.add_code(c); ++total}
    bool invalid(double ratio) { return selection.max_vote() < ratio*total; }
}

class ref_seq {
public:
    ref_seq(const t_bseq *pseq) {};
    void append(char *pseg, int len) {
        memmove(txt_buf + post, pseg, len);
        post += len;
    };
    void prepend(char *pseg, int len) {
        pre = pre - len;
        memmove(txt_buf + pre, pseg, len);
    };
    bool contained(int pos) { return pos >= pre && pos < post; }
    unsigned length() { return end - beg; }
    seq_accessor get_accessor(bool forward, int pos) {
        assert(contained(pos));
        return seq_accessor(txt_buf + pos, forward, 
                forward ? post-pos : pos-pre);
    }
    /* *
     * build (rebuild) seedmap for reference sequence
     * */
    unsigned void get_seedmap(hash_table &seedmap, t_seed sd_pat) {
        unsigned nseed = end - beg - N_SEQ_WORD;
        for (unsigned i = 0; i < nseed; ++i) {
            unsigned sd = dna_seq::seed_at(bin_buf, i);
            // there are a lot of 'AAAAAAAAAAAAAAAA' segments, ignore them
            if (sd & sd_pat) seedmap[sd & sd_pat].push_back(i);
        }
        return nseed;
    };
    // seedmap will be invalidated once evolved
    bool evolve() {
    }
    void elect(int pos, edit *pedit, int nedit, bool forward) {
        if (forward) {
            std::deque<vote_box>::iterator it = consensus.begin();
            it += (beg + pos);
        }
        for (int i = 0; i < nedit; ++i) {
            if (pedit->op == DELETE) {

            }
            if (pedit->op == MATCH)
        }
    }

    int beg;        // origin of current iteration
    int end;        // end of current iteration
    int pre;        // extension before beg
    int post;       // extension after end

    char txt_buf[3*MAX_SEQ_LEN];
    unsigned char bin_buf[4+MAX_SEQ_LEN/N_SEQ_BYTE];
    std::deque<vote_box> consensus; 
};

#endif
