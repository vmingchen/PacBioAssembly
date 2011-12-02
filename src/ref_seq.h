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
#include	"seq_aligner.h"
#include	"common.h"

template <class Iter> 
void apply_edits(edit *pedit, int nedit, Iter it) {
    for (int i = 0; i < nedit; ++i) {
        if (pedit->op == DELETE) {
            it->ignore();
            ++it;
        } else if (pedit->op == MATCH) {
            it->select(pedit->val);
            ++it;
            ++pedit;
        } else if (pedit->op == INSERT) {
            it->supply(pedit->val);
            ++pedit;
        }
    }
}

class base_vote {
private:
    unsigned short acgt[4];
public:
    base_vote() { reset(); }
    base_vote(char c) { reset(); add_char(c); }
    void add_char(char  c) { ++acgt[C2I(c)]; }
    void add_code(int c) { ++acgt[c]; }
    void reset() { acgt[0] = acgt[1] = acgt[2] = acgt[3] = 0; }
    int max_vote() { 
        return std::max(acgt[0], std::max(acgt[1], 
                    std::max(acgt[2], acgt[3])));
    }
    char winner() {
        int mv = max_vote();
        return mv == acgt[0] ? 'A' 
            : (mv == acgt[1] ? 'C' : (mv == acgt[2] ? 'G' : 'T'));
    }
}; 

class vote_box {
public:
    vote_box(char c) : selection(c), total(1) {}
    base_vote selection;
    base_vote suppliment;
    int total; 
    void select(char c) { selection.add_char(c); ++total; }
    void ignore() { ++total; }
    void supply(char c) { suppliment.add_char(c); ++total; }
    bool is_valid(double ratio) { return selection.max_vote() > ratio*total; }
    bool has_supply(double ratio) { return suppliment.max_vote() > ratio*total; } 
    char get_vote() { return selection.winner(); }
    char get_supply() { return suppliment.winner(); }
};

class ref_seq {
public:
    ref_seq(const t_bseq *pseq) {
        beg = pre = MAX_SEQ_LEN;
        end = post = beg + dna_seq::bin2text(pseq, txt_buf+beg, MAX_SEQ_LEN);
        char *p = txt_buf + beg;
        for (int i = beg; i < end; ++i)
            consensus.push_back(vote_box(*p++));
    };
    void append(char *pseg, int len) {
        memmove(txt_buf + post, pseg, len);
        post += len;
        for (int i = 0; i < len; ++i) {
            consensus.push_back(vote_box(*pseg++));
        }
    }
    void prepend(char *pseg, int len) {
        pre = pre - len;
        memmove(txt_buf + pre, pseg, len);
        pseg = pseg + len - 1;
        for (int i = 0; i < len; ++i) {
            consensus.push_front(vote_box(*pseg--));
        }
    }
    bool contained(int pos) { return pos+beg >= pre && pos+beg < post; }
    unsigned length() { return end - beg; }
    seq_accessor get_accessor(int pos, bool forward) {
        assert(contained(pos));
        return seq_accessor(txt_buf + beg + pos, forward, 
                forward ? post-beg-pos : pos+beg-pre);
    }
    /* *
     * build (rebuild) seedmap for reference sequence
     * */
    unsigned get_seedmap(hash_table &seedmap, t_seed sd_pat) {
        unsigned nseed = end - beg - N_SEQ_WORD;
        char *ptext = txt_buf + beg;
        for (unsigned i = 0; i < nseed; ++i) {
            unsigned sd = dna_seq::encode(ptext++);
            // there are a lot of 'AAAAAAAAAAAAAAAA' segments, ignore them
            if (sd & sd_pat) seedmap[sd & sd_pat].push_back(i);
        }
        return nseed;
    };
    // seedmap will be invalidated once evolved
    void evolve() {
        pre = beg = MAX_SEQ_LEN;
        std::deque<vote_box>::iterator it = consensus.begin();
        char *p = txt_buf + beg;
        while (it != consensus.end()) {
            if (it->is_valid(0.5)) *p++ = it->get_vote();
            if (it->has_supply(0.5)) *p++ = it->get_supply();
            ++it;
        }
        post = end = (p - (txt_buf + beg));
        consensus.clear();
        p = txt_buf + beg;
        for (int i = beg; i < end; ++i) 
            consensus.push_back(vote_box(*p++));
    }
    void elect(int pos, edit *pedit, int nedit, bool forward) {
        if (forward) {
            std::deque<vote_box>::iterator it = consensus.begin();
            it += (beg + pos);
            apply_edits(pedit, nedit, it);
        } else {
            std::deque<vote_box>::reverse_iterator it = consensus.rbegin();
            it += (post - pos);
            apply_edits(pedit, nedit, it);
        }
    }
private:
    int beg;        // origin of current iteration
    int end;        // end of current iteration
    int pre;        // extension before beg
    int post;       // extension after end

    char txt_buf[3*MAX_SEQ_LEN];
    unsigned char bin_buf[4+MAX_SEQ_LEN/N_SEQ_BYTE];
    std::deque<vote_box> consensus; 
};

#endif
