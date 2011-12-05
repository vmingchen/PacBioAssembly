/*
 * ===========================================================================
 *
 *       Filename:  ref_seq.h
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/29/2011 12:19:14 AM
 *
 *    Description:  sequence reference
 *
 *       Revision:  none
 *
 *
 * ===========================================================================
 */

#ifndef REF_SEQ_H
#define REF_SEQ_H

#include	<list>

#include	"dna_seq.h"
#include	"seq_aligner.h"
#include	"common.h"

// the first edit cannot be INSERT
template <class Iter> 
void apply_edits(edit *pedit, int nedit, Iter it, bool forward) {
    for (int i = 0; i < nedit; ++i) {
        if (pedit->op == DELETE) {
            it->ignore();
            ++it;
        } else if (pedit->op == MATCH) {
            it->select(pedit->val);
            ++it;
        } else if (pedit->op == INSERT) {
            if (forward) --it;
            it->supply(pedit->val);
            if (forward) ++it;
        }
        ++pedit;
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
    void absorb(base_vote& other) {
        acgt[0] += other.acgt[0];
        acgt[1] += other.acgt[1];
        acgt[2] += other.acgt[2];
        acgt[3] += other.acgt[3];
        other.reset();
    }
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
    vote_box() : total(0) {}
    vote_box(char c) : selection(c), total(1) {}
    base_vote selection;
    base_vote suppliment;
    int total; 
    void select(char c) { selection.add_char(c); ++total; }
    void ignore() { ++total; }
    void supply(char c) { suppliment.add_char(c); }
    void split(vote_box *other) {
        other->selection = suppliment;
        other->total = total;
        suppliment.reset();
    }
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
    bool try_align(t_aligner *paligner, int pos, seq_accessor *pac_seg) {
        bool forward = pac_seg->is_forward();
        seq_accessor ac_ref = get_accessor(pos, forward);
        // don't mistake the order of the two parameters
        // pac_seg now behave like a reference
        if (paligner->align(&ac_ref, pac_seg) < 0) return false;
        if (paligner->matlen_a < OVERLAP_MIN) return false;
        elect(pos, paligner->edits, paligner->nedit, forward);
        if (paligner->matlen_a == ac_ref.length()) {
            int add_len = pac_seg->length() - paligner->matlen_b;
            if (forward) {
                append(pac_seg->pt(paligner->matlen_b), add_len);
            } else {
                prepend(pac_seg->pt(pac_seg->length()-1), add_len);
            }
        }
        return true;
    }
    seq_accessor get_accessor(int pos, bool forward) {
        assert(contained(pos));
        return seq_accessor(txt_buf + beg + pos, forward, 
                forward ? post-beg-pos : pos+beg-pre+1);
    }
    /**
     * build (rebuild) seedmap for reference sequence
     **/
    unsigned get_seedmap(hash_table &seedmap, t_seed sd_pat) {
        int len = end - beg;
        int nmax = len - N_SEQ_WORD;
        int nhead = std::min(nmax, MAX_READ_LEN);
        seedmap.clear();
        char *ptext = txt_buf + beg;
        for (int i = 0; i < nhead; ++i) {
            unsigned sd = dna_seq::encode(ptext++);
            // there are a lot of 'AAAAAAAAAAAAAAAA' segments, ignore them
            if (sd & sd_pat) seedmap[sd & sd_pat].push_back(i);
        }

        int ntail = std::min(len-MAX_READ_LEN-N_SEQ_WORD, MAX_READ_LEN);
        ptext = txt_buf + end - N_SEQ_WORD;
        for (int i = 0; i < ntail; ++i) {
            unsigned sd = dna_seq::encode(ptext--);
            if (sd & sd_pat) seedmap[sd & sd_pat].push_back(len-i-N_SEQ_WORD);
        }

        return nhead + (ntail < 0 ? 0 : ntail);
    };
    // seedmap will be invalidated once evolved
    void evolve() {
        end = pre = beg = MAX_SEQ_LEN;
        char *p = txt_buf + beg;
        vote_box vb;
        std::list<vote_box>::iterator prev;
        std::list<vote_box>::iterator next;
        std::list<vote_box>::iterator cur = consensus.begin();
        while (cur != consensus.end()) {
            if (cur->has_supply(0.5)) {     // insert
                cur->split(&vb);
                next = cur;
                ++next;
                if (next == consensus.end()) 
                    consensus.push_back(vb);
                else
                    consensus.insert(next, vb);
            }
            if (cur->is_valid(0.5)) {       // match 
                *p++ = cur->get_vote(); 
                ++end; 
                ++cur;
            } else {                        // delete
                if (cur != consensus.begin()) { // has prev
                    prev = cur;
                    --prev;
                    prev->suppliment.absorb(cur->selection);
                } 
                cur = consensus.erase(cur);
            }
        }
        post = end;
    }
    // pos should be contained
    void elect(int pos, edit *pedit, int nedit, bool forward) {
        if (forward) {
            std::list<vote_box>::iterator it = consensus.begin();
            advance(it, pos + beg - pre);
            apply_edits(pedit, nedit, it, forward);
        } else {
            std::list<vote_box>::reverse_iterator it = consensus.rbegin();
            advance(it, post - beg - pos - 1);
            apply_edits(pedit, nedit, it, forward);
        }
    }
private:
    int beg;        // origin of current iteration
    int end;        // end of current iteration
    int pre;        // extension before beg
    int post;       // extension after end

    char txt_buf[3*MAX_SEQ_LEN];
    unsigned char bin_buf[4+MAX_SEQ_LEN/N_SEQ_BYTE];
    std::list<vote_box> consensus; 
};

#endif
