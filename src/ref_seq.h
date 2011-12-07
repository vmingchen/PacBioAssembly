/*
 * ===========================================================================
 *
 *       Filename:  ref_seq.h
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/29/2011 12:19:14 AM
 *
 *    Description:  Three classes are defined in this file: base_vote,
 *    vote_box, ref_seq. 
 *
 *       Revision:  none
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

/**
 * Vote of DNA base class. It represents opinions from different segments
 * overlapped at the same place. 
 * */
class base_vote {
private:
    unsigned short acgt[4];
public:
    /**
     * Constructor of empty vote.
     * */
    base_vote() { reset(); }

    /**
     * Constructor of initial n votes on base c. 
     */
    base_vote(char c, int n = 1) { reset(); add_char(c, n); }

    /**
     * Add a preference on a base specified by its character [ACGT]. 
     **/
    void add_char(char  c) { ++acgt[C2I(c)]; }
    void add_char(int c, int n) { acgt[C2I(c)] += n; }

    /**
     * Add a preference on a base specified by its code [0123]. 
     **/
    void add_code(int c) { ++acgt[c]; }

    /**
     * Reset the vote to be empty. 
     **/
    void reset() { acgt[0] = acgt[1] = acgt[2] = acgt[3] = 0; }

    /**
     * Add all information all another vote and then reset it. 
     **/
    void absorb(base_vote& other) {
        acgt[0] += other.acgt[0];
        acgt[1] += other.acgt[1];
        acgt[2] += other.acgt[2];
        acgt[3] += other.acgt[3];
        other.reset();
    }

    /**
     * The maximum number of votes. 
     **/
    int max_vote() { 
        return std::max(acgt[0], std::max(acgt[1], 
                    std::max(acgt[2], acgt[3])));
    }
    
    /**
     * Base that with the maximum vote. 
     **/
    char winner() {
        int mv = max_vote();
        return mv == acgt[0] ? 'A' 
            : (mv == acgt[1] ? 'C' : (mv == acgt[2] ? 'G' : 'T'));
    }
}; 

/**
 * A vote box at a particular place in the reference. It represents the
 * reference's base value at one place. Every vote box actually has two
 * base_vote, one for base value at that place, another for a possible
 * suppliment right after that place. 
 **/
class vote_box {
public:
    /**
     * Default constructor for an empty vote box.
     **/
    vote_box() : total(0) {}

    /**
     * Constructor with an initial n votes on c. 
     **/
    vote_box(char c, int n=1) : selection(c, n), total(1) {}

    /** 
     * Vote for base value. 
     **/
    base_vote selection;

    /**
     * Vote for possible suppliment. 
     **/
    base_vote suppliment;

    /**
     * How many segments have vote at this particular box (place). 
     **/
    int total; 
    
    /**
     * Add vote. This happens when there is a match or mismatch between the
     * reference and the segment. 
     **/
    void select(char c) { selection.add_char(c); ++total; }

    /**
     * Ignore vote. This happens when the segment want to delete this base
     * in the reference.
     **/
    void ignore() { ++total; }

    /**
     * Provide a suppliment. This happens when the segment want to insert
     * another base right after. 
     **/
    void supply(char c) { suppliment.add_char(c); }

    /**
     * Split suppliment to make it an separated vote_box.
     **/
    void split(vote_box *other) {
        other->selection = suppliment;
        other->total = total;
        suppliment.reset();
    }

    /**
     * Check if the base with maximum vote is large than the ratio of total
     * votes, which make it an effective result. 
     **/
    bool is_valid(double ratio) { return selection.max_vote() > ratio*total; }

    /**
     * Check if there is any effective suppliment. 
     * */
    bool has_supply(double ratio) { return suppliment.max_vote() > ratio*total; } 

    /**
     * Get the vote result. Call is_valid to make sure its return value makes
     * sense. 
     **/
    char get_vote() { return selection.winner(); }

    /**
     * Get the suppliment. Call has_supply to make sure its return value makes
     * sense. 
     **/
    char get_supply() { return suppliment.winner(); }
};


/**
 * Reference sequence. DNA reads (segments) will aligned against it. Once
 * aligned, the segment will express its opinion of the real base at
 * overlapped places. If the read is aligned at the two boundaries of the
 * reference, its unmatched part will be appended or prepened to the
 * reference, so that the reference can grow. However, if the reference is in a
 * locked state, the above-mentioned grow mechanism is disabled. The
 * reference is stable within one iteration. Even its state will be
 * changed during the alignment, but the changes are not disclosed to
 * outside. After iteration, its state can be updated by 'evolve'. 
 **/
class ref_seq {
public:
    /**
     * Constructor of reference with binary sequence.
     * */
    ref_seq(const t_bseq *pseq, bool lk = false) : locked(lk) {
        beg = pre = MAX_SEQ_LEN;
        end = post = beg + dna_seq::bin2text(pseq, txt_buf+beg, MAX_SEQ_LEN);
        char *p = txt_buf + beg;
        for (int i = beg; i < end; ++i)
            consensus.push_back(vote_box(*p++));
    };

    /**
     * Constructor of reference with text sequence. 
     * */
    ref_seq(const char *ptxt, int len, bool l, int w = 1) : locked(l) {
        beg = pre = MAX_SEQ_LEN;
        end = post = beg + len;
        strncpy(txt_buf + beg, ptxt, len);
        char *p = txt_buf + beg;
        for (int i = beg; i < end; ++i)
            consensus.push_back(vote_box(*p++, w));
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
    
    /**
     * Checked if there is anything at pos of reference. 
     * */
    bool contained(int pos) { return pos+beg >= pre && pos+beg < post; }

    /**
     * Return the length of reference. 
     * */
    unsigned length() { return end - beg; }

    /*
     * Try to align pac_seg against the reference starting from pos.
     * Return true on success. The details of the alignment are available in
     * paligner. 
     */
    bool try_align(t_aligner *paligner, int pos, seq_accessor *pac_seg) {
        bool forward = pac_seg->is_forward();
        seq_accessor ac_ref = get_accessor(pos, forward);
        // don't mistake the order of the two parameters
        // pac_seg now behave like a reference
        if (paligner->align(&ac_ref, pac_seg) < 0) return false;
        if (paligner->matlen_a < OVERLAP_MIN) return false;
        if (locked) return true;
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

    /**
     * Get accessor into the reference. 
     * */
    seq_accessor get_accessor(int pos, bool forward) {
        assert(contained(pos));
        return seq_accessor(txt_buf + beg + pos, forward, 
                forward ? post-beg-pos : pos+beg-pre+1);
    }

    /**
     * Build (rebuild) seedmap for reference sequence. 
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

    /**
     * Refresh the reference to reflect updated state. The previous
     * seedmap will be invalidated once evolved. 
     * */
    void evolve() {
        if (locked) return ;
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
    bool locked;    // prevent from vote and grow

    char txt_buf[3*MAX_SEQ_LEN];
    unsigned char bin_buf[4+MAX_SEQ_LEN/N_SEQ_BYTE];
    std::list<vote_box> consensus; 
};

#endif
