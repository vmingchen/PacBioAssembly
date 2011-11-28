/*
 * ===========================================================================
 *
 *       Filename:  seq_aligner.h
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/10/2011 11:03:08 PM
 *
 *    Description:  perform hurestic local alignment
 *
 *       Revision:  none
 *
 *
 * ===========================================================================
 */

#ifndef SEQ_ALIGNER_H
#define SEQ_ALIGNER_H

#include	<string.h>
#include	<math.h>
#include	<limits.h>
#include	<vector>
#include	"dna_seq.h"
#include	"common.h"

// ratio is the maximum bias between segment and reference
//#define DEBUG_ALIGNER
#define R 0.3
#define MAXN 2000
#define MAXM (int)(MAXN*R)

enum OP {
    MATCH = 1,
    INSERT,
    DELETE
};

typedef struct {
    int score;        
    int parent;      
} state;

typedef struct {
    short op;
    char val;
} edit;

class seq_aligner {
public:
    seq_aligner() {};
    int seg_len;                // max possible length of match in segment
    int ref_len;                // max possible length of match in ref
    int max_dst;                // max distance allowed
    int seg_ml;                 // length of match in seg
    int ref_ml;                 // length of match in ref
    state mat[MAXN][MAXM];      // DP matrix
    int nedit;                  // number of edits
    edit edits[MAXN + MAXM];    // edits of transform pseg to pref[beg:end]
    int align(seq_accessor *pseg, seq_accessor *pref) {
        // work out parameters
        if (pref->length() >= pseg->length()) { 
            seg_len = pseg->length();
            max_dst = 1 + (int)(seg_len * R);
            ref_len = MIN(pref->length(), seg_len + max_dst);
        } else {
            ref_len = pref->length();
            max_dst = 1 + (int)(ref_len * R);
            seg_len = MIN(pseg->length(), ref_len + max_dst);
        }

        if (seg_len >= MAXN || max_dst >= MAXM) {
            LOG("segment too long: %d\n", seg_len);
            return -1;
        }

        init_cell();

        if (!search(pseg, pref)) return -1;

        goal_cell();
        if (ref_ml < ref_len*(1-R)) return -1;
        nedit = 0;
        find_path(seg_ml, ref_ml, pref);

#ifdef DEBUG_ALIGNER
        print_matrix(pseg, pref);
        printf("(%d, %d)\n", seg_ml, ref_ml);
        printf("seg_len: %d, ref_len: %d\n", seg_len, ref_len);
#endif

        return ref_ml;
    };
    int get_score(int i, int j) {
        return mat[i][j-i+max_dst].score;
    };
    void set_score(int i, int j, int v) {
        mat[i][j-i+max_dst].score = v;
    };
    int get_parent(int i, int j) {
        return mat[i][j-i+max_dst].parent;
    }
    void set_parent(int i, int j, int p) {
        mat[i][j-i+max_dst].parent = p;
    }
    int final_score() {
        return get_score(seg_ml, ref_ml);
    }
private:
    int match(char c, char d) { return c == d ? 5 : -4; };
    int indel(char c) { return -4; };
    // translate index into rectangle matrix to index into diagonal stripe
    void init_cell() {
        for (int i=0; i<=seg_len; ++i) {
            int beg = MAX(1, i - max_dst);
            int end = MIN(ref_len, i + max_dst);
            for (int j=beg; j<=end; ++j) {
                set_score(i, j, -4*MAX(i, j));
                set_parent(i, j, i>j ? DELETE : (i==j ? MATCH : INSERT));
            }
        }
        set_score(0, 0, 0);
        set_parent(0, 0, 0);
    };
    bool search(seq_accessor *pseg, seq_accessor *pref) {
        pseg->reset(0);     // start from the first
        int best_score = 0;
        for (int i=1; i<=seg_len; ++i) {
            char c = pseg->next();
            int beg = MAX(1, i - max_dst);
            int end = MIN(ref_len, i + max_dst);
            pref->reset(beg-1);     // start from the beg-th element
            for (int j=beg; j<=end; ++j) {
                char d = pref->next();
                int t; 
                int score = get_score(i, j);
                if ((t = get_score(i-1, j-1) + match(c, d)) > score) {
                    set_score(i, j, t);
                    set_parent(i, j, MATCH);
                    score = t;
                }
                if ((t = get_score(i, j-1) + indel(c)) > score) {
                    set_score(i, j, t);
                    set_parent(i, j, INSERT);
                    score = t;
                }
                if ((t = get_score(i-1, j) + indel(d)) > score) {
                    set_score(i, j, t);
                    set_parent(i, j, DELETE);
                    score = t;
                }
                best_score = MAX(score, best_score);
            }
            // early failure 
            if (i > 10 && best_score < 3*i) return false;
        }
        return true;
    };
    void goal_cell() {
        seg_ml = ref_ml = 0;
        int best_score = get_score(0, 0);
        for (int i=1; i<=seg_len; ++i) {
            int beg = MAX(1, i - max_dst);
            int end = MIN(ref_len, i + max_dst);
            for (int j=beg; j<=end; ++j) {
                if (get_score(i, j) > best_score) {
                    seg_ml = i;
                    ref_ml = j;
                    best_score = get_score(i, j);
                }
            }
        }
    };
    void find_path(int i, int j, seq_accessor *pref) {
        int p = get_parent(i, j);
        if (p == MATCH) {
            find_path(i-1, j-1, pref);
            edits[nedit].op = MATCH;
            edits[nedit].val = pref->at(j-1);
            ++nedit;
        }
        if (p == INSERT) {
            find_path(i, j-1, pref);
            edits[nedit].op = INSERT;
            edits[nedit].val = pref->at(j-1);
            ++nedit;
        } 
        if (p == DELETE) {
            find_path(i-1, j, pref);
            edits[nedit].op = DELETE;
            ++nedit;
        }
    };
    void print_matrix(seq_accessor *pseg, seq_accessor *pref) {
        printf(" \t \t");
        for (int j = 1; j <= ref_len; ++j) {
            printf("%c\t", pref->at(j-1));
        }
        printf("\n");
        for (int i = 0; i <= seg_len; ++i) {
            printf("%c\t", i == 0 ? ' ' : pseg->at(i-1));
            for (int j = 0; j <= ref_len; ++j) {
                if (j < i-max_dst || j > i+max_dst) 
                    printf("%d\t", -1);
                else
                    printf("%d\t", get_score(i, j));
            }
            printf("\n");
        }
    }
};

#endif
