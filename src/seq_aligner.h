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
#include	<algorithm>
#include	"dna_seq.h"
#include	"common.h"

// ratio is the maximum bias between segment and reference
//#define DEBUG_ALIGNER
#define R 0.4
#define MAXN 20000
#define MAXM (int)(MAXN*R)

enum OP {
    MATCH = 1,
    INSERT,
    DELETE
};

typedef struct {
    int cost;        
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
            ref_len = std::min(pref->length(), seg_len + max_dst);
        } else {
            ref_len = pref->length();
            max_dst = 1 + (int)(ref_len * R);
            seg_len = std::min(pseg->length(), ref_len + max_dst);
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
    int get_cost(int i, int j) { return mat[i][j-i+max_dst].cost; };
    void set_cost(int i, int j, int v) { mat[i][j-i+max_dst].cost = v; };
    int get_parent(int i, int j) { return mat[i][j-i+max_dst].parent; }
    void set_parent(int i, int j, int p) { mat[i][j-i+max_dst].parent = p;}
    int final_cost() { return get_cost(seg_ml, ref_ml); }
private:
    int match(char c, char d) { return c != d; };
    int indel(char c) { return 1; };
    // translate index into rectangle matrix to index into diagonal stripe
    void init_cell() {
        for (int i=0; i<=seg_len; ++i) {
            int beg = std::max(1, i - max_dst);
            int end = std::min(ref_len, i + max_dst);
            for (int j=beg; j<=end; ++j) {
                set_cost(i, j, std::max(i, j));
                set_parent(i, j, i>j ? DELETE : (i==j ? MATCH : INSERT));
            } 
        }
        set_cost(0, 0, 0);
        set_parent(0, 0, 0);
    };
    bool search(seq_accessor *pseg, seq_accessor *pref) {
        pseg->reset(0);     // start from the first
        for (int i=1; i<=seg_len; ++i) {
            int best_cost = seg_len;
            char c = pseg->next();
            int beg = std::max(1, i - max_dst);
            int end = std::min(ref_len, i + max_dst);
            pref->reset(beg-1);     // start from the beg-th element
            for (int j=beg; j<=end; ++j) {
                char d = pref->next();
                int t; 
                int cost = get_cost(i, j);
                if ((t = get_cost(i-1, j-1) + match(c, d)) < cost) {
                    set_cost(i, j, t);
                    set_parent(i, j, MATCH);
                    cost = t;
                }
                if (i-j<max_dst && (t=get_cost(i,j-1)+indel(c)) < cost) {
                    set_cost(i, j, t);
                    set_parent(i, j, INSERT);
                    cost = t;
                }
                if (j-i<max_dst && (t=get_cost(i-1,j)+indel(d)) < cost) {
                    set_cost(i, j, t);
                    set_parent(i, j, DELETE);
                    cost = t;
                }
                best_cost = std::min(cost, best_cost);
            }
#ifdef DEBUG_ALIGNER
            LOG("i = %d, best_cost = %d\n", i, best_cost);
#endif
            // early failure 
            if (i > 10 && get_cost(i, i) > i*R) {
//                printf("failed at %d\n", i);
                return false;
            }
        }
        return true;
    };
    void goal_cell() {
        seg_ml = seg_len;
        ref_ml = ref_len;
        while (seg_ml > ref_len 
                && get_cost(seg_ml-1, ref_ml) <= get_cost(seg_ml, ref_ml))
            --seg_ml;
        while (ref_ml > seg_len 
                && get_cost(seg_ml, ref_ml-1) <= get_cost(seg_ml, ref_ml))
            --ref_ml;
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
                    printf("%d\t", get_cost(i, j));
            }
            printf("\n");
        }
    }
};

#endif
