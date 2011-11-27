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
} cell;

typedef struct {
    short op;
    char val;
} edit;

class seq_aligner {
public:
    seq_aligner() {};
    cell mat[MAXN][MAXM];       // DP matrix
    edit edits[MAXN + MAXM];    // edits of transform pseg to pref[beg:end]
    int l;                      // radius of diagonal stripe
    int nedit;                  // number of edits
    int align(seq_accessor *pseg, seq_accessor *pref) {
        int n = pseg->length();
        l = 1 + (int)(n * R);               // largest distance allowed
        int m = MIN(pref->length(), n + l);

        if (n >= MAXN || m >= MAXM) {
            LOG("segment too long: %d\n", n);
            return -1;
        }

        init_cell(n+1, m+1);

        if (!search(pseg, pref, n, m, l)) return -1;

        int ti, tj;
        goal_cell(&ti, &tj, n+1, m+1);
        if (ti < n*(1-R)) return -1;
        nedit = 0;
        find_path(ti, tj, pref);

        return tj;
    };
    int get_score(int i, int j) {
        return mat[i][j-i+l].score;
    };
    void set_score(int i, int j, int v) {
        mat[i][j-i+l].score = v;
    };
    int get_parent(int i, int j) {
        return mat[i][j-i+l].parent;
    }
    void set_parent(int i, int j, int p) {
        mat[i][j-i+l].parent = p;
    }
private:
    int match(char c, char d) { return c == d ? 5 : -4; };
    int indel(char c) { return -4; };
    // translate index into rectangle matrix to index into diagonal stripe
    void init_cell(int n, int m) {
        memset(mat, 0, sizeof(cell)*MAXN*MAXM);
//        for (int i=0; i<n; ++i) { 
//            for (int j=0; j<m; ++j) { 
//                mat[i][j].parent = 0;
//                mat[i][j].score = 0;
//            }
//        }
    };
    bool search(seq_accessor *pseg, seq_accessor *pref, int n, int m, int l) {
        pseg->reset(0);     // start from the first
        for (int i=1; i<=n; ++i) {
            char c = pseg->next();
            int beg = MAX(1, i - l);
            int end = MIN(m, i + l);
            int best_score = 0;
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
    void goal_cell(int *pti, int *ptj, int n, int m) {
        *pti = *ptj = 1;
        for (int i=1; i<=n; ++i) {
            for (int j = 1; j<=m; ++j) {
                if (mat[i][j].score > mat[*pti][*ptj].score) {
                    *pti = i;
                    *ptj = j;
                }
            }
        }
        *ptj += *pti -l;
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
};

#endif
