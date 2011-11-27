/*
 * ===========================================================================
 *
 *       Filename:  local_aligner.h
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

#ifndef LOCAL_ALIGNER_H
#define LOCAL_ALIGNER_H

#include	<string.h>
#include	<math.h>
#include	<limits.h>
#include	<vector>
#include	"dna.h"

// ratio is the maximum bias between segment and reference
#define R 0.3
#define MAXN 20000
#define MAXM (MAXN*R)

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

class local_aligner {
public:
    local_aligner() {};
    cell mat[MAXN][MAXM];       // DP matrix
    edit edits[MAXN + MAXM];    // edits of transform pseg to pref[beg:end]
    int nedit;      // number of edits
    int align(dna_accessor *pseg, dna_accessor *pref) {
        int n = pseg->length();
        int l = 1 + (int)(n * R);               // largest distance allowed
        int m = MIN(pref->length(), n + l);

        if (n >= MAXN || m >= MAXM) {
            LOG("segment too long: %d\n", n);
            return -1;
        }

        init_cell(n+1, m+1);

        if (!search(pseg, pref, n, m, l)) return -1;

        int ti, tj;
        goal_cell(&ti, &tj, n, m);
        if (ti < n*(1-ratio)) return -1;
        nedit = 0;
        find_path(ti, tj, pref);

        return tj + ti - l;
    };
private:
    int match(char c, char d) { return c == d ? 5 : -4; };
    int indel(char c) { return -4; };
    void init_cell(int n, int m) {
        memset(mat, 0, sizeof(cell)*MAXN*MAXM);
//        for (int i=0; i<n; ++i) { 
//            for (int j=0; j<m; ++j) { 
//                mat[i][j].parent = 0;
//                mat[i][j].cost = 0;
//            }
//        }
    };
    bool search(dna_accessor *pseg, dna_accessor *pref, int n, int m, int l) {
        pseg->reset(0);     // start from the first
        for (int i=1; i<=n; ++i) {
            char c = pseg->next();
            int beg = MAX(1, i - l);
            int end = MIN(m, i + l);
            int best_score = 0;
            pref->reset(beg-1);     // start from the beg-th element
            for (int k=beg; k<=end; ++k) {
                char d = pref->next();
                int j = k - i + l + 1;  // index into DP matrix
                int t; 
                if ((t = mat[i-1][j-1].score + match(c, d)) > mat[i][j].score) {
                    mat[i][j].score = t;
                    mat[i][j].parent = MATCH;
                }
                if ((t = mat[i][j-1].score + indel(c)) > mat[i][j].score) {
                    mat[i][j].score = t;
                    mat[i][j].parent = INSERT;
                }
                if ((t = mat[i-1][j].score + indel(d)) > mat[i][j].score) {
                    mat[i][j].score = t;
                    mat[i][j].parent = DELETE;
                }
                best_score = MAX(mat[i][j].score, best_score);
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
                if (mat[i][j].cost > mat[*pti][*ptj].cost) {
                    *pti = i;
                    *ptj = j;
                }
            }
        }
    };
    void find_path(int i, int j, dna_accessor *pref) {
        if (mat[i][j].parent == 0) return;
        if (mat[i][j].parent == MATCH) {
            find_path(i-1, j-1);
            edits[nedit].op = MATCH;
            edits[nedit].val = pref->at(j);
            ++nedit;
        }
        if (mat[i][j].parent == INSERT) {
            find_path(i, j-1);
            edits[nedit].op = INSERT;
            edits[nedit].val = pref->at(j);
            ++nedit;
        } 
        if (mat[i][j].parent == DELETE) {
            find_path(i-1, j);
            edits[nedit].op = DELETE;
            ++nedit
        }
    };
};

#endif
