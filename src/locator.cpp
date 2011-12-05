/*
 * ===========================================================================
 *
 *       Filename:  locator.cpp
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  12/04/2011 06:24:56 PM
 *
 *    Description:  locate sequence's position in contig
 *
 *       Revision:  none
 *
 * ===========================================================================
 */


#include	<stdlib.h>
#include	<stdio.h>
#include	<string.h>
#include	"common.h"
#include	"dna_seq.h"
#include	"seq_aligner.h"

#define MAX_LOC_LEN 40000
#define MAX_LOC_DIFF 6000

char contig[MAX_SEQ_LEN];
char sequence[MAX_SEQ_LEN];
hash_table seedmap(1<<23);
unsigned seed_pattern;
seq_aligner<MAX_LOC_LEN, MAX_LOC_DIFF> *paligner;

typedef std::list<int>::iterator list_it;

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  main
 *  Description:  
 * ===========================================================================
 */
    int
main ( int argc, char *argv[] )
{ 
    if (argc <= 2) {
        fprintf(stderr, "usage: locator contig_file seed < seq_file\n");
        return EXIT_FAILURE;
    }

    FILE *fp = fopen(argv[1], "r");
    fscanf(fp, "%s", contig);

    char str_pat[20];
    for (int i = 0; i < strlen(argv[2]); ++i)
        str_pat[i] = (argv[2][i] == '1' ? 'T' : 'A');
    seed_pattern = dna_seq::encode(str_pat);

    // convert N to A
    seq_accessor ac_contig(contig, true, strlen(contig));
    char *p = contig;
    for (int i = 0; i < ac_contig.length(); ++i) 
        if (*p == 'N') *p = 'A';

    for (int i = 0; i < ac_contig.length(); ++i) {
        int sd = dna_seq::encode(contig+i);
        if (sd & seed_pattern) 
            seedmap[sd & seed_pattern].push_back(i);
    }

    paligner = new seq_aligner<MAX_LOC_LEN, MAX_LOC_DIFF>(0.15);
    int nseq = 0;
    while (scanf("%s", sequence) != EOF) {
        int len = strlen(sequence);
        if (len < 500) continue;
        bool found = false;
        for (int j = 0; j < 50 && !found; ++j) {
            int seed = dna_seq::encode(sequence+j) & seed_pattern;
            sm_it sit = seedmap.find(seed);
            if (sit == seedmap.end()) continue;
            seq_accessor ac_seg(sequence+j, true, len - j);
            for (list_it it = sit->second.begin(); it != sit->second.end(); ++it) {
                seq_accessor ac_ref(contig + *it, true, 
                        ac_contig.length() - *it);
                if (paligner->align(&ac_seg, &ac_ref) > 0) {
                    found = true;
                    printf("%d\t%d\t%d\t%d\t%d\n", nseq, *it, 
                            paligner->final_cost(), len - j,
                            paligner->get_cost(len -j , len - j));
                    break;
                }
            }
        }
        if (!(++nseq & 0xFFF)) LOG("%d sequences processed\n", nseq);
    }
    LOG("totally %d sequences processed\n", nseq);

    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
