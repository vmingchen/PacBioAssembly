/*
 * ===========================================================================
 *
 *       Filename:  stat_hash.c
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/07/2011 01:36:14 PM
 *
 *    Description:  generate sequence hash using base statistics
 *
 *       Revision:  none
 *
 * ===========================================================================
 */


#include	<stdlib.h>
#include	<stdio.h>

#define QUANTIZE(a) (((a>>4) > 0xff) ? 0xff : ((a>>4) & 0xff))
#define HASH(a, b, c, d) ((a << 24) | (b << 16) | (c << 8) | d)

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  main
 *  Description:  
 * ===========================================================================
 */
    int
main ( int argc, char *argv[] )
{ 
    unsigned hash = 0; 
    unsigned a = 0, c = 0, g = 0, t = 0;
    char i; 
    while ((i = getchar()) != EOF) {
        if (i == '\n') {
            hash = HASH(QUANTIZE(a), QUANTIZE(c), QUANTIZE(g), QUANTIZE(t));
            printf("%08x\n", hash);
            a = c = g = t = 0;
        }
        if (i == 'A') ++a;
        else if (i == 'C') ++c;
        else if (i == 'G') ++g;
        else if (i == 'T') ++t;
        else if (i != 0x0a) { printf("Error character: %02x\n", i); } 
    }
    hash = HASH(QUANTIZE(a), QUANTIZE(c), QUANTIZE(g), QUANTIZE(t));
    printf("%08x\n", hash);

    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
