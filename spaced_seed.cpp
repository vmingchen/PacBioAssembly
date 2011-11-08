/*
 * ===========================================================================
 *
 *       Filename:  spaced_seed.cpp
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/07/2011 03:03:04 PM
 *
 *    Description:  
 *
 *       Revision:  none
 *
 * ===========================================================================
 */

#include	<stdlib.h>
#include	<map>

#define MAX_PAT_LEN 16
#define TR(x) ((x == 'A') ? 0 : ((X == 'C') ? 1 : (X == 'G' ? 2 : 3)))

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  parse_pattern
 *  Description:  
 * ===========================================================================
 */
    unsigned
parse_pattern ( const char *pat )
{
    unsigned pattern = 0;
    for (int i = 0; i < strlen(pat) && i < MAX_PAT_LEN; ++i) {
        pattern = (pat[i] == '1') ? ((pattern << 2) | 0x3) : (pattern << 2);
    }
    return pattern;
}		/* -----  end of function parse_pattern  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  seed
 *  Description:  
 * ===========================================================================
 */
    void
seed ( const char *pstr )
{
    return ;
}		/* -----  end of function seed  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  main
 *  Description:  
 * ===========================================================================
 */
    int
main ( int argc, char *argv[] )
{ 
    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
