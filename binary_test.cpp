/*
 * ===========================================================================
 *
 *       Filename:  binary_test.cpp
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/09/2011 11:18:52 PM
 *
 *    Description:  test binary_parser
 *
 *       Revision:  none
 *
 * ===========================================================================
 */

#include	<stdlib.h>
#include	<stdio.h>
#include	<string.h>
#include	"binary_parser.h"

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  main
 *  Description:  
 * ===========================================================================
 */
    int
main ( int argc, char *argv[] )
{ 
//    const char *pdna = "GGTCGCAATGAAGAAGAGGCGCTGGAAATCTATTTCCTGAGTTGCGCGCGGCTATTGAGATTGACATTTTATGGAGACCG";
    const size_t N = 20000;
    char pdna[N];
    unsigned char binary[4*N];
    char text[N];

//    while (fgets(pdna, N, stdin)) { 
    while (scanf("%s", pdna) != EOF) { 
        binary_parser::text2bin(pdna, binary, N);
        binary_parser::bin2text(binary, text, N);
        printf("%s\n", text);
    }

//    printf("%d\n", strcmp(pdna, text));
//    printf("%s\n%s\n", pdna, text);
    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
