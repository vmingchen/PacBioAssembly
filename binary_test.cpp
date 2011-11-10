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
    const size_t N = 20000;
    char pdna[N];
    unsigned char binary[4*N];
    char text[N];
    size_t blen, tlen = 0;
    FILE *fp = NULL;

    if (argc < 3) {
        fprintf(stderr, "usage: %s option filename\n"
                "options: 0 do in-memory checks for DNA text from stdin\n"
                "         1 read text input from stdin write to binary file\n"
                "         2 read binary file write text to stdout\n", argv[0]);
        return EXIT_FAILURE;
    }

    if (argv[1][0] == '0') {
        while (scanf("%s", pdna) != EOF) { 
            blen = binary_parser::text2bin(pdna, binary, N);
            tlen = binary_parser::bin2text(binary, text, N);
            if (strcmp(pdna, text) != 0) { 
                printf("Error:%s\n%s\n", pdna, text);
                return EXIT_FAILURE;
            }
        }
    }

    if (argv[1][0] == '1') { 
        fp = fopen(argv[2], "wb");
        while (scanf("%s", pdna) != EOF) { 
            blen = binary_parser::text2bin(pdna, binary, N);
            assert(blen == fwrite(binary, 1, blen, fp));
            tlen = binary_parser::bin2text(binary, text, N);
//            printf("%s\n", text);
        }
        fclose(fp);
    }

    if (argv[1][0] == '2') { 
        fp = fopen(argv[2], "rb");
        unsigned char *p = binary;
        while (fread(p, sizeof(unsigned), 1, fp) != 0) {
            tlen = *((unsigned*)p);
            blen = sizeof(unsigned) + (tlen+4-1)/4;
            fread(p + sizeof(unsigned), 1, blen - sizeof(unsigned), fp);
            assert(tlen == binary_parser::bin2text(p, text, N));
            printf("%s\n", text);
        }
        fclose(fp);
    }

    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
