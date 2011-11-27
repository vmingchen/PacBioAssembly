/*
 * ===========================================================================
 *
 *       Filename:  quality.cpp
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/16/2011 01:20:42 PM
 *
 *    Description:  generate mean quality for a DNA sequence
 *
 *       Revision:  none
 *
 * ===========================================================================
 */

#include	<iostream>
#include	<string.h>
#include	<stdlib.h>

const size_t MAXLEN = 200000;

using namespace std;
/* 
 * ===  FUNCTION  ============================================================
 *         Name:  main
 *  Description:  
 * ===========================================================================
 */
    int
main ( int argc, char *argv[] )
{ 
    char line[MAXLEN];
    while (cin.getline(line, MAXLEN)) {
        unsigned quality = 0;
        size_t len = strlen(line);
        for (size_t i = 0; i < len; i++) {
            quality += line[i];
        }
        cout << (quality/len) << endl;
    }
    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
