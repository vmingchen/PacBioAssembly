/*
 * ===========================================================================
 *
 *       Filename:  visual_align.cpp
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/28/2011 01:34:45 PM
 *
 *    Description:  visualize alignment
 *
 *       Revision:  none
 *
 * ===========================================================================
 */

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<string>
#include	<vector>
#include	<iostream>

#include	"dna_seq.h"
#include	"seq_aligner.h"
#include	"common.h"


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
    std::string ref_str;
    std::string seg_str;

    seq_aligner *paligner = new seq_aligner();
    while (cin >> ref_str >> seg_str) {
        seq_accessor ref((char*)ref_str.c_str(), true, ref_str.length());
        seq_accessor seg((char*)seg_str.c_str(), true, seg_str.length());
        if (paligner->align(&seg, &ref) <= 0) {
            cerr << "cannot align" << endl;
            cerr << ref_str << endl << seg_str << endl;
        }
        cout << paligner->final_cost() << endl;
        vector<char> aref;
        vector<char> aseg;
        int iref = 0;
        int iseg = 0;
        for (int i = 0; i < paligner->nedit; ++i) {
            if (paligner->edits[i].op == MATCH) {
                aref.push_back(ref_str[iref++]);
                aseg.push_back(seg_str[iseg++]);
            } else if (paligner->edits[i].op == INSERT) {
                aseg.push_back('-');
                aref.push_back(ref_str[iref++]);
            } else {
                aref.push_back('-');
                aseg.push_back(seg_str[iseg++]);
            }
        }
        
        for (vector<char>::iterator it = aref.begin(); it != aref.end(); ++it)
            cout << *it;
        cout << endl;
        for (vector<char>::iterator it = aseg.begin(); it != aseg.end(); ++it)
            cout << *it;
        cout << endl;

    }
    delete paligner;
    
    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
