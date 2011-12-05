/*
 * ===========================================================================
 *
 *       Filename:  common.h
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/14/2011 06:08:52 PM
 *
 *    Description:  
 *
 *       Revision:  none
 *
 *
 * ===========================================================================
 */
#ifndef COMMON_H
#define COMMON_H

#include	<ext/hash_map>
#include	<algorithm>
#include	<list>

#define DBG

#ifdef DBG
#define LOG(...) fprintf(stderr, __VA_ARGS__)
#else
#define LOG(...) 
#endif

//! max length of genome allowed
#define MAX_SEQ_LEN 800000
//! max length of segment reads processed
#define MAX_READ_LEN 20000
//! max difference (distance) allowed between overlapped reads
#define MAX_DIFF_LEN 6000
//! max ratio of difference (distance)
#define MAXR 0.3
//! min length of aligned region to justify overlap
#define OVERLAP_MIN 32

typedef unsigned t_seed;
typedef unsigned char t_bseq;
typedef __gnu_cxx::hash_map< unsigned, std::list<int> > hash_table;
typedef __gnu_cxx::hash_map< unsigned, std::list<int> >::iterator sm_it;

#endif
