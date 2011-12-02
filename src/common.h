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

#define MAX_SEQ_LEN 100000

typedef unsigned t_seed;
typedef unsigned char t_bseq;
typedef __gnu_cxx::hash_map< unsigned, std::list<int> > hash_table;

#endif
