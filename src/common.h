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

#define DBG

#ifdef DBG
#define LOG(...) fprintf(stderr, __VA_ARGS__)
#else
#define LOG(...) 
#endif

#define MAX(x, y) ((x > y) ? (x) : (y))
#define MIN(x, y) ((x < y) ? (x) : (y))

#define MAX_SEQ_LEN 100000
typedef unsigned char t_bseq;

#endif
