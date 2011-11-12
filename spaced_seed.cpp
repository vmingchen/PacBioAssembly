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
#include    <sys/types.h>
#include    <sys/stat.h>
#include    <sys/mman.h>
#include    <fcntl.h>
#include    <unistd.h>
#include	<assert.h>
#include	<string.h>
#include	<stdio.h>
#include	<vector>
#include	<deque>
#include	<map>
#include	"binary_parser.h"

// number of sequnce in a word
#define N_SEQ_WORD 16
// number of sequnce in a byte
#define N_SEQ_BYTE 4

#define N_TRIAL 20
#define MAX_PAT_LEN 16
#define TR(x) ((x == 'A') ? 0 : ((X == 'C') ? 1 : (X == 'G' ? 2 : 3)))
#define handle_error(msg) do { perror(msg); exit(EXIT_FAILURE); } while (0)

struct Vote {
    unsigned char A, C, G, T;
    Vote() : A(0), C(0), G(0), T(0) {};
    void add(int c) {
        if (c == 0) 
            A = (A == 255 ? 255 : A+1);
        else if (c == 1)
            C = (C == 255 ? 255 : C+1);
        else if (c == 2)
            G = (G == 255 ? 255 : G+1);
        else
            T = (T == 255 ? 255 : T+1);
    };
}; 

// spaced seed
unsigned seed = 0;

// buf for binary DNA sequence
unsigned char *buf = NULL;
// indices for binary DNA sequence
std::list<size_t> indices;  

// length of reference DNA sequence
size_t ref_len;
unsigned char *ref_bin = NULL;
char *ref_txt = NULL;

// seedmap for reference sequence
std::map< unsigned, list<size_t> > seedmap;
// vote against reference sequence
std::deque<Vote> votes;


/* 
 * ===  FUNCTION  ============================================================
 *         Name:  set_ref
 *  Description:  set reference using binary (type 1) or text (type 2)
 * ===========================================================================
 */
    bool
set_ref ( void *new_ref, int type )
{
    if (type == 1) {

    } else (type == 2) {

    }
    return ;
}		/* -----  end of function set_ref  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  parse_pattern
 *  Description:  parse spaced seed pattern to a mask
 * ===========================================================================
 */
    unsigned
parse_pattern ( const char *pat )
{
    unsigned pattern = 0;
    for (size_t i = 0; i < strlen(pat) && i < MAX_PAT_LEN; ++i) {
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
 *         Name:  build_seedmap
 *  Description:  build (rebuild) seedmap for reference sequence
 * ===========================================================================
 */
    void
build_seedmap (  )
{
    unsigned *pseg = (unsigned*)(ref_bin + sizeof(unsigned)); 
    unsigned nseed = ref_len - N_SEQ_WORD;

    seedmap.clear();
    for (size_t i = 0; i < nseed; ++i) {
        size_t j = ((i & 0xf) << 1);
        size_t tmp = (*pseg << j) & (*(pseg+1) >> (32-j));
        seedmap[(seed & tmp)].insert(i);
    }
}		/* -----  end of function build_seedmap  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  try_align
 *  Description:  
 * ===========================================================================
 */
    void
try_align ( unsigned char *seq, std::list<size_t> &list, size_t pos, int dir)
{

//            if (it != indices.end()) {
//                list<size_t>::iterator beg = it->second.begin();
//                list<size_t>::iterator end = it->second.end();
//                while (beg != end) {
//                    if (align(seq, *beg, j)) break;
//                    ++beg;
//                }
//                break;
//            }
    return ;
}		/* -----  end of function try_align  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  open_binary
 *  Description:  open binary sequence file, read into buf, build index of all
 *  sequences into indices, and return the index of the longest sequence. 
 * ===========================================================================
 */
    size_t
open_binary ( const char *fname, void *buf, std::list<size_t> &indices )
{
    struct stat fst;
    size_t len;
    size_t i_max_len = 0, max_len = 0;
    int fd; 
    
    if ((fd = open(fname, O_RDONLY)) == -1)
        handle_error("open");

    if (fstat(fd, &fst) == -1) 
        handle_error("fstat");

    len = fst.st_size;
    if ((buf = mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, 0)) == MAP_FAILED) 
        handle_error("mmap");

    for (size_t i = 0; i < len; ) {
        indices.insert(i);
        size_t seq_len = *((unsigned*)(buf + i));
        if (seq_len > max_len) {
            max_len = seq_len;
            i_max_len = i;
        }
        i += sizeof(unsigned) + (seq_len + 4 - 1)/4;
    }

    return i_max_len;
}		/* -----  end of function open_binary  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  main
 *  Description:  
 * ===========================================================================
 */
    int
main ( int argc, char *argv[] )
{ 
    size_t i_max_len, max_len, seq_len;
    unsigned char *seq;
    unsigned *pseg;

    typedef std::list<size_t>::iterator lit;
    typedef std::map< unsigned, list<size_t> >::iterator sm_it;

    if (argc < 3) 
        handle_error("usage: spaced_seed bin seed\n");

    // read binary sequence file and build index for all DNA sequences
    i_max_len = (unsigned char*)open_binary(argv[1], buf, indices);
    // set the longest DNA sequence as initial reference


    // parse spaced seed
    seed = parse_pattern(argv[2]);

    // build seedmap for reference sequences (the longest sequence)
    max_len = *((unsigned*)buf[i_max_len]);
    pseg = (unsigned*)(buf[i_max_len] + sizeof(unsigned)); 
    for (size_t i = 0; i < max_len; ++i) {
        size_t j = ((i & 0xf) << 1);
        size_t tmp = (*pseg << j) & (*(pseg+1) >> (32-j));
        seedmap[(seed & tmp)].insert(i);
    }

    // find repeat
    for (lit it = indices.begin(); it != indices.end(); ) {
        seq = buf[*it];
        pseg = (unsigned*)(seq + sizeof(unsigned));
        seq_len = (*((unsigned*)seq) >> 4); // seq length in 32-word
        size_t tmp = max(seq_len, N_TRIAL);
        bool found = false;
        // number of trial 
        for (size_t j = 0; j < max(seq_len, N_TRIAL); ++j) {
            sm_it it = seedmap.find(pseg[j]);
            if (it != indices.end() && try_align(seq, it->second, j, 1)) { 
                found = true;
                break;
            }
            it = seedmap.find(pseg[seq_len-j-1]);
            if (it != indices.end() && try_align(seq, it->second, 
                        seq_len-j-1, -1)) { 
                found = true;
                break;
            }
        }
        it = found ? indices.erase(it) : it+1;
    }

    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
