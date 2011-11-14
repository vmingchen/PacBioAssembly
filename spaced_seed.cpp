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
#include	<math.h>
#include	<vector>
#include	<deque>
#include	<list>
#include	<hash_map>
#include	"binary_parser.h"

#define DBG

// number of sequnce in a word
#define N_SEQ_WORD 16
// number of sequnce in a byte
#define N_SEQ_BYTE 4

#define N_SEGMENT 50
#define N_TRIAL 10
#define MAX_PAT_LEN 16
#define MAX(x, y) (x > y ? x : y)
#define MIN(x, y) (x < y ? x : y)
#define TR(x) ((x == 'A') ? 0 : ((x == 'C') ? 1 : (x == 'G' ? 2 : 3)))
#define handle_error(msg) do { perror(msg); exit(EXIT_FAILURE); } while (0)
#define get_seq_len(x) (*((unsigned *)x))

#ifdef DBG
#define LOG(...) fprintf(stderr, __VA_ARGS__)
size_t _ntrials = 0;
size_t _nmatches = 0;
#else
#define LOG(...) 
#endif

enum Op {
    MATCH = 0, INSERT, DELETE 
};

struct Vote {
    unsigned char A, C, G, T;
    Vote() : A(0), C(0), G(0), T(0) {};
    Vote(char c) { add(TR(c)); };
    char get() {
        if (A > C && A > G && A > T)
            return 'A';
        else if (C > A && C > G && C > T)
            return 'C';
        else if (G > A && G > C && G > T)
            return 'G';
        else
            return 'T';
    }
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

typedef unsigned char t_bseq;

// spaced seed
unsigned seed = 0;

// buf for binary DNA sequence
t_bseq *buf = NULL;
// indices for binary DNA sequence
std::list<size_t> indices;  

// length of reference DNA sequence
size_t ref_len;
t_bseq *ref_bin = NULL;
char *ref_txt = NULL;

// seedmap for reference sequence
__gnu_cxx::hash_map< unsigned, std::list<size_t> > seedmap(1<<20);
// vote against reference sequence
std::deque<Vote> votes;
size_t ref_beg;
size_t ref_end;

typedef std::list<size_t>::iterator ls_it;
typedef __gnu_cxx::hash_map< unsigned, std::list<size_t> >::iterator sm_it;


/* 
 * ===  FUNCTION  ============================================================
 *         Name:  print substr
 *  Description:  
 * ===========================================================================
 */
    void
print_substr ( const char *pseq, int beg, int end )
{
    char fmt[20];

    sprintf(fmt, "%%.%ds\n", end-beg);
    printf(fmt, pseq+beg);
}		/* -----  end of function print substr  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  set_ref
 *  Description:  set reference using binary (type 1), 
 *  text (type 2), 
 *  or votes (type 3)
 * ===========================================================================
 */
    void
set_ref ( void *new_ref, int type )
{
    if (type == 1) {
        if (ref_bin != NULL) 
            free(ref_bin);
        ref_bin = (t_bseq*)new_ref;
        ref_len = get_seq_len(ref_bin);
        if (ref_txt != NULL)
            free(ref_txt);
        if ((ref_txt = (char*)malloc(ref_len + 1)) == NULL)
            handle_error("fail to alloc memory for ref_txt");
        assert(binary_parser::bin2text(ref_bin, ref_txt, ref_len+1) == ref_len); 
    } else if (type == 2) {
        if (ref_txt != NULL)
            free(ref_txt);
        ref_txt = (char*)new_ref;
        ref_len = strlen(ref_txt);
        if (ref_bin != NULL)
            free(ref_bin);
        size_t blen = (ref_len+4-1)/4 + sizeof(unsigned);
        if ((ref_bin = (t_bseq*)malloc(blen)) == NULL)
            handle_error("fail to alloc memory for ref_bin");
        assert(binary_parser::text2bin(ref_txt, ref_bin, blen) == blen);
    } else if (type == 3) {
        char *ptxt = (char*)malloc(votes.size()+1);
        if (ptxt == NULL)
            handle_error("fail to alloc memory for string from ptxt");
        ptxt[votes.size()] = '\0';
        for (size_t i = 0; i < votes.size(); ++i)
            ptxt[i] = votes[i].get();
        ref_beg = 0;
        set_ref(ptxt, 2);
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
 *         Name:  build_seedmap
 *  Description:  build (rebuild) seedmap for reference sequence
 * ===========================================================================
 */
    void
build_seedmap (  )
{
    unsigned tseg = 0;
    unsigned nseed = ref_len - N_SEQ_WORD;

    LOG("nseed: %d\n", nseed);
    seedmap.clear();

    for (size_t i = 0; i < nseed; ++i) {
        tseg = 0;
        unsigned char *t = (unsigned char*)&tseg;
        char *p = ref_txt + i;
        for (size_t j = 0; j < N_SEQ_WORD; ++j, ++p) {
            *t = (*t << 2) | TR(*p);
            if ((j & 0x3) == 0x3) ++t;
        }
//        LOG("%08x\n", tseg);
        if (seed & tseg) seedmap[seed & tseg].push_back(i);
    }

//    for (size_t i = N_SEQ_WORD; i < ref_len; ++i) {
//        size_t j = ((i & 0xf) << 1);
//        size_t tmp = (*pseg << j) | (*(pseg+1) >> (32-j));
//        size_t key = seed & tmp;
//        LOG("%08x\n", tmp);
//        // there are a lot of 'AAAAAAAAAAAAAAAA' segments, ignore them
//        if (key) seedmap[key].push_back(i);
//        if ((i & 0xf) == 0xf) ++pseg;
//    }
}		/* -----  end of function build_seedmap  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  match_point
 *  Description:  
 * ===========================================================================
 */
    int
match_point ( unsigned sv, int base, int dist )
{
    int mp = 0;
    int diff = ref_len;
    int nd;
    sm_it it = seedmap.find(sv & seed);
    if (it == seedmap.end()) 
        return -1;
    for (ls_it it1 = it->second.begin(); it1 != it->second.end(); ++it1) {
        if ((nd = abs(*it1 - base - dist)) < diff) {
            diff = nd;
            mp = *it1;
        }
    }
    return mp;
}		/* -----  end of function match_point  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  align_seg
 *  Description:  
 * ===========================================================================
 */
    bool
align_seg ( const char *txt, size_t sb, size_t se, size_t rb, size_t re )
{
    sb = MIN(sb, se);
    se = MAX(sb, se);
    rb = MIN(rb, re);
    re = MAX(rb, re);

    int slen = se - sb + 1;
    int rlen = re - rb + 1;

    const char *ptxt = txt + sb;
    const char *pref = ref_txt + rb;

#ifdef DBG
    LOG("segment(%d, %d), reference(%d, %d)\n", sb, slen, rb, rlen);
    print_substr(txt, sb, se);
    print_substr(ref_txt, rb, re);
#endif

    /* 
    std::vector< std::vector<int> > m(slen, std::vector<int>(rlen, ref_len)); 

    for (int i = 0; i < slen; ++i) m[i][0] = i;
    for (int j = 0; j < rlen; ++j) m[0][j] = j;

    for (int i = 1; i < slen; ++i) {
        for (int j = 1; j < rlen; ++j) {
            int vm = m[i-1][j-1] + (1 - (ptxt[i-1] == pref[j-1]));
            int vd = m[i][j-1] + 1;
            int vi = m[i-1][j] + 1;
            if (vm < vd && vm < vi) {
                m[i][j] = vm;
            } 
            if (vd < vm && vd < vi) {
                m[i][j] = vd;
            }
            if (vi < vd && vi < vm) {
                m[i][j] = vi;
            }
        }
    }

    if (m[slen-1][rlen-1] > 0.3*rlen) 
        return false;

    int x = slen - 1;
    int y = rlen - 1;
    while (x > 0 && y > 0) {
        if (x > 0 && y > 0 && m[x][y] == 
                (m[x-1][y-1] + 1 - (ptxt[x-1] == pref[y-1]))) { 
            votes[ref_beg + rb + y - 1].add(ptxt[x-1]);
            --x;
            --y;
        } else if (x > 0 && m[x][y] == (m[x-1][y] + 1)) { 
            --x;
        } else {
            --y;
        }
    }
    */

    return true;
}		/* -----  end of function align_seg  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  align
 *  Description:  
 * ===========================================================================
 */
    bool
align ( t_bseq *seq, int sb, int rb, int dir )
{
    int seq_len = get_seq_len(seq);
    int seq_len_in_word = seq_len / N_SEQ_WORD;
    unsigned *pseg = (unsigned*)(seq + sizeof(unsigned));
    char *ptxt = NULL;
    bool match = false;
    int se = sb, re;

    assert((ptxt = (char*)malloc(seq_len + 1)) != NULL);
    assert(binary_parser::bin2text(seq, ptxt, seq_len+1) == seq_len);

#ifdef DBG
    unsigned *pref = (unsigned*)(ref_bin + sizeof(unsigned));
    LOG("seg: %08x, ref: %08x%08x\n", pseg[sb]&seed, pref[rb>>4], pref[(rb>>4)+1]);
#endif

    // find segment of seq, that match segment of reference
    int dist = 0;
    do {
        se += dir;
        re = match_point(pseg[se], rb, (dist<<4)*dir);
        if (re != -1 && align_seg(ptxt, sb<<4, se<<4, rb, re)) { 
            sb = se; 
            rb = re;
            match = true;
        }
    } while ((++dist<N_SEGMENT || match) && se > 0 && se < seq_len_in_word);

    /*
    if (match) {
        if (dir == 1) { 
            // prefix of seq match suffix of reference
            if ((seq_len-(sb<<4)) > (ref_len-rb)) {  
                align_seg(ptxt, sb<<4, ref_len-rb+sb, rb, ref_len);
                for (int i = ref_len-rb+sb; i < seq_len; ++i)
                    votes.push_back(Vote(ptxt[i]));
            } else {
                align_seg(ptxt, sb<<4, seq_len, rb, rb+seq_len-(sb<<4));
            }
        } else { 
            // suffix of seq match prefix of reference
            if (sb > rb) {
                align_seg(ptxt, seq_len-rb, seq_len, 0, sb+1);
                for (int i = 0; i < seq_len-rb; ++i)
                    votes.push_front(Vote(ptxt[i]));
                ref_beg += seq_len-rb;
            } else {
                align_seg(ptxt, 0, (sb+1)<<4, rb - ((sb+1)<<4), rb);
            }
        }
    }
    */

    free(ptxt);
    return match;
}		/* -----  end of function align  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  try_align
 *  Description:  
 * ===========================================================================
 */
    bool
try_align ( t_bseq *seq, std::list<size_t> &cand, size_t pos, int dir)
{
#ifdef DBG
    ++_ntrials;
#endif
    for (ls_it it = cand.begin(); it != cand.end(); ++it) {
        if (align(seq, pos, *it, dir)) 
            return true;
    }
    return false;
}		/* -----  end of function try_align  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  open_binary
 *  Description:  open binary sequence file, read into buf, build index of all
 *  sequences into indices, and return the index of the longest sequence. 
 * ===========================================================================
 */
    size_t
open_binary ( const char *fname, std::list<size_t> &indices )
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
    if ((buf = (t_bseq*)mmap(NULL, len, PROT_READ, MAP_PRIVATE, fd, 0)) == MAP_FAILED) 
        handle_error("mmap");

    for (size_t i = 0; i < len; ) {
        indices.push_back(i);
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
    size_t i_max_len, seq_len;
    t_bseq *seq;
    unsigned *pseg;

    if (argc < 3) 
        handle_error("usage: spaced_seed bin seed\n");

    // read binary sequence file and build index for all DNA sequences
    i_max_len = open_binary(argv[1], indices);
    // set the longest DNA sequence as initial reference
    set_ref(buf+i_max_len, 1);
    LOG("i_max_len: %d\n", i_max_len);
    LOG("indices: size %d\n", indices.size());
    LOG("ref_len: %d\n", ref_len);

    // parse spaced seed
    seed = parse_pattern(argv[2]);
    LOG("seed: %x\n", seed);

    // build seedmap for reference sequences (the longest sequence)
    build_seedmap();
    LOG("seedmap done\n");

    // find repeat
    int count = 0;
    for (ls_it it = indices.begin(); it != indices.end(); ++it) {
        if (*it == i_max_len) continue;
        seq = buf + *it;
        pseg = (unsigned*)(seq + sizeof(unsigned));
        seq_len = (*((unsigned*)seq) >> 4); // seq length in 32-word
        if (seq_len == 0) {
            LOG("segment too short! length: %d\n", get_seq_len(seq));
            continue;
        }
        bool found = false;
        // number of trial 
        for (size_t j = 0; j < MIN(seq_len, N_TRIAL); ++j) {
            sm_it it = seedmap.find(pseg[j] & seed);
//            LOG("%08x\n", pseg[j]);
            if (it != seedmap.end() && try_align(seq, it->second, j, 1)) { 
#ifdef DBG
                ++_nmatches;
#endif
                found = true;
                break;
            }
            it = seedmap.find(pseg[seq_len-j-1] & seed);
            if (it != seedmap.end() && try_align(seq, it->second, 
                        seq_len-j-1, -1)) { 
#ifdef DBG
                ++_nmatches;
#endif
                found = true;
                break;
            }
        }
        if (found) LOG("found %d\n", count+1);
//        else ++it;
//        it = found ? indices.erase(it) : ++it;
        if (!(++count & 0xFFF)) LOG("%d sequences processed\n", count);
    }

#ifdef DBG
    LOG("#trials: %d\n", _ntrials);
    LOG("#matches: %d\n", _nmatches);
#endif

    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
