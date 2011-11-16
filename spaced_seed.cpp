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
#include	<ext/hash_map>

#include	"binary_parser.h"
#include	"common.h"

#define MAX_SEQ_LEN 100000
#define STRONG 3
#define N_SEGMENT 100
#define N_TRIAL 20
#define MAX_PAT_LEN N_SEQ_WORD
#define handle_error(msg) do { perror(msg); exit(EXIT_FAILURE); } while (0)
#define get_seq_len(x) (*((unsigned *)x))

#ifdef DBG
size_t _ntrials = 0;
size_t _nmatches = 0;
#endif

enum Op {
    MATCH = 0, INSERT, DELETE 
};

struct Vote {
    unsigned char A, C, G, T;
    Vote() : A(0), C(0), G(0), T(0) {};
    Vote(char c) { add(C2I(c)); };
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

struct Reference {
    Reference(const t_bseq *pseq) {};
    bool align(int rb, int re, const char *pseg, int sb, int se) {
    };
    void append(char *pseg, int len) {
    };
    void prepend(char *pseg, int len) {
    };
    bool contained(int pos) {
        return pos >= pre && pos < post;
    };
    int beg, end, pre, post; 
    char _buf[3*MAX_SEQ_LEN];
    std::deque<Vote> consensus; 
};

// spaced seed
unsigned seed = 0;

// buf for binary DNA sequence
t_bseq *buf = NULL;
// indices for binary DNA sequence
std::list<int> indices;  

// length of reference DNA sequence
size_t ref_len;
t_bseq *ref_bin = NULL;
char *ref_txt = NULL;

// seedmap for reference sequence
__gnu_cxx::hash_map< unsigned, std::list<int> > seedmap(1<<20);
// vote against reference sequence
std::deque<Vote> votes;
size_t ref_beg;
size_t ref_end;

typedef std::list<int>::iterator ls_it;
typedef __gnu_cxx::hash_map< unsigned, std::list<int> >::iterator sm_it;

inline void swap(int &a, int &b) { int t = a; a = b; b = t; }
inline unsigned get_seq_len(const t_bseq *x) { return  *((unsigned *)x); }


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
    if (beg > end) swap(beg, end);
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
    char dnapat[MAX_PAT_LEN+1] = "AAAAAAAAAAAAAAAA";
    size_t len = MIN(strlen(pat), MAX_PAT_LEN);

#ifdef DBG
    if (len < MAX_PAT_LEN) 
        LOG("WARNING: pattern is shorter than %d\n", MAX_PAT_LEN);
#endif

    for (size_t i = 0; i < len; ++i) 
        dnapat[i] = (pat[i] == '1' ? 'T' : 'A');
    return binary_parser::encode(dnapat);
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
    unsigned nseed = ref_len - N_SEQ_WORD;

    LOG("nseed: %d\n", nseed);
    seedmap.clear();

    for (size_t i = 0; i < nseed; ++i) {
        unsigned tseg = binary_parser::encode(ref_txt + i);
        // there are a lot of 'AAAAAAAAAAAAAAAA' segments, ignore them
        if (seed & tseg) seedmap[seed & tseg].push_back(i);
    }
}		/* -----  end of function build_seedmap  ----- */


/* 
 * ===  FUNCTION  ============================================================
 *         Name:  filter_seq
 *  Description:  return true if pa and pb are not likely to be similar
 * ===========================================================================
 */
    bool
filter_seq ( const char *pa, int la, const char *pb, int lb )
{
    int stats[2][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}};

    for (int i = 0; i < la; ++i, ++pa) ++stats[0][C2I(*pa)];
    for (int i = 0; i < lb; ++i, ++pb) ++stats[1][C2I(*pb)];

    int diff = abs(stats[0][0] - stats[1][0]) 
        + abs(stats[0][1] - stats[1][1]) 
        + abs(stats[0][2] - stats[1][2]) 
        + abs(stats[0][3] - stats[1][3]);

    return (diff*4) > (la + lb);
}		/* -----  end of function filter_seq  ----- */

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

    if (abs(dist) > 800) return -1;

    sm_it it = seedmap.find(sv & seed);
#ifdef NDBG
    char seg[N_SEQ_WORD];
    binary_parser::decode(sv & seed, seg);
    LOG("%.16s\n", seg);
#endif
    if (it == seedmap.end()) 
        return -1;
    for (ls_it it1 = it->second.begin(); it1 != it->second.end(); ++it1) {
//        LOG("offset: %d, base: %d, dist: %d\n", *it1, base, dist);
        int nd = abs(*it1 - base - dist);
        if (8*nd <= abs(dist) && nd < diff) {
            diff = nd;
            mp = *it1;
        }
    }
//    LOG("selected: %d\n", mp);
    return diff == ref_len ? -1 : mp;
}		/* -----  end of function match_point  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  align_seg
 *  Description:  
 * ===========================================================================
 */
    bool
align_seg ( const char *txt, int sb, int se, int rb, int re, int dir)
{
#ifdef NDBG
    sb = sb + ((1 - dir) << 3);
    se = se + ((1 + dir) << 3);
    rb = rb + ((1 - dir) << 3);
    re = re + ((1 + dir) << 3);
#endif

    if (dir == -1) {
        sb += 16;
        se += 16;
        rb += 16;
        re += 16;
    }

    if (se < sb) swap(sb, se);
    if (re < rb) swap(rb, re);

    int slen = se - sb;
    int rlen = re - rb;

    const char *ptxt = txt + sb;
    const char *pref = ref_txt + rb;

#ifdef DBG
//    LOG("segment(%d, %d), reference(%d, %d)\n", sb, slen, rb, rlen);
//    print_substr(txt, sb, se);
//    print_substr(ref_txt, rb, re);
#endif

//    if (filter_seq(ptxt, slen, pref, rlen)) return false;

//    slen += 1;
//    rlen += 1;
//    std::vector< std::vector<int> > m(slen, std::vector<int>(rlen, ref_len)); 
//
//    for (int i = 0; i < slen; ++i) m[i][0] = i;
//    for (int j = 0; j < rlen; ++j) m[0][j] = j;
//
//    for (int i = 1; i < slen; ++i) {
//        for (int j = 1; j < rlen; ++j) {
//            int vm = m[i-1][j-1] + (1 - (ptxt[i-1] == pref[j-1]));
//            int vd = m[i][j-1] + 1;
//            int vi = m[i-1][j] + 1;
//            if (vm < vd && vm < vi) {
//                m[i][j] = vm;
//            } 
//            if (vd < vm && vd < vi) {
//                m[i][j] = vd;
//            }
//            if (vi < vd && vi < vm) {
//                m[i][j] = vi;
//            }
//        }
//    }
//
//    if (m[slen-1][rlen-1] > 0.3*rlen) 
//        return false;

    /* 
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
 *  Description:  ap_s, aligned point in segment, against ap_r in reference
 * ===========================================================================
 */
    bool
align ( t_bseq *seq, int ap_s, int ap_r, int dir )
{
    int seq_len = get_seq_len(seq);
    int seq_len_in_word = seq_len / N_SEQ_WORD;
    unsigned *pseg = (unsigned*)(seq + sizeof(unsigned));
    char *ptxt = NULL;
    bool mismatch = false;
    size_t n_matched_seg = 0;       
    size_t l_matched_seg = 0;       // length of aligned segments in word
    int sb, rb, se, re; 

    assert((ptxt = (char*)malloc(seq_len + 1)) != NULL);
    assert(binary_parser::bin2text(seq, ptxt, seq_len+1) == seq_len);

    // find segment of seq, that match segment of reference
    se = sb = ap_s;
    rb = ap_r;
    int dist = 0;
    do {
        ++dist;
        se += dir;
        re = match_point(pseg[se], rb, (dist<<4)*dir);
#ifdef NDBG
        if (re != -1) { 
            print_substr(ptxt, (sb<<4), (se<<4));
            print_substr(ref_txt, rb, re);
        }
#endif
//        if (re != -1 && align_seg(ptxt, sb<<4, se<<4, rb, re, dir)) { 
        if (re != -1) { 
            sb = se; 
            rb = re;
            dist = 0;
            l_matched_seg += abs(se - sb);
            ++n_matched_seg;
        } 

        if (dist > N_SEGMENT) {  
            if (8*l_matched_seg < abs(se - ap_s)) 
                mismatch = true;
            else
                sb = se;
        }
    } while (!mismatch && se > 0 && se < seq_len_in_word);

    if (!mismatch && n_matched_seg > STRONG) {
//        if (dir == 1) { 
//            // prefix of seq match suffix of reference
//            if ((seq_len-(sb<<4)) > (ref_len-rb) && n_matched_seg > STRONG) {  
//                align_seg(ptxt, sb<<4, ref_len-rb+(sb<<4), rb, ref_len, dir);
////                for (int i = ref_len-rb+sb; i < seq_len; ++i)
////                    votes.push_back(Vote(ptxt[i]));
////                print_substr(ptxt, ref_len-rb+(sb<<4), seq_len);
//            } else {
//                align_seg(ptxt, sb<<4, seq_len, rb, rb+seq_len-(sb<<4), dir);
//            }
//        } else { 
//            // suffix of seq match prefix of reference
//            if ((sb<<4) > rb && n_matched_seg > STRONG) {
//                align_seg(ptxt, (sb<<4)-rb, (sb<<4)+16, 0, rb+16, dir);
////                for (int i = 0; i < seq_len-rb; ++i)
////                    votes.push_front(Vote(ptxt[i]));
////                print_substr(ptxt, 0, (sb<<4)-rb);
//                ref_beg += seq_len-rb;
//            } else {
//                align_seg(ptxt, 0, (sb<<4)+16, rb-(sb<<4), rb+16, dir);
//            }
//        }

#ifdef DBG
        print_substr(ptxt, ap_s, (sb<<4));
        print_substr(ref_txt, ap_r, rb);
        LOG("nseg: %d\n", n_matched_seg);
#endif
    }  

    free(ptxt);
    return !mismatch && n_matched_seg > STRONG;
}		/* -----  end of function align  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  try_align
 *  Description:  pos is index into seq in units of 16
 * ===========================================================================
 */
    inline bool
try_align ( t_bseq *seq, size_t pos, int dir)
{
    unsigned *pseg = (unsigned*)(seq + sizeof(unsigned));
    sm_it sit = seedmap.find(pseg[pos] & seed);
    if (sit == seedmap.end()) return false;

#ifdef DBG
    ++_ntrials;
#endif
    ls_it it = sit->second.begin();
    ls_it end = sit->second.end();
    for (; it != end; ++it) {
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
open_binary ( const char *fname, std::list<int> &indices )
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

    if (argc < 3) 
        handle_error("usage: spaced_seed bin seed\n");

    // read binary sequence file and build index for all DNA sequences
    i_max_len = open_binary(argv[1], indices);
    // set the longest DNA sequence as initial reference
    
    i_max_len = *indices.begin();
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
        seq_len = (*((unsigned*)seq) >> 4); // seq length in 32-word
        if (seq_len == 0) {
            LOG("segment too short! length: %d\n", get_seq_len(seq));
            continue;
        }
        // number of trial 
        for (size_t j = 0; j < MIN(seq_len, N_TRIAL); ++j) {
            if (try_align(seq, j, 1) || try_align(seq, seq_len-j-1, -1)) {
#ifdef DBG
                LOG("found %d\n", count+1);
                ++_nmatches;
#endif
                break;
            }
        }
//        it = found ? indices.erase(it) : ++it;
        if (!(++count & 0xFFF)) LOG("%d sequences processed\n", count);
    }

#ifdef DBG
    LOG("#trials: %d\n", _ntrials);
    LOG("#matches: %d\n", _nmatches);
#endif

    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
