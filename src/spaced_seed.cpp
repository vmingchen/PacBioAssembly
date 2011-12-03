/*
 * ===========================================================================
 *
 *       Filename:  spaced_seed.cpp
 *         Author:  Ming Chen, brianchenming@gmail.com
 *        Created:  11/07/2011 03:03:04 PM
 *
 *    Description:  perform spaced seed
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
#include	<fstream>

#include	"dna_seq.h"
#include	"seq_aligner.h"
#include	"common.h"
#include	"ref_seq.h"

#define STRONG 3
#define N_SEGMENT 100
#define N_TRIAL 20
#define MAX_PAT_LEN N_SEQ_WORD
#define handle_error(msg) do { perror(msg); exit(EXIT_FAILURE); } while (0)

#ifdef DBG
size_t _ntrials = 0;
size_t _nmatches = 0;
size_t _nfound = 0;
#endif

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

char tmp[MAX_SEQ_LEN];
seq_aligner *paligner = NULL;
ref_seq *pref = NULL;

// seedmap for reference sequence
hash_table seedmap(1<<20);

typedef std::list<int>::iterator ls_it;
typedef __gnu_cxx::hash_map< unsigned, std::list<int> >::iterator sm_it;

inline unsigned get_seq_len(const t_bseq *x) { return *((unsigned *)x); }


/* 
 * ===  FUNCTION  ============================================================
 *         Name:  print_seq
 *  Description:  
 * ===========================================================================
 */
    void
print_seq ( seq_accessor *pac, int length )
{
    assert(pac->length() < length);
    for (int i = 0; i < length; ++i) 
        putc(pac.next());
    putc('\n');
}		/* -----  end of function print_seq  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  select_ref
 *  Description:  
 * ===========================================================================
 */
    int
select_ref ( const char *fname )
{
    std::ifstream fin(fname);
    int qualtiy;
    int best_idx, best_val = 100, best_len;
    for (ls_it it = indices.begin(); fin >> qualtiy; ++it) {
        int seq_len = get_seq_len(buf+*it);
        if (seq_len < 2000) continue;
        if (qualtiy < best_val 
                || (qualtiy == best_val && seq_len > best_len)) {
            best_val = qualtiy;
            best_idx = *it;
            best_len = seq_len;
        }
    }
    return best_idx;
}		/* -----  end of function select_ref  ----- */

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
    size_t len = std::min(strlen(pat), (size_t)MAX_PAT_LEN);

#ifdef DBG
    if (len < MAX_PAT_LEN) 
        LOG("WARNING: pattern is shorter than %d\n", MAX_PAT_LEN);
#endif

    for (size_t i = 0; i < len; ++i) 
        dnapat[i] = (pat[i] == '1' ? 'T' : 'A');
    return dna_seq::encode(dnapat);
}		/* -----  end of function parse_pattern  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  init
 *  Description:  init pref, paligner, seed, seedmap
 * ===========================================================================
 */
    void
init ( t_bseq *bseq, const char *ptn_str )
{
    // set reference
    pref = new ref_seq(bseq);
    assert(pref != NULL);
    LOG("ref_len: %d\n", pref->length());

    // instantiate aligner
    paligner = new seq_aligner(0.25);
    assert(paligner != NULL);

    // parse spaced seed
    seed = parse_pattern(ptn_str);
    LOG("seed: %x\n", seed);

    // build seedmap for reference sequences (the longest sequence)
    LOG("seedmap (size: %d) done!\n", pref->get_seedmap(seedmap, seed));
}		/* -----  end of function init  ----- */

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
    dna_seq::decode(sv & seed, seg);
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
 *         Name:  try_align
 *  Description:  pos is index into seq 
 * ===========================================================================
 */
    inline bool
try_align ( t_bseq *seq, size_t pos, int dir)
{
    sm_it sit = seedmap.find(dna_seq::seed_at(seq, pos) & seed);
    if (sit == seedmap.end()) return false;

#ifdef DBG
    ++_ntrials;
#endif

    int seq_len = get_seq_len(seq);
    dna_seq::bin2text(seq, tmp, seq_len+1); 
    bool forward = dir == 1;
    int s_offset = forward ? pos : pos+16-1;
    int s_len = forward ? seq_len - s_offset : s_offset + 1;
    seq_accessor pseq_acsr(tmp+s_offset, forward, s_len);

    if (s_len < 200) return false;

    ls_it it = sit->second.begin();
    ls_it end = sit->second.end();
    for (; it != end; ++it) {
        int r_offset = forward ? (*it) : (*it)+16-1;
        if (pref->try_align(paligner, r_offset, &pseq_acsr)) { 
#ifdef DBG
            fprintf(stderr, "ref_ml = %d, seg_ml = %d\n", 
                    paligner->ref_ml, paligner->seg_ml);
            char fmt[100];
            if (_nfound < 20) {
                seq_accessor ac_ref = pref->get_accessor(r_offset, forward);
                print_seq(&ac_ref, paligner->seg_ml);
            }

            if (forward && _nfound < 20) {


                print_substr(pref->get_accessor(r_offset, forward).pt(0),
                        0, paligner->ref_ml);

                sprintf(fmt, "%%.%ds\n", paligner->ref_ml);
                printf(fmt, ref_txt+r_offset);

                sprintf(fmt, "%%.%ds\n", paligner->seg_ml);
                printf(fmt, tmp+s_offset);

                fflush(stdout);
                ++_nfound;
            }
            if (!forward && _nfound < 20) {
                sprintf(fmt, "%%.%ds\n", paligner->ref_ml);
                printf(fmt, ref_txt+r_offset-paligner->ref_ml+1);
                sprintf(fmt, "%%.%ds\n", paligner->seg_ml);
                printf(fmt, tmp+s_offset-paligner->seg_ml+1);
                ++_nfound;
            }
#endif
            return true;
        }
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

    // set the best DNA sequence as initial reference
//    int iref = select_ref("src/quality.in");
    ls_it it = indices.begin();
    srand( (unsigned)time(0) );
    for (int i = rand() % indices.size(); i > 0; --i)
        ++it;
    while (get_seq_len(buf + *it) > 5000 || get_seq_len(buf + *it) < 2000)
        ++it;
    init(buf + *it, argv[2]);

    LOG("i_max_len: %d\n", i_max_len);
    LOG("indices: size %d\n", indices.size());
    // find repeat
    int count = 0;
    for (it = indices.begin(); it != indices.end(); ++it) {
//        if (*it == i_max_len) continue;
        seq = buf + *it;
        seq_len = get_seq_len(seq);
        if (seq_len < N_TRIAL+16-1) {
//            LOG("segment too short! length: %d\n", get_seq_len(seq));
            continue;
        }

        // number of trial 
        for (size_t j = 0; j < N_TRIAL; ++j) {
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
