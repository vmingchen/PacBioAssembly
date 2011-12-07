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
#define SEQ_THRESHOLD 500
#define N_SEGMENT 100
#define N_TRIAL 50
#define MAX_PAT_LEN N_SEQ_WORD
#define handle_error(msg) do { perror(msg); exit(EXIT_FAILURE); } while (0)

#ifdef DBG
size_t _ntrials = 0;
size_t _nfound = 0;
#endif

const char *usage_str = "usage: %s [options] bin seedfile\n"
    "options: [-f:r:d:m:t:lh]\n"
    "   -h          Get help and usage.\n"
    "   -f file     Use the string from file as starting reference.\n"
    "               Without this option the problem will select a random\n"
    "               segments as starting reference instead.\n"
    "   -r ratio    Ratio of difference (0.3 by default) allowed.\n"
    "   -d dumpfile Dump matched segments.\n"
    "   -m nround   Maximum number of round of iteration.\n"
    "   -t ntrials  Number of seeding trial for each segment.\n"
    "   -l          Lock reference during iteration.\n";

class seq_index {
public:
    seq_index(int i, int o) : id(i), offset(o) {};
    int id;             // the id-th sequence
    int offset;         // offset into binary file 
};

// spaced seed
unsigned seed = 0;

// buf for binary DNA sequence
t_bseq *buf = NULL;
// indices for binary DNA sequence
std::list<seq_index> indices;  

t_aligner *paligner = NULL;
ref_seq *pref = NULL;

// information of active segment
t_bseq *seg_bin; 
char seg_txt[MAX_SEQ_LEN];
int seg_len;
int seg_id = -1;

// seedmap for reference sequence
hash_table seedmap(1<<20);
std::vector<unsigned> seeds; 

// max number of iteration round
int max_round = INT_MAX;
int max_trial = 32;

FILE *fpdump = NULL;
FILE *fpref = NULL;

typedef std::list<int>::iterator list_it;
typedef std::list<seq_index>::iterator index_it;

inline unsigned get_seq_len(const t_bseq *x) { return *((unsigned *)x); }

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  set_active_seg
 *  Description:  set current segment of interest
 * ===========================================================================
 */
    void
set_active_seg ( seq_index &idx )
{
    if (idx.id != seg_id) {
        seg_id = idx.id;
        seg_bin = buf + idx.offset;
        seg_len = get_seq_len(seg_bin);
        dna_seq::bin2text(seg_bin, seg_txt, seg_len+1); 
    }
}		/* -----  end of function set_active_seg  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  dump_seq
 *  Description:  
 * ===========================================================================
 */
    void
dump_seq ( FILE *fp, seq_accessor *pac, int length )
{
    assert(pac->length() >= length);
    for (int i = 0; i < length; ++i) 
        fputc(pac->next(), fp);
    fputc('\n', fp);
}		/* -----  end of function dump_seq  ----- */

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
    for (index_it it = indices.begin(); fin >> qualtiy; ++it) {
        int seq_len = get_seq_len(buf + it->offset);
        if (seq_len < 2000) continue;
        if (qualtiy < best_val 
                || (qualtiy == best_val && seq_len > best_len)) {
            best_val = qualtiy;
            best_idx = it->offset;
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
init ( FILE *fp, const char *seed_file, double ratio, bool l )
{
    char tmp[MAX_SEQ_LEN];

    // init random function with a seed
    srand( (unsigned)time(0) );

    // set reference
    if (fp) {     // from file
        if (fgets(tmp, MAX_SEQ_LEN, fp) == NULL)
            handle_error("failed to open ref_file");
        pref = new ref_seq(tmp, strlen(tmp), l);
        fclose(fp);
    } else {                        // from random-selected segment
        index_it it = indices.begin();
        advance(it, rand() % indices.size());
        pref = new ref_seq(buf + it->offset, l);
        LOG("%d selected as the initial reference.\n", it->id);
    }
    assert(pref != NULL);
    LOG("ref_len: %d\n", pref->length());

    // instantiate aligner
    paligner = new t_aligner(ratio);
    assert(paligner != NULL);

    // parse spaced seed
    fp = fopen(seed_file, "r");
    if (fp == NULL) 
        handle_error("failed to open seedfile");

    char ptn_str[1024];
    while (fgets(ptn_str, 1024, fp) != NULL) {
        ptn_str[strlen(ptn_str)-1] = '\0';  // get rid of new-line 
        seeds.push_back(parse_pattern(ptn_str));
        LOG("seed %s: %08x\n", ptn_str, seeds.back());
    } 
    fclose(fp);
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
 *         Name:  try_align
 *  Description:  try align segment to reference from postion at pos in
 *  direction of dir. Return true if aligned. 
 * ===========================================================================
 */
    inline bool
try_align ( seq_index &idx, size_t pos, int dir)
{
    t_bseq *seq = buf + idx.offset;
    sm_it sit = seedmap.find(dna_seq::seed_at(seq, pos) & seed);
    if (sit == seedmap.end()) return false;

#ifdef DBG
    ++_ntrials;
#endif

    set_active_seg(idx);

    bool forward = dir == 1;
    int s_offset = forward ? pos : pos+16-1;
    int s_len = forward ? seg_len - s_offset : s_offset + 1;
    seq_accessor ac_seg(seg_txt+s_offset, forward, s_len);

    // too short to justify overlap
    if (s_len < OVERLAP_MIN) return false;       

    list_it it = sit->second.begin();
    list_it end = sit->second.end();
    for (; it != end; ++it) {
        int r_offset = forward ? (*it) : (*it)+16-1;
        if (pref->try_align(paligner, r_offset, &ac_seg)) { 
            if (fpdump) { 
                seq_accessor ac_ref = pref->get_accessor(r_offset, forward);
                dump_seq(fpdump, &ac_ref, paligner->matlen_a);
                ac_seg.reset(0);
                dump_seq(fpdump, &ac_seg, paligner->matlen_b); 
                fflush(fpdump);
            }
            return true;
        }
    }

    return false;
}		/* -----  end of function try_align  ----- */

/* 
 * ===  FUNCTION  ============================================================
 *         Name:  open_binary
 *  Description:  open binary sequence file, mmap to buf, build index of all
 *  segments into indices, and return the index of the longest segment.
 *  Segments too long or too short is ignored. 
 * ===========================================================================
 */
    size_t
open_binary ( const char *fname, std::list<seq_index> &indices )
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

    if (close(fd) != 0)
        handle_error("close");

    int i = 0;
    for (size_t offset = 0; offset < len; ) {
        size_t seq_len = *((unsigned*)(buf + offset));
        // make sure segments in indices are not too short or too long
        if (seq_len > SEQ_THRESHOLD && seq_len < MAX_READ_LEN) { 
            indices.push_back(seq_index(i++, offset));
        }
        if (seq_len > max_len) {
            max_len = seq_len;
            i_max_len = offset;
        }
        offset += sizeof(unsigned) + (seq_len + 4 - 1)/4;
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
    size_t i_max_len;
    int opt;
    double ratio = MAXR;
    bool locked = false;
    char ref_file[PATH_MAX] = {0};

    if (argc < 3) {
        fprintf(stderr, usage_str, argv[0]);
        return EXIT_FAILURE;
    }

    while ((opt = getopt(argc, argv, "f:r:d:m:t:lh")) != -1) {
        switch (opt) {
            case 'h':
                fprintf(stdout, usage_str, argv[0]);
                return EXIT_SUCCESS;
            case 'f':
                if ((fpref = fopen(optarg, "r")) == NULL)
                    handle_error("failed to read ref_file"); 
                break;
            case 'd':
                if ((fpdump = fopen(optarg, "w")) == NULL) 
                    handle_error("failed to create dump file");
                break;
            case 'r':
                ratio = atof(optarg);
                break;
            case 'l':
                locked = true;
                break;
            case 'm':
                max_round = atoi(optarg);
                break;
            case 't':
                max_trial = atoi(optarg);
                break;
            default: /*  '?' */
                fprintf(stderr, usage_str, argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    // read binary sequence file and build index for all DNA sequences
    i_max_len = open_binary(argv[optind], indices);
    LOG("indices: size %d\n", indices.size());
    LOG("number of seeding trial: %d\n", max_trial);
//    LOG("i_max_len: %d\n", i_max_len);

    // set the best DNA sequence as initial reference
//    int iref = select_ref("src/quality.in");
    // pick up a random segment as the starting reference
    init(fpref, argv[optind+1], ratio, locked);

    for (int nround = 1; nround <= max_round; ++nround) { 
        LOG("--------------- round %d ---------\n", nround);
        LOG("seed: %08x\n", (seed = seeds[rand() % seeds.size()]));
        LOG("seedmap size: %d\n", pref->get_seedmap(seedmap, seed));
        LOG("reference length: %d\n", pref->length());
        int nmatches = 0;
        int count = 0;
        index_it it = indices.begin();
        while (it != indices.end()) {
            bool found = 0;
            unsigned slen = get_seq_len(buf + it->offset);
            // number of trial 
            for (size_t j = 0; j < max_trial; ++j) {
                // try both forward and backward
                if (try_align(*it, j, 1) || try_align(*it, slen-j-16, -1)) {
                    found = 1;
#ifdef DBG
                    LOG("found %d at cost %d:\tref_ml=%d,\tseg_ml=%d\n",
                            it->id, paligner->final_cost(), 
                            paligner->matlen_a, paligner->matlen_b);
#endif
                    ++nmatches;
                    break;
                }
            }
            it = found ? indices.erase(it) : ++it;
            if (!(++count & 0xFFFF)) LOG("%d sequences processed\n", count);
        }
#ifdef DBG
        LOG("#trials: %d\n", _ntrials);
        LOG("#matches: %d\n", nmatches);
#endif
        if (nmatches == 0) break;
        pref->evolve();
        // print out consensus
        seq_accessor ac_ref = pref->get_accessor(0, true); 
        dump_seq(stdout, &ac_ref, pref->length());
    }

    return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
