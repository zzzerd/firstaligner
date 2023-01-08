#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "klib/kseq.h"
#include "struct_ref_read.h"
#include <stdint.h>
#include "idx.c"
#include "klib/ksw.h"
#define MAX_LOC 30
#define GAP_CHAR '-'
#define readkmernum 79
KSEQ_INIT(gzFile, gzread)
typedef struct {
    /* data */
    int min;
    char* cig_1;
    int ref_loc;
} Alg_read;

typedef struct {
    /* data */
    int start;
    int end;

} Less_max;
typedef struct {
    /* data */
    int start;
    int end;

} Over_max;
typedef struct {
    unsigned char k_mer_read[22];
    int loction_read;
}Read_kmer;
int g(int l, int m);
void rle(char* cig);
Alg_read hanming(char* seq_read_00, char* seqs, Less_max* less_max, Over_max* over_max, int lessloction, int overloction, W_yuanzu* w_yuanzu, u_int64_t bases, int* n_cigar_,
    uint32_t** cigar_);
void aligner(int argc, char* argv[], unsigned char** p_yuanzu, W_yuanzu* w_yuanzu, char* seqs);