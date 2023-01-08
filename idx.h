#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "klib/kseq.h"
#include "struct_ref_read.h"
#include <stdint.h>
#define lengthkmer 32
u_int64_t kmernumber(unsigned char test[]);
void merge(Yuanzu arr[], uint64_t start, uint64_t mid, uint64_t end, uint64_t length);
void mergeSort(Yuanzu arr[], uint64_t start, uint64_t end, uint64_t length);
void ref(unsigned char** p_yuanzu, W_yuanzu* w_yuanzu, Yuanzu* yuanzu, char* seqs, uint64_t k_mers);