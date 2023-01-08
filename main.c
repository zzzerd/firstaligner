#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
//#include "kseq.h"
//#include "struct_ref_read.h"
//#include "idx.h"
#include "aln.h"

// STEP 1: declare the type of file handler and the read() function
unsigned char table_zm[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
u_int64_t bidui() {
    u_int64_t count = 0;
    FILE* f_sam = fopen("result.sam", "r");
    FILE* f_sam_1 = fopen("single_dat.sam", "r");
    while ((fgetc(f_sam_1)) == '@') //跳过头部
    {
        /* code */
        while (fgetc(f_sam_1) == '\n') {
            /* code */
        }
    }
    while (fgetc(f_sam_1) != EOF || fgetc(f_sam) != EOF) {
        int count_sam = 0;
        while (count_sam == 3) {

            while (fgetc(f_sam_1) == '\t') {
                /* code */
            }
            while (fgetc(f_sam) == '\t') {
                /* code */
                count_sam++;
            }
        }
        int i = fgetc(f_sam_1);
        int j = fgetc(f_sam);
        if (i == j) {
            count++;
        }
        while (fgetc(f_sam_1) == '\n') {
            /* code */
        }
        while (fgetc(f_sam) == '\n') {
            /* code */
        }
    }
    printf("%d", count);
    return count;
    fclose(f_sam);
    fclose(f_sam_1);
}
KSEQ_INIT(gzFile, gzread)
int main(int argc, char* argv[]) {
    gzFile fp;
    kseq_t* seq;
    uint64_t bases = 0;
    //char* seqs;
    int l;
    //int c;
    uint64_t k_mers;
    Yuanzu* yuanzu = NULL;
    unsigned char** p_yuanzu = NULL;
    W_yuanzu* w_yuanzu = NULL;
    //char getzhi;
    u_int64_t count = 0;
    /*     long k_mers_1 = 0;//使用哈希函数后的剩余量
        int flag ; //k-mer哈希后是否留下
        unsigned char test[32];
        long i=0;
        long chushi = 0;
        long chushi2 = k_mers_1-1;
        unsigned char test_22[22];
        long p_k_mer = 0;//p数组k-mer数量 */

    if (argc == 1) {
        fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
        return 1;
    }
    /*    while ((c = getopt(argc, argv, "a:b:")) >= 0) {
           switch (c) {
               case 'a': ref; break;
               case 'b': aligner; break;

           }
       } */
    fp = gzopen(argv[1], "r"); // STEP 2: open the file handler
    seq = kseq_init(fp);       // STEP 3: initialize seq
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        bases += strlen(seq->seq.s);
        //seqs = seq->seq.s;
        k_mers = bases - 32 + 1;
        yuanzu = (Yuanzu*)malloc(k_mers * sizeof(Yuanzu));
        for (uint64_t bases_1 = 0; bases_1 < bases; bases_1++) {
            unsigned char table_1[4] = { 0, 1, 2, 3 };
            if (table_zm[seq->seq.s[bases_1]] == 4) {
                /* code */

                srand((unsigned)time(0));
                seq->seq.s[bases_1] == table_1[rand() % 4];
            }
            else {
                seq->seq.s[bases_1] = table_zm[seq->seq.s[bases_1]];
            }
        }

        ref(p_yuanzu, w_yuanzu, yuanzu, seq->seq.s, k_mers);
        aligner(argc, argv, p_yuanzu, w_yuanzu, seq->seq.s);
        count = bidui();
    }
    kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(fp);       // STEP 6: close the file handler

    return 0;
}