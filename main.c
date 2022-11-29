#include <zlib.h>  
#include <stdio.h>
#include <string.h>  
#include <time.h>
#include "klib/kseq.h"  
#include "struct_ref_read.h"
#include "idx.c"
// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)


/* //归并排序
void merge(Yuanzu arr[], long start, long mid, long end, long length) {
	Yuanzu result[length];
	long k = 0;
	long i = start;
	long j = mid + 1;
	while (i <= mid && j <= end) {
		if (strcmp(arr[i].k_mer, arr[j].k_mer)< 0)
    {
			result[k++] = arr[i++];
    }
		else{
			result[k++] = arr[j++];
        }
	}
	if (i == mid + 1) {
		while(j <= end)
			result[k++] = arr[j++];
	}
	if (j == end + 1) {
		while (i <= mid)
			result[k++] = arr[i++];
	}
	for (j = 0, i = start ; j < k; i++, j++) {
		arr[i] = result[j];
	}
}
 
void mergeSort(Yuanzu arr[], long start, long end, long length) {
	if (start >= end)
		return;
	long mid = ( start + end ) / 2;
	mergeSort(arr, start, mid, mid+1);
	mergeSort(arr, mid + 1, end, mid+1);
	merge(arr, start, mid, end, length);
}
 */
int main(int argc, char *argv[])
{
    gzFile fp;
    kseq_t *seq;
    long bases = 0;
    char *seqs;
    int l;
    long k_mers = bases-32+1;
    Yuanzu *yuanzu;
    yuanzu = (Yuanzu*)malloc(k_mers*sizeof(Yuanzu));
    P_yuanzu *p_yuanzu;
    W_yuanzu *w_yuanzu;
    long p_k_mer = 0;
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
    fp = gzopen(argv[1], "r"); // STEP 2: open the file handler
    seq = kseq_init(fp); // STEP 3: initialize seq
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        bases += strlen(seq->seq.s);
        seqs=seq->seq.s;
    }
   for(long bases_1=0; bases_1 < bases; bases_1++)
   {

          if (table[bases_1]==4)
          {
            /* code */
             unsigned char table_1[4]= {0, 1, 2 , 3};
             srand((unsigned)time(0));
             seqs[bases_1] == table_1[rand()%4];
          }          
          else 
          {
               seqs[bases_1] = table[bases_1];
          }
          
    }
    p_k_mer = ref(p_yuanzu, w_yuanzu, yuanzu, seqs, k_mers);
    kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(fp); // STEP 6: close the file handler
    
    return 0;
}