#include <zlib.h>  
#include <stdio.h>
#include <string.h>  
#include <time.h>
#include "klib/kseq.h"  

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)
unsigned char table[4] = {'0','1','2','3'};
typedef struct 
{
    unsigned char k_mer[32];
    int Loction;
    /* data */
}Yuanzu;
//归并排序
void printList(int arr[], int len) {
	int i;
	for (i = 0; i < len; i++) {
		printf("%d\t", arr[i]);
	}
}
void merge(int arr[], int start, int mid, int end, int length) {
	int result[length];
	int k = 0;
	int i = start;
	int j = mid + 1;
	while (i <= mid && j <= end) {
		if (arr[i] < arr[j]){
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
 
void mergeSort(int arr[], int start, int end, int length) {
	if (start >= end)
		return;
	int mid = ( start + end ) / 2;
	mergeSort(arr, start, mid, length);
	mergeSort(arr, mid + 1, end, length);
	merge(arr, start, mid, end, length);
}

int main(int argc, char *argv[])
{
    gzFile fp;
    kseq_t *seq;
    long bases = 0;
    char *seqs;
    int l;
    if (argc == 1) {
        fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
        return 1;
    }
    fp = gzopen(argv[1], "r"); // STEP 2: open the file handler
    seq = kseq_init(fp); // STEP 3: initialize seq
    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
        //printf("name: %s\n", seq->name.s);
        //if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
        //printf("seq: %s\n", seq->seq.s);
        //if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
        bases += strlen(seq->seq.s);
        seqs[bases]=seq->seq.s;
    }
   for(long bases_1=0; bases_1 < bases; bases_1++){
       if (seqs[bases_1] == 'A')
          seqs[bases_1] == table[0];
        else if (seqs[bases_1] == 'C')
          seqs[bases_1] == table[1];
        else if (seqs[bases_1] == 'G')
          seqs[bases_1] == table[2];
        else if (seqs[bases_1] == 'T')
          seqs[bases_1] == table[3];
        else 
          srand((unsigned)time(0));
          seqs[bases_1] == table[rand()%4];
   }
    long k_mers = bases-32+1;
    Yuanzu yuanzu[k_mers];
    long k_mers_1 = 0;//使用哈希函数后的剩余量
    int flag ; //k-mer哈希后是否留下
    for(long i=0; i<= k_mers; i++){
        for(int j=i; j<i+32; j++){
             unsigned char test[32];
             test[j] = seqs[j];
             int number=atoi(test);
             flag = 0 ;
             if(number % 7 == 5 || number % 17 == 7 || number % 19 == 13 || number % 53 == 31 || number % 71 == 47 ){
                for(int k=i; k<i+32; k++){
                   yuanzu[i].k_mer[k] = seqs[k];
                }
                flag = 1 ;
             }  //哈希函数           
        }
        if(flag == 1){
                yuanzu[i].Loction = i;
                k_mers_1++ ;
        }
        
    }
    kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(fp); // STEP 6: close the file handler
    return 0;
}