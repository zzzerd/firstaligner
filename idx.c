#include <zlib.h>  
#include <stdio.h>
#include <string.h>  
#include <time.h>
#include "klib/kseq.h"  
#include "struct_ref_read.h"
KSEQ_INIT(gzFile, gzread)
//归并排序
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
int ref(P_yuanzu *p_yuanzu, W_yuanzu *w_yuanzu, Yuanzu *yuanzu,char *seqs, long k_mers)
{
    long k_mers_1 = 0;//使用哈希函数后的剩余量
    int flag ; //k-mer哈希后是否留下
    unsigned char test[32];
    long i=0;
    long chushi = 0;
    long chushi2 = k_mers_1-1;
    unsigned char test_22[22];
    long p_k_mer = 0;//p数组k-mer数量
    
    
    for( i=0; i< k_mers; i++)
    {
        int j_1 =0 ;
        for(int j=i; j<i+32; j++)
        {                    
                  test[j_1] = seqs[j];
                  j_1++;
             
        }
             int number=atoi(test);
             flag = 0 ;
             if(number % 7 == 5 || number % 17 == 7 || number % 19 == 13 || number % 53 == 31 || number % 71 == 47 )
             {
                int k_1 =0;
                for(int k=i; k<i+32; k++)
                {
                    yuanzu[k_mers_1].k_mer[k_1] = seqs[k];
                    k_1++;
                }
                   
                flag = 1 ;
             }  //哈希函数           
              if(flag == 1)
              {
                yuanzu[k_mers_1].Loction = i;
                k_mers_1++ ;
              }
        
    }
  
    mergeSort(yuanzu, chushi, chushi2, k_mers_1);//排序基因组k-mer
    long w_loction=0;
    p_yuanzu = (P_yuanzu*)malloc(k_mers_1*sizeof(P_yuanzu));
    w_yuanzu = (W_yuanzu*)malloc(k_mers_1*sizeof(W_yuanzu));
/*     P_yuanzu p_yuanzu[k_mers_1];
    W_yuanzu w_yuanzu[k_mers_1]; */
    while (w_loction<k_mers_1)
    {
      /* code */
      int w_dingwei = 0 ; //用于w数组的k-mer
      strncpy(p_yuanzu[p_k_mer].k_mer22, yuanzu[w_loction].k_mer, 22);
      p_yuanzu[p_k_mer].start=w_loction;
      for (int i_dingwei = 22; i_dingwei < 32; i_dingwei++)
      {
        /* code */
        w_yuanzu[w_loction].k_mer10[w_dingwei] = yuanzu[w_loction].k_mer[i_dingwei];
        w_dingwei++;
      }
      w_yuanzu[w_loction].Loction1 = yuanzu->Loction;
      w_loction++;
      for (; strcmp(p_yuanzu[p_k_mer].k_mer22, yuanzu[w_loction].k_mer)== 0; w_loction++)
      {
        /* code */
        int w_dingwei2 = 0 ;
        for (int i_dingwei = 22; i_dingwei < 32; i_dingwei++)
        {
        /* code */
        w_yuanzu[w_loction].k_mer10[w_dingwei2] = yuanzu[w_loction].k_mer[i_dingwei];
        w_dingwei2++;
        }
      }
      p_yuanzu[p_k_mer].end = w_loction-1;
      p_k_mer++;
    }//p数组，w数组
    return p_k_mer;
}