#include <zlib.h>  
#include <stdio.h>
#include <string.h>  
#include <time.h>
#include "klib/kseq.h"  
#define READ 1000000
// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)
unsigned char table[4] = {'0','1','2','3'};
typedef struct 
{
    unsigned char k_mer[32];
    long Loction;
    /* data */
}Yuanzu;
typedef struct 
{
   unsigned char k_mer22[22];
   long start;
   long end;
  /* data */
}P_yuanzu;
typedef struct 
{
    unsigned char k_mer10[10];
    long Loction1;
    /* data */
}W_yuanzu;
//汉明距离
int hanming(P_yuanzu canshu[], unsigned char canshu2)
{

}
//汉明距离
//NW
struct Unit {
    int W1;   // 是否往上回溯一格
    int W2;   // 是否往左上回溯一格
    int W3;   // 是否往左回溯一格
    float X;
    float Y;
    float M;
    float O;      // 得分矩阵第(i, j)这个单元的分值，即序列s(1,...,i)与序列r(1,...,j)比对的最高得分
};
typedef struct Unit *pUnit;
void strUpper(char *s) {
    while (*s != '\0') {
        if (*s >= 'a' && *s <= 'z') {
            *s -= 32;
        }
        s++;
    }
}
float max2(float a, float b) {
    return a > b ? a : b;
}

float max3(float a, float b, float c) {
    float f = a > b ? a : b;
    return f > c ? f : c;
}
// 替换矩阵：match分值为5，mismatch分值为-4
// 数组下标是两个字符的ascii码减去65之后的和
float FMatrix[] = {
    5, 0, -4, 0, 5, 0, -4, 0, -4, 0,
    0, 0, 5, 0, 0, 0, 0, 0, 0, -4,
    0, -4, 0, 0, 0, -4, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 5
};

float getFScore(char a, char b) {
    return FMatrix[a + b - 'A' - 'A'];
}

void printAlign(pUnit** a, const int i, const int j, char* s, char* r, char* saln, char* raln, int n) {
    int k;
    pUnit p = a[i][j];
    if (! (i || j)) {   // 到矩阵单元(0, 0)才算结束，这代表初始的两个字符串的每个字符都被比对到了
        for (k = n - 1; k >= 0; k--)
            printf("%c", saln[k]);
        printf("\n");
        for (k = n - 1; k >= 0; k--)
            printf("%c", raln[k]);
        printf("\n\n");
        return;
    }
    if (p->W1) {    // 向上回溯一格
        saln[n] = s[i - 1];
        raln[n] = GAP_CHAR;
        printAlign(a, i - 1, j, s, r, saln, raln, n + 1);
    }
    if (p->W2) {    // 向左上回溯一格
        saln[n] = s[i - 1];
        raln[n] = r[j - 1];
        printAlign(a, i - 1, j - 1, s, r, saln, raln, n + 1);
    }
    if (p->W3) {    // 向左回溯一格
        saln[n] = GAP_CHAR;
        raln[n] = r[j - 1];
        printAlign(a, i, j - 1, s, r, saln, raln, n + 1);
    }
}

void align(char *s, char *r) {
    int i, j;
    int m = strlen(s);
    int n = strlen(r);
    float d = -7;     // 对第一个空位的罚分
    float e = -2;     // 第二个及以后空位的罚分
    pUnit **aUnit;
    char* salign;
    char* ralign;
    float f;
    // 初始化
    if ((aUnit = (pUnit **) malloc(sizeof(pUnit*) * (m + 1))) == NULL) {
        fputs("Error: Out of space!\n", stderr);
        exit(1);
    }
    for (i = 0; i <= m; i++) {
        if ((aUnit[i] = (pUnit *) malloc(sizeof(pUnit) * (n + 1))) == NULL) {
            fputs("Error: Out of space!\n", stderr);
            exit(1);     
        }
        for (j = 0; j <= n; j++) {
            if ((aUnit[i][j] = (pUnit) malloc(sizeof(struct Unit))) == NULL) {
                fputs("Error: Out of space!\n", stderr);
                exit(1);     
            }
            aUnit[i][j]->W1 = 0;
            aUnit[i][j]->W2 = 0;
            aUnit[i][j]->W3 = 0;
        }
    }
    aUnit[0][0]->X = d;
    aUnit[0][0]->Y = d;
    aUnit[0][0]->M = 0;
    aUnit[0][0]->O = max3(aUnit[0][0]->X, aUnit[0][0]->Y, aUnit[0][0]->M);
    for (i = 1; i <= m; i++) {
        aUnit[i][0]->X = d + (i - 1) * e;
        aUnit[i][0]->Y = 2 * d + (i - 1) * e;
        aUnit[i][0]->M = d + (i - 1) * e;
        aUnit[i][0]->O = max3(aUnit[i][0]->X, aUnit[i][0]->Y, aUnit[i][0]->M);
        aUnit[i][0]->W1 = 1;
    }
    for (j = 1; j <= n; j++) {
        aUnit[0][j]->X = 2 * d + (j - 1) * e;
        aUnit[0][j]->Y = d + (j - 1) * e;
        aUnit[0][j]->M = d + (j - 1) * e;
        aUnit[0][j]->O = max3(aUnit[0][j]->X, aUnit[0][j]->Y, aUnit[0][j]->M);
        aUnit[0][j]->W3 = 1;
    }
    // 将字符串都变成大写
    strUpper(s);
    strUpper(r);
    // 动态规划算法计算得分矩阵每个单元的分值
    for (i = 1; i <= m; i++) {
        for (j = 1; j <= n; j++) {
            aUnit[i][j]->X = max2(aUnit[i - 1][j]->X + e, aUnit[i - 1][j]->M + d);
            aUnit[i][j]->Y = max2(aUnit[i][j - 1]->Y + e, aUnit[i][j - 1]->M + d);
            f = getFScore(s[i - 1], r[j - 1]);
            aUnit[i][j]->M = max3(aUnit[i - 1][j - 1]->X + f, aUnit[i - 1][j - 1]->Y + f, aUnit[i - 1][j - 1]->M + f);
            aUnit[i][j]->O = max3(aUnit[i][j]->X, aUnit[i][j]->Y, aUnit[i][j]->M);
            if (aUnit[i][j]->O == aUnit[i][j]->X) aUnit[i][j]->W1 = 1;
            if (aUnit[i][j]->O == aUnit[i][j]->M) aUnit[i][j]->W2 = 1;
            if (aUnit[i][j]->O == aUnit[i][j]->Y) aUnit[i][j]->W3 = 1;
        }
    }
/*
    // 打印得分矩阵
    for (i = 0; i <= m; i++) {
        for (j = 0; j <= n; j++)
            printf("%f ", aUnit[i][j]->O);
        printf("\n");
    }
*/
    printf("max score: %f\n", aUnit[m][n]->O);
    // 打印最优比对结果，如果有多个，全部打印
    // 递归法
    if ((salign = (char*) malloc(sizeof(char) * (m + n + 1))) == NULL) {
        fputs("Error: Out of space!\n", stderr);
        exit(1);
    }
    if ((ralign = (char*) malloc(sizeof(char) * (m + n + 1))) == NULL) {
        fputs("Error: Out of space!\n", stderr);
        exit(1);
    }
    printAlign(aUnit, m, n, s, r, salign, ralign, 0);
    // 释放内存
    free(salign);
    free(ralign);
    for (i = 0; i <= m; i++) {
        for (j = 0; j <= n; j++) {
            free(aUnit[i][j]);
        }
        free(aUnit[i]);
    }
    free(aUnit);
}
//NW
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

int main(int argc, char *argv[])
{
    gzFile fp;
    kseq_t *seq;
    long bases = 0;
    char *seqs;
    int l;
    long k_mers = bases-32+1;
    Yuanzu yuanzu[k_mers];
    long k_mers_1 = 0;//使用哈希函数后的剩余量
    int flag ; //k-mer哈希后是否留下
    unsigned char test[32];
    long i=0;
    long chushi = 0;
    long chushi2 = k_mers_1-1;
    unsigned char test_22[22];
    long p_k_mer = 0;//p数组k-mer数量
    
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
        seqs=seq->seq.s;
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
    P_yuanzu p_yuanzu[k_mers_1];
    W_yuanzu w_yuanzu[k_mers_1];
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
    
    kseq_destroy(seq); // STEP 5: destroy seq
    gzclose(fp); // STEP 6: close the file handler
    //read部分
    FILE *fp_sam;
    int i, a[6],b[6]; //i为数组元素下标
    fp_sam = fopen("a.txt", "w");
    if (fp_sam == NULL)
    {
        printf("Can't open the file!\n");
    }
    /* for (i = 0; i < 6; i++)
    {
        fprintf(fp, "%d ", a[i]);
    } */
    
    gzFile fp_raed;
    kseq_t *seq_read;
    int l_read;
    if (argc == 1) 
    {
    fprintf(stderr, "Usage: %s <in.fasta>\n", argv[0]);
    return 1;
    }
    
    fp_raed = gzopen(argv[2], "r");
    seq_read = kseq_init(fp_raed);
    
    while ((l_read = kseq_read(seq_read)) >= 0) 
    {
/*        printf("name: %s\n", seq->name.s);
       if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
       printf("seq: %s\n", seq->seq.s);
       if (seq->qual.l) printf("qual: %s\n", seq->qual.s); */
      char *read_name = seq_read->name;
      fprintf(fp_raed, "%s"+/t, read_name);
    }
    fclose(fp_sam);
    kseq_destroy(seq_read);
    gzclose(fp_raed);
    return 0;
}