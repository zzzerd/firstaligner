//read部分
#include <zlib.h>  
#include <stdio.h>
#include <string.h>  
#include <time.h>
#include "klib/kseq.h"  
#include "struct_ref_read.h"
#define MAX_LOC 30
#define GAP_CHAR '-'
KSEQ_INIT(gzFile, gzread)
typedef struct 
{
   /* data */
   int min;
   char *cig_1;
   int ref_loc;
}Alg_read;

typedef struct 
{
  /* data */
  int start;
  int end;

}Less_max;
typedef struct 
{
  /* data */
  int start;
  int end;

}Over_max;
struct Unit {
    int W1;   // 是否往上回溯一格
    int W2;   // 是否往左上回溯一格
    int W3;   // 是否往左回溯一格
    int X;
    int Y;
    int M;
    int O;      // 得分矩阵第(i, j)这个单元的分值，即序列s(1,...,i)与序列r(1,...,j)比对的最高得分
};
typedef struct Unit *pUnit;

int g(int l , int m)
{
   if(l > (m+2))
    return 1 ;
    else
    return 0 ;
}
void rle(char *cig)//游程编码
{
 int count = 0;
 char samechar;
 int length = sizeof(cig);
 char *cig_rle;
 cig_rle = (char*)malloc(length*sizeof(char));
 char temp = cig[0];
 int loction = 0;
 int rle_lo = 0; 
 samechar = temp;
 while (loction < 100)
 {
   /* code */

   if (temp == samechar)
   {
      count++;
      loction++;
      temp = cig[loction];
      continue;
   }
   else
   {
      cig_rle[rle_lo] = count;
      rle_lo++;
      cig_rle[rle_lo] = samechar;
      rle_lo++;
      samechar = temp;
      count = 0;
   }
 }
 memset(cig,'\0',sizeof(cig));
 strcpy(cig,cig_rle);
}
int max2(int a, int b) {
    return a > b ? a : b;
}
int max3(int a, int b, int c) {
    int f = a > b ? a : b;
    return f > c ? f : c;
}
//match分值为5，mismatch分值为-4
int getFScore(char a, char b) {
    if(a == b)
    return 5 ;
    else
    return -4;
}
void printAlign(pUnit** a, const int i, const int j, char* s, char* r, char* saln, char* raln, int n, char *cig_n) {
   // int k;
    pUnit p = a[i][j];
/*     if (! (i || j)) {   // 到矩阵单元(0, 0)才算结束，这代表初始的两个字符串的每个字符都被比对到了
        for (k = n - 1; k >= 0; k--)
            printf("%c", saln[k]);
        printf("\n");
        for (k = n - 1; k >= 0; k--)
            printf("%c", raln[k]);
        printf("\n\n");
        return;
    } */
    if (p->W1) {    // 向上回溯一格
        saln[n] = s[i - 1];
        raln[n] = GAP_CHAR;
        cig_n[n] = 'I';
        printAlign(a, i - 1, j, s, r, saln, raln, n + 1, cig_n);
    }
    if (p->W2) {    // 向左上回溯一格
        saln[n] = s[i - 1];
        raln[n] = r[j - 1];
        cig_n[n] = 'M';
        printAlign(a, i - 1, j - 1, s, r, saln, raln, n + 1, cig_n);
    }
    if (p->W3) {    // 向左回溯一格
        saln[n] = GAP_CHAR;
        raln[n] = r[j - 1];
        cig_n[n] = 'D';
        printAlign(a, i, j - 1, s, r, saln, raln, n + 1, cig_n);
    }
}

void align(char *s, char *r, int *max, char *cig_n) {
    int i, j;
    int m = strlen(s);
    int n = strlen(r);
    int d = -7;     // 对第一个空位的罚分
    int e = -2;     // 第二个及以后空位的罚分
    pUnit **aUnit;
    char* salign;
    char* ralign;
    int f;
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
    //printf("max score: %f\n", aUnit[m][n]->O);
    // 打印最优比对结果，如果有多个，全部打印
    // 递归法
    max[0] = aUnit[m][n]->O;
    if ((salign = (char*) malloc(sizeof(char) * (m + n + 1))) == NULL) {
        fputs("Error: Out of space!\n", stderr);
        exit(1);
    }
    if ((ralign = (char*) malloc(sizeof(char) * (m + n + 1))) == NULL) {
        fputs("Error: Out of space!\n", stderr);
        exit(1);
    }
    cig_n = (char*) malloc(sizeof(char) * (m + n + 1));
    printAlign(aUnit, m, n, s, r, salign, ralign, 0, cig_n);

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

Alg_read NW (char *seq_read_00, char *seqs, Less_max *less_max, Over_max *over_max, int lessloction, int overloction, W_yuanzu *w_yuanzu, long bases)
{
  
  char *r;
  char *cig_n;
  Alg_read alg;
   alg.min = 0;
   int l = 0;
   alg.ref_loc=0;
   int flag_ov = 0 ;
   int max_1 = 0;
   int max[1];
   for(int i =0; i<=lessloction; i++)
   {
    int start = less_max[i].start;
    int end = less_max[i].end;
     for( ; start <= end ;start++)
     {
          int alg_loc =w_yuanzu[start].Loction1;
          int length = bases - alg_loc;
          r =(char *)malloc(length*sizeof(char));
          long j = 0;
          for(long i = alg_loc; i<bases; i++)
          {
               r[j] = seqs[alg_loc];
               j++;
          }
          align(r, seq_read_00, max, cig_n);
          if(max_1< max[0])
          {
            max_1 = max[0];
            rle(cig_n);
            alg.ref_loc = alg_loc;
            alg.cig_1 = cig_n;
          }
     }

   }
   for(int i =0; i<=overloction; i++)
   {
     int start = over_max[i].start;
     int end = over_max[i].end;
     for( ; start <= end ;start++)
     {
          int alg_loc =w_yuanzu[start].Loction1;
          int length = bases - alg_loc;
          r =(char *)malloc(length*sizeof(char));
          long j = 0;
          for(long i = alg_loc; i<bases; i++)
          {
               r[j] = seqs[alg_loc];
               j++;
          }
          align(r, seq_read_00, max, cig_n);
          if(max_1< max[0])
          {
            max_1 = max[0];
            rle(cig_n);
            alg.ref_loc = alg_loc;
            alg.cig_1 = cig_n;
          }
     }
   }
   return alg;
}
Alg_read hanming( char *seq_read_00, char *seqs, Less_max *less_max, Over_max *over_max, int lessloction, int overloction, W_yuanzu *w_yuanzu, long bases)
{
   Alg_read alg;
   alg.min = 0;
   int l = 0;
   alg.ref_loc=0;
   int flag_ov = 0 ;
   long base = bases;
   alg.cig_1 = (char *)malloc(100*sizeof(char));
   for(int i =0; i<=lessloction; i++)
   {
    int start = less_max[i].start;
    int end = less_max[i].end;
    char cig[100];
    int mismatch =0;
    
    for( ; start <= end ;start++)
    {
       int alg_loc =w_yuanzu[start].Loction1;
       for (int j = 0; j < 100; j++)
       {
         /* code */
         if (seq_read_00[j]==seqs[alg_loc])
         {
            /* code */
            cig[j] = '=';
         }
         else
         {
            cig[j] = 'X';
            mismatch++;
         }
         alg_loc++;
       }
       if (alg.min > mismatch)
       {
         /* code */
         alg.min = mismatch;
         rle(cig);
         memcpy(alg.cig_1,cig,sizeof(alg.cig_1));
         alg.ref_loc = w_yuanzu[start].Loction1;
       }
       
    }
    l++;
    if(g(l,alg.min))
      {
      flag_ov = 1;
      break;
      }
   }
   while (flag_ov == 0)
   {
      /* code */
      for(int i =0; i<=overloction; i++)
      {
          int start = over_max[i].start;
          int end = over_max[i].end;
          char cig[100];
          int mismatch =0;
         for( ; start <= end ;start++)
        {
             int alg_loc =w_yuanzu[start].Loction1;
             for (int j = 0; j < 100; j++)
         {
            /* code */
            if (seq_read_00[j]==seqs[alg_loc])
          {
            /* code */
            cig[j]= '=';
           }
           else
          {
            cig[j]= 'X';
            mismatch++;
         }
           alg_loc++;
         }
       if (alg.min > mismatch)
       {
         /* code */
         alg.min = mismatch;
         rle(cig);
         
         alg.ref_loc = w_yuanzu[start].Loction1;
       }
       
       }
       l++;
       if(g(l,alg.min))
      {
      break;
      }
     }
     flag_ov = 1;
   }
   
   if (alg.min>50)
   {
      /* code */
      //NW函数
    return  NW(seq_read_00, seqs, less_max, over_max, lessloction, overloction, w_yuanzu, base);
   }
   else
    return alg;
}
void aligner (int argc, char *argv[], long p_k_mer, P_yuanzu *p_yuanzu, W_yuanzu *w_yuanzu, char *seqs)
{
    FILE *fp_sam;
    long p_k_mer1 = p_k_mer;
    /* int i, a[6],b[6]; //i为数组元素下标 */
    fp_sam = fopen("result.sam", "w");
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
    int length;
    if (argc == 1) 
    {
    fprintf(stderr, "Usage: %s <in.fasta>\n", argv[0]);
    
    }
    
    fp_raed = gzopen(argv[2], "r");
    seq_read = kseq_init(fp_raed);
    
    while ((l_read = kseq_read(seq_read)) >= 0) 
    {
/*        printf("name: %s\n", seq->name.s);
       if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
       printf("seq: %s\n", seq->seq.s);
       if (seq->qual.l) printf("qual: %s\n", seq->qual.s); */
      fprintf(fp_sam, "%s\t", seq_read->name.s);
      char *seq_read_seq = seq_read->seq.s;
      length = strlen(seq_read->seq.s);
      char seq_read_00 [100];
      Less_max *less_max;
      Over_max *over_max;
      int lessloction = 0 ;
      int overloction = 0 ;
      char *cigar;
      Alg_read alg;
      for (int seq_swich = 0; seq_swich < 100; seq_swich++)
      {
        /* code */
        if (table[seq_swich]==4)
          {
            /* code */
             unsigned char table_1[4]= {0, 1, 2 , 3};
             srand((unsigned)time(0));
             seq_read_seq[seq_swich] == table_1[rand()%4];
          }          
          else 
          {
               seq_read_seq[seq_swich] = table[seq_swich];
          }
      }
      int k_mer_order = 0;
      Read_kmer read_kmer[79];
      int i_dingwei = 0;
         while (i_dingwei<88)//不重叠的k_mer
         {
            /* code */
            unsigned char test_read[22];
            int flag_read = 0;
            for(int j_dingwei=0; j_dingwei<22; j_dingwei++)
            {
                test_read[j_dingwei] = seq_read_00[i_dingwei];
                i_dingwei++;
            }
            i_dingwei++;
            if (k_mer_order == 0)
            {
                /* code */
              int number=atoi(test_read);
              if(number % 7 == 5 || number % 17 == 7 || number % 19 == 13 || number % 53 == 31 || number % 71 == 47 )
              {strcpy(read_kmer[k_mer_order].k_mer_read, test_read); 
              read_kmer[k_mer_order].loction_read = i_dingwei-22;
              k_mer_order++;
              }
            }
             for (int k = 0; k < k_mer_order; k++)
             {
                /* code */
                if (strcmp(read_kmer[k].k_mer_read, test_read) == 0)
                {
                    /* code */
                    flag_read = 1;
                    

                }
                
             }
             if (flag_read == 0)
             {
                /* code */
                int number=atoi(test_read);
               if(number % 7 == 5 || number % 17 == 7 || number % 19 == 13 || number % 53 == 31 || number % 71 == 47 )
                {
                strcpy(read_kmer[k_mer_order].k_mer_read, test_read);
                read_kmer[k_mer_order].loction_read = i_dingwei-22;
                k_mer_order++;
                }
             }
             
         }

          for (int  kk = 1; kk < 79; kk++)//重叠k_mer
          {
            /* code */
               int j_2 =0 ;
               unsigned char test_read[22];
               int flag_read = 0;
               for(int j=kk; j<kk+22; j++)
               {                    
                  test_read[j_2] = seq_read_00[kk] ;
                  j_2++;
               }
             for (int k = 0; k < k_mer_order; k++)
             {
                /* code */
                if (strcmp(read_kmer[k].k_mer_read, test_read) == 0)
                {
                    /* code */
                    flag_read = 1;
                    

                }
              }
              if (flag_read == 0)
             {
                /* code */
               int number=atoi(test_read);
              if(number % 7 == 5 || number % 17 == 7 || number % 19 == 13 || number % 53 == 31 || number % 71 == 47 )
                {
                strcpy(read_kmer[k_mer_order].k_mer_read, test_read);
                read_kmer[k_mer_order].loction_read = kk;
                k_mer_order++;
                }
             }
          }
          less_max = (Less_max*)malloc(k_mer_order*sizeof(Less_max));
          over_max = (Over_max*)malloc(k_mer_order*sizeof(Over_max));

          for(int kmer_dingwei = 0; kmer_dingwei<k_mer_order; kmer_dingwei++)
          {         
             
             for (long p_dingwei = 0; p_dingwei < p_k_mer1; p_dingwei++)//在p数组中查找与k_mer相同的
            {
            /* code */
                      if(strcmp(read_kmer[kmer_dingwei].k_mer_read, p_yuanzu[p_dingwei].k_mer22) == 0)
                      {
                         int read_loctions = p_yuanzu[p_dingwei].end - p_yuanzu[p_dingwei].start + 1;
                         if (read_loctions<= MAX_LOC)
                         {
                            less_max[lessloction].start = p_yuanzu[p_dingwei].start;
                            less_max[lessloction].end = p_yuanzu[p_dingwei].end;
                            lessloction++;
                         }
                         else
                         {
                           over_max[overloction].start = p_yuanzu[p_dingwei].start;
                           over_max[overloction].end = p_yuanzu[p_dingwei].end;
                           overloction++;
                         }
                      }
            }

          }
       alg = hanming(seq_read_00, seqs, less_max, over_max, lessloction, overloction, w_yuanzu, length);
       fprintf(fp_sam, "%c\t", '*');
       fprintf(fp_sam, "%c\t", '*');
       fprintf(fp_sam, "%d\t", alg.ref_loc);
       fprintf(fp_sam, "%c\t", '*');
       fprintf(fp_sam, "%s\t", alg.cig_1);
       fprintf(fp_sam, "%c\t", '*');
       fprintf(fp_sam, "%c\t", '*');  
       fprintf(fp_sam, "%c\t", '*');
       fprintf(fp_sam, "%c\t", '*'); 
       fprintf(fp_sam, "%s\t", seq_read_seq); 
       fprintf(fp_sam, "%c\n", '*'); 

    }
    fclose(fp_sam);
    kseq_destroy(seq_read);
    gzclose(fp_raed);
}