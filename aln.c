//read部分
#include <zlib.h>  
#include <stdio.h>
#include <string.h>  
#include <time.h>
#include "klib/kseq.h"  
#include "struct_ref_read.h"
#define MAX_LOC 30
KSEQ_INIT(gzFile, gzread)
struct 
{
   /* data */
   int min;
   char cig_1[100];
   int ref_loc;
}Alg_read;

struct 
{
  /* data */
  int start;
  int end;

}Less_max;
struct 
{
  /* data */
  int start;
  int end;

}Over_max;
int g(int l , int m)
{
   if(l > (m+2))
    return 1 ;
    else
    return 0 ;
}
void rle()//游程编码
{
 int count = 0;
 char samechar;
 char temp = cig[0];
 int loction = 0;
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
      ;
      ;
      samechar = temp;
      count = 0;
   }
 }
 
}
Alg_read hanming( char *seq_read_00, char *seqs, Less_max *less_max, Over_max *over_max, int lessloction, int overloction, W_yuanzu *w_yuanzu)
{
   Alg_read alg;
   alg.min = 0;
   int l = 0;
   alg.ref_loc=0;
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
            cig[j]= "=";
         }
         else
         {
            cig[j]= "X";
            mismatch++;
         }
         alg_loc++;
       }
       if (min > mismatch)
       {
         /* code */
         min = mismatch;
         memcpy(alg.cig_1,cig,sizeof(cig_1));
         ref_loc = w_yuanzu[start].Loction1;
       }
       
    }
    l++;
    if(g(l,min))
      break;
   }
   if (min>50)
   {
      /* code */
      //NW函数
    return  ;
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
      char *cigar;
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
              strcpy(read_kmer[k_mer_order].k_mer_read, test_read); 
              read_kmer[k_mer_order].loction_read = i_dingwei-22;
              k_mer_order++;
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
                strcpy(read_kmer[k_mer_order].k_mer_read, test_read);
                read_kmer[k_mer_order].loction_read = i_dingwei-22;
                k_mer_order++;
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
                strcpy(read_kmer[k_mer_order].k_mer_read, test_read);
                read_kmer[k_mer_order].loction_read = kk;
                k_mer_order++;
             }
          }
          less_max = (Less_max*)malloc(k_mer_order*sizeof(Less_max));
          Over_max = (Over_max*)malloc(k_mer_order*sizeof(Over_max));
          int lessloction = 0 ;
          int overloction = 0 ;
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

         
      
      

    }
    fclose(fp_sam);
    kseq_destroy(seq_read);
    gzclose(fp_raed);
}