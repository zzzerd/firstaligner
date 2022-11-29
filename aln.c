//read部分
#include <zlib.h>  
#include <stdio.h>
#include <string.h>  
#include <time.h>
#include "klib/kseq.h"  
#include "struct_ref_read.h"
KSEQ_INIT(gzFile, gzread)


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
      unsigned char seq_read_00 [100];
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
          for(int kmer_dingwei = 0; kmer_dingwei<k_mer_order; kmer_dingwei++)
          {         
             for (long p_dingwei = 0; p_dingwei < p_k_mer1; p_dingwei++)//在p数组中查找与k_mer相同的
            {
            /* code */
                      if(strcmp(read_kmer[kmer_dingwei].k_mer_read, p_yuanzu[p_dingwei].k_mer22) == 0)
                      {
                         
                      }
            }
          }
         
      
      

    }
    fclose(fp_sam);
    kseq_destroy(seq_read);
    gzclose(fp_raed);
}