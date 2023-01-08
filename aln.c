// read部分

/* struct Unit {
  int W1; // 是否往上回溯一格
  int W2; // 是否往左上回溯一格
  int W3; // 是否往左回溯一格
  int X;
  int Y;
  int M;
  int O; // 得分矩阵第(i, j)这个单元的分值，即序列s(1,...,i)与序列r(1,...,j)比对的最高得分
};
typedef struct Unit* pUnit; */
#include "aln.h"
int g(int l, int m) {
  if (l > (m + 2))
    return 1;
  else
    return 0;
}
void rle(char* cig) //游程编码
{
  int count = 0;
  char samechar;
  int length = 100;
  char* cig_rle = NULL;
  cig_rle = (char*)malloc(length * sizeof(char));
  char temp = cig[0];
  int loction = 0;
  int rle_lo = 0;
  samechar = temp;
  while (loction < 100) {
    /* code */

    if (temp == samechar) {
      count++;
      loction++;
      temp = cig[loction];
      continue;
    }
    else {
      cig_rle[rle_lo] = count;
      rle_lo++;
      cig_rle[rle_lo] = samechar;
      rle_lo++;
      samechar = temp;
      count = 0;
    }
  }
  memset(cig, '\0', 100);
  strncpy(cig, cig_rle, 100);
}

Alg_read hanming(unsigned char* seq_read_00, char* seqs, Less_max* less_max, Over_max* over_max, int lessloction, int overloction, W_yuanzu* w_yuanzu, int bases, int* n_cigar_,
  uint32_t** cigar_) {
  Alg_read alg;
  alg.min = 0;
  int l = 0;
  alg.ref_loc = 0;
  int flag_ov = 0;
  int base = bases;
  alg.cig_1 = (u_int32_t*)malloc(100 * sizeof(u_int32_t));
  for (int i = 0; i <= lessloction; i++) {
    int start = less_max[i].start;
    int end = less_max[i].end;
    char cig[100];
    int mismatch = 0;

    for (; start <= end; start++) {
      int alg_loc = w_yuanzu[start].Loction1;
      for (int j = 0; j < 100; j++) {
        /* code */
        if (seq_read_00[j] == seqs[alg_loc]) {
          /* code */
          cig[j] = '=';
        }
        else {
          cig[j] = 'X';
          mismatch++;
        }
        alg_loc++;
      }
      if (alg.min > mismatch) {
        /* code */
        alg.min = mismatch;
        rle(cig);
        memcpy(alg.cig_1, cig, sizeof(*alg.cig_1));
        alg.ref_loc = w_yuanzu[start].Loction1;
      }
    }
    l++;
    if (g(l, alg.min)) {
      flag_ov = 1;
      break;
    }
  }
  while (flag_ov == 0) {
    /* code */
    for (int i = 0; i <= overloction; i++) {
      int start = over_max[i].start;
      int end = over_max[i].end;
      char cig[100];
      int mismatch = 0;
      for (; start <= end; start++) {
        int alg_loc = w_yuanzu[start].Loction1;
        for (int j = 0; j < 100; j++) {
          /* code */
          if (seq_read_00[j] == seqs[alg_loc]) {
            /* code */
            cig[j] = '=';
          }
          else {
            cig[j] = 'X';
            mismatch++;
          }
          alg_loc++;
        }
        if (alg.min > mismatch) {
          /* code */
          alg.min = mismatch;
          rle(cig);

          alg.ref_loc = w_yuanzu[start].Loction1;
        }
      }
      l++;
      if (g(l, alg.min)) {
        break;
      }
    }
    flag_ov = 1;
  }

  if (alg.min > 50) {
    /* code */
    // NW函数
    /* return NW(seq_read_00, seqs, less_max, over_max, lessloction, overloction, w_yuanzu, base); */
    int sa = 1, sb = 3;
    int k;
    int gapo = 5, gape = 2;
    int8_t mat[16];

    for (int i = k = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j)
        mat[k++] = i == j ? sa : -sb;
    }
    ksw_global(100, seq_read_00, bases, seqs, 4, mat, gapo, gape, 100, n_cigar_, cigar_);
    alg.ref_loc = *n_cigar_;
    alg.cig_1 = *cigar_;
    return alg;
  }
  else
    return alg;
}
KSEQ_INIT(gzFile, gzread)
void aligner(int argc, char* argv[], unsigned char** p_yuanzu, W_yuanzu* w_yuanzu, char* seqs) {

  FILE* fp_sam = NULL;
  /* uint64_t p_k_mer1 = p_k_mer; */
  /* int i, a[6],b[6]; //i为数组元素下标 */
  fp_sam = fopen("result.sam", "w");
  if (fp_sam == NULL) {
    printf("Can't open the file!\n");
  }
  /* for (i = 0; i < 6; i++)
  {
      fprintf(fp, "%d ", a[i]);
  } */
  gzFile fp_raed;
  kseq_t* seq_read = NULL;
  int l_read;
  int length;
  if (argc == 1) {
    fprintf(stderr, "Usage: %s <in.fasta>\n", argv[0]);
  }

  fp_raed = gzopen(argv[2], "r");
  seq_read = kseq_init(fp_raed);

  while ((l_read = kseq_read(seq_read)) >= 0) {
    /*        printf("name: %s\n", seq->name.s);
           if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
           printf("seq: %s\n", seq->seq.s);
           if (seq->qual.l) printf("qual: %s\n", seq->qual.s); */
    fprintf(fp_sam, "%s\t", seq_read->name.s);
    unsigned char* seq_read_seq = seq_read->seq.s;
    length = strlen(seq_read->seq.s);
    Less_max* less_max = NULL;
    Over_max* over_max = NULL;
    int lessloction = 0;
    int overloction = 0;
    char* cigar = NULL;
    Alg_read alg;
    int* n_cigar_ = NULL;
    uint32_t** cigar_ = NULL;
    for (int seq_swich = 0; seq_swich < 100; seq_swich++) {
      /* code */
      if (table_zm[seq_read_seq[seq_swich]] == 4) {
        /* code */
        unsigned char table_1[4] = { 0, 1, 2, 3 };
        srand((unsigned)time(0));
        seq_read_seq[seq_swich] == table_1[rand() % 4];
      }
      else {
        seq_read_seq[seq_swich] = table_zm[seq_read_seq[seq_swich]];
      }
    }
    int k_mer_order = 0;
    Read_kmer read_kmer[readkmernum];
    int i_dingwei = 0;
    while (i_dingwei < 88) //不重叠的k_mer
    {
      /* code */
      unsigned char test_read[22];
      int flag_read = 0;
      for (int j_dingwei = 0; j_dingwei < 22; j_dingwei++) {
        test_read[j_dingwei] = seq_read_seq[i_dingwei];
        i_dingwei++;
      }
      i_dingwei++;
      if (k_mer_order == 0) {
        /* code */
        uint64_t number = kmernumber(test_read, 22);
        if (number % 7 == 5 || number % 17 == 7 || number % 19 == 13 || number % 53 == 31 || number % 71 == 47) {
          strncpy(read_kmer[k_mer_order].k_mer_read, test_read, 22);
          read_kmer[k_mer_order].loction_read = i_dingwei - 22;
          k_mer_order++;
        }
      }
      for (int k = 0; k < k_mer_order; k++) {
        /* code */
        if (strncmp(read_kmer[k].k_mer_read, test_read, 22) == 0) {
          /* code */
          flag_read = 1;
        }
      }
      if (flag_read == 0) {
        /* code */
        uint64_t number = kmernumber(test_read, 22);
        if (number % 7 == 5 || number % 17 == 7 || number % 19 == 13 || number % 53 == 31 || number % 71 == 47) {
          strncpy(read_kmer[k_mer_order].k_mer_read, test_read, 22);
          read_kmer[k_mer_order].loction_read = i_dingwei - 22;
          k_mer_order++;
        }
      }
    }

    for (int kk = 1; kk < 79; kk++) //重叠k_mer
    {
      /* code */
      int j_2 = 0;
      unsigned char test_read[22];
      int flag_read = 0;
      for (int j = kk; j < kk + 22; j++) {
        test_read[j_2] = seq_read_seq[kk];
        j_2++;
      }
      for (int k = 0; k < k_mer_order; k++) {
        /* code */
        if (strncmp(read_kmer[k].k_mer_read, test_read, 22) == 0) {
          /* code */
          flag_read = 1;
        }
      }
      if (flag_read == 0) {
        /* code */
        uint64_t number = kmernumber(test_read, 22);
        if (number % 7 == 5 || number % 17 == 7 || number % 19 == 13 || number % 53 == 31 || number % 71 == 47) {
          strncpy(read_kmer[k_mer_order].k_mer_read, test_read, 22);
          read_kmer[k_mer_order].loction_read = kk;
          k_mer_order++;
        }
      }
    }
    less_max = (Less_max*)malloc(k_mer_order * sizeof(Less_max));
    over_max = (Over_max*)malloc(k_mer_order * sizeof(Over_max));

    for (int kmer_dingwei = 0; kmer_dingwei < k_mer_order; kmer_dingwei++) {
      u_int64_t tmp = kmernumber(read_kmer[kmer_dingwei].k_mer_read, 22);
      if (p_yuanzu[tmp][0] != '\0') {
        if (p_yuanzu[tmp][0] <= MAX_LOC) {
          /* code */
          less_max[lessloction].start = p_yuanzu[tmp][1];
          less_max[lessloction].end = p_yuanzu[tmp][0] + p_yuanzu[tmp][1] - 1;
          lessloction++;
        }
        else {
          over_max[lessloction].start = p_yuanzu[tmp][1];
          over_max[lessloction].end = p_yuanzu[tmp][0] + p_yuanzu[tmp][1] - 1;
          overloction++;
        }

      }


    }
    alg = hanming(seq_read_seq, seqs, less_max, over_max, lessloction, overloction, w_yuanzu, length, n_cigar_,
      cigar_);
    fprintf(fp_sam, "%c\t", '*');
    fprintf(fp_sam, "%c\t", '*');
    fprintf(fp_sam, "%d\t", alg.ref_loc);
    fprintf(fp_sam, "%c\t", '*');
    fprintf(fp_sam, "%s\t", (char*)alg.cig_1);
    fprintf(fp_sam, "%c\t", '*');
    fprintf(fp_sam, "%c\t", '*');
    fprintf(fp_sam, "%c\t", '*');
    fprintf(fp_sam, "%c\t", '*');
    fprintf(fp_sam, "%s\t", seq_read_seq);
    fprintf(fp_sam, "%c\n", '*');
    free(over_max);
    free(less_max);
  }

  fclose(fp_sam);
  kseq_destroy(seq_read);
  gzclose(fp_raed);
}