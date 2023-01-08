#include "idx.h"
#define pair_lt(a, b) (kmernumber((a).k_mer,32) < kmernumber((b).k_mer,32) )
//归并排序
u_int64_t kmernumber(unsigned char test[], int length) {
  u_int64_t tmp = 0;
  int i = length - 1;
  for (; i >= 0; i--) {
    /* code */
    tmp <<= 2;
    tmp |= test[i];
  }
  return tmp;

}
/* void merge(Yuanzu* arr, uint64_t start, uint64_t mid, uint64_t end, uint64_t length) {
  Yuanzu result[length];
  uint64_t k = 0;
  uint64_t i = start;
  uint64_t j = mid + 1;
  while (i <= mid && j <= end) {
    if (strncmp(arr[i].k_mer, arr[j].k_mer, 32) < 0) {
      result[k++] = arr[i++];
    }
    else {
      result[k++] = arr[j++];
    }
  }
  if (i == mid + 1) {
    while (j <= end)
      result[k++] = arr[j++];
  }
  if (j == end + 1) {
    while (i <= mid)
      result[k++] = arr[i++];
  }
  for (j = 0, i = start; j < k; i++, j++) {
    arr[i] = result[j];
  }
}

void mergeSort(Yuanzu* arr, uint64_t start, uint64_t end, uint64_t length) {
  if (start >= end)
    return;
  uint64_t mid = (start + end) / 2;
  mergeSort(arr, start, mid, mid + 1);
  mergeSort(arr, mid + 1, end, mid + 1);
  merge(arr, start, mid, end, length);
} */
KSORT_INIT(pair, Yuanzu, pair_lt)
KSORT_INIT_GENERIC(long)
void ref(unsigned char** p_yuanzu, W_yuanzu* w_yuanzu, Yuanzu* yuanzu, char* seqs, uint64_t k_mers) {
  uint64_t k_mers_1 = 0; //使用哈希函数后的剩余量
  int flag;              // k-mer哈希后是否留下
  unsigned char test[lengthkmer];
  uint64_t i = 0;
  /*  uint64_t chushi = 0;
   uint64_t chushi2 = 0; */
  unsigned char test_22[22];
  uint64_t p_k_mer = 0; // p数组k-mer数量

  for (i = 0; i < k_mers; i++) {
    uint64_t j_1 = 0;
    for (uint64_t j = i; j < i + 32; j++) {
      test[j_1++] = seqs[j];
      /*  j_1++; */
    }
    u_int64_t number = kmernumber(test, 32);
    flag = 0;
    if (number % 7 == 5 || number % 17 == 7 || number % 19 == 13 || number % 53 == 31 || number % 71 == 47) {
      int k_1 = 0;
      for (uint64_t k = i; k < i + 32; k++) {
        yuanzu[k_mers_1].k_mer[k_1] = seqs[k];
        k_1++;
      }

      flag = 1;
    } //哈希函数
    if (flag == 1) {
      yuanzu[k_mers_1].Loction = i;
      k_mers_1++;
    }
  }
  //chushi2 = k_mers_1 - 1;
  //mergeSort(yuanzu, chushi, chushi2, k_mers_1); //排序基因组k-mer
  ks_mergesort(pair, k_mers_1, yuanzu, 0);
  uint64_t w_loction = 0;
  p_yuanzu = (unsigned char**)malloc(0x3FFFFF * sizeof(unsigned char));
  for (unsigned char i = 0;i <= 0x3FFFFF;i++) {
    p_yuanzu[i] = (unsigned char*)malloc(sizeof(unsigned char) * 2);
  }
  w_yuanzu = (W_yuanzu*)malloc(k_mers_1 * sizeof(W_yuanzu));
  /*     P_yuanzu p_yuanzu[k_mers_1];
      W_yuanzu w_yuanzu[k_mers_1]; */
  while (w_loction < k_mers_1) {
    /* code */
    uint64_t w_dingwei = 0; //用于w数组的k-mer
    unsigned char k_mer22[22];
    strncpy(k_mer22, yuanzu[w_loction].k_mer, 22);
    uint64_t loc = kmernumber(k_mer22, 22);
    u_int64_t now = w_loction;
    p_yuanzu[loc][1] = w_loction;
    for (int i_dingwei = 22; i_dingwei < 32; i_dingwei++) {
      /* code */
      w_yuanzu[w_loction].k_mer10[w_dingwei] = yuanzu[w_loction].k_mer[i_dingwei];
      w_dingwei++;
    }
    w_yuanzu[w_loction].Loction1 = yuanzu->Loction;
    w_loction++;
    unsigned char yuanzu_k_mer22[22];
    for (; strncmp(k_mer22, strncpy(yuanzu_k_mer22, yuanzu[w_loction].k_mer, 22), 22) == 0; w_loction++) {
      /* code */
      int w_dingwei2 = 0;
      for (int i_dingwei = 22; i_dingwei < 32; i_dingwei++) {
        /* code */
        w_yuanzu[w_loction].k_mer10[w_dingwei2] = yuanzu[w_loction].k_mer[i_dingwei];
        w_dingwei2++;
      }
    }
    p_yuanzu[loc][0] = w_loction - now;
  } // p数组，w数组
  free(yuanzu);
}