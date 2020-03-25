#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#include "mt.h"

#define ITER 1000000
#define GAMMA 19
#define K 100
#define N 100
#define skip 9

#define randombp() (1 + (randomMT() / 1431655765))

static int POS[151] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
   21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,
   45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,
   69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,
   93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,
   113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,
   131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,
   149,150};

int shuffle (const void * a, const void * b) {
  return (randomMT() < 2147483648) ? -1 : 1;
}

void print (char * seq) {
  for (int i = 0 ; i < K ; i++) {
    fprintf(stdout, "%d", seq[i]);
  }
  fprintf(stdout, "\n");
}


int main(int argc, char **argv) {

//  char **ignore;
//
  size_t E    = atoi(argv[1]);
//  double prob = strtod(argv[2], ignore); // p
//
  if (E == 0) {
    fprintf(stderr, "argument error\n");
    exit(EXIT_FAILURE);
  }

//  const double prob = 0.01;
  const double mu   = 0.06;

//  const unsigned long int p = (prob * 4294967295);
  const unsigned long int m = (mu   * 4294967295);

  // Set the random seed.
  seedMT(123);

  char dup[N+1][K] = {0};
  char * read = dup[0];

  int  err[N+1] = {0};      // Total errors.
  int  str[N+1] = {0};      // Streak score.

  int  has_seed[N+1] = {0}; // Seeded duplicates.

  int total_case_1 = 0;
  int total_case_2 = 0;

  // Run the simulation.
  for (size_t iter = 0 ; iter < ITER ; iter++) {

    // Erase seed info.
    bzero(has_seed, (N+1) * sizeof(int));
    
    // Erase read and duplicates.
    bzero(dup, (N+1) * K);
    bzero(err, (N+1) * sizeof(int));
    bzero(str, (N+1) * sizeof(int));

    // Get error positions in the read.
    qsort(POS, K, sizeof(int), shuffle);

    // Introduce E errors in the read.
    for (int e = 0 ; e < E ; e++) {
      read[POS[e]] = randombp();
    }

    for (int i = 0 ; i < K ; i++) {
      if (read[i] != 0) {
        str[0] = (i % (skip+1)) - skip;
        err[0]++;
      }
      else {
        str[0]++;
      }
      for (int n = 1 ; n < N+1; n++) {
        if (randomMT() < m) dup[n][i] = randombp();
        if (dup[n][i] != read[i]) {
          str[n] = (i % (skip+1)) - skip;
          err[n]++;
        }
        else {
          str[n]++;
        }
      }
      // Check seeds.
      for (int n = 0 ; n < N+1 ; n++) {
        if (str[n] >= GAMMA) has_seed[n] = 1;
      }
    }

    int there_is_a_better_hit = 0;
    for (int n = 1 ; n < N+1 ; n++) {
      if (err[n] < E && has_seed[n]) {
        there_is_a_better_hit = 1;
//           for (int n = 0 ; n < N+1 ; n++) {
//             print(dup[n]);
//           }
//           fprintf(stdout, "---\n");
        break;
      }
    }

    int has_false_hit = 0;
    for (int n = 0 ; n < N+1 ; n++) {
      if (has_seed[n]){
        has_false_hit = 1;
        break;
      }
    }

    total_case_1 += !has_seed[0] && has_false_hit;
    total_case_2 += has_seed[0] && there_is_a_better_hit;

  }

  fprintf(stdout, "Case 1: %f Case 2: %f\n",
      total_case_1 / (float) ITER, total_case_2 / (float) ITER);

}
