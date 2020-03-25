#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#include "mt.h"

#define GAMMA 19
#define K 100
#define N 100
#define skip 9

// Random number from 1 to 3 included.
#define randombp() (1 + (randomMT() / 1431655765))

void print (char * seq) {
  for (int i = 0 ; i < K ; i++) {
    fprintf(stdout, "%d", seq[i]);
  }
  fprintf(stdout, "\n");
}


int main(int argc, char **argv) {

  const double mu   = 0.06;
  const unsigned long int m = (mu * 4294967295);

  FILE * mutfile = fopen(argv[1], "r");
  if (mutfile == NULL) {
     fprintf(stderr, "cannot open file %s\n", argv[1]);
     exit(EXIT_FAILURE);
  }

  // Set the random seed.
  seedMT(123);

  char dup[N+1][K] = {0};
  char * read = dup[0];

  int  err[N+1] = {0};      // Total errors.
  int  str[N+1] = {0};      // Streak score.

  int  has_seed[N+1] = {0}; // Seeded duplicates.

  int total_case_1 = 0;
  int total_case_2 = 0;


  size_t sz = 151;
  ssize_t len;
  char * mut = malloc(sz);

  // Run the simulation.
  int ITER = 0;
  while ((len = getline(&mut, &sz, mutfile)) != -1) {

    ITER++;

    // Erase seed info.
    bzero(has_seed, (N+1) * sizeof(int));
    
    // Erase read and duplicates.
    bzero(dup, (N+1) * K);
    bzero(err, (N+1) * sizeof(int));
    bzero(str, (N+1) * sizeof(int));

    // Introduce errors in the read.
    bzero(read, K);
    int E = 0; // Number of errors.
    for (int i = 0 ; i < K ; i++) {
      if (mut[i] == '1') {
         read[i] = randombp();
         E++;
      }
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
