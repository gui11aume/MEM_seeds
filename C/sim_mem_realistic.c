#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#include "mt.h"

#define GAMMA 19
#define K 100
#define N 1

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

  int  falls[N+1] = {0};    // Which threads fall.
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

    int dominant = 0;

    for (int i = 0 ; i < K ; i++) {
      // See which threads fall.
      falls[0] = read[i] != 0;
      for (int n = 1 ; n < N+1; n++) {
        if (randomMT() < m) dup[n][i] = randombp();
        falls[n] = dup[n][i] != read[i];
      }
      // Update.
      int update_dominant = 0;
      if (falls[dominant]) {
        update_dominant = 1;
        if (str[dominant] >= GAMMA) {
          // It's seed time.
          int another_takes_over = 0;
          for (int n = 0 ; n < N+1 ; n++) {
            if (str[n] == str[dominant] && !falls[n]) {
              another_takes_over = 1;
              break;
            }
          }
          // If another thead takes over, there
          // is nothing to do. Otherwise, we have
          // a seed (strict or shared).
          if (!another_takes_over) {
            for (int n = 0 ; n < N+1 ; n++) {
              if (str[n] == str[dominant]) has_seed[n] = 1;
            }
          }
        }
      }
      // Update streaks and errors.
      for (int n = 0 ; n < N+1 ; n++) {
        if (falls[n]) {
          str[n] = 0;
          err[n]++;
        }
        else {
          str[n]++;
        }
      }
      // If required, update dominant.
      if (update_dominant) {
        dominant = 0;
        for (int n = 1 ; n < N+1 ; n++) {
          if (str[n] > str[dominant]) dominant = n;
        }
      }
    }

    // Final wrap up.
    if (str[dominant] >= GAMMA) {
      for (int n = 0 ; n < N+1 ; n++) {
        if (str[n] == str[dominant]) has_seed[n] = 1;
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
