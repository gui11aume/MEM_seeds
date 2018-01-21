#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>

#include "mt.h"

// Global variables.
unsigned int K;
unsigned int GAMMA;
unsigned int NDUP;


void
update_dup
(
         unsigned int * dup,
   const unsigned long int mt_match
)
{
   for (int j = 0 ; j < NDUP ; j++)
      dup[j] = randomMT() < mt_match ? dup[j] + 1 : 0;
   return;
}


int
there_is_a_seed
(
   const unsigned int   target,
   const unsigned int * dup
)
{
   if (target < GAMMA) return 0;
   for (int j = 0 ; j < NDUP ; j++)
      if (dup[j] >= target) return 0;
   return 1;
}


int
simulate_read
(
   unsigned long int   mt_error,
   unsigned long int   mt_match1,
   unsigned long int   mt_match2,
   unsigned int      * dup
)
{
   int read_error;

   // Initialize counters.
   bzero(dup, NDUP * sizeof(int));
   int target = 0;

   for (int i = 0 ; i < K ; i++) {

      read_error = randomMT() < mt_error;
      update_dup(dup, read_error ? mt_match2 : mt_match1);

      // The code of 'there_is_a_seed()' works for both strict and
      // shared seeds, but we need to update 'target' before when
      // there is no error, and after when there is one.
      if (!read_error) target++;
      if (there_is_a_seed(target, dup)) return 1;
      if (read_error) target = 0;

   }

   // Check if there is a shared seed at the end of the read.
   if (target >= GAMMA && !read_error)
      return there_is_a_seed(target+1, dup);

   return 0;

}



int main(int argc, char **argv) {

                K     = atoi(argv[1]);            // Read size.
                GAMMA = atoi(argv[2]);            // Seed length.
   double       prob  = strtod(argv[3], NULL);    // Read error.
   double       kappa = strtod(argv[4], NULL);    // Divergence.
                NDUP  = atoi(argv[5]);            // Duplicate number.
   long int     ITER  = atoi(argv[6]);

   // Check arguments.
   int arguments_OK = K > 0 && GAMMA > 0 && ITER > 0 &&
         (prob >= 0 && prob <= 1.0) && (kappa >= 0 && kappa <= 1.0);
   if (!arguments_OK) {
      fprintf(stderr, "argument error\n");
      exit(EXIT_FAILURE);
   }

   // Convert probabilities to integers for Mersenne twister.
   const unsigned long int mt_error = (prob * 4294967295);
   const unsigned long int mt_match1 = ((1-kappa) * 4294967295);
   const unsigned long int mt_match2 = (kappa/3 * 4294967295);

   // Create data record for target and duplicates.
   unsigned int *dup = malloc(NDUP * sizeof(int));
   if (dup == NULL) {
      fprintf(stderr, "memory error\n");
      exit(EXIT_FAILURE);
   }

   // Set the random seed.
   seedMT(123);

   // Count the reads without seed.
   double total = ITER;

   // Run the simulation.
   for (long int iter = 0 ; iter < ITER ; iter++) {
      total -= simulate_read(mt_error, mt_match1, mt_match2, dup);
   }

   free(dup);
   fprintf(stdout, "%d\t%f\n", K, total / ITER);

}
