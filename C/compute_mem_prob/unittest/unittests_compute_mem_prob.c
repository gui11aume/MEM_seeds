#include "unittest.h"
#include "compute_mem_prob.c"

void
test_trunc_pol_mult
(void)
{

   test_assert(0);

   return;

}

// Test cases for export.
const test_case_t test_cases_compute_mem_prob[] = {
   {"compute_mem_prob/kpoly_mult",      test_trunc_pol_mult},
   {NULL, NULL},
};
