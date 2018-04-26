#include "unittest.h"
#include "compute_mem_prob.c"

void
test_set_params_mem_prob
(void)
{

   int success;

   // Case 1. //
   success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert(success);

   test_assert(G == 17);
   test_assert(K == 50);
   test_assert(P == 0.01);
   test_assert(U == 0.05);

   test_assert(HIGH == 50);
   test_assert(KSZ == sizeof(trunc_pol_t) + 51*sizeof(double));

   for (int i = 0 ; i < MAXN ; i++) {
      test_assert(ARRAY[i] == NULL);
   }

   // Case 2. //
   set_params_mem_prob(50, 20, 0.99, 0.95);
   test_assert(success);

   test_assert(G == 50);
   test_assert(K == 20);
   test_assert(P == 0.99);
   test_assert(U == 0.95);

   test_assert(HIGH == 50);
   test_assert(KSZ == sizeof(trunc_pol_t) + 21*sizeof(double));

   for (int i = 0 ; i < MAXN ; i++) {
      test_assert(ARRAY[i] == NULL);
   }

   return;

}


void
test_error_set_params_mem_prob
(void)
{

   int success;

   // Case 1. //
   redirect_stderr();
   success = set_params_mem_prob(17, 50, 0.00, 0.05);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 2. //
   redirect_stderr();
   success = set_params_mem_prob(17, 50, 0.01, 0.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 3. //
   redirect_stderr();
   success = set_params_mem_prob(17, 50, 1.00, 0.05);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 4. //
   redirect_stderr();
   success = set_params_mem_prob(17, 50, 0.01, 1.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Test memory error //
   set_alloc_failure_rate_to(1.0);
   redirect_stderr();
   success = set_params_mem_prob(17, 50, 0.01, 0.05);
   unredirect_stderr();
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   return;

}


void
test_trunc_pol_mult
(void)
{

   test_assert(0);

   return;

}

// Test cases for export.
const test_case_t test_cases_compute_mem_prob[] = {
   {"set_params_mem_prob",       test_set_params_mem_prob},
   {"error_set_params_mem_prob", test_error_set_params_mem_prob},
   {"trunc_pol_mult",            test_trunc_pol_mult},
   {NULL, NULL},
};
