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

   clean_mem_prob();

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

   // Case 5. //
   redirect_stderr();
   success = set_params_mem_prob(0, 50, 0.01, 1.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 6. //
   redirect_stderr();
   success = set_params_mem_prob(17, 0, 0.01, 1.00);
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
test_uninitialized_error
(void)
{

   // Do not call 'set_params_mem_prob()'.
   redirect_stderr();
   double x = compute_mem_prob(5, 20);
   unredirect_stderr();
   test_assert_stderr("[compute_mem_prob] error in function `comp");
   test_assert(x != x);

   return;

}


void
test_new_zero_trunc_pol
(void)
{

   size_t ksz = 50;
   int success = set_params_mem_prob(17, ksz, 0.01, 0.05);
   test_assert(success);

   trunc_pol_t *a = new_zero_trunc_pol();
   test_assert_critical(a);

   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   free(a);
   clean_mem_prob();

   return;

}


void
test_trunc_pol_mult
(void)
{
   
   // Test multiplications between zero polynomials.

   size_t ksz = 50;
   set_params_mem_prob(17, ksz, 0.01, 0.05);
   trunc_pol_t *a = new_zero_trunc_pol();

   test_assert_critical(a != NULL);

   test_assert(trunc_pol_mult(a, NULL, NULL) == NULL);

   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   trunc_pol_t *b = new_zero_trunc_pol();
   test_assert_critical(b != NULL);

   test_assert(trunc_pol_mult(a, b, NULL) == NULL);

   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   test_assert(trunc_pol_mult(a, NULL, b) == NULL);

   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(a->coeff[i] == 0);
   }

   trunc_pol_t *c = new_zero_trunc_pol();
   test_assert_critical(c != NULL);

   // Returns 'a' if arguments are not NULL.
   test_assert(trunc_pol_mult(a, b, c) == a);

   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      test_assert(a->coeff[i] == 0);
   }


   // Test multiplications between monomials (b = 5z and c = z^2).
   b->mono.deg = 1; b->mono.coeff = 5; b->coeff[1] = 5;
   c->mono.deg = 2; c->mono.coeff = 1; c->coeff[2] = 1;

   test_assert(trunc_pol_mult(a, b, c) == a);
   test_assert(a->mono.deg == 3);
   test_assert(a->mono.coeff == 5);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = i == 3 ? 5 : 0;
      test_assert(a->coeff[i] == target);
   }


   // Test multiplications between a monomial and a
   // polynomial (b = 5z and c = z^2 + 2z^3).
   b->mono.deg = 1; b->mono.coeff = 5; b->coeff[1] = 5;
   c->mono.deg = 0; c->mono.coeff = 0; c->coeff[2] = 1; c->coeff[3] = 2;

   test_assert(trunc_pol_mult(a, b, c) == a);
   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 10;
      test_assert(a->coeff[i] == target);
   }

   // Test symmetry.
   test_assert(trunc_pol_mult(a, c, b) == a);
   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 10;
      test_assert(a->coeff[i] == target);
   }

   // Test multiplications between two polynomials
   // (b = 5z + 3z^2 and c = z^2 + 2z^3).
   b->mono.deg = 0; b->mono.coeff = 0; b->coeff[1] = 5; b->coeff[2] = 3;
   c->mono.deg = 0; c->mono.coeff = 0; c->coeff[2] = 1; c->coeff[3] = 2;

   test_assert(trunc_pol_mult(a, b, c) == a);
   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 13;
      if (i == 5) target = 6;
      test_assert(a->coeff[i] == target);
   }

   // Test symmetry.
   test_assert(trunc_pol_mult(a, c, b) == a);
   test_assert(a->mono.deg == 0);
   test_assert(a->mono.coeff == 0);
   for (int i = 0 ; i <= ksz ; i++) {
      double target = 0;
      if (i == 3) target = 5;
      if (i == 4) target = 13;
      if (i == 5) target = 6;
      test_assert(a->coeff[i] == target);
   }

   free(a);
   free(b);
   free(c);
   clean_mem_prob();

   return;

}


void
test_new_trunc_pol_A
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert(success);

   trunc_pol_t *A = new_trunc_pol_A(17, 2, NO);
   test_assert_critical(A != NULL);
   test_assert(A->mono.deg == 0);
   test_assert(A->mono.coeff == 0);
   test_assert(A->coeff[0] == 0);
   double omega = .01 * pow(1-.05/3,2);
   for (int i = 1 ; i <= 17 ; i++) {
      double target = omega * pow(.99, i-1) * (1-pow(1-pow(.95,i-1),2)); 
      test_assert(fabs(A->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= 50 ; i++) {
      test_assert(A->coeff[i] == 0);
   }

   trunc_pol_t *_A = new_trunc_pol_A(50, 2, YES);
   test_assert_critical(_A != NULL);
   test_assert(_A->mono.deg == 0);
   test_assert(_A->mono.coeff == 0);
   test_assert(_A->coeff[0] == 0);
   double _omega = .01 * (1-pow(1-.05/3,2));
   for (int i = 1 ; i <= 17 ; i++) {
      double target = _omega * pow(.99, i-1) * (1-pow(1-pow(.95,i-1),2)); 
      test_assert(fabs(_A->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= 50 ; i++) {
      double alpha_i_sq = pow(1-pow(.95,i-1) * .05/3,2);
      double target = 0.01 * pow(.99, i-1) * (1-alpha_i_sq);
      test_assert(fabs(_A->coeff[i]-target) < 1e-9);
   }

   free(A);
   free(_A);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_A
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   redirect_stderr();
   new_trunc_pol_A(0, 2, YES);
   unredirect_stderr();

   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_A(51, 2, YES);
   unredirect_stderr();

   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   clean_mem_prob();

}



// Test cases for export.
const test_case_t test_cases_compute_mem_prob[] = {
   {"set_params_mem_prob",       test_set_params_mem_prob},
   {"error_set_params_mem_prob", test_error_set_params_mem_prob},
   {"uninitialized_error",       test_uninitialized_error},
   {"new_zero_trunc_pol",        test_new_zero_trunc_pol},
   {"trunc_pol_mult",            test_trunc_pol_mult},
   {"new_trunc_pol_A",           test_new_trunc_pol_A},
   {"error_new_trunc_pol_A",     test_error_new_trunc_pol_A},
   {NULL, NULL},
};
