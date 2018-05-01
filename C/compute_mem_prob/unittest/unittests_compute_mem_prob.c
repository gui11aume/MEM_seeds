#include "unittest.h"
#include "compute_mem_prob.c"

void
test_set_params_mem_prob
(void)
{

   int success;

   // Case 1.
   success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   test_assert(G == 17);
   test_assert(K == 50);
   test_assert(P == 0.01);
   test_assert(U == 0.05);

   test_assert(HIGH == 50);
   test_assert(KSZ == sizeof(trunc_pol_t) + 51*sizeof(double));

   for (int i = 0 ; i < MAXN ; i++) {
      test_assert(ARRAY[i] == NULL);
   }

   // Case 2.
   set_params_mem_prob(50, 20, 0.99, 0.95);
   test_assert_critical(success);

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

   // Case 1.
   redirect_stderr();
   success = set_params_mem_prob(17, 50, 0.00, 0.05);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 2.
   redirect_stderr();
   success = set_params_mem_prob(17, 50, 0.01, 0.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 3.
   redirect_stderr();
   success = set_params_mem_prob(17, 50, 1.00, 0.05);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 4.
   redirect_stderr();
   success = set_params_mem_prob(17, 50, 0.01, 1.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 5.
   redirect_stderr();
   success = set_params_mem_prob(0, 50, 0.01, 1.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Case 6.
   redirect_stderr();
   success = set_params_mem_prob(17, 0, 0.01, 1.00);
   unredirect_stderr();
   test_assert(!success);
   test_assert_stderr("[compute_mem_prob] error in function `set_");

   // Test memory error.
   set_alloc_failure_rate_to(1.0);
   redirect_stderr();
   success = set_params_mem_prob(17, 50, 0.01, 0.05);
   unredirect_stderr();
   reset_alloc();
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
   const double omega = .01 * pow(1-.05/3,2);
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
   const double _omega = .01 * (1-pow(1-.05/3,2));
   for (int i = 1 ; i <= 17 ; i++) {
      double target = _omega * pow(.99, i-1) * (1-pow(1-pow(.95,i-1),2)); 
      test_assert(fabs(_A->coeff[i]-target) < 1e-9);
   }
   for (int i = 18 ; i <= 50 ; i++) {
      double alpha_i_sq = pow(1-pow(.95,i-1) * .05/3,2);
      double target = 0.01 * pow(.99, i-1) * (1-alpha_i_sq);
      test_assert(fabs(_A->coeff[i]-target) < 1e-9);
   }

   // Test special case N = 0.
   trunc_pol_t *A0 = new_trunc_pol_A(50, 0, NO);
   test_assert_critical(A0 != NULL);
   test_assert(A0->mono.deg == 1);
   test_assert(A0->mono.coeff == .01);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(A0->coeff[i] == (i == 1 ? 0.01 : 0));
   }

   trunc_pol_t *_A0 = new_trunc_pol_A(50, 0, YES);
   test_assert_critical(_A0 != NULL);
   test_assert(_A0->mono.deg == 0);
   test_assert(_A0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(_A0->coeff[i] == 0);
   }

   free(A);
   free(_A);
   free(A0);
   free(_A0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_A
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *A;

   redirect_stderr();
   A = new_trunc_pol_A(0, 0, NO);
   unredirect_stderr();

   test_assert(A == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   A = new_trunc_pol_A(0, 2, NO);
   unredirect_stderr();

   test_assert(A == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_A(51, 2, NO);
   unredirect_stderr();

   test_assert(A == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   new_trunc_pol_A(18, 2, NO);
   unredirect_stderr();
   reset_alloc();

   test_assert(A == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_B
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert(success);

   trunc_pol_t *B = new_trunc_pol_B(50, 2, NO);
   test_assert_critical(B != NULL);
   test_assert(B->mono.deg == 0);
   test_assert(B->mono.coeff == 0);
   test_assert(B->coeff[0] == 0);
   const double omega = .01 * pow(1-.05/3,2);
   const double denom = 1-pow(1-.05/3,2);
   for (int i = 1 ; i <= 50 ; i++) {
      double alpha_i_sq = pow(1-pow(.95,i-1) * .05/3,2);
      double target = omega * pow(.99, i-1) * (1-alpha_i_sq) / denom;
      test_assert(fabs(B->coeff[i]-target) < 1e-9);
   }

   trunc_pol_t *_B = new_trunc_pol_B(50, 2, YES);
   test_assert_critical(_B != NULL);
   test_assert(_B->mono.deg == 0);
   test_assert(_B->mono.coeff == 0);
   test_assert(_B->coeff[0] == 0);
   const double _omega = .01 * (1-pow(1-.05/3,2));
   for (int i = 1 ; i <= 50 ; i++) {
      double alpha_i_sq = pow(1-pow(.95,i-1) * .05/3,2);
      double target = _omega * pow(.99, i-1) * (1-alpha_i_sq) / denom; 
      test_assert(fabs(_B->coeff[i]-target) < 1e-9);
   }

   // Test the special case N = 0.
   trunc_pol_t *B0 = new_trunc_pol_B(50, 0, NO);
   test_assert_critical(B0 != NULL);
   test_assert(B0->mono.deg == 0);
   test_assert(B0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(B0->coeff[i] == 0);
   }

   trunc_pol_t *_B0 = new_trunc_pol_B(50, 0, YES);
   test_assert_critical(_B0 != NULL);
   test_assert(_B0->mono.deg == 0);
   test_assert(_B0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(_B0->coeff[i] == 0);
   }

   free(B);
   free(_B);
   free(B0);
   free(_B0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_B
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *B;

   redirect_stderr();
   B = new_trunc_pol_B(0, 0, NO);
   unredirect_stderr();

   test_assert(B == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   B = new_trunc_pol_B(0, 2, NO);
   unredirect_stderr();

   test_assert(B == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_B(51, 2, NO);
   unredirect_stderr();

   test_assert(B == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   new_trunc_pol_B(18, 2, NO);
   unredirect_stderr();
   reset_alloc();

   test_assert(B == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_C
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a C polynomial of degree 10 with N = 2.
   trunc_pol_t *C = new_trunc_pol_C(10, 2, NO);
   test_assert_critical(C != NULL);
   test_assert(C->mono.deg == 0);
   test_assert(C->mono.coeff == 0);
   test_assert(C->coeff[0] == 0);
   const size_t j = 7; // G-10
   const double omega = .01 * pow(1-.05/3,2);
   const double denom = pow(1-pow(1-.05,j)*.05/3,2) - \
      pow(1-pow(1-.05,j-1)*.05/3,2) - pow(1-pow(1-.05,j),2) + \
      pow(1-(1-.05+.05*.05/3)*pow(1-.05,j-1),2);
   for (int i = 1 ; i <= 10 ; i++) {
      double num =  pow(1-pow(1-.05,j)*.05/3,2) - \
         pow(1-pow(1-.05,j-1)*.05/3,2) - \
         pow(1-pow(1-.05,j)*.05/3-pow(1-.05,i+j-1)*(1-.05/3),2) + \
         pow(1-pow(1-.05,j-1)*.05/3-pow(1-.05,i+j-1)*(1-.05/3),2);
      double target = omega * pow(.99, i-1) * num / denom;
      test_assert(fabs(C->coeff[i]-target) < 1e-9);
   }
   for (int i = 11 ; i <= 50 ; i++) {
      test_assert(C->coeff[i] == 0);
   }

   trunc_pol_t *_C = new_trunc_pol_C(10, 2, YES);
   test_assert_critical(_C != NULL);
   test_assert(_C->mono.deg == 0);
   test_assert(_C->mono.coeff == 0);
   test_assert(_C->coeff[0] == 0);
   const double _omega = .01 * (1 - pow(1-.05/3,2));
   for (int i = 1 ; i <= 10 ; i++) {
      double num =  pow(1-pow(1-.05,j)*.05/3,2) - \
         pow(1-pow(1-.05,j-1)*.05/3,2) - \
         pow(1-pow(1-.05,j)*.05/3-pow(1-.05,i+j-1)*(1-.05/3),2) + \
         pow(1-pow(1-.05,j-1)*.05/3-pow(1-.05,i+j-1)*(1-.05/3),2);
      double target = _omega * pow(.99, i-1) * num / denom;
      test_assert(fabs(_C->coeff[i]-target) < 1e-9);
   }
   for (int i = 11 ; i <= 50 ; i++) {
      test_assert(_C->coeff[i] == 0);
   }

   // Test the special cases N = 0 and N = 1
   trunc_pol_t *C0 = new_trunc_pol_C(50, 0, NO);
   trunc_pol_t *_C0 = new_trunc_pol_C(50, 0, YES);
   trunc_pol_t *C1 = new_trunc_pol_C(50, 1, NO);
   trunc_pol_t *_C1 = new_trunc_pol_C(50, 1, YES);
   test_assert_critical(C0 != NULL);
   test_assert_critical(_C0 != NULL);
   test_assert_critical(C1 != NULL);
   test_assert_critical(_C0 != NULL);
   test_assert(C0->mono.deg == 0);
   test_assert(C0->mono.coeff == 0);
   test_assert(_C0->mono.deg == 0);
   test_assert(_C0->mono.coeff == 0);
   test_assert(C1->mono.deg == 0);
   test_assert(C1->mono.coeff == 0);
   test_assert(_C1->mono.deg == 0);
   test_assert(_C1->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(C0->coeff[i] == 0);
      test_assert(_C0->coeff[i] == 0);
      test_assert(C1->coeff[i] == 0);
      test_assert(_C1->coeff[i] == 0);
   }

   free(C);
   free(_C);
   free(C0);
   free(_C0);
   free(C1);
   free(_C1);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_C
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *C;

   redirect_stderr();
   C = new_trunc_pol_C(0, 0, NO);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   C = new_trunc_pol_C(0, 2, NO);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_C(51, 2, NO);
   unredirect_stderr();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   new_trunc_pol_C(18, 2, NO);
   unredirect_stderr();
   reset_alloc();

   test_assert(C == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_D
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a D polynomial of degree 10 with N = 2.
   trunc_pol_t *D = new_trunc_pol_D(10, 2, NO);
   test_assert_critical(D != NULL);
   test_assert(D->mono.deg == 0);
   test_assert(D->mono.coeff == 0);
   test_assert(D->coeff[0] == 0);
   const double omega = .01 * pow(1-.05/3,2);
   for (int i = 1 ; i <= 10 ; i++) {
      double target = omega * pow(.99, i-1);
      test_assert(fabs(D->coeff[i]-target) < 1e-9);
   }
   for (int i = 11 ; i <= 50 ; i++) {
      test_assert(D->coeff[i] == 0);
   }

   trunc_pol_t *_D = new_trunc_pol_D(10, 2, YES);
   test_assert_critical(_D != NULL);
   test_assert(_D->mono.deg == 0);
   test_assert(_D->mono.coeff == 0);
   test_assert(_D->coeff[0] == 0);
   const double _omega = .01 * (1 - pow(1-.05/3,2));
   for (int i = 1 ; i <= 10 ; i++) {
      double target = _omega * pow(.99, i-1);
      test_assert(fabs(_D->coeff[i]-target) < 1e-9);
   }
   for (int i = 11 ; i <= 50 ; i++) {
      test_assert(_D->coeff[i] == 0);
   }

   // Test the special cases N = 0.
   trunc_pol_t *D0 = new_trunc_pol_D(50, 0, NO);
   trunc_pol_t *_D0 = new_trunc_pol_D(50, 0, YES);
   test_assert_critical(D0 != NULL);
   test_assert_critical(_D0 != NULL);
   test_assert(D0->mono.deg == 0);
   test_assert(D0->mono.coeff == 0);
   test_assert(_D0->mono.deg == 0);
   test_assert(_D0->mono.coeff == 0);

   test_assert(D0->coeff[0] == 0);
   test_assert(_D0->coeff[0] == 0);
   for (int i = 1 ; i <= 50 ; i++) {
      double target = pow(1-.01,i-1) * 0.01;
      test_assert(fabs(D0->coeff[i]-target) < 1e-9);
      test_assert(_D0->coeff[i] == 0);
   }

   free(D);
   free(_D);
   free(D0);
   free(_D0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_D
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *D;

   redirect_stderr();
   D = new_trunc_pol_D(0, 0, NO);
   unredirect_stderr();

   test_assert(D == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   D = new_trunc_pol_D(0, 2, NO);
   unredirect_stderr();

   test_assert(D == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_D(51, 2, NO);
   unredirect_stderr();

   test_assert(D == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   new_trunc_pol_D(18, 2, NO);
   unredirect_stderr();
   reset_alloc();

   test_assert(D == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_u
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a u polynomial of degree 10 with N = 2.
   trunc_pol_t *u = new_trunc_pol_u(10, 2);
   test_assert_critical(u != NULL);

   const double target = pow(1-.01,10) * \
      (pow(1-pow(1-.05,10),2) - pow(1-pow(1-.05,9),2));
   test_assert(u->mono.deg == 10);
   test_assert(fabs(u->mono.coeff-target) < 1e-9);
   for (int i = 0 ; i <= 50 ; i++) {
      if (i == 10)
         test_assert(fabs(u->coeff[i]-target) < 1e-9);
      else
         test_assert(u->coeff[i] == 0);
   }

   // Test the special case N = 0.
   trunc_pol_t *u0 = new_trunc_pol_u(1, 0);
   test_assert_critical(u0 != NULL);
   test_assert(u0->mono.deg == 1);
   test_assert(fabs(u0->mono.coeff-.99) < 1e-9);
   for (int i = 0 ; i <= 50 ; i++) {
      if (i == 1)
         test_assert(fabs(u0->coeff[i]-.99) < 1e-9);
      else
         test_assert(u0->coeff[i] == 0);
   }
   
   free(u);
   free(u0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_u
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *u;

   redirect_stderr();
   u = new_trunc_pol_u(0, 0);
   unredirect_stderr();

   test_assert(u == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   u = new_trunc_pol_u(0, 2);
   unredirect_stderr();

   test_assert(u == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_u(18, 2);
   unredirect_stderr();

   test_assert(u == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   new_trunc_pol_u(10, 2);
   unredirect_stderr();
   reset_alloc();

   test_assert(u == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_new_trunc_pol_v
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // Test a u polynomial of degree 10 with N = 2.
   trunc_pol_t *v = new_trunc_pol_v(10, 2);
   test_assert_critical(v != NULL);

   const double denom = 1-pow(1-.05/3,2);
   const double num = pow(1-pow(1-.05,10)*.05/3,2) - \
      pow(1-pow(1-.05,9)*.05/3,2) - pow(1-pow(1-.05,10),2) + \
      pow(1-(1-.05+.05*.05/3)*pow(1-.05,9),2);
   const double target = pow(1-.01,10) * num / denom;
   test_assert(v->mono.deg == 10);
   test_assert(fabs(v->mono.coeff-target) < 1e-9);
   for (int i = 0 ; i <= 50 ; i++) {
      if (i == 10)
         test_assert(fabs(v->coeff[i]-target) < 1e-9);
      else
         test_assert(v->coeff[i] == 0);
   }

   // Test the special case N = 0.
   trunc_pol_t *v0 = new_trunc_pol_v(1, 0);
   test_assert_critical(v0 != NULL);
   test_assert(v0->mono.deg == 0);
   test_assert(v0->mono.coeff == 0);
   for (int i = 0 ; i <= 50 ; i++) {
      test_assert(v0->coeff[i] == 0);
   }
   
   free(v);
   free(v0);
   clean_mem_prob();

}


void
test_error_new_trunc_pol_v
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   trunc_pol_t *v;

   redirect_stderr();
   v = new_trunc_pol_v(0, 0);
   unredirect_stderr();

   test_assert(v == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   v = new_trunc_pol_v(0, 2);
   unredirect_stderr();

   test_assert(v == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   redirect_stderr();
   new_trunc_pol_v(51, 2);
   unredirect_stderr();

   test_assert(v == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_t");

   set_alloc_failure_rate_to(1);
   redirect_stderr();
   new_trunc_pol_v(10, 2);
   unredirect_stderr();
   reset_alloc();

   test_assert(v == NULL);
   test_assert_stderr("[compute_mem_prob] error in function `new_z");

   clean_mem_prob();

}


void
test_compute_mem_prob
(void)
{

   int success = set_params_mem_prob(17, 50, 0.01, 0.05);
   test_assert_critical(success);

   // The first terms can be computed directly.
   for (int i = 0 ; i < 17 ; i++) {
      test_assert(fabs(compute_mem_prob(2,i)-1) < 1e-9);
   }
   double target_17 = 1-pow(.99,17);
   test_assert(fabs(compute_mem_prob(2,17)-target_17) < 1e-9);

   double target_18;
   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,2);
   test_assert(fabs(compute_mem_prob(2,18)-target_18) < 1e-9);

   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,3);
   test_assert(fabs(compute_mem_prob(3,18)-target_18) < 1e-9);

   target_18 = 1-pow(.99,18) - \
      2*.01*pow(.99,17) * pow(1-pow(.95,17)*.05/3,4);
   test_assert(fabs(compute_mem_prob(4,18)-target_18) < 1e-9);


   success = set_params_mem_prob(20, 50, 0.02, 0.05);
   test_assert_critical(success);

   // The first terms can be computed directly.
   for (int i = 0 ; i < 20 ; i++) {
      test_assert(fabs(compute_mem_prob(0,i)-1) < 1e-9);
      test_assert(fabs(compute_mem_prob(1,i)-1) < 1e-9);
      test_assert(fabs(compute_mem_prob(2,i)-1) < 1e-9);
   }

   const double target_20 = 1-pow(.98,20);
   test_assert(fabs(compute_mem_prob(0,20)-target_20) < 1e-9);
   test_assert(fabs(compute_mem_prob(1,20)-target_20) < 1e-9);
   test_assert(fabs(compute_mem_prob(2,20)-target_20) < 1e-9);
   test_assert(fabs(compute_mem_prob(3,20)-target_20) < 1e-9);

   double target_21;

   // Special case N = 0.
   target_21 = 1-pow(.98,21) - 2*.02*pow(.98,20);
   test_assert(fabs(compute_mem_prob(0,21)-target_21) < 1e-9);

   // Special case N = 1.
   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * (1-pow(.95,20)*.05/3);
   test_assert(fabs(compute_mem_prob(1,21)-target_21) < 1e-9);

   // Cases N > 1.
   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,2);
   test_assert(fabs(compute_mem_prob(2,21)-target_21) < 1e-9);

   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,3);
   test_assert(fabs(compute_mem_prob(3,21)-target_21) < 1e-9);

   target_21 = 1-pow(.98,21) - \
      2*.02*pow(.98,20) * pow(1-pow(.95,20)*.05/3,4);
   test_assert(fabs(compute_mem_prob(4,21)-target_21) < 1e-9);

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
   {"new_trunc_pol_B",           test_new_trunc_pol_B},
   {"error_new_trunc_pol_B",     test_error_new_trunc_pol_B},
   {"new_trunc_pol_C",           test_new_trunc_pol_C},
   {"new_error_trunc_pol_C",     test_error_new_trunc_pol_C},
   {"new_trunc_pol_D",           test_new_trunc_pol_D},
   {"new_error_trunc_pol_D",     test_error_new_trunc_pol_D},
   {"new_trunc_pol_u",           test_new_trunc_pol_u},
   {"new_error_trunc_pol_u",     test_error_new_trunc_pol_u},
   {"new_trunc_pol_v",           test_new_trunc_pol_v},
   {"new_error_trunc_pol_v",     test_error_new_trunc_pol_v},
//   {"compute_mem_prob",          test_compute_mem_prob},
   {NULL, NULL},
};
