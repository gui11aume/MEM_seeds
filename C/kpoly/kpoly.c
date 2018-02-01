#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>


typedef unsigned int uint;

typedef struct kpoly_t kpoly_t;
typedef struct mat_t mat_t;
typedef struct mono_t mono_t;

struct mono_t {
  size_t deg;
  double coeff;
};

struct kpoly_t {
  const size_t k;            // Max degree of the polynomial.
        mono_t mono;         // Monomer (if applicable).
        double coeff[];      // Terms of the polynomial.
};

struct mat_t {
  const size_t    dim;       // Column / row number.
  const size_t    k;         // Max degree of the polynomials.
        kpoly_t * term[];    // Terms of the matrix.
};



void critical_error (char *, int) __attribute__ ((noreturn));


void print_kpoly (kpoly_t *P) {
   if (P == NULL) {
      fprintf(stderr, "0\n");
      return;
   }
   fprintf(stderr, "%f", P->coeff[0]);
   fprintf(stderr, " + %fz", P->coeff[1]);
   for (int i = 2 ; i <= P->k ; i++) {
      fprintf(stderr, " + %fz^%d", P->coeff[i], i);
   }
   fprintf(stderr, "\n");
}


void
critical_error
(
   char * msg,
   int    lineno
)
{
   fprintf(stderr, "error line %d: %s\n", lineno, msg);
   abort();
}


void
destroy_mat
(
   mat_t *mat
)
{

   size_t dim = mat->dim;
   for (int i = 0 ; i < dim*dim ; i++) free(mat->term[i]);
   free(mat);

}


kpoly_t *
new_zero_kpoly
(
   size_t k
)
{

   // Initialize to zero.
   kpoly_t *new = calloc(1, sizeof(kpoly_t) + (k+1) * sizeof(double));

   if (new == NULL) critical_error("memory error", __LINE__);

   // The size of the k-polynomial is set upon
   // creation and must never change afterwards.
   *(size_t *)&new->k = k;

   return new;

}


mat_t *
new_null_matrix
(
   size_t k,
   size_t dim
)
{

   // Initialize to zero.
   mat_t *new = calloc(1, sizeof(mat_t) + dim*dim * sizeof(kpoly_t *));

   if (new == NULL) critical_error("memory error", __LINE__);

   // The size of the k-polynomials and the dimension
   // are set upon creation and must never change afterwards.
   *(size_t *)&new->dim = dim;
   *(size_t *)&new->k = k;

   return new;

}


kpoly_t *
new_kpoly_A
(
   size_t k,    // Max degree.
   size_t n,    // Degree of A.
   double p,    // Probability of error.
   double mu,   // Divergence rate.
   uint   N     // Number of duplicates.
)
{

   if (n > k || n == 0) critical_error("parameter error", __LINE__);

   kpoly_t *new = new_zero_kpoly(k);

   double cst = p * (1-pow(1-mu/3.0, N));
   double q_power_m = 1.0;
   for (int m = 1 ; m <= n ; m++) {
      new->coeff[m] = cst * q_power_m;
      q_power_m *= (1-p);
   }

   return new;

}


kpoly_t *
new_kpoly_tilde_A
(
   size_t k,    // Max degree.
   size_t n,    // Degree of A.
   double p,    // Probability of error.
   double mu,   // Divergence rate.
   uint   N     // Number of duplicates.
)
{

   if (n > k || n == 0) critical_error("parameter error", __LINE__);

   kpoly_t *new = new_zero_kpoly(k);

   double cst = p * (1-pow(1-mu/3.0, N));
   double one_minus_mu_power_m = 1.0;
   double q_power_m = 1.0;
   for (int m = 1 ; m <= n ; m++) {
      double Cm = 1 - pow(1-one_minus_mu_power_m, N);
      new->coeff[m] = cst * Cm * q_power_m;
      one_minus_mu_power_m *= (1-mu);
      q_power_m *= (1-p);
   }

   return new;

}


kpoly_t *
new_kpoly_B
(
   size_t k,    // Max degree.
   size_t n,    // Degree of A.
   double p,    // Probability of error.
   double mu,   // Divergence rate.
   uint   N     // Number of duplicates.
)
{

   if (n > k || n == 0) critical_error("parameter error", __LINE__);

   kpoly_t *new = new_zero_kpoly(k);

   double cst = p * pow(1-mu/3.0, N);
   double q_power_m = 1.0;
   for (int m = 1 ; m <= n ; m++) {
      new->coeff[m] = cst * q_power_m;
      q_power_m *= (1-p);
   }

   return new;

}


kpoly_t *
new_kpoly_tilde_B
(
   size_t k,    // Max degree.
   size_t n,    // Degree of A.
   double p,    // Probability of error.
   double mu,   // Divergence rate.
   uint   N     // Number of duplicates.
)
{

   if (n > k || n == 0) critical_error("parameter error", __LINE__);

   kpoly_t *new = new_zero_kpoly(k);

   double cst = p * pow(1-mu/3.0, N);
   double q_power_m = 1.0;
   double one_minus_mu_power_m = 1.0;
   for (int m = 1 ; m <= n ; m++) {
      double Cm = 1 - pow(1-one_minus_mu_power_m, N);
      new->coeff[m] = cst * Cm * q_power_m;
      one_minus_mu_power_m *= (1-mu);
      q_power_m *= (1-p);
   }

   return new;

}


kpoly_t *
new_kpoly_r
(
   size_t k,    // Max degree.
   size_t n,    // Degree of A.
   double p,    // Probability of error.
   double mu,   // Divergence rate.
   uint   N     // Number of duplicates.
)
{

   if (n > k || n == 0) critical_error("parameter error", __LINE__);

   kpoly_t *new = new_zero_kpoly(k);

   double thisC = pow(1-pow(1-mu,n), N);
   double thatC = pow(1-pow(1-mu,n-1), N);
   new->mono.deg = n;
   new->mono.coeff = (thisC - thatC) * pow(1-p,n);
   new->coeff[n] = (thisC - thatC) * pow(1-p,n);

   return new;

}


kpoly_t *
new_kpoly_F
(
   size_t k,    // Max degree.
   size_t n,    // Degree of A.
   double p,    // Probability of error.
   double mu,   // Divergence rate.
   uint   N     // Number of duplicates.
)
{

   // Note that 'n' can be 0 for this polynomial.
   if (n > k) critical_error("parameter error", __LINE__);

   kpoly_t *new = new_zero_kpoly(k);

   double q_power_m = 1.0;
   for (int m = 0 ; m <= n ; m++) {
      new->coeff[m] = q_power_m;
      q_power_m *= (1-p);
   }

   return new;

}


kpoly_t *
new_kpoly_tilde_F
(
   size_t k,    // Max degree.
   size_t n,    // Degree of A.
   double p,    // Probability of error.
   double mu,   // Divergence rate.
   uint   N     // Number of duplicates.
)
{

   if (n > k || n == 0) critical_error("parameter error", __LINE__);

   kpoly_t *new = new_zero_kpoly(k);

   double q_power_m = 1.0;
   double one_minus_mu_power_m = 1.0;
   for (int m = 0 ; m <= n ; m++) {
      double Cm = 1 - pow(1-one_minus_mu_power_m, N);
      new->coeff[m] = Cm * q_power_m;
      one_minus_mu_power_m *= (1-mu);
      q_power_m *= (1-p);
   }

   return new;

}


mat_t *
new_matrix_M
(
   size_t k,
   size_t gamma,
   double p,
   double mu,
   uint   N
)
{

   size_t dim = gamma+3;
   mat_t *M = new_null_matrix(k, dim);

   // First row.
   M->term[0*dim+2] = new_zero_kpoly(k);
   M->term[0*dim+2]->coeff[0] = 1.0;

   // Second row.
   M->term[1*dim+1] = new_kpoly_tilde_A(k, k, p, mu, N);
   M->term[1*dim+2] = new_kpoly_tilde_B(k, k, p, mu, N);
   for (int i = 3 ; i < dim-1 ; i++)
      M->term[1*dim+i] = new_kpoly_r(k, i-2, p, mu, N);
   M->term[1*dim+dim-1] = new_kpoly_tilde_F(k, k, p, mu, N);

   // Third row.
   M->term[2*dim+1] = new_kpoly_tilde_A(k, k, p, mu, N);
   M->term[2*dim+2] = new_kpoly_tilde_B(k, gamma, p, mu, N);
   for (int i = 3 ; i < dim-1 ; i++)
      M->term[2*dim+i] = new_kpoly_r(k, i-2, p, mu, N);
   M->term[2*dim+dim-1] = new_kpoly_tilde_F(k, gamma-1, p, mu, N);

   // Middle rows.
   for (int j = 3 ; j < dim-1 ; j++) {
      M->term[j*dim+1] = new_kpoly_A(k, gamma+2-j, p, mu, N);
      M->term[j*dim+2] = new_kpoly_B(k, gamma+2-j, p, mu, N);
      M->term[j*dim+dim-1] = new_kpoly_F(k, gamma+1-j, p, mu, N);
   }

   // Last row is null.

   return M;

}



void do_shortcut (kpoly_t *dest,const  kpoly_t *b, int deg, double coeff) {
   for (int i = deg ; i <= dest->k ; i++) {
      dest->coeff[i] = coeff * b->coeff[i-deg];
   }   
}



kpoly_t *
kpoly_mult
(
         kpoly_t * dest,
   const kpoly_t * a,
   const kpoly_t * b
)
{

   // If any of the two k-polynomials is zero,
   // set 'dest' to zero and return 'NULL'.
   if (a == NULL || b == NULL) {
      bzero(dest->coeff, (dest->k+1) * sizeof(double));
      return NULL;
   }

   // All k-polynomials must have the same 'k'.
   if (a->k != dest->k || b->k != dest->k)
      critical_error("incongruent k-polynomials", __LINE__);

   if (a->mono.deg) {
      // If a is a monomial, use a shortcut.
      bzero(dest->coeff, (dest->k+1) * sizeof(double));
		do_shortcut(dest, b, a->mono.deg, a->mono.coeff);
//      for (int i = a->mono.deg ; i <= dest->k ; i++) {
//         dest->coeff[i] = a->mono.coeff * b->coeff[i-a->mono.deg];
//      }
   }
   else {
      // Standard convolution product.
      for (int i = 0 ; i <= dest->k ; i++) {
         dest->coeff[i] = a->coeff[0] * b->coeff[i];
         for (int j = 1 ; j <= i ; j++) {
            dest->coeff[i] += a->coeff[j] * b->coeff[i-j];
         }
      }
   }

/*
   for (int i = 0 ; i <= dest->k ; i++) {
      dest->coeff[i] = a->coeff[0] * b->coeff[i];
      for (int j = 1 ; j <= i ; j++) {
         dest->coeff[i] += a->coeff[j] * b->coeff[i-j];
      }
   }
*/

   return dest;

}


void
kpoly_update_add
(
         kpoly_t * dest,
   const kpoly_t * a
)
{

   // No update needed is second k-polynomial is zero..
   if (a == NULL) return;

   // All k-polynomials must have the same 'k'.
   if (a->k != dest->k)
      critical_error("incongruent k-polynomials", __LINE__);

   for (int i = 0 ; i <= dest->k ; i++) {
      dest->coeff[i] += a->coeff[i];
   }

}


mat_t *
matrix_mult
(
         mat_t * dest,
   const mat_t * a,
   const mat_t * b
)
{

   if (a->dim != dest->dim || b->dim != dest->dim)
      critical_error("incongruent matrices", __LINE__);

   size_t dim = dest->dim;

   kpoly_t *temp = new_zero_kpoly(dest->k);

   for (int i = 0 ; i < dim ; i++) {
   for (int j = 0 ; j < dim ; j++) {
      // Initialize destination to zero.
      if (dest->term[i*dim+j] == NULL)
         dest->term[i*dim+j] = new_zero_kpoly(dest->k);
      else
         bzero(dest->term[i*dim+j]->coeff, (dest->k+1) * sizeof(double));

      // Matrix multiplication.
      for (int k = 0 ; k < dim ; k++) {
         kpoly_update_add(dest->term[i*dim+j],
               kpoly_mult(temp, a->term[i*dim+k], b->term[k*dim+j]));
      }
   }
   }

   free(temp);

   return dest;

}


int main(void) {

   uint gamma = 17;

   uint N = 5;

   kpoly_t *w = new_zero_kpoly(100);
   mat_t *M = new_matrix_M(100, gamma, .01, .05, N);

   mat_t *powM = new_null_matrix(100, gamma+3);
   matrix_mult(powM, M, M);

   kpoly_update_add(w, powM->term[gamma+2]);
   // print_kpoly(w);

   for (int i = 0 ; i < 16 ; i++) {
      mat_t *tmpM = new_null_matrix(100, gamma+3); // FIXME //
      matrix_mult(tmpM, M, powM);
      kpoly_update_add(w, tmpM->term[gamma+2]);
      destroy_mat(powM);
      powM = tmpM;
   }

   destroy_mat(powM);
   destroy_mat(M);

   print_kpoly(w);

   free(w);

}
