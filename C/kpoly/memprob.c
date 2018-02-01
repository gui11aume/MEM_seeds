#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>


//  TYPE DECLARATIONS  //

typedef unsigned int  uint;   // Just a shortcut.

// Struct declarations and definitions //

typedef struct kpoly_t  kpoly_t;
typedef struct mat_t    mat_t;
typedef struct mono_t   mono_t;

struct mono_t {
   size_t deg;
   double coeff;
};

struct kpoly_t {
   mono_t mono;              // Monomer (if applicable).
   double coeff[];           // Terms of the polynomial.
};

struct mat_t {
   const size_t    dim;      // Column / row number.
         kpoly_t * term[];   // Terms of the matrix.
};


// GLOBAL VIARABLES //

double    P;     // Probability of a read error.
double    U;     // Divergence rate between duplicates.
size_t    K;     // Max degree of the polynomials (read size).

size_t    KSZ;   // Size of the 'kpoly_t' struct.

kpoly_t * TEMP = NULL;        // Matrix multipliciation.
kpoly_t * ARRAY[1024] = {0};  // Store the results.

int       ERRNO; // Error codes.



// FUNCTION DEFINITIONS //

int
memprob_init
(
   double p,
   double u,
   size_t k
)
{

   P = p;
   U = u;
   K = k; 

   KSZ = sizeof(kpoly_t) + (K+1) * sizeof(double);

   for (int i = 0 ; i < 1024 ; i++) free(ARRAY[i]);
   bzero(ARRAY, 1024 * sizeof(kpoly_t *));

   free(TEMP);
   TEMP = calloc(1, KSZ);

   ERRNO = TEMP == NULL ? __LINE__ : 0;

   // Return 1 on error, 0 on success.
   return TEMP != NULL;

}



kpoly_t *
new_zero_kpoly (void)
{
   return calloc(1, KSZ);
}



kpoly_t *
new_kpoly_A
(
   size_t n,    // Degree of A.
   uint   N     // Number of duplicates.
)
{

   if (n > K || n == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();

   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      double cst = P * (1-pow(1-U/3.0, N));
      double q_power_m = 1.0;
      for (int m = 1 ; m <= n ; m++) {
         new->coeff[m] = cst * q_power_m;
         q_power_m *= (1-P);
      }
   }

   return new;

}


kpoly_t *
new_kpoly_tilde_A
(
   size_t n,    // Degree of A.
   uint   N     // Number of duplicates.
)
{

   if (n > K || n == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();

   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      double cst = P * (1-pow(1-U/3.0, N));
      double one_minus_mu_power_m = 1.0;
      double q_power_m = 1.0;
      for (int m = 1 ; m <= n ; m++) {
         double Cm = 1 - pow(1-one_minus_mu_power_m, N);
         new->coeff[m] = cst * Cm * q_power_m;
         one_minus_mu_power_m *= (1-U);
         q_power_m *= (1-P);
      }
   }

   return new;

}


kpoly_t *
new_kpoly_B
(
   size_t n,    // Degree of A.
   uint   N     // Number of duplicates.
)
{

   if (n > K || n == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();

   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      double cst = P * pow(1-U/3.0, N);
      double q_power_m = 1.0;
      for (int m = 1 ; m <= n ; m++) {
         new->coeff[m] = cst * q_power_m;
         q_power_m *= (1-P);
      }
   }

   return new;

}


kpoly_t *
new_kpoly_tilde_B
(
   size_t n,    // Degree of A.
   uint   N     // Number of duplicates.
)
{

   if (n > K || n == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();

   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      double cst = P * pow(1-U/3.0, N);
      double q_power_m = 1.0;
      double one_minus_mu_power_m = 1.0;
      for (int m = 1 ; m <= n ; m++) {
         double Cm = 1 - pow(1-one_minus_mu_power_m, N);
         new->coeff[m] = cst * Cm * q_power_m;
         one_minus_mu_power_m *= (1-U);
         q_power_m *= (1-P);
      }
   }

   return new;

}


kpoly_t *
new_kpoly_r
(
   size_t n,    // Degree of A.
   uint   N     // Number of duplicates.
)
{

   if (n > K || n == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();

   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      double thisC = pow(1-pow(1-U,n), N);
      double thatC = pow(1-pow(1-U,n-1), N);
      new->coeff[n] = (thisC - thatC) * pow(1-P,n);
      // Polynomial is a monomial.
      new->mono.deg = n;
      new->mono.coeff = (thisC - thatC) * pow(1-P,n);
   }

   return new;

}


kpoly_t *
new_kpoly_F
(
   size_t n,    // Degree of A.
   uint   N     // Number of duplicates.
)
{

   // Note that 'n' can be 0 for F polynomials.
   if (n > K) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();
   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      double q_power_m = 1.0;
      for (int m = 0 ; m <= n ; m++) {
         new->coeff[m] = q_power_m;
         q_power_m *= (1-P);
      }
   }

   return new;

}


kpoly_t *
new_kpoly_tilde_F
(
   size_t n,    // Degree of A.
   uint   N     // Number of duplicates.
)
{

   if (n > K) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();
   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      double q_power_m = 1.0;
      double one_minus_mu_power_m = 1.0;
      for (int m = 0 ; m <= n ; m++) {
         double Cm = 1 - pow(1-one_minus_mu_power_m, N);
         new->coeff[m] = Cm * q_power_m;
         one_minus_mu_power_m *= (1-U);
         q_power_m *= (1-P);
      }
   }

   return new;

}



mat_t *
new_null_matrix
(
   size_t dim
)
// Create a matrix where all kpoly are NULL.
{

   // Initialize to zero.
   mat_t *new = calloc(1, sizeof(mat_t) + dim*dim * sizeof(kpoly_t *));

   // The dimension is set upon creation
   // and must never change afterwards.
   if (new != NULL) *(size_t *)&new->dim = dim;

   return new;

}


void
destroy_mat
(
   mat_t *mat
)
{

   for (int i = 0 ; i < (mat->dim)*(mat->dim) ; i++)
      free(mat->term[i]);
   free(mat);

}



mat_t *
new_zero_matrix
(
   size_t dim
)
// Create a matrix where all terms
// are 'kpoly_t' struct set to zero.
{

   mat_t *new = new_null_matrix(dim);

   if (new == NULL) {
      ERRNO = __LINE__;
      return NULL;
   }

   for (int i = 0 ; i < dim*dim ; i++) {
      new->term[i] = new_zero_kpoly();
      if (new->term[i] == NULL) {
         ERRNO = __LINE__;
         destroy_mat(new);
         return NULL;
      }
   }

   return new;

}



mat_t *
new_matrix_M
(
   size_t gamma,
   uint   N
)
{

   size_t dim = gamma+3;
   mat_t *M = new_null_matrix(dim);

   if (M == NULL) {
      ERRNO = __LINE__;
   }
   else {
      // First row.
      M->term[0*dim+2] = new_zero_kpoly();
      M->term[0*dim+2]->coeff[0] = 1.0;

      // Second row.
      M->term[1*dim+1] = new_kpoly_tilde_A(K, N);
      M->term[1*dim+2] = new_kpoly_tilde_B(K, N);
      for (int i = 3 ; i < dim-1 ; i++)
         M->term[1*dim+i] = new_kpoly_r(i-2, N);
      M->term[1*dim+dim-1] = new_kpoly_tilde_F(K, N);

      // Third row.
      M->term[2*dim+1] = new_kpoly_tilde_A(K, N);
      M->term[2*dim+2] = new_kpoly_tilde_B(gamma, N);
      for (int i = 3 ; i < dim-1 ; i++)
         M->term[2*dim+i] = new_kpoly_r(i-2, N);
      M->term[2*dim+dim-1] = new_kpoly_tilde_F(gamma-1, N);

      // Middle rows.
      for (int j = 3 ; j < dim-1 ; j++) {
         M->term[j*dim+1] = new_kpoly_A(gamma+2-j, N);
         M->term[j*dim+2] = new_kpoly_B(gamma+2-j, N);
         M->term[j*dim+dim-1] = new_kpoly_F(gamma+1-j, N);
      }

      // Last row is null.
   }

   return M;

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
      bzero(dest, KSZ);
      return NULL;
   }

   if (b->mono.deg) {
      // If 'b' is a monomial, use a shortcut.
      bzero(dest, KSZ);
      for (int i = b->mono.deg ; i <= K ; i++) {
         dest->coeff[i] = b->mono.coeff * a->coeff[i-b->mono.deg];
      }
   }
   else {
      // Standard convolution product.
      for (int i = 0 ; i <= K ; i++) {
         dest->coeff[i] = a->coeff[0] * b->coeff[i];
         for (int j = 1 ; j <= i ; j++) {
            dest->coeff[i] += a->coeff[j] * b->coeff[i-j];
         }
      }
   }

   return dest;

}



void
kpoly_update_add
(
         kpoly_t * dest,
   const kpoly_t * a
)
{

   // No update when adding zero.
   if (a == NULL) return;

   for (int i = 0 ; i <= K ; i++) {
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

   if (a->dim != dest->dim || b->dim != dest->dim) {
      ERRNO = __LINE__;
      return NULL;
   }

   size_t dim = dest->dim;

   for (int i = 0 ; i < dim ; i++) {
   for (int j = 0 ; j < dim ; j++) {
      bzero(dest->term[i*dim+j], KSZ);
      for (int m = 0 ; m < dim ; m++) {
         kpoly_update_add(dest->term[i*dim+j],
               kpoly_mult(TEMP, a->term[i*dim+m], b->term[m*dim+j]));
      }
   }
   }

   return dest;

}

// FIXME //
void print_kpoly (kpoly_t *p) {
   if (p == NULL) {
      fprintf(stderr, "0\n");
      return;
   }
   fprintf(stderr, "%f", p->coeff[0]);
   fprintf(stderr, " + %fz", p->coeff[1]);
   for (int i = 2 ; i <= K ; i++) {
      fprintf(stderr, " + %fz^%d", p->coeff[i], i);
   }
   fprintf(stderr, "\n");
}


int main(void) {

   memprob_init(0.01, 0.05, 100);

   uint gamma = 17;
   uint N = 5;

   kpoly_t *w = new_zero_kpoly();
   mat_t *M = new_matrix_M(gamma, N);

   mat_t *powM1 = new_zero_matrix(gamma+3);
   mat_t *powM2 = new_zero_matrix(gamma+3);

   if (powM1 == NULL || powM2 == NULL) {
      ERRNO = __LINE__;
      return 1; // FIXME //
   }

   matrix_mult(powM1, M, M);
   kpoly_update_add(w, powM1->term[gamma+2]);


   for (int i = 0 ; i < 8 ; i++) {
      matrix_mult(powM2, powM1, M);
      kpoly_update_add(w, powM2->term[gamma+2]);
      matrix_mult(powM1, powM2, M);
      kpoly_update_add(w, powM1->term[gamma+2]);
   }

   destroy_mat(powM1);
   destroy_mat(powM2);
   destroy_mat(M);

   print_kpoly(w);

   free(w);

}
