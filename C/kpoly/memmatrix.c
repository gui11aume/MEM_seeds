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


// GLOBAL VARIABLES / CONSTANTS //

size_t    G;       // Minimum size of MEM seeds.
size_t    N;       // Number of duplicate sequences.
size_t    K;       // Max degree of the polynomials (read size).

size_t    HIGH;    // Proxy for infinite degree polynomials.

double    P;       // Probability of a read error.
double    U;       // Divergence rate between duplicates.
double    OMEGA;   // Probability of double-down arrow.
double   _OMEGA;   // Probability of down arrow.

size_t    KSZ;     // Size of the 'kpoly_t' struct.

kpoly_t * TEMP = NULL;        // Matrix multipliciation.
kpoly_t * ARRAY[1024] = {0};  // Store the results (indexed by N).

int       ERRNO; // Error codes.


// MACROS //
#define YES 1
#define NO  0

// Creation of a new 'kpoly_t' is just a call to 'calloc()'.
#define new_zero_kpoly() calloc(1, KSZ)

// Prob one of m altnerative threads survives i steps.
#define xi(i,m) ( 1.0 - pow( 1.0 - pow(1.0-U,(i)), (m) ))

// Calculation intermediates (one index).
#define aN(i) pow( 1.0 - pow(1.0-U,(i)) * U/3.0, N )
#define gN(i) pow( 1.0 - pow(1.0-U,(i)), N )
#define dN(i) pow( 1.0 - (1.0 - U + U*U/3.0) * pow(1.0-U,(i)), N )

// Calculation intermediates (two indices).
#define bN(j,i) pow( 1.0 - pow(1.0-U,(j))*U/3.0 + \
                           pow(1-U,(i))*(1-U/3.0), N)


void print_kpoly (kpoly_t *p);


// FUNCTION DEFINITIONS //

int
memp_init
(
   size_t g,
   size_t n,
   size_t k,
   double p,
   double u
)
// Initialize the global variables from user-defined values.
{

   G = g;  // MEM size.
   N = n;  // Number duplicates.
   K = k;  // Read size.
   P = p;  // Sequencing error.
   U = u;  // Divergence rate.

   // Set 'HIGH' to the greater of 'K' or 'G'.
   HIGH = K > G ? K : G;

   // All 'kpoly_t' must be congruent.
   KSZ = sizeof(kpoly_t) + (K+1) * sizeof(double);

   // Initialize symbolic constants.
    OMEGA = P*pow(1.0-U/3, N);
   _OMEGA = P*(1-pow(1.0-U/3, N));

   // Clean previous values (if any).
   for (int i = 0 ; i < 1024 ; i++) free(ARRAY[i]);
   bzero(ARRAY, 1024 * sizeof(kpoly_t *));

   // Allocate or reallocate 'TEMP'.
   free(TEMP);
   TEMP = calloc(1, KSZ);

   // If (re)allocation worked, clear 'ERRNO'.
   ERRNO = TEMP == NULL ? __LINE__ : 0;

   // Return 1 on error, 0 on success.
   return TEMP != NULL;

}


kpoly_t *
new_kpoly_A
(
   const size_t deg,  // Degree of polynomial D.
   const int    tilde // Return D or tilde D.
)
{

   if (deg > K || deg == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();

   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      // See definition of polynomial A.
      const double cst = tilde ? _OMEGA : OMEGA;
      double pow_of_q = 1.0;
      for (int i = 1 ; i <= deg ; i++) {
         new->coeff[i] = cst * (xi(i-1,N)) * pow_of_q;
         pow_of_q *= (1.0-P);
      }
   }

   return new; // NULL in case of failure.

}


kpoly_t *
new_kpoly_B
(
   const size_t deg,  // Degree of polynomial B.
   const int    tilde // Return D or tilde B.
)
{

   if (deg > K || deg == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();

   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      // See definition of polynomial B.
      const double cst = tilde ? _OMEGA : OMEGA;
      const double denom = 1.0 - pow(1-U/3.0, N);
      double pow_of_q = 1.0-P;
      new->coeff[1] = cst;
      for (int i = 2 ; i <= deg ; i++) {
         double numer = aN(i-1) - aN(0);
         new->coeff[i] = cst * numer / denom * pow_of_q;
         pow_of_q *= (1.0-P);
      }
   }

   return new; // NULL in case of failure.

}


kpoly_t *
new_kpoly_C
(
   const size_t deg,  // Degree of polynomial C.
   const int    tilde // Return C or tilde C.
)
{

   if (deg > K || deg == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();

   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      // See definition of polynomial C.
      const int j = G - deg;
      const double denom_j = aN(j) - aN(j-1) - gN(j) + dN(j-1);
      const double cst = tilde ? _OMEGA : OMEGA;
      double pow_of_q = 1.0;
      for (int i = 1 ; i <= deg ; i++) {
         double numer = aN(j) - aN(j-1) - bN(j,i+j-1) + bN(j-1,i+j-1);
         new->coeff[i] = cst * numer / denom_j * pow_of_q;
         pow_of_q *= (1.0-P);
      }
   }

   return new; // NULL in case of failure.

}


kpoly_t *
new_kpoly_D
(
   const size_t deg,  // Degree of polynomial D.
   const int    tilde // Return D or tilde D.
)
{

   if (deg > K || deg == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();

   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      // See definition of polynomial D.
      const double cst = tilde ? _OMEGA : OMEGA;
      double pow_of_q = 1.0;
      for (int i = 1 ; i <= deg ; i++) {
         new->coeff[i] = cst * pow_of_q;
         pow_of_q *= (1.0-P);
      }
   }

   return new; // NULL in case of failure.

}


kpoly_t *
new_kpoly_u
(
   const size_t deg  // Degree of polynomial u.
)
{

   if (deg > K || deg >= G || deg == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();

   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      // See definition of polynomial u.
      new->mono.deg = deg;
      new->mono.coeff = (xi(deg-1,N) - xi(deg,N)) * pow(1.0-P,deg);
   }

   return new; // NULL in case of failure.

}


kpoly_t *
new_kpoly_v
(
   const size_t deg  // Degree of polynomial v.
)
{

   if (deg > K || deg >= G || deg == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();

   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      // See definition of polynomial v.
      new->mono.deg = deg;
      double numer =  aN(deg) - aN(deg-1) - gN(deg) + dN(deg-1);
      double denom = 1.0 - pow(1-U/3.0, N);
      new->mono.coeff = numer / denom * pow(1.0-P, deg);
   }

   return new; // NULL in case of failure.

}


kpoly_t *
new_kpoly_w
(
   const size_t deg  // Degree of polynomial w.
)
{

   if (deg > K || deg >= G || deg == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();

   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      // See definition of polynomial w.
      new->mono.deg = deg;
      double numer = gN(deg) - dN(deg-1);
      double denom = 1.0 - pow(1-U/3.0, N);
      new->mono.coeff = numer / denom * pow(1.0-P, deg);
   }

   return new; // NULL in case of failure.

}



kpoly_t *
new_kpoly_y
(
   const int j, 
   const int i
)
{

   const size_t deg = i;

   if (deg > K || deg >= G || deg == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();

   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      // See definition of polynomial y.
      new->mono.deg = deg;
      double numer = bN(j,j+i) - bN(j,j+i-1) - bN(j-1,i+j) + bN(j-1,j+i-1);
      double denom_j = aN(j) - aN(j-1) - gN(j) + dN(j-1);
      new->mono.coeff = numer / denom_j * pow(1.0-P, deg);
   }

   return new; // NULL in case of failure.

}


kpoly_t *
new_kpoly_T_down
(
   void
)
{

   kpoly_t *new = new_zero_kpoly();
   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      double pow_of_q = 1.0;
      for (int i = 0 ; i <= K ; i++) {
         new->coeff[i] = pow_of_q;
         pow_of_q *= (1-P);
      }
   }

   return new;

}


kpoly_t *
new_kpoly_T_double_down
(
   void
)
{

   kpoly_t *new = new_zero_kpoly();
   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      double pow_of_q = 1.0;
      for (int i = 0 ; i <= G-1 ; i++) {
         new->coeff[i] = xi(i,N) * pow_of_q;
         pow_of_q *= (1-P);
      }
   }

   return new;

}


kpoly_t *
new_kpoly_T_up
(
   size_t deg
)
{

   if (deg > K || deg >= G || deg == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();
   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      double pow_of_q = 1.0;
      for (int i = 0 ; i <= deg ; i++) {
         new->coeff[i] = pow_of_q;
         pow_of_q *= (1-P);
      }
   }

   return new;

}


kpoly_t *
new_kpoly_T_sim
(
   size_t deg
)
{

   if (deg > K || deg >= G || deg == 0) {
      ERRNO = __LINE__;
      return NULL;
   }

   kpoly_t *new = new_zero_kpoly();
   if (new == NULL) {
      ERRNO = __LINE__;
   }
   else {
      const int j = G-1 - deg;
      const double denom_j = aN(j) - aN(j-1) - gN(j) + dN(j-1);
      double pow_of_q = 1.0;
      for (int i = 0 ; i <= deg ; i++) {
         double numer = aN(j) - aN(j-1) - bN(j,i+j) + bN(j-1,i+j);
         new->coeff[i] = numer / denom_j * pow_of_q;
         pow_of_q *= (1.0-P);
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

   for (int i = 0 ; i < (mat->dim)*(mat->dim) ; i++) free(mat->term[i]);
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
   void
)
{

   size_t dim = 2*G+2;
   mat_t *M = new_null_matrix(dim);

   if (M == NULL) {
      ERRNO = __LINE__;
   }
   else {
      // First row.
      M->term[0*dim+1] = new_zero_kpoly();
      M->term[0*dim+1]->coeff[0] = 1.0;

      // Second row.
      M->term[1*dim+1] = new_kpoly_A(G, NO);
      M->term[1*dim+2] = new_kpoly_A(HIGH, YES);
      for (int j = 1 ; j <= G-1 ; j++)
         M->term[1*dim+j+G+1] = new_kpoly_u(j);
      M->term[1*dim+dim-1] = new_kpoly_T_double_down();

      // Third row.
      M->term[2*dim+1] = new_kpoly_B(HIGH, NO);
      M->term[2*dim+2] = new_kpoly_B(HIGH, YES);
      for (int j = 1 ; j <= G-1 ; j++)
         M->term[2*dim+j+2] = new_kpoly_v(j);
      for (int j = 1 ; j <= G-1 ; j++)
         M->term[2*dim+j+G+1] = new_kpoly_w(j);
      M->term[1*dim+dim-1] = new_kpoly_T_down();

      // First series of middle rows.
      for (int j = 1 ; j <= G-1 ; j++) {
         M->term[(j+2)*dim+1] = new_kpoly_C(G-j, NO);
         M->term[(j+2)*dim+2] = new_kpoly_C(G-j, YES);
         for (int i = 1 ; i <= G-j-1 ; i++)
            M->term[(j+2)*dim+G+2+i] = new_kpoly_y(j,i);
         M->term[(j+2)*dim+dim-1] = new_kpoly_T_sim(G-j-1);
      }

      // Second series of middle rows.
      for (int j = 1 ; j <= G-1 ; j++) {
         M->term[(j+G+1)*dim+1] = new_kpoly_D(G-j, NO);
         M->term[(j+G+1)*dim+2] = new_kpoly_D(G-j, YES);
         M->term[j*dim+dim-1] = new_kpoly_T_up(G-j-1);
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

   if (b->mono.deg || b->mono.coeff) {
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
      fprintf(stderr, " + %.10fz^%d", p->coeff[i], i);
   }
   fprintf(stderr, "\n");
}


int main(void) {

   memp_init(17, 2, 100, 0.01, 0.05);

   kpoly_t *w = new_zero_kpoly();
   mat_t *M = new_matrix_M();

   mat_t *powM1 = new_zero_matrix(2*G+2);
   mat_t *powM2 = new_zero_matrix(2*G+2);

   if (powM1 == NULL || powM2 == NULL) {
      ERRNO = __LINE__;
      return 1; // FIXME //
   }

   matrix_mult(powM1, M, M);
   kpoly_update_add(w, powM1->term[2*G+1]);

   for (int i = 0 ; i < 2 ; i++) {
      matrix_mult(powM2, powM1, M);
      kpoly_update_add(w, powM2->term[2*G+1]);
      matrix_mult(powM1, powM2, M);
      kpoly_update_add(w, powM1->term[2*G+1]);
   }

   destroy_mat(powM1);
   destroy_mat(powM2);
   destroy_mat(M);

   print_kpoly(w);

   free(w);

}
