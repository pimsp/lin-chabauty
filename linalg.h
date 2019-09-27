#include<pari/pari.h>

GEN GetFq1(GEN T); // Returns 1 as an element of FQ
GEN Z2Fq(GEN n, GEN T); // The map from Z to F_q
long ZX_is0mod(GEN x,GEN p); // Check if an element of Z or Z[X] is zero modulo p
GEN FpXM_red(GEN,GEN); // Reduce a matrix modulo p
GEN FpXM_add(GEN,GEN,GEN); // Add two matrices using FpX_add on every coefficient. Coefficients may be polynomials
GEN FpXM_sub(GEN,GEN,GEN); // Same but with FpX_sub
GEN FqV_Fq_mul(GEN,GEN,GEN,GEN); // Multiplies a vector in Fq with an element of Fq
GEN FqM_Fq_mul(GEN,GEN,GEN,GEN); // Same, but a matrix with a scalar
GEN ZXM_Z_mul(GEN,GEN); // Same, but a matrix with coeffs in Z[X] with an element of Z
GEN RandVec_1(GEN A,GEN pe); // Get a random linear combination of the columns of A, with all coefficients being 1 or -1 (and the first one being 1)

// The functions with padic in their name work modulo p^e, where e is also a parameter)

GEN RandVec_padic(GEN,GEN,GEN,GEN); // Get a random linear combination of the columns of A
GEN matkerpadic(GEN,GEN,GEN,long); // Get a matrix representing the kernel of the matrix A
GEN matkerpadic_hint(GEN,GEN,GEN,long,GEN,ulong); // The same, but this time knowing the dimension of the kernel 
GEN mateqnpadic(GEN,GEN,GEN,long); // The `equation matrix`: take the transpose, then the kernel, then the transpose
GEN matimagepadic(GEN,GEN,GEN,long); // The image
GEN matF(GEN,GEN,GEN,long); // Unclear
GEN mat2col(GEN); // Go from a mxn matrix to a column vector of length mn
GEN col2mat(GEN,long,long); // Go from a column vector of length mn to a mxn matrix
GEN M2ABCD(GEN,GEN); // Unclear
GEN M2ABCD_1block(GEN,ulong,ulong,GEN); // Unclear
GEN VecSmallCompl(GEN,ulong); // Unclear
GEN FqM_MinorCompl(GEN,GEN,GEN); // Unclear
