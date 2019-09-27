#include<pari/pari.h>

/*
This and pic.c implement the Jacobian itself. As a Pari object, a Jacobian is
represented as an array of length lgJ, containing many of its useful parameters.

The construction seems to happen in hyper.c, more specifically in HyperInit.
*/

#define lgJ 15 // The number of elements in the Jacobian

GEN Jgetf(GEN J); // The equation for the curve of which this is the Jacobian
long Jgetg(GEN J); // The genus of the curve
long Jgetd0(GEN J); // The degree of D_0
GEN JgetT(GEN J); // Unclear
GEN Jgetp(GEN J); // The prime p we're looking at
long Jgete(GEN J); // The exponent e s.t. we're looking modulo p^e
GEN Jgetpe(GEN J); // p^e
GEN JgetFrob(GEN J); // Unclear
GEN JgetV(GEN J); // V = V_2 = L(2D_0)
GEN JgetKV(GEN J); // An `equation matrix` for V
GEN JgetW0(GEN J); // W0 = WD_0 = L(D_0), representing the trivial element in the Jacobian
GEN JgetZ(GEN J); // The random points on the curve
GEN JgetFrobCyc(GEN J); // Unclear
GEN JgetV3(GEN J); // V_3 = L(3D_0)
GEN JgetKV3(GEN J); // An `equation matrix` for V_3
void JgetTpe(GEN J, GEN* T, GEN* pe, GEN* p, long* e); // Unclear

GEN PicRed(GEN J, ulong e); // Calculata J modulo p^e

GEN ZpXQ_FrobMat(GEN T, GEN p, long e, GEN pe); // Unclar
GEN Frob(GEN x, GEN FrobMat, GEN T, GEN pe); // Unclear

GEN DivMul(GEN f, GEN W, GEN T, GEN pe); // Given W, a matrix of functions, and f, a function, calculate the matrix corresponding to W*f
GEN DivAdd0(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEN pe, ulong excess); // Given M_A, M_B rep L(A), L(B), return L(A+B). d is the dimension of the result. There are some conditions, in 2.1.2 of the paper.
GEN DivAdd1(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEN pe, ulong excess); // Unclear what the difference is with the previous
GEN DivAdd(GEN WA, GEN WB, ulong d, GEN T, GEN p, long e, GEN pe, ulong excess); // Unclear what the difference is with the previous
GEN DivSub(GEN WA, GEN WB, GEN KV, ulong d, GEN T, GEN p, long e, GEN pe, ulong nIGS); // Same, but L(A-B) instead
GEN PicChord(GEN J, GEN WA, GEN WB, long flag); // Given x,y in the Jacobian, return -(x+y). Unclear what flag does
GEN PicAdd(GEN J, GEN WA, GEN WB); // Given x,y in the Jacobian, return x+y
GEN PicSub(GEN J, GEN WA, GEN WB); // Given x,y in the Jacobian, return y-x
GEN PicNeg(GEN J, GEN W); // Given x in the Jacobian, return -x
GEN PicMul(GEN J, GEN W, GEN n, long flag); // Given x in the Jacobian, n in Z, return n*x. Unclear what flag does. 2 seems to be default
GEN PicFrob(GEN J, GEN W); // Unclear
GEN PicFrobPoly(GEN J, GEN W, GEN F); // Unclear
long PicEq(GEN J, GEN WA, GEN WB); // Given x,y in the Jacobiann return x=y
long PicIsZero(GEN J, GEN W); // Given x in the Jacobian, return x=0
GEN PicRand0(GEN J); // Generate a random element in the Jacobian

GEN PicChart(GEN J, GEN W, ulong P0); // Unclear
