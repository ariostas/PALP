#pragma once

#include <vector>

#define WZinput (1) /* WZ-input (in progress)  */

#define W_Nmax (POLY_Dmax + 1)

#if (POLY_Dmax < 7)
using Pint = int;
#else
using Pint = long;
#endif /* type of coefficients in PolyCoeffList */

typedef struct {
  int n;
  std::vector<int> e;
  std::vector<Pint> c;
  int A;
} PoCoLi; /* e=exp */

void AllocPoCoLi(PoCoLi *P);                    /* allocate e[P.A] and c[P.A] */
void Free_PoCoLi(PoCoLi *P);                    /* free P.e and P.c */
void Poly_Sum(PoCoLi *A, PoCoLi *B, PoCoLi *S); /* S = A+B */
void Poly_Dif(PoCoLi *A, PoCoLi *B, PoCoLi *D); /* D = A-B */
void PolyProd(PoCoLi *A, PoCoLi *B, PoCoLi *AB);              /* AB = A*B */
int BottomUpQuot(PoCoLi *N, PoCoLi *D, PoCoLi *Q, PoCoLi *R); /* Q*D = N-R */
void PolyCopy(PoCoLi *X, PoCoLi *Y);                          /*  Y = X   */
void PrintPoCoLi(PoCoLi *P, FILE *out = stdout);
void UnitPoly(PoCoLi *P);
void Init1_xN(PoCoLi *P, int N); /* 1 - x^N */
void PoincarePoly(int N, int *w, int d, PoCoLi *PP, PoCoLi *Naux, PoCoLi *Raux,
                  FILE *out = stdout);

int IsDigit(char c);

typedef struct {
  int d, N, z[POLY_Dmax][W_Nmax], m[POLY_Dmax], M, r, R; /* Ref */
  Long w[W_Nmax], B[W_Nmax][POLY_Dmax], A[W_Nmax], rI[POLY_Dmax];
  PolyPointList *P;
} /* Eq: Ei.c=Ai Ei.a[]=Bi[]} */
/* 0<=A+B*x  r=sum(w)/d  rI=IP(r*P)  n=(r,rI) */ Weight;

typedef struct {
  int D, E, sts;
  Pint h[POLY_Dmax][POLY_Dmax];
} VaHo;

/* AmbiPointList is defined identically in Nef.h and LG.cpp.  Use Nef.h's tag.
 */
typedef struct AmbiPointList_ AmbiPointList;

/* AmbiLatticeBasis is used by nef.c; the full definition must match LG.cpp. */
typedef struct AmbiLatticeBasis_ {
  Long x[POLY_Dmax][W_Nmax];
  int N, n;
} AmbiLatticeBasis;

int Read_W_PP(Weight *, PolyPointList *, FILE *out = stdout);
int Read_Weight(Weight *_W);
void WeightLatticeBasis(Weight *_w, AmbiLatticeBasis *_B);
void WeightMakePoints(Weight *_W, AmbiPointList *_P);
int ChangeToTrianBasis(AmbiPointList *_AP, AmbiLatticeBasis *_B,
                       PolyPointList *_PP);
int Trans_Check(Weight);
void LGO_VaHo(Weight *, VaHo *, FILE *out = stdout);
void Write_Weight(Weight *_W, FILE *out = stdout);
void Write_WH(Weight *_W, BaHo *_BH, VaHo *_VH, int rc, int tc,
              PolyPointList *_P, VertexNumList *_V, EqList *_E,
              FILE *out = stdout);
void Make_Poly_Points(Weight *_W_in, PolyPointList *_PP, FILE *out = stdout);
