#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
/*
These are include files that should exist in your C library.
*/

/*  ============	basic choice of PARAMETERS	      ============  */

#define	                Long            long
#define                 LLong           long long
/*
For reflexive polytopes in 4 or less dimensions, everything should work with
Long set to 32-bit-integer and LLong set to 64 bits.
Many applications will even work with LLong at 32 bits.
For higher dimensional or complicated non-reflexive polytopes it may be
necessary to set even Long to 64 bits.
*/

#ifndef POLY_Dmax       /* You can set POLY_Dmax at compilation:
                           use -D POLY_Dmax=<value> in the C flags  */
#define   		POLY_Dmax	6	/* max dim of polytope	    */
#endif
/*
POLY_Dmax should be set to the dimension of the polytopes that are analysed.
While the programs still work if POLY_Dmax is set to a higher value, they may
be considerably slowed down.
*/

#if	(POLY_Dmax <= 3)
#define   		POINT_Nmax	40	/* max number of points	    */
#define   		VERT_Nmax	16	/* max number of vertices   */
#define   		FACE_Nmax	30      /* max number of faces      */
#define	                SYM_Nmax	88	/* cube: 2^D*D! plus extra  */

#elif	(POLY_Dmax == 4)
#define   		POINT_Nmax	700	/* max number of points	    */
#define 		VERT_Nmax	64      /* max number of vertices   */
#define   		FACE_Nmax	824	/* max number of faces      */
#define	                SYM_Nmax	1200

#else
#define   		POINT_Nmax	2000000 
#define   		VERT_Nmax	64	/* !! use optimal value !!  */
#define   		FACE_Nmax	10000	/* max number of faces      */
#define                 SYM_Nmax        46080   /* symmetry (P_1)^6: 2^6*6! */
#define			EQUA_Nmax	1280    /* up to 20000 without alloc */
#endif

#ifndef			EQUA_Nmax			/* default setting */
#define                 EQUA_Nmax       VERT_Nmax
#endif
/*
POINT_Nmax, VERT_Nmax and FACE_Nmax denote the maximal numbers of points,
vertices and faces, respectively.
SYM_Nmax is the maximal number of symmetries of a polytope, i.e. the order of
the finite subgroup S of the group GL(n,Z) of lattice automorphisms that leaves
a polytope invariant.
EQUA_Nmax denotes the maximal number of facets (given by equations) of a
polytope. By duality this is just the number of vertices of the dual polytope,
so it makes sense to have the default setting EQUA_Nmax = VERT_Nmax.
In applications not related to reflexive polytopes or in large dimensions a
larger value may be useful. While CPU-time is almost independent of EQUA_Nmax,
it strongly depends on VERT_Nmax/32 (i.e. use 32, 64, 96, ...).
Our settings for dimensions less than or equal to 4 are such that they work
for any reflexive polytope.
*/

#define                 AMBI_Dmax       (5 * POLY_Dmax)	/* default setting */
/*
If a polytope is determined by a combined weight system it is first realised
by an embeddeding in an ambient space of dimension (Poly-dim + number of
weight systems). AMBI_Dmax is the maximal dimension of this ambient space.
*/


#define                 FIB_Nmax	3000 /*NOW: 27/5/11 default setting*/
/*
Given a polytope P* it is possible to analyze the IP simplices among its 
points. These simplices are given in terms of weight relations among points 
of P*. FIB_Nmax is the maximal number of allowed relations.
*/


#define  CD2F_Nmax               FACE_Nmax
/*
Max number of codimension 2 faces.
*/


#define GL_Long		Long
/*
Uses W_to_GLZ like in Rat.c
*/

#define MAXLD (26)

/*
Used in the handling of large lists of weight systems (cf. C5stats)
*/

extern FILE *inFILE, *outFILE;
/*
Ascii-files for input and output. If not given in the parameter list they
default to stdin and stdout, respectively.
*/


/*  ==========         Global typedefs           		==========  */

typedef struct {int n, np; Long x[POINT_Nmax][POLY_Dmax];}   PolyPointList;
/*
A list (not necessarily complete) of lattice points of a polytope.
P.x[i][j] is the j'th coordinate of the i'th lattice point.
P.n is the dimension of the polytope and P.np the number of points in the list.
*/

typedef struct {int v[VERT_Nmax]; int nv;}                   VertexNumList;
/*
The list of vertices of a polytope, referring to some PolyPointList P.
The j'th coordinate of the i'th vertex is then given by P.x[V.v[i]][j].
V.nv is the number of vertices of P.
*/

typedef struct {Long a[POLY_Dmax], c;}                       Equation;
/*
This structure determines an equation of the type ax+c=0, explicitly:
sum_{i=1}^n E.a[i] x_i + E.c = 0.
*/

typedef struct {int ne; Equation e[EQUA_Nmax];}		     EqList;
/*
A list of equations; EL.ne is the number of equations in the list.
*/

typedef struct {EqList B; Long W[AMBI_Dmax][AMBI_Dmax], d[AMBI_Dmax]; 
  int nw, N, z[POLY_Dmax][AMBI_Dmax], m[POLY_Dmax], nz, index;}       CWS;
/*
Combined weight system: W[i][j] and d[i] are the j'th weight and the "degree"
of the i'th weight system, respectively; nw is the number of weight systems,
N is the dimension of the ambient space.
z[i][j]/m[i] are the phases of nz symmetries for polytopes on sublattices.
B describes the ambient space coordinate hyperplanes in terms of the new
(non-redundant) coordinates.
*/

typedef Long PairMat[EQUA_Nmax][VERT_Nmax];
/*
The matrix whose entries are the pairings av+c between the vertices v and
the equations (a,c).
*/

typedef struct {int mp, mv, np, nv, n, pic, cor, h22, h1[POLY_Dmax-1];}
                                                                     BaHo;
/*
This structure is related to Batyrev's formulas for Hodge numbers.
n     ... dimension of the polytope
pic   ... Picard number
cor   ... sum of correction terms
h1[i] ... Hodge number h_{1i}
h22   ... Hodge number h_{22} (if n = 5)
mp, mv, np, nv denote the numbers of points/vertices in the M and N lattices,
repectively.
*/

typedef struct {
  Long W[FIB_Nmax][VERT_Nmax]; 
  int nw, PS, ZS, nv, f[VERT_Nmax],r[VERT_Nmax],nf,nz[FIB_Nmax], n0[FIB_Nmax],
    Z[FIB_Nmax][VERT_Nmax], M[FIB_Nmax];
  GL_Long G[VERT_Nmax][POLY_Dmax][POLY_Dmax];
  PolyPointList *P;
}                                                                     FibW;
/*
This list is an extension of the PolyPointList with the combined weight system.
W[i][j] is the j'th weight; nw is the number of weight systems.
*/

typedef struct{
  long n_nonIP, n_IP_nonRef, n_ref, // numbers of WS of certain types
    max_w, nr_max_w, //maximum weight in the reflexive/non-reflexive cases
    nr_n_w[MAXLD], n_w[MAXLD]; //numbers of weights of given [ld]
  int nr_max_mp, nr_max_mv, nr_max_nv, max_mp, max_mv, max_np, max_nv,
    max_h22, max_h1[POLY_Dmax-1], //max values of certain entries of BH
    min_chi, max_chi, max_nf[POLY_Dmax+1]; //range for chi, max facet numbers
}                                                                     C5stats;
/* 
statistics on large lists of weight systems, cf. classification of 4fold weights
*/

/*  ==========         I/O functions (from Coord.c)		==========  */

int  Read_CWS_PP(CWS *C, PolyPointList *P);
/*
Reads either a CWS or a PolyPointList.
If *C is read, the PolyPointList *P determined by *C is calculated, otherwise
C->nw is set to 0 to indicate that no weight has been read.
CWS-input consists of a single line of the form
d1 w11 w12 ... d2 w21 w22 ...,
whereas PolyPointList-input begins with a line
#columns #lines
followed by #lines further lines. It reads P->x as a matrix such that
either P->n = #columns and P->np = #lines or vice versa (the result is
unique because of P->np > P->n).
*/

int  Read_CWS(CWS *_CW, PolyPointList *_P);
/*
 Reads CWS input *C, the PolyPointList *P determined by *C is calculated.
*/

int  Read_PP(PolyPointList *_P);
/*
Reads the PolyPointList input *P
*/

void Print_PPL(PolyPointList *P, const char *comment);
void Print_VL(PolyPointList *P, VertexNumList *V, const char *comment);
void Print_EL(EqList *EL, int *n, int suppress_c, const char *comment);
void Print_Matrix(Long Matrix[][VERT_Nmax], int n_lines, int n_columns,
		  const char *comment);
/*
Each of these routines prints a matrix in the format
#columns #lines  *comment
line_0
...
line_{#lines-1}.
With Print_PPL and Print_VL, points/vertices are displayed as column vectors
if there's enough space and as row vectors otherwise.
Print_EL always displays equations in line format.
If *suppress_c is zero line_i reads
EL->e[i].a[0] ... EL->e[i].a[*n-1]  EL->e[i].c,
otherwise the last entry EL->e[i].c is suppressed so that the
resulting output can be used as input for Read_CWS_PP.
*/

void Print_CWH(CWS *C, BaHo *BH);
/*
Writes a single line that reproduces *C (if C->nw isn't 0, i.e. if the
input was of CWS type), information on the numbers of points and
vertices of the polytope and its dual, and the Hodge numbers of the
corresponding Calabi-Yau variety.
*C is reproduced in the format
d1 w11 w12 ... d2 w21 w22 ...
Information on the polytope is given in the form
M:mp mv N:np nv
for reflexive polytopes.
Here mp and mv are the numbers of lattice points and vertices of the
polytope, respectively, and np and nv denote the corresponding numbers
for the dual polytope.
If a polytope is not reflexive, "N:np nv" is replaced by "F:ne" (the
number of facets/equations).
Hodge number information is given in the format
H: h11 h12 ... h1(n-2) [chi],
where the h1i are the corresponding Hodge numbers and chi is the Euler
number. This output is suppressed for polytopes that are not reflexive.
As an example, the complete output for the quintic threefold reads
5 1 1 1 1 1 M:126 5 N:6 5 H:1,101 [-200].
*/

void Initialize_C5S(C5stats *_C5S, int n);
void Update_C5S(BaHo *_BH, int *nf, Long *W, C5stats *_C5S);
void Print_C5S(C5stats *_C5S);
/*
Routines for handling the structure C5stats
*/

/*  ==========              From Polynf.c                        ========== */

int  Make_Poly_Sym_NF(PolyPointList *P, VertexNumList *VNL, EqList *EL,
		      int *SymNum, int V_perm[][VERT_Nmax],
		      Long NF[POLY_Dmax][VERT_Nmax], int t, int S, int N);
/*
Given *P, *VNL and *EL, the following objects are determined:
the number *SymNum of GL(n,Z)-symmetries of the polytope,
the *SymNum vertex permutations V_perm realising these symmetries,
the normal form coordinates NF of the vertices,
the number of symmetries of the vertex pairing matrix
    (this number is the return value of Make_Poly_Sym_NF).
If t/S/N are non-zero, the output of the corresponding options of poly.x
is displayed.
*/

void IP_Simplex_Decomp(Long CM[][POLY_Dmax], int p, int d,
        int *nw, Long W[][VERT_Nmax], int Wmax, int codim);
/*
Given the matrix CM of coordinates of p points in Z^d, the list W[i] of *nw
weight systems corresponding to IP-simplices spanned by the points in CM is
created.
If codim!=0 only the IP-simplices with dimension > 1 and codimension
between 1 and codim are computed.
It is assumed that p<=VERT_Nmax and that W can hold at least Wmax sets of
coefficients.
*/

void IP_Simplices(PolyPointList *P, int nv, int PS, int VS, int CD);
/*
Realizes the -P,-V,-Z, and fibration options of poly (the results of this
routine are displayed as output; *P is not modified).
*/

int  Sublattice_Basis(int d, int p, Long *P[],     /* return index=det(D) */
	Long Z[][VERT_Nmax], Long *M, int *r, Long G[][POLY_Dmax], Long *D);
/*
Given a vector P[] of pointers at p points in N=Z^d that generate a 
(sub)lattice N' of the same dimension d, the following data are determined:
D[i] with 0 <= i < d  such that the lattice quotient N/N' is the product of 
cyclic groups Z_{D[i]} with D[i] dividing D[i+1], and a GL(d,Z) matrix G 
corresponding to a base change P->GxP such that the i'th new coordinate of 
each of the lattice points is divisible by D[i].
If p<=VERT_Nmax the program also computes *r coefficient vectors Z[i] for 
linear combinations of the points on P that are M[i]-fold multiples of
primitive lattice vectors, where M[i]=D[d-i] for i<*r.
If p>VERT_Nmax it is asserted that the index of the lattice quotient is 1.
*/

void Make_Poly_UTriang(PolyPointList *P);
/*
A coordinate change is performed that makes the matrix P->x upper triangular,
with minimal entries above the diagonal.
*/

void Make_ANF(PolyPointList *P, VertexNumList *V, EqList*E, 
	      Long ANF[][VERT_Nmax]);
/*
Given *P, *V and *E, the affine normal form ANF (i.e., a normal form
that also works for non-reflexive polytopes), is computed.
*/

int SimpUnimod(PolyPointList *P, VertexNumList *V, EqList *E, int vol);
/*
If vol is 0, the return value is 1 if all facets are simplicial, 0 otherwise
If vol is not 0, the return value is 1 if all facets are unimoular
(i.e. of volume 1) and 0 otherwise.
*/

int ConifoldSing(PolyPointList *P, VertexNumList *V, EqList *E,
		 PolyPointList *dP, EqList *dE, int CYorFANO);
/*
Realizes the -C1 or -C2 options of poly for CYorFANO being 1 or 2, respectively.
*/

int  Fano5d(PolyPointList *, VertexNumList *, EqList *);
/*
Realizes the -U5 option of poly.
*/

void Einstein_Metric(CWS *CW,PolyPointList *P,VertexNumList *V,EqList *E);
/*
Realizes the -E option of poly.
*/

int  Divisibility_Index(PolyPointList *P, VertexNumList *V);
/*
Returns the largest integer g for which *P is a g-fold multiple of some
other polytope.
*/

Long LatVol_Barycent(PolyPointList *P, VertexNumList *V, Long *B, Long *N);
/*
Given *P and *V, the coordinates of the barycenter of *P are computed (with 
the i'th coordinate as B[i] / *N) and the lattice volume of *P is returned.
*/

void IPs_degD(PolyPointList *P, VertexNumList *V, EqList *E, int l);
/*
 *P is interpreted as the origin and the first level of a Gorenstein cone. 
The points of the cone up to level l are computed and displayed together with 
information on the type of face of the cone they represent (option -B# of poly).
*/

void Make_Facet(PolyPointList *P, VertexNumList *V, EqList *E, int e, 
		Long vertices_of_facet[POLY_Dmax][VERT_Nmax], int *nv_of_facet);
/*
The e'th facet of *P is determined as a (P->n-1)-dimensional polytope:
*nv_of_facet vertices represented by vertices_of_facet.
*/

/*  ==========     General purpose functions from Vertex.c   	==========  */

void swap(int *i,int *j);
/*
Swaps *i and *j.
*/

void Sort_VL(VertexNumList *V);
/*
Sorts the entries _V->v[i] in ascending order.
*/

Long Eval_Eq_on_V(Equation *E, Long *V, int n);
/*
Evaluates E on V, i.e. calculates \sum_{i=0}^{n-1} E->a[i] * V[i] + E->c.
*/

int  Span_Check(EqList *EL, EqList *HL, int *n);
/*
Returns 1 if every equation of *HL is contained in *EL and 0 otherwise.
*n is the dimension.
*/

int  Vec_Greater_Than(Long *X, Long *Y, int n);
/*
Returns 1 if *X > *Y in the sense that X[i] > Y[i] for the first i where
X[i] and Y[i] differ, returns 0 if *X < *Y and gives an error message if
X[i] equals Y[i] for all i in {0,...n-1}.
*/

int Vec_is_zero(Long *X, int n);
/*
Returns 1 if X[i]==0 for 0<=i<n; returns 0 otherwise.
*/

void Swap_Vecs(Long *X, Long *Y, int n);
/*
Exchanges the n-dimensional vectors X and Y.
*/

Equation EEV_To_Equation(Equation *E1, Equation *E2, Long *V, int n);
/*
Returns the equation describing the span of the vector V and the intersection
of the hyperplanes corresponding to E1 and E2; n is the dimension.
*/

void Make_VEPM(PolyPointList *P, VertexNumList *VNL, EqList *EL, PairMat PM);
/*
Calculates the matrix of pairings between the vertices in VNL and the
equations in EL.
*/

int EL_to_PPL(EqList *EL, PolyPointList *DP, int *n);
/*
Converts *EL to the incomplete PolyPointList *DP corresponding to the dual
polytope; *n is the dimension. Returns 1 if all equations of *EL are at
distance 1 from the origin and 0 otherwise.
*/

int VNL_to_DEL(PolyPointList *P, VertexNumList *V, EqList *DE);
/*
Converts *V, which refers to *P, into the list *DE of equations of the
dual polytope (assuming reflexivity).
Returns 0 if _V->nv exceeds EQUA_Nmax and 1 otherwise.
*/

int Transpose_PM(PairMat PM, PairMat DPM, int nv, int ne);
/*
Transposes PM into DPM; returns 1 if the dimensions nv, ne are within their
limits and 0 otherwise.
*/


/*  ==========   Polytope analysis functions (from Vertex.c)    ==========  */

int  Find_Equations(PolyPointList *P, VertexNumList *VNL, EqList *EL);
/*
For the polytope determined by P, *VNL and *EL are calculated.
*VNL is the complete list of vertices of P.
*EL is the complete list of equations determining the facets of P.
Find_Equations returns 1 if P has IP property (i.e., it has the
origin in its interior) and 0 otherwise.
*/

int  IP_Check(PolyPointList *P, VertexNumList *VNL, EqList *EL);
/*
Same as Find_Equations, but returns immediately without
calculating *VNL and *EL if P does not have the IP property.
*/

int  Ref_Check(PolyPointList *P, VertexNumList *VNL, EqList *EL);
/*
Returns 1 if P is reflexive and 0 otherwise.
Only in the reflexive case *VNL and *EL are calculated.
*/

void Make_Dual_Poly(PolyPointList *P, VertexNumList *VNL, EqList *EL,
		    PolyPointList *DP);
/*
Given P, VNL and EL for a reflexive polytope, the complete list *DP
of lattice points of the dual polytope is determined.
*/

void Complete_Poly(Long VPM[][VERT_Nmax],EqList *E,int nv,PolyPointList *P);
/*
Given the vertex pairing matrix VPM, the EqList *E and the number nv of
vertices, the complete list of lattice points *P is determined.
*/

void RC_Calc_BaHo(PolyPointList *P, VertexNumList *VNL, EqList *EL,
		  PolyPointList *DP, BaHo *BH);
/*
Given *P, *VNL, *EL and *DP (points of dual polytope) as input, the elements
of *BH are calculated. *P must be reflexive; *P and *DP must be complete.
*/


/*  ======  typedefs and functions (from Vertex.c) related to INCIs  ====  */

#define                 INT_Nbits            32
#define                 LONG_LONG_Nbits      64
/*
These numbers should be set to the actual numbers of bits occupied by the
structures "unsigned int" and "unsigned long long" in your version of C.
If they are set to lower values, everything still works but may be
considerably slowed down.
*/

#if (VERT_Nmax <= INT_Nbits)
typedef		        unsigned int            INCI;
#elif (VERT_Nmax <= LONG_LONG_Nbits)
typedef		        unsigned long long	INCI;
#else
#define I_NUI     ((VERT_Nmax-1)/INT_Nbits+1)
typedef struct {unsigned int ui[I_NUI];}   INCI;
#endif
/*
An INCI encodes the incidence relations between a face and a list of
vertices as a bit pattern (1 if a vertex lies on the face, 0 otherwise).
Depending on the allowed number VERT_Nmax of vertices, a single "unsigned int"
or "unsigned long long" may be sufficient.
If VERT_Nmax is larger than the number of bits in a "long long integer", an
array of unsigned integers is used to simulate an integer type of the required
size.
*/

typedef struct {int nf[POLY_Dmax+1];			  /* #(faces)[dim]  */
 	INCI v[POLY_Dmax+1][FACE_Nmax]; 		  /*  vertex info   */
 	INCI f[POLY_Dmax+1][FACE_Nmax]; 		  /* V-on-dual info */
 	Long nip[POLY_Dmax+1][FACE_Nmax];		   /* #IPs on face  */
 	Long dip[POLY_Dmax+1][FACE_Nmax];} 	FaceInfo;  /* #IPs on dual  */
/*
nf[i] denotes the number of faces of dimension i
   (the number of faces of dimension n-i-1 of the dual polytope).
v[i][j] encodes the incidence relation of the j'th dim-i face with the vertices
nip[i][j] is the number of interior points of the j'th dim-i face.
f[i][j] and dip[i][j] give the same informations for the dual (n-i-1
   dimensional) faces, with f[i][j] referring to the dual vertices.
*/

#if (VERT_Nmax <= LONG_LONG_Nbits)
#define INCI_M2(x)     ((x) % 2)              /* value of first bit      */
#define	INCI_AND(x,y)  ((x) & (y))            /* bitwise logical and     */
#define	INCI_OR(x,y)   ((x) | (y))            /* bitwise logical or      */
#define	INCI_XOR(x,y)  ((x) ^ (y))            /* bitwise exclusive or    */
#define	INCI_EQ(x,y)   ((x) == (y))           /* check on equality       */
#define INCI_LE(x,y)   INCI_EQ(INCI_OR(x,y),y)/* bitwise less or equal */
#define INCI_EQ_0(x)   INCI_EQ(x,INCI_0())    /* check if all bits = 0   */
#define INCI_0()       (0)                    /* set all bits to 0       */
#define INCI_1()       (1)                    /* set only first bit to 1 */
#define INCI_D2(x)     ((x) / 2)              /* shift by one bit        */
#define INCI_PN(x,y)   (2 * (x) + !(y))       /* shift and set first bit */
/*
For an INCI defined as a single unsigned (long long) integer whose bits are
regarded as representing incidences, these are useful definitions.
INCI_PN is particularly useful when a new vertex is added: if x represents
an equation E w.r.t. some vertex list and y is the result of evaluating E
on some new vertex V, then INCI_PN(x,y) represents x w.r.t. the vertex list
enhanced by V.
*/

#else
#define INCI_M2(x)      ((x).ui[0] % 2)
INCI INCI_AND(INCI x, INCI y);
INCI INCI_OR(INCI x, INCI y);
INCI INCI_XOR(INCI x, INCI y);
int  INCI_EQ(INCI x, INCI y);
int  INCI_LE(INCI x, INCI y);
int  INCI_EQ_0(INCI x);
INCI INCI_0();
INCI INCI_1();
INCI INCI_D2(INCI x);
INCI INCI_PN(INCI x, Long y);
#endif
/*
If we need more bits than can be represented by a single unsigned long long,
these routines are designed to simulate the above definitions.
*/

int  INCI_abs(INCI X);
/*
Returns the number of bits of X whose value is 1.
*/

int  Print_INCI(INCI X);
/*
Prints X as a pattern of 0's and 1's, omitting the 0's after the last 1.
*/

INCI Eq_To_INCI(Equation *E, PolyPointList *P, VertexNumList *VNL);
/*
Converts *E to an INCI.
*/

void Make_Incidence(PolyPointList *P, VertexNumList *VNL, EqList *EL,
                    FaceInfo *FI);
/*
Creates the structure FaceInfo *FI from *P, *VNL and *EL.
*/

void Print_FaceInfo(int n, FaceInfo *FI);
/*
Displays the information contained in the FaceInfo *FI.
*/

int QuickAnalysis(PolyPointList *_P, BaHo *_BH, FaceInfo *_FI);
/*
Fast computation of FaceInfo and Hodge numbers.
*/
