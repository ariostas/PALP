/* =========================================================== */
/* ===                                                     === */
/* ===                  M o r i . h                        === */
/* ===                                                     === */
/* ===	Authors: Maximilian Kreuzer, Nils-Ole Walliser	   === */
/* ===	Last update: 19/03/12                              === */
/* ===                                                     === */
/* =========================================================== */


/* ======================================================== */
/* =========            D E F I N I T I O N s     ========= */

/*** local for Moricone.c ***/
#define Inci64		unsigned long long

/* ====================================================== */
/* =========           T Y P E D E F s          ========= */

/*** from mori.c ***/
typedef struct {
  	int FilterFlag, g, m, P, K, i, t, c, d, a, b, D, H, I, M, Read_HyperSurfCounter;
} MORI_Flags;
/*
List of flags that correspond to the options of mori.x -h for
more info plus Read_HyperSurfCounter. The latter controls if an 
hypersurface class has been inserted (option -H). Avoids to type 
the hyper class each time that HyperSurfSingula() is called. 
HyperSurfSingular() is called for each triangulation of the polytope.
*/

/*** from Moricone.c ***/
typedef struct {
	int d, v, n, nmax;
	Inci64 *I;
} triang;
/*
INCIDENCE structure needed for Triangulation and SR-ideal data 
d: dimension of the lattice polytope
v: number of digits in an INCIDENCE
n: number of INCIDENCES
nmax: max number of n
I: INCIDENCE, (Inci64: unsigned long long)
*/

/* ====================================================== */
/* =========         P R O T O T Y P E s        ========= */

/*** from MoriCone.c ***/
void TriList_to_MoriList(PolyPointList *_P, FibW *F, MORI_Flags * _Flag);
/*
Having read the Polytope matrix, asks for the triangulation and computes
SR-Ideal and Mori cone of the ambient space. The user has to specify the
number of triangulations, the number of simplices and insert the
simplices in INCI format.
 */

void HyperSurfDivisorsQ(PolyPointList *_P,VertexNumList*V,EqList *E,
		MORI_Flags *Flag);
/*
Hypersurface divisors Q(charges) permutes the N-lattice points of
non-intersecting divisors to the end of the PPL *_P, calls IP_Simplex_Fiber
and then prints the charges = linear relations.
 */

void Subdivide(PolyPointList *P,int v,Inci64 I[],int p,Inci64 *T,int *t,
		MORI_Flags *_Flag, FibW *F);
/*
Triangulate and call InterSectionRing. Incidences of ni=*t bounding
equations on I[]; subdivide -> T[*t]; assuming 4 dimensions and at most 3
non-vertices, etc.; triangulate and call InterSectionRing(T,t,P,p...).
For each T[] if simplicial add <=3 pts by hand else use GKZ for <= 3d secondary
fan.
*/

void GKZsubdivide(Inci64 *F,int f,PolyPointList *P,int p,int *Tp,int *ntp,
		int nPS, MORI_Flags *_Flag, FibW *_F );
/*
List maximal 2ndary fans of facets with descending dimensions for all compatible
max triangulations make (induced) triangulation of facets. Depending in dim of
the 2ndary fan triangulates circuits calling either Triang1dSFan(), Triang2dSFan()
or Triang2dSFan(). If dim (2ndary fan) > 3 then exit(0);
*/

void InterSectionRing(Inci64 *Tri,int *t,PolyPointList *P, int p,
		MORI_Flags *_Flag, FibW *F);
/*
Print triangulation and call StanleyReisner(SR,T). Call HyperSurfaceSingular(P,T,SR...).
Call Print_Mori(P,p,t,Tri).
*/

void StanleyReisner(triang *SR,triang *T);
/*
Determine and print the SR ideal.
*/

void DivClassBasis(int SF,PolyPointList *P,int v,char *D,char *B);
/*
Find basis for divisor classes, i.e. integral basis of intersection ring.
Simply trying to find toric divisors that span the lattice of DivClasses i.e. find
subdeterminants=volumes=1 for elimination.
*/


/*  ======  typedefs and functions related to INCIDENCEs  ======  */

Inci64 makeN(int N);
void putN(int N,Inci64 *I); /* make INCIDENCE */
void setN(int N,Inci64 *I); /* make INCIDENCE */
int  getN(int N,Inci64 I);  /* read INCIDENCE */

void prnI(int N,Inci64 I); /* print INCIDENCE */
void fprI(int N,Inci64 I); /* print INCIDENCE to file */

int Inci64_LE(Inci64 A,Inci64 B);
int Inci64_LT(Inci64 A,Inci64 B);

void PRNtriang(triang *SR,const char *c);
/*
print (Inci64) triangulation
*/

/*** from SingularInput.c ***/
void HyperSurfSingular(PolyPointList *P,triang *T, triang *SR ,MORI_Flags *_Flag ,
		FibW *F, int *cp);
/*
Interface between the Mori and SINGULAR. Given the polytope P and the linear relations
among their points, the triangulation T and the SR ideal it generates the the input
for SINGULAR.
*/


/*======== 	dependences from other modules  ======== */

/*** from Polynf.c ***/
void IP_Simplex_Fiber(Long PM[][POLY_Dmax], int p, int d, /* need PM[i]!=0 */ FibW *F, int Wmax, int CD);
/*
It analyzes the IP simplices among (1 <= codim <= CD) points of P*. 
Given the matrix CM of coordinates of p points in Z^d, the list W[i] of *nw
weight systems corresponding to IP-simplices spanned by the points in CM is
created. If codim!=0 only the IP-simplices with dimension > 1 and codimension
between 1 and codim are computed. It is assumed that p<=VERT_Nmax and that W 
can hold at least Wmax sets of coefficients.
*/

Long SimplexVolume(Long *V[POLY_Dmax+1],int d);

