#define Nef_Max 	500000
#define NP_Max          500000
#define W_Nmax (POLY_Dmax+1)
#define MAXSTRING	100

#undef	WRITE_CWS
#define	WRITE_CWS

#define Pos_Max   (POLY_Dmax + 2)
#define FIB_POINT_Nmax  VERT_Nmax

typedef struct {
  Long W[FIB_Nmax][FIB_POINT_Nmax];
  Long VM[FIB_POINT_Nmax][POLY_Dmax];
  int nw;
  int nv;
  int d;
  int Wmax;
} LInfo;

struct Poset_Element {
  int num, dim;
};

struct Interval {
  int min, max;
};

typedef struct Interval Interval;

typedef struct {
  struct Interval *L;
  int n;
} Interval_List;

typedef struct Poset_Element Poset_Element;

typedef struct {
  struct Poset_Element x, y;
} Poset;

typedef struct {
  struct Poset_Element *L;
  int n;
} Poset_Element_List;

typedef struct {
  int nface[Pos_Max];
  int dim;
  INCI edge[Pos_Max][FACE_Nmax];
} Cone;

typedef struct {
  Long S[2*Pos_Max];
} SPoly;

typedef struct {
  Long B[Pos_Max][Pos_Max];
} BPoly;

typedef struct {
  int E[4*(Pos_Max)][4*(Pos_Max)];
} EPoly;

typedef struct {
    Long x[POINT_Nmax][W_Nmax];
    int N, np;
} AmbiPointList;

typedef struct {
  int n;
  int nv;
  int codim;
  int S[Nef_Max][VERT_Nmax];
  int DirProduct[Nef_Max];
  int Proj[Nef_Max];
  int DProj[Nef_Max];
} PartList;

typedef struct {
  int n;
  int nv;
  int S[Nef_Max][VERT_Nmax];
} Part;

typedef struct {
  int n, y, w, p, t, S, Lv, Lp, N, u, d, g, VP, B, T, H, dd, gd, 
    noconvex, Msum, Sym, V, Rv, Test, Sort, Dir, Proj, f, G;		
} Flags;

typedef struct {
	int noconvex, Sym, Test, Sort;		
} NEF_Flags;

struct Vector {
  Long x[POLY_Dmax];
};

typedef struct Vector Vector ;

typedef struct {
  struct Vector *L;
  int n;
  Long np, NP_max;						} DYN_PPL;

void part_nef(PolyPointList *, VertexNumList *, EqList *, PartList *,
	      int *, NEF_Flags *);

void Make_E_Poly(FILE *, CWS *, PolyPointList *, VertexNumList *, EqList *,
		 int *, Flags *, int *);

void Mink_WPCICY(AmbiPointList * _AP_1, AmbiPointList * _AP_2,
		 AmbiPointList * _AP);

int  IsDigit(char);

int  IntSqrt(int q);

void Die(char *);

void Print_CWS_Zinfo(CWS *);

void AnalyseGorensteinCone(CWS *_CW,  PolyPointList *_P, VertexNumList *_V, 
			   EqList *_E, int *_codim, Flags * _F);
