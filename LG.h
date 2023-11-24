#define	 WZinput 	(1)	/* WZ-input (in progress)  */

#define  W_Nmax		(POLY_Dmax+1)

#define  min(a,b)  	(((a)<(b)) ? (a) : (b))
#define  max(a,b)  	(((a)>(b)) ? (a) : (b))

#define	 Pint		int	   /* type of coefficients in PolyCoeffList */
typedef struct {int n; int *e; Pint *c; int A;}		  PoCoLi;  /* e=exp */

void AllocPoCoLi(PoCoLi *P);		      /* allocate e[P.A] and c[P.A] */
void Free_PoCoLi(PoCoLi *P);				/* free P.e and P.c */
void Poly_Sum(PoCoLi *A,PoCoLi *B,PoCoLi *S);			 /* S = A+B */
void Poly_Dif(PoCoLi *A,PoCoLi *B,PoCoLi *D);			 /* D = A-B */
void PolyProd(PoCoLi *A,PoCoLi *B,PoCoLi *AB);			/* AB = A*B */
int  BottomUpQuot(PoCoLi *N,PoCoLi *D,PoCoLi *Q,PoCoLi *R);    /* Q*D = N-R */
void PolyCopy(PoCoLi *X,PoCoLi *Y);				/*  Y = X   */
void PrintPoCoLi(PoCoLi *P);
void UnitPoly(PoCoLi *P);
void Init1_xN(PoCoLi *P,int N);					 /* 1 - x^N */
void PoincarePoly(int N,int *w,int d, PoCoLi *PP, PoCoLi *Naux, PoCoLi *Raux);

int  IsDigit(char c);

typedef struct {int d, N, z[POLY_Dmax][W_Nmax], m[POLY_Dmax], M, r, R;/* Ref */
		Long w[W_Nmax], B[W_Nmax][POLY_Dmax], A[W_Nmax], rI[POLY_Dmax];
		PolyPointList *P;} /* Eq: Ei.c=Ai Ei.a[]=Bi[]} */
		/* 0<=A+B*x  r=sum(w)/d  rI=IP(r*P)  n=(r,rI) */	Weight;

typedef struct {int D,E,sts; int h[POLY_Dmax][POLY_Dmax];}		VaHo;

int Read_W_PP(Weight *, PolyPointList *);
int Trans_Check(Weight);
void LGO_VaHo(Weight *,VaHo *);
void Write_Weight(Weight *_W);
void Write_WH(Weight *_W, BaHo *_BH, VaHo *_VH, int rc, int tc,
	      PolyPointList *_P, VertexNumList *_V, EqList *_E);
void Make_Poly_Points(Weight *_W_in, PolyPointList *_PP);

Long V_to_G_GI(Long *V,int d, Long G[][POLY_Dmax],Long GI[][POLY_Dmax]);
