/*  ======================================================================  */
/*  ==========                                                  ==========  */
/*  ==========  test-              N E F . C                    ==========  */
/*  ==========                                                  ==========  */
/*  ==========  Modules: IO, Rat, Points, Vertex, Convex, Hodge ==========  */
/*  ==========                                                  ==========  */
/*  ======================================================================  */

#include "Global.h"
#include "Nef.h"
#include "LG.h"

/*  ==========            l o c a l  T Y P E D E F s            ==========  */

typedef struct { Long x[POLY_Dmax][W_Nmax]; int N, n; }	     AmbiLatticeBasis;
typedef struct { Long x[AMBI_Dmax][AMBI_Dmax]; int n, N; }   CWLatticeBasis;
typedef struct { Long P[POINT_Nmax]; Long n; }               Pstat;

/*  ==========          l o c a l  P R O T O T Y P E s          ==========  */

int  ReadCwsPp(CWS *_CW, PolyPointList *_P, int codim, int index);

void Compute_E_Poly(EPoly *_EP,  
		    PolyPointList * _P_D, VertexNumList * _V_D, EqList * _E_D,
		    PolyPointList * _P_N, VertexNumList * _V_N, EqList * _E_N,
		    int *_codim, Flags * _F, time_t *_Tstart, clock_t *_Cstart);

void Sort_PPL(PolyPointList *_P, VertexNumList *_V);

void NormTriangularBasis(AmbiLatticeBasis * _B);

void WeightLatticeBasis(Weight * _w, AmbiLatticeBasis * _B);

void WeightMakePoints(Weight * _W, AmbiPointList * _P);

void MakeRefWeights(int N, int from_d, int to_d);

void ChangeToTrianBasis(AmbiPointList *, AmbiLatticeBasis *,
			PolyPointList *);
int IN_WEIGHT(Weight *, CWS *, int *, PolyPointList *, Flags *, int);

void OUT_CWS(CWS *, int *, int *);

void Print_VP(PolyPointList *, VertexNumList *, int, int, Pstat *);

void Print_Pstat(Pstat *, int, int, int);

void Die(char *);

//=== NOW INIZIO ===//
void Print_Nefinfo(PartList *_PTL, /* Flags *_F,*/ time_t *_Tstart, clock_t *_Cstart);
//=== NOW FINE ===//

FILE *inFILE, *outFILE;

#define OSL (31)  /* opt_string's length */

void  PrintNefUsage(char *c){
  int i;
  char *opt_string[OSL]={
    "Options: -h        prints this information",
    "         -f or -   use as filter; otherwise parameters denote I/O files",
    "         -N        input is in N-lattice (default is M)",
    "         -H        gives full list of Hodge numbers",
    "         -Lv       prints L vector of Vertices (in N-lattice)",
    "         -Lp       prints L vector of Points (in N-lattice)",
    "         -p        prints only partitions, no Hodge numbers",
    "         -D        calculates also direct products",
    "         -P        calculates also projections",
    "         -t        full time info",
    "         -cCODIM   codimension (default = 2)",
    "         -Fcodim   fibrations up to codim (default = 2)",
    "         -y        prints poly/CWS in M lattice if it has nef-partitions",
    "         -S        information about #points calculated in S-Poly",
    "         -T        checks Serre-duality          ",
    "         -s        don't remove symmetric nef-partitions   ",
    "         -n        prints polytope only if it has nef-partitions",
    "         -v        prints vertices and #points of input polytope in one",
    "                   line; with -u, -l the output is limited by #points:",
    "             -uPOINTS  ... upper limit of #points (default = POINT_Nmax)",
    "             -lPOINTS  ... lower limit of #points (default = 0)",
    "         -m        starts with [d  w1 w2 ... wk d=d_1 d_2 (Minkowski sum)",
    "         -R        prints vertices of input if not reflexive",
    "         -V        prints vertices of N-lattice polytope",
    "         -Q        only direct products (up to lattice Quotient)",
    "         -gNUMBER  prints points of Gorenstein polytope in N-lattice",
    "         -dNUMBER  prints points of Gorenstein polytope in M-lattice",
    "               if NUMBER = 0 ... no            0/1 info",
    "               if NUMBER = 1 ... no redundant  0/1 info (=default)",
    "               if NUMBER = 2 ... full          0/1 info",
    "         -G        Gorenstein cone: input <-> support polytope"
  };
  puts("");
  printf("This is '%s':  calculate Hodge numbers of nef-partitions\n",c);
  printf("Usage:   %s <Options>\n", c);
  puts("");
  for (i=0;i<OSL;i++) puts(opt_string[i]);
  exit(0);
}

int READ_INT(int *n, int narg, char *fn[], char *s){
  char *a;

  if((fn[*n][2]==0) && (narg>*n+1)) a=fn[++*n]; else a=&fn[*n][2];
  if(*a == 0) PrintNefUsage(fn[0]); 
  if(!IsDigit(*a)){ printf("after %s there must be digit(s)!\n",s); exit(0);} 
  return atoi(a);
}

int Make_Mirror(EPoly *_EP, int h[][POLY_Dmax], int D, int dim);

int main(int narg, char *fn[])
{
    Flags F;
    int N = 0, n = 0, FilterFlag = 0, VPmax = POINT_Nmax-1, VPmin = 0;
    CWS CW;
    Weight W;
    int D[2], codim=2;
    Pstat *_PS;
    VertexNumList *_V;
    EqList *_E;
    PolyPointList *_P;

    _P = (PolyPointList *) malloc(sizeof(PolyPointList));
    if (_P == NULL) Die("Unable to allocate space for _P");
    _V = (VertexNumList *) malloc(sizeof(VertexNumList));
    if (_V == NULL) Die("Unable to alloc space for VertexNumList _V");
    _E = (EqList *) malloc(sizeof(EqList));
    if (_E == NULL) Die("Unable to alloc space for EqList _E");
    _PS = (Pstat *) malloc(sizeof(Pstat));
    if (_PS == NULL) Die("Unable to alloc space for Pstat _PS");
 
    F.p = 0;
    F.Lv = 0;
    F.Lp = 0;
    F.t = 0;
    F.S = 0;
    F.T = 0;
    F.N = 0;
    F.H = 0;
    F.Msum = 0;
    F.Sym = 1;
    F.Rv = 0;
    F.V = 0;
    F.Test = 0;
    F.Sort = 1;
    F.Dir = 0;
    F.B = 1;
    F.VP = 0;
    F.Proj = 0;
    F.w = 0;
    F.f = 0;
    F.n = 0;
    F.g = 0; F.gd=1;
    F.d = 0; F.dd=1;
    F.noconvex=0;
    F.y=0;
    F.G = 0;

    while (narg > ++n)
	if (fn[n][0] != '-')
	    break;
	else {
	    if ((fn[n][1] == '?') || (fn[n][1] == 'h')) {
		PrintNefUsage(fn[0]);
		return 0;
	    }
	    else if (fn[n][1] == 'x'){
	        PrintNefUsage(fn[0]);
		return 0;
	    }
	    else if ((fn[n][1] == 'f') || (fn[n][1] == 0))
		FilterFlag = 1;
	    else if (fn[n][1] == 'S')
		F.S = 1;
	    else if (fn[n][1] == 'H')
		F.H = 1;
	    else if (fn[n][1] == 't')
		F.t = 1;
	    else if (fn[n][1] == 'W')
		F.w = 1;
	    else if (fn[n][1] == 's')
		F.Sym = 0;
	    else if (fn[n][1] == 'T')
		F.T = 1;
	    else if (fn[n][1] == 'N')
		F.N = 1;
	    else if (fn[n][1] == 'v')
		F.VP = 1;
	    else if (fn[n][1] == 'n'){
		F.n = 1;
		F.p = 1;
	    }
	    else if (fn[n][1] == 'u')
	        VPmax = READ_INT(&n, narg, fn, "u");
	    else if (fn[n][1] == 'l')
	        VPmin = READ_INT(&n, narg, fn, "l");
	    else if (fn[n][1] == 'm')
		F.Msum = 1;
	    else if (fn[n][1] == 'g'){
		F.p = 1;
		F.g = 1;
		if(IsDigit(fn[n][2]))
		  F.gd = atoi(&fn[n][2]);
	    }
	    else if (fn[n][1] == 'd'){
		F.d = 1;
		F.p = 1;
		if(IsDigit(fn[n][2]))
		  F.dd = atoi(&fn[n][2]);
	    }
	    else if (fn[n][1] == 'R')
		F.Rv = 1;
	    else if (fn[n][1] == 'V')
		F.V = 1;
	    else if (fn[n][1] == 'L'){
	        if (fn[n][2] == 'p')
		    F.Lp = 1;
		if (fn[n][2] == 'v')
		    F.Lv = 1;
	    }
	    else if (fn[n][1] == 'k')
	        F.noconvex=1;
	    else if (fn[n][1] == 'p')
		F.p = 1;
	    else if (fn[n][1] == 'y')
		F.y = 1;
	    else if (fn[n][1] == 'D')
		F.Dir = 1;
            else if (fn[n][1] == 'Q')
                F.Dir = 2;
	    else if (fn[n][1] == 'P')
		F.Proj = 1;
	    else if (fn[n][1] == 'c')
	       codim = READ_INT(&n, narg, fn, "c");
	    else if (fn[n][1] == 'F'){
	        F.f = 2;
		if(fn[n][2]!=0) F.f = READ_INT(&n, narg, fn, "F");
	    }
	    else if (fn[n][1] == 'G')
	      F.G = 1;
	    else {
	      printf("Unknown option '-%c'; use -h for help\n",fn[n][1]); 
	      exit(0);}
	}
    n--;
    if (FilterFlag) {
	inFILE = NULL;
	outFILE = stdout;
    } 
    else {
	if (narg > ++n)
	    inFILE = fopen(fn[n], "r");
	else
	    inFILE = stdin;
	if (inFILE == NULL) {
	    printf("Input file %s not found!\n", fn[n]);
	    exit(0);
	}
	if (narg > ++n)
	    outFILE = fopen(fn[n], "w");
	else
	    outFILE = stdout;
    }	
    while (IN_WEIGHT(&W, &CW, D, _P, &F, codim)) {
      /* _P is the M-lattice polytope */
      if (F.G) AnalyseGorensteinCone(&CW,_P,_V,_E,&codim,&F);
      else if (Ref_Check(_P, _V, _E)){
	int nv=_V->nv, ne=_E->ne;
	Long PM[EQUA_Nmax][VERT_Nmax];
     
	Make_VEPM(_P,_V,_E, PM);
	Complete_Poly(PM,_E,_V->nv,_P);
	Find_Equations(_P,_V,_E);
	/*Print_PPL(_P,"vorher");*/
	Sort_VL(_V);
	Sort_PPL(_P, _V);
	/*Print_PPL(_P,"nachher");*/
	assert(nv==_V->nv && ne==_E->ne);
	if (!F.VP){
#ifdef  WRITE_CWS
          OUT_CWS(&CW, D, &F.Msum);
#endif
	  if (POLY_Dmax  < (_P->n + codim - 1)){
	    printf("Please increase POLY_Dmax to at least %d = %d + %d - 1\n",
		   (_P->n + codim - 1), _P->n, codim);
	    printf("(%s requires POLY_Dmax >= dim N + codim - 1)\n",
		   fn[0]);
	    exit(0);	  }
	  Make_E_Poly(outFILE, &CW, _P, _V, _E, &codim, &F, &D[0]);	}
	else{
	  N++;
	  Print_VP(_P, _V, VPmax, VPmin, _PS);	}       }
      else{
	if ((F.Rv == 1) || ((F.V == 1)&&(F.N == 1))){
	  Find_Equations(_P,_V,_E);
	  Print_VL(_P, _V, "Vertices of input polytope:");      }    }}
    if (F.VP){
      assert(VPmax < POINT_Nmax); assert(VPmax >= VPmin); 
      Print_Pstat(_PS, N, VPmax, VPmin);    }
    free(_E); free(_V); free(_P); free(_PS);
    return 0;
}

void Print_Pstat(Pstat *_PS, int N, int VPmax, int VPmin){
  int i;
  
  fprintf(outFILE,"\n\n%d  of  %d\n\n",(int) _PS->n,(int) N);
  for(i=VPmin; i<=VPmax; i++) 
    if (_PS->P[i] != 0)
      fprintf(outFILE,"%4d# %4d\n",(int) i, (int) _PS->P[i]);
}

void Print_VP(PolyPointList *_P, VertexNumList *_V, int VPmax, int VPmin, 
	Pstat *_PS){
  int i,j;

  if((_P->np <= VPmax) && (_P->np >= VPmin)){
    _PS->P[_P->np] ++; _PS->n ++;
    if(_V->nv>20){
      fprintf(outFILE,"%d %d P:%d E",_V->nv,_P->n,_P->np);
      for(i=0;i<_V->nv;i++) {
	for(j=0;j<_P->n;j++) fprintf(outFILE,"%d ",(int) _P->x[_V->v[i]][j]); 
	if(i!=(_V->nv-1)) fprintf(outFILE,"E");}}
    else {
      fprintf(outFILE,"%d %d P:%d E",_P->n,_V->nv,_P->np);
      for(i=0;i<_P->n;i++) {
	for(j=0;j<_V->nv;j++) fprintf(outFILE," %4d",(int) _P->x[_V->v[j]][i]);
	if(i!=(_P->n-1)) fprintf(outFILE,"E");}}
    fprintf(outFILE,"\n");
  }
}
         
void Make_Sort(Long * _N, Long * _n, Long * _B_num, Long * _c_num,
	       Long * _c_less)
{
  Long N, n;

  n = *_n;
  N = *_N;
  while (n != 0) {
    if (_c_less[n - 1] == N) {
      _B_num[N + n - 1] = _c_num[n - 1];
      n--;
    } else {
      _B_num[N + n - 1] = _B_num[N - 1];
      N--;
    }
  }
  *_N += *_n;
  *_n = 0;
}

Long Make_Bi_section(DYN_PPL *_MP, Long *_N_List, Long *_N, Long *x)
{
  Long max_pos, min_pos, pos, diff;
  int d;

  max_pos = *_N;
  min_pos = -1;
  while ((max_pos - min_pos) > 1) {
    pos = (max_pos + min_pos) / 2;
    diff = 0;
    d = _MP->n;
    while ((d != 0) && (diff == 0)) {
      diff = (x[d - 1] - _MP->L[_N_List[pos]].x[d - 1]);
      d--;
    }
    if (diff > 0)
      min_pos = pos;
    else if (diff < 0)
      max_pos = pos;
    else
      max_pos = -1;
  }
  return max_pos;
}

void Mink_WPCICY(AmbiPointList * _AP_1, AmbiPointList * _AP_2,
		 AmbiPointList * _AP)
{
    int l;
    Long *_x, *_B_num, *_c_num, *_c_less, n = 0, N = 0, Num,
	num, P_max, p_max, j, k, m;
    DYN_PPL B;

    P_max = _AP_1->np * _AP_2->np;
    p_max = ((Long) IntSqrt(P_max));

    _x = (Long *) calloc(_AP_1->N, sizeof(Long));
    B.L = (Vector *) calloc(P_max, sizeof(Vector));
    _B_num = (Long *) calloc(P_max, sizeof(Long));
    _c_num = (Long *) calloc(p_max, sizeof(Long));
    _c_less = (Long *) calloc(p_max, sizeof(Long));

    assert((B.L != NULL) && (_B_num != NULL) && (_c_num != NULL) &&
	   (_c_less != NULL) && (_x != NULL));

    B.n = _AP_1->N; B.NP_max = P_max;
    for (j = 0; j < _AP_1->np; j++)
	for (k = 0; k < _AP_2->np; k++) {
	    if (n == p_max)
		Make_Sort(&N, &n, _B_num, _c_num, _c_less);
	    for (l = 0; l < _AP_1->N; l++)
		_x[l] = _AP_1->x[j][l] + _AP_2->x[k][l];
	    Num = Make_Bi_section(&B, _B_num, &N, _x);
	    if (Num >= 0) {
		num = Make_Bi_section(&B, _c_num, &n, _x);
		if (num >= 0) {
		    assert((n + N) < P_max);
		    for (l = 0; l < _AP_1->N; l++)
			B.L[(n + N)].x[l] = _x[l];
		    for (m = n; m > num; m--) {
			_c_num[m] = _c_num[m - 1];
			_c_less[m] = _c_less[m - 1];
		    }
		    _c_num[num] = N + n;
		    _c_less[num] = Num;
		    n++;
		}
	    }
	}
    assert((n + N) <= POINT_Nmax);
    B.np = n + N;
    _AP->np = B.np;
    _AP->N = _AP_1->N;
    for (j = 0; j < (N + n); j++)
	for (l = 0; l < _AP_1->N; l++)
	    _AP->x[j][l] = B.L[j].x[l];

    free(_x);
    free(B.L);
    free(_B_num);
    free(_c_num);
    free(_c_less);
}

void Make_Poly_WPCICY(Weight * _W, int *_D, PolyPointList * _PP)
{
    AmbiLatticeBasis B;

    AmbiPointList *_AP_1 = (AmbiPointList *) malloc(sizeof(AmbiPointList));
    AmbiPointList *_AP_2 = (AmbiPointList *) malloc(sizeof(AmbiPointList));
    AmbiPointList *_AP = (AmbiPointList *) malloc(sizeof(AmbiPointList));
    assert((_AP_1 != NULL) && (_AP_2 != NULL) && (_AP != NULL));


    WeightLatticeBasis(_W, &B);

    _W->d = _D[0];
    WeightMakePoints(_W, _AP_1);
    assert(POINT_Nmax >= _AP_1->np);

    _W->d = _D[1];
    WeightMakePoints(_W, _AP_2);
    assert(POINT_Nmax >= _AP_2->np);

    _W->d += _D[0];

    Mink_WPCICY(_AP_1, _AP_2, _AP);
    ChangeToTrianBasis(_AP, &B, _PP);

    free(_AP);
    free(_AP_1);
    free(_AP_2);
    return;
}

void Make_CW_WPCICY(Weight * _W, CWS * _CW)
{
    int i;

    _CW->nw = 1;
    _CW->N = _W->N;
    _CW->d[0] = _W->d;
    for (i = 0; i < _CW->N; i++)
	_CW->W[0][i] = _W->w[i];
    _CW->nz = 0;
}

int Read_WPCICY(Weight * _W, int *_D)
     /* read "d" and "w_i" till sum=d or non-digit */
{
    char c;
    int long nl, sum, FilterFlag = (inFILE == NULL);
    if (inFILE == stdin)
	printf("type degrees and weights [d  w1 w2 ... wk d=d_1 d_2]: ");
    else if (FilterFlag) inFILE = stdin;
    c = fgetc(inFILE);
    if (!IsDigit(c))
	return 0;
    ungetc(c, inFILE);
    fscanf(inFILE, "%ld", &nl);
    _W->d = nl;
    sum = nl;

    for (_W->N = 0; _W->N < W_Nmax; sum -= nl) {
	assert(_W->N < W_Nmax);
	while (' ' == (c = fgetc(inFILE)));
	ungetc(c, inFILE);
	if (IsDigit(c)) {
	    fscanf(inFILE, "%ld", &nl);
	    _W->w[(_W->N)++] = nl;
	} else
	    break;
    }
    while (!IsDigit(c = fgetc(inFILE)));
    ungetc(c, inFILE);
    fscanf(inFILE, "%d", &_D[1]);
    while (' ' == (c = fgetc(inFILE)));
    ungetc(c, inFILE);
    fscanf(inFILE, "%d", &_D[0]);

    while (fgetc(inFILE) - '\n')
	if (feof(inFILE))
	    return 0;		/* read to EOL */

    if (_W->N > POLY_Dmax) {
      printf("Please increase POLY_Dmax ");
      printf("(POLY_Dmax >= number of weights is required)\n");
      exit(0);}

    assert((_D[0] + _D[1]) == _W->d);

    if (_W->N < 2) {
	puts("I need at least 2 weights!");
	exit(0);
    }
    if (FilterFlag) inFILE = NULL;
    return 1;
}

int Make_WPCICY(Weight * _W, CWS * _CW, int *_D, PolyPointList * _P)
{
    int r = Read_WPCICY(_W, _D);

    Make_Poly_WPCICY(_W, _D, _P);
    Make_CW_WPCICY(_W, _CW);
    return r;
}

void Make_RGC_Points(CWS *Cin, PolyPointList *_P);

int IN_WEIGHT(Weight * _W, CWS * _CW, int *_D, PolyPointList * _P,
	      Flags *_F, int codim){
  if (_F->Msum) return Make_WPCICY(_W, _CW, _D, _P);
  if (_F->G) return ReadCwsPp(_CW, _P , 1, codim);
  return ReadCwsPp(_CW, _P, codim, 1);
}

void OUT_CWS(CWS * _W, int *_D, int *_M_Flag)
{
    int i, j;

    for (i = 0; i < _W->nw; i++) {
	fprintf(outFILE, "%d ", (int) _W->d[i]);
	for (j = 0; j < _W->N; j++)
	    fprintf(outFILE, "%d ", (int) _W->W[i][j]);
	if (i + 1 < _W->nw)
	    fprintf(outFILE, " ");
    }
    if (*_M_Flag)
	fprintf(outFILE, "d=%d %d ", (int) _D[1], _D[0]);
    Print_CWS_Zinfo(_W);
}
