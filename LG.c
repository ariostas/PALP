#include "Global.h"
#include "Rat.h"
#include "LG.h"

#define SHOW_b01_TWIST	(0)
#define Tout		(0)
#define	DET_WARN_ONLY	(0)		     /* continue if group has det!=1 */

#define	COEFF_Nmax	(d*D+2*N)

#define	ABBREV_POLY_PRINT	(4)	  /* !=0 => #(leading/trailing terms */

#define NO_COORD_IMPROVEMENT		/* switch off weight permutation */
#undef	W_PERM_CODE
int  Is_Gen_CY(int index, PolyPointList *P)
{    int i; Long *IP=P->x[P->np]; VertexNumList V;     /* IP = IP(index * P) */
     EqList *E = (EqList *) malloc(sizeof(EqList)); assert(E!=NULL);
     Find_Equations(P,&V,E); for(i=0;i<E->ne;i++) E->e[i].c *= index;
     for(i=0;i<E->ne;i++) if(1!=Eval_Eq_on_V(&E->e[i],IP,P->n)) break;
     free(E); return i==E->ne;
}
void WZerror(char *c)
{    printf("Format error %s in Read_WZeight\n",c);exit(0);
}    
int  auxString2SInt(char *c,int *n)
{    int j=0,neg=0; *n=0;     if(c[0]=='-')if(('0'<c[1])&&(c[1]<='9')) neg=j=1;
     while(('0'<=c[j])&&(c[j]<='9')) *n = 10*(*n)+c[j++]-'0';
     if(neg)*n *=-1;
     if(j) while(c[j]==' ') j++;
     return j;
}
int  Read_WZ_PP(Weight *WZ)     /* read "d w_i" [ or "w_i d" if last=max ] */
{    int i,j,k,a,n,d,shift=1, I[W_Nmax+2], *nz=&WZ->M; 
     int FilterFlag=(inFILE==NULL); char C, c[999],b=' ';
     Long BM[W_Nmax][W_Nmax], *B[W_Nmax], Wa[POLY_Dmax], Za[POLY_Dmax], 
	F[W_Nmax], G[POLY_Dmax][POLY_Dmax], GI[POLY_Dmax][POLY_Dmax],X;
     if(FilterFlag) inFILE=stdin; else
     if(inFILE==stdin) printf("type degree and weights  [d  w1 w2 ...]: "); 
     C=fgetc(inFILE); if( !IsDigit(C) ) return 0; ungetc(C,inFILE);	 *nz=0;
     fscanf(inFILE,"%d",I);
     for(i=1; i<W_Nmax+2; i++)
     {	while(' '==(C=fgetc(inFILE))); ungetc(C,inFILE);
	if(IsDigit(C)) fscanf(inFILE,"%d",&I[i]); else break;
     }	WZ->N=i-1; if(WZ->N>W_Nmax) { puts("Increase POLY_Dmax"); exit(0); }
     for(i=0;i<=WZ->N;i++) assert(I[i]>0);
     if(I[WZ->N]>I[0]) { WZ->d=I[WZ->N]; shift=0; } else WZ->d=I[0];
     for(i=0;i<WZ->N;i++) {WZ->w[i]=I[i+shift]; assert(WZ->w[i]<WZ->d);} 
     WZ->r=0; for(i=0;i<WZ->N;i++) WZ->r+=WZ->w[i]; 
     if(WZ->r%WZ->d)WZ->r=0; else WZ->r/=WZ->d;
     for(n=0;n<999;n++)					/* read /Z*: * * * */
     {  c[n]=fgetc(inFILE); 
	if(feof(inFILE)) {if(FilterFlag) inFILE=NULL; return 0;} 
	if(c[n]=='\n') break;
     }  if(n==999) {puts("Out of space in Read_WZeight");exit(0);} i=0;
     while(c[i]==b)i++;
     if((c[i]=='=')&&(c[i+1]=='d'))i+=2;
     while(c[i]==b)i++; 
     while(i<n)
     {  int j, k, s; if((c[i]!='/')||(c[i+1]!='Z')) break;
        i+=2; assert(*nz<POLY_Dmax);		 WZ->m[*nz]=0;
	while((i<n)&&IsDigit(c[i])) WZ->m[*nz]=10*WZ->m[*nz]+c[i++]-'0';
	j=0;
	if(WZ->m[*nz]<=0) WZerror("Order not positive");
	    if(c[i+j]!=':') WZerror(":"); else {i+=j+1; while(c[i]==b)i++;}
            for(k=0;k<WZ->N;k++)
            {   if((j=auxString2SInt(&c[i],&WZ->z[*nz][k]))) {if(j) i+=j; 
                else WZerror("missing");}
            }   s=0;
            for(k=0;k<WZ->N;k++) s+=WZ->z[*nz][k]; 
            /* if(s % WZ->m[*nz])  WZerror("det!=1") */;
            (WZ->M)++;
     }	/* Eq_i::{c=Ai a[]=Bi[]}  0<=A+B*x  r=sum(w)/d  rI=IP(r*P)  n=(r,rI) */

     a=WZ->N; for(i=0;i<a;i++)B[i]=BM[i]; assert( 1 == W_to_GLZ(WZ->w,&a,B) );
     if(1==WZ->r) for(i=0;i<a;i++) WZ->A[i]=F[i]=1;	/* Z.(X-F) \cong 0 */
     else {for(i=0;i<a;i++) { WZ->A[i]=WZ->d*B[0][i]; F[i]=0;} 
	for(j=0;j<WZ->M;j++) {Long s=0; for(n=0;n<a;n++) s+=WZ->z[j][n]; 
	if(s%WZ->m[j]){puts("det=1 required if Gorenstein index r>1:");
#if	DET_WARN_ONLY
	WZ->M=0; }
#else
	    Write_Weight(WZ); exit(0);}
#endif
     }}	d=a-1; n=0; shift=1;
     for(i=0;i<a;i++)for(j=0;j<d;j++)WZ->B[i][j]=B[j+1][i]; /* LATTICE BASIS */
     for(i=a-1;i>0;i--) { X = 1 - WZ->r * WZ->A[i]; 
	for(j=i;j<d;j++) X -= WZ->rI[j]*WZ->B[i][j]; 
	assert(0==(X%WZ->B[i][i-1])); WZ->rI[i-1]=X/WZ->B[i][i-1];}
     for(i=0;i<*nz;i++)				    /* divide by symmetries: */
     {	Long g=WZ->m[i]; X=0; for(j=0;j<d;j++){Za[j]=0; 
	for(k=0;k<a;k++)Za[j]+=WZ->z[i][k]*WZ->B[k][j];
	g=NNgcd(Za[j],g);}
if(Tout){printf("WZ->m[%d]=%d g=%ld  ",i,WZ->m[i],g);if(g==WZ->m[i])puts("");}
	if(g==WZ->m[i]) continue;
	if(g>1) for(k=0;k<a;k++) Za[j]/=g; 
	V_to_G_GI(Za,d,G,GI); for(k=0;k<a;k++)
	{   int l; for(l=0;l<d;l++) Wa[l]=WZ->B[k][l]; for(l=0;l<d;l++) 
	    {	WZ->B[k][l]=0; for(j=0;j<d;j++) WZ->B[k][l] += Wa[j]*G[l][j];
	    }		/* shift A to A'=A+B.x with x=dx=Z_ia*(F_a-A_a) */
	    X+=WZ->z[i][k]*(WZ->A[k]-F[k]); 
	}   X%=WZ->m[i];/* IP = r*A+B.rI, B'=B.G^T => rI' = GI^T . rI + r dx */
	if(X<-WZ->m[i]/2)X+=WZ->m[i]; else if(X>WZ->m[i]/2)X-=WZ->m[i];
	if(Tout) printf("X=%ld\n",X);
	for(k=0;k<a;k++) WZ->A[k]-=X*WZ->B[k][0];
	for(j=0;j<d;j++) {Za[j]=WZ->rI[j]; WZ->rI[j]=0;} WZ->rI[0] =X*WZ->r;
	for(j=0;j<d;j++) for(k=0;k<d;k++) WZ->rI[j]+=Za[k]*GI[k][j];
	g=WZ->m[i]/g; for(k=0;k<a;k++) WZ->B[k][0]*=g;    /* goto sublattice */
	assert((WZ->rI[0]%g)==0); WZ->rI[0]/=g; shift*=g;
     }	X=1; for(j=0;j<*nz;j++) X*=WZ->m[j]; assert(X%shift==0); shift=X/shift;

if(Tout){for(j=0;j<a;j++)printf("%2ld ", WZ->A[j]); printf("=A rI=");
  for(j=0;j<d;j++) printf("%2ld ", WZ->rI[j]);
  puts("");
  for(k=0;k<a;k++){for(j=0;j<d;j++)printf("%2ld ",WZ->B[k][j]);puts("=B");} 
  for(j=0;j<a;j++){X=0;for(k=0;k<d;k++)X+=WZ->B[j][k]*WZ->rI[k];
    printf("%2ld ",X);}
  puts("=B.rI");Write_Weight(WZ);} 

     {	PairMat VPM; VertexNumList V; PolyPointList *P=WZ->P;
	EqList *E=(EqList *)malloc(sizeof(EqList));assert(E!=NULL);
	for(i=0;i<a;i++){E->e[i].c=WZ->A[i]; VPM[i][0]=WZ->d/WZ->w[i];
	    for(j=0;j<d;j++)E->e[i].a[j]=WZ->B[i][j];} E->ne=a; 
	WZ->P->n=d; WZ->P->np=0; Complete_Poly(VPM,E,1,WZ->P);
/*	    assert(POINT_Nmax>WZ->P->np); / * Print_PPL(WZ->P,""); obsolete */
/*	    for(i=0;i<d;i++) WZ->P->x[WZ->P->np][i]=WZ->rI[i]; / * obsolete */
	for(i=0;i<d;i++) if(WZ->rI[i]) break;
	if(i<d) /* rI!=0 => make P0=0 */
	{   for(i=1;i<P->np;i++) for(j=0;j<d;j++) P->x[i][j]-=P->x[0][j];
	    for(i=0;i<a;i++)for(j=0;j<d;j++) WZ->A[i]+=WZ->B[i][j]*P->x[0][j];
	    for(i=0;i<d;i++) {WZ->rI[i]-=WZ->r*P->x[0][i]; P->x[0][i]=0;}
	}	
     	Find_Equations(WZ->P,&V,E); for(i=0;i<E->ne;i++) E->e[i].c *= WZ->r;
     	for(i=0;i<E->ne;i++) if(1!=Eval_Eq_on_V(&E->e[i],WZ->rI,d)) break;
     	WZ->R = (i==E->ne);  free(E);
     }	if(shift!=1) 				/* normalize sublattice data */
     {	Long *PB[W_Nmax], Z[POLY_Dmax][VERT_Nmax], M[POLY_Dmax], D[POLY_Dmax]; 
	for(i=0;i<a;i++) PB[i]=WZ->B[i]; 
	Sublattice_Basis(d,a,PB,Z,M,&k,G,D);
	for(i=0;i<k;i++){X=0;for(j=0;j<a;j++)X+=Z[i][j];assert(0==(X%M[i]));}
	for(i=0;i<k;i++){WZ->m[i]=M[i];	for(j=0;j<a;j++)WZ->z[i][j]=Z[i][j];}
	WZ->M=k;	/* printf("shift=%d   rank=%d \n",shift,k); */
     }	X=0; for(j=0;j<a;j++) X+=WZ->w[j]*WZ->A[j];		    /* X=w.A */
     for(i=0;i<WZ->M;i++) {Long A=0,M=WZ->m[i],cx,cm,g;		    /* A=Z.A */
	for(j=0;j<a;j++) A+=WZ->z[i][j]*WZ->A[j];
	if(0==(A%=M)) continue; 
	g=Egcd(X%M,M,&cx,&cm); assert(0==(A%g));     /* Z_i -= cx * A/g * w_ */
	cx*=A/g; if(cx>0)cx-=M; 
	for(j=0;j<a;j++) {WZ->z[i][j]-=cx*WZ->w[j]; WZ->z[i][j]%=M;}
     }
     X=0;for(j=0;j<a;j++)X+=WZ->w[j]*WZ->A[j];assert(X==WZ->d);/* TEST w.A=d */
     for(j=0;j<a;j++){
       X=WZ->r*WZ->A[j];
       for(k=0;k<d;k++) X+=WZ->B[j][k]*WZ->rI[k];
       assert(X==1);}	      /* TEST:  r*A+B.rI==IP */
     for(i=0;i<WZ->M;i++) {X=0; for(j=0;j<a;j++)
	X+=WZ->z[i][j]*WZ->A[j];assert(0==X%WZ->m[i]);}/* TEST: A*z[]%m[]==0 */
     
     if(FilterFlag) inFILE=NULL;
return 1;
}


#undef	TEST_PP				    /*  print some diagnostic info  */
#define	TEST_PD				   /*  check Poincare duality of PP */

extern FILE *inFILE, *outFILE;

typedef struct {Long x[POINT_Nmax][W_Nmax]; int N, np;}   AmbiPointList;
typedef struct {Long x[POLY_Dmax][W_Nmax]; int N, n;}     AmbiLatticeBasis;

Long Wperm_to_GLZ(Long *W, int *d, Long **G, int *P);
void WeightLatticeBasis(Weight *_w, AmbiLatticeBasis *_B);	
void WeightMakePoints(Weight *_W , AmbiPointList *_P);
int  ChangeToTrianBasis(AmbiPointList*,	AmbiLatticeBasis *, PolyPointList *);

/*  ==========  	  I/O functions:                	==========  */

int  IsDigit(char c) { return (('0'<=c) && (c<='9'));}
void Write_Weight(Weight *W)
{    int n; fprintf(outFILE,"%d",(int) W->d);
     for(n=0; n < W->N; n++) fprintf(outFILE," %d",(int) W->w[n]); 
#if	WZinput
     {  int i,j; if(W->M>0)printf(" "); for(i=0;i<W->M;i++) 
     	{  fprintf(outFILE,"/Z%d: ",W->m[i]);
           for(j=0;j<W->N;j++) fprintf(outFILE,"%d ",W->z[i][j]);
     	} 
     }
#endif
     fprintf(outFILE,"\n");     /* puts("  = d  w_1 ... w_N"); */
}
int  Read_Weight(Weight *_W)     /* read "d w_i" [ or "w_i d" if last=max ] */
{    char c; int i, shift=1, I[W_Nmax+2], FilterFlag=(inFILE==NULL); 
     if(FilterFlag) inFILE=stdin; else
     if(inFILE==stdin) printf("type degree and weights  [d  w1 w2 ...]: "); 
     c=fgetc(inFILE); if( !IsDigit(c) ) return 0; ungetc(c,inFILE); 
     fscanf(inFILE,"%d",I);
     for(i=1; i<W_Nmax+2; i++)
     {	while(' '==(c=fgetc(inFILE))); ungetc(c,inFILE);
	if(IsDigit(c)) fscanf(inFILE,"%d",&I[i]); else break;
     }  while(fgetc(inFILE )-'\n') if(feof(inFILE))return 0; /* read to EOL */
     _W->N=i-1; if(_W->N>W_Nmax) { puts("Increase POLY_Dmax"); exit(0); }
     for(i=0;i<=_W->N;i++) assert(I[i]>0);
     if(I[_W->N]>I[0]) { _W->d=I[_W->N]; shift=0; } else _W->d=I[0];
     for(i=0;i<_W->N;i++) {_W->w[i]=I[i+shift]; assert(_W->w[i]<_W->d);} 
     if(FilterFlag) inFILE=NULL;
     return 1;
}
void Write_WH(Weight *_W, BaHo *_BH, VaHo *_VH, int rc, int tc,
	      PolyPointList *_P, VertexNumList *_V, EqList *_E){
  int i, j;
#if	WZinput
  fprintf(outFILE,"%d ",(int)_W->d);
  for(i=0;i<_W->N;i++) fprintf(outFILE,"%d ",(int)_W->w[i]);
  for(i=0;i<_W->M;i++){fprintf(outFILE,"/Z%d: ",(int)_W->m[i]);
    for(j=0;j<_W->N;j++)fprintf(outFILE,"%d ",(int)_W->z[i][j]);}
#else  
  for(i=0;i<_W->N;i++) fprintf(outFILE,"%d ",(int)_W->w[i]); 
  fprintf(outFILE,"%d=d ",(int)_W->d);
#endif
  fprintf(outFILE,"M:%d %d ",_P->np,_V->nv);  /* PolyData */
  if(rc) fprintf(outFILE,"N:%d %d ",_BH->np,_E->ne); 
  else fprintf(outFILE,"F:%d ",_E->ne); 		 /* END of PolyData */
  if(tc&&(_P->n!=1+_VH->D)) {
    int D=0, d=_W->d;		/* LG non-geometric */
    fprintf(outFILE,"LG: "); for(i=0;i<_W->N;i++) D+=d-2*_W->w[i]; 
    if(D%d) {int g=Fgcd(D,d); fprintf(outFILE,"c/3=%d/%d\n",D/g,d/g);}
    else { int r=_W->N-(D/=d); if((r%2)||(r<=2)) r=0;assert(D==_VH->D);
      for(i=0;i<=D;i++) {fprintf(outFILE,"%sH%d:",i?" ":"",i);
        for(j=0;j<=D-i;j++)fprintf(outFILE,"%s%d",j?",":"",(int)_VH->h[i][j]);}
#if     WZinput
      /* if(r) fputs(_W->R ? " RefGC": " NOrgc",outFILE); */
      if(r) { if(_W->R) fprintf(outFILE," RefI%d",_W->r);
	else fputs(" nonRG",outFILE); }
#endif 
      fputs("\n",outFILE);     }
    return ;
  }
  if((!tc)&&(!rc)){fprintf(outFILE,"non-transversal\n");return;}
  if(_P->n>2){
    if(tc&&rc) for(i=1;i<_P->n-1;i++) if(_VH->h[1][i]!=_BH->h1[i]) 
      { puts("Vafa!=Batyrev"); for(i=1;i<_P->n-1;i++) printf(
	"V[1,%d]=%d  B[1,%d]=%d\n",i,_VH->h[1][i],i,_BH->h1[i]);exit(0);}
    if(tc)      {	
      fprintf(outFILE,"V:%d",_VH->h[1][1]); 
      for(i=2;i<_P->n-1;i++) fprintf(outFILE,",%d",_VH->h[1][i]);  
      if(_P->n > 5){
	for (j=2; 2*j <= _P->n-1; j++){
	  fprintf(outFILE,";%d",_VH->h[j][j]);
	  for(i=j+1;i<_P->n-j;i++) fprintf(outFILE,",%d",_VH->h[j][i]);}}
    }
    else if(rc)     {	
      fprintf(outFILE,"H:%d",_BH->h1[1]); 
      for(i=2;i<_P->n-1;i++) fprintf(outFILE,",%d",_BH->h1[i]);     }
    if(6<_P->n) fprintf(outFILE," [???]\n"); 
    else if(tc||rc) /* Euler number */ {	
      int chi=0, *ho=_BH->h1; if(tc) ho=_VH->h[1];
      if(_P->n==3) chi=4+ho[1];
      if(_P->n==4) chi=2*(ho[1]-ho[2]);
      if(_P->n==5) chi=48+6*(ho[1]-ho[2]+ho[3]);
      if(_P->n==6) chi=24*(ho[1]-ho[2]+ho[3]-ho[4]);
      fprintf(outFILE," [%d]\n",chi);     }	
    else fprintf(outFILE,"\n");} 
  else fprintf(outFILE,"\n");
}

void TEST_LatticeBasis(AmbiLatticeBasis *_B)      /* print AmbiLatticeBasis */
{    int p,a; for(p=0;p<_B->n;p++){
	for(a=0;a<_B->N;a++) printf(" %3d",(int) _B->x[p][a]);
	puts("");     }	
}
#ifdef	TEST
void TEST_WeightMakePoints(AmbiPointList *_P)
{    int i,j; static int MaxPoNum; if(_P->np>MaxPoNum) MaxPoNum = _P->np;
     if(_P->np>20)
     {	for(i=0;i<_P->np;i++) 
	{  for(j=0;j<_P->N;j++) printf("%d ",_P->x[i][j]); puts("");}
     }
     else
     {	for(i=0;i<_P->N;i++) 
	{  for(j=0;j<_P->np;j++) printf("%4d",_P->x[j][i]); puts("");}
     }
     printf("PointNum=%d [max=%d]\n",_P->np,MaxPoNum);
}
#endif

void Ambi_2_Lattice(Long *A,AmbiLatticeBasis *B,Long *P)
{    int i, p=B->n; while(p--) 
     {   int a=p + B->N - B->n; P[p]=A[a] /* -Ai[a]*/ ;
	 for(i=p+1; i<B->n; i++) P[p] -= P[i]*B->x[i][a];
	 if(P[p] % B->x[p][a]) { puts("Error in BasisChange!"); exit(0);}
	 else P[p] /= B->x[p][a];
     } 
}
void Make_Poly_Points(Weight *_W_in, PolyPointList *_PP)
{    AmbiLatticeBasis B; Weight *_W=_W_in; int /* index=0,*/ nip;
     AmbiPointList * _AP = (AmbiPointList *) malloc(sizeof (AmbiPointList));
#ifndef NO_COORD_IMPROVEMENT		/* ==== Perm Coord Improvement ==== */
     Weight Waux=(*_W); _W=&Waux; 
     {	int i,n=_W->N,pi[AMBI_Dmax]; Long M[AMBI_Dmax][AMBI_Dmax],
	    *G[AMBI_Dmax]; for(i=0;i<n;i++) G[i]=M[i];
	Wperm_to_GLZ(_W->w,&n,G,pi);for(i=0;i<n;i++)Waux.w[i]=_W_in->w[pi[i]];
     }
#endif				       /* = End of Perm Coord Improvement = */
     if(_AP==NULL) {printf("Unable to allocate space for _AP\n"); exit(0);}
     WeightLatticeBasis(_W, &B);	/* TEST_LatticeBasis(&B); */
     WeightMakePoints(_W, _AP);	
#ifdef TEST
     puts("\nWeights:");     Write_Weight(_W);  
     puts("AmbiPoints:");    TEST_WeightMakePoints(_AP); 	
     puts("Basis:");	     TEST_LatticeBasis(&B);
     printf("POINT_Num=%d\n",_AP->np);	
#endif
     assert(_AP->np <= POINT_Nmax);   
     nip=ChangeToTrianBasis(_AP, &B, _PP);
#ifdef	TEST
     Print_PPL(_PP,"PolyPoints:");
#endif
#ifdef GEN_CY
     {	int i; for(i=0;i<_W->N;i++) index+=_W->w[i]; if(index%_W->d)index=0;}
     index/=_W->d; if(index>1)
     {  Long A[AMBI_Dmax];int i;for(i=0;i<B.N;i++)A[i]=1-index*_AP->x[nip][i];
	assert(_PP->np<POINT_Nmax); /* need one more to store IP for GEN_CY */
	Ambi_2_Lattice(A,&B,_PP->x[_PP->np]);
     }
#endif
     free(_AP); return;
}
#if	(WZinput) 
int Read_W_PP(Weight *W, PolyPointList *P){ W->P=P; return Read_WZ_PP(W); }
#else
int Read_W_PP(Weight *_W, PolyPointList *_PP){
    if (!Read_Weight(_W)) return 0; Make_Poly_Points(_W,_PP); return 1;
}
#endif

/*  ==========	   WeightLatticeBasis(Weight *_w, AmbiLatticeBasis *_B)
 *   
 *   M[0]:= (-n1, n0, 0, ...) / gcd(n0,n1);
 *   M[i]:= (0,...,0,g/ng,0,...)- (ni/ng) * egcd.Vout(n0,...,n(i-1),0,...);
 *   	    with g=gcd(n0,...,n[i-1]); ng=gcd(g,ni);
 *									    */
#ifdef TEST
void NormTriangularBasis(AmbiLatticeBasis *_B)
{    int p=_B->n-1, a=_B->N-1;	while(p--) 
     {	int pi=_B->n; 	while(0==_B->x[p][--a]);
	while(--pi > p)
	{   int ai=a+1, r = _B->x[pi][a] / _B->x[p][a];
/**/	    if((_B->x[pi][a] % _B->x[p][a]) < 0) r--;
/**/	    if(r) while(ai--) _B->x[pi][ai] -= r * _B->x[p][ai];
	}
     }
}
void OldWeightLatticeBasis(Weight *_w, AmbiLatticeBasis *_B)	
{    int i, j; 
     Long *V, Vin[W_Nmax], Vout[W_Nmax]; 
     for(i=0; i<_w->N; i++) Vin[i]=_w->w[i];
     V=_B->x[0]; *Vout=Fgcd(Vin[0],Vin[1]); 
     V[0] = -Vin[1] / *Vout; V[1] = *Vin / *Vout; _B->N=_w->N; _B->n=_w->N-1;
     for(j=2; j<_w->N; j++) V[j]=0;
     for(i=2; i<_w->N; i++) {
	int g=REgcd(Vin, &i, Vout), gn=Fgcd(g,_w->w[i]); V=_B->x[i-1]; 
        V[i] = g/gn; g = _w->w[i]/gn;
	for(j=0;j<i;j++) V[j]=-Vout[j]*g; 
	j++; while(j<_w->N) V[j++]=0;
     }
     NormTriangularBasis(_B); /* TEST_LatticeBasis(_B); */
}
#endif
void WeightLatticeBasis(Weight *_w, AmbiLatticeBasis *_B)	
{    Long *B[W_Nmax], E[W_Nmax]; int i; _B->N=_w->N; _B->n=_w->N-1;
     B[0]=E; for(i=0;i<_B->n;i++) B[i+1]=_B->x[i]; W_to_GLZ(_w->w,&(_w->N),B);
}

/*  ==========		WeightMakePoints(Weight *_W, AmbiPointList *_P)
 *
 *	MakePoints: 	X[i]-IP[i] = x[j] * B[j][i] :: i-1<=j<B.n,
 *	     		with  \sum X_i w_i=d  and  0<=X[i].
 *									    */
void WeightMakePoints(Weight *_W, AmbiPointList *_P)
{    int i=_W->N-1, X_max[W_Nmax], d_Rest[W_Nmax], X[W_Nmax+1]={0};
     _P->np=0; _P->N=_W->N; d_Rest[i]=_W->d; X_max[i]=d_Rest[i]/_W->w[i];
     while(i<_W->N) 
     {  if(X[i]>X_max[i]) X[++i]++;
	else if(i) 
	{   d_Rest[i-1]=d_Rest[i]-X[i]*_W->w[i];  X[--i]=0; 
	    X_max[i]=d_Rest[i]/_W->w[i];
	}
	else
	{   int j; 
	    Long *Y; 
	    *d_Rest=d_Rest[1]-X[1]*_W->w[1];
	    if( !(*d_Rest % *(_W->w)) ) 
	    {	if(POINT_Nmax <= _P->np){puts("increase POINT_Nmax");exit(0);}
		Y=_P->x[(_P->np)++];
		*Y= *d_Rest/ *(_W->w); for(j=1;j<_W->N;j++) Y[j]=X[j];
	    }
	    X[1]++; i=1;
	}
     }
}

/*  ==========	 ChangeToTrianBasis(AmbiPointList *, Basis *, PolyPointList *)
 *
 *	Assuming that AmbiLatticeBasis is of triangular form, i.e.:
 *	_B->x[p][a] = 0 for a < p + CoDim  ... a < W_Nmax and p < POLY_Dmax
 *	this solves AP=B.PP for AP (to be precise we must take AP.x-IP, 
 *				    where IP is the unique interior point.)
 *	DP=AP.x-IP; DP[a]=PP[p]*B[p][a] => 
 *		PP[p]*B[p][p+1] = DP[p+1] - sum_{a>p+1} PP[p]*B[p][p+1]
 *									    */
int  ChangeToTrianBasis(AmbiPointList *_AP,
			AmbiLatticeBasis *_B, PolyPointList *_PP)
{    int n, ipcount=0, nIP=0;
     if( _AP->N - _B->N) {puts("Dimensions don't match!"); exit(0);}
     else { _PP->n = _B->n; _PP->np = _AP->np; }
     for(n=0; n<_AP->np; n++)
     {	int i, p=1; for(i=0; i < _B->N; i++) if( ! _AP->x[n][i] ) p=0;
	if(p) { ipcount++; nIP=n; }
     }
#ifdef	TEST
     if(ipcount-1) puts("*** Wrong IP-count: take any point as origin! ***");
#endif
     for(n=0; n<_AP->np; n++)
     {	int i, p=_B->n; 
	while(p--) 
	{   int a=p + _B->N - _B->n; _PP->x[n][p]=_AP->x[n][a]-_AP->x[nIP][a];
	    for(i=p+1; i<_B->n; i++) _PP->x[n][p] -= _PP->x[n][i]*_B->x[i][a];
	    if(_PP->x[n][p] % _B->x[p][a]) 
	    { puts("Error in BasisChange!"); exit(0);}
	    else _PP->x[n][p] /= _B->x[p][a];
	}
     }	return nIP;
}
/*  =============		make weights		===========  */
int IfRefWWrite(Weight *W, PolyPointList *P)
{    VertexNumList V; EqList E; Make_Poly_Points(W, P); 
     if(Ref_Check(P,&V,&E)) {Write_Weight(W); fflush(stdout); return 1;} 
     else return 0;
}
void Rec_RefWeights(Weight *W, PolyPointList *P, int g, int sum, int *npp, 
	int *nrp, int n)
{    int wmax=W->d/(W->N-n+1); wmax=min(wmax,W->w[n+1]); 
     wmax=min(wmax,sum-n);
     if(n) for(W->w[n]=wmax;(n+1)*W->w[n]>=sum;W->w[n]--)
	Rec_RefWeights(W,P,Fgcd(g,W->w[n]),sum-W->w[n],npp,nrp,n-1);
     else if(1==Fgcd(g,W->w[0]=sum)) {(*npp)++;if(IfRefWWrite(W,P))(*nrp)++;};
}
void MakeRefWeights(int N, int from_d, int to_d)
{    int npp=0, nrp=0; Weight W; 
     PolyPointList *P = (PolyPointList *) malloc (sizeof(PolyPointList)); 
     W.N=N; assert((N<=W_Nmax)&&(N<POLY_Dmax+2)); assert(P!=NULL);
     for(W.d=from_d;W.d<=to_d;W.d++)
        for(W.w[N-1]=W.d/2; W.d <= N*W.w[N-1]; W.w[N-1]--)
	Rec_RefWeights(&W,P,Fgcd(W.d,W.w[W.N-1]),W.d-W.w[W.N-1],&npp,&nrp,N-2);
     fprintf(outFILE,"#primepartitions=%d #refpolys=%d\n",npp,nrp);exit(0);
}

/*  =============	Landau-Ginzburg-Calculations:		===========  */

void PrintPoCoLi(PoCoLi *P)
{    int i; for(i=0;i<P->n;i++){ 
#if	ABBREV_POLY_PRINT
	int a=ABBREV_POLY_PRINT;
	if((a<=i)&&(i<P->n-a)) {if(i==a) printf("+ ... ");} else
#endif
	{ int co=abs(P->c[i]); printf("%s ", i ? 
	((P->c[i] > 0) ? "+" : "-" ) : "");if(!i&&(P->c[0]<0))printf("- ");
	if(co!=1) printf("%d ",co);
	printf("x^%d ",(int)P->e[i]); }}
     puts(i ? "" : " 0");
}
void UnitPoly(PoCoLi *P) { P->n=1; P->e[0]=0; P->c[0]=1; }
void Add_Mono_2_Poly(int e, Pint c, PoCoLi *P)		   /* use bisection */
{    int m=-1, M=P->n;
     if(P->n) if(e <= P->e[P->n-1]) while(M-m>1)
     {	int pos=(M+m)/2; if(e > P->e[pos]) m=pos; else 
	if (e < P->e[pos]) M=pos; else { P->c[pos]+=c; return; }/* e exists */
     }	if(P->A<=P->n){printf("n=A=%d\n",P->A); fflush(0);
       assert(P->n < P->A);}			/* check #(coeff.) of Poly. */
     for(m=P->n++;M<m;m--) {			/* insert new exponent at M */
	P->c[m]=P->c[m-1]; P->e[m]=P->e[m-1]; } P->e[m]=e;P->c[m]=c;
#ifdef	TEST_PP
     {static int M=1; if(P->n/1000000 > M) {printf("#c=%dM ",++M);fflush(0);}}
#endif
}
void Init1_xN(PoCoLi *P,int N)					 /* 1 - x^N */
{    UnitPoly(P); Add_Mono_2_Poly(N,-1,P);
}
void Remove_Zeros(PoCoLi *AB)
{    int i, s=0; for(i=0; i<AB->n; i++)	       /* s=shift: eliminatet zeros */
     {	if(0==AB->c[i]) s++; else 
	if(s) { AB->c[i-s]=AB->c[i];  AB->e[i-s]=AB->e[i]; }
     }	AB->n -= s;
}
void PolyCopy(PoCoLi *X,PoCoLi *Y)
{    assert(X->n<=Y->A); for(Y->n=0;Y->n<X->n;Y->n++){
	Y->e[Y->n]=X->e[Y->n]; Y->c[Y->n]=X->c[Y->n];}
}
int  BottomUpQuot(PoCoLi *N,PoCoLi *D,PoCoLi *Q,PoCoLi *R)	/* Q*D = N-R */
{    int i, c, Npos, e, E, dD, dN;  Q->n=R->n=0;	      /* return R==0 */
     assert(D->n>0); if(N->n==0) return 1; assert(N->n>0);  /* assume D != 0 */
     dD=D->e[D->n-1]-D->e[0]; dN=N->e[N->n-1]-N->e[0];  
     e=N->e[0]-D->e[0]; if((e<0)||(dD>dN)) { PolyCopy(N,R); return 0; }
     for(Npos=0;Npos<N->n;Npos++)
     {	if(N->e[Npos]<=e+dD)Add_Mono_2_Poly(N->e[Npos],N->c[Npos],R);
	else break;
     }	     E = N->e[N->n-1] - D->e[D->n-1];
     while(0==(R->c[0]%D->c[0]))
     {	Add_Mono_2_Poly(e,c=R->c[0]/D->c[0],Q);		      /* next Q.c[] */
	for(i=0;i<D->n;i++)Add_Mono_2_Poly(e+D->e[i],-c*D->c[i],R); /* R-cQ */
	Remove_Zeros(R); if((Npos==N->n)&&(0==R->n)) return 1; /* finished? */
        e = (R->n) ?  R->e[0]-D->e[0] : N->e[Npos]-D->e[0];   /* next Q.e[] */
     	while(Npos<N->n) { if(N->e[Npos]<=e+dD)
	    {   Add_Mono_2_Poly(N->e[Npos],N->c[Npos],R); Npos++;} else break;
	Remove_Zeros(R); if(e > E) return 0;
     	}
     }	while(Npos<N->n){Add_Mono_2_Poly(N->e[Npos],N->c[Npos],R); Npos++;}
     return 0;
}
void PolyProd(PoCoLi *A,PoCoLi *B,PoCoLi *AB)			/* AB = A*B */
{    int i, j, s=A->n+B->n-1; AB->n=0; assert(s<AB->A);
     for(i=0;i<s;i++)
     {	int m = (i<B->n) ? 0 : i+1-B->n,   M = (i<A->n) ? i+1 : A->n; 
	for(j=m;j<M;j++)
	    Add_Mono_2_Poly(A->e[j]+B->e[i-j],A->c[j]*B->c[i-j],AB);
     }	Remove_Zeros(AB);
}
void Poly_Sum(PoCoLi *A,PoCoLi *B,PoCoLi *S)			/* S = A+B */
{    int a, b=0; S->n=0; for(a=0;a<A->n;a++) { int q=1;
	while(b<B->n) if(A->e[a]<=B->e[b]) break; else {Pint s=B->c[b++];
	    if(s) {if(S->n>=S->A){printf("S.n>%d in S=A+B\n",S->n);exit(0);} 
		S->c[S->n]=s; S->e[S->n++]=B->e[b-1];}}
	if(b<B->n)if(B->e[b]==A->e[a]){Pint s=A->c[a]+B->c[b++]; q=0; 
	    if(s) {assert(S->n<S->A); S->e[S->n]=A->e[a]; S->c[S->n++]=s;} }
	if(q) if(A->c[a]) {assert(S->n<S->A); S->e[S->n]=A->e[a]; 
	    S->c[S->n++]=A->c[a];} }
     while(b<B->n) if(B->c[b]) {assert(S->n<S->A);
	    S->e[S->n]=B->e[b]; S->c[S->n++]=B->c[b++];} else b++; 
}
void Poly_Dif(PoCoLi *A,PoCoLi *B,PoCoLi *D)			/* D = A-B */
{    int a, b=0; D->n=0; for(a=0;a<A->n;a++) { int q=1;
	while(b<B->n) if(A->e[a]<=B->e[b]) break; else {Pint d=-B->c[b++];
	    if(d) {assert(D->n<D->A); D->c[D->n]=d; D->e[D->n++]=B->e[b-1];}}
	if(b<B->n)if(B->e[b]==A->e[a]){Pint d=A->c[a]-B->c[b++]; q=0; 
	    if(d) {assert(D->n<D->A); D->e[D->n]=A->e[a]; D->c[D->n++]=d;} }
	if(q) if(A->c[a]) {assert(D->n<D->A); D->e[D->n]=A->e[a]; 
	    D->c[D->n++]=A->c[a];} }
     while(b<B->n) if(B->c[b]) {assert(D->n<D->A);
	    D->e[D->n]=B->e[b]; D->c[D->n++]=-B->c[b++];} else b++; 
}
void AllocPoCoLi(PoCoLi *P)			/* allocate e[A] and c[A] */
{    assert(0<P->A); assert( NULL != ( P->e = (int *) malloc( P->A * 
	(sizeof(int)+sizeof(Pint) ) ) ));  P->c = (Pint *) & P->e[P->A];
  // printf("AllocPoCoLi: P->A = %d\n", P->A);
}
void Free_PoCoLi(PoCoLi *P) { free(P->e); }	/* free P.e and P.c */
void PoincarePoly(int N, int *w, int d, PoCoLi *P, PoCoLi *Z,PoCoLi *R)
{    int i, e[2]; Pint c[2]; PoCoLi B, *In=Z,*Out=P,*aux;
     if(N==0) {UnitPoly(P); return;} 
     B.e=e; B.c=c; B.A=B.n=2; e[0]=0; c[0]=1; c[1]=-1; Init1_xN(In,d - w[0]);
     for(i=1;i<N;i++)
     {	e[1]=d-w[i]; PolyProd(In,&B,Out); aux=Out; Out=In; In=aux;
     }  /* printf("\nN =");PrintPoCoLi(Z); */
     while(i--)
     {	e[1]=w[i]; assert(BottomUpQuot(In,&B,Out,R)); aux=Out; Out=In; In=aux;
     }	/* printf("Q =");PrintPoCoLi(P); */	assert((R->n)==0); 
#ifdef	TEST_PD
     	{  int M=P->n-1, I=(M+1)/2, E=P->e[P->n-1];
	for(i=0;i<I;i++) {assert(P->e[i]==E-P->e[M-i]);
	assert(P->c[i]==P->c[M-i]); }}
#endif
}


/* ======	I N D E X  &  T R A C E / b01  computation of VaHo	==== */

void Print_VaHo(VaHo *V)
{    int i,j,D=V->D; for(i=0;i<=D;i++){fprintf(outFILE,"H%d*: ",i);
	for(j=0;j<=D;j++) fprintf(outFILE,"%d ",V->h[i][j]);}    
} 
int  DoHodgeTest(VaHo *V)/*[holo] Poincare duality, Hodge duality, sum rule */
{    int i,j,D=V->D,X;for(i=D/2;i<=D;i++){if(V->h[i][0]!=V->h[D-i][0])return 0;
	for(j=0;j<=D;j++)if(V->h[i][j]!=V->h[D-i][D-j]) return 0;}  V->E=X=0;
     for(i=0;i<=D;i++)for(j=i+1;j<=D;j++)if(V->h[i][j]!=V->h[j][i]) return 0;
     for(i=0;i<=D;i++)for(j=0;j<=D;j++){int x=2*i-D; x=x*x*V->h[i][j];
	if((i+j)%2) {V->E -= V->h[i][j];X-=x;} else {V->E += V->h[i][j];X+=x;}}
     assert(3*X==D*V->E);			return (V->h[0][0]==1);
}
int  Hodge_Test(VaHo *V) /* [holo] Poincare duality, Hodge duality, sum rule */
{    if(DoHodgeTest(V)) return 1; else {Print_VaHo(V); return 0;}
}
/*   do {;} while(Multiloop(N,I,j,J));  <-->  do forall 0<=I[j]<N[j], 0<=j<J */
int  Init_Multiloop(int *N, int *I, int *j, int *J)
{    long X=1; for(*j=0;*j<*J;(*j)++) {I[*j]=0;assert(N[*j]>0);} assert(0<=*J);
     if((*J)==0) {N[0]=1; I[0]=0;} 		     /*  safe also for J=0   */
     for(*j=0;*j<*J;(*j)++) X*=N[*j];
     *j=0; return X; /* return order=prod(N) */
}
int  Multiloop(int *N,int *I,int *j,int *J)/* need j=I[]=0 => Init_Multiloop */
{    assert((*j)==0);     if(++(I[0])<N[0]) return 1; 	while(N[*j]<=I[*j])
     {	I[(*j)++]=0; if(*J<=*j) return 0; I[*j]++;} (*j)=0; return 1;
}

int  Count_b01(Weight *W)
{    int b=0,c=0,l,M[POLY_Dmax],I[POLY_Dmax],eq[W_Nmax],j,J=W->M,N=W->N,vac=0;
     for(j=0;j<J;j++) M[j]=W->m[j];
     for(l=0;l<W->d;l++)      /* Z_d==GSO */
     {	c=Init_Multiloop(M,I,&j,&J); do { 	     /* BEGIN generate group */
	int i,can=1, sum=0, k;           /* can(didate) for b01-contribution */
        for(i=0;can&&(i<N);i++){Rat th=rR(l*W->w[i],W->d); /* th=0/q -> eq[] */
        for(k=0;k<J;k++){th=rS(th,rR(I[k]*W->z[k][i],M[k]));} th.N%=th.D;
        assert(th.N>=0); 
        if(th.N) {th=rP(th,rR(W->d,W->w[i])); if(th.N!=th.D) can=0;
          else {eq[i]=1;sum+=W->d-2*W->w[i];}} 
	else eq[i]=0;}					/* eq=(th_i==q_i); */
       	if(sum&&(sum != W->d)) can=0;	            	/* i.e. charge=1? */
        for(k=0;can&&(k<J);k++){int gph=0; /* check inv. under G_k: gph(ase) */
          for(i=0;i<N;i++) {if(eq[i]) gph+=W->z[k][i];  gph%=M[k];}
          if(gph) can=0; }
        if(can) {if(sum) b++; else vac++;}		
#if	(SHOW_b01_TWIST)
	if(can&&sum){ if(b==1)puts("");printf("b01=%d: j^%d, I[]=",b,l);
	for(k=0;k<J;k++)printf(" %d",I[k]);puts("");}
/* 	printf("M[%d]=",J);for(i=0;i<J;i++)printf("%d ",M[i]);  	     */
/* 	printf(" I=");for(i=0;i<J;i++)printf("%d ",I[i]);printf(" j=%d\n",j);*/
#endif
     	assert(0<(c--));} while(Multiloop(M,I,&j,&J)); assert(c==0);  /* END */
     }	assert(vac==1); return b;
}
void Fast_c9_VaHo(Weight *W,VaHo *V)/* -->> V.D<=3 <<--  via index and trace */
{    int i,j,k,l,c, b01=0;        /* ... from "proced" in lgotwist.c */
     int fac, wort, nvar, expo, ns=W->M, N=W->N, R, mask[W_Nmax+1], n=N;
     Long zsum1=0, zsum2=0, mo=W->d, ng, ngb, kgV=mo, omega[POLY_Dmax]; 
     int *wo[3],*woG,*woA,*woS, WM=1<<W_Nmax,I[POLY_Dmax]; Rat prod; V->D=0; 
     for(i=0;i<W->N;i++) V->D+= W->d-2*W->w[i]; 
     assert(V->D%W->d==0); V->D/=W->d; assert((0<V->D)&&(V->D<=3));
     V->h[0][0]=V->h[0][V->D]=V->h[V->D][0]=V->h[V->D][V->D]=1;
     if(V->D==1){assert(Count_b01(W)==1); return;}
     if(V->D==2){b01=Count_b01(W);
	V->h[0][1]=V->h[1][0]=V->h[2][1]=V->h[1][2]=b01;
	assert((b01==0)||(b01==2)); V->h[1][1]=(b01) ? 4 : 20; return;}
     wo[0] = woG = (int *) malloc(WM*3*sizeof(int)); 
     assert(woG!=NULL); wo[1]=woS=&woG[WM]; wo[2]=woA=&woS[WM];
     for(i=0;i<=N;i++) mask[i]=1<<i;
     for(i=0;i<WM;i++) woG[i]=woS[i]=woA[i]=0;
     for(k=0;k<ns;k++) {kgV *= W->m[k]/Fgcd(kgV,W->m[k]); mo*=W->m[k];}
     for(k=0;k<ns;k++) omega[k]=kgV/W->m[k];
     R=kgV/W->d;
     for(l=0;l<W->d;l++)					 /* Z_d==GSO */
     {  c=Init_Multiloop(W->m,I,&j,&ns); do {        /* BEGIN generate group */
	nvar=0; for(i=0;i<N;i++) { expo=l*W->w[i]*R; for(k=0;k<ns;k++)expo+=
	  I[k]*omega[k]*W->z[k][i];if(!(expo%kgV))nvar+=mask[i];}woS[nvar]++;
        assert(0<(c--));}while(Multiloop(W->m,I,&j,&ns));assert(c==0);/* END */
     }  assert(woS[mask[n]-1]==1);	 /* over=overcount=wo[mask[n]-1][1]; */
     for(i=0;i<mask[n];i++) { fac = 1;
       for(j=0;j<n;j++) if ((mask[j] & i) == mask[j])  fac*=-1;
       if(woS[i]) for(k=0;k<mask[n];k++) if(woS[k]){ /* woS = sum of contrib */
            wort = i & k;
            wo[fac+1][wort]+=woS[i]*woS[k]; }  }
     for(i=mask[n]-1;i>=0;i--) if (woG[i]+woA[i]) {
       prod=rI(1);
       for (j=0;j<n;j++) if ((mask[j]&i)==mask[j])
          prod = rP(prod,rR(W->w[j]-W->d,W->w[j]));
       if (prod.D != 1){ fprintf(outFILE,
		"\nDenominator != 1 in Fast_c9_VaHo (LG.c)\n");exit(0);}
       zsum1+= woG[i]*prod.N; zsum2+= woA[i]*prod.N;}
     ng = -(zsum1/mo+2)/2; ngb = (zsum2/mo-2)/2;		    /* /over */
     assert((zsum1==-2*(ng+1)*mo) && (zsum2 == 2*(ngb+1)*mo));	    /* *over */
     free(woG);		/* woS[word] = #group elements :: survivors==word */
			/* woG / woA = contributions to ng / ngb */
     if(ng==ngb){b01=Count_b01(W); assert((b01==0)||(b01==1)||(b01==3));}
     for(i=1;i<=2;i++){
	V->h[i][i]=ngb-2*b01; V->h[i][3-i]=ng-2*b01;
	V->h[0][i]=V->h[3][i]=b01; V->h[i][0]=V->h[i][3]=b01;
	V->h[0][0]=V->h[3][3]=V->h[3][0]=V->h[0][3]=1;
     }
}

int  WIndex_HTrace(Weight *W, int *WI, int *T)/* T=sum(Hij), return over=H00 */
{    int i,j,k,l,c,vacnum;  /* ... from "proced" in lgotwist.c */
     int fac, wort, nvar, expo, ns=W->M, N=W->N, R, mask[W_Nmax+1], n=N;
     Long zsum1=0, zsum2=0, mo=W->d, kgV=mo, omega[POLY_Dmax]; 
     int *wo[3],*woG,*woA,*woS, WM=1<<W_Nmax,I[POLY_Dmax]; Rat prod;
     wo[0] = woG = (int *) malloc(WM*3*sizeof(int)); 
     assert(woG!=NULL); wo[1]=woS=&woG[WM]; wo[2]=woA=&woS[WM];
     for(i=0;i<=N;i++) mask[i]=1<<i;
     for(i=0;i<WM;i++) woG[i]=woS[i]=woA[i]=0;
     for(k=0;k<ns;k++) {kgV *= W->m[k]/Fgcd(kgV,W->m[k]); mo*=W->m[k];}
     for(k=0;k<ns;k++) omega[k]=kgV/W->m[k];
     R=kgV/W->d;
     for(l=0;l<W->d;l++)					 /* Z_d==GSO */
     {  c=Init_Multiloop(W->m,I,&j,&ns); do {        /* BEGIN generate group */
	nvar=0; for(i=0;i<N;i++) { expo=l*W->w[i]*R; for(k=0;k<ns;k++)expo+=
	  I[k]*omega[k]*W->z[k][i];if(!(expo%kgV))nvar+=mask[i];}woS[nvar]++;
        assert(0<(c--));}while(Multiloop(W->m,I,&j,&ns));assert(c==0);/* END */
     }  vacnum=woS[mask[n]-1];		 /* over=overcount=wo[mask[n]-1][1]; */
     for(i=0;i<mask[n];i++) { fac = 1;
       for(j=0;j<n;j++) if ((mask[j] & i) == mask[j])  fac*=-1;
       if(woS[i]) for(k=0;k<mask[n];k++) if(woS[k]){ /* woS = sum of contrib */
            wort = i & k;
            wo[fac+1][wort]+=woS[i]*woS[k]; }  }
     for(i=mask[n]-1;i>=0;i--) if (woG[i]+woA[i]) {
       prod=rI(1);
       for (j=0;j<n;j++) if ((mask[j]&i)==mask[j])
          prod = rP(prod,rR(W->w[j]-W->d,W->w[j]));
       if (prod.D != 1){ fprintf(outFILE,
		"\nDenominator != 1 in Fast_c9_VaHo (LG.c)\n");exit(0);}
       zsum1+= woG[i]*prod.N; zsum2+= woA[i]*prod.N;}
     free(woG);		/* woS[word] = #group elements :: survivors==word */
			/* woG / woA = contributions to ng / ngb */
     assert(zsum1%mo==0); assert(zsum2%mo==0); assert(vacnum==1);
      *WI=(zsum2+zsum1)/mo; *T=(zsum2-zsum1)/mo; return vacnum;
}
/*   Z_d is taken care of by the integral charge condition.		     *
 *   x[j]=lcm(phase denom.), go over G=\prod M[j], denom=\prod(1-t^(w*x))    *
 *   numerator: go over group G and project: this is a critical part that    *
 *   can have extremely many terms			     *
 *   U=charge unit for th[i], X=phase unit for (non-GSO) group projection    *
 *   QL=sum(th_i-q_i)+Q_ring; dQ=QR-QL=sum(1-2th_i);			     *
 * 	proj: g|h> = (-1)^(Kg*Kh) \e(g,h)(\det g_h)/(det g)|h>, thus	     *
 *	K=\e=0 => ph[j] is the sum of the phases on invariant fields.	     *
 */	
int Test_BottomUpQuot(PoCoLi *Num,PoCoLi *Den,PoCoLi *Quo,PoCoLi *Rem)
{   int i=BottomUpQuot(Num,Den,Quo,Rem); if(Rem->n)assert(i==0);
    if(i)return 1;
    printf("Num="); PrintPoCoLi(Num);
    printf("Den="); PrintPoCoLi(Den); printf("Quo=");PrintPoCoLi(Quo);
    printf("Rem="); PrintPoCoLi(Rem);
/* printf("N=%d D=%d Q=%d R=%d\n",Num.A,Den.A,Quo.A,Rem.A);exit(0);fflush(0);*/
    {	PoCoLi A; A.A=10000;AllocPoCoLi(&A);PolyProd(Den,Quo,&A);
	Poly_Dif(&A,Num,Den);Poly_Sum(Den,Rem,&A); assert(A.n==0);
	Free_PoCoLi(&A);	puts("BottomUpQuot: Test o.k.");
    }	return 0;
}

void pff(char *c){puts(c);fflush(0);}

typedef	struct {int X, n, *d, **mt;} /* mt[i][j]=mobius(j,i) */	MobiusData;

void MakeMobius(MobiusData *M,int X)   /* d[i] divisors, mt[i][j] 0<=j<=i<=n */
{    int i,j=0,n=0,*d; 
     for(i=1;i<X/i;i++) if(X%i==0)n++;
     M->n=n=2*n+(X==i*i);	/* #(Div(X)) */
     M->d = d = (int *)malloc(((n*(n+3))/2) * sizeof(int) + n * sizeof(int*));
     assert(d != NULL); M->mt = (int **) &d[(n*(n+3))/2]; M->mt[0]=&(d[n]);
     for(i=1;i<n;i++)M->mt[i]=&M->mt[i-1][i];
     for(i=1;i<=X/i;i++) if(X%i==0){d[j]=i;d[n-(++j)]=X/i;}
     for(i=0;i<n;i++){int k; M->mt[i][i]=1;
	for(j=i-1;0<=j;j--) {M->mt[i][j]=0; for(k=j+1;k<=i;k++) 
	    if(d[i]%d[k]==0) if(d[k]%d[j]==0) M->mt[i][j] -= M->mt[i][k];}}
#ifdef	PRINT_MOBIUS_FUNCTION
	printf("Div(%d)=",X);for(i=0;i<n;i++)printf("%d ",d[i]);
	for(i=0;i<n;i++){printf("m[%d]=",d[i]);
	for(j=0;j<=i;j++)printf("%d ",M->mt[i][j]);} ;puts("");
#endif
}
void FreeMobius(MobiusData *M){free(M->d);}
void Calc_VaHo(Weight *W,VaHo *V);
void PoincarePoly(int N, int *w, int d, PoCoLi *P, PoCoLi *Z,PoCoLi *R);
void Aux_Phase_Poly(PoCoLi *P,int w,int d, int r, int s, int x);
int  Index_Trace_Test(VaHo *V,int WI,int T)
{    int i,j,D=V->D;for(i=0;i<=D;i++)for(j=0;j<=D;j++){T-=V->h[i][j];
	if((i+j)%2) WI += V->h[i][j]; else WI -= V->h[i][j];}
	assert((T==0)&&(WI==0)); return 1;
}
/* Thru 4-folds it is sufficient to know boundary (i.e. H{0i} and H{di} *)   *
 * contributions, Witten index and trace for the twisted sectors; the	     *
 * untwisted sector can be reconstructed from "WIndex_HTrace(W,&WI,&T)".     *
 * goint twice over the group is not too costly if the projection is reduced *
 * to the effectively acting subgroup for each twisted sector		     */
void LGO_VaHo(Weight *W,VaHo *V)	      
{    int i,d=W->d,D=0; for(i=0;i<W->N;i++) D+=d-2*W->w[i];assert(!(D%d));D/=d;
  /* printf("In LGO_VaHo: W->N: %d W->M: %d D: %d\n", W->N, W->M, D); */
  /* From palp-2.0 to palp-2.11 the next lines were
   * if((V->sts)||(D>3)) {if(W->M==0) W->m[0]=1;} 
   * else {if(D<=3) {Fast_c9_VaHo(W,V);return;}
   *   if(W->M==0){Calc_VaHo(W,V);return;}}
   * but that makes no sense since it's equivalent to
   * if((V->sts)||(D>3)) {if(W->M==0) W->m[0]=1;} 
   * else {Fast_c9_VaHo(W,V);return;} 
   * therefore:                                        */
     if ((V->sts == 0) && (W->M == 0)) {  /* V->sts ... lg-flag set to 2,
                                            W->M ... # of Z_n symmetries */
       if (D<=3) {Fast_c9_VaHo(W,V);return;}
       else {Calc_VaHo(W,V);return;} }
     if(W->M==0) W->m[0]=1;
     {	Long dQ,QL,q; PoCoLi *S,*SO,*SN;/* dQ=sum(1-2th)=QR-QL; QL=sum(th-w) */
	int j,k,s, J=W->M, M[POLY_Dmax],I[POLY_Dmax],X=W->m[0],x[POLY_Dmax];
	int WI,T; Long ph[POLY_Dmax], th[W_Nmax], U=X*d/Fgcd(X,d), rd=U/d, G=1;
	PoCoLi Den,Num,Quo,Rem; MobiusData MX; MakeMobius(&MX,X);
	Den.A = (1<<W->N);	AllocPoCoLi(&Den);
	Quo.A = X*d*D+2*W->N; 	AllocPoCoLi(&Quo);
	Num.A=Rem.A= Quo.A;	AllocPoCoLi(&Num); AllocPoCoLi(&Rem);
	V->D=D; for(i=0;i<=D;i++) for(j=0;j<=D;j++) V->h[i][j]=0;
	for(j=0;j<J;j++) {G*=(M[j]=W->m[j]); assert(*M%M[j]==0); x[j]=X/M[j];}
	S = (PoCoLi *) malloc( 2*X*sizeof(PoCoLi) ); assert(S!=NULL); 
  	for(k=0;k<W->d;k++)					 /* Z_d==GSO */
     	{   int v,a=Init_Multiloop(W->m,I,&v,&J); do {	  /*BEGIN gen TWISTS */
		int n[W_Nmax], N=0,h[POLY_Dmax+1],hn=0;
		QL=dQ=0; for(j=0;j<J;j++) ph[j]=0; for(i=0;i<=D;i++) h[i]=0;
		for(i=0;i<W->N;i++) {th[i]=k*W->w[i]*rd; 
		    for(j=0;j<J;j++) th[i]+=I[j]*W->z[j][i]*(U/W->m[j]);
		    assert(th[i]>=0); th[i]%=U; 
		    if(th[i]) {QL+=th[i]-W->w[i]*rd; dQ+=U-2*th[i];}
		    else{n[N++]=i;for(j=0;j<J;j++)ph[j]+=x[j]*W->z[j][i];}
	    	}   assert((dQ%U)==0); assert((QL%rd)==0); assert(QL+dQ>=0);
	    	dQ/=U; QL/=rd; for(j=0;j<J;j++)ph[j]%=X; /* Q & ph of |h>_gs */

		/* int pntw=0; if(100*k+10*I[0]+I[1]==202)pntw=1; */

	/* U==lcm(X,d)=d*rd -> th[i]/U, X==lcm(M[])=x[j]*M[j] -> ph[j]/X |g> */
	/* QL/d+dQ	G=\prod M[j]=|group|   n[i], i<N = untw.fields	     */
	/* twist= gso^k*g_j^I[j]     r/oz=phase(\Ph_n)  s/oz=phase of vacuum */
		if(N==0) {i=1; for(j=0;j<J;j++) if(ph[j]) i=0;    
		     if(i&&(QL%d==0)) hn=h[QL/d]=1; }/* trivial ring sectors */
	/*  s1+l r1 \in X\IZ, g=ar+bX => l=(l1 X - a1 s1)/g1,		     */
	/*  A g1+B g2=g, if((s1 a1 g2/g - s2 a2 g1/g) % X ==0) then	     */
	/*  l=l' X/g +(B X-a1 s1)/g1 = l' X/g - (A X +a2 s2)/g2 else break;  */
	/*  r'=g, s'=(A X + a2 s2)/g2, a'=1;  assert(s' % g' ==0);	     */
		if(N==1) {int w=W->w[*n], r=0,g=X; Long A,B; s=0; 
		    for(j=0;j<J;j++) {			/* cyclic projection */
			int rj=(W->z[j][*n]*x[j])%X, sj=ph[j]%X; assert(0<=rj);
			if(rj==0) { if(sj==0) continue; else {g=0;break;} }
			g=Egcd(rj,X,&A,&B); sj=(sj*A)%X; 
				assert(rj>0);assert(g>0);assert(sj%g==0);
			if(sj%g) {g=0;break;}
			rj=g; g=Egcd(r,rj,&A,&B); assert(g>0);
			if((s*rj-sj*r)% (g*X)) {g=0;break;}
			r=g; s=(A*X+sj)/rj; if(s%g) {g=0;break;}
			}	assert((j<J)==(g==0));
		    if(g==0) Quo.n=0; else if(g==X) Init1_xN(&Num,d-w);   
		    else {/* assert(g>1); */ Aux_Phase_Poly(&Num,w,d,r,s,X/g);}
		    if(g) {Init1_xN(&Den,w*(X/g));
		    	BottomUpQuot(&Num,&Den,&Quo,&Rem);
		    	for(i=0;i<Quo.n;i++) if(( q=(QL+Quo.e[i]) )%d==0)
			    {h[q/d]+=Quo.c[i]; hn=1;}
		    }}
	if(N>1)		/* N>1: S[] sector -> SN[] new sector PoCoLi */
	{  int u,b,l,L[POLY_Dmax],w[W_Nmax],o[W_Nmax],io[W_Nmax],mm=1;
	for(i=0;i<N;i++) for(j=i+1;j<N;j++)/* sort invariant weights */
	    if(W->w[n[j]]<W->w[n[i]]) swap(&n[i],&n[j]); 
	for(i=0;i<N;i++) {w[i]=W->w[n[i]]; o[i]=1; 
	    for(j=0;j<J;j++){ u=NNgcd(W->z[j][n[i]],W->m[j]); u=W->m[j]/u;
		o[i]=u*o[i]/Fgcd(o[i],u);} mm=mm*o[i]/Fgcd(mm,o[i]);}
	/* phase(g.s.) = gs/mm; mm=lcm(o[i]); phase(X(n(i))=z[i]/o[i];     */

						  /*  BEGIN GROUP PROJECTION */
	    if(mm>1){int *mo, m, ego=1; /* effective group order */
		PoCoLi *A,*B,*C,Ax; Ax.A=Quo.A; AllocPoCoLi(&Ax); Num.n=0;
		for(j=0;j<J;j++){int g=W->z[j][n[0]]; 
		    for(i=1;i<N;i++) g=NNgcd(g,W->z[j][n[i]]);
		    M[j] = W->m[j] / NNgcd(g,W->m[j]);}
		for(m=0;m<MX.n;m++)if(mm<=MX.d[m])break;assert(mm==MX.d[m]);
	    	mo=MX.mt[m]; l=0; for(j=0;j<m;j++) 
		if(mo[j]==1)l++; else if(mo[j]==-1)l++; else assert(mo[j]==0);

		b=ego=Init_Multiloop(M,L,&u,&J); do {	 /* BEGIN make group */
		int gs=0,z[W_Nmax];for(i=0;i<N;i++){z[i]=0;   /* g.s. phases */
		    for(j=0;j<J;j++) z[i]+=L[j]*((W->z[j][n[i]]*o[i])/W->m[j]);
			io[i]=mm/o[i]; gs+=io[i]*z[i];} gs%=mm; assert(gs>=0);
		SO=S; SN=&S[mm];	/* for each projection group element */
		for(s=0;s<mm;s++){int ds=gs+s; S[s].A=2*mm;	  /* init SO */
		    AllocPoCoLi(&S[s]); if( ds % io[0] ) S[s].n=0;
		    else if(*o>1)Aux_Phase_Poly(&S[s],w[0],d,z[0],-ds/io[0],
			 o[0]);  else Init1_xN(&S[s],d-w[0]);}
		for(i=1;i<N;i++){PoCoLi *Saux=SO,*ac; int t; for(t=0;t<mm;t++){
		    A=&Ax; B=&Quo; C=&Rem; C->n=0;
		    for(s=0;s<mm;s++)if(s%io[i]==0){
			if(o[i]>1) Aux_Phase_Poly(A,w[i],d,z[i],-s/io[i],o[i]);
			else Init1_xN(A,d-w[i]);
			PolyProd(&SO[(mm+t-s)%mm],A,B); Poly_Sum(B,C,A);
			ac=A;A=C;C=ac;}
		    SN[t].A=C->A; AllocPoCoLi(&SN[t]); PolyCopy(C,&SN[t]);
		    }	for(t=0;t<mm;t++) Free_PoCoLi(&SO[t]);
		    SO=SN; SN=Saux;
		}  			/* S[s] finished */

		/* if(pntw)for(s=0;s<mm;s++){
		printf("gs=%d S[%d]= ",gs,s);PrintPoCoLi(&SO[s]);} */

		Poly_Sum(&SO[0],&Num,&Rem); B=&Num; A=&Rem;
		for(j=0;j<m;j++) if(mo[j]==1) {Poly_Sum(A,&SO[MX.d[j]],B);
		    C=A;A=B;B=C;} else if(mo[j]==-1)
		    {Poly_Dif(A,&SO[MX.d[j]],B); C=A;A=B;B=C;}
		assert((A==&Num)==(l%2)); 	if(l%2==0) PolyCopy(A,&Num);
		for(s=0;s<mm;s++){Free_PoCoLi(&SO[s]);}
		assert(0<(b--));} while(Multiloop(M,L,&u,&J)); assert(b==0);

		for(i=0;i<Num.n;i++){assert(0==Num.c[i]%ego); Num.c[i]/=ego;}

		/* if(pntw){printf("Num[%d:%d,%d]=",k,I[0],I[1]); 
		   PrintPoCoLi(&Num);fflush(0);} */

		A=&Num; B=&Quo;
		if(Num.n) for(i=0;i<N;i++){Init1_xN(&Den,w[i]*o[i]);
		    BottomUpQuot(A,&Den,B,&Rem); C=A;A=B;B=C;}
		for(i=0;i<A->n;i++) if(( q=(QL+A->e[i]) )%d==0)
			    {h[q/d]+=A->c[i]; hn=1;}
		
	/* if(pntw){printf("Quo[%d:%d,%d]=",k,I[0],I[1]); 
		PrintPoCoLi(A);fflush(0);} */

	    	Free_PoCoLi(&Ax);		    /*  END GROUP PROJECTION */
	    }
	    else
	    { 	PoincarePoly(N,w,d,&Quo,&Num,&Den);
	    	for(i=0;i<Quo.n;i++) if(( q=(QL+Quo.e[i]) )%d==0)
		    {h[q/d]+=Quo.c[i]; hn=1;} /*END of no projection */
	    }
	}						/* N>1 CASE FINISHED */

		if(hn) for(i=0;i<=D;i++) if(h[i]) V->h[D-i][i+dQ]+=h[i];
	if(V->sts)	if(hn){fprintf(outFILE,"sec[%d",k);
	    	    for(j=0;j<J;j++) fprintf(outFILE,"%s%d",j?",":":",I[j]); 
	    	    fputs("]",outFILE);
	fprintf(outFILE," th=%2ld",th[0]);for(i=1;i<W->N;i++)
	fprintf(outFILE," %2ld",th[i]); /*fprintf(outFILE,"/%d ",U);*/ 
	/*fprintf(outFILE," %d*ph=",U);for(i=0;i<J;i++)printf("%d ",ph[i]); */
	fprintf(outFILE,"  QL=%2ld/%d dQ=%2ld ",QL,d,dQ);
	/*fprintf(outFILE,"N=%d ",N);*/
	    	for(i=0;i<=D;i++) if(h[i]) fprintf(outFILE," q%d%ld+=%d",
		    i,i+dQ,h[i]); fputs("\n",outFILE);}
							/* if(!cont)exit(0); */
	    assert(0<(a--));
	    } while(Multiloop(W->m,I,&v,&J));assert(a==0); /* END gen TWISTS */
     	}   free(S); FreeMobius(&MX);	Free_PoCoLi(&Rem); Free_PoCoLi(&Quo);
	assert(Hodge_Test(V)); Free_PoCoLi(&Den); Free_PoCoLi(&Num);
	assert(1==WIndex_HTrace(W,&WI,&T)); assert(Index_Trace_Test(V,WI,T));
     	if(V->sts) printf("WittenIndex=%d, Trace=%d\n",WI,T);
     }
}
/* T=\x*t: (1-T^(d-w))/(1-T^w)=(1-T^(d-w))*(1+T^w+T^2w+...+T^(x-1)w)/(1-T^wx)*
 * phase(\x)=r/x. All terms in numerator with phase s/x, p=0,...,x-1 =>	     *
 * exp=p mod x in numerator  1+...+T^(x-1)w-T^(d-w)-T^d-...-T^(d+w(x-2))     */
void Aux_Phase_Poly(PoCoLi *P,int w,int d, int r, int s, int x)
{    int i, nocancel = (d%w) || (x*w<d) || ((r*d/w)%x), 
	 negmin = nocancel ? d-w : x*w,	negmax=d+x*w-2*w, 
	 posmax = nocancel ? (x-1)*w : d-2*w;	assert(x>1); P->n=0;
     for(i=0;i<=posmax;i+=w) if((r*(i/w)-s) % x == 0) Add_Mono_2_Poly(i,1,P);
     for(i=negmin;i<=negmax;i+=w)if((r*(i-d)/w-s)%x==0)Add_Mono_2_Poly(i,-1,P);
}

void Calc_VaHo(Weight *W,VaHo *V) {    
     int i, j, k, D=0, w[W_Nmax], N=W->N, d=W->d;
     PoCoLi A,B,C, *Z=&A, *R=&B, *P=&C; 
     for(i=0;i<N;i++) { w[i]=W->w[i]; D+=d-2*w[i]; }
     if(D % d) {V->D=0; return;} else V->D = (D/=d);
     for(i=0;i<N;i++)for(j=i+1;j<N;j++)if(w[j]<w[i]) swap(&w[i],&w[j]);
#ifndef	COEFF_Nmax
#define	COEFF_Nmax	(d*D+2*N)		/* 22999000 for W_Nmax == 6 */
#endif
     P->A=Z->A=COEFF_Nmax; R->A=w[N-1]+1; //printf("%d %d %d", P->A, Z->A, R->A); 	
     AllocPoCoLi(P); AllocPoCoLi(Z); AllocPoCoLi(R);  
     PoincarePoly(N,w,d,P,Z,R);		
     assert(D*d==P->e[P->n-1]);
     {  int n=P->n, cM=0; long long sum=0,num=1,den=1; /* check sum = P(0) */
	for(i=0;i<n;i++) {
	  Pint co=P->c[i]; assert(co>0); if(co>cM) cM=co;
	  sum+=co; }
	for(i=0;i<W->N;i++){	num*=W->d-W->w[i]; den*=W->w[i];
	   {long long g=LFgcd(num,den); num/=g; den/=g; }} 
#ifdef	TEST_PP
	if(P->n<99) {printf("PP =");PrintPoCoLi(P);} else printf(
	   "#(Exp,Co)=%d  Exp<=%d  Coeff<=%d  sum=%lld\n",n,P->e[n-1],cM,sum);
#endif
	if(den==1) assert(num==sum);
	else { printf("sum=%lld  test=%lld/%lld \n",sum,num,den);}
	}
     for(i=0;i<=D;i++) for(j=0;j<=D;j++) V->h[i][j] = 0;
     for(i=0;i<P->n;i++)if(P->e[i] % d==0) {j=P->e[i]/d; V->h[D-j][j]=P->c[i];}
#ifdef	TEST_PP
     for(i=0;i<=D;i++){for(j=0;j<=D;j++)printf("%6d ",(int)V->h[i][j]);
	puts("= V-pri");}fflush(0);
#endif
     for(k=1;k<d;k++)					/* k-twisted sector */
     {	int th[W_Nmax], Tmt, et=0, eT=0, Deff=0, n=0;
	for(i=0;i<N;i++)
	{   if((th[i]=((long long) k * W->w[i]) % (long long)d))
	    {	et+=th[i]-W->w[i];eT+=W->d-th[i]-W->w[i];}
	    else { Deff+=W->d-2*W->w[i]; w[n++]=W->w[i];}    /* effective w */
	}
	if(0 == ((Tmt=eT-et) % d))
	{   for(i=0;i<n;i++)for(j=i+1;j<n;j++)if(w[j]<w[i]) swap(&w[j],&w[j]);
	    PoincarePoly(n,w,d,P,Z,R); assert(Deff==P->e[P->n-1]); Tmt/=d;
	    for(i=0;i<P->n;i++) 
	    {	j=et+P->e[i]; 
		if(j % d==0) {j/=d; if(j+Tmt>=0)V->h[D-j][j+Tmt]+=P->c[i];}
	    }
	}
     }	Free_PoCoLi(P); Free_PoCoLi(Z); Free_PoCoLi(R); 
#ifdef	TEST_PP
     for(i=0;i<=D;i++){for(j=0;j<=D;j++)printf("%6d ",V->h[i][j]);puts("= V");}
#endif
}


/*  =============	     Trans_Check(Weight W)		==========  */
/*  OLD version ignores
 *  orbifold projection */

int OLDsymcheck(int sum, int link, int *mon, Weight W, int *mask){
/* symcheck returns 1 if there is a monomial mon in
 * the variables indicated by link whose total weight is sum and which
 * transforms under the k'th symmetry with a phase sum[k];
 * if X_i occurs in the monomial then mon[i] is set to 1.                   */
   int i, j, newsum, newlink;

   for (i=0;i<W.N;i++) if (link&mask[i]){
      if (!(sum % W.w[i])){
        *mon=*mon|mask[i]; 
        return 1;};   
      newlink=link-mask[i];
      for (j=0;j*W.w[i]<sum;j++){
        newsum=sum-j*W.w[i];
        if (OLDsymcheck(newsum,newlink,mon,W,mask)) {
          if (j) *mon=*mon|mask[i];
          return 1; };   };   };
   return 0;
}
void OLDInit_Trans_Check(int *mask, int **targets, int **mighty)
{    int i, j, maxPN=1<<W_Nmax; mask[0]=1;		  /* maxPN=2^W_Nmax */
     for(i=1;i<=W_Nmax;i++)mask[i]=2*mask[i-1]; /* mask={1,2,4,8,16,32,...} */

     assert(maxPN==mask[W_Nmax]);/* maximum number of pointers at one point */
     assert( NULL != (*targets = (int *) malloc(maxPN * sizeof(int))) );
     assert( NULL != (*mighty  = (int *) malloc(maxPN * sizeof(int))) );

     for (i=0;i<W_Nmax;i++) for (j=mask[i];j<mask[i+1];j++)/* make mighty[] */
       (*mighty)[j]=(*mighty)[j-mask[i]]+1;		       
}

int OLDTrans_Check(Weight W){/* returns 1 if non-degenerate potential exists */
   int i, j, k, mon;   			/* j represents the link!!!        */
   static int mask[1+W_Nmax], *targets, *mighty; assert(W.N<=W_Nmax);
   if(targets==NULL) OLDInit_Trans_Check(mask,&targets,&mighty);
   targets[0]=0;
   for (i=0;i<W.N;i++) for (j=mask[i];j<mask[i+1];j++){
     targets[j]=0;
     for (k=0;k<=i;k++) targets[j]=targets[j]|targets[j&~mask[k]];
     for (k=0;(k<W.N)&&(mighty[targets[j]]<mighty[j]);k++) 
       if (!(mask[k]&targets[j])) if (OLDsymcheck(W.d-W.w[k],j,&mon,W,mask)) 
         targets[j]=targets[j]|mask[k];
     if (mighty[targets[j]]<mighty[j]) return 0;}
   return 1;
}

typedef struct {int d, m[POLY_Dmax];}            symlist;

int symcheck(symlist sum, int link, Weight W, int *mask){
/* symcheck returns 1 if there is a monomial in
 * the variables indicated by link whose total weight is sum.d and which
 * transforms under the k'th symmetry with a phase sum.m[k];    */

   int i, j, k, check, expo, newlink;
   symlist newsum;

   for (i=0;i<W.N;i++) if (link&mask[i]){
     if (!(sum.d % W.w[i])){
       expo=sum.d/W.w[i];
       check=1;
       for (j=0;(j<W.M)&&check;j++)
         if ((expo*W.z[j][i]-sum.m[j])% W.m[j]) check=0;
       if (check) return 1;}
     newlink=link-mask[i];
     for (j=0;j*W.w[i]<sum.d;j++){
       newsum.d=sum.d-j*W.w[i];
       for (k=0;k<W.M;k++) newsum.m[k]=sum.m[k]-j*W.z[k][i];
       if (symcheck(newsum,newlink,W,mask)) return 1;  }   }
   return 0;
}

void Init_Trans_Check(int *mask, int **targets, int **mighty)
{    int i, j, maxPN=1<<W_Nmax; mask[0]=1;                /* maxPN=2^W_Nmax */
     for(i=1;i<=W_Nmax;i++)mask[i]=2*mask[i-1]; /* mask={1,2,4,8,16,32,...} */

     assert(maxPN==mask[W_Nmax]);/* maximum number of pointers at one point */
     assert( NULL != (*targets = (int *) malloc(maxPN * sizeof(int))) );
     assert( NULL != (*mighty  = (int *) malloc(maxPN * sizeof(int))) );

     for (i=0;i<W_Nmax;i++) for (j=mask[i];j<mask[i+1];j++)/* make mighty[] */
       (*mighty)[j]=(*mighty)[j-mask[i]]+1;
}

/*   targets[j] bin-encodes the gradients that can be nonzero with variables
     bin-encoded in j;  transversality <==> mighty[j] <= mighty[targets[j]] */

int Trans_Check(Weight W){  /* returns 1 if non-degenerate potential exists */
   int i, j, k, l;                      /* j represents the link!!!        */
   symlist dw;
   static int mask[1+W_Nmax], *targets, *mighty;
   assert(W.N<=W_Nmax);
   if(targets==NULL) Init_Trans_Check(mask,&targets,&mighty);
   targets[0]=0;
   for (i=0;i<W.N;i++) for (j=mask[i];j<mask[i+1];j++){
     targets[j]=0;	/* union of targets with one variable less -> */
     for (k=0;k<=i;k++) targets[j]=targets[j]|targets[j&~mask[k]]; 
     for (k=0;(k<W.N)&&(mighty[targets[j]]<mighty[j]);k++)
       if (!(mask[k]&targets[j])) {	      /* targets=lower bounds only! */
         dw.d=W.d-W.w[k];				/* gradient degrees */
         for (l=0;l<W.M;l++) dw.m[l]=W.m[l]-W.z[l][k]; 	/* gradient phases */
         if (symcheck(dw,j,W,mask)) targets[j]=targets[j]|mask[k];}
     if (mighty[targets[j]]<mighty[j]) return 0;}
   return 1;
}
