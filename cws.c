#include "Global.h"
#include "LG.h"
#include "Rat.h"

FILE *inFILE, *outFILE;

#define	Only_IP_CWS 			1
#define TRANS_INFO_FOR_IP_WEIGHTS 	0

#define NFmax  10			/* maximal number of WS-files */

#define SIMPLEX_POINT_Nmax		50

typedef struct {  int  u[NFmax];  int nu; } 	CWS_type;

#define OSL (24)  /* opt_string's length */

void  PrintCWSUsage(char *c){
  int i;
  char *opt_string[OSL]={
"Options: -h        print this information",
"         -f        use as filter; otherwise parameters denote I/O files",
"         -w# [L H] make IP weight systems for #-dimensional polytopes.",
"                   For #>4 the lowest and highest degrees L<=H are required.",
"             -r/-t make reflexive/transversal weight systems (optional).",
"         -m# L H   make transverse weight systems for #-dimensional polytopes",
"                   (fast). The lowest and highest degrees L<=H are required.",
"         -c#       make combined weight systems for #-dimensional polytopes.",
"                   For #<=4 all relevant combinations are made by default,",
"                   otherwise the following option is required:",
"             -n[#] followed by the names wf_1 ... wf_# of weight files",
"                   currently #=2,3 are implemented.",
"              [-t] followed by # numbers n_i specifies the CWS-type, i.e.",
"                   the numbers n_i of weights to be selected from wf_i.",
"                   Currently all cases with n_i<=2 are implemented.",
"          -i       compute the polytope data M:p v [F:f] N:p [v] for all IP",
"                   CWS, where p and v denote the numbers of lattice points",
"                   and vertices of a dual pair of IP polytopes; an entry  ",
"                   F:f and no v for N indicates a non-reflexive `dual pair'.",
"          -d#      compute basic IP weight systems for #-dimensional",
"                   reflexive Gorenstein cones;",
"              -r#  specifies the index as #/2.",
"          -2       adjoin a weight of 1/2 to the input weight system.",
"          -N       make CWS for PPL in N lattice."
};
  printf("This is `%s': ", c);
  puts("create weight systems and combined weight systems.");
  printf("Usage:   %s -<options>; ", c);
  puts("the first option must be -w, -c, -i, -d or -h.");
  for (i=0;i<OSL;i++) puts(opt_string[i]);
  exit(0);
}

void Die(char *comment){ printf("\n%s\n",comment); exit(0);}

int Read_Weight(Weight *);

int READ_Weight(Weight *_W, FILE *INFILE){
  inFILE=INFILE;
  return Read_Weight(_W);
}

int READ_CWS_PP(CWS *_CW, PolyPointList *_P, FILE *INFILE){
  inFILE=INFILE;
  return Read_CWS_PP(_CW, _P);
}

void Print_CWS(CWS *_W);
void Conv(int narg, char* fn[]);
void SimplexPointCount(int narg, char* fn[]);
void Init_IP_Weights(int narg, char* fn[]);
void Init_moon_Weights(int narg, char* fn[]);
void Init_IP_CWS(int narg, char* fn[]);
void IP_Poly_Data(int narg, char* fn[]);
void Make_CWS_Points(CWS *, PolyPointList *);
void Npoly2cws(int narg, char* fn[]);
void RgcWeights(int narg, char* fn[]);
void AddHalf(void);
void PrintCWSextUsage(char *c)
{    printf( "This is `%s': -x gives undocumented extensions:\n",c);
       	puts("              -ip    printf PolyPointList");
       	puts("              -id    printf dual PolyPointList");
       	puts("              -p#    [infile1] [infile2] makes cartesian product"
     );	puts("                     of Vertices. # dimensions are identified.");
       	puts("              -S     count simplex points for weight system");
       	puts("              -L     count using LattE (-> count redcheck cdd)");
}

int  main (int narg, char* fn[])
{    inFILE=stdin; outFILE=stdout;
     if(narg==1) {printf("\nFor help type `%s -h'\n\n",fn[0]);exit(0);}
     if((fn[1][0]!='-')||(fn[1][1]=='h')) PrintCWSUsage(fn[0]);
     else if(fn[1][1]=='w') Init_IP_Weights(narg, fn);
     else if(fn[1][1]=='m') Init_moon_Weights(narg, fn);
     else if(fn[1][1]=='c') Init_IP_CWS(narg, fn);
     else if(fn[1][1]=='i') IP_Poly_Data(narg, fn);
     else if(fn[1][1]=='N') Npoly2cws(narg, fn);
     else if(fn[1][1]=='p') Conv(narg, fn);
     else if(fn[1][1]=='x') PrintCWSextUsage(fn[0]);
     else if(fn[1][1]=='S') SimplexPointCount(narg, fn);
     else if(fn[1][1]=='L') SimplexPointCount(narg, fn);
     else if(fn[1][1]=='d') RgcWeights(narg, fn);
     else if(fn[1][1]=='2') AddHalf();
     else printf("Unknown option '-%c'; use -h for help\n",fn[1][1]);
     return 0;
}


#define  WDIM      800000

typedef struct {
  int d, r2, allow11;
  Long wnum, winum, candnum;
  Long x[POLY_Dmax+1][POLY_Dmax];
  EqList q[POLY_Dmax];
  INCI qI[POLY_Dmax][EQUA_Nmax];
  int f0[POLY_Dmax];
  Equation wli[WDIM];
}                                                      RgcClassData;

int RgcWeicomp(Equation w1, Equation w2, int d) {
  /* w2-w1, i.e. pos for w1<w2,neg for w1>w2  */
  int i=d-1;
  if (w1.c-w2.c) return w1.c-w2.c;
  while((i)&&(w1.a[i]==w2.a[i])) i--;  
  return w2.a[i]-w1.a[i];
}

void RgcInsertat(Equation ww, int position, RgcClassData *X){
  int i;
  for(i=X->wnum-1; i>=position; i--) X->wli[i+1]=X->wli[i];
  X->wli[position]=ww;
  X->wnum++;
}

void RgcAddweight(Equation wn, RgcClassData *X){
  int i, j, p, n0, n1, k;
  if (X->wnum>=WDIM) {
    if(X->wnum>WDIM) return; 
    else {X->wnum++; printf("WDIM too small!\n");fflush(0);return;}}
  for(i=0;i<X->d;i++) if (2*wn.a[i]==-wn.c) return;
  for(i=0;i<X->d;i++) if (wn.a[i]==-wn.c) return;
  if (!X->allow11) for(i=0;i<X->d-1;i++) 
		     for(j=i+1;j<X->d;j++) if (wn.a[i]+wn.a[j]==-wn.c) return;
  X->candnum++;
  for(i=0;i<X->d-1;i++) for(p=i+1;p<X->d;p++) if (wn.a[i]>wn.a[p]) {
	k=wn.a[i]; wn.a[i]=wn.a[p]; wn.a[p]=k;}  /* make n0<=n1<=...<=n# */
  if (X->wnum) {
    i = RgcWeicomp(wn,X->wli[n0=0],X->d);
    if (!i) return;
    if (i>0) {RgcInsertat(wn,0,X); return;}
    i = RgcWeicomp(wn,X->wli[n1=X->wnum-1],X->d);
    if (!i) return; 
    if (i<0) {RgcInsertat(wn,X->wnum,X); return;}
    while(n1>n0+1) {
      p=(n0+n1)/2; 
      i=RgcWeicomp(wn,X->wli[p],X->d);
      if(!i) return;
      if(i>0) n1=p; 
      else n0=p;} 
    RgcInsertat(wn,n1,X);}
  else RgcInsertat(wn,0,X);
}

void PrintPoint(int n, RgcClassData *X){
  int i;
  assert(n<POLY_Dmax);
  for (i=0;i<n;i++) printf(" ");
  printf("X: ");
  for (i=0;i<X->d;i++) printf("%d ",(int) X->x[n][i]); 
  printf("\n");
}

void PrintQ(int n, RgcClassData *X){
  int i, j;
  assert(n<POLY_Dmax);
  for (i=0;i<n;i++) printf(" ");
  printf("q: ne=%d\n",X->q[n].ne);
  for (j=0;j<X->q[n].ne;j++){
    printf("%d  ", (int) X->q[n].e[j].c);
    for (i=0; i<X->d; i++) printf("%d ", (int) X->q[n].e[j].a[i]);
    printf("\n");}
}

void PrintEquation(Equation *q, int d/*, char *c, int j*/){
  int i;
  printf("%d  ", (int) -q->c);
  for (i=0;i<d;i++) printf("%d ", (int) q->a[i]);
  /*printf("  %s  np=%d\n", c, j);*/
}

int LastPointForbidden(int n, RgcClassData *X){
  int l;
  Long *y = X->x[n];
  Long ysum=0, ymax=0;
  assert(n<X->d);
  for (l=0;l<X->d;l++){
    ysum += y[l];
    if (y[l]>ymax) ymax=y[l];}
  if (ysum < 2) return 1;
  if (ysum == 2) if ((!X->allow11) ||(ymax == 2))	return 1;
  if (X->r2 == 2) if (ymax < 2) return 1;
  if (X->r2 == 1) if (ymax < 3) return 1;
  for (l=X->f0[n-1];l<X->d-1;l++) if (y[l]<y[l+1]) return 1;
  return 0;
}

void ComputeQ0(RgcClassData *X){
  int i, j;
  X->q[0].ne = X->d;
  for (i=0; i<X->q[0].ne; i++){
    for (j=0; j<X->d; j++) X->q[0].e[i].a[j] = 0;
    if (X->r2%2) {X->q[0].e[i].a[i] = X->r2; X->q[0].e[i].c = -2;}
    else {X->q[0].e[i].a[i] = X->r2/2; X->q[0].e[i].c = -1;}}
  X->f0[0] = 0;
  X->qI[0][0] = INCI_1();
  for (i=1; i<X->q[0].ne; i++) X->qI[0][i] = INCI_PN(X->qI[0][i-1],1); 
}

int IsRedundant(INCI newINCI, INCI *qINew, int ne){
  int i;
  for (i=0;i<ne;i++) if (INCI_LE(qINew[i], newINCI)) return 1;
  return 0;
}

void ComputeQ(int n, RgcClassData *X){
  /* q[n] from q[n-1], x[n] */
  int i, j, k;
  Long *y = X->x[n];
  Long yqOld[EQUA_Nmax];
  EqList *qOld = &X->q[n-1], *qNew = &X->q[n];
  INCI *qIOld = X->qI[n-1], *qINew = X->qI[n];
  INCI newINCI;
  assert(n<X->d);
  for (i=X->d-1;(i>=X->f0[n-1])&&(y[i]==0);i--) ;
  X->f0[n] = ++i;
  qNew->ne = 0;
  for (i=0;i<qOld->ne;i++) 
    if (!(yqOld[i] = Eval_Eq_on_V(&(qOld->e[i]),y,X->d))){
      qNew->e[qNew->ne]=qOld->e[i];
      qINew[qNew->ne]=qIOld[i];
      (qNew->ne)++;}
  for (i=0;i<qOld->ne-1;i++) 
    for (j=i+1;j<qOld->ne;j++) 
      if (yqOld[i]*yqOld[j]<0)
	if (INCI_abs(newINCI=INCI_OR(qIOld[i],qIOld[j]))<=n+1)
	  if (!IsRedundant(newINCI, qINew, qNew->ne)){
	    for (k=qNew->ne-1; k>=0; k--) if (INCI_LE(newINCI,qINew[k])){
		qINew[k] = qINew[qNew->ne-1];
		qNew->e[k] = qNew->e[qNew->ne-1];
		qNew->ne--;}
	    assert(qNew->ne < EQUA_Nmax-1);
	    qINew[qNew->ne] = newINCI;
	    qNew->e[qNew->ne] = EEV_To_Equation(&qOld->e[i],&qOld->e[j],y,X->d);
	    if (qNew->e[qNew->ne].c > 0){
	      for (k=0;k<X->d;k++) qNew->e[qNew->ne].a[k] *= -1;
	      qNew->e[qNew->ne].c *= -1;}
	    qNew->ne++;}
}

Long Flcm(Long a, Long b){
  return (a*b)/Fgcd(a,b);
}

void Cancel(Equation *q, int d){
  Long gcd = -q->c;
  int j;
  for (j=0;j<d;j++) gcd = Fgcd(gcd,q->a[j]);
  if (gcd > 1){
    for (j=0;j<d;j++) q->a[j] /= gcd;
    q->c /= gcd;}
}

int ComputeAndAddAverageWeight(Equation *q, int n, RgcClassData *X){
  int i, j;
  if (X->q[n].ne < X->d - n) return 0;
  q->c = -1;
  for (i=0;i<X->q[n].ne;i++) {
    if (X->q[n].e[i].c >= 0) {PrintQ(n,X); exit(0);}
    q->c = -Flcm(-q->c,-X->q[n].e[i].c);}
  for (j=0;j<X->d;j++){
    q->a[j] = 0;
    for (i=0;i<X->q[n].ne;i++) 
      q->a[j]+=X->q[n].e[i].a[j]*(q->c/X->q[n].e[i].c);
    if (q->a[j]<=0){ assert(q->a[j]==0); return 0;}    }
  q->c *= X->q[n].ne;
  Cancel(q, X->d);
  RgcAddweight(*q, X);
  return 1;
}

void ComputeAndAddLastQ(RgcClassData *X){
  /* q[d-1] from q[d-2], x[d-1] */
  int i, d = X->d;
  Equation q;
  Equation *q0 = &X->q[d-2].e[0], *q1 = &X->q[d-2].e[1];
  Long *y = X->x[d-1];
  Long yq0 = Eval_Eq_on_V(q0,y,d), yq1;
  if (LastPointForbidden(d-1, X)) return;
  if (!yq0) return;
  yq1 = Eval_Eq_on_V(q1,y,d);
  if (yq0 < 0){
    if (yq1 <= 0) return; 
    yq0 *= -1;}
  else {
    if (yq1 >= 0) return; 
    yq1 *= -1;}
  q.c = yq0 * q1->c + yq1 * q0->c;
  for (i=0;i<d;i++) q.a[i] = yq0 * q1->a[i] + yq1 * q0->a[i];
  Cancel(&q, d);
  RgcAddweight(q, X);
}

void RecConstructRgcWeights(int n, RgcClassData *X){
  /* we have q[n-1], x[n] */
  int k, l;
  Equation q;
  Long yq[POLY_Dmax];
  Long *y = X->x[n+1];
  if (n == 0) ComputeQ0(X);
  else if (LastPointForbidden(n, X)) return;
  else ComputeQ(n, X);
  if (!ComputeAndAddAverageWeight(&q, n, X)) return;
  if (n >= X->d-1) return;
  /* Examine all integer points of simplex:                               */
  for (k=0;k<X->d-1;k++) {y[k]=0; yq[k]=0;}  /* sets k=d-1; important!    */
  y[k]=-1;           /* starting point just outside                       */
  yq[k]=-q.a[k];
  while(k>=0){
    y[k]++; 
    yq[k]+=q.a[k];
    for (l=k+1;l<X->d;l++) yq[l]=yq[k];
    if (n==X->d-2){assert(X->q[X->d-2].ne==2); ComputeAndAddLastQ(X);}
    else RecConstructRgcWeights(n+1, X);
    for(k=X->d-1;(k>=0 ? (yq[k]+q.a[k]>=-q.c) : 0);k--) y[k]=0;}
         /* sets k to the highest value where y[k] didn't exceed max value; 
            resets the following max values to min values                 */
}

void AddPointToPoly(Long *y, PolyPointList *P){
  int i;
  for (i=0;i<P->n;i++) P->x[P->np][i] = y[i];
  P->np++;
}

int WsIpCheck(Equation *q, int d){
  int k,l;
  PolyPointList *P = (PolyPointList *) malloc(sizeof(PolyPointList));
  assert (P != NULL);
  VertexNumList V;
  EqList E;
  Long y[POLY_Dmax];
  Long yq[POLY_Dmax];
  P->n=d;
  P->np=0;
  for (k=0;k<d;k++) {y[k]=0; yq[k]=0;}
  k = d-1;
  AddPointToPoly(y, P);
  y[k]=-1;           /* starting point just outside                       */
  yq[k]=-q->a[k];
  while(k>=0){
    y[k]++; 
    yq[k]+=q->a[k];
    for (l=k+1;l<d;l++) yq[l]=yq[k];
    if (yq[k]==-q->c) AddPointToPoly(y, P);
    for(k=d-1;(k>=0 ? (yq[k]+q->a[k]>-q->c) : 0);k--) y[k]=0;}
         /* sets k to the highest value where y[k] didn't exceed max value; 
            resets the following max values to min values                 */
  if (P->np <= d) {free(P); return 0;}
  Find_Equations(P, &V, &E);
  if (E.ne < d) {free(P); return 0;}
  for (k=0;k<d;k++) y[k]=1;
  for (k=0;k<E.ne;k++) 
    if (Eval_Eq_on_V(&(E.e[k]), y, d) <= 0) if (!E.e[k].c) {free(P); return 0;}
  k=P->np-1;
  free(P);   
  return k;
}

void RgcWeights(int narg, char* fn[])
{    int i, j, d, n=1, r2=2;
     char *c=&fn[1][2];
     RgcClassData *X = (RgcClassData *) malloc(sizeof(RgcClassData));
     if(narg>2) if(c[0]==0) c=fn[++n];
     if(!IsDigit(c[0])){ puts("-d must be followed by a number");exit(0);}
     if(POLY_Dmax<(d=atoi(c))){printf("Increase POLY_Dmax to %d\n",d);exit(0);}
     if(narg>++n){
       if((fn[n][0]!='-')||(fn[n][1]!='r')){
	 printf("the second option has to be of the type -r\n");exit(0);}
       c=&fn[n][2]; r2=atoi(c);}
     X->d=d; X->r2=r2; X->wnum=0; X->winum=0; X->candnum=0; X->allow11=0;
     RecConstructRgcWeights(0,X);
     if (X->wnum <= WDIM){
       for (i=0;i<X->wnum;i++) {
	 j=WsIpCheck(&X->wli[i], d); 
	 if (j) {
	   PrintEquation(&X->wli[i], X->d);
	   printf("  np=%d\n",j);
	   X->winum++;
	   /*else PrintEquation(&X->wli[i], X->d, "n");*/}}}
     printf("#ip=%ld, #cand=%ld(%ld)\n", 
	    X->winum,  X->wnum,  X->candnum);
}

int IsNextDigit(void);

void AddHalf(void)
{
  int IN[AMBI_Dmax*(AMBI_Dmax+1)]; 
  int i, j, n = 0;

  inFILE=stdin; 
  n=0;
  while (1){
    for(i=0;i<AMBI_Dmax*(AMBI_Dmax+1);i++)    { 
      char c;
      while(' '==(c=fgetc(inFILE )));
      ungetc(c,inFILE);    /* read blanks */
      if(IsNextDigit()) fscanf(inFILE,"%d",&IN[i]); else break;  }
    if (i==0) break;
    n++;
    while('\n'!=fgetc(inFILE));      /* read to end of line */
    for (j=0; j<i; j++) printf("%d ", IN[j]);
    printf("%d\n", IN[0]/2);}
  printf("#ws=%d\n",n);
}







void Make_IP_Weights(int d, int Dmin, int Dmax,int rFlag,int tFlag);
void MakeMoonWeights(int d, int Dmin, int Dmax);
void Make_34_Weights(int d, int tFlag);
void Init_IP_Weights(int narg, char* fn[])
{    int n=1,d,L=0,H=0,rf=0,tf=0; char *c=&fn[1][2];
     if(narg>2) if(c[0]==0) c=fn[++n];
     if(!IsDigit(c[0])){ puts("-w must be followed by a number");exit(0);}
     if(POLY_Dmax<(d=atoi(c))){printf("Increase POLY_Dmax to %d\n",d);exit(0);}
     if(++n<narg) if((fn[n][0]!='-')&&(IsDigit(fn[n][0])))
     {	L=atoi(fn[n]); assert(++n<narg); assert(IsDigit(fn[n][0]));
	H=atoi(fn[n]); assert(L<=H); n++;
     }	n--;
     while(++n<narg) if(fn[n][0]=='-')
     {  if(fn[n][1]=='r')rf=1; else if(fn[n][1]=='t')tf=1; 
	else {printf("Illegal option %s\n",fn[n]);exit(0);}
     }
     if(++n<narg) {printf("Want %s as output file?\n",fn[n]);exit(0);}
     if(n<narg) {puts("Too many arguments");exit(0);}
     if(H) Make_IP_Weights(d,L,H,rf,tf); else Make_34_Weights(d,tf);
     {  ;
     }
}

void Init_moon_Weights(int narg, char* fn[]){
  int n=1,d,L=0,H=0; char *c=&fn[1][2];
  if(narg>2) if(c[0]==0) c=fn[++n];
  if(!IsDigit(c[0])){ puts("-m must be followed by a number");exit(0);}
  if(POLY_Dmax<(d=atoi(c))){printf("Increase POLY_Dmax to %d\n",d);exit(0);}
  if(++n<narg) if((fn[n][0]!='-')&&(IsDigit(fn[n][0]))) {
      L=atoi(fn[n]); assert(++n<narg); assert(IsDigit(fn[n][0]));
      H=atoi(fn[n]); assert(L<=H); n++;     }
  n--;
  while(++n<narg) if(fn[n][0]=='-')     {
      printf("Illegal option %s\n",fn[n]); exit(0);}     
  if(++n<narg) {printf("Want %s as output file?\n",fn[n]); exit(0);}
  if(n<narg) {puts("Too many arguments");exit(0);}
  if(H) MakeMoonWeights(d+1,L,H);
  else  {puts("Please give lowest and highest d");exit(0);}
}

int  VP_2_CWS(Long *V[], int d, int v, CWS *W);
void Npoly2cws(int narg, char* fn[])
{    int n=2; CWS W; EqList E; VertexNumList V; Long *X[VERT_Nmax]; FILE *OF;
     PolyPointList *P=(PolyPointList *)malloc(sizeof(PolyPointList)); 
     assert(P!=NULL); assert(!strcmp(fn[1],"-N")); inFILE=stdin;outFILE=stdout;
     if(narg>2) {if(fn[2][0]=='-'){assert(fn[2][1]=='f');inFILE=NULL;} else
     {	inFILE=fopen(fn[2],"r"); assert(NULL!=inFILE);
        if(narg>3) {outFILE=fopen(fn[3],"w"); assert(NULL!=outFILE);}
     }}	OF=outFILE; while(Read_CWS_PP(&W,P))
     {  if(W.N) Die("Only PPL-input in Npoly2cws!");
	if(!IP_Check(P,&V,&E)) Die("Not IP!");
	Sort_VL(&V);
	for(n=0;n<V.nv;n++) X[n]=P->x[V.v[n]];
	if(VP_2_CWS(X,P->n,V.nv,&W)) {Print_CWS(&W);fprintf(outFILE,"\n");}
	else {outFILE=stderr; Print_PPL(P,"CWS not found"); outFILE=OF;}
     }
}

void Make_IP_CWS(int narg, char* fn[]);
void Make_34_CWS(int d);
void Init_IP_CWS(int narg, char* fn[])
{    int d,n=1,nop=0; char *c=&fn[1][2]; if(narg>2) if(c[0]==0) c=fn[++n];
     if(!IsDigit(c[0])){ puts("-c must be followed by a number");exit(0);}
     if(POLY_Dmax<(d=atoi(c))){printf("Increase POLY_Dmax to %d\n",d);exit(0);}
     if(++n<narg) if((fn[n][0]=='-')&&(fn[n][1]=='n')) nop=1;
     if(nop) Make_IP_CWS(narg,fn); else if(d<=4) Make_34_CWS(d); else Die(
       "`-c#' has to be followed by `-n' and the weight file names for dim>4");
}

/*  ==========             ALL  IP  WEIGHTS  in  d <= 4     	==========  */

#define  INFO      0

#define  lcm(a,b)  ((a)*(b)/NNgcd((a),(b)))
typedef struct {int n[W_Nmax+1];} weights;
typedef Rat ratmat[W_Nmax][W_Nmax];
typedef Rat ratvec[W_Nmax];
typedef struct { int wnum, N, points[W_Nmax][W_Nmax], nsubsets[W_Nmax-1],
	subsets[W_Nmax-1][10][W_Nmax];  weights wli[WDIM]; }	WSaux;

int  weicomp(weights w1,weights w2,int *_N) 
                                /* w2-w1, i.e. pos for w1<w2,neg for w1>w2  */
{    int i=*_N;
     while((i)&&(w1.n[i]==w2.n[i])) i--;
     return w2.n[i]-w1.n[i];
}
void insertat(WSaux *X,weights ww, int position)
{    int i, j;
     for(i=X->wnum-1;i>=position;i--)
       for(j=0;j<=X->N;j++) X->wli[i+1].n[j]=X->wli[i].n[j];
     for(j=0;j<=X->N;j++) X->wli[position].n[j]=ww.n[j];
     X->wnum++;
}

void addweight(WSaux *X,weights wn)
{    int i, p, n0, n1, k;
     if (X->wnum>=WDIM) {
       if(X->wnum>WDIM) return; 
       else {X->wnum++; printf("WDIM too small!\n");fflush(0);return;}}
     for(i=0;i<X->N-1;i++) for(p=i+1;p<X->N;p++) if (wn.n[i]>wn.n[p]) 
       {k=wn.n[i]; wn.n[i]=wn.n[p]; wn.n[p]=k;}  /* make n0<=n1<=...<=n# */
     if (X->wnum) {if ((i = weicomp(wn,X->wli[n0=0],&X->N)))
               {if (i>0) {insertat(X,wn,0); return;}} else return;
            if ((i = weicomp(wn,X->wli[n1=X->wnum-1],&X->N)))
               {if (i<0) {insertat(X,wn,X->wnum); return;}} else return;
      while(n1>n0+1) {p=(n0+n1)/2; i=weicomp(wn,X->wli[p],&X->N);
                  if(i) {if(i>0) n1=p; else n0=p;} else return;}
      insertat(X,wn,n1);}
     else insertat(X,wn,0);
}
int checkwrite(WSaux *X,weights ws){
  int i;
  for (i=0;i<X->N;i++) if (!ws.n[i]) return 0;
  addweight(X,ws);
  return 1;
}
weights testweisys(WSaux *X,int npoints){
  weights tws;
  ratmat rm;
  ratvec newboundwei, boundwei[10], rattws, rs;
  int i, ii, j, k, nboundwei=0, New, rankrm, one[W_Nmax];
  Long minnbw, maxnbw;
  for (k=0;k<X->nsubsets[npoints-2];k++   /* alle 0-systeme */ ){
    rankrm=0;
    for (i=0;i<npoints;i++){
      for (j=0;j<npoints;j++){
        rm[i][j]=rI(X->points[i][X->subsets[npoints-2][k][j]]);}
      rs[i]=rI(1);   }
    for (i=0;i<npoints;i++){
      one[i]=-1;
      for (j=0;(j<npoints)&&(one[i]<0);j++) if (rm[i][j].N) {
        one[i]=j;
        rankrm++;   }
      if (one[i]>=0){
        for (j=one[i]+1;j<npoints;j++) /* normalize i'th line */
          rm[i][j]=rQ(rm[i][j],rm[i][one[i]]);
        rs[i]=rQ(rs[i],rm[i][one[i]]);
        rm[i][one[i]]=rI(1);
        for (ii=0;ii<npoints;ii++) if (ii-i){ /* subtract multiple of i'th
                    line from ii'th line */
          for (j=0;j<npoints;j++) 
            if (j!=one[i]) rm[ii][j]=rD(rm[ii][j],rP(rm[ii][one[i]],rm[i][j]));
          rs[ii]=rD(rs[ii],rP(rm[ii][one[i]],rs[i]));
          rm[ii][one[i]]=rI(0);  }   }    }
    for (j=0;j<X->N;j++) newboundwei[j]=rI(0);
    if (rankrm==npoints) for (i=0;i<npoints;i++) 
      newboundwei[X->subsets[npoints-2][k][one[i]]]=rs[i];
    minnbw=0; maxnbw=0; 
    for (j=0;j<X->N;j++) {
      minnbw=min(newboundwei[j].N,minnbw);
      maxnbw=max(newboundwei[j].N,maxnbw); }
    New=((minnbw>=0)&&(maxnbw>0));
    for (i=0;New&&(i<nboundwei);i++) {
      New=0; 
      for (j=0;j<X->N;j++) if (rD(newboundwei[j],boundwei[i][j]).N) New=1; }
    if (New) {
      for (j=0;j<X->N;j++) boundwei[nboundwei][j]=newboundwei[j]; 
      nboundwei++;    }   }
  for (j=0;j<X->N;j++) {
    rattws[j]=rI(0);
    for (i=0;i<nboundwei;i++) rattws[j]=rS(rattws[j],boundwei[i][j]);
    rattws[j]=rQ(rattws[j],rI(max(nboundwei,1)));}
  tws.n[X->N]=1;
  for (j=0;j<X->N;j++) tws.n[X->N]=lcm(tws.n[X->N],rattws[j].D);
  for (j=0;j<X->N;j++) tws.n[j]=rP(rI(tws.n[X->N]),rattws[j]).N;
  return tws;
}
void createweights(WSaux *X,int npoints){
  int x0, x1, x2, x3, x4, sum, maxx;  weights tws; tws=testweisys(X,npoints);
  if(checkwrite(X,tws)) if (npoints<X->N){
    for (x0=0;x0*tws.n[0]<tws.n[X->N];x0++)
    for (x1=0;x0*tws.n[0]+x1*tws.n[1]<tws.n[X->N];x1++)
    for (x2=0;x0*tws.n[0]+x1*tws.n[1]+x2*tws.n[2]<tws.n[X->N];x2++)
/* #if (N>3) */    for (x3=0; (x3==0) || ((X->N>3) &&
	(x0*tws.n[0]+x1*tws.n[1]+x2*tws.n[2]+x3*tws.n[3]<tws.n[X->N])) ;x3++)
/* #if (N>4) */    for (x4=0; (x4==0) || ((X->N>4) && (x0*tws.n[0]+
	x1*tws.n[1]+x2*tws.n[2]+x3*tws.n[3]+x4*tws.n[4]<tws.n[X->N]));x4++)
    { sum=0; maxx=0;
      X->points[npoints][0]=x0; sum+=x0; maxx=max(maxx,x0);
      X->points[npoints][1]=x1; sum+=x1; maxx=max(maxx,x1);
      X->points[npoints][2]=x2; sum+=x2; maxx=max(maxx,x2);
  if (X->N>3)    
      {X->points[npoints][3]=x3; sum+=x3; maxx=max(maxx,x3);}
  if (X->N>4)    
      {X->points[npoints][4]=x4; sum+=x4; maxx=max(maxx,x4);}
      if ((sum>2)&&(maxx>1)) createweights(X,npoints+1);
      /* if (npoints<3) {printf("%d",npoints); fflush(0);}*/}  }
}
void makesubsets(WSaux *X){
  int i, p0,p1,p2,p3,p4;
  for (i=0;i<X->N-1;i++) X->nsubsets[i]=0;
  for (p0=0;p0<X->N-1;p0++) for (p1=p0+1;p1<X->N;p1++) {
    X->subsets[0][X->nsubsets[0]][0]=p0; 
    X->subsets[0][X->nsubsets[0]][1]=p1;
    X->nsubsets[0]++;
    for (p2=p1+1;p2<X->N;p2++) {
      X->subsets[1][X->nsubsets[1]][0]=p0; 
      X->subsets[1][X->nsubsets[1]][1]=p1;
      X->subsets[1][X->nsubsets[1]][2]=p2;
      X->nsubsets[1]++;
      for (p3=p2+1;p3<X->N;p3++) {
        X->subsets[2][X->nsubsets[2]][0]=p0; 
        X->subsets[2][X->nsubsets[2]][1]=p1;
        X->subsets[2][X->nsubsets[2]][2]=p2;
        X->subsets[2][X->nsubsets[2]][3]=p3;
        X->nsubsets[2]++;
        for (p4=p3+1;p4<X->N;p4++) {
          X->subsets[3][X->nsubsets[3]][0]=p0; 
          X->subsets[3][X->nsubsets[3]][1]=p1;
          X->subsets[3][X->nsubsets[3]][2]=p2;
          X->subsets[3][X->nsubsets[3]][3]=p3;
          X->subsets[3][X->nsubsets[3]][4]=p4;
          X->nsubsets[3]++;  }  }  }  } 
}
void WRITE_Weight(Weight *_W);
void Make_34_Weights(int d, int tFlag)
{    int i, Info=0; WSaux *X = (WSaux *) malloc(sizeof(WSaux)); 
     PolyPointList *P = (PolyPointList *) malloc(sizeof (PolyPointList));
     assert(P!=NULL); assert(X!=NULL); X->wnum=1; X->N=d+1; assert(d<=4); 
     makesubsets(X);	for (i=0;i<X->N;i++) X->points[0][i]=1;
     X->wli[0].n[X->N]=X->N;	for (i=0;i<X->N;i++) X->wli[0].n[i]=1;
  if (X->N>4){ X->points[1][0]=4; for (i=1;i<X->N;i++) X->points[1][i]=0;
    createweights(X,2); if(Info){printf("Did (4,0,0,0,0)\n");fflush(0);} 
    X->points[1][0]=3; X->points[1][1]=1; for(i=2;i<X->N;i++)X->points[1][i]=0;
    createweights(X,2);  if(Info){printf("Did (3,1,0,0,0)\n"); fflush(0);}
    X->points[1][0]=2; X->points[1][1]=2; for(i=2;i<X->N;i++)X->points[1][i]=0;
    createweights(X,2);  if(Info){printf("Did (2,2,0,0,0)\n"); fflush(0);}
    X->points[1][0]=2; X->points[1][1]=1; X->points[1][2]=1;
      for (i=3;i<X->N;i++) X->points[1][i]=0;
    createweights(X,2);  if(Info){printf("Did (2,1,1,0,0)\n"); fflush(0);}}
  if (X->N>3){ X->points[1][0]=3; for (i=1;i<X->N;i++) X->points[1][i]=0;
    createweights(X,2);  if(Info){printf("Did (3,0,0,0,0)\n"); fflush(0);}
    X->points[1][0]=2; X->points[1][1]=1; for(i=2;i<X->N;i++)X->points[1][i]=0;
    createweights(X,2);  if(Info){printf("Did (2,1,0,0,0)\n"); fflush(0);}}         X->points[1][0]=2; for (i=1;i<X->N;i++) X->points[1][i]=0;
  createweights(X,2);  if(Info){printf("Did (2,0,0,0,0)\n"); fflush(0);}
  Info=0; for (i=0;i<X->wnum;i++){ int j; Weight W; VertexNumList V; EqList E;
    W.d=X->wli[i].n[W.N=X->N]; for(j=0;j<X->N;j++)W.w[j]=X->wli[i].n[j]; W.M=0;
    Make_Poly_Points(&W, P); if(Ref_Check(P,&V,&E))
    { int t=Trans_Check(W); char c[5]="  rt"; if(t||!tFlag) {c[3]=(t)?'t':0; 
	if(Info++) puts("");
	WRITE_Weight(&W); printf("%s",c);}}
  } fprintf(outFILE,"  #=%d  #cand=%d\n",Info,X->wnum); free(X);
}
/*  ==========       End of  ALL  IP  WEIGHTS  in  d <= 4    	==========  */
/*  ==========  	      MAKE WEIGHTS d>4:                	==========  */

void WRITE_Weight(Weight *_W)
{    int n; fprintf(outFILE,"%d ",(int) _W->d);
     for(n=0; n < _W->N; n++) fprintf(outFILE," %d",(int) _W->w[n]); 
}

int IfIpWWrite(Weight *W, PolyPointList *P, int *rFlag, int *tFlag)
{    VertexNumList V; EqList E; int r=1,i=-1; Make_Poly_Points(W, P); 
     if(IP_Check(P,&V,&E)){
       while(r && (++i < E.ne)) if(E.e[i].c != 1) r = 0;
       if(*tFlag && Trans_Check(*W)){
	 WRITE_Weight(W); if(r) fprintf(outFILE," r"); 
	 fprintf(outFILE,"\n");fflush(stdout); return 1;}
       if(*rFlag && r){Write_Weight(W); fflush(stdout); return 1;} 
       if(!*tFlag && !*rFlag)
	 {WRITE_Weight(W); if(r) fprintf(outFILE," r"); 
#if(TRANS_INFO_FOR_IP_WEIGHTS)
	   if(Trans_Check(*W)) fprintf(outFILE,"%st", r ? "" : " "); 
#endif
	   fprintf(outFILE,"\n"); fflush(stdout); return 1;}
       return 0;
     } 
     else return 0;
}
void Rec_IpWeights(Weight *W, PolyPointList *P, int g, int sum, int *npp, 
	int *nrp, int n, int *rFlag, int *tFlag)
{    int wmax=W->d/(W->N-n+1); wmax=min(wmax,W->w[n+1]); 
     wmax=min(wmax,sum-n);
     if(n) for(W->w[n]=wmax;(n+1)*W->w[n]>=sum;W->w[n]--)
       Rec_IpWeights(W,P,Fgcd(g,W->w[n]),sum-W->w[n],npp,nrp,n-1,rFlag,tFlag);
     else if(1==Fgcd(g,W->w[0]=sum)) {
       (*npp)++;if(IfIpWWrite(W,P,rFlag,tFlag))(*nrp)++;};
}
void MakeIpWeights(int N, int from_d, int to_d, int *rFlag, int *tFlag)
{    int npp=0, nrp=0; Weight W;  
     PolyPointList *P = (PolyPointList *) malloc (sizeof(PolyPointList)); 
     assert((N<=W_Nmax)&&(N<POLY_Dmax+2)); assert(P!=NULL);W.N=N; W.M=0;
     for(W.d=from_d;W.d<=to_d;W.d++)
        for(W.w[N-1]=W.d/2; W.d <= N*W.w[N-1]; W.w[N-1]--)
	Rec_IpWeights(&W,P,Fgcd(W.d,W.w[W.N-1]), W.d-W.w[W.N-1],&npp,&nrp,
		      N-2,rFlag,tFlag);
     if(*rFlag) fprintf(outFILE,"#primepartitions=%d #refpolys=%d\n",npp,nrp);
    if(*tFlag)fprintf(outFILE,"#primepartitions=%d #transpolys=%d\n",npp,nrp);
     if(!*rFlag && !*tFlag)
       fprintf(outFILE,"#primepartitions=%d #IPpolys=%d\n",npp,nrp);
     exit(0);
}
void Make_IP_Weights(int d, int Dmin, int Dmax, int rFlag, int tFlag)
{
  MakeIpWeights(d+1, Dmin, Dmax, &rFlag, &tFlag);
}


#define MOONSHINE_CRITERIA (0)
/* if set to (1), weights are written according to the criteria
   formulated in IfMoonWWrite, otherwise if they satisfy Trans_Check */

#if(MOONSHINE_CRITERIA)

int IfMoonWWrite(Weight *W, Rat *dmwow, long *nsmallchi, int *stf){
  /* Write the weights determined by the algorithm;
     if desired, put in extra conditions like the ones below 
     (here: chi should be a small multiple of 24)                   */
  int i, D;
  Rat chi = rI(0);
  Rat auxrat;
  for (D=1; D<=W->d; D++)
    if (!(W->d % D)){
      /* sf = Strangefun(W->d / D); */
      auxrat = rR(stf[W->d / D], W->d);
      for(i=0; i<W->N; i++)
	if (!(D * W->w[i] % W->d)) {
	  auxrat = rP(auxrat, dmwow[i]);
	  auxrat.N *= -1;}
      chi = rS(auxrat, chi);}
  if (chi.D != 1) return 0;
  /* exit if the weight doesn't satisfy certain desired criteria */
  if (chi.N % 24) return 0;
  if (abs(chi.N) > 96) return 1;  
  (*nsmallchi)++; 
  WRITE_Weight(W);
  printf("  chi=%ld\n", chi.N);
  fflush(0);
  return 1;
}

int BlaRout(int quot){/* compute |{(i,j): gcd(i,j,quot)=1}| */
  int sf = 0, bla, qob;
  int i,j;
  for (i=1;i<=quot;i++){
    bla = Fgcd(i,quot);
    qob = quot/bla;
    for (j=1;j<=bla;j++)
      if (Fgcd(j,bla) == 1) sf += qob;}
  return sf;
}

#endif

void RecMoonWeights(Weight *W, int g, int sum, long *npp, long *nintPP1,
		    long *nintchi, int n, long long *PP1N, long long *PP1D
#if(MOONSHINE_CRITERIA)
		    , long *nsmallchi, Rat*dmwow, int *stf
#endif
		    ){
  int wmax=W->d/(W->N-n+1);
  wmax=min(wmax,W->w[n+1]); 
  wmax=min(wmax,sum-n);
  if(n)
    for(W->w[n]=wmax; (n+1)*W->w[n]>=sum; W->w[n]--){
#if(MOONSHINE_CRITERIA)
      dmwow[n] = rR(W->d - W->w[n], W->w[n]);
#endif
      /*PP1[n] = rP(dmwow[n], PP1[n+1]);*/
      PP1N[n] = PP1N[n+1] * (long long) (W->d - W->w[n]); assert(PP1N[n]>0);
      PP1D[n] = PP1D[n+1] * (long long) W->w[n];
      RecMoonWeights(W, Fgcd(g,W->w[n]), sum-W->w[n], npp, nintPP1, nintchi,
		     n-1, PP1N, PP1D
#if(MOONSHINE_CRITERIA)
		     , nsmallchi, dmwow, stf
#endif
		     );}
  else if (1 == Fgcd(g,sum)){
    W->w[0]=sum;
    /*PP1[n] = rP(dmwow[n], PP1[n+1]);*/
    PP1N[n] = PP1N[n+1] * (long long) (W->d - W->w[n]);
    PP1D[n] = PP1D[n+1] * (long long) W->w[n];
    (*npp)++;
    /*if (PP1->D != 1) return;*/
    if (PP1N[0] % PP1D[0]) return;
    if (PP1N[n]/PP1N[n+1] != W->d - W->w[n]) return;
    (*nintPP1)++;
#if(MOONSHINE_CRITERIA)
    dmwow[n] = rR(W->d - W->w[n], W->w[n]);
    if(IfMoonWWrite(W, dmwow, nsmallchi, stf)) (*nintchi)++;
#else
    if(Trans_Check(*W)){
      WRITE_Weight(W); fprintf(outFILE,"\n");fflush(stdout); (*nintchi)++;}
#endif
  }
}

void MakeMoonWeights(int N, int from_d, int to_d){
  /* Fast routine for creating candidates for transverse weight systems;
     conditions: Poincare-polynomial at t=1 integer, chi integer */
  long int npp=0, nintPP1=0, nintchi=0;
  Weight W;
  /* Rat PP1[W_Nmax]; */   /* Poincare polynomial evaluated at t=1 */
  long long PP1N[W_Nmax], PP1D[W_Nmax];
#if(MOONSHINE_CRITERIA)
  long int nsmallchi=0, i;
  Rat dmwow[W_Nmax]; /* dmwow[i] = rR(W.d-W.w[i], W.w[i]); "d minus w over w" */
  int *stf = (int *) malloc((to_d +1) * sizeof(int));
  for (i=1; i<=to_d; i++) stf[i] = BlaRout(i);
#endif
  assert((N<=W_Nmax)&&(N<POLY_Dmax+2));
  W.N=N;
  W.M=0;
  for(W.d=from_d; W.d<=to_d; W.d++){
#if(MOONSHINE_CRITERIA)
    /* if results are rare, create output to show that the program is still
       doing something  */
    if (W.d>607) printf("d=%d\n",W.d);
    fflush(0);
#endif
    for(W.w[N-1]=W.d/2; W.d <= N*W.w[N-1]; W.w[N-1]--){
      /*PP1[N-1] = rR(W.d-W.w[W.N-1], W.w[W.N-1]);*/
      PP1N[N-1] = W.d-W.w[W.N-1];
      PP1D[N-1] = W.w[W.N-1];
#if(MOONSHINE_CRITERIA)
     dmwow[N-1] = rR(W.d-W.w[W.N-1], W.w[W.N-1]);
#endif
     RecMoonWeights(&W, Fgcd(W.d, W.w[W.N-1]), W.d-W.w[W.N-1], &npp, &nintPP1,
		     &nintchi, N-2, PP1N, PP1D
#if(MOONSHINE_CRITERIA)
		     , &nsmallchi, dmwow, stf
#endif
		     );}}
#if(MOONSHINE_CRITERIA)
  fprintf(outFILE,"#partitions=%ld #intPP1=%ld #intchi=%ld #smallchi=%ld\n",
	  npp, nintPP1, nintchi, nsmallchi);
#else
  fprintf(outFILE,"#partitions=%ld #intPP1=%ld #trans=%ld\n",
	  npp, nintPP1, nintchi);
#endif
  exit(0);
}


	/* ----------   LG/transversal stuff  ------------ */
#define  ALLOWHALF	(1)		      /* i.e. trivial LG potentials */
#define  CHAT		(0)			/* 3 ... for positive c_1 */
#define  TWDIM  16384      /* 16384  8192  4096  dimension of weight-buffer */
#define  mod(a,b)  ((a)%(b))

typedef int T_weight[AMBI_Dmax+2];			/* NM::AMBI_Dmax */
typedef struct {int n,d,wnum,jmax;T_weight wei, wli[TWDIM];} 	T_aux;

void T_Chon(int,int,int,int,T_aux*); /* i, {-fermat,0=closed,+open}, nmax, g */
void T_Addweight(T_weight,T_aux *X);
int PPT_Check(T_weight nli,T_aux *X);

void Make_Trans_Weights(int n,int dmin,int dmax /*,int rFlag */)
{    int i,j,inc=1; T_aux X; X.n=n; X.wnum=0;
     outFILE=stdout; assert(n<=AMBI_Dmax); X.wei[0]=n;
     if(CHAT){ assert(CHAT==3); if(!(n%2)) {inc++; dmin+=(dmin%2);}}
      for(X.d=dmin;X.d<=dmax;X.d+=inc) 
      { X.wei[n+1]=X.d; X.wnum=0;
	if(ALLOWHALF) X.jmax=X.d/2; else X.jmax=(X.d-1)/2; 
	if(CHAT) T_Chon(1,-1,(X.d*(n- CHAT))/2-n+1,X.d,&X);  
	else T_Chon(1,-1, X.d -n+1,X.d,&X);
        for(i=0;i<X.wnum;i++) 
        {   Weight W; W.N=n; W.d=X.d; for(j=0;j<n;j++) W.w[j]=X.wli[i][j+1];
            W.M=0; if(Trans_Check(W)) Write_Weight(&W);       
        }
        fflush(outFILE);
      }
}
/* T_Chon chooses 0<wei[i]<=nm=d(n-3)/2-n+i-wei[1]-...-wei[i-1] such that a *
 * pointer structur with at most one unresolved pointer(urp) pointing to the*
 * right at a time is not ruled out. If (urp!=0) then the number of the un- *
 * resolved pointer is urp; in addition put all fermats to the left (urp<0).*
 * i has to point at l>=urp                                                 */
/* let j run; check mod(d||(d-n),j); if (upr) check if upr is resolved by j;*/
void T_Chon(int i, int urp, int nm, int g,T_aux *X)
{    int res, j, l=0, ip=i+1, jm=min(nm,X->jmax);
     if (i<X->n) for(j=(i==X->n-1) ? (1+nm-jm) : 1;j<=jm;j++)
		{X->wei[i]=j;/*next step*/
        if(urp<0) {if(X->d%j) {res=1; /* i.e. not ferm; res=0 -> resolved */
                      for(l=1;(l<i)&&res;l++) res=(X->d-X->wei[l])%j;
                      if (res) T_Chon(ip,i,nm-j+1,Fgcd(g,j),X);
                      else T_Chon(ip,0,nm-j+1,Fgcd(g,j),X);}
                 else{if((i==1)||(j>=X->wei[i-1]))
			T_Chon(ip,-1,nm-j+1,Fgcd(g,j),X);}
        } else  /* continue;} */
        if(X->d%j) {                   /* now there can be no more fermat */
        l=max(urp,1); for(res=1;(l<i)&&res;l++) res=(X->d-X->wei[l])%j;
        if(urp){if(!res){if((X->d-j)%X->wei[urp]) 
			T_Chon(ip,urp,nm-j+1,Fgcd(g,j),X);
                            else T_Chon(ip,0,nm-j+1,Fgcd(g,j),X);}}
        else {if(res) T_Chon(ip,i,nm-j+1,Fgcd(g,j),X);
              else T_Chon(ip,0,nm-j+1,Fgcd(g,j),X);} }}
     else {if(X->jmax<nm) return;                              /* last step */
        if(urp>0) if((X->d-nm)%X->wei[urp]) return;    /* pointer resolved? */
        if((res=(X->d%nm))||(urp<0)){ 
            for(l=1;(l<i)&&res;l++) res=(X->d-X->wei[l])%nm;
        if(!res){X->wei[i]=nm;if(Fgcd(g,nm)==1)PPT_Check(X->wei,X);}} }
}
/*  ppcheck checks whether the formal poincare polynomial is a polynomial  */
int PPT_Check(T_weight nli,T_aux *X)
{    int i=0, n, t, tt, j,d=X->d; assert(d==nli[nli[0]]);
     for (i=1;i<=nli[0];i++){
        n=1; tt=nli[i];
        for (j=i+1;j<=nli[0];j++) if (!(nli[j]%tt)) n++;
        for (j=1;j<=nli[0];j++) if (!((d-nli[j])%tt)) n--;
        if (n>0) return 0;
        for (t=2;t*t<=nli[i];t++)
           if (!mod(nli[i],t)) {
              n=1;
              for (j=i+1;j<=nli[0];j++) if (!mod(nli[j],t)) n++;
              for (j=1;j<=nli[0];j++) if (!mod(d-nli[j],t)) n--;
              if (n>0) return 0;
              n=1; tt=nli[i]/t;
              for (j=i+1;j<=nli[0];j++) if (!mod(nli[j],tt)) n++;
              for (j=1;j<=nli[0];j++) if (!mod(d-nli[j],tt)) n--;
              if (n>0) return 0;};};
     T_Addweight(nli,X); return 1;
}
int  T_Weicomp(T_weight w1,T_weight w2)/* w2-w1,i.e.pos if w1<w2,neg if w1>w2*/
{    int i=1; while((i<=(*w1)) && (w1[i]==w2[i])) i++; return w2[i]-w1[i];
}
void T_Insertat(T_weight ww, int position,T_aux *X)
{    int i, j;
     for(i=X->wnum-1;i>=position;i--)
         for(j=0;j<X->wli[i][0]+2;j++) X->wli[i+1][j]=X->wli[i][j];
     for(j=0;j<ww[0]+2;j++) X->wli[position][j]=ww[j];
     X->wnum++;
}
void T_Addweight(T_weight win,T_aux *X)
{    int i, p, n0, n1; T_weight wn; for(i=0;i<*win+2;i++) wn[i]=win[i];
     for(i=1;i<wn[0];i++) for(p=i+1;p<=wn[0];p++)
             if (wn[i]>wn[p]) swap(&wn[i],&wn[p]);  /* make n0<=n1<=...<=n# */
     if (X->wnum) {if ((i = T_Weicomp(wn,X->wli[n0=0])))
                   {if (i>0) {T_Insertat(wn,0,X); return;}} else return;
                if ((i = T_Weicomp(wn,X->wli[n1=X->wnum-1])))
                   {if (i<0) {T_Insertat(wn,X->wnum,X); return;}} else return;
        while(n1>n0+1) {p=(n0+n1)/2; i=T_Weicomp(wn,X->wli[p]);
        if(i) {if(i>0) n1=p; else n0=p;} else return;}
        T_Insertat(wn,n1,X);}
     else T_Insertat(wn,0,X);
}
/*  ==========  	  End of MAKE WEIGHTS d>4:		==========  */




/*  ==========       	    ALL  CWS  in  d <= 4		==========  */
typedef struct {int d, w[2];}   wei2; 
typedef struct {int d, w[3];}   wei3; 
typedef struct {int d, w[4];}   wei4;
const wei2 W2={2,{1,1}};
const wei3 W3[ 3]={{3,{1,1,1}},   {4,{1,1,2}},      {6,{1,2,3}}};
const wei4 W4[95]={
{4, {1,1,1,1}},	 {5,{1,1,1,2}},	{6,{1,1,2,2}},	{6,{1,1,1,3}},	{7,{1,1,2,3}},
{8, {1,2,2,3}},	 {8,{1,1,2,4}},	{9,{1,2,3,3}},	{9,{1,1,3,4}},	{10,{1,2,3,4}},
{10,{1,2,2,5}},	 {10,{1,1,3,5}},{11,{1,2,3,5}},	{12,{2,3,3,4}},	{12,{1,3,4,4}},
{12,{2,2,3,5}},	 {12,{1,2,4,5}},{12,{1,2,3,6}},	{12,{1,1,4,6}},	{13,{1,3,4,5}},
{14,{2,3,4,5}},	 {14,{2,2,3,7}},{14,{1,2,4,7}},	{15,{3,3,4,5}},	{15,{2,3,5,5}},
{15,{1,3,5,6}},	 {15,{1,3,4,7}},{15,{1,2,5,7}},	{16,{1,4,5,6}},	{16,{2,3,4,7}},
{16,{1,3,4,8}},	 {16,{1,2,5,8}},{17,{2,3,5,7}},	{18,{3,4,5,6}},	{18,{1,4,6,7}},
{18,{2,3,5,8}},	 {18,{2,3,4,9}},{18,{1,3,5,9}},	{18,{1,2,6,9}},	{19,{3,4,5,7}},
{20,{2,5,6,7}},	 {20,{3,4,5,8}},{20,{2,4,5,9}},	{20,{2,3,5,10}},{20,{1,4,5,10}},
{21,{3,5,6,7}},	 {21,{1,5,7,8}},{21,{2,3,7,9}},	{21,{1,3,7,10}},{22,{2,4,5,11}},
{22,{1,4,6,11}}, {22,{1,3,7,11}},{24,{3,6,7,8}},{24,{4,5,6,9}},	{24,{1,6,8,9}},
{24,{3,4,7,10}},{24,{2,3,8,11}},{24,{3,4,5,12}},{24,{2,3,7,12}},{24,{1,3,8,12}},
{25,{4,5,7,9}},	{26,{2,5,6,13}},{26,{1,5,7,13}},{26,{2,3,8,13}},{27,{5,6,7,9}},
{27,{2,5,9,11}},{28,{4,6,7,11}},{28,{3,4,7,14}},{28,{1,4,9,14}},{30,{5,6,8,11}},
{30,{3,4,10,13}}, {30,{4,5,6,15}}, {30,{2,6,7,15}}, {30,{1,6,8,15}},
{30,{2,3,10,15}}, {30,{1,4,10,15}},{32,{4,5,7,16}}, {32,{2,5,9,16}},
{33,{3,5,11,14}}, {34,{4,6,7,17}}, {34,{3,4,10,17}},{36,{7,8,9,12}},
{36,{3,4,11,18}}, {36,{1,5,12,18}},{38,{5,6,8,19}}, {38,{3,5,11,19}},
{40,{5,7,8,20}},  {42,{3,4,14,21}},{42,{2,5,14,21}},{42,{1,6,14,21}},
{44,{4,5,13,22}}, {48,{3,5,16,24}},{50,{7,8,10,25}},{54,{4,5,18,27}},
{66,{5,6,22,33}}};

void MakeSelections(FILE *, FILE *, int);
void Make2CWS(FILE *, FILE *, int, int);
void RW_TO_CWS(CWS *, Weight *, int, int, int, int);
void W_TO_CWS(CWS *, Weight *, int, int, int, int);
void PRINT_CWS(CWS *);
void Make_111_CWS(FILE **, int *);
void Make_nno_CWS(FILE **, int, int);

void STtmp(FILE *w2FILE, FILE *w3FILE, FILE *w4FILE)
{
  int i,j;

  for(i=0; i<95; i++){fprintf(w4FILE, "%d  ",W4[i].d);
    for(j=0; j<4; j++) fprintf(w4FILE, "%d ",W4[i].w[j]);
    fprintf(w4FILE,"\n");
  }
  for(i=0; i<3; i++){fprintf(w3FILE, "%d  ",W3[i].d);
    for(j=0; j<3; j++) fprintf(w3FILE, "%d ",W3[i].w[j]);
    fprintf(w3FILE,"\n");
  }
  fprintf(w2FILE, "%d  ",W2.d); for(j=0;j<2;j++)fprintf(w2FILE, "%d ",W2.w[j]);
  fprintf(w2FILE,"\n");
  rewind(w2FILE); rewind(w3FILE); rewind(w4FILE);
}

void mkold2(char *outfile, FILE *INFILE1, FILE *INFILE2, int u, int ef)
{
  
  FILE *AUXFILE1, *AUXFILE2;
  
  if((AUXFILE1 = tmpfile()) == NULL) 
	Die("Unable to open tmpfile for read/write");
  if((AUXFILE2 = tmpfile()) == NULL) 
	Die("Unable to open tmpfile for read/write");

  MakeSelections(INFILE1, AUXFILE1, u);
  MakeSelections(INFILE2, AUXFILE2, u);

  if(strcmp(outfile,""))
    if((outFILE = fopen(outfile, "w")) == NULL){
      printf("\nUnable to open file %s for write\n",outfile); exit(0);}
 
  Make2CWS(AUXFILE1, AUXFILE2, u, ef);
  fclose(AUXFILE1); fclose(AUXFILE2);
  if(strcmp(outfile,""))
    fclose(outFILE);
}

void mk2xxx(char *outfile, int n)
{
  int i,j,d=0;
  CWS CW; CW.nz=0;

  CW.N=2*n; CW.nw=n; 
  for(i=0; i<n; i++){
    CW.d[i]=2; 
    for(j=0; j<CW.N; j++){
      if((j==d) || (j==d+1)) CW.W[i][j]=1;
      else CW.W[i][j]=0;
    }
    d=d+2; 
  }
  if(strcmp(outfile,""))
    if((outFILE = fopen(outfile, "w")) == NULL){
      printf("\nUnable to open file %s for write\n",outfile); exit(0);}
  PRINT_CWS(&CW);
  if(strcmp(outfile,""))
    fclose(outFILE);
}

void mk3u3u3(char *outfile, FILE *INFILE)
{
 
  FILE *AUXFILE[3] = {NULL};
  int i, u = 1, eq[2]; 

  for(i = 0; i < 3; i++){
    if((AUXFILE[i] = tmpfile()) == NULL) 
	Die("Unable to open tmpfile for read/write");
    MakeSelections(INFILE, AUXFILE[i], u);
  }
  if(strcmp(outfile,""))
    if((outFILE = fopen(outfile, "w")) == NULL){
      printf("\nUnable to open file %s for write\n",outfile); exit(0);}
 
  eq[0] = eq[1] = 1;
  Make_111_CWS(AUXFILE, eq);
  for(i = 0; i < 3; i++) fclose(AUXFILE[i]); 
  if(strcmp(outfile,""))
    fclose(outFILE);
}

void mkold_nno(char *outfile, FILE *INFILE1, FILE *INFILE2, FILE *INFILE3, 
	int u, int eq)
{
  FILE *AUXFILE[3]; 
  
  if((AUXFILE[0]=tmpfile())==NULL)Die("Unable to open tmpfile for read/write");
  if((AUXFILE[1]=tmpfile())==NULL)Die("Unable to open tmpfile for read/write");

  MakeSelections(INFILE1, AUXFILE[0], u);
  MakeSelections(INFILE2, AUXFILE[1], u);

  AUXFILE[2] = INFILE3;
  if(strcmp(outfile,""))
    if((outFILE = fopen(outfile, "w")) == NULL){
      printf("\nUnable to open file %s for write\n",outfile); exit(0);}
  
  Make_nno_CWS(AUXFILE, u, eq);
  fclose(AUXFILE[0]); fclose(AUXFILE[1]); 
  if(strcmp(outfile,""))
    fclose(outFILE);
}

void Make_34_CWS(int d)
{    
  FILE *w2FILE, *w3FILE, *w4FILE; 
  int u, ef;
  char *outfile;

  assert(d<=4); /*puts("Implement Make_34_CWS");*/
  if((w2FILE = tmpfile())==NULL) Die("Unable to open tmpfile for read/write");
  if((w3FILE = tmpfile())==NULL) Die("Unable to open tmpfile for read/write");
  if((w4FILE = tmpfile())==NULL) Die("Unable to open tmpfile for read/write");
  
  STtmp(w2FILE, w3FILE, w4FILE);
  
  if(d == 3){outfile="";
  /*outfile="3u3"*/; u = 1; ef = 1; mkold2(outfile, w3FILE, w3FILE, u, ef);
  /*outfile="2x3"*/; u = 0; ef = 0; mkold2(outfile, w2FILE, w3FILE, u, ef);
  /*outfile="2x2x2"*/; mk2xxx(outfile, atoi("3"));
  }
  if(d == 4){ outfile="";
    /*outfile="4uu4";*/ u = 2; ef = 1; mkold2(outfile, w4FILE, w4FILE, u, ef);
    /*outfile="3u4"; */ u = 1; ef = 0; mkold2(outfile, w3FILE, w4FILE, u, ef);
    /*outfile="3x3"; */ u = 0; ef = 1; mkold2(outfile, w3FILE, w3FILE, u, ef);
    /*outfile="2x4"; */ u = 0; ef = 0; mkold2(outfile, w2FILE, w4FILE, u, ef);
    /*outfile="3u3x2";*/ u=1;ef=1;mkold_nno(outfile,w3FILE,w3FILE,w2FILE,u,ef);
    /*outfile="3x2x2";*/ u=0;ef=0;mkold_nno(outfile,w3FILE,w2FILE,w2FILE,u,ef);
    /*outfile="3u3u3";*/ 		mk3u3u3(outfile, w3FILE);
    /*outfile="2x2x2x2";*/ 		mk2xxx(outfile, atoi("4"));
  }
  fclose(w2FILE);fclose(w3FILE);fclose(w4FILE);
}
/*  ==========          End of  ALL  CWS  in  d <= 4    	==========  */

/*  ==========  	        MAKE CWS d>4               	==========  */

void Print_CWS_Zinfo(CWS *CW);
void Print_CWS(CWS *_W)
{
  int i, j;
  
  for (i = 0; i < _W->nw; i++) {
    fprintf(outFILE, "%d ", (int) _W->d[i]);
    for (j = 0; j < _W->N; j++)
      fprintf(outFILE, "%d ", (int) _W->W[i][j]);
    if (i + 1 < _W->nw)
      fprintf(outFILE, " ");
  }  Print_CWS_Zinfo(_W);
}

void print_W(Weight *_s, Weight *_W, FILE *auxFILE){

  int i,j=0;

  fprintf(auxFILE, "%d ", (int) _W->d);
  for(i = 0; i < _s->N; i++)
    fprintf(auxFILE, "%d ", (int) _W->w[_s->w[i]]);
  for(i = 0; i < _W->N; i++)
    if(i != _s->w[j])
      fprintf(auxFILE, "%d ", (int) _W->w[i]);
    else
      if(j < (_s->N - 1))
        j++;
  fprintf(auxFILE, "\n");
}

void next_n(Weight *_s, Weight *_W, int *_n, FILE *auxFILE){
  int i;

  if(_s->N == *_n)
    print_W(_s, _W, auxFILE);
  else
    if(_s->w[_s->N - 1] != (_W->N - 1)){
      if (_W->w[_s->w[_s->N - 1] + 1] == _W->w[_s->w[_s->N - 1]]){
	_s->w[_s->N] = _s->w[_s->N - 1] + 1;
	_s->N ++;
	next_n(_s, _W, _n, auxFILE);
	_s->N --;
      }
      for(i = _s->w[_s->N - 1] + 1; i < _W->N; i++){
	if(_W->w[i] > _W->w[i - 1]){
	  _s->w[_s->N] =  i;
	  _s->N ++;
	  next_n(_s, _W, _n, auxFILE);
	  _s->N --;
	}
      }
    }
}

void Select_n_of_W(Weight *_W, int n, FILE *auxFILE){

  int i;
  Weight s;
  
  if(n == 0){
    fprintf(auxFILE, "%d ", (int) _W->d);
    for(i = 0; i < _W->N; i++)
      fprintf(auxFILE, "%d ", (int) _W->w[i]);
    fprintf(auxFILE, "\n");
  }
  else{
    for(i = 1; i < _W->N; i++)
      if(_W->w[i] < _W->w[i-1])
	Die("Weights must be sortet: W_1 <= W_2 <= .... <=W_N!");
    s.N = 1; s.d = _W->d;
    s.w[0] = 0;
    next_n(&s, _W, &n, auxFILE);
    for(i = 1; i < _W->N; i++)
      if(_W->w[i] > _W->w[i-1]){
	s.w[0] = i;
	next_n(&s, _W, &n, auxFILE);
      }
  }
}

void PRINT_CWS(CWS *CW){
#if (!Only_IP_CWS)
  {
    Print_CWS(CW);
    fprintf(outFILE,"\n");
  }
#else
  {
    PolyPointList *P, *DP; EqList E; VertexNumList V;
    P = (PolyPointList *) malloc(sizeof(PolyPointList));
    if (P == NULL) Die("Unable to allocate space for P");
    DP = (PolyPointList *) malloc(sizeof(PolyPointList));
    if (DP == NULL) Die("Unable to allocate space for DP");
    CW->index = 1;
    Make_CWS_Points(CW, P);
    if (IP_Check(P,&V,&E)){
      int r=1,i=-1;
      Print_CWS(CW);
      while(r && (++i < E.ne))
	if(E.e[i].c != 1) r = 0;
      Make_Dual_Poly(P, &V, &E, DP);
      fprintf(outFILE," M:%d %d",P->np, V.nv);
      if(r) fprintf(outFILE," N:%d %d",DP->np, E.ne);
      else  fprintf(outFILE," F:%d N:%d", E.ne,DP->np);
      assert(IP_Check(DP,&V,&E));
      fprintf(outFILE,"\n");
    }
    free(DP); free(P); 
  }
#endif  
}

void W_TO_CWS(CWS *CW, Weight *_W, int Nf, int Nr, int Nb, int u)
{  /* d (Nf x 0) w_1 ... w_u  (Nb x 0)  w_(u+1) ... w_N  (Nr x 0) */
  int i, j, d = 0;
    CW->d[CW->nw]=_W->d; CW->N=(Nf+Nb+Nr+_W->N);
  if(CW->N > AMBI_Dmax) Die("increase AMBI_Dmax !");
  for(i=0; i < Nf; i++) CW->W[CW->nw][i]=0;
  for(i=0; i < _W->N; i++){
    CW->W[CW->nw][i+Nf+d]=_W->w[i];
    if(u && (i == (u - 1))){ 
      for(j=0; j < Nb; j++)
	CW->W[CW->nw][i+j+Nf+1]=0;
      d=Nb;
    }
  }
  for(i=0; i < Nr; i++) CW->W[CW->nw][Nf+d+_W->N+i]=0;
  CW->nw++; CW->nz=0;
}

void RW_TO_CWS(CWS *CW, Weight *_W, int Nf, int Nr, int Nb, int u)
{  /* d (Nf x 0) w_N ... w_(u+1)  (Nb x 0)  w_u ... w_1  (Nr x 0) */
  int i, j, d = 0;  
  CW->d[CW->nw]=_W->d; CW->N=(Nf+Nb+Nr+_W->N);
  if(CW->N > AMBI_Dmax) Die("increase AMBI_Dmax !");
  for(i=0; i < Nf; i++) CW->W[CW->nw][i]=0;
  for(i=0; i < _W->N; i++){
    CW->W[CW->nw][i+Nf+d]=_W->w[_W->N - i - 1];
    if(u && (i == (u - 1))){ 
      for(j=0; j < Nb; j++)
        CW->W[CW->nw][i+j+Nf+1]=0;
      d=Nb;
    }
  }
  for(i=0; i < Nr; i++) CW->W[CW->nw][Nf+d+_W->N+i]=0;
  CW->nw++; CW->nz=0;
}

void SWAP(Long *_A, Long *_B)
{
  Long C;
  C = *_A; *_A = *_B; *_B = C;
}

void scan_dim(int nF, char *infile[], int D[])
{
  int i, j = 0;
  FILE *INfile[NFmax];
  Weight W;

  for(i = 0; i < nF; i++){
    if((INfile[i] = fopen(infile[i], "r"))== NULL){
      printf("\nUnable to open file %s for read\n",infile[i]);exit(0);}
    j = 0; while(READ_Weight(&W, INfile[i]))if (j++) break;
    D[i] = W.N - 1;
    fclose(INfile[i]);
  }
}

int Wcomp(Weight *_W1, Weight *_W2){

  int i,j;
  int A[POLY_Dmax+1], B[POLY_Dmax+1];

  if (_W1->N > (POLY_Dmax+1)) Die("increase POLY_Dmax!");
  if (_W1->N != _W2->N) Die("N1 != N2 in Wcomp!");
      
  for(i = 0; i < _W1->N; i++){
    A[i] = i;
    B[i] = i;
  }
  for(i = 0; i < _W1->N - 1; ++i)
    for(j = _W1->N - 1; j > i; --j){
      if(_W1->w[A[j-1]] >  _W1->w[A[j]])
        swap(&A[j-1], &A[j]);
      if(_W2->w[B[j-1]] >  _W2->w[B[j]])
        swap(&B[j-1], &B[j]);
    }
  for(i = 0; i < _W1->N; i++){
    if(_W1->w[A[i]] >  _W2->w[B[i]])
      return 1;
    if(_W1->w[A[i]] <  _W2->w[B[i]])
      return -1;
  }
  return 0;
}

void Make_nno_CWS(FILE *AUXFILE[], int u, int ef)
{
  int n=0, l[2]; 
  Weight W[3];
  CWS CW;

  l[0]=0; while (READ_Weight(&W[0], AUXFILE[0])){
    l[1]=0; l[0]++; while (READ_Weight(&W[1], AUXFILE[1])){
      l[1]++; while (READ_Weight(&W[2], AUXFILE[2])){
	if((l[0] <= l[1]) || !ef){ CW.nw=0;
	  RW_TO_CWS(&CW, &W[0], n, (W[1].N + W[2].N - u), n, n);
	  W_TO_CWS(&CW, &W[1], (W[0].N - u), W[2].N, n, n);
	  W_TO_CWS(&CW, &W[2], (W[0].N + W[1].N - u), n, n, n);
	  PRINT_CWS(&CW);
	  if(u == 2)
	    if((W[0].w[0] != W[0].w[1]) && (W[1].w[0] != W[1].w[1])){
	      SWAP(&W[1].w[0], &W[1].w[1]); CW.nw=0;
	      RW_TO_CWS(&CW, &W[0], n, (W[1].N + W[2].N - u), n, n);
	      W_TO_CWS(&CW, &W[1], (W[0].N - u), W[2].N, n, n);
	      W_TO_CWS(&CW, &W[2], (W[0].N + W[1].N), n, n, n);
	      PRINT_CWS(&CW);
	    }
	}
      }
      rewind(AUXFILE[2]);
    }
    rewind(AUXFILE[1]);
  }
  rewind(AUXFILE[0]); 
}

void Make_111_CWS(FILE *AUXFILE[], int ef[])
{
  int u = 1, n=0, l[3]; 
  Weight W[3]; 
  CWS CW;

  l[0]=0; while (READ_Weight(&W[0], AUXFILE[0])){
    l[1]=0; l[0]++; while (READ_Weight(&W[1], AUXFILE[1])){
      l[2]=0; l[1]++; while (READ_Weight(&W[2], AUXFILE[2])){
	l[2]++; 
	if(((l[0] <= l[1]) || !ef[0]) && ((l[1] <= l[2]) || !ef[1])){
	  CW.nw=0;
	  RW_TO_CWS(&CW, &W[0], n, (W[1].N + W[2].N - 2*u), n, n);
	  W_TO_CWS(&CW, &W[1], (W[0].N - u), (W[2].N - u), n, n);
	  W_TO_CWS(&CW, &W[2], (W[0].N - u), n, (W[1].N - u), u);
	  PRINT_CWS(&CW);
	}
      }
      rewind(AUXFILE[2]);
    }
    rewind(AUXFILE[1]);
  }
  rewind(AUXFILE[0]);
}

void Make_222_CWS(FILE *AUXFILE[], int ef[])
{
  int u = 2, n=0, l[3]; 
  Weight W[3]; 
  CWS CW;

  l[0]=0; while (READ_Weight(&W[0], AUXFILE[0])){
    l[1]=0; l[0]++; while (READ_Weight(&W[1], AUXFILE[1])){
      l[2]=0; l[1]++; while (READ_Weight(&W[2], AUXFILE[2])){
	l[2]++; 
	if(((l[0] <= l[1]) || !ef[0]) && ((l[1] <= l[2]) || !ef[1])){
	  CW.nw=0;
	  RW_TO_CWS(&CW, &W[0], n, (W[1].N + W[2].N - 2*u), n, n);
	  W_TO_CWS(&CW, &W[1], (W[0].N - u), (W[2].N - u), n, n);
	  W_TO_CWS(&CW, &W[2], (W[0].N - u), n, (W[1].N - u), u);
	  PRINT_CWS(&CW);
	  if(W[0].w[0] == W[0].w[1]){
	    if((W[1].w[0] != W[1].w[1]) && (W[2].w[0] != W[2].w[1])){
	      SWAP(&W[1].w[0], &W[1].w[1]); CW.nw=0;
	      RW_TO_CWS(&CW, &W[0], n, (W[1].N + W[2].N - 2*u), n, n);
	      W_TO_CWS(&CW, &W[1], (W[0].N - u), (W[2].N - u), n, n);
	      W_TO_CWS(&CW, &W[2], (W[0].N - u), n, (W[1].N - u), u);
	      PRINT_CWS(&CW);
	    }
	  }
	  else{
	    if((W[1].w[0] != W[1].w[1]) || (W[2].w[0] != W[2].w[1])){
	      SWAP(&W[0].w[0], &W[0].w[1]); CW.nw=0;
	      RW_TO_CWS(&CW, &W[0], n, (W[1].N + W[2].N - 2*u), n, n);
	      W_TO_CWS(&CW, &W[1], (W[0].N - u), (W[2].N - u), n, n);
	      W_TO_CWS(&CW, &W[2], (W[0].N - u), n, (W[1].N - u), u);
	      PRINT_CWS(&CW);
	    }
	    if((W[1].w[0] != W[1].w[1]) && (W[2].w[0] != W[2].w[1])){
	      SWAP(&W[1].w[0], &W[1].w[1]); CW.nw=0;
	      RW_TO_CWS(&CW, &W[0], n, (W[1].N + W[2].N - 2*u), n, n);
	      W_TO_CWS(&CW, &W[1], (W[0].N - u), (W[2].N - u), n, n);
	      W_TO_CWS(&CW, &W[2], (W[0].N - u), n, (W[1].N - u), u);
	      PRINT_CWS(&CW);
	    }
	  }
	}
      }
      rewind(AUXFILE[2]);
    }
    rewind(AUXFILE[1]);
  }
  rewind(AUXFILE[0]);
}

void Make_221_CWS(FILE *AUXFILE[], int ef)
{
  int n=0, U = 2, u = 1, i=0, l[2]; 
  Weight W[3]; 
  CWS CW;

  l[0]=0; while (READ_Weight(&W[0], AUXFILE[0])){
    l[1]=0; l[0]++; while (READ_Weight(&W[1], AUXFILE[1])){
      l[1]++; while (READ_Weight(&W[2], AUXFILE[2])){
	if((l[0] <= l[1]) || !ef){
	  i=0; CW.nw=0;
	  RW_TO_CWS(&CW, &W[0], n, (W[1].N + W[2].N - U -u ), n, n);
	  W_TO_CWS(&CW, &W[1], (W[0].N - U), (W[2].N - u ), n, n);
	  W_TO_CWS(&CW, &W[2], (W[0].N - U + i), n, (W[1].N - U + 1 - i), u);
	  PRINT_CWS(&CW);
	  if((W[0].w[0] != W[0].w[1]) && (W[1].w[0] != W[1].w[1])){
	    SWAP(&W[1].w[0], &W[1].w[1]); CW.nw=0;
	    RW_TO_CWS(&CW, &W[0], n, (W[1].N + W[2].N - U -u ), n, n);
	    W_TO_CWS(&CW, &W[1], (W[0].N - U), (W[2].N - u ), n, n);
	    W_TO_CWS(&CW, &W[2], (W[0].N - U + i), n, (W[1].N - U + 1 - i), u);
	    PRINT_CWS(&CW);
	  }
	  if((W[0].w[0] != W[0].w[1]) || (W[1].w[0] != W[1].w[1])){
	    i=1; CW.nw=0;
	    RW_TO_CWS(&CW, &W[0], n, (W[1].N + W[2].N - U -u ), n, n);
	    W_TO_CWS(&CW, &W[1], (W[0].N - U), (W[2].N - u ), n, n);
	    W_TO_CWS(&CW, &W[2], (W[0].N - U + i), n, (W[1].N - U + 1 - i), u);
	    PRINT_CWS(&CW);
	    if((W[0].w[0] != W[0].w[1]) && (W[1].w[0] != W[1].w[1])){
	      SWAP(&W[1].w[0], &W[1].w[1]); CW.nw=0;
	      RW_TO_CWS(&CW, &W[0], n, (W[1].N + W[2].N - U -u ), n, n);
	      W_TO_CWS(&CW, &W[1], (W[0].N - U), (W[2].N - u ), n, n);
	      W_TO_CWS(&CW, &W[2], (W[0].N - U + i), n,(W[1].N - U + 1 - i),u);
	      PRINT_CWS(&CW);
	    }
	  }
	}
      }
      rewind(AUXFILE[2]);
    }
    rewind(AUXFILE[1]);
  }
  rewind(AUXFILE[0]); 
}

void Make_211_CWS(FILE *AUXFILE[], int ef)
{
  int u = 1, n=0, l[3]; 
  Weight W[3]; 
  CWS CW;

  while (READ_Weight(&W[0], AUXFILE[0])){
    l[1]=0; while (READ_Weight(&W[1], AUXFILE[1])){
      l[2]=0; l[1]++; while (READ_Weight(&W[2], AUXFILE[2])){
	l[2]++; 
	if(((l[1] <= l[2]) || !ef)){ CW.nw=0;
	  RW_TO_CWS(&CW, &W[0], n, (W[1].N + W[2].N - 2*u), n, n);
	  W_TO_CWS(&CW, &W[1], (W[0].N - u), (W[2].N - u), n, n);
	  W_TO_CWS(&CW, &W[2], (W[0].N - 2*u), n, W[1].N, u);
	  PRINT_CWS(&CW);
	  if((W[0].w[0] != W[0].w[1]) && (W[1].w[0] !=  W[2].w[0])){
	    SWAP(&W[0].w[0], &W[0].w[1]); CW.nw=0;
	    RW_TO_CWS(&CW, &W[0], n, (W[1].N + W[2].N - 2*u), n, n);
	    W_TO_CWS(&CW, &W[1], (W[0].N - u), (W[2].N - u), n, n);
	    W_TO_CWS(&CW, &W[2], (W[0].N - 2*u), n, W[1].N, u);
	    PRINT_CWS(&CW);
	  }
	}
      }
      rewind(AUXFILE[2]);
    }
    rewind(AUXFILE[1]);
  }
  rewind(AUXFILE[0]);
}

void Make2CWS(FILE *AUXFILE1, FILE *AUXFILE2, int u, int ef){

  Weight W1, W2;
  int n=0, l[2]; 
  CWS CW;

  if (u > 2) Die("for u > 2 no support !");
    l[0]=0;while (READ_Weight(&W1, AUXFILE1)){
      l[1]=0; l[0]++;
      while (READ_Weight(&W2, AUXFILE2)){
	l[1]++;
	if((l[1] >= l[0]) || !ef){CW.nw=0;
	  RW_TO_CWS(&CW, &W1, n, (W2.N-u), n, n);
	  W_TO_CWS(&CW, &W2, (W1.N - u), n, n, n);
	  PRINT_CWS(&CW);
	  if(u == 2)
	    if((W1.w[0] != W1.w[1]) && (W2.w[0] != W2.w[1])){
	      SWAP(&W2.w[0], &W2.w[1]); CW.nw=0;
	      RW_TO_CWS(&CW, &W1, n, (W2.N - u), n, n);
	      W_TO_CWS(&CW, &W2, (W1.N - u), n, n, n);
	      PRINT_CWS(&CW);
	    }
	}
      }
      rewind(AUXFILE2);
    }
    rewind(AUXFILE1);
}

void MakeSelections(FILE *INFILE, FILE *AUXFILE, int u)
{
  Weight W;

  while (READ_Weight(&W, INFILE))
    Select_n_of_W(&W, u, AUXFILE);
  rewind(AUXFILE); rewind(INFILE);
}

void PrintCWSTypes(void)
{
  const char B[]="         ";
  printf("\nThe following types are available:\n\n");
  printf("#infiles = 2 (need no -t option):\n");
  printf("%s-c# -n2 [intile1] [infile2] (-t 0 0)\n",B);
  printf("%s-c# -n2 [intile1] [infile2] (-t 1 1)\n",B);
  printf("%s-c# -n2 [intile1] [infile2] (-t 2 2)\n",B);
  printf("#infiles = 3:\n");
  printf("%s-c# -n3 [intile1] [infile2] [infile3] -t n n 0\n",B);
  printf("%s-c# -n3 [intile1] [infile2] [infile3] -t 1 1 1\n",B);
  printf("%s-c# -n3 [intile1] [infile2] [infile3] -t 2 1 1\n",B);
  printf("%s-c# -n3 [intile1] [infile2] [infile3] -t 2 2 1\n",B);
  printf("%s-c# -n3 [intile1] [infile2] [infile3] -t 2 2 2\n",B);
  exit(0);
}

void Make_IP_CWS(int narg, char* fn[])
{
  FILE *INFILE[NFmax] = {NULL}, *AUXFILE[NFmax] = {NULL};
  char *infile[NFmax] = {NULL}, *outfile = NULL, *a;
  int n = 0, d = 0, u = -1, nF = 0, i, D[NFmax];
  CWS_type t;

  t.nu = 0;
  for (i=0; i<NFmax; i++) t.u[i]=0;

  while (narg > ++n) {
    if(fn[n][0] != '-') 
      break;
    if(fn[n][1] == 'c'){
      if((fn[n][2]==0) && (narg>n+1)) a=fn[++n]; else a=&fn[n][2];
      if(!IsDigit(*a)) Die("after -c there must be a digit!");
      d = atoi(a);
    }
    if(fn[n][1] == 'n'){ 
      if((fn[n][2]==0) && (narg>n+1)) a=fn[++n]; else a=&fn[n][2];
      if(!IsDigit(*a)) Die("after -n there must be a digit!");
      nF = atoi(a);
      n++; break;
    }
    if(fn[n][1] == 't') PrintCWSTypes();
  }
  for(i = 0; i < nF; i++){
    if((n >= narg)||(fn[n][0] == '-'))
	{printf("#infiles = %d < %d!\n",i,nF); exit(0);}
    infile[i] = fn[n];
    n++;
  }
  if((narg > n) && (fn[n][0] != '-')) Die("too many infiles!");
  n--; t.nu=0;
  while (narg > ++n) {
    if(fn[n][0] != '-') break;
    if(fn[n][1] == 't'){ 
      if((fn[n][2]==0) && (narg>n+1)) a=fn[++n]; else a=&fn[n][2];
      if(*a == 0) PrintCWSTypes();  
      if(!IsDigit(*a)) Die("after -t there must be digit(s)!");  
      t.u[0] = atoi(a); t.nu = 1;
      while ((narg > ++n)  && (t.nu < nF)) {
	if(fn[n][0] == '-') break;
	a=fn[n]; 
	if(!IsDigit(*a)) Die("after -t there must be digit(s)!");
	assert(t.nu < NFmax); t.u[t.nu] = atoi(a); t.nu++;
      }
      n--;
      if(narg > ++n) outfile = fn[n]; 
    }
  }
  if(nF == 0) Die("there is no -n#infiles!"); 
  if(d == 0) Die("No dimensoin specified!");
  if(t.nu && (t.nu != nF)) 
    Die("if input is -nN -t k_1,...,k_n then N must be equal to n!");
  if(outfile == NULL) outFILE = stdout;
  else 
    if((outFILE = fopen(outfile, "w")) == NULL){
      printf("\nUnable to open file %s for write\n",fn[n]); exit(0);}
  scan_dim(nF, infile, D);
  for(i = 0; i < nF; i++){
    if((AUXFILE[i] = tmpfile()) == NULL) 
	Die("Unable to open tmpfile to read/write");
    if((INFILE[i] = fopen(infile[i], "r"))==NULL)
	Die("Unable to open infile to read");
  }
  if(nF == 2){
    if(!t.nu) assert((u = D[0] + D[1] -d) >= 0); 
    else{
      if(t.u[0] != t.u[1]) 
	Die("if input is -n2 -t k_1 k_2 then k_1 must be equal to k_2!");
      if(t.u[0] != (D[0] + D[1] - d))
	Die("wrong DIM -cDIM or wrong TYPES -t TYPE1 TYPE2");
      u = t.u[0];
    }
    for(i = 0; i < nF; i++)
      MakeSelections(INFILE[i], AUXFILE[i], u);
    Make2CWS(AUXFILE[0], AUXFILE[1], u, !strcmp(infile[0], infile[1]));
  }
  if(nF == 3){
    if(!t.nu) Die("with nNUMBER and NUMBER>2 I need -t TYPE1 TYPE2 TYPE3!");
    if(((t.u[0] == 1) && (t.u[1] == 1) && (t.u[2] == 1)) ||
       ((t.u[0] == 2) && (t.u[1] == 2) && (t.u[2] == 2))){ int eq[2];
       if((D[0] + D[1] + D[2] -2*t.u[0] - d) != 0)
	Die("wrong DIM -cDIM or wrong TYPES -t TYPE1 TYPE2 TYPE3");
      eq[0] = 0; eq[1] = 0;
      if(!strcmp(infile[0], infile[1])) eq[0] = 1;
      if(!strcmp(infile[1], infile[2])) eq[1] = 1;
      if(((!eq[0]) && (!eq[1])) && (!strcmp(infile[0], infile[2]))){
	eq[0] = 1;
	MakeSelections(INFILE[0], AUXFILE[0], t.u[0]);
	MakeSelections(INFILE[2], AUXFILE[1], t.u[1]);
	MakeSelections(INFILE[1], AUXFILE[2], t.u[2]);
      }
      else
	for(i = 0; i < nF; i++)
	  MakeSelections(INFILE[i], AUXFILE[i], t.u[i]);
      if(t.u[0] == 1)
	Make_111_CWS(AUXFILE, eq);
      else 
	Make_222_CWS(AUXFILE, eq);
    }
    else if(((t.u[0] == 2) && (t.u[1] == 1) && (t.u[2] == 1)) ||
	    ((t.u[0] == 2) && (t.u[1] == 2) && (t.u[2] == 1))){
      if((D[0] + D[1] + D[2] -(t.u[1] + t.u[2]) - d) != 0)
	Die("wrong DIM -cDIM or wrong TYPES -t TYPE1 TYPE2 TYPE3");
      for(i = 0; i < nF; i++)
	MakeSelections(INFILE[i], AUXFILE[i], t.u[i]);
      if(t.u[1] == 2)
	Make_221_CWS(AUXFILE, !strcmp(infile[1], infile[2]));
      else
	Make_211_CWS(AUXFILE, !strcmp(infile[0], infile[1]));
    }
    else if((t.u[0] == t.u[1]) && (t.u[2] == 0)){
      if((D[0] + D[1] + D[2] - t.u[0] - d) != 0)
	Die("wrong DIM -cDIM or wrong TYPES -t TYPE1 TYPE2 TYPE3");
      for(i = 0; i < nF; i++)
	MakeSelections(INFILE[i], AUXFILE[i], t.u[i]);
      Make_nno_CWS(AUXFILE, t.u[0], !strcmp(infile[0], infile[1]));
    }
    else PrintCWSTypes();
  }
  for(i = 0; i < nF; i++){fclose(INFILE[i]);fclose(AUXFILE[i]);}
}

/*  ==========  	      POLY DATA:                	==========  */

void FileRW(char *file, char *m, FILE *rwFILE){

  if((rwFILE = fopen(file, m)) == NULL){
    printf("\n\nUnable to open file %s for %s!\n\n",file,m);
    exit(0);
  }
}

void IP_Poly_Data(int narg, char* fn[])
{
  int r = 1, i, n = 0, p=0, d=0;
  CWS CW;
  PolyPointList *_P, *_DP;
  VertexNumList *_V;
  EqList *_E;

  inFILE=stdin; outFILE=stdout; /*puts("IP_Poly_Data: to be done");*/

  _P = (PolyPointList *) malloc(sizeof(PolyPointList));
  if (_P == NULL) Die("Unable to allocate space for _P");
  _DP = (PolyPointList *) malloc(sizeof(PolyPointList));
  if (_DP == NULL) Die("Unable to allocate space for _DP");
  _V = (VertexNumList *) malloc(sizeof(VertexNumList));
  if (_V == NULL) Die("Unable to alloc space for VertexNumList _V");
  _E = (EqList *) malloc(sizeof(EqList));
  if (_E == NULL) Die("Unable to alloc space for EqList _E");
  
  


  while ((narg > ++n) && (fn[n][0] == '-')){
    if (fn[n][1] == 'i'){
      if (fn[n][2] == 'p') p=1;
      if (fn[n][2] == 'd') d=1;
    }
    if ((fn[n][1] == 'f') || (fn[n][1] == 0))
      inFILE=NULL;
  }
  n--;
  if (narg > ++n){
    if((inFILE = fopen(fn[n], "r")) == NULL){
      printf("\nUnable to open file %s for read\n",fn[n]);
      exit(0);
    }
  }
  if (narg > ++n){
    if((outFILE = fopen(fn[n], "w")) == NULL){
      printf("\nUnable to open file %s for write\n",fn[n]);
      exit(0);
    }
  }
  while (Read_CWS_PP(&CW, _P))
    if (IP_Check(_P,_V,_E)){
      r=1; i=-1;
      while(r && (++i < _E->ne))
	if(_E->e[i].c != 1)
	  r = 0;
      Make_Dual_Poly(_P, _V, _E, _DP);
      if((!p) && (!d)){
	Print_CWS(&CW);
	fprintf(outFILE," M:%d %d",_P->np, _V->nv);
	if(r) fprintf(outFILE," N:%d %d",_DP->np, _E->ne);
	else  fprintf(outFILE," F:%d N:%d", _E->ne,_DP->np);
      }
      if(p) Print_PPL(_P,"");
      if(d) Print_PPL(_DP,"");	
      assert(IP_Check(_DP,_V,_E));
      fprintf(outFILE,"\n");
    }
}
/*  ==========  	      END POLY DATA                	==========  */

/*  ==========  	         Convex Hull:                	==========  */

int Remove_Identical_Points(PolyPointList *);

int ConvHull(PolyPointList *P1, PolyPointList *P2, PolyPointList *P, 
	     VertexNumList *V1, int x)
{
  int i, j;
  VertexNumList *V2;
  EqList *E1, *E2;
  
  E1 = (EqList *) malloc(sizeof(EqList));
  if (E1 == NULL) Die("Unable to alloc space for EqList E1");
  E2 = (EqList *) malloc(sizeof(EqList));
  if (E2 == NULL) Die("Unable to alloc space for EqList E2");
  V2 = (VertexNumList *) malloc(sizeof(VertexNumList));
  if (V2 == NULL) Die("Unable to alloc space for VertexNumList V2");

  Find_Equations(P1, V1, E1); Find_Equations(P2, V2, E2);
  
  if((V1->nv+V2->nv) > VERT_Nmax) Die("increase VERT_Nmax!");
  if((V1->nv+V2->nv) > POINT_Nmax) Die("increase POINT_Nmax!");
  if((P1->n+x) > POLY_Dmax) Die("increase POLY_Dmax!");
  if(x<0) Die("if input is -p#, # must be less than dim of first poly");

  P->np=V1->nv;
  for(i=0; i<V1->nv; i++){
    for(j=P1->n;j<(P2->n+x);j++)P->x[i][j]=0;
    for(j=0;j<P1->n;j++) P->x[i][j]=P1->x[V1->v[i]][j];
  }
  for(i=0; i<V2->nv; i++){
    for(j=0;j<x;j++)P->x[P->np][j]=0;
    for(j=0;j<P2->n;j++) P->x[P->np][x+j]=P2->x[V2->v[i]][j];
    P->np++;
  }
  P->n=P2->n+x;
  if(x==0) Remove_Identical_Points(P);
  i=Find_Equations(P, V1, E1);
  Sort_VL(V1);
  free(E1);free(E2);free(V2);
  return i;
}
void Conv(int narg, char* fn[])
{
  FILE *INFILE[2];
  int n=0, x=0, nF=2, i;
  char *infile[2] = {NULL}, *outfile = NULL, *a;
  PolyPointList *P[2], *PP;
  CWS *CW[2];
  VertexNumList *V;

  V = (VertexNumList *) malloc(sizeof(VertexNumList));
  if (V == NULL) Die("Unable to alloc space for VertexNumList V");
  PP = (PolyPointList *) malloc(sizeof(PolyPointList));
  if (PP == NULL) Die("Unable to allocate space for PolyPointList");

  for(i = 0; i < nF; i++){
    P[i] = (PolyPointList *) malloc(sizeof(PolyPointList));
    if (P[i] == NULL) Die("Unable to allocate space for PolyPointList");
    CW[i] = (CWS *) malloc(sizeof(CWS));
    if (CW[i] == NULL) Die("Unable to allocate space for CWS");
  }

  while ((narg > ++n) && (fn[n][0] == '-')){
    if(fn[n][1] == 'p'){
      if(fn[n][2]!=0){ a=&fn[n][2];
      if(!IsDigit(*a)) Die("after -c there must be a digit!");
      x=atoi(a);}
    }
  }
  for(i = 0; i < nF; i++){
    if((n >= narg)||(fn[n][0] == '-'))
	{printf("#infiles = %d < %d!\n",i,nF); exit(0);}
    infile[i] = fn[n];
    n++;
  }
  if(narg > n)
    outfile = fn[n];
  for(i = 0; i < nF; i++)
    if((INFILE[i] = fopen(infile[i], "r"))==NULL)
      Die("Unable to open infile to read");
  if(outfile == NULL) outFILE = stdout;
  else 
    if((outFILE = fopen(outfile, "w")) == NULL){
      printf("\nUnable to open file %s for write\n",fn[n]); exit(0);}
  while(READ_CWS_PP(CW[0], P[0], INFILE[0])){
    while(READ_CWS_PP(CW[1], P[1], INFILE[1]))
      if(ConvHull(P[0], P[1], PP, V, (P[0]->n-x))) Print_VL(PP, V, "Vertices of P");
    rewind(INFILE[1]);
  }
  for(i = 0; i < nF; i++){free(P[i]); free(CW[i]);} free(PP);free(V);  
}
/*  ==========  	      END of Convex Hull		==========  */

/*uses latte instead of Aux_Complete_Poly for counting points*/
void td_Print_EL(EqList *_E, int *n, int suppress_c, const char *comment){
  int i,j;
  char command[100];
  sprintf(command,"rm zzL.tmp");
  system(command);
  outFILE=fopen("zzL.tmp","w");
  fprintf(outFILE,"%d %d  %s\n",_E->ne,(*n)+1,comment);
  for(i=0;i<_E->ne;i++) {
    if (!suppress_c) fprintf(outFILE,"%d",(int) _E->e[i].c);
    for(j=0;j<*n;j++) fprintf(outFILE," %3d",(int) _E->e[i].a[j]); 
    fprintf(outFILE,"\n");}
    fclose(outFILE);
}
Long NP_use_lat(EqList *_E, PolyPointList *_P)
{
    int tmp;
    char command[100];
    sprintf(command,"count zzL.tmp | grep '*' | awk '{print $7}' > zzL.tmp1");
    
    td_Print_EL(_E,&_P->n,0,"");
    system(command);outFILE=fopen("zzL.tmp1","r");
    while((fscanf(outFILE,"%d",&tmp))!=EOF);
    fclose(outFILE);
    return tmp;
}
Long L_Point_Count(Weight *W,PolyPointList *P,VertexNumList *V,EqList *E){
  int i,j,d; 
  Long *G[W_Nmax], GM[W_Nmax][VERT_Nmax];
  d=W->N; for(i=0;i<d;i++) G[i]=GM[i]; W_to_GLZ(W->w,&d, G);
  P->n=d-1; P->np=d; for(i=0;i<d-1;i++)for(j=0;j<d;j++)P->x[j][i]=GM[i+1][j]; 
  Find_Equations(P,V,E); return NP_use_lat(E,P);
}

int Read_Weight(Weight *);
Long Poly_Point_Count(PolyPointList *P,VertexNumList *V,EqList *E);

Long W_Point_Count(Weight *W,PolyPointList *P,VertexNumList *V,EqList *E){
  int i,j,d; 
  Long *G[W_Nmax], GM[W_Nmax][VERT_Nmax];
  d=W->N; for(i=0;i<d;i++) G[i]=GM[i]; W_to_GLZ(W->w,&d, G);
  P->n=d-1; P->np=d; for(i=0;i<d-1;i++)for(j=0;j<d;j++)P->x[j][i]=GM[i+1][j]; 
  return Poly_Point_Count(P,V,E); 
  /*	char c[50]="#points="; sprintf(&c[8],"%d",P->np);  
	if(P->np<20) Print_PPL(P,c); else printf("%s\n",c); */
}
void SimplexPointCount(int narg, char* fn[])
{ Weight W; VertexNumList V; EqList *E = (EqList *) malloc(sizeof(EqList));
  int L; PolyPointList *P = (PolyPointList *) malloc(sizeof(PolyPointList));
  if ((E == NULL)||(P == NULL)) Die("Unable to allocate space for E or P");
  assert(narg>1); if(fn[1][2]=='f') inFILE=NULL; W.M=0; L=(fn[1][1]=='L');
  while(Read_Weight(&W)){int n; Long np= L ? L_Point_Count(&W,P,&V,E) :
    W_Point_Count(&W,P,&V,E); 
    if(np<=SIMPLEX_POINT_Nmax){printf("%d",W.d); 
      for(n=0;n<W.N;n++) printf(" %ld",W.w[n]);
      printf(" N:%ld\n",np);fflush(0);}}
  free(P); free(E); }
