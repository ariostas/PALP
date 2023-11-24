/*  ======================================================================  */
/*  ==========                                                  ==========  */
/*  ==========                P O L Y . C                       ==========  */
/*  ==========                                                  ==========  */
/*  ======================================================================  */

/*
 *	Coord   get coordinates by reading them or converting weight input
 *	Rat	rational functions
 *	Vertex	computes Vertices and Faces
 *      Polynf	normal form and symmetries
 */

#include "Global.h"
#include "LG.h"
#define OSL (42)  /* opt_string's length */


FILE *inFILE, *outFILE;

void  PrintUsage(char *c){
  int i;
  char *opt_string[OSL]={
  "h  print this information           ",
  "f  use as filter                    ",
  "g  general output:                  ",
  "   P reflexive: numbers of (dual) points/vertices, Hodge numbers ",
  "   P not reflexive: numbers of points, vertices, equations    ",
  "p  points of P                      ",
  "v  vertices of P                    ",
  "e  equations of P/vertices of P-dual",
  "m  pairing matrix between vertices and equations                  ",
  "d  points of P-dual (only if P reflexive)          ",
  "a  all of the above except h,f      ",
  "l  LG-`Hodge numbers' from single weight input                   ",
  "r  ignore non-reflexive input       ",
  "o  ignore non-IP input       ",
  "q  quick version of og, incompatible with other options",   
  "Q  like 'q', with statistics on weights for 5d classification",   
  "D  dual polytope as input (ref only)",
  "n  do not complete polytope or calculate Hodge numbers        ",
  "i  incidence information            ",
  "s  check for span property (only if P from CWS)           ",
  "I  check for IP property            ",
  "S  number of symmetries             ",
  "T  upper triangular form	       ",
  "N  normal form                      ",
  "t  traced normal form computation   ",
  "V  IP simplices among vertices of P*",
  "P  IP simplices among points of P* (with 1<=codim<=# when # is set)",
  "Z  lattice quotients for IP simplices",
  "#  #=1,2,3  fibers spanned by IP simplices with codim<=#        ",
  "## ##=11,22,33,(12,23): all (fibered) fibers with specified codim(s) ",
  "   when combined: ### = (##)#       ",
  "A  affine normal form",
  "B  Barycenter and lattice volume [# ... points at deg #]",
  "F  print all facets",
  "G  Gorenstein: divisible by I>1",
  "L  like 'l' with Hodge data for twisted sectors",
  "U  simplicial facets in N-lattice",
  "U1 Fano (simplicial and unimodular facets in N-lattice)",
  "U5 5d fano from reflexive 4d projections (M lattice)",
  "C1 conifold CY (unimodular or square 2-faces)",
  "C2 conifold FANO (divisible by 2 & basic 2 faces)",
  "E  symmetries related to Einstein-Kaehler Metrics"};

  puts("");
  printf("This is '%s':  computing data of a polytope P\n",c);
  printf("Usage:   %s [-<Option-string>] [in-file [out-file]]\n", c);
  puts("");
  printf("Options (concatenate any number of them into <Option-string>):\n");
  for (i=0;i<OSL;i++) puts(opt_string[i]);
  puts("");
  puts("Input:    degrees and weights `d1 w11 w12 ... d2 w21 w22 ...'");
  puts("          or `d np' or `np d' (d=Dimension, np=#[points]) and");
  puts("              (after newline) np*d coordinates");
  puts("Output:   as specified by options");
}

typedef	struct 	{int p[SYM_Nmax][VERT_Nmax];}		VPermList;

int main (int narg, char* fn[]){
  int n=0, k, FilterFlag=0, lg=0, s=0, i=0, m=0, p=0, v=0, e=0, d=0, t=0, z=0,
    S=0, N=0, I=0, r=0, nc=0, g=0, D=0, IP, R, Tr, T=0, PS=0, VS=0, CD=0, ZS=1,
    A=0, B=0, G=0, F=0, U=0, dd=0, Einstein=0, o=0, q=0, Q=0;
  char c; 
  CWS *CW=(CWS *) malloc(sizeof(CWS));
  Weight W;
  VertexNumList V;
  EqList *E = (EqList *) malloc(sizeof(EqList));
  EqList *DE = (EqList *) malloc(sizeof(EqList));
  BaHo BH;
  VaHo VH; 
  PolyPointList *_P = (PolyPointList *) malloc(sizeof(PolyPointList)),
               *_DP = (PolyPointList *) malloc(sizeof(PolyPointList));
  FaceInfo *FI=NULL;
  PairMat *PM = (PairMat *) malloc(sizeof(PairMat)),
         *DPM = (PairMat *) malloc(sizeof(PairMat));
  C5stats C5S;

  if((CW==NULL)||(E==NULL)||(_P==NULL)||(DE==NULL)||(_DP==NULL)
	||(PM==NULL)||(DPM==NULL)) {
    puts("Allocation failure: Reduce dimensions!"); exit(0);}
  CW->nw=0;

  while(narg > ++n) {
    if(fn[n][0]!='-') break;
    k=0;
    while ((c=fn[n][++k])!='\0'){
      if(c=='A') A=1;
      else if(c=='B') B=1;
      else if(c=='F') F=1;
      else if(c=='G') G=1;
      else if(c=='L') lg=2;
      else if(c=='U') U=1;
      else if(c=='C') U=2;
      else if(c=='E') Einstein=1;
      else if(c=='h') { PrintUsage(fn[0]); exit(0);}
      else if(c=='f') FilterFlag=1;
      else if(c=='g') g=1;
      else if(c=='l') lg=1; 
      else if(c=='s') s=1;
      else if(c=='i') i=1;
      else if(c=='I') I=1;
      else if(c=='m') m=1;
      else if(c=='p') p=1;
      else if(c=='v') v=1;
      else if(c=='e') e=1;
      else if(c=='d') d=1;
      else if(c=='t') t=1;
      else if(c=='S') S=1;
      else if(c=='N') N=1;
      else if(c=='T') T=1;
      else if(c=='r') r=1;
      else if(c=='o') o=1;
      else if(c=='q') q=1;
      else if(c=='Q') Q=1;
      else if(c=='n') nc=1;
      else if(c=='P') PS=1;
      else if(c=='V') VS=1;
      else if(c=='Z') ZS=-1;
      else if(c=='z') z=1;
      else if(c=='D') D=1;
      else if(('0'<=c)&&(c<='9')) CD=10*CD+c-'0';
      else if(c=='a') {g=1; m=1; p=1; v=1; e=1; d=1; }
      else {printf("Unknown option '-%c'; use -h for help\n",c); exit(0);}}}
  n--;

  if(s+i+I+m+p+v+e+d+t+S+N+T+PS+VS+CD+G+A+B+F+z==0) g=1;
  VH.sts=(lg==2);
  if ((U==1)&&(CD==1)&&(s+i+I+m+p+v+e+d+t+S+N+T+PS+VS+G+A+B+F+z==0)) g=1;
  if(g+lg+p+d+PS+CD==0) nc=1; /* don't need completion of points */
  if((T==1)&&(B+U+lg+g+s+i+I+m+p+v+e+d+t+S+N+PS+VS+(1-ZS)==0)){ puts(
    "\n-T: Please specify desired output, e.g. via -v or -p \n");exit(0);}
  if(FilterFlag) {inFILE=NULL; outFILE=stdout;}
  else {
    if (narg > ++n)  inFILE=fopen(fn[n],"r");
    else inFILE=stdin;
    if (inFILE==NULL){printf("Input file %s not found!\n",fn[n]);exit(0);}
    if (narg > ++n) outFILE=fopen(fn[n],"w");
    else outFILE=stdout;     }	
  if(U) {dd=CD; CD=0; if((U==2)||(dd==5)) nc=0;}
  if(i){
    FI=(FaceInfo *) malloc(sizeof(FaceInfo));
    if (FI==NULL) {puts("Unable to allocate space for FaceInfo FI"); exit(0);}}
  if(Q) Initialize_C5S(&C5S, POLY_Dmax); // Initialize statistics
  if(Einstein) Einstein_Metric(CW,_P,&V,E);
  while(lg ? Read_W_PP(&W,_P) : Read_CWS_PP(CW,_P)) {
    if(q||Q) {
      FaceInfo FI;
      if(!QuickAnalysis(_P, &BH, &FI)) {if(Q) C5S.n_nonIP++; continue;} //non-IP
      Print_CWH(CW, &BH);
      if(Q) Update_C5S(&BH, FI.nf, CW->W[0], &C5S);
      continue;}
    if(T) {
      if (CW->nw) {
	puts("Please do not use weight input with the -T option"); continue;}
      else Make_Poly_UTriang(_P);}
    R=0;
    if ((IP=Find_Equations(_P,&V,E))){
      if (D){
	int k;
	*_DP=*_P;
	R=EL_to_PPL(E, _P, &_P->n);
	VNL_to_DEL(_DP, &V, E);
	for (k=0;k<_P->np;k++) V.v[k]=k;
	V.nv=_P->np;      }
      else R=EL_to_PPL(E, _DP, &_P->n);}
    else if (o&&!IP) continue;
    if (D&&!R) {fprintf(outFILE,"Input not reflexive!\n"); continue;}
    if (r&&!R) continue;
    Sort_VL(&V);
    Make_VEPM(_P,&V,E,*PM);
    if(!nc) {
      if(D||!(lg||CW->nw)) Complete_Poly(*PM,E,V.nv,_P);
      if(R&&!(D&&(lg||CW->nw))) {
	if(0==Transpose_PM(*PM, *DPM, V.nv, E->ne))
	{   fprintf(stderr,"Transpose_PM failed because #eq=%d > VERT_Nmax\n",
	    E->ne);exit(0);}
	VNL_to_DEL(_P, &V, DE);
	Complete_Poly(*DPM,DE,E->ne,_DP);}}
    if (U==1) {
      if(dd==5) {if(!Fano5d(_P,&V,E)) continue;}
      else if(!SimpUnimod(_P,&V,E,dd)) continue; }
    if (U==2) {
      VNL_to_DEL(_P,&V,DE); 
      if(!ConifoldSing(_P,&V,E,_DP,DE,dd)) continue;}
    if(g){
      if(!R||nc) {BH.mp=_P->np; BH.mv=V.nv; BH.np=0; BH.nv=E->ne;}
      else RC_Calc_BaHo(_P,&V,E,_DP,&BH);
      if(lg) {
	if ((Tr=Trans_Check(W))) LGO_VaHo(&W,&VH);
	Write_WH(&W, &BH, &VH, R, Tr, _P, &V, E); }
      else Print_CWH(CW, &BH); }
    if(s&&CW->nw) if(!Span_Check(E,&(CW->B),&_P->n))
      fprintf(outFILE,"No Span\n");
    if(I && !IP) fprintf(outFILE,"No IP\n");
    if(p) Print_PPL(_P,"Points of P");
    if(v) Print_VL(_P, &V, "Vertices of P");
    if(e) Print_EL(E, &_P->n, R,
          (R ? "Vertices of P-dual <-> Equations of P" : "Equations of P"));
    if(i){Make_Incidence(_P,&V,E,FI); Print_FaceInfo(_P->n,FI);}
    if(m) Print_Matrix(*PM, E->ne, V.nv,
		       "Pairing matrix of vertices and equations of P");
    if(d&&(_DP->np>E->ne)) Print_PPL(_DP, "Points of P-dual");
    if(S||N||t){
      int SymNum /*, VPMSymNum*/; Long NF[POLY_Dmax][VERT_Nmax]; 
      VPermList *VP = (VPermList*) malloc(sizeof(VPermList)); 
      assert(VP!=NULL);
      /* VPMSymNum=*/ Make_Poly_Sym_NF(_P, &V, E, &SymNum, VP->p, NF, t, S, N);
      free(VP);}
    if(R&&(PS||VS||CD)) IP_Simplices(_DP, (!D)*E->ne, PS*ZS, VS*ZS, CD);
    if(G) {
      char divi[99]; Long g=Divisibility_Index(_P,&V); 
      if(g>1){sprintf(divi,"divisible by factor=%ld",g); Print_VL(_P,&V,divi);}}
    if(B) {
      Long vB[POLY_Dmax],Z; int j;
      Long vol=LatVol_Barycent(_P,&V,vB,&Z); printf("vol=%ld, baricent=(",vol); 
      for(j=0;j<_P->n;j++) printf("%s%ld",j?",":"",vB[j]);
      printf(")/%ld\n",Z); if(CD) IPs_degD(_P,&V,E,CD);}
    if(F) { 
      int j, cc;
      Long VM[POLY_Dmax][VERT_Nmax];
      for (j=0; j<E->ne; j++){
	Make_Facet(_P, &V, E, j, VM, &cc);
	Print_Matrix(VM,_P->n-1,cc,"");}  }
    if(A) {
      Long ANF[POLY_Dmax][VERT_Nmax];
      Make_ANF(_P,&V,E,ANF); 
      Print_Matrix(ANF, _P->n, V.nv,"Affine normal form");}
    fflush(outFILE);     }
  if(Q) Print_C5S(&C5S);
  return 0;
}
