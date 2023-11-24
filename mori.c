/* =========================================================== */
/* ===                                                     === */
/* ===                  m o r i . c                        === */
/* ===                                                     === */
/* ===	Authors: Maximilian Kreuzer, Nils-Ole Walliser	   === */
/* ===	Last update: 19/04/12                              === */
/* ===                                                     === */
/* =========================================================== */	


/* ======================================================== */
/* =========            H E A D E R s             ========= */

#include "Global.h"
#include "LG.h"
#include "Mori.h"

/*==========================================================*/

FILE *inFILE, *outFILE;

void  PrintUsage(char *c){
  printf("This is ``%s'':  star triangulations of a polytope P* in N\n",c);
  printf("                     Mori cone of the corresponding toric ambient spaces\n");
  printf("                     intersection rings of embedded (CY) hypersurfaces\n");
  printf("Usage:   %s [-<Option-string>] [in-file [out-file]]\n", c);
  printf("Options (concatenate any number of them into <Option-string>):\n");

  printf("    -h      print this information \n");
  printf("    -f      use as filter\n");
  printf("    -g      general output: triangulation and Stanley-Reisner ideal\n");
  printf("    -I      incidence information of the facets (ignoring IPs of facets)\n");
  printf("    -m      Mori generators of the ambient space\n");
  printf("    -P      IP-simplices among points of P* (ignoring IPs of facets)\n");
  printf("    -K      points of P* in Kreuzer polynomial form\n");
  printf("    -b      arithmetic genera and Euler number\n");
  printf("    -i      intersection ring\n");
  printf("    -c      Chern classes of the (CY) hypersurface\n");
  printf("    -t      triple intersection numbers\n");
  printf("    -d      topological information on toric divisors & del Pezzo conditions\n");
  printf("    -a      all of the above except h, f, I and K\n");
  printf("    -D      lattice polytope points of P* as input (default CWS)\n");
  printf("    -H      arbitrary (also non-CY) hypersurface `H = c1*D1 + c2*D2 + ...'\n");
  printf("            input: coefficients `c1 c2 ...'\n");
  printf("    -M      manual input of triangulation\n");
  puts("Input: 1) standard: degrees and weights `d1 w11 w12 ... d2 w21 w22 ...'");
  puts("       2) alternative (use -D): `d np' or `np d' (d=Dimension, np=#[points])");
  puts("                                and (after newline) np*d coordinates");
  puts("Output:   as specified by options");
}

int main (int narg, char* fn[]){


  int  n=0, i, k;

  /* flags */
  MORI_Flags Flag;
  
       Flag.FilterFlag = 0; // filter
       Flag.g = 0; // -g: general output
       Flag.m = 0; // -m: Mori cone
       Flag.P = 0; // -P: IP simplices
       Flag.K = 0; // -K: Newton polynomial
       Flag.i = 0; // -i: intersection ring
       Flag.t = 0; // -t: triple intersection number
       Flag.c = 0; // -c: Chern classes
       Flag.d = 0; // -d; del Pezzo
       Flag.a = 0; // -a: all of the above except h,f and K
       Flag.b = 0; // -b: Hodge number of toric div
       Flag.D = 0; // -D: dual poly as input
       Flag.H = 0; // -H: arbitrary hypersurface
       Flag.I = 0; // -I: incidence information
       Flag.M = 0; // -M: allows to insert a triangulation
       Flag.Read_HyperSurfCounter = 0; // see Mori.h for description
  char c;

  CWS *CW=(CWS *) malloc(sizeof(CWS));

  VertexNumList V;
  EqList *E = (EqList *) malloc(sizeof(EqList));
  EqList *DE = (EqList *) malloc(sizeof(EqList));
  
  PolyPointList *_P = (PolyPointList *) malloc(sizeof(PolyPointList)),
               *_DP = (PolyPointList *) malloc(sizeof(PolyPointList));

  PairMat *PM = (PairMat *) malloc(sizeof(PairMat)),
         *DPM = (PairMat *) malloc(sizeof(PairMat));

  if((CW==NULL)||(E==NULL)||(_P==NULL)||(DE==NULL)||(_DP==NULL)||(PM==NULL)||(DPM==NULL)){
    puts("Allocation failure: Reduce dimensions!");
    exit(0);
  }
  CW->nw=0;

  while(narg > ++n) {
    if(fn[n][0]!='-') break;
    k=0;
    while ((c=fn[n][++k])!='\0'){
      if(c=='h') { PrintUsage(fn[0]); exit(0);}
      if(c=='f') Flag.FilterFlag=1;
      if(c=='g') Flag.g=1;
      if(c=='m') Flag.m=1;
      if(c=='P') Flag.P=1;
      if(c=='K') Flag.K=1;
      if(c=='i') Flag.i=1;
      if(c=='t') Flag.t=1;
      if(c=='c') Flag.c=1;
      if(c=='d') Flag.d=1;
      if(c=='a') Flag.a=1;
      if(c=='b') Flag.b=1;
      if(c=='D') Flag.D=1;
      if(c=='H') Flag.H=1;
      if(c=='I') Flag.I=1;
      if(c=='M') Flag.M=1; 
      }
  }
  n--;

  /*if ((Flag.M)&&(!Flag.D)){ 
    puts("-M works only when combined with -D!");
    exit(0);}*/
  if(Flag.g + Flag.m + Flag.P + Flag.K + Flag.i + Flag.t + Flag.c + Flag.d
     + Flag.a + Flag.b + Flag.H + Flag.I ==0) Flag.g=1;

  if(Flag.a){
	  Flag.g=1;
	  Flag.m=1;
	  Flag.P=1;
	  //Flag.K=1;
	  Flag.i=1;
	  Flag.t=1;
	  Flag.c=1;
	  Flag.d=1;
	  Flag.b=1;
  }

  if(Flag.H==1 && 
     (Flag.g + Flag.m + Flag.P + Flag.K + Flag.i + Flag.t + Flag.c + 
           Flag.d + Flag.a + Flag.b + Flag.I ==0 )){
	Flag.b=1;
	//Flag.g=1;	
	}

  if(Flag.FilterFlag) {inFILE=NULL; outFILE=stdout;}

  else {
    if (narg > ++n)  inFILE=fopen(fn[n],"r");
    else inFILE=stdin;

    if (inFILE==NULL){printf("Input file %s not found!\n",fn[n]);exit(0);}

    if (narg > ++n) outFILE=fopen(fn[n],"w");
    else outFILE=stdout;
  }
  
  while((Flag.D ? Read_PP(_P) : Read_CWS(CW,_P))) {
    if (!Ref_Check(_P,&V,E)){
      fprintf(outFILE,"Input not reflexive!\n");
      continue;    }
    if (Flag.D == 0){ /* dualize: _P should become the N-polytope! */
      assert(EL_to_PPL(E, _P, &_P->n));
      assert(Ref_Check(_P,&V,E));}
    Sort_VL(&V);
    if (!(Flag.D&&Flag.M)){
      Make_VEPM(_P,&V,E,*PM);
      Complete_Poly(*PM,E,V.nv,_P);
      for (i=V.nv; i<_P->np-1; i++)
	if(Vec_is_zero(_P->x[i],_P->n)) {
	  Swap_Vecs(_P->x[i],_P->x[_P->np-1],_P->n);
	  break;      }   }
    else {
      for (i=0; i<_P->np; i++)
	if(Vec_is_zero(_P->x[i],_P->n)) {
	  Swap_Vecs(_P->x[i],_P->x[_P->np-1],_P->n);
	  break;      }
      if (i==_P->np){
	for (k=0;k<_P->n;k++) _P->x[_P->np][k]=0;
	_P->np++;   } }
    if(Flag.M){
      if (POLY_Dmax  < (_P->np - _P->n)){
	printf("Please increase POLY_Dmax to at least %d = %d - %d - 1\n",
	       (_P->np - _P->n -1), _P->np, _P->n);
	printf("(%s -M requires POLY_Dmax >= #(points) - dim N -1)\n",
	       fn[0]);
	exit(0);      }   }
    HyperSurfDivisorsQ(_P,&V,E,&Flag);
    fflush(outFILE);  }
  return 0;
}

