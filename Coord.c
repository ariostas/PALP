#include "Global.h"
#include "Rat.h"

#undef	TEST_Wbase
#undef  USE_Old_Wbase
#define NO_COORD_IMPROVEMENT		/* switch off weight permutation */

typedef struct {Long x[AMBI_Dmax][AMBI_Dmax]; int n, N;}     CWLatticeBasis;

void Make_CWS_Points(CWS *_C, PolyPointList *_P);
void Make_RGC_Points(CWS *Cin, PolyPointList *_P);
void CWS_to_PermCWS(CWS *Cin, CWS *C, int *pi);

/*  ==========  	  I/O functions:                	==========  */

int  IsNextDigit(void){
  char c; c=fgetc(inFILE); ungetc(c,inFILE);
  if(c=='0') return -1;
  if((c<'0') || ('9'<c)) return 0; else return 1;
}

void Print_PPL(PolyPointList *_P, const char *comment){
  int i,j;
  if(_P->np>20){
    fprintf(outFILE,"%d %d  %s\n",_P->np,_P->n,comment);
    for(i=0;i<_P->np;i++) {
      for(j=0;j<_P->n;j++) fprintf(outFILE,"%d ",(int) _P->x[i][j]); 
      fprintf(outFILE,"\n");}}
  else {
    fprintf(outFILE,"%d %d  %s\n",_P->n,_P->np,comment);
    for(i=0;i<_P->n;i++) {
      for(j=0;j<_P->np;j++) fprintf(outFILE," %4d",(int) _P->x[j][i]); 
      fprintf(outFILE,"\n");}}
}

void Print_VL(PolyPointList *_P, VertexNumList *_V, const char *comment){
  int i,j;
  if(_V->nv>20){
    fprintf(outFILE,"%d %d  %s\n",_V->nv,_P->n,comment);
    for(i=0;i<_V->nv;i++) {
      for(j=0;j<_P->n;j++) fprintf(outFILE,"%d ",(int) _P->x[_V->v[i]][j]); 
      fprintf(outFILE,"\n");}}
  else {
    fprintf(outFILE,"%d %d  %s\n",_P->n,_V->nv,comment);
    for(i=0;i<_P->n;i++) {
      for(j=0;j<_V->nv;j++) fprintf(outFILE," %4d",(int) _P->x[_V->v[j]][i]); 
      fprintf(outFILE,"\n");}}
}

void Print_EL(EqList *_E, int *n, int suppress_c, const char *comment){
  int i,j;
  fprintf(outFILE,"%d %d  %s\n",_E->ne,*n,comment);
  for(i=0;i<_E->ne;i++) {
    for(j=0;j<*n;j++) fprintf(outFILE," %3d",(int) _E->e[i].a[j]); 
    if (!suppress_c) fprintf(outFILE," %5d",(int) _E->e[i].c);
    fprintf(outFILE,"\n");}
}

void Print_Matrix(Long Matrix[][VERT_Nmax], int n_lines, int n_columns, 
		  const char *comment){
  int i,j;
  fprintf(outFILE,"%d %d  %s\n",n_lines, n_columns, comment);
  for(i=0;i<n_lines;i++) {
    for(j=0;j<n_columns;j++) fprintf(outFILE," %3d",(int) Matrix[i][j]); 
    fprintf(outFILE,"\n");}
}
int  auxString2Int(char *c,int *n)
{    int j=0; *n=0;  while(('0'<=c[j])&&(c[j]<='9')) *n = 10*(*n)+c[j++]-'0';
     if(j) while(c[j]==' ') j++;
     return j;
}
void CWSZerror(char *c)
{    printf("Format error %s in Read_CWS_Zinfo\n",c);exit(0);
}
void Print_CWS_Zinfo(CWS *CW)
{    int i,j; if(CW->nw) for(i=0;i<CW->nz;i++) 
     {	fprintf(outFILE,"/Z%d: ",CW->m[i]);
	for(j=0;j<CW->N;j++) fprintf(outFILE,"%d ",CW->z[i][j]);
     }	
}
int  Read_CWS_Zinfo(FILE *inFILE,CWS *CW)		      /* return !EOF */
{    int *nz=&CW->nz; int i=0,n; char c[999],b=' '; *nz=0; 
     for(n=0;n<999;n++)
     {	c[n]=fgetc(inFILE); if(feof(inFILE)) return 0; if(c[n]=='\n') break;
     }  if(n==999) {puts("Out of space in Read_CWS_Zinfo");exit(0);}
     while(c[i]==b)i++;
     if((c[i]=='=')&&(c[i+1]=='d'))i+=2;
     while(c[i]==b)i++; 
     /*	printf("Read_CWS_Zinfo n=%d -> ",n);
	{int I; for(I=0;I<n;I++)printf("%c",c[I]);puts("");} */
     while(i<n)
     {	int j, k, s; if((c[i]!='/')||(c[i+1]!='Z')) return 1;
        i+=2; assert(*nz<POLY_Dmax);
        if((j=auxString2Int(&c[i],&CW->m[*nz])))
	{   if(c[i+j]!=':') CWSZerror(":"); else {i+=j+1; while(c[i]==b)i++;}
	    for(k=0;k<CW->N;k++)
	    {	if((j=auxString2Int(&c[i],&CW->z[*nz][k]))) 
			{if((j)) i+=j; else CWSZerror("missing");} /* ????? */
	    }   s=0;
	    for(k=0;k<CW->N;k++) s+=CW->z[*nz][k]; 
	    /*	if(s % CW->m[*nz])  CWSZerror("det!=1"); */
	    (CW->nz)++;
        }   else CWSZerror("Order");
     }	return 1;
}

void checkDimension(int polyDim, int codim, int index){
  if (index == 1){ 
    if (POLY_Dmax  < (polyDim + codim - 1)){
      if (codim > 1){
	printf("Please increase POLY_Dmax to at least %d = %d + %d - 1\n",
	       (polyDim + codim - 1), polyDim, codim);
	printf("(POLY_Dmax >= dim N + codim - 1 is required)\n");}
      else printf("Please increase POLY_Dmax to at least %d\n", polyDim);
      exit(0);  }}
  else if (POLY_Dmax  < polyDim +1){
    printf("Please increase POLY_Dmax to at least %d = %d + 1\n",
	   (polyDim + 1), polyDim);
    printf("(option -G requires POLY_Dmax >= dim(cone) = dim(support) + 1)\n"
	   );
    exit(0);}
}
 
int  ReadCwsPp(CWS *_CW, PolyPointList *_P, int codim, int index)
/*   _P is always an M-lattice polytope
     codim = 1, index = 1: CY hypersurface 
     codim > 1, index = 1: _P reflexive, CICY with codimension codim
     codim = 1, index > 1: _P a Gorenstein polytope with index index, 
                           determines a reflexive Gorenstein cone       */
{    int i, j, FilterFlag=(inFILE==NULL); 
     int IN[AMBI_Dmax*(AMBI_Dmax+1)], S; 
     static int InputOK;
     _CW->nw=_CW->N=_CW->nz=0;
     _CW->index = index;
     if(FilterFlag) inFILE=stdin; 
     else if(inFILE==stdin)     {
       puts("Degrees and weights  `d1 w11 w12 ... d2 w21 w22 ...'");
       puts("  or `#lines #columns' (= `PolyDim #Points' or `#Points PolyDim'):");}
     for(i=0;i<AMBI_Dmax*(AMBI_Dmax+1);i++)
     { char c;
       while(' '==(c=fgetc(inFILE )));
       ungetc(c,inFILE);    /* read blanks */
       if(IsNextDigit()) fscanf(inFILE,"%d",&IN[i]); else break;
     }
     if(i==0) { if(!InputOK){puts("-h gives you help\n"); exit(0);}
	else return 0; } InputOK++;
     if(i==1) { puts("Error in INPUT: need at least 2 numbers!"); exit(0);}
     if(i==2) 					      /* READ PolyPointList */
     {	int tr=0;
       while('\n'!=fgetc(inFILE));      /* read to end of line */
       if(IN[0]==IN[1]){ 
	 puts("The number of points must be larger than the dimension!");
	 exit(0);}
       if(IN[0]>IN[1]) {tr=IN[0];IN[0]=IN[1];IN[1]=tr;} tr=!tr;
	checkDimension(IN[0], codim, index);
 	if(IN[1]>POINT_Nmax) {puts("Please increase POINT_Nmax"); exit(0);}
 	if((inFILE==stdin)&&!FilterFlag)  
	  printf("Type the %d coordinates as %s=%d lines with %s=%d columns:\n",
		IN[0]*IN[1], tr ? "dim" : "#pts", tr ? IN[0] : IN[1],
		tr ? "#pts" : "dim", tr ? IN[1] : IN[0]);
	_P->n=IN[0]; _P->np=IN[1];
	/* allow all numbers in one string (distributed over lines) or
	 * matrix blocks with trailing comments	*/
	if(tr) for(i=0;i<IN[0];i++) for(j=0;j<IN[1];j++)
	  {   int X; fscanf(inFILE,"%d",&X); _P->x[j][i]=X;	}
	else for(i=0;i<IN[1];i++) for(j=0;j<IN[0];j++)
	  {   int X; fscanf(inFILE,"%d",&X); _P->x[i][j]=X;	}
	/* Finish_Poly_Points(_P); */
	while(fgetc(inFILE )-'\n') if(feof(inFILE)) return 0;/* read to EOL */
     	if(FilterFlag) inFILE=NULL;
	return 1;
     }  /* End of reading PolyPointList */
     assert(i!=3);
     S=IN[i-1]; 
     for(j=0;j<i-1;j++) if((IN[j]==0)||(S<IN[j])) break; 
     if(j==i-1)				 /* Single Weights with: "w1 ... d" */
     {	_CW->nw=1; _CW->d[0]=S; _CW->N=i-1; 
        assert(_CW->N <= POLY_Dmax+1);		      /* Increase POLY_Dmax */ 
	for(j=0;j<_CW->N;j++) {_CW->W[0][j]=IN[j]; assert(IN[j]>0);} 
	goto MAP;
     }
     _CW->d[0] = IN[0];                        /* ASSIGN (COMBINED) WEIGHTS */
     S = IN[0] * index; 		     
     for(j=1; j<i; j++)
     {	if(S<IN[j]) 
	{   if(S) { puts("Error in INPUT: degree vs. weights!"); exit(0);}
	    else  
	    {	_CW->d[++(_CW->nw)]=IN[j];
	        S = IN[j] * index; 
		if(1==_CW->nw) 
		{   if(i%(1+_CW->N)) {puts("INPUT error: numbers?"); exit(0);}
		}   if(j%(1+_CW->N)) {puts("INPUT error: degrees?"); exit(0);}
	    }
	}
	else
	{   _CW->W[_CW->nw][j - 1 - (1 + _CW->N)*_CW->nw] = IN[j]; S-=IN[j];
	    if(0==_CW->nw) if(++(_CW->N)>AMBI_Dmax) 
	    {   puts("Increase AMBI_Dmax!"); exit(0);
	    }
	}
     }
     ++(_CW->nw);
     checkDimension(_CW->N - _CW->nw, codim, index);

MAP: for(i=0;i<_CW->nw;i++)			/* check consistency of CWS */
     {	Long sum = _CW->d[i] * index, *w = _CW->W[i];
	for(j=0;j<_CW->N;j++) { assert(w[j]>=0); sum-=w[j]; }
	if((sum)&&(index==1)){ 
	  printf("Use option -l for (single) WeightSystems with ");
	  printf("d!=\\sum(w)\n(only Read_Weight makes the correct ");
	  puts("PolyPointList in that case)"); exit(0);}
     }
     _CW->nz=0;
     if(!Read_CWS_Zinfo(inFILE,_CW))			   /* read Z to EOL */
     {	if(!InputOK) puts("-h gives you help\n"); return 0;}
     Make_CWS_Points(_CW, _P);
     if(FilterFlag) inFILE=NULL;
     return 1;
}

int  Read_CWS_PP(CWS *_CW, PolyPointList *_P){
  return ReadCwsPp(_CW, _P, 1, 1);
}

int  Read_PP(PolyPointList *_P)
{    int i, j, FilterFlag=(inFILE==NULL);
     int IN[AMBI_Dmax*(AMBI_Dmax+1)];
     static int InputOK;
     /* _CW->nw=_CW->N=_CW->nz=0; */

     if(FilterFlag) inFILE=stdin;
     else if(inFILE==stdin)     {
       printf("`#lines #columns' (= `PolyDim #Points' or `#Points PolyDim'):\n");
     };
     for(i=0;i<AMBI_Dmax*(AMBI_Dmax+1);i++)
     { char c;
       while(' '==(c=fgetc(inFILE )));
       ungetc(c,inFILE);    /* read blanks */
       if(IsNextDigit()) fscanf(inFILE,"%d",&IN[i]); else break;
     }
     if(i==0) { if(!InputOK){puts("-h gives you help\n"); exit(0);}
	else return 0; } InputOK++;
     if(i==1) { puts("Error in INPUT: need at least 2 numbers!"); exit(0);}
     if(i==2) 					      /* READ PolyPointList */
     {	int tr=0;
       while('\n'!=fgetc(inFILE));      /* read to end of line */
       if(IN[0]==IN[1]){ 
	 puts("The number of points must be larger than the dimension!");
	 exit(0);}
        if(IN[0]>IN[1]) {tr=IN[0];IN[0]=IN[1];IN[1]=tr;} tr=!tr;
	if(IN[0]>POLY_Dmax) {puts("increase POLY_Dmax!"); exit(0);}
 	if(IN[1]>POINT_Nmax) {puts("increase POINT_Nmax!"); exit(0);}
 	if((inFILE==stdin)&&!FilterFlag)
	  printf("Type the %d coordinates as %s=%d lines with %s=%d columns:\n",
		IN[0]*IN[1], tr ? "dim" : "#pts", tr ? IN[0] : IN[1],
		tr ? "#pts" : "dim", tr ? IN[1] : IN[0]);
	_P->n=IN[0]; _P->np=IN[1];
/* allow all numbers in one string (distributed over lines) or
 * matrix blocks with trailing comments	*/
	if(tr) for(i=0;i<IN[0];i++) for(j=0;j<IN[1];j++)
	  {   int X; fscanf(inFILE,"%d",&X); _P->x[j][i]=X;	}
	else for(i=0;i<IN[1];i++) for(j=0;j<IN[0];j++)
	  {   int X; fscanf(inFILE,"%d",&X); _P->x[i][j]=X;	}
	/* Finish_Poly_Points(_P); */
	while(fgetc(inFILE )-'\n') if(feof(inFILE)) return 0;/* read to EOL */
     	if(FilterFlag) inFILE=NULL;
	return 1;
     }
     if(i>2){puts("Error: expected input format is matrix of polytope points!");exit(0);}
     if(FilterFlag) inFILE=NULL;
     return 1;
}

int  Read_CWS(CWS *_CW, PolyPointList *_P)
{    int i, j, FilterFlag=(inFILE==NULL);
     int IN[AMBI_Dmax*(AMBI_Dmax+1)], S;
     static int InputOK;
     _CW->nw=_CW->N=_CW->nz=0;
     _CW->index = 1;

     if(FilterFlag) inFILE=stdin;
     else if(inFILE==stdin)     {
       printf("Degrees and weights  `d1 w11 w12 ... d2 w21 w22 ...':\n");
     };
     for(i=0;i<AMBI_Dmax*(AMBI_Dmax+1);i++)
     { char c;
       while(' '==(c=fgetc(inFILE )));
       ungetc(c,inFILE);    /* read blanks */
       if(IsNextDigit()) fscanf(inFILE,"%d",&IN[i]); else break;
     }
     if(i==0) { if(!InputOK){puts("-h gives you help\n"); exit(0);}
	else return 0; } InputOK++;
     if(i==1) { puts("Error in INPUT: need at least 2 numbers!"); exit(0);}
     if(i==2) {puts("Error: expected input format is CWS!"); exit(0);}
     assert(i!=3);
     S=IN[i-1]; for(j=0;j<i-1;j++) if(S<IN[j]) break;
     if(j==i-1)				 /* Single Weights with: "w1 ... d" */
     {	_CW->nw=1; _CW->d[0]=S; _CW->N=i-1;
        assert(_CW->N <= POLY_Dmax+1);		      /* Increase POLY_Dmax */
	for(j=0;j<_CW->N;j++) {_CW->W[0][j]=IN[j]; assert(IN[j]>0);}
	goto MAP;
     }
     S=_CW->d[_CW->nw]=IN[0]; 			          /* ASSIGN WEIGHTS */
     for(j=1; j<i; j++)
     {	if(S<IN[j])
	{   if(S) { puts("Error in INPUT: degree vs. weights!"); exit(0);}
	    else
	    {	S=_CW->d[++(_CW->nw)]=IN[j];
		if(1==_CW->nw)
		{   if(i%(1+_CW->N)) {puts("INPUT error: numbers?"); exit(0);}
		}   if(j%(1+_CW->N)) {puts("INPUT error: degrees?"); exit(0);}
	    }
	}
	else
	{   _CW->W[_CW->nw][j - 1 - (1 + _CW->N)*_CW->nw] = IN[j]; S-=IN[j];
	    if(0==_CW->nw) if(++(_CW->N)>AMBI_Dmax)
	    {   puts("Increase AMBI_Dmax!"); exit(0);
	    }
	}
     }
     ++(_CW->nw);
     if(_CW->N-_CW->nw > POLY_Dmax){
       printf("Please increase POLY_Dmax to at least %d\n", _CW->N-_CW->nw);
       exit(0);} /* increase POLY_Dmax */

MAP: for(i=0;i<_CW->nw;i++)			/* check consistency of CWS */
     {	Long sum=_CW->d[i], *w=_CW->W[i];
	for(j=0;j<_CW->N;j++) { assert(w[j]>=0); sum-=w[j]; }
	if(sum){ /*printf("Use poly.x -w for (single) WeightSystems with ");
		   printf("d!=\\sum(w)\n(only Read_Weight makes the correct ");
		   puts("PolyPointList in that case)"); exit(0);*/
	  printf("cannot handle (single) WeightSystems with ");
	  printf("d!=\\sum(w)\n");
	  exit(0);}}
     _CW->nz=0;
     if(!Read_CWS_Zinfo(inFILE,_CW))			   /* read Z to EOL */
     {	if(!InputOK) puts("-h gives you help\n"); return 0;}
     Make_CWS_Points(_CW, _P);				   /* now make POLY */
     if(FilterFlag) inFILE=NULL;
     return 1;
}

void Print_CWH(CWS *_W, BaHo *_BH){
  int i, j;
  for(i=0;i<_W->nw;i++)     {	
    fprintf(outFILE,"%d ",(int) _W->d[i]);
    for(j=0;j<_W->N;j++) fprintf(outFILE,"%d ",(int) _W->W[i][j]);
    if(i+1<_W->nw) fprintf(outFILE," ");     }
  Print_CWS_Zinfo(_W);
  if(_BH->np)     {
    fprintf(outFILE,"M:%d %d N:%d %d",_BH->mp,_BH->mv,_BH->np,_BH->nv);
    if (_BH->n == 3) fprintf(outFILE," Pic:%d Cor:%d",_BH->pic, _BH->cor); 
    if (_BH->n > 3) {
      fprintf(outFILE," H:%d",_BH->h1[1]); 
      for(i=2;i<_BH->n-1;i++) fprintf(outFILE,",%d",_BH->h1[i]);
      if(_BH->n==4) fprintf(outFILE," [%d]",2*(_BH->h1[1]-_BH->h1[2])); 
      if(_BH->n==5) 
	fprintf(outFILE," [%d]",48+6*(_BH->h1[1]-_BH->h1[2]+_BH->h1[3]));}}
  else if(_BH->mp)
    fprintf(outFILE,"M:%d %d F:%d",_BH->mp,_BH->mv,_BH->nv);
  else fprintf(outFILE,"V:%d F:%d",_BH->mv,_BH->nv);
  fprintf(outFILE,"\n");
}

/*  ==========  	      END of I/O functions		==========  */


/*  ==========	   Solve W-Eq. -> triangular LatticeBasis
 *   
 *   B[0]:= (-n1, n0, 0, ...) / GCD(n0,n1);
 *   B[i]:= (0,...,0,g/G,0,...)- (ni/G) * ExtGCD.K(n0,...,n(i-1),0,...);
 *   	    with g=GCD(n0,...,n[i-1]); G=GCD(g,ni);
 *
 *   CWS by iteration: NextWeight=W[]*B; NewB=Basis(NextW); B[i+1]=B[i]*NewB;
 *									    */

void PrintBasis(CWLatticeBasis *_B)
{    int i,j; 		     puts("Basis:");
     for(i=0;i<_B->n;i++)
     {	for(j=0;j<_B->N;j++) fprintf(outFILE,"%6d ",(int) _B->x[i][j]); 
	puts("");
     }	puts("End of Basis  - -");
}

void Orig_Solve_Next_WEq(Long *NW, CWLatticeBasis *_B)
{    int i, j, P=0, p[AMBI_Dmax]; Long W[AMBI_Dmax], G;	_B->n=_B->N-1; 
     for(i=0;i<_B->N;i++) 
     {	for(j=0;j<_B->n;j++) _B->x[j][i]=0; 		      /* init B.x=0 */
	if(NW[i]) {p[P]=i; W[P++]=NW[i];}		/* non-zero weights */
     }
     if(P<2) puts("need two non-zero weights in  >>Solve_Next_WEq<<");
     for(i=0;i<p[0];i++) _B->x[i][i]=1;
     while((++i)<p[1]) _B->x[i-1][i]=1;
     G=Fgcd(W[0],W[1]); if(W[0]/G<0) G=-G;
     _B->x[i-1][p[0]]=-W[1]/G; _B->x[i-1][p[1]]=W[0]/G;
     j=2; while(++i<_B->N)
     {	if(NW[i])
	{   int k; Long *X=_B->x[i-1], K[AMBI_Dmax], g=REgcd(W, &j, K); 
	    G=Fgcd(g,NW[i]); if(g/G<0) G=-G; X[i]= g/G; g=W[j]/G;
	    for(k=0;k<j;k++) X[p[k]]=-K[k]*g;
	    j++;
	}
	else _B->x[i-1][i]=1;
     }
}
void Solve_Next_WEq(Long *NW, CWLatticeBasis *_B)
{    Long W[AMBI_Dmax], *X[AMBI_Dmax], GLZ[AMBI_Dmax][AMBI_Dmax];
     int i, j, P=0, p[AMBI_Dmax]; _B->n=_B->N-1;
#ifdef TEST_Wbase
     Orig_Solve_Next_WEq(NW,_B); PrintBasis(_B);
#endif
     for(i=0;i<_B->N;i++) 
     {	for(j=0;j<_B->n;j++) _B->x[j][i]=0;		      /* init B.x=0 */
	if(NW[i]) {p[P]=i; X[P]=GLZ[P]; W[P++]=NW[i];}	/* non-zero weights */
     }	if(P>1) W_to_GLZ(W,&P,X);		/* P>1, compute GLZ */
     else {/* printf("P=%d W[0]=%d for W_to_GLZ\n",P,W[0]);exit(0);*/assert(P);
	for(i=0;i<p[0];i++)_B->x[i][i]=1;
	while((++i)<_B->N)_B->x[i-1][i]=1;
	return; }
     for(i=1;i<P;i++) if(X[i][i]<0) for(j=0;j<=i;j++) X[i][j] *= -1;
     for(i=0;i<p[0];i++) _B->x[i][i]=1;
     while((++i)<p[1]) _B->x[i-1][i]=1;
     _B->x[i-1][p[0]]=X[1][0]; _B->x[i-1][p[1]]=X[1][1]; j=2;
     while(++i<_B->N)
     {	if(NW[i])
	{   int k; Long *B=_B->x[i-1]; for(k=0;k<=j;k++) B[p[k]]=X[j][k]; j++;
	}
	else _B->x[i-1][i]=1;
     }
#ifdef TEST_Wbase
	/* for(i=0;i<_B->N;i++)printf("%d ",NW[i]);puts("=NW");for(i=0;i<P;i++)
	{for(j=0;j<P;j++)printf("%4d ",GLZ[i][j]);puts("=GL");} */
	printf("New version: "); PrintBasis(_B);
#endif
#ifdef USE_Old_Wbase
	Orig_Solve_Next_WEq(NW,_B);
#endif
}
void Make_CWS_Basis(CWS *_C, CWLatticeBasis *_B){    
  int i,j,k,l; Long W[AMBI_Dmax]; CWLatticeBasis NB, AuxB; 
  AuxB.N=_B->N=_C->N;
  Solve_Next_WEq(_C->W[0],_B);
  for(i=1;i<_C->nw;i++){ 
    for(j=0;j<_B->n;j++){   
      W[j]=0; 
      for(k=0;k<_B->N;k++) W[j]+=_C->W[i][k]*_B->x[j][k]; }
    NB.N=_B->n; 
    Solve_Next_WEq(W,&NB);
    AuxB.n=NB.n;
    for(j=0;j<AuxB.N;j++) for(k=0;k<AuxB.n;k++) {   
      AuxB.x[k][j]=0; 
      for(l=0;l<NB.N;l++) AuxB.x[k][j] += NB.x[k][l] * _B->x[l][j]; }
    *_B=AuxB;     }
  /* init coordinate hyperplanes */
  _C->B.ne=_B->N;	  /* == transposed of Basis */
  for(i=0;i<_B->N;i++) {
    for(j=0;j<_B->n;j++) _C->B.e[i].a[j]=_B->x[j][i];
    _C->B.e[i].c=1;}
}

void Poly_To_Ambi(CWLatticeBasis *_B, Long *x, Long *X)
{    int i, j;
     for(i=0;i<_B->N;i++)
     {	X[i]=1; for(j=0;j<_B->n;j++) X[i] += _B->x[j][i] * x[j];
     }
}
/*   B.x[i][A] is lower triangular: zero for i<p and Amin[p] <= A <Amin[p+1]
 *   X^A=x^p B_p^A+1, hence  for A \in [ Amin[p],Amin[p+1] ) and all k:
 *   x^iB_i^A \in [0,d_k/w_k^A] - 1 - \sum_{l<i} x^l B_l^A
 *									     */

void QuotZ_2_SublatG(Long Z[][VERT_Nmax],int *zm,Long *M,int *d,
     Long G[][POLY_Dmax]); 		      /* normalize and diagonalize Z */
void CWS_2_SublatZ(CWS *C,CWLatticeBasis *B,			       /* in */
	int *m, Long *M, Long G[POLY_Dmax][POLY_Dmax])		      /* out */
{    int i,j,k,d=C->N-C->nw; Long Z[POLY_Dmax][VERT_Nmax]; *m=C->nz; 
     for(i=0;i<C->nz;i++){ M[i]=C->m[i]; for(j=0;j<d;j++)
     {	Z[i][j]=0; for(k=0;k<C->N;k++) Z[i][j]+=C->z[i][k]*B->x[j][k];}
     }	/* PrintBasis(B); */	QuotZ_2_SublatG(Z,m,M,&d,G);
}
void Reduce_PPL_2_Sublat(PolyPointList *P,int *nm,Long *M,Long G[][POLY_Dmax])
{    int i,j,n=0,N; for(N=0;N<P->np;N++) 
     {	Long X[POLY_Dmax]; for(i=0;i<P->n;i++)
	{   X[i]=0; for(j=0;j<P->n;j++) X[i]+=G[i][j]*P->x[N][j];
	}   for(i=0;i<*nm;i++) if(X[i]%M[i]) break; if(i<*nm) continue;
	for(i=0;i<*nm;i++) P->x[n][i]=X[i]/M[i]; 
	for(i=*nm;i<P->n;i++) P->x[n][i]=X[i];
	n++;
     }	P->np=n; /* Print_PPL(P,"");printf("nm=%d\n",*nm); */
}
Long PD_Floor(Long N,Long D)	/*  assuming PosDenom  D>0:  F <= N/D < F+1  */
{    Long F=N/D; return (F*D>N) ? F-1 : F; 
}

void Old_Make_CWS_Points(CWS *Cin, PolyPointList *_P)
{    int i, j, Amin[POLY_Dmax+1]; Long *x=_P->x[_P->np=0], xmin[POLY_Dmax], 
	xmax[POLY_Dmax], Xmax[AMBI_Dmax], xaux[POLY_Dmax], L, R; CWS *_C=Cin;
     CWLatticeBasis B; Long G[POLY_Dmax][POLY_Dmax],M[POLY_Dmax];int m=Cin->nz;
#ifndef NO_COORD_IMPROVEMENT		/* ==== Perm Coord Improvement ==== */
     int pi[AMBI_Dmax]; CWS Caux; _C=&Caux; CWS_to_PermCWS(Cin,_C, pi);
#endif				       /* = End of Perm Coord Improvement = */
     Make_CWS_Basis(_C, &B);	  		 /* make `triangular' Basis */
#ifndef NO_COORD_IMPROVEMENT		/* ==== Perm Coord Improvement ==== */
     assert(_C->N==_C->B.ne);
     for(i=0;i<Cin->N;i++) Cin->B.e[i]=_C->B.e[pi[i]]; Cin->B.ne=_C->N;
#endif				       /* = End of Perm Coord Improvement = */
     i=_P->n=B.n; Amin[0]=0; Amin[B.n]=j=B.N;	      /* inversion structure */
     while(--i) {while(!B.x[i-1][--j]); Amin[i]=++j;}
     for(i=0;i<B.N;i++) 				     /* compute Xmax */
     {	Xmax[i]=0; for(j=0;j<_C->nw;j++) if(_C->W[j][i])
	{   L=_C->d[j]/_C->W[j][i]; 
	    if(Xmax[i]) {if(L<Xmax[i]) Xmax[i]=L;}
	    else Xmax[i]=L;
	}
     }
     j=B.n-1; i=Amin[j+1]-1; R=B.x[j][i];     /* compute xmin[j] and xmax[j] */
     xmin[j]=-PD_Floor(1,R); xmax[j]=PD_Floor(Xmax[i]-1,R);   /* since R > 0 */
/** /printf("R=%2d:  %2d <= x[%d=&%d] <= %2d\n",R,xmin[j],j,i,xmax[j]);/ **/
     while((i--)>Amin[j]) 	/* if(R=B.x[B.n-1][i]):  R!=0  => new limits */
     {	Long Low=-1, Upp=Low+Xmax[i]; 		  R=B.x[B.n-1][i];
	if(R>0) {if(xmax[j]>(L=PD_Floor(Upp,R)))   xmax[j]=L;
		 if(xmin[j]<(L=-PD_Floor(-Low,R))) xmin[j]=L;}
	else	{if(xmax[j]>(L=PD_Floor(-Low,-R))) xmax[j]=L;
		 if(xmin[j]<(L=-PD_Floor(Upp,-R))) xmin[j]=L;}
/** /printf("R=%2d:  %2d <= x[%d=&%d] <= %2d\n",R,xmin[j],j,i,xmax[j]);/ **/
     }				   /* this completes the limits for x[B.n-1] */
     x[j]=xmin[j]; while(j<B.n)
     {	int k; if(x[j]>xmax[j]) {if(B.n==(++j)) break; else x[j]++;}
	else		/* compute limits[j-1] and initialize x[j-1] */
	{   Long Low=-1, Upp=Xmax[i=Amin[j--]-1]; int RangeFlag=0;
	    for(k=j+1;k<B.n;k++) Low-=x[k]*B.x[k][i];	   /* compute offset */
	    Upp+=Low; R=B.x[j][i]; 
	    xmin[j]=-PD_Floor(-Low,R); xmax[j]=PD_Floor(Upp,R);
/** /printf("R=%2d:  %2d <= x[%d=&%d] <= %2d\n",R,xmin[j],j,i,xmax[j]);/ **/
	    while((i--)>Amin[j]) if((R=B.x[j][i]))    /* R!=0  => new limits */
     	    {	Low=-1; Upp=Xmax[i];
		for(k=j+1;k<B.n;k++) Low-=x[k]*B.x[k][i];
		Upp+=Low; R=B.x[j][i]; 
		if(R>0) {if(xmax[j]>(L=PD_Floor(Upp,R)))   xmax[j]=L;
		 	 if(xmin[j]<(L=-PD_Floor(-Low,R))) xmin[j]=L;}
		else	{if(xmax[j]>(L=PD_Floor(-Low,-R))) xmax[j]=L;
			 if(xmin[j]<(L=-PD_Floor(Upp,-R))) xmin[j]=L;}
/** /printf("R=%2d:  %2d <= x[%d=&%d] <= %2d\n",R,xmin[j],j,i,xmax[j]);/ **/
	    }			     	/* completes limits for x[] except */
	    else			/*   when R=0 and X out of Range:  */
	    {	Long X=1; for(k=j+1;k<B.n;k++) X+=x[k]*B.x[k][i];
		if((X<0)||(X>Xmax[i])) RangeFlag=1;
	    }
	    if(RangeFlag) ++x[++j]; else x[j]=xmin[j];            
	    if(j==0)
	    {   while(x[0]<=xmax[0]) 
	    	{   Long *y; 
		    if((++_P->np) < POINT_Nmax) 
		    {	y=(_P->x[_P->np]); for(k=0;k<B.n;k++) y[k]=x[k]; 
		    }
		    else if(_P->np == POINT_Nmax)
		    {	y=xaux; for(k=0;k<B.n;k++) y[k]=x[k];
		    }
		    else {puts("Increase POINT_Nmax");exit(0);}
		    x=y; ++x[0];
	    	}
	        x[j=1]++;
	    }
	}
     }	if(m)CWS_2_SublatZ(_C,&B,&m,M,G);if(m)Reduce_PPL_2_Sublat(_P,&m,M,G);
}

int Compute_X0(int N, CWS *_C, Long *X0){
  /* finds a first X = X0 such that sum w_i X_i = d for every weight w in _C,
     this becomes the origin after the change of coordinates from X to x */
  int j;
  Long Xmax=0; 
  if (!N){
    for(j=0;j<_C->nw;j++) 
      if(_C->W[j][0]){
	if (_C->d[j]%_C->W[j][N]) return 0;
	if(Xmax) {if (Xmax != _C->d[j]/_C->W[j][N]) return 0;}
	else Xmax = _C->d[j]/_C->W[j][N];   }
      else if(_C->d[j]) return 0;
    X0[0] = Xmax;
    return 1;}
  for(j=0;j<_C->nw;j++) if(_C->W[j][N]!=0){
      Long L=_C->d[j]/_C->W[j][N]; 
      if(Xmax) {if(L<Xmax) Xmax = L;}
      else Xmax = L;   }
  for (X0[N]=0; X0[N] <= Xmax; X0[N]++){
    if (Compute_X0(N-1, _C, X0)){
      for(j=0; j<_C->nw; j++) _C->d[j] += X0[N] * _C->W[j][N]; /* reset Cin*/
      return 1;}
    for(j=0; j<_C->nw; j++) _C->d[j] -= _C->W[j][N];}
  for(j=0; j<_C->nw; j++) _C->d[j] += (Xmax + 1) * _C->W[j][N]; /* reset Cin */
  return 0;
}

void Make_CWS_Points(CWS *Cin, PolyPointList *_P)
{    int i, j, Amin[POLY_Dmax+1], m=Cin->nz; 
     Long *x=_P->x[_P->np=0], xmin[POLY_Dmax], 
       xmax[POLY_Dmax], Xmax[AMBI_Dmax], X0[AMBI_Dmax], xaux[POLY_Dmax], L, R; 
     CWS *_C=Cin;
     CWLatticeBasis B; Long G[POLY_Dmax][POLY_Dmax],M[POLY_Dmax];
#ifndef NO_COORD_IMPROVEMENT		/* ==== Perm Coord Improvement ==== */
     int pi[AMBI_Dmax]; CWS Caux; _C=&Caux; CWS_to_PermCWS(Cin,_C, pi);
#endif				       /* = End of Perm Coord Improvement = */
     Make_CWS_Basis(_C, &B);	  		 /* make `triangular' Basis */
#ifndef NO_COORD_IMPROVEMENT		/* ==== Perm Coord Improvement ==== */
     assert(_C->N==_C->B.ne);
     for(i=0;i<Cin->N;i++) Cin->B.e[i]=_C->B.e[pi[i]]; Cin->B.ne=_C->N;
#endif				       /* = End of Perm Coord Improvement = */
     if (Cin->index == 1) for (i=0; i< Cin->N; i++) X0[i] = 1;
     else if(!Compute_X0(Cin->N - 1, Cin, X0)){_P->n=0; puts("no X0!");return;}
     /* X0 is the reference point in X-space that is transformed to the 
	origin of x-space for _P    */
     /*printf("\nX0: ");for (i=0;i<Cin->N;i++) printf("%ld ", X0[i]);puts("");
     for(i=0;i<Cin->nw;i++)     {	
       fprintf(outFILE,"%d ",(int) Cin->d[i]);
       for(j=0;j<Cin->N;j++) fprintf(outFILE,"%d ",(int) Cin->W[i][j]);
       if(i+1<Cin->nw) fprintf(outFILE," ");     }*/
     i=_P->n=B.n; Amin[0]=0; Amin[B.n]=j=B.N;	      /* inversion structure */
     while(--i) {while(!B.x[i-1][--j]); Amin[i]=++j;}
     for(i=0;i<B.N;i++) 				     /* compute Xmax */
     {	Xmax[i]=0; for(j=0;j<_C->nw;j++) if(_C->W[j][i])
	{   L=_C->d[j]/_C->W[j][i]; 
	    if(Xmax[i]) {if(L<Xmax[i]) Xmax[i]=L;}
	    else Xmax[i]=L;
	}
     }
     j=B.n-1; i=Amin[j+1]-1; R=B.x[j][i];     /* compute xmin[j] and xmax[j] */
     xmin[j] = -PD_Floor(X0[i],R); 
     xmax[j] = PD_Floor(Xmax[i]-X0[i],R);     /* since R > 0 */
/** /printf("R=%2d:  %2d <= x[%d=&%d] <= %2d\n",R,xmin[j],j,i,xmax[j]);/ **/
     while((i--)>Amin[j]) 	/* if(R=B.x[B.n-1][i]):  R!=0  => new limits */
     {	Long Low=-X0[i], Upp=Low+Xmax[i]; 		  R=B.x[B.n-1][i];
	if(R>0) {if(xmax[j]>(L=PD_Floor(Upp,R)))   xmax[j]=L;
		 if(xmin[j]<(L=-PD_Floor(-Low,R))) xmin[j]=L;}
	else	{if(xmax[j]>(L=PD_Floor(-Low,-R))) xmax[j]=L;
		 if(xmin[j]<(L=-PD_Floor(Upp,-R))) xmin[j]=L;}
/** /printf("R=%2d:  %2d <= x[%d=&%d] <= %2d\n",R,xmin[j],j,i,xmax[j]);/ **/
     }				   /* this completes the limits for x[B.n-1] */
     x[j]=xmin[j]; 
     while(j<B.n)
     {	int k; if(x[j]>xmax[j]) {if(B.n==(++j)) break; else x[j]++;}
	else		/* compute limits[j-1] and initialize x[j-1] */
	{   Long Upp=Xmax[i=Amin[j--]-1], Low=-X0[i]; int RangeFlag=0;
	    for(k=j+1;k<B.n;k++) Low-=x[k]*B.x[k][i];	   /* compute offset */
	    Upp+=Low; R=B.x[j][i]; 
	    xmin[j]=-PD_Floor(-Low,R); xmax[j]=PD_Floor(Upp,R);
/** /printf("R=%2d:  %2d <= x[%d=&%d] <= %2d\n",R,xmin[j],j,i,xmax[j]);/ **/
	    while((i--)>Amin[j]) if((R=B.x[j][i]))    /* R!=0  => new limits */
     	    {	Low=-X0[i]; Upp=Xmax[i];
		for(k=j+1;k<B.n;k++) Low-=x[k]*B.x[k][i];
		Upp+=Low; R=B.x[j][i]; 
		if(R>0) {if(xmax[j]>(L=PD_Floor(Upp,R)))   xmax[j]=L;
		 	 if(xmin[j]<(L=-PD_Floor(-Low,R))) xmin[j]=L;}
		else	{if(xmax[j]>(L=PD_Floor(-Low,-R))) xmax[j]=L;
			 if(xmin[j]<(L=-PD_Floor(Upp,-R))) xmin[j]=L;}
/** /printf("R=%2d:  %2d <= x[%d=&%d] <= %2d\n",R,xmin[j],j,i,xmax[j]);/ **/
	    }			     	/* completes limits for x[] except */
	    else			/*   when R=0 and X out of Range:  */
	    {	Long X=1; for(k=j+1;k<B.n;k++) X+=x[k]*B.x[k][i];
		if((X<0)||(X>Xmax[i])) RangeFlag=1;
	    }
	    if(RangeFlag) ++x[++j]; else x[j]=xmin[j];            
	    if(j==0)
	    {   while(x[0]<=xmax[0]) 
	    	{   Long *y; 
		    if((++_P->np) < POINT_Nmax) 
		    {	y=(_P->x[_P->np]); for(k=0;k<B.n;k++) y[k]=x[k]; 
		    }
		    else if(_P->np == POINT_Nmax)
		    {	y=xaux; for(k=0;k<B.n;k++) y[k]=x[k];
		    }
		    else {puts("Increase POINT_Nmax");exit(0);}
		    x=y; ++x[0];
	    	}
	        x[j=1]++;
	    }
	}
     }	if(m)CWS_2_SublatZ(_C,&B,&m,M,G);if(m)Reduce_PPL_2_Sublat(_P,&m,M,G);
}

/*  ==========    Coordinate improvement via CWS-Permutations   ==========  */

typedef struct {Long s, g, *N, **G; int *p;}		Tri_GLZ_MPaux;  
Long Tri_GLZ_Norm(int *d, Long **S)		      /* sum-norm = default */
{    Long norm=0; int i,j; for(i=0;i<*d;i++)for(j=0;j<*d;j++)
#ifdef	MAX_NORM						/* max(|*|) */
     {	Long x=S[i][j]; if(x<0) x=-x; if(x>norm) norm=x; } return norm;
#else								/* sum(|*|) */
     {	Long x=S[i][j]; norm += (x>0) ? x : -x ; } return norm; 
#endif
}
void Tri_GLZ_Basis_Perm(int *d,int *pi, /* int *pinv, */ Tri_GLZ_MPaux * AP)
{    Long g, norm, N[AMBI_Dmax], M[AMBI_Dmax][AMBI_Dmax], *S[AMBI_Dmax];  
     int i, j; /* if(pi[0]>pi[1]) return; *//* eliminate equiv. permut. */
     for(i=0;i<*d;i++) { N[i]=AP->N[pi[i]]; S[i]=M[i];}
     g=W_to_GLZ(N,d,S); norm=Tri_GLZ_Norm(d,S);
     if(AP->s) assert(g==AP->g); else AP->g=g;
#ifdef	TEST_MIN_GLZ
     {	Long err=0, max=0, tsum=0, pos; for(i=0;i<*d;i++)
	{   for(j=0;j<*d;j++) {pos=S[i][j]; if(pos<0)pos=-pos; tsum+=pos; 
	    if(max<pos)max=pos;} if(err==0) { 
		for(j=0;j<*d;j++) err+=N[j]*S[i][j]; if(i==0)err-=g;}
	}   assert((norm==tsum)||(norm==max)); 
	printf("max=%lld, sum=%lld\n",(long long) max, (long long) tsum);
	if(err) { for(i=0;i<*d;i++) printf("%d",pi[i]); 
	printf("  g=%d  norm=%d",(int)g,(int)norm); puts(""); 
	for(i=0;i<*d;i++){printf("S[%d]=",i);for(j=0;j<*d;j++)printf(" %5d",
	S[i][j]);if(!i)printf("   g=%d  norm=%d",(int)g,(int)norm);puts("");}}
     }
#endif
     if((0 == AP->s) || (norm < AP->s))		   /* init or improve AP->G */
     {	for(i=0;i<*d;i++)for(j=0;j<*d;j++)AP->G[i][pi[j]]=S[i][j]; AP->s=norm;
	for(i=0;i<*d;i++) AP->p[i]=pi[i];
     }
}

#ifndef NO_COORD_IMPROVEMENT	
				/* improved by permutation pi=pinv^-1 */
Long Wperm_to_GLZ(Long *W, int *d, Long **G, int *P)
{    Tri_GLZ_MPaux AS; int i, j, pi[AMBI_Dmax], pinv[AMBI_Dmax];
     AS.s=0; AS.N=W; AS.G=G; AS.p=P;
     Map_Permut(d,pi, /* pinv,*/ (Tri_GLZ_Basis_Perm),(void *) &AS);
     for(i=0;i<*d;i++)for(j=0;j<*d;j++)G[i][j]=AS.G[i][j];
#ifdef	TEST_MIN_GLZ
     	for(i=0;i<*d;i++){printf("  G[%d]=",i);for(j=0;j<*d;j++)printf(" %5d",
	G[i][j]);if(!i)printf("  g=%d norm=%d",(int)AS.g,(int)AS.s);puts("");}
#endif
     return AS.g;
}
void CWS_to_PermCWS(CWS *Cin, CWS *C, int *pi)
{    int i, j, N=C->N=Cin->N, n=C->nw=Cin->nw, l=0, A[AMBI_Dmax]; 
     Long *X, *G[AMBI_Dmax], M[AMBI_Dmax][AMBI_Dmax], W[AMBI_Dmax]; 
     for(j=1;j<n;j++) if(Cin->d[j]>Cin->d[l]) l=j; X=Cin->W[l]; 
     j=0; for(i=0;i<N;i++) if(X[i]) {A[j]=i; W[j]=labs(X[i]); G[j]=M[j]; j++;}
     Wperm_to_GLZ(W,&j,G,pi); for(i=0;i<j;i++) pi[i]=A[pi[i]];
     for(i=0;i<N;i++) if(X[i]==0) pi[j++]=i;
     C->d[0]=Cin->d[l]; for(i=0;i<N;i++) C->W[0][i]=Cin->W[l][pi[i]];
     for(j=1;j<n;j++) {int L=j-(j<=l); C->d[j]=Cin->d[L];
	for(i=0;i<N;i++) C->W[j][i]=Cin->W[L][pi[i]];}
     C->nz=Cin->nz; for(j=0;j<C->nz;j++){C->m[j]=Cin->m[j];
	for(i=0;i<N;i++) C->z[j][i]=Cin->z[j][pi[i]];}
}
#endif

/*  ==========        For WS for 5d-polytopes, Sep 2017         ==========  */
int int_ld(Long w){int i=-1; while (w) {w /= 2; i++;} return i;}

void Initialize_C5S(C5stats *_C5S, int n){
  int k;
  if (n < 5) {puts("Option '-Q' requires POLY_Dmax > 4!"); exit(0);};
  _C5S->n_nonIP = 0;
  _C5S->n_IP_nonRef = 0;
  _C5S->n_ref = 0;
  _C5S->nr_max_mp = 0;
  _C5S->nr_max_mv =  0;
  _C5S->nr_max_nv = 0;
  _C5S->nr_max_w = 0;
  for (k=0; k<MAXLD; k++) {_C5S->n_w[k] = 0; _C5S->nr_n_w[k] = 0;}
  _C5S->max_mp = 0;
  _C5S->max_mv = 0;
  _C5S->max_np = 0;
  _C5S->max_nv = 0;
  _C5S->max_h22 = 0;
  _C5S->max_w = 0;
  for (k=1; k<n-1; k++) _C5S->max_h1[k] = 0;
  for (k=0; k<=n; k++) _C5S->max_nf[k] = 0;
  _C5S->max_chi = -100000000;
  _C5S->min_chi = 100000000; 
}

void Update_C5S(BaHo *_BH, int *nf, Long *W, C5stats *_C5S){
  assert(POLY_Dmax>4);
  if (_BH->np) { // reflexive case
    int i, chi = 48+6*(_BH->h1[1]-_BH->h1[2]+_BH->h1[3]), ld = int_ld(W[5]);
    assert(0<=ld); assert(ld<MAXLD);
    if (_BH->mp > _C5S->max_mp) _C5S->max_mp = _BH->mp;
    if (_BH->mv > _C5S->max_mv) _C5S->max_mv = _BH->mv;
    if (_BH->np > _C5S->max_np) _C5S->max_np = _BH->np;
    if (_BH->nv > _C5S->max_nv) _C5S->max_nv = _BH->nv;
    if (_BH->h22 > _C5S->max_h22) _C5S->max_h22 = _BH->h22;
    for (i=1; i<_BH->n-1; i++)
      if (_BH->h1[i] > _C5S->max_h1[i]) _C5S->max_h1[i] = _BH->h1[i];
    for (i=0; i<=_BH->n; i++)
      if (nf[i] > _C5S->max_nf[i]) _C5S->max_nf[i] = nf[i];
    if (chi > _C5S->max_chi) _C5S->max_chi = chi;
    if (chi < _C5S->min_chi) _C5S->min_chi = chi;
    _C5S->n_ref++;
    _C5S->n_w[ld]++;
    if(W[5] > _C5S->max_w) _C5S->max_w = W[5];}
  else { // IP but non-reflexive case
    if (_BH->mp > _C5S->nr_max_mp) _C5S->nr_max_mp = _BH->mp;
    if (_BH->mv > _C5S->nr_max_mv)  _C5S->nr_max_mv = _BH->mv;
    if (_BH->nv > _C5S->nr_max_nv) _C5S->nr_max_nv = _BH->nv;
    _C5S->n_IP_nonRef++;
    _C5S->nr_n_w[int_ld(W[5])]++;
    if(W[5] > _C5S->nr_max_w) _C5S->nr_max_w = W[5];}
}
    
void Print_C5S(C5stats *_C5S){
  int i;
  assert(POLY_Dmax>4);
  printf("non-IP: #=%ld\n", _C5S->n_nonIP);
  printf("IP, non-reflexive: #=%ld, max_mp=%d, max_mv=%d, max_nv=%d, max_w=%ld\n",
	 _C5S->n_IP_nonRef, _C5S->nr_max_mp, _C5S->nr_max_mv, _C5S->nr_max_nv,
	 _C5S->nr_max_w);
  printf("  #(w5) of given ld: ");
  for (i=0; i<MAXLD; i++) printf(" %d:%ld", i, _C5S->nr_n_w[i]);
  puts("");
  printf("reflexive: #=%ld, max_mp=%d, max_mv=%d, max_np=%d, max_nv=%d, max_w=%ld\n",
	 _C5S->n_ref, _C5S->max_mp, _C5S->max_mv, _C5S->max_np, _C5S->max_nv,
	 _C5S->max_w);
  printf("  #(w5) of given ld: ");
  for (i=0; i<MAXLD; i++) printf(" %d:%ld", i, _C5S->n_w[i]);
  puts("");
  printf("  max #(faces): %d %d %d %d %d\n", _C5S->max_nf[0], _C5S->max_nf[1],
	 _C5S->max_nf[2], _C5S->max_nf[3], _C5S->max_nf[4]);
  printf("  h11<=%d, h12<=%d, h13<=%d, h22<=%d, %d<=chi<=%d\n",
	 _C5S->max_h1[1], _C5S->max_h1[2], _C5S->max_h1[3],
	 _C5S->max_h22, _C5S->min_chi, _C5S->max_chi);
}

