#include "Global.h"
#include "Subpoly.h"
#include "Rat.h"

#define  subl_int                LLong
#define  UnAided_IP_CHECK        (POLY_Dmax>4)
#define  SIMPLE_CTH	         (0)
#define	 INCOMPLETE_SL_REDUCTION (0)
#define  TEST_Aided_IP_CHECK     (0)
#ifndef	 CEQ_Nmax	
#define  CEQ_Nmax                EQUA_Nmax
#endif

#define  IMPROVE_SL_COORD	(0)	      /* 0=no  1=SL(old)  2=GL(new) */
#define  IMPROVE_SL_REGCD	(0)

typedef struct {int nk, k[VERT_Nmax];} KeepList;
typedef struct {int ne; Equation e[CEQ_Nmax];}              CEqList;

/*  ==========	Auxiliary routines from other modules           ==========  */

Long CompareEq(Equation *X, Equation *Y, int n);
int  IsGoodCEq(Equation *_E, PolyPointList *_P, VertexNumList *_V);
int  Finish_IP_Check(PolyPointList *_P, VertexNumList *_V, EqList *_F,
                    CEqList *_CEq, INCI *F_I, INCI *CEq_I); /* from Vertex.c */
int  Improve_Coords(PolyPointList *_P,VertexNumList *_V);  /* from Polynf.c */

int  New_Improve_Coords(PolyPointList *_P,VertexNumList *_V)
{    switch(IMPROVE_SL_COORD)
     {	case 0: return 1; case 1: return Improve_Coords(_P,_V);
	case 2: Make_Poly_UTriang(_P); return 1; default: assert(0);}
}

/*  ======================================================================  */
/*  ==========		     			  		==========  */
/*  ==========		S U B P O L Y S        			==========  */
/*  ==========							==========  */ 
/*  ======================================================================  */

int  IsNewEq(Equation *E, CEqList *_C, EqList *_F, int *_n){
  int i=_C->ne; 
  while(i--) if(!CompareEq(E,&_C->e[i],*_n)) return 0;
  i=_F->ne; 
  while(i--) if(!CompareEq(E,&_F->e[i],*_n)) return 0;
  return 1;
}

int Irrel(INCI *_INCI, INCI *INCI_List, int *n_irrel){
  int i;
  for (i=0;i<*n_irrel;i++) if (INCI_EQ(*_INCI, INCI_List[i])) return 1;
  return 0;
}

void INCI_To_VertexNumList(INCI X, VertexNumList *_V, int n){
  _V->nv=0;  
  while(!INCI_EQ_0(X)) {n--; if(INCI_M2(X)) _V->v[_V->nv++]=n; X=INCI_D2(X);} 
}

void Remove_INCI_From_List(INCI New_INCI, INCI *INCI_List, int *n_INCI){
  int i;
  for (i=*n_INCI-1;i>=0;i--) {
    if(INCI_LE(INCI_List[i],New_INCI)) INCI_List[i]=INCI_List[--(*n_INCI)];}
}    

void Add_INCI_To_List(INCI New_INCI, INCI *INCI_List, int *n_INCI){
  int i;
  for (i=*n_INCI-1;i>=0;i--) {
    if(INCI_LE(New_INCI,INCI_List[i])) return;
    if(INCI_LE(INCI_List[i],New_INCI)) INCI_List[i]=INCI_List[--(*n_INCI)];}
  assert(*n_INCI < CD2F_Nmax);
  INCI_List[(*n_INCI)++]=New_INCI;
}
/*
An INCI_List corresponds to a list of facets of some polytope, which might 
be the original polytope or one of its faces. This implies that no two
INCIs x,y of the INCI_List fulfill INCI_LE(x,y). 
Consequently, Add_INCI_To_List enhances INCI_List by New_INCI only if there is
no y in INCI_List such that INCI_LE(New_INCI,y) is fulfilled; 
at the same time it removes any y with INCI_LE(y,New_INCI) from INCI_List.
Remove_INCI_From_List removes any y with INCI_LE(y,New_INCI) from INCI_List.
*/

INCI List_Complete(int n, int n_polys, int n_facets, 
		   INCI *polys, INCI *facets);

INCI Poly_Complete(int n, int *n_polyfacets, INCI *poly, INCI *polyfacets){
  int j, k, n_cd2_faces=0;
  INCI cd2_faces[CD2F_Nmax];
  INCI X=INCI_0();
  for (j=0;j<*n_polyfacets;j++) X=INCI_OR(X,polyfacets[j]);
  if (n==1) return INCI_XOR(X,*poly);
  for (j=0;j<*n_polyfacets;j++) for (k=j+1;k<*n_polyfacets;k++)
    Add_INCI_To_List(INCI_AND(polyfacets[j],polyfacets[k]),
		     cd2_faces,&n_cd2_faces);
  return INCI_OR(INCI_XOR(X,*poly),
	   List_Complete(n-1,*n_polyfacets,n_cd2_faces,polyfacets,cd2_faces));
}
 
INCI List_Complete(int n, int n_polys, int n_facets, 
		   INCI *polys, INCI *facets){
  int i, j, n_polyfacets;
  INCI polyfacets[CEQ_Nmax];
  INCI X=INCI_0();
  for (i=0;i<n_polys;i++){
    n_polyfacets=0;
    for (j=0;j<n_facets;j++) if (INCI_LE(facets[j],polys[i]))
      polyfacets[n_polyfacets++]=facets[j];
    X=INCI_OR(X,Poly_Complete(n, &n_polyfacets, &(polys[i]), polyfacets));}
  return X;
}

void FE_Close_the_Hole(PolyPointList *_P, VertexNumList *_V, EqList *_E, 
           CEqList *_CEq, int n_old_v, INCI *CEq_INCI, INCI Hole_Vert_INCI){
  static PolyPointList P;
  VertexNumList Hole_Verts;
  EqList BVE;
  int i,j,n_new_Eq=0;
  INCI_To_VertexNumList(Hole_Vert_INCI, &Hole_Verts, n_old_v); 
  if(Hole_Verts.nv<_P->n){
    printf("Hole_Verts.nv=%d\n",Hole_Verts.nv);
    puts("Hole_Verts:"); 
    Print_INCI(Hole_Vert_INCI);
    puts("\nCEq_INCI:"); 
    for(i=0;i<_CEq->ne;i++) {Print_INCI(CEq_INCI[i]); puts("");}
    Print_VL(_P,_V,"");
    Print_EL((EqList *) _E,&_P->n,0,"");
    Print_EL((EqList *) _CEq,&_P->n, 0,"");
    exit(0);}
  P.n=_P->n;
  P.np=Hole_Verts.nv;
  for (i=0;i<Hole_Verts.nv;i++) for (j=0;j<_P->n;j++)
    P.x[i][j]=_P->x[Hole_Verts.v[i]][j];
  Find_Equations(&P,&Hole_Verts,&BVE);
  assert(BVE.ne);
  for (i=0;i<BVE.ne;i++) if(IsGoodCEq(&BVE.e[i],_P,_V)) 
    if(IsNewEq(&BVE.e[i],_CEq,_E,&_P->n)){
      n_new_Eq++;
      assert(_CEq->ne<CEQ_Nmax); 
      CEq_INCI[_CEq->ne]=Eq_To_INCI(&BVE.e[i],_P,_V);
      _CEq->e[_CEq->ne++]=BVE.e[i];}
  assert(n_new_Eq);  
}

void Close_the_Hole(PolyPointList *_P, VertexNumList *_V, EqList *_E, 
       CEqList *_CEq, int old_ne, int n_old_v, int n_Hole_Faces, 
       INCI *E_INCI, INCI *CEq_INCI, INCI Hole_Verts, INCI *Hole_Faces){

  INCI Bad_Vert_INCI, cd2_Faces[CD2F_Nmax];
  int n_cd2_Faces=n_Hole_Faces, np=_P->np, i, j;

  /* puts("CTH");
     Print_PPL(_P);
     fprintf(outFILE,"_E:\n");
     for (i=0;i<_E->ne;i++){
     for (j=0;j<_P->n;j++) fprintf(outFILE,"%d ", _E->e[i].a[j]);
     fprintf(outFILE," %d\n", _E->e[i].c);}
     fprintf(outFILE,"n_irrel=%d, n_old_v=%d\n", n_irrel, n_old_v);
     fprintf(outFILE,"E_INCI: ");
     for (i=0;i<_E->ne;i++) Print_INCI(E_INCI[i]); 
     fprintf(outFILE,"\nHole_Faces: ");
     for (i=0;i<n_Hole_Faces;i++) Print_INCI(Hole_Faces[i]);
     fprintf(outFILE,"\n"); 
     fflush(0); */
  
  for (i=0;i<n_cd2_Faces;i++) cd2_Faces[i]=Hole_Faces[i];
  _P->np=n_old_v;
  while(n_Hole_Faces--){
    int not_found=1;
    VertexNumList Other_Verts;
    Equation Eq, E1, E2;
    /* Connect the codim 2 faces along the hole with other vertices
       bounding the hole; calculate corresponding equations and 
       check if they are good and new */
    
    for (j=0;(j<_E->ne)&&not_found;j++) 
      if (INCI_LE(Hole_Faces[n_Hole_Faces],E_INCI[j])){ 
	E1=_E->e[j]; not_found=0;}
    if (not_found) 
      fprintf(outFILE,"Equation1 not found in Close_The_Hole!\n");
    not_found=1;
    for (j=_E->ne;(j<old_ne)&&not_found;j++) 
      if (INCI_LE(Hole_Faces[n_Hole_Faces],E_INCI[j])){ 
	E2=_E->e[j]; not_found=0;}
    if (not_found) 
      fprintf(outFILE,"Equation2 not found in Close_The_Hole!\n");
    INCI_To_VertexNumList(INCI_XOR(Hole_Faces[n_Hole_Faces],Hole_Verts),
			  &Other_Verts, n_old_v);
    for (i=0;i<Other_Verts.nv;i++){
      Eq=EEV_To_Equation(&E1,&E2,_P->x[Other_Verts.v[i]],_P->n); 
      if (IsGoodCEq(&Eq,_P,_V)) {
	CEq_INCI[_CEq->ne]=Eq_To_INCI(&Eq,_P,_V);
	if (!Irrel(&(CEq_INCI[_CEq->ne]),E_INCI,&_E->ne)){
	  assert(_CEq->ne<CEQ_Nmax);
	  _CEq->e[_CEq->ne]=Eq;
	  Remove_INCI_From_List(CEq_INCI[_CEq->ne],
				Hole_Faces,&n_Hole_Faces); 
	  _CEq->ne++;
	  break;}}}}
  
  /* for (i=0;i<_CEq->ne;i++){
     fprintf(outFILE,"_CEq->e[%d]: ",i);
     for (j=0;j<_P->n;j++) fprintf(outFILE," %d",_CEq->e[i].a[j]); 
     fprintf(outFILE,"  %d\n",_CEq->e[i].c);}
     for (i=0;i<_E->ne;i++){
     fprintf(outFILE,"_E->e[%d]: ",i);
     for (j=0;j<_P->n;j++) fprintf(outFILE," %d",_E->e[i].a[j]); 
     fprintf(outFILE,"  %d\n",_E->e[i].c);}
     fflush(0); */
  
  for (i=0;i<_CEq->ne;i++) for (j=0;j<i;j++) {
    INCI New_Face=INCI_AND(CEq_INCI[i],CEq_INCI[j]);
    int k;
    if (INCI_abs(New_Face)<_P->n-1) continue;
    for (k=0;k<_CEq->ne;k++) if (INCI_LE(New_Face,CEq_INCI[k])) 
      if((k!=i)&&(k!=j)) break;
    if (k<_CEq->ne) continue;
    for (k=0;k<_E->ne;k++) if (INCI_LE(New_Face,E_INCI[k])) break;
    if (k<_E->ne) continue;
    assert(n_cd2_Faces<CD2F_Nmax);
    cd2_Faces[n_cd2_Faces++]=New_Face; }
  Bad_Vert_INCI=List_Complete(_P->n-1, _CEq->ne, n_cd2_Faces, 
		   CEq_INCI, cd2_Faces);
  if (!INCI_EQ_0(Bad_Vert_INCI)) 
    FE_Close_the_Hole(_P, _V, _E, _CEq, n_old_v, CEq_INCI, Bad_Vert_INCI);
  _P->np=np; 
}

int  Aided_IP_Check(PolyPointList *_P, VertexNumList *_V, EqList *_E, 
		    int n_irrel, int n_old_v){
     
  /* The first n_old_v entries of _P are the old vertices; 
     the first n_irrel entries of _E are the irrelevant old facets (which 
     remain facets), the remaining facets are the relevant facets */

  int i, j, n_Hole_Faces=0, old_ne=_E->ne;
  INCI CEq_INCI[CEQ_Nmax], Hole_Verts;
  INCI E_INCI[EQUA_Nmax];
  INCI Hole_Faces[CD2F_Nmax];
  CEqList CEq; 

  if (n_irrel==0) { 
    fprintf(outFILE, "n_irrel=0 in Aided_IP_Check!"); exit(0); }

  /* Create E_INCI: Incidences between _E and old vertices; 
     Intersect relevant and irrelevant facets to get the codim 2 faces
     that bound the hole; 
     Create Hole_Verts as the union of Hole_Faces */

  for (i=0;i<_E->ne;i++){
    E_INCI[i]=INCI_0();
    for (j=0;j<n_old_v;j++) E_INCI[i]=
      INCI_PN(E_INCI[i],Eval_Eq_on_V(&(_E->e[i]),_P->x[j],_P->n));} 

  for (i=0;i<n_irrel;i++) for (j=n_irrel;j<_E->ne;j++) {
    INCI New_Face=INCI_AND(E_INCI[i],E_INCI[j]);
    int k;
    if (INCI_abs(New_Face)<_P->n-1) continue;
    for (k=0;k<_E->ne;k++) if (INCI_LE(New_Face,E_INCI[k])) 
      if((k!=i)&&(k!=j)) break;
    if (k<_E->ne) continue;
    assert(n_Hole_Faces<CD2F_Nmax);
    Hole_Faces[n_Hole_Faces++]=New_Face;}

  /* puts("AIP:");
  Print_PPL(_P); 
  Print_EL(_E,&_P->n,0);
  puts("E_INCI:");
  for(i=0;i<_E->ne;i++) {Print_INCI(E_INCI[i]); puts("");}
  printf("n_irrel=%d\n", n_irrel); */

  Hole_Verts=INCI_0();
  for (i=n_irrel;i<_E->ne;i++) Hole_Verts=INCI_OR(Hole_Verts, E_INCI[i]); 
  for (i=0;i<n_old_v;i++) _V->v[i]=i;
  _V->nv=n_old_v;
  _E->ne=n_irrel;

  if (n_irrel==1) {
    /* Find vertex of maximal distance from the irrelevant facet and make
       Candidate equations by connecting it with the Hole_Faces */
    Long maxdist;
    int newvert=n_old_v;
    if (n_old_v==_P->np) return 0;
    maxdist=Eval_Eq_on_V(&(_E->e[0]),_P->x[n_old_v],_P->n);
    for (i=n_old_v+1; i<_P->np; i++){
      Long dist=Eval_Eq_on_V(&(_E->e[0]),_P->x[i],_P->n);
      if (dist>maxdist) {
	newvert=i;
	maxdist=dist;
	continue;}
      if (dist==maxdist) if (Vec_Greater_Than(_P->x[i],_P->x[newvert],_P->n)){
	newvert=i;
	maxdist=dist;}}
    _V->v[_V->nv++]=newvert;

    for (i=0;i<n_Hole_Faces;i++){
      Long dist;
      int not_found=1;
      Equation E1, E2;
      for (j=0;(j<n_irrel)&&not_found;j++) 
	if (INCI_LE(Hole_Faces[i],E_INCI[j])){ E1=_E->e[j]; not_found=0;}
      if (not_found) 
	fprintf(outFILE,"Equation1 not found in Aided_IP_Check!\n");
      not_found=1;
      for (j=n_irrel;(j<old_ne)&&not_found;j++) 
	if (INCI_LE(Hole_Faces[i],E_INCI[j])){ E2=_E->e[j]; not_found=0;}
      if (not_found) 
	fprintf(outFILE,"Equation2 not found in Aided_IP_Check!\n");
      CEq.e[i]=EEV_To_Equation(&E1,&E2,_P->x[newvert],_P->n); 
      for (j=0;j<n_old_v;j++) 
	if((dist=Eval_Eq_on_V(&(CEq.e[i]),_P->x[j],_P->n))){
	  if (dist<0){
	    int k;
	    CEq.e[i].c=-CEq.e[i].c; 
	    for(k=0;k<_P->n;k++) CEq.e[i].a[k]=-CEq.e[i].a[k]; 	} 
	  break;}}

    E_INCI[0]=INCI_PN(E_INCI[0],1);
    for (i=0;i<n_Hole_Faces;i++) CEq_INCI[i]=INCI_PN(Hole_Faces[i],0);

    CEq.ne=n_Hole_Faces;   }

  else {
    CEq.ne=0; 
#if SIMPLE_CTH
    FE_Close_the_Hole(_P,_V,_E,&CEq,n_old_v,CEq_INCI,Hole_Verts);}
#else
    Close_the_Hole(_P,_V,_E,&CEq,old_ne,n_old_v,n_Hole_Faces,
		   E_INCI, CEq_INCI, Hole_Verts, Hole_Faces);}
#endif

  /*Print_PPL(_P);
  Print_VL(_P,_V,"Verts:"); 
  for (i=0;i<CEq.ne;i++){
    fprintf(outFILE,"CEq.e[%d]: ",i);
    for (j=0;j<_P->n;j++) fprintf(outFILE," %d",CEq.e[i].a[j]); 
    fprintf(outFILE,"  %d\n",CEq.e[i].c);}
  fprintf(outFILE,"CEq_INCI: ");
  for (i=0;i<CEq.ne;i++) Print_INCI(CEq_INCI[i]);
  fprintf(outFILE,"\n"); 
  for (i=0;i<_E->ne;i++){
    fprintf(outFILE,"_E->e[%d]: ",i);
    for (j=0;j<_P->n;j++) fprintf(outFILE," %d",_E->e[i].a[j]); 
    fprintf(outFILE,"  %d\n",_E->e[i].c);}
  fprintf(outFILE,"E_INCI: ");
  for (i=0;i<_E->ne;i++) Print_INCI(E_INCI[i]); 
  fflush(0); */

  return Finish_IP_Check(_P, _V, _E, &CEq, E_INCI, CEq_INCI);
}

void Make_All_Subpolys(PolyPointList *_P, EqList *_E, VertexNumList *_V,
		       KeepList *_KL, NF_List *_NFL);
     
int Relevant(Long *X, int *n, EqList *_E, int *n_irrel){
  int i;
  for (i=0;i<*n_irrel;i++) 
    if (!Eval_Eq_on_V(&(_E->e[i]),X,*n)) return 0;
  for (i=0;i<*n;i++) if (X[i]) return 1;
  return 0;
}

int kept(int i, KeepList *_KL){
  int j;
  for(j=0;j<_KL->nk;j++) if(_KL->k[j]==i) return j;
  return -1;
}

void Drop_and_Keep(PolyPointList *_P, VertexNumList *_V, EqList *_E, 
		   KeepList *_KL, int drop_num,  NF_List *_NFL){
    
  int i, j, n_irrel=0, IP;
  int *new2old =(int *) malloc(sizeof(int)*POINT_Nmax);
  Long drop_point[POLY_Dmax];
  VertexNumList red_V;
  EqList *new_E=(EqList *) malloc(sizeof(EqList));
  PolyPointList *red_P = (PolyPointList *) malloc(sizeof (PolyPointList));

  if(red_P==NULL) {puts("Unable to allocate space for red_P"); exit(0);}
  if(new2old==NULL) {puts("Unable to allocate space for new2old"); exit(0);}

  for (j=0;j<_P->n;j++) drop_point[j]=_P->x[_V->v[drop_num]][j];

  /* Create new_E: Same as *_E, but irrelevant facets first */

  new_E->ne=_E->ne;
  for (i=0;i<_E->ne;i++){ 
    if (Eval_Eq_on_V(&(_E->e[i]),drop_point,_P->n)){
      for (j=0;j<_P->n;j++) new_E->e[n_irrel].a[j]=_E->e[i].a[j];
      new_E->e[n_irrel++].c=_E->e[i].c;}
    else {
      for (j=0;j<_P->n;j++) new_E->e[_E->ne+n_irrel-i-1].a[j]=_E->e[i].a[j];
      new_E->e[_E->ne+n_irrel-i-1].c=_E->e[i].c;}}

  /* Create red_P: old vertices, points that are not on irrelevant facets */
  red_P->np=0;
  red_P->n=_P->n;
  for (i=0;i<_V->nv;i++) if (i!=drop_num){
    for (j=0;j<_P->n;j++) red_P->x[red_P->np][j]=_P->x[_V->v[i]][j];
    new2old[red_P->np++]=_V->v[i];}
  for (i=0;i<_P->np;i++) 
    if ((i!=_V->v[drop_num])&&Relevant(_P->x[i],&_P->n,new_E,&n_irrel)){
      for (j=0;j<_P->n;j++) red_P->x[red_P->np][j]=_P->x[i][j];
      new2old[red_P->np++]=i;}

  /* Print_PPL(_P); 
  fprintf(outFILE,"new_E:\n");
  for (i=0;i<new_E->ne;i++){
    for (j=0;j<_P->n;j++) fprintf(outFILE,"%d ", new_E->e[i].a[j]);
    fprintf(outFILE," %d\n", new_E->e[i].c);}
  fprintf(outFILE,"\n"); 
  Print_PPL(red_P);
  fprintf(outFILE,"\n");
  printf("drop_point: ");
  for (i=0;i<_P->n;i++) printf("%d ", (int) drop_point[i]);
  puts("");*/

  _NFL->nIP++;

#if UnAided_IP_CHECK
  IP=IP_Check(red_P,&red_V,new_E);
#else
  IP=Aided_IP_Check(red_P,&red_V,new_E,n_irrel,_V->nv-1);
#endif

  /* puts("After AIP (in D&K):");
  Print_VL(red_P,&red_V);
  Print_EL(new_E,&red_P->n,0); */

  free(red_P);

  if(IP){

    VertexNumList new_V;
#if TEST_Aided_IP_CHECK
    VertexNumList test_V;
    EqList test_E;
#endif

    /* Create new_V from red_V */
    new_V.nv=red_V.nv;
    for (i=0;i<red_V.nv;i++) if((new_V.v[i]=new2old[red_V.v[i]])==_P->np-1) 
      new_V.v[i]=_V->v[drop_num];

    /* Drop the vertex _V->v[drop_num]: */
    _P->np--;
    for (j=0;j<_P->n;j++) _P->x[_V->v[drop_num]][j]=_P->x[_P->np][j]; 
    j=kept(_P->np,_KL);
    if (j>=0) _KL->k[j]=_V->v[drop_num];

#if TEST_Aided_IP_CHECK
  if(!IP_Check(_P,&test_V,&test_E)||(test_V.nv!=red_V.nv)||
      (test_E.ne!=new_E->ne)){
    int k;
    fprintf(outFILE,"_V: ");
    for (i=0;i<_V->nv;i++) fprintf(outFILE,"%d ",_V->v[i]);
    fprintf(outFILE,"\n"); 
    fprintf(outFILE,"_E:\n");
    for (i=0;i<_E->ne;i++){
      for (k=0;k<_P->n;k++) fprintf(outFILE,"%d ", (int) _E->e[i].a[k]);
      fprintf(outFILE," %d\n", (int) _E->e[i].c);}
    fprintf(outFILE,"\n");
    fprintf(outFILE,"drop_num: %d\n",drop_num);
    fprintf(outFILE,"drop_point: ");
    for (i=0;i<_P->n;i++) fprintf(outFILE,"%d ", (int) drop_point[i]);
    fprintf(outFILE,"\n");
    fprintf(outFILE,"_KL: ");
    for (i=0;i<_KL->nk;i++) fprintf(outFILE,"%d ", _KL->k[i]);
    fprintf(outFILE,"\n");

    fprintf(outFILE,"red_V: ");
    for (i=0;i<red_V.nv;i++) fprintf(outFILE,"%d ",red_V.v[i]);
    fprintf(outFILE,"\n"); 
    fprintf(outFILE,"new_E:\n");
    for (i=0;i<new_E->ne;i++){
      for (k=0;k<_P->n;k++) fprintf(outFILE,"%d ", (int) new_E->e[i].a[k]);
      fprintf(outFILE," %d\n", new_E->e[i].c);}
    fprintf(outFILE,"\n");
    fprintf(outFILE,"n_irrel: %d\n",n_irrel);
    fprintf(outFILE,"new2old:");
    for (i=0;i<red_P->np;i++) fprintf(outFILE," %d", new2old[i]);
    fprintf(outFILE,"\n");

    fprintf(outFILE,"test_V: ");
    for (i=0;i<test_V.nv;i++) fprintf(outFILE,"%d ",test_V.v[i]);
    fprintf(outFILE,"\n"); 
    fprintf(outFILE,"test_E:\n");
    for (i=0;i<test_E.ne;i++){
      for (k=0;k<_P->n;k++) fprintf(outFILE,"%d ", (int) test_E.e[i].a[k]);
      fprintf(outFILE," %d\n", (int) test_E.e[i].c);}
    fprintf(outFILE,"\n");

    exit(0);} 
#endif

    Make_All_Subpolys(_P, new_E, &new_V, _KL, _NFL);
    
    /* Reconstruct _P  */
    if (j>=0) _KL->k[j]=_P->np;
    for (j=0;j<_P->n;j++){
      _P->x[_P->np][j]=_P->x[_V->v[drop_num]][j];
      _P->x[_V->v[drop_num]][j]=drop_point[j];}
    _P->np++;}
  
  /* Attach a keep label to the vertex _V->v[drop_num]: */
  _KL->k[_KL->nk]=_V->v[drop_num];
  _KL->nk++; 
  free(new2old); free(new_E);
}

  

void Start_Make_All_Subpolys(PolyPointList *_P, NF_List *_NFL){
  int i,j;
  VertexNumList V, new_V; 
  EqList E, new_E; 
  KeepList KL; 
  Long drop_point[POLY_Dmax];
  _NFL->nIP++; 

  if (_NFL->rf) _NFL->rf=_NFL->rd;
  _NFL->rd=0;

  if(!IP_Check(_P,&V,&E)) {
    fprintf(outFILE,"IP_check negative!\n"); return;}

  if(_NFL->Nmin==0) {PairMat PM; Make_VEPM(_P,&V,&E,PM); /* complete poly if */
    Complete_Poly(PM, &E, V.nv, _P); _NFL->Nmin=_P->np;} /* non-CWS input   */

  /* fprintf(outFILE,"StartMASP:\n");
  for (i=0;i<V.nv;i++){
    for (j=0;j<_P->n;j++) fprintf(outFILE,"%d ",(int) _P->x[V.v[i]][j]);
    fprintf(outFILE,"\n");}
  for (i=0;i<E.ne;i++){
    for (j=0;j<_P->n;j++) fprintf(outFILE,"%d ", (int) E.e[i].a[j]);
    fprintf(outFILE," %d\n", (int) E.e[i].c);}
  fflush(0);*/

  KL.nk=0;

  if(_NFL->kf){
    fprintf(outFILE,"The vertices are:\n");
    for (i=0;i<V.nv;i++){
      fprintf(outFILE,"%d.  ",i);
      for (j=0;j<_P->n;j++) fprintf(outFILE,"%d ",(int) _P->x[V.v[i]][j]);
      fprintf(outFILE,"\n"); }
    fprintf(outFILE,"How many of them do you want to keep?\n");
    fscanf(inFILE,"%d",&(KL.nk));
    fprintf(outFILE,"Which %d of them do you want to keep?\n",KL.nk);
    for (i=0; i<KL.nk;i++){
      fscanf(inFILE,"%d",&j);
      KL.k[i]=V.v[j]; } 
    fprintf(outFILE,"Keeping\n");
    for (i=0; i<KL.nk;i++) {
      for(j=0;j<_P->n;j++) fprintf(outFILE,"%d ",(int) _P->x[KL.k[i]][j]);
    fprintf(outFILE,"\n"); }}

  else for(i=0;i<V.nv;i++) {
    /* Try dropping the vertex V.v[i]: */
    _P->np--;
    for (j=0;j<_P->n;j++){
      drop_point[j]=_P->x[V.v[i]][j];
      _P->x[V.v[i]][j]=_P->x[_P->np][j]; }
    if (!IP_Check(_P,&new_V,&new_E)) KL.k[KL.nk++]=V.v[i];
    /* Reconstruct _P  */
    for (j=0;j<_P->n;j++){
      _P->x[_P->np][j]=_P->x[V.v[i]][j];
      _P->x[V.v[i]][j]=drop_point[j];}
    _P->np++; }

  Make_All_Subpolys(_P, &E, &V, &KL, _NFL);
}


void IREgcd(Long *vec_in, int *d, Long *vec_out){
  int i, j, n_rem_perm, n_perm=1, a, b, perm_j, inv_perm[POLY_Dmax];
  Long perm_vec_in[POLY_Dmax], perm_vec_out[POLY_Dmax], abs_sum, max_abs_sum=0;
  for (i=1;i<=*d;i++) n_perm*=i;
  for (i=0;i<n_perm;i++){
    b=i;
    n_rem_perm=n_perm;
    for (j=0;j<*d;j++) inv_perm[j]=-1;
    for (j=*d;j>0;j--){
      n_rem_perm/=j;
      a=b/n_rem_perm;
      b=b%n_rem_perm;
      perm_j=*d;
      while (a>=0){
	perm_j--;
	if (inv_perm[perm_j]<0) a--;}
      perm_vec_in[j-1]=vec_in[perm_j];
      inv_perm[perm_j]=j-1;}
    REgcd(perm_vec_in, d, perm_vec_out);
    abs_sum=0;
    for (j=0;j<*d;j++) abs_sum+=labs(perm_vec_out[j]);
    if ((i==0)||(abs_sum<max_abs_sum)){
      max_abs_sum=abs_sum;
      for (j=0;j<*d;j++) vec_out[j]=perm_vec_out[inv_perm[j]];    }  }
}


int Make_RedVec(Long n, Long *badFE, Long *RedVec){
  int i=0, j;
  Long RV[POLY_Dmax], FE[POLY_Dmax], dist=0;
  for (j=0;j<n;j++) if (badFE[j]) FE[i++]=badFE[j];
  if(i==1) RV[0]=1;
  else IREgcd(FE,&i,RV);
  for (j=0;j<i;j++) dist+=FE[j]*RV[j];
  if ((dist!=1)&&(dist!=-1)) {
    fprintf(outFILE,"FE: ");
    for (j=0;j<i;j++) fprintf(outFILE,"%d ",(int) FE[j]);
    fprintf(outFILE,"\n");     fprintf(outFILE,"RV: ");
    for (j=0;j<i;j++) fprintf(outFILE,"%d ",(int) RV[j]);
    fprintf(outFILE,"\n");
    fprintf(outFILE,"dist= %d\n", (int) dist);
    return 0;}
  i=0;
  for (j=0;j<n;j++) if (badFE[j]) RedVec[j]=-dist*RV[i++];
  else RedVec[j]=0;
  return 1;
}

void Reduce_Poly(PolyPointList *_P, EqList *_E, KeepList *_KL, 
		 NF_List *_NFL, int badfacet){

  /* reduce _P to the sublattice where badfacet has distance 1;
     call Make_All_Subpolys for the reduced polyhedron new_P;   */

  VertexNumList V;
  PolyPointList *new_P= (PolyPointList *) malloc(sizeof (PolyPointList));
  KeepList new_KL;
  int i, j;
  Long RedVec[POLY_Dmax];
  Equation BE=_E->e[badfacet];
  Rat VPM_checksum_old, VPM_checksum_new;

  /*Generate new_P (in old coordinates!): */

  new_P->np=0;
  new_KL.nk=0;
  new_P->n=_P->n;
  for (i=0;i<_P->np;i++){
    Long dist=0;
    for (j=0;j<_P->n;j++) dist-=BE.a[j]*_P->x[i][j];  
    if (!(dist%BE.c)){
      for (j=0;j<_P->n;j++) 
	new_P->x[new_P->np][j]=_P->x[i][j];
      if (kept(i,_KL)>=0) new_KL.k[new_KL.nk++]=new_P->np;
      new_P->np++; } }

  /*   Print_PPL(new_P); 
  fprintf(outFILE,"kept: ");
  for (j=0;j<new_KL.nk;j++) fprintf(outFILE,"%d ",new_KL.k[j]);
  fprintf(outFILE,"\n\n"); */

  if(new_KL.nk!=_KL->nk) {free(new_P); return;}
  _NFL->nIP++; 

  if(!IP_Check(new_P,&V,_E)) {free(new_P); return;}

  VPM_checksum_old=rI(0);
  for (i=0;i<_E->ne;i++) {
    Long s=0;
    for (j=0;j<V.nv;j++) 
     s+=Eval_Eq_on_V(&(_E->e[i]),new_P->x[V.v[j]], new_P->n);
    VPM_checksum_old=rS(VPM_checksum_old, rQ(rI(s),rI(_E->e[i].c)));}

  /* Generate RedVec: */
  if (!Make_RedVec(_P->n, BE.a, RedVec)){
    Print_PPL(_P,"");
    fprintf(outFILE,"Bad Facet: ");
    for (j=0;j<_P->n;j++) fprintf(outFILE,"%d ",(int) BE.a[j]);
    fprintf(outFILE," %d\n",(int) BE.c);
    fprintf(outFILE,"kept: ");
    for (j=0;j<_KL->nk;j++) fprintf(outFILE,"%d ",_KL->k[j]);
    fprintf(outFILE,"\n"); 
    fprintf(outFILE,"RedVec: ");
    for (j=0;j<_P->n;j++) fprintf(outFILE,"%d ",(int) RedVec[j]);
    fprintf(outFILE,"\n"); 
    exit(0);}

  /*Coordinates of new_P: */

  for (i=0;i<new_P->np;i++){
    Long dist=0;
    for (j=0;j<new_P->n;j++) dist-=BE.a[j]*new_P->x[i][j];  
    if (!(dist%BE.c)){
      for (j=0;j<new_P->n;j++) 
	new_P->x[i][j]-=(dist-dist/BE.c)*RedVec[j]; } }

  /*   Print_PPL(new_P); 
  fprintf(outFILE,"\n\n"); */

  if(!New_Improve_Coords(new_P,&V))   {    
    Print_PPL(_P,"");
    fprintf(outFILE,"Bad Facet: ");
    for (j=0;j<_P->n;j++) fprintf(outFILE,"%d ",(int) BE.a[j]);
    fprintf(outFILE," %d\n",(int) BE.c);
    fprintf(outFILE,"kept: ");
    for (j=0;j<_KL->nk;j++) fprintf(outFILE,"%d ",_KL->k[j]);
    fprintf(outFILE,"\n"); 
    fprintf(outFILE,"RedVec: ");
    for (j=0;j<_P->n;j++) fprintf(outFILE,"%d ",(int) RedVec[j]);
    fprintf(outFILE,"\n"); 
    puts("WARNING: INCOMPLETE SL REDUCTION"); fflush(stdout);
#if	INCOMPLETE_SL_REDUCTION
    free(new_P); return;
#else
    exit(0);
#endif
  }

  if(!IP_Check(new_P,&V,_E)) {
    fprintf(outFILE,"Trouble in Reduce_Poly!"); exit(0);}

  VPM_checksum_new=rI(0);
  for (i=0;i<_E->ne;i++) {
    Long s=0;
    for (j=0;j<V.nv;j++) 
     s+=Eval_Eq_on_V(&(_E->e[i]),new_P->x[V.v[j]],new_P->n);
    VPM_checksum_new=rS(VPM_checksum_new, rQ(rI(s),rI(_E->e[i].c)));}
  if(VPM_checksum_new.D*VPM_checksum_old.N-
     VPM_checksum_new.N*VPM_checksum_old.D) {
         fprintf(outFILE,"Checksums don't match in Reduce_Poly!"); exit(0); }

  Make_All_Subpolys(new_P,_E,&V,&new_KL, _NFL); 
  free(new_P);
}


void Make_All_Subpolys(PolyPointList *_P, EqList *_E, VertexNumList *_V,
		       KeepList *_KL, NF_List *_NFL){

  /* Creates all subpolyhedra of an IP-polyhedron given by _P, _E, _V.
     The basic structure is:
     Search for bad facet;
     If (reflexive) if (!Add_NF_to_List) return;
     Drop_and_Keep vertices (on bad facet / everywhere);
     If (bad facet exists) Reduce_Poly;  */

  int i, j, badfacet=-1, nbadvert=0, maxnkept=-1;
  KeepList new_KL=*_KL;

  /* fprintf(outFILE,"MASP:\n");
     Print_PPL(_P); 
     Print_VL(_P,_V);
     Print_EL(_E,&_P->n,0);
     fflush(0);*/
  
  if (_V->nv>_NFL->VN) _NFL->VN=_V->nv;
  if (_E->ne>_NFL->FN) _NFL->FN=_E->ne;

  /* Choose bad facet with largest number of kept vertices: */

  for (i=0;i<_E->ne;i++) if (_E->e[i].c>1){
    int nkept=0;
    for (j=0;j<_KL->nk;j++) 
      if (!Eval_Eq_on_V(&(_E->e[i]),_P->x[_KL->k[j]],_P->n)) nkept++;
    if (nkept>maxnkept){
      maxnkept=nkept; 
      badfacet=i; }}

  if  (!_NFL->rf) if (badfacet==-1){	/* not recover mode && reflexive */
    if (_NFL->of<=-2) {			/* flags -o0  (break on reflexive) */
      if (_NFL->rd) {			/*    or -oc  (recover missing)	*/
        if(_NFL->of==-2){Add_NF_to_List(_P,_V,_E,_NFL); return;} /* flag -o0 */
        else {if(!Add_NF_to_List(_P,_V,_E,_NFL)) return;}	 /* flag -oc */
        }}
    else if(!Add_NF_to_List(_P,_V,_E,_NFL)) return;}		 /* standard */

  if(0 < _NFL->of) if(_NFL->of <= _NFL->rd) return; /* limit rd */

  if (badfacet!=-1) {
    /* Find vertices on bad facet: */
    nbadvert=0;
    for (i=0;i<_V->nv;i++) 
      if (!Eval_Eq_on_V(&(_E->e[badfacet]),_P->x[_V->v[i]],_P->n)){
	j=_V->v[nbadvert];
	_V->v[nbadvert]=_V->v[i];
	_V->v[i]=j;
	nbadvert++; } }

  if (_NFL->rf&&(_NFL->rf<_NFL->rd)) fprintf(outFILE,"_NFL->rf<_NFL->rd\n");
  if (_NFL->rf==_NFL->rd) _NFL->rf=0;
  if (_NFL->rf) for (i=0;i<_NFL->b[_NFL->rd];i++) 
    if((kept(_V->v[i],&new_KL)<0)) new_KL.k[new_KL.nk++]=_V->v[i];
  _NFL->rd++;
  
  for(i=0;i<((badfacet==-1)?_V->nv:nbadvert);i++) 
    if((kept(_V->v[i],&new_KL)<0)){
      _NFL->b[_NFL->rd-1]=i;
      Drop_and_Keep(_P, _V, _E, &new_KL, i, _NFL);    }
  
  /* If a bad facet is kept, reduce the poly: */
  if((_NFL->of==0)||(_NFL->of==-3)) if (badfacet!=-1) {/* -oc recover[of=-3] */
    int currentSL=_NFL->SL; 
    _NFL->SL=1;
    _NFL->b[_NFL->rd-1]=nbadvert;
    Reduce_Poly(_P, _E, &new_KL, _NFL, badfacet); 
    _NFL->SL=currentSL;   }

  _NFL->rd--;
}
 
void Ascii_to_Binary(CWS *W, PolyPointList *P, 
  char *dbin, char *polyi, char *polyo){
  NF_List *_NFL=(NF_List *) malloc(sizeof(NF_List)); 
  VertexNumList V; 
  EqList F;
  assert(_NFL!=NULL);
  if(!(*polyo)) {
    puts("You have to specify an output file via -po in -a-mode!\n"); 
    printf("For more help type use option `-h'\n");
    exit(0);}
  _NFL->of=0; _NFL->rf=0;	
  _NFL->iname=polyi; _NFL->oname=polyo; _NFL->dbname=dbin; 
  Init_NF_List(_NFL);
  _NFL->SL=0;
  
  while(Read_CWS_PP(W,P))    {
    _NFL->hc = 0;
    _NFL->V= _NFL->F= _NFL->VN= _NFL->FN = _NFL->Xnuc = _NFL->Xdif = 0;	
    _NFL->Nmin=P->np; _NFL->Nmax=0; 
    if(_NFL->d==0) _NFL->d=P->n;
    else if(_NFL->d-P->n) {puts("different dim!"); exit(0);}
    if (!IP_Check(P,&V,&F)){
      puts("IP_Check failed in Ascii_to_Binary!\n"); exit(0);}
    if (Add_NF_to_List(P,&V,&F,_NFL)) if (outFILE!=stdout){
      int i, j;
      for(i=0;i<W->nw;i++)     {	
	fprintf(outFILE,"%d ",(int) W->d[i]);
	for(j=0;j<W->N;j++) fprintf(outFILE,"%d ",(int) W->W[i][j]);
	if(i+1<W->nw) fprintf(outFILE," "); else fprintf(outFILE,"\n");   } 
      fflush(0);}  }
  Write_List_2_File(polyo,_NFL); 
  free(_NFL);
}

void Do_the_Classification(CWS *W, PolyPointList *P, /* char *fn, */
  int oFlag, int rFlag, int kFlag, char *polyi, char *polyo, char *dbin) {

  /* static int nw; */
  NF_List *_NFL=(NF_List *) malloc(sizeof(NF_List)); 
  time_t W_SAVE_TIME=time(NULL); 
  if(!(*polyo)) { 
    puts("You have to specify an output file via -po!\n"); 
    printf("For more help use option '-h'\n");
    exit(0);}
  assert(_NFL!=NULL); _NFL->of=oFlag; _NFL->rf=rFlag; _NFL->kf=kFlag; 
  _NFL->iname=polyi; _NFL->oname=polyo; _NFL->dbname=dbin; 
  Init_NF_List(_NFL);
  rFlag=0; /* now used as "read flag" */

  while(Read_CWS_PP(W,P)) {		   /* make subpolys */
    if(W->nw>0) _NFL->Nmin=P->np; else _NFL->Nmin=0 /* :: Complete_Poly() */ ;
    _NFL->Nmax=0; _NFL->SL = _NFL->hc = 0;
    _NFL->V= _NFL->F= _NFL->VN= _NFL->FN = _NFL->Xnuc = _NFL->Xdif = 0;	
    if(_NFL->d==0) _NFL->d=P->n;
    else if(_NFL->d-P->n) {puts("different dim!"); exit(0);}
    if(rFlag) {Read_File_2_List(polyo,_NFL); rFlag=0;}
    Start_Make_All_Subpolys(P, _NFL); 
    Print_Weight_Info(W,_NFL);
    if((WRITE_DIM <= P->n) && (MIN_NEW <= _NFL->NP))
    if((int)difftime(time(NULL),W_SAVE_TIME) > MIN_W_SAVE_TIME)   {   
      Write_List_2_File(polyo,_NFL); 
      rFlag=1; 
      _NFL->SAVE=W_SAVE_TIME=time(NULL);    }  }

  if(rFlag==0) Write_List_2_File(polyo,_NFL); 
  _NFL->TIME=time(NULL); 
  fputs(ctime(&_NFL->TIME),stdout);
  free(_NFL);
}


/*  ======================================================================  */
/*  ==========		     			  		==========  */
/*  ==========		S U B L A T T I C E S  			==========  */
/*  ==========							==========  */
/*  ======================================================================  */

void Write_VPM(FILE *F,int *v, int *f, subl_int x[][VERT_Nmax]){    
  int i,j; 
  fprintf(F,"%d %d\n",*v,*f);
  for(i=0;i<*f;i++)     {   
    for(j=0;j<*v-1;j++) fprintf(F,"%lld ",(long long) x[i][j]);
    fprintf(F,"%lld\n",(long long) x[i][j]);  }
}

void Write_M(FILE *F,int *v, int *f, subl_int x[][POLY_Dmax])
{    int i,j; fprintf(F,"%d %d\n",*v,*f);
     for(i=0;i<*f;i++)
     {   for(j=0;j<*v-1;j++) fprintf(F,"%lld ",(long long) x[i][j]);
         fprintf(F,"%lld\n",(long long) x[i][j]);
     }
}


void Make_All_Sublat(NF_List *_L, int n, int v, subl_int diag[POLY_Dmax], 
		     subl_int u[][VERT_Nmax], char *mFlag, PolyPointList *_P){
  /* create all inequivalent decompositions diag=s*t  into upper
     triangular matrices s,t;
     t*(first lines of u) becomes the poly P;
     choose the elements of t (columns rising, line # falling),
     calculate corresponding element of s at the same time;
     finally calculate P and add to list  */

   
  VertexNumList V; 
  EqList F;
  subl_int s[POLY_Dmax][POLY_Dmax], t[POLY_Dmax][POLY_Dmax], ts_lincol;
  int i,j,k, lin=0, col=0, err=0;

  _L->d=_P->n=n;
  _P->np=V.nv=v;
  for(i=0;i<n;i++) for (j=0;j<n;j++) t[i][j]=(i>=j)-1;
  while (col>=0){
    /* fprintf(outFILE,"t[%d][%d]=%lld\n",lin, col, t[lin][col]); */
    if (lin==col){
      do 
	t[lin][col]++; 
      while ((diag[lin]%t[lin][col])&&(t[lin][col]<=diag[lin]));
      if (t[lin][col]>diag[lin]){
	if (col==0) return;
	t[lin][col]=0;
	lin=0;
	col--;}
      else{ 
	s[lin][col]=diag[lin]/t[lin][col];
	if (col==0) lin=col=1; 
	else lin--; } }
    else{
      /* (t*s)[lin][col]=sum_j t[lin][j]*s[j][col]=0 */
      ts_lincol=0;
      for (j=col;j>lin;j--) ts_lincol+=t[lin][j]*s[j][col];
      do{
	t[lin][col]++;
	ts_lincol+=s[col][col];}
      while ((ts_lincol%t[lin][lin])&&(t[lin][col]<t[lin][lin]));
      if (t[lin][col]>=t[lin][lin]){
	t[lin][col]=-1;
	lin++;}
      else{
	s[lin][col]=-ts_lincol/t[lin][lin];
	if (lin) lin--;
	else if (col<n-1) {col++; lin=col;}
	else{
	  for (i=0;i<n;i++) for (j=0;j<v;j++){
	    /* _P->x[j][i]=s[i][i]*u[i][j]; 
	       for (k=i+1;k<n;k++) _P->x[j][i]+=s[i][k]*u[k][j]; } */
	    subl_int Pji=s[i][i]*u[i][j]; 
	    for (k=i+1;k<n;k++) Pji+=s[i][k]*u[k][j]; 
	    _P->x[j][i]=Pji;
	    if (labs(Pji)>20000000) err=1;}
	  _P->np=v;
	  if (err){
	    Print_PPL(_P,"");
	    fprintf(outFILE,"u: "); Write_VPM(outFILE,&v,&v,u);
	    fprintf(outFILE,"s: "); Write_M(outFILE,&n,&n,s);
	    fprintf(outFILE,"t: "); Write_M(outFILE,&n,&n,t);}
	  for(i=0;i<v;i++) V.v[i]=i;
	  New_Improve_Coords(_P,&V);

	  if (!IP_Check(_P,&V,&F)){
	    puts("IP_Check failed in Make_All_Sublat!\n"); exit(0);}
	  for (i=0;i<F.ne;i++) if(F.e[i].c!=1) {
	    fprintf(outFILE,"Not reflexive in Make_All_Sublat!\n"); 
	    exit(0);}
	  if (*mFlag != 'r') Add_NF_to_List(_P,&V,&F,_L);  
	  else if (Poly_Max_check(_P,&V,&F)) if(Add_NF_to_List(_P,&V,&F,_L)){
		Print_PPL(_P,"");		}}}}}
}

subl_int SI_abs(subl_int x){return (x>0 ? x : -x);}

int Line_Recomb(int i, int line, int v, int f, subl_int x[][VERT_Nmax]){
  int j, k, changes=0; 
  subl_int y11, y12, y21, y22, g, new_;	
  if (SI_abs(x[i][i])==1) return 0;
  for (j=line;j<f;j++){
    /* take a combination of i'th and j'th line such that x[i][i] becomes
       gcd(x[i][i],x[j][i]) */
    if (!x[j][i]) continue;
    if (!x[i][i]||(x[j][i]%x[i][i])){
      changes=1;
      if ((g=LEgcd(x[i][i],x[j][i],&y11,&y12))<0) {y11*=-1; y12*=-1; g*=-1;}
      y21=-x[j][i]/g;
      y22=x[i][i]/g;
      for (k=i;k<v;k++){
	new_=y11*x[i][k]+y12*x[j][k];
	x[j][k]=y21*x[i][k]+y22*x[j][k];
	x[i][k]=new_; } }  }
  return changes;
}

void Col_Recomb(int i, int col, int v, int f, subl_int x[][VERT_Nmax], 
		subl_int u[][VERT_Nmax] ){
    /* take a combination of i'th and col'th column such that x[i][i] becomes
       gcd(x[i][i],x[i][col]); change u accordingly */
  int k; 
  subl_int y11, y12, y21, y22, g, new_;	/* "new" conflicts with g++ => new_ */
  if ((g=LEgcd(x[i][i],x[i][col],&y11,&y12))<0) {y11*=-1; y12*=-1; g*=-1;}
  y21=-x[i][col]/g;
  y22=x[i][i]/g;
  for (k=i;k<f;k++){
    new_=y11*x[k][i]+y12*x[k][col];
    x[k][col]=y21*x[k][i]+y22*x[k][col];
    x[k][i]=new_; }  
  for (k=0;k<v;k++){
    new_=y22*u[i][k]-y21*u[col][k];
    u[col][k]=-y12*u[i][k]+y11*u[col][k];
    u[i][k]=new_; }  
}

void lin_col_swap(int i, int lin, int col, int v, int f, 
		  subl_int x[][VERT_Nmax], subl_int u[][VERT_Nmax] ){
  /* swap i'th and lin'th line, i'th and col'th column, change u accordingly */
  int k;
  subl_int swap_vec[VERT_Nmax];
  if (lin!=i) for (k=i;k<v;k++){
    swap_vec[k]=x[lin][k];
    x[lin][k]=x[i][k];
    x[i][k]=swap_vec[k];}
  if (col!=i){
    for (k=i;k<f;k++){
      swap_vec[k]=x[k][col];
      x[k][col]=x[k][i];
      x[k][i]=swap_vec[k];}
    for (k=0;k<v;k++){
      swap_vec[k]=u[col][k];
      u[col][k]=u[i][k];
      u[i][k]=swap_vec[k];}}
}


void MakePolyOnSublat(NF_List *_L, subl_int x[VERT_Nmax][VERT_Nmax],
		      int v, int f, int *max_order, char *mFlag, 
		       PolyPointList *_P){

  /* Decompose the VPM x as x=w*diag*u, where w and u are SL(Z);
     the first lines of u are the coordinates on the coarsest lattices */    

  int i, j, k, col, lin;
  subl_int diag[POLY_Dmax];
  subl_int fac, order=1;
  subl_int u[VERT_Nmax][VERT_Nmax];

  for (i=0;i<v;i++) for (j=0;j<v;j++) u[i][j]=(i==j);
 
  for (i=0;i<VERT_Nmax;i++) {
    /* put gcd of remaining matrix elements at x[i][i], make zeroes at
       x[i][j] and x[j][i] for all j unequal to i   */
    int chosen_lin=0, chosen_col=-1;
    subl_int min_entry, min_col_max=0, min_lin_gcd=0, min_lin_max=0, rem_gcd=0;

    for (lin=i;lin<f;lin++){
      subl_int lin_gcd=0, lin_max=0;
      for (col=i;col<v;col++) if(x[lin][col]){
	if (!lin_gcd) lin_gcd=SI_abs(x[lin][col]);
	else lin_gcd=NNgcd(lin_gcd,x[lin][col]);
	if (SI_abs(x[lin][col])>lin_max) lin_max=SI_abs(x[lin][col]);}
      if (lin_max) if ((!min_lin_gcd)||(lin_gcd<min_lin_gcd)||
		       ((lin_gcd==min_lin_gcd)&&(lin_max<min_lin_max))){
	chosen_lin=lin;
	min_lin_gcd=lin_gcd;
	min_lin_max=lin_max;}
      rem_gcd=NNgcd(rem_gcd, lin_gcd);}

    if (!min_lin_max) break;

    min_entry=min_lin_max+1;
    for (col=i;col<v;col++) 
      if(x[chosen_lin][col]) if(SI_abs(x[chosen_lin][col])<=min_entry){
	subl_int col_max=0;
	for (lin=i;lin<f;lin++) if (SI_abs(x[lin][col])>col_max) 
	  col_max=SI_abs(x[lin][col]);
	if (!col_max){printf("col_max==0!!!"); exit(0);}
	if ((chosen_col==-1)||(SI_abs(x[chosen_lin][col])<min_entry)||
	    ((SI_abs(x[chosen_lin][col])==min_entry)&&(col_max<min_col_max))){
	  chosen_col=col;
	  min_col_max=col_max;}}
    
    lin_col_swap(i, chosen_lin, chosen_col, v, f, x, u);

    /* put gcd of remaining matrix elements at x[i][i] */
    while(SI_abs(x[i][i])>rem_gcd){
      for (col=i+1;col<v;col++) if ((x[i][col]%x[i][i])) 
	Col_Recomb(i,col,v,f,x,u);
      if (Line_Recomb(i,i+1,v,f,x)) continue;
      if (SI_abs(x[i][i])<=rem_gcd) break;
      for (col=i+1;col<v;col++){
	for(j=i+1;(j<f)&&(!(x[j][col]%x[i][i]));j++);
	if (j<f){
	  /* set i'th element of j'th line to zero, then add j'th line to
	     i'th line, change w accordingly */
	  fac=x[j][i]/x[i][i];
	  for (k=i;k<v;k++){
	    x[j][k]-=fac*x[i][k];
	    x[i][k]+=x[j][k];}
	  break;}}}

    /* make zeroes */
    for (col=i+1;col<v;col++) {
      fac=x[i][col]/x[i][i];
      for (k=i;k<f;k++) x[k][col]-=fac*x[k][i];
      for (k=0;k<v;k++) u[i][k]+=fac*u[col][k];    }
    for (k=i+1;k<v;k++) if (x[i][k]){
      printf("error in MakePolyOnSublat!!!\n"); exit(0);}}

  for (j=0;j<i;j++) order*=(diag[j]=SI_abs(x[j][j]));
  if (order>*max_order) *max_order=order;
  if (i>POLY_Dmax) {printf("diag has %d entries!!!\n", i); exit(0);}
  if (order>1) Make_All_Sublat(_L, i, v, diag, u, mFlag, _P);
}

void uc_nf_to_P(PolyPointList *_P, int *MS, int *d, int *v, int *nuc, 
                unsigned char *uc){    
  Long tNF[POLY_Dmax][VERT_Nmax]; int i, j; 
  UCnf2vNF(d, v, nuc, uc, tNF, MS); 
  _P->n=*d; _P->np=*v; (*MS)%=4;
  for(i=0;i<*v;i++) for(j=0;j<*d;j++) _P->x[i][j]=tNF[j][i];
}

void Find_Sublat_Polys(char mFlag, char *dbin, char *polyi, char *polyo, 
		       PolyPointList *_P){	
  NF_List *_NFL=(NF_List *) malloc(sizeof(NF_List)); 
  VertexNumList Vnl; EqList Fel;
  subl_int x[VERT_Nmax][VERT_Nmax], y[VERT_Nmax][VERT_Nmax];
  time_t Tstart;
  int v, nu, i, j, k, max_order=1;

  assert(_NFL!=NULL);
  if(!(*polyo)) {
    puts("You have to specify an output file via -po in -sm-mode."); 
    printf("For more help use option `-h'\n");
    exit(0);}
  _NFL->of=0; _NFL->rf=0;
  _NFL->iname=polyi; _NFL->oname=polyo; _NFL->dbname=dbin;
  Init_NF_List(_NFL);
  _NFL->SL=0;

  if (*dbin){
    DataBase *DB=&(_NFL->DB);
    char *dbname = (char *) malloc(1+strlen(dbin)+File_Ext_NCmax), *fx;
    unsigned char uc_poly[NUC_Nmax];
    int MS;
    
    strcpy(dbname,dbin);
    strcat(dbname,".");
    fx=&dbname[strlen(dbin)+1]; 
    
    printf("Reading DB-files, calculating sublattices:\n"); fflush(0);
    
    /* read the DB-files and calculate Sublattices */
    for (v=2;v<=DB->nVmax;v++) if(DB->nNUC[v]){
      FILE *dbfile;
      char ext[4]={'v',0,0,0};
      ext[1]='0' + v / 10; ext[2]='0' + v % 10;
      Tstart=time(NULL); 
      strcpy(fx,ext); 
      dbfile=fopen(dbname,"rb"); 
      assert(dbfile!=NULL); 
      if (!mFlag) {printf("Reading %s\n", dbname); fflush(0);}
      for (nu=0;nu<=DB->NUCmax;nu++) for (i=0; i<DB->NFnum[v][nu]; i++){
	for (j=0; j<nu; j++) uc_poly[j]=fgetc(dbfile);
	uc_nf_to_P(_P, &MS, &(_NFL->d), &v, &nu, uc_poly);
	if (MS>3) {printf("MS=%d!!!\n",MS); exit(0);}
	assert(IP_Check(_P,&Vnl,&Fel));
	assert(v==Vnl.nv);
	/* compute VPM */
	for(j=0;j<Fel.ne;j++) for(k=0;k<Vnl.nv;k++)
	  x[j][k]=Eval_Eq_on_V(&(Fel.e[j]),_P->x[Vnl.v[k]],_P->n)-1;
	for(j=0;j<Fel.ne;j++) for(k=0;k<Vnl.nv;k++) y[k][j]=x[j][k];
	if (MS!=2) 
	  MakePolyOnSublat(_NFL, x, Vnl.nv, Fel.ne, &max_order, &mFlag, _P); 
	if (MS>1) 
	  MakePolyOnSublat(_NFL, y, Fel.ne, Vnl.nv, &max_order, &mFlag, _P); }

      if(ferror(dbfile)) {printf("File error in %s\n",dbname); exit(0);}
      fclose(dbfile);
      printf(" %dp (%ds)\n", (int)_NFL->NP, (int) difftime(time(NULL),Tstart)); 
      fflush(0);  }}
  
  else { 
    CWS W;
    while(Read_CWS_PP(&W,_P)){
      assert(IP_Check(_P,&Vnl,&Fel));
      /* compute VPM */
      for(j=0;j<Fel.ne;j++) for(i=0;i<Vnl.nv;i++)
	x[j][i]=Eval_Eq_on_V(&(Fel.e[j]),_P->x[Vnl.v[i]],_P->n)-1;
      MakePolyOnSublat(_NFL, x, Vnl.nv, Fel.ne, &max_order, &mFlag, _P); 
      if (!mFlag) Print_Weight_Info(&W,_NFL);}}
    
  printf("max_order=%d\n", max_order);
  Write_List_2_File(polyo,_NFL); 
  _NFL->TIME=time(NULL); fputs(ctime(&_NFL->TIME),stdout);
  free(_NFL);
}

 

/*  ======================================================================  */
/*  ==========		     			  		==========  */
/*  ==========		   I R ,   V I R   &   M A X		==========  */
/*  ==========							==========  */
/*  ======================================================================  */

int irred(PolyPointList *_P)
{
  int i,j,drop_point[POLY_Dmax];
  VertexNumList V,NV; EqList E;

  for (j=0;j<_P->n;j++) drop_point[j]=0;  /* just to silence the compiler */
  IP_Check(_P,&V,&E);
  for(i=0;i<V.nv;i++){
    
    /* Dropping the vertex V.v[i]: */
    _P->np--;
    for (j=0;j<_P->n;j++){
      drop_point[j]=_P->x[V.v[i]][j];
      _P->x[V.v[i]][j]=_P->x[_P->np][j]; }
    if (IP_Check(_P,&NV,&E)) return 0;
    
    /* Reconstruct _P  */
    for (j=0;j<_P->n;j++){
      _P->x[_P->np][j]=_P->x[V.v[i]][j];
      _P->x[V.v[i]][j]=drop_point[j];}
    _P->np++; }
  return 1;
}
    
void DPircheck(CWS *_W, PolyPointList *_P){
  int i, j, k; 
  EqList E,DE;
  VertexNumList V,DV;
  EqList *B=&_W->B; 
  PolyPointList *_RDP= (PolyPointList *) malloc(sizeof (PolyPointList));
  PolyPointList *_PD = (PolyPointList *) malloc(sizeof (PolyPointList));
  if(_RDP==NULL) {puts("Unable to allocate space for _P"); exit(0);}
  if(_PD==NULL) {printf("Unable to allocate _PD\n"); exit(0);}
  
  IP_Check(_P,&V,&E);
  Make_Dual_Poly(_P,&V,&E,_PD);
  _RDP->n=_P->n;
  _RDP->np=B->ne;
  for (i=0;i<B->ne;i++) for (j=0;j<_P->n;j++) _RDP->x[i][j]=B->e[i].a[j];
  IP_Check(_RDP,&DV,&DE);
  _RDP->np=0;
  for (i=0;i<_PD->np;i++){
    int Good=1;
    for (k=0;k<DE.ne;k++){
      int Pairing=DE.e[k].c;
      for (j=0;j<_P->n;j++) Pairing+=DE.e[k].a[j]*_PD->x[i][j];
      if (Pairing<0) Good=0;}
    if(Good){
      for (j=0;j<_P->n;j++) _RDP->x[_RDP->np][j]=_PD->x[i][j];
      _RDP->np++;}}
  /*       Print_PPL(_RDP);   */
  if (irred(_RDP)){
    /* Print_PPL(_RDP); */
    for(i=0;i<_W->nw;i++){
      fprintf(outFILE,"%d ",(int) _W->d[i]);
      for(j=0;j<_W->N;j++) fprintf(outFILE,"%d ",(int) _W->W[i][j]);
      if(i+1<_W->nw) fprintf(outFILE," ");     }
    fprintf(outFILE,"\n"); }
  free(_PD); free(_RDP); 
} 


int virred(PolyPointList *_P, EqList *B)
{
  int i,j,k,drop_point[POLY_Dmax],equal;
  VertexNumList V; EqList E;
  assert(Ref_Check(_P,&V,&E));
  for (j=0;j<_P->n;j++) drop_point[j]=0;  /* just to silence the compiler */
  /* for(i=0;(i<V.nv);i++) {
    for (j=0;(j<_P->n);j++) fprintf(outFILE,"%d ",(int) _P->x[V.v[i]][j]);
    fprintf(outFILE,"\n");}
  fprintf(outFILE,"_CL->B:\n");
  for(k=0;k<B->ne;k++){
    for (j=0;(j<_P->n);j++) fprintf(outFILE,"%d ",(int) _CL->B[k][j]);
    fprintf(outFILE,"\n");} */
  for(k=0;k<B->ne;k++){
    equal=0;
    for(i=0;(i<V.nv)&&(!equal);i++){
      equal=1;
      for (j=0;(j<_P->n)&&equal;j++) 
	if (_P->x[V.v[i]][j]-B->e[k].a[j]) equal=0; }
    /* if (equal) break;}*/
    i--;
    if (!equal){ 
      fprintf(outFILE,"Vertex not found!\n"); 
      for(i=0;(i<V.nv);i++) {
	for (j=0;(j<_P->n);j++) fprintf(outFILE,"%d ",(int) _P->x[V.v[i]][j]);
	fprintf(outFILE,"\n");}
      for (j=0;j<_P->n;j++) fprintf(outFILE,"%d ",(int) _P->x[V.v[i]][j]);
      fprintf(outFILE,"\n");
      for (j=0;j<_P->n;j++) fprintf(outFILE,"%d ",(int) B->e[k].a[j]);
      fprintf(outFILE,"\n");
      return 1;}
    /* Dropping the vertex V.v[i]: */
    _P->np--;
    for (j=0;j<_P->n;j++){
      drop_point[j]=_P->x[V.v[i]][j];
      _P->x[V.v[i]][j]=_P->x[_P->np][j]; }
    if (IP_Check(_P,&V,&E)) return 0;
    
    /* Reconstruct _P  */
    for (j=0;j<_P->n;j++){
      _P->x[_P->np][j]=_P->x[V.v[i]][j];
      _P->x[V.v[i]][j]=drop_point[j];}
    _P->np++; 
    assert(Ref_Check(_P,&V,&E));}
  return 1;
}
    
void DPvircheck(CWS *_W, PolyPointList *_P) {
  int i, j; 
  EqList E; EqList *B=&_W->B;
  VertexNumList V;
  PolyPointList *_PD = (PolyPointList *) malloc(sizeof (PolyPointList));
  if(_PD==NULL) {printf("Unable to allocate _PD\n"); exit(0);}

     Ref_Check(_P,&V,&E);
     Make_Dual_Poly(_P,&V,&E,_PD);
     if (virred(_PD,B)){
       /* Print_PPL(_PD); */
       for(i=0;i<_W->nw;i++){
	 fprintf(outFILE,"%d ",(int) _W->d[i]);
	 for(j=0;j<_W->N;j++) fprintf(outFILE,"%d ",(int) _W->W[i][j]);
	 if(i+1<_W->nw) fprintf(outFILE," ");     }
       fprintf(outFILE,"\n"); }
     free(_PD); 
}

int Find_Ref_Subpoly(PolyPointList *_P, EqList *_E, VertexNumList *_V,
		       KeepList *_KL, int *rd);

int Find_RSP_Drop_and_Keep(PolyPointList *_P, VertexNumList *_V, EqList *_E, 
		   KeepList *_KL, int drop_num, int *rd){
    
  int i, j, n_irrel=0, IP/*, test_IP*/;
  int *new2old =(int *) malloc(sizeof(int)*POINT_Nmax);
  Long drop_point[POLY_Dmax];
  VertexNumList red_V/*, test_V*/;
  EqList new_E/*, test_E*/;
  PolyPointList *red_P = (PolyPointList *) malloc(sizeof (PolyPointList));

  if(red_P==NULL) {puts("Unable to allocate space for red_P"); exit(0);}
  if(new2old==NULL) {puts("Unable to allocate space for new2old"); exit(0);}

  for (j=0;j<_P->n;j++) drop_point[j]=_P->x[_V->v[drop_num]][j];

  /* Create new_E: Same as *_E, but irrelevant facets first */

  new_E.ne=_E->ne;
  for (i=0;i<_E->ne;i++){ 
    if (Eval_Eq_on_V(&(_E->e[i]),drop_point,_P->n)){
      for (j=0;j<_P->n;j++) new_E.e[n_irrel].a[j]=_E->e[i].a[j];
      new_E.e[n_irrel++].c=_E->e[i].c;}
    else {
      for (j=0;j<_P->n;j++) new_E.e[_E->ne+n_irrel-i-1].a[j]=_E->e[i].a[j];
      new_E.e[_E->ne+n_irrel-i-1].c=_E->e[i].c;}}

  /* Create red_P: old vertices, points that are not on irrelevant facets */
  red_P->np=0;
  red_P->n=_P->n;
  for (i=0;i<_V->nv;i++) if (i!=drop_num){
    for (j=0;j<_P->n;j++) red_P->x[red_P->np][j]=_P->x[_V->v[i]][j];
    new2old[red_P->np++]=_V->v[i];}
  for (i=0;i<_P->np;i++) 
    if ((i!=_V->v[drop_num])&&Relevant(_P->x[i],&_P->n,&new_E,&n_irrel)){
      for (j=0;j<_P->n;j++) red_P->x[red_P->np][j]=_P->x[i][j];
      new2old[red_P->np++]=i;}

  IP=Aided_IP_Check(red_P,&red_V,&new_E,n_irrel,_V->nv-1);

  /*  test_IP=IP_Check(red_P,&test_V,&test_E);
  printf("%d %d %d\n", _P->np, red_P->np, *rd); fflush(0);
  printf("%d %d %d\n", test_V.nv, test_E.ne, test_IP);
  printf("%d %d \n", red_V.nv, new_E.ne); fflush(0);
  if (IP) if ((!test_IP)||(test_V.nv!=red_V.nv)||
      (test_E.ne!=new_E.ne)){
    int k;
    fprintf(outFILE,"red_P:\n");fflush(0);
    Print_PPL(red_P);
    fprintf(outFILE,"\n"); 
    fprintf(outFILE,"red_V: ");
    for (i=0;i<red_V.nv;i++) fprintf(outFILE,"%d ",red_V.v[i]);
    fprintf(outFILE,"\n"); 
    fprintf(outFILE,"new_E:\n");
    for (i=0;i<new_E.ne;i++){
      for (k=0;k<_P->n;k++) fprintf(outFILE,"%d ", (int) new_E.e[i].a[k]);
      fprintf(outFILE," %d\n", (int) new_E.e[i].c);}
    fprintf(outFILE,"\n");
    fprintf(outFILE,"n_irrel: %d\n",n_irrel);
    fprintf(outFILE,"new2old:");
    for (i=0;i<red_P->np;i++) fprintf(outFILE," %d", new2old[i]);
    fprintf(outFILE,"\n");
    fprintf(outFILE,"test_V: ");
    for (i=0;i<test_V.nv;i++) fprintf(outFILE,"%d ",test_V.v[i]);
    fprintf(outFILE,"\n"); 
    fprintf(outFILE,"test_E:\n");
    for (i=0;i<test_E.ne;i++){
      for (k=0;k<_P->n;k++) fprintf(outFILE,"%d ", (int) test_E.e[i].a[k]);
      fprintf(outFILE," %d\n", (int) test_E.e[i].c);}
    fprintf(outFILE,"\n");
    exit(0);} */


  free(red_P);

  if(IP){
    VertexNumList new_V;

    /* Create new_V from red_V */

    new_V.nv=red_V.nv;
    for (i=0;i<red_V.nv;i++) if((new_V.v[i]=new2old[red_V.v[i]])==_P->np-1) 
      new_V.v[i]=_V->v[drop_num];

    /* Drop the vertex _V->v[drop_num]: */
    _P->np--;
    for (j=0;j<_P->n;j++) _P->x[_V->v[drop_num]][j]=_P->x[_P->np][j]; 
    j=kept(_P->np,_KL);
    if (j>=0) _KL->k[j]=_V->v[drop_num];

    if (Find_Ref_Subpoly(_P, &new_E, &new_V, _KL, rd)) {
      free(new2old); 
      return 1;}
    
    /* Reconstruct _P  */
    if (j>=0) _KL->k[j]=_P->np;
    for (j=0;j<_P->n;j++){
      _P->x[_P->np][j]=_P->x[_V->v[drop_num]][j];
      _P->x[_V->v[drop_num]][j]=drop_point[j];}
    _P->np++;}
  
  /* Attach a keep label to the vertex _V->v[drop_num]: */
  _KL->k[_KL->nk]=_V->v[drop_num];
  _KL->nk++; 
  free(new2old);
  return 0;
}

  
int Start_Find_Ref_Subpoly(PolyPointList *_P){
  int i,j;
  VertexNumList V, new_V; 
  EqList E, new_E; 
  KeepList KL; 
  Long drop_point[POLY_Dmax];
  int rd=0;

  if(!IP_Check(_P,&V,&E)) {
    fprintf(outFILE,"IP_check negative!\n"); exit(0);}

  KL.nk=0;
  for (j=0;j<_P->n;j++) drop_point[j]=0;  /* just to silence the compiler */

  for(i=0;i<V.nv;i++) {
    
    /* Try dropping the vertex V.v[i]: */
    _P->np--;
    for (j=0;j<_P->n;j++){
      drop_point[j]=_P->x[V.v[i]][j];
      _P->x[V.v[i]][j]=_P->x[_P->np][j]; }
    
    if (!IP_Check(_P,&new_V,&new_E)) KL.k[KL.nk++]=V.v[i];
    
    /* Reconstruct _P  */

    for (j=0;j<_P->n;j++){
      _P->x[_P->np][j]=_P->x[V.v[i]][j];
      _P->x[V.v[i]][j]=drop_point[j];}
    _P->np++; }

  return Find_Ref_Subpoly(_P, &E, &V, &KL, &rd);
}


int Find_Ref_Subpoly(PolyPointList *_P, EqList *_E, VertexNumList *_V,
		       KeepList *_KL, int *rd){

  /* Creates all subpolyhedra of an IP-polyhedron given by _P, _E, _V.
     The basic structure is:
     Search for bad facet;
     If (reflexive) if (!Add_NF_to_List) return;
     Drop_and_Keep vertices (on bad facet / everywhere);
     If (bad facet exists) Reduce_Poly;  */

  int i, j, badfacet=-1, nbadvert=0, maxnkept=-1;
  KeepList new_KL=*_KL;

  /* Choose bad facet with largest number of kept vertices: */

  for (i=0;i<_E->ne;i++) if (_E->e[i].c>1){
    int nkept=0;
    for (j=0;j<_KL->nk;j++) 
      if (!Eval_Eq_on_V(&(_E->e[i]),_P->x[_KL->k[j]],_P->n)) nkept++;
    if (nkept>maxnkept){
      maxnkept=nkept; 
      badfacet=i; }}

  if  (*rd>0) if (badfacet==-1) return 1; 

  if (badfacet!=-1) {
    /* Find vertices on bad facet: */
    nbadvert=0;
    for (i=0;i<_V->nv;i++) 
      if (!Eval_Eq_on_V(&(_E->e[badfacet]),_P->x[_V->v[i]],_P->n)){
	j=_V->v[nbadvert];
	_V->v[nbadvert]=_V->v[i];
	_V->v[i]=j;
	nbadvert++; } }

  (*rd)++;

  for(i=0;i<((badfacet==-1)?_V->nv:nbadvert);i++) 
    if((kept(_V->v[i],&new_KL)<0))
      if (Find_RSP_Drop_and_Keep(_P, _V, _E, &new_KL, i, rd)) return 1;

  (*rd)--;
  return 0;
}
 
void Max_check(CWS *_W, PolyPointList *_P) {
  int i, j; 
  EqList E;
  VertexNumList V;
  IP_Check(_P,&V,&E);
  if (Poly_Max_check(_P,&V,&E)){
    for(i=0;i<_W->nw;i++){
      fprintf(outFILE,"%d ",(int) _W->d[i]);
      for(j=0;j<_W->N;j++) fprintf(outFILE,"%d ",(int) _W->W[i][j]);
      if(i+1<_W->nw) fprintf(outFILE," ");     }
    fprintf(outFILE,"\n"); }
}

int  Poly_Max_check(PolyPointList *_P, VertexNumList *_V, EqList *_E){
  int rm;
  PolyPointList *_PD = (PolyPointList *) malloc(sizeof(PolyPointList));
  assert(_PD!=NULL);
  Make_Dual_Poly(_P,_V,_E,_PD);
  rm=!Start_Find_Ref_Subpoly(_PD);
  free(_PD);
  return rm;
}

int  Poly_Min_check(PolyPointList *_P, VertexNumList *_V, EqList *_E){
  PairMat PM; 
  Make_VEPM(_P,_V,_E,PM); 
  Complete_Poly(PM, _E, _V->nv, _P);
  return !Start_Find_Ref_Subpoly(_P);
}

void Overall_check(CWS *_W, PolyPointList *_P) {
  int i, j, k, span, lpm=0, vm=0, r=0; 
  EqList E, DE;
  VertexNumList V, DV;
  EqList *B=&_W->B;
  PolyPointList *_RDP= (PolyPointList *) malloc(sizeof (PolyPointList));
  PolyPointList *_PD = (PolyPointList *) malloc(sizeof (PolyPointList));
  PolyPointList *_PD2 = (PolyPointList *) malloc(sizeof (PolyPointList));
  if(_RDP==NULL) {puts("Unable to allocate space for _P"); exit(0);}
  if(_PD==NULL) {printf("Unable to allocate _PD\n"); exit(0);}
  if(_PD2==NULL) {printf("Unable to allocate _PD\n"); exit(0);}

  if ((_P->n<5) ? !IP_Check(_P,&V,&E) : !Ref_Check(_P,&V,&E)) {
    free(_PD); free(_PD2); free(_RDP); return;}
  
  /* assert(Ref_Equations(&E)); */

  for(i=0;i<_W->nw;i++){
    fprintf(outFILE,"%d ",(int) _W->d[i]);
    for(j=0;j<_W->N;j++) fprintf(outFILE,"%d ",(int) _W->W[i][j]);
    if(i+1<_W->nw) fprintf(outFILE," ");     }
  fflush(0);

  span=Span_Check(&E, B, &_P->n);
  
  Make_Dual_Poly(_P,&V,&E,_PD);
  *_PD2=*_PD;

  _RDP->n=_P->n;
  _RDP->np=B->ne;
  for (i=0;i<B->ne;i++) for (j=0;j<_P->n;j++) _RDP->x[i][j]=B->e[i].a[j];
  IP_Check(_RDP,&DV,&DE);
  _RDP->np=0;
  for (i=0;i<_PD->np;i++){
    int Good=1;
    for (k=0;k<DE.ne;k++){
      int Pairing=DE.e[k].c;
      for (j=0;j<_P->n;j++) Pairing+=DE.e[k].a[j]*_PD->x[i][j];
      if (Pairing<0) Good=0;}
    if(Good){
      for (j=0;j<_P->n;j++) _RDP->x[_RDP->np][j]=_PD->x[i][j];
      _RDP->np++;}}
  if (irred(_RDP)) lpm=1;

  if (!Start_Find_Ref_Subpoly(_PD)) r=1;

  if (span) if (virred(_PD2,B)) vm=1;
  
  if ((!span&&vm)||(!lpm&&vm)||(r!=vm)) 
    fprintf(outFILE,"span:%d lpm:%d vm:%d r:%d\n", span, lpm, vm, r);
  else{
    if(r) fprintf(outFILE,"r");
    if(lpm) fprintf(outFILE,"l");
    if(span) fprintf(outFILE,"s");
    if(!r&&!lpm&&!span) fprintf(outFILE,"-");
    fprintf(outFILE,"\n");}
  fflush(0);
  free(_PD); free(_PD2); free(_RDP); 
}
