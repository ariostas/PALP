#include "Global.h"
#include "Nef.h"
#include "Rat.h"

/*   ===============	    Typedefs and Headers	===================  */
#define min(a,b) (((a)<(b)) ? (a):(b))
#define max(a,b) (((a)>(b)) ? (a):(b))

void Sort_PPL(PolyPointList *_P, VertexNumList *_V);
void part_nef(PolyPointList *, VertexNumList *, EqList *, PartList *,
	      int *, NEF_Flags *);

void IP_Fiber_Data(PolyPointList *, PolyPointList *,int nv,
		   Long G[VERT_Nmax][POLY_Dmax][POLY_Dmax],int fd[VERT_Nmax],
		   int *nf,int CD);

/*   ===============	End of Typedefs and Headers	===================  */

/*   ===============	Begin of DYNamical Complete     ===================  */

void Print_DYN_PPL( DYN_PPL *_MP, const char *comment){
  int i,j;
  if(_MP->np > 20){
    fprintf(outFILE,"%d %d  %s\n", (int) _MP->np,_MP->n, comment);
    for(i = 0; i < _MP->np; i++) {
      for(j = 0; j < _MP->n; j++) 
	fprintf(outFILE,"%d ",(int) _MP->L[i].x[j]); 
      fprintf(outFILE,"\n");
    }
  }
  else {
    fprintf(outFILE,"%d %d  %s\n",_MP->n, (int)_MP->np, comment);
    for(i = 0;i < _MP->n; i++) {
      for(j = 0; j < _MP->np; j++) 
	fprintf(outFILE," %4d",(int) _MP->L[j].x[i]); 
      fprintf(outFILE,"\n");
    }
  }
}

void DYNadd_for_completion(Long *yDen, Long Den, EqList *_E, DYN_PPL *_CP){
  int i,n=_CP->n;
  Long yold[POLY_Dmax];

  if(Den>1) for(i=0;i<n;i++) {
    if(yDen[i]%Den) return;
    yold[i]=yDen[i]/Den;}
  else for(i=0;i<n;i++) yold[i]=yDen[i];
  for (i=0;i<_E->ne;i++) 
    if (Eval_Eq_on_V(&(_E->e[i]), yold, n) < 0) return;
  if(!(_CP->np < _CP->NP_max)){
    _CP->NP_max += 1000000;
    if((_CP->L=(Vector *)realloc(_CP->L,_CP->NP_max * sizeof(Vector))) == NULL)
      Die("Unable to realloc space for _CP->L");
  }
  for(i=0;i<n;i++) _CP->L[_CP->np].x[i]=yold[i];
  _CP->np++;
}

void DYNComplete_Poly(Long VPM[][VERT_Nmax], EqList *_E, int nv, DYN_PPL *_CP){
  int i,j,k,l,InsPoint,rank=0,n=_CP->n;
  Long MaxDist[EQUA_Nmax], InvMat[POLY_Dmax][POLY_Dmax], Den=1;
  Long yDen[POLY_Dmax];
  int OrdFac[VERT_Nmax], 
    BasFac[POLY_Dmax], one[POLY_Dmax], position[POLY_Dmax];
  LRat ind[POLY_Dmax][POLY_Dmax], x[POLY_Dmax], y[POLY_Dmax], f, 
    PInvMat[POLY_Dmax][POLY_Dmax];

  _CP->np=0;

  /* Calculate maximal distances from facets of Delta^* (Vertices of Delta) */

  for (i=0;i<_E->ne;i++) {  
    MaxDist[i]=0;
    for (j=0;j<nv;j++) 
    if (MaxDist[i]<VPM[i][j]) MaxDist[i]=VPM[i][j];}

  /* Order facets of Delta^* (Vertices of Delta) w.r.t. MaxDist   */

  OrdFac[0]=0;
  for (i=1;i<_E->ne;i++){
    InsPoint=i; 
    while (InsPoint&&(MaxDist[i]<MaxDist[OrdFac[InsPoint-1]])) InsPoint--;
    for (j=i;j>InsPoint;j--) OrdFac[j]=OrdFac[j-1];
    OrdFac[InsPoint]=i; }

  /* Find first POLY_Dmax linearly independent facets + Inverse Matrix */

  for (i=0;i<n;i++) for (j=0;j<n;j++) PInvMat[i][j]=LrI(0);
  for (i=0;i<n;i++) PInvMat[i][i]=LrI(1);
  i=0;
  while (rank<n){
    for (j=0;j<n;j++) x[j]=LrI(_E->e[OrdFac[i]].a[j]);
    for (j=0;j<n;j++) y[j]=LrI(0);
    y[rank]=LrI(1);
    for (j=0;j<rank;j++) {
      f=x[one[j]];
      for (k=0;k<n;k++) {
        x[k]=LrD(x[k],LrP(f,ind[j][k])); 
        y[k]=LrD(y[k],LrP(f,PInvMat[j][k]));  } }
    one[rank]=-1;
    for (l=0;(l<n)&&(one[rank]==-1);l++) if (x[l].N) one[rank]=l;
    if(one[rank]>-1){
      for (k=0;k<n;k++) {
        ind[rank][k]=LrQ(x[k],x[one[rank]]);
        PInvMat[rank][k]=LrQ(y[k],x[one[rank]]); }
      for (j=0;j<rank;j++) {
        f=ind[j][one[rank]];
        for (k=0;k<n;k++)         {
          ind[j][k]=LrD(ind[j][k],LrP(ind[rank][k],f));   
          PInvMat[j][k]=LrD(PInvMat[j][k],LrP(PInvMat[rank][k],f));  }     }
      BasFac[rank]=OrdFac[i];
      rank++; }  
    i++; }
  for (i=0;i<n;i++) for (j=0;j<n;j++) 
    Den=(Den/LFgcd(Den,PInvMat[i][j].D))*PInvMat[i][j].D;
  for (i=0;i<n;i++) for (j=0;j<n;j++) 
    InvMat[one[i]][j]=(Den/PInvMat[i][j].D)*PInvMat[i][j].N;

  for (i=0;i<n;i++){
    for (j=0;j<n;j++) {
      long long s=0;
      for(k=0;k<n;k++) s+=((long long) (InvMat[k][i]))*
			 ((long long) (_E->e[BasFac[j]].a[k]));
      if (s!=Den*(i==j)) {
	puts("something wrong in Make_Dual_Poly");
	exit(0);}}} 

  /* Examine all integer points of parallelogram:                         */
  /* The basic structure of the algorithm is:
  for (k=0;k<n-1;k++) position[k]=-1;      / * sets k=n-1; important!      *
  position[n-1]=-2;  / * starting point just outside the parallelogram     *
  while(k>=0){
    position[k]++;
    DO AT position;
    for(k=n-1;((position[k]==MaxDist[BasFac[k]]-1)&&(k>=0));k--) 
       position[k]=-1;  }
         / * sets k to the highest value where pos.[k] wasn't the max value; 
            resets the following max values to min values                 */
  /* Quantities linear in position can be changed with every change of
     position (here: yDen)                                                */

  for(i=0;i<n;i++) yDen[i]=0;
  for (k=0;k<n-1;k++) {   /* sets k=n-1; important!   */
    position[k]=-_E->e[BasFac[k]].c;   
    for(i=0;i<n;i++) yDen[i]-=_E->e[BasFac[k]].c*InvMat[i][k]; }
  position[n-1]=-_E->e[BasFac[n-1]].c-1;
  for(i=0;i<n;i++) yDen[i]-=(_E->e[BasFac[k]].c+1)*InvMat[i][n-1];
  while(k>=0){
    position[k]++;
    for(i=0;i<n;i++) yDen[i]+=InvMat[i][k];
    DYNadd_for_completion(yDen, Den, _E, _CP);
    for(k=n-1;(k>=0);k--){
      if (position[k]!=MaxDist[BasFac[k]]-_E->e[BasFac[k]].c) break;
      position[k]=-_E->e[BasFac[k]].c;
      for (i=0;i<n;i++) yDen[i]-=MaxDist[BasFac[k]]*InvMat[i][k]; }}
}

void DYNMake_VEPM(DYN_PPL *_P, VertexNumList *_V, EqList *_E, 
	       Long PM[][VERT_Nmax]){
  int i, j;
  for (i=0;i<_E->ne;i++) for (j=0;j<_V->nv;j++) 
    PM[i][j]=Eval_Eq_on_V(&_E->e[i],_P->L[_V->v[j]].x,_P->n);
}

/*   ===============	End of DYNamical Complete	===================  */
/*   ===============	Begin of FIBRATIONS             ===================  */
void PRINT_APL(AmbiPointList *_AP, const char *comment){
  int i,j;

  fprintf(outFILE,"%d %d  %s\n", _AP->N, _AP->np, comment); 
  for(i = 0; i < _AP->N; i++){	
    for(j = 0; j < _AP->np; j++) 
      fprintf(outFILE,(_AP->np>20) ? " %2d" : " %4d", (int) _AP->x[j][i]);
    fprintf(outFILE,"\n");
  }
}

void PRINT_MATRIX(Long *M, int l, int c, int C){
  int i,j;
  for(j = 0; j < l; j++){
    for(i = 0; i < c; i++)
      fprintf(outFILE,(c>20) ? " %3d" : " %4d", (int) *(M+i+C*j));
    fprintf(outFILE,"\n"); 
  }
}

void PRINT_TMATRIX(Long *M, int l, int c, int C){
  int i,j;

  for(i = 0; i < c; i++){	
      for(j = 0; j < l; j++) 
	fprintf(outFILE,(l>20) ? " %3d" : " %4d", (int) *(M+i+C*j));
      fprintf(outFILE,"\n");
    }
}

void PRINT_PL(PolyPointList *_P, const char *comment){

  fprintf(outFILE,"%d %d  %s\n", _P->n, _P->np, comment); 
  PRINT_TMATRIX(&_P->x[0][0], _P->np, _P->n, POLY_Dmax);
}

void PRINT_GORE(PolyPointList *_P, int codim, int n, const char *comment){

  PolyPointList *_P_AUX; 
  VertexNumList *_V_AUX; 
  EqList *_E_AUX;
  int i, j, Z=0;
  
  _P_AUX = (PolyPointList *) malloc(sizeof(PolyPointList));
  if (_P_AUX == NULL) Die("Unable to alloc space for PolyPointLis _P_AUX");
  _V_AUX = (VertexNumList *) malloc(sizeof(VertexNumList));
  if (_V_AUX == NULL) Die("Unable to alloc space for VertexNumList _V_AUX");
  _E_AUX = (EqList *) malloc(sizeof(EqList));
  if (_E_AUX == NULL) Die("Unable to alloc space for EqList _E_AUX");

  _P_AUX->n = _P->n - codim + 1;

 /* troubles if IP==0 not in last position:
  _P_AUX->np = _P->np - codim + 1;
  for(i = 0; i < _P_AUX->np; i++)
    for(j = 0; j < _P_AUX->n; j++)
      _P_AUX->x[i][j] = _P->x[i][j + codim - 1];	** replaced by:	*/
  _P_AUX->np = _P->np - codim;
  for(i = 0; i < _P->np; i++){ int z=0; j=0;
    while((j < _P_AUX->n)&&(!z)) {if(_P->x[i][j + codim - 1]) z=1; j++;};
    if(z){ for(j = 0; j < _P_AUX->n; j++)
           _P_AUX->x[i-Z][j] = _P->x[i][j + codim - 1]; }
    else  Z++;}
  for(j = 0; j < _P_AUX->n; j++) _P_AUX->x[_P_AUX->np][j] = 0;
  _P_AUX->np++;
				/* Print_PPL(_P_AUX,"Ref_Check input"); */
  assert(Ref_Check(_P_AUX, _V_AUX, _E_AUX));
  /* Find_Equations(_P_AUX, _V_AUX, _E_AUX);  ...  redundant */
  if(n == 0){
    Sort_PPL(_P_AUX, _V_AUX);
    fprintf(outFILE,"%d %d %s (nv=%d)\n",_P_AUX->n,_P_AUX->np,comment,_V_AUX->nv);
    for(i = 0; i < _P_AUX->n; i++){
      for(j = 0; j < _P_AUX->np; j++)
	fprintf(outFILE,(_P_AUX->np>20) ? " %3d" : " %4d", (int) _P_AUX->x[j][i]);
      fprintf(outFILE,"\n");
    } /*PRINT_PL(_P_AUX, comment);*/
  }
  else{
    Sort_VL(_V_AUX);
    Sort_PPL(_P, _V_AUX);
    if(n == 1) 
      fprintf(outFILE,"%d %d %s (nv=%d)\n",_P->n,_P->np,comment,_V_AUX->nv );
    if(n == 2){
      int o;
      fprintf(outFILE,"%d %d %s (nv=%d)\n",_P->n+1,_P->np,comment,_V_AUX->nv );
      for(j = 0; j < _P->np; j++){
	o = 1; i = 0;
	while(o && (i < codim - 1)){
	  if(_P->x[j][i] == 1) 
	    o = 0;
	  i++;
	}
	fprintf(outFILE,(_P->np>20) ? " %3d" : " %4d", o);
      }
      fprintf(outFILE,"\n");  
    }
    for(i = 0; i < _P->n; i++){
      for(j = 0; j < _P->np; j++)
	fprintf(outFILE,(_P->np>20) ? " %3d" : " %4d", (int) _P->x[j][i]);
      fprintf(outFILE,"\n");
    }
  }
  free(_P_AUX);free(_V_AUX);free(_E_AUX);
}

Long G_x_P(Long *Gi, Long *V, int *d){     

  Long x=0; int j; 
 
  for(j=0;j<*d;j++) 
    x+=Gi[j]*V[j]; 
  return x;
}

void PRINT_Fibrations(VertexNumList *_V, PolyPointList *_P, Flags *_F
		      /* ,PartList *_PTL */  ){
 
  Long G[VERT_Nmax][POLY_Dmax][POLY_Dmax];
  int s[VERT_Nmax], CD=_F->f, n, c, i, j, fib, nf, nv, np, dim[VERT_Nmax];
  PolyPointList *_P_AUX; 
  VertexNumList *_V_AUX; 
  EqList *_E_AUX;
  char C[VERT_Nmax];

  _P_AUX = (PolyPointList *) malloc(sizeof(PolyPointList));
  if (_P_AUX == NULL) Die("Unable to alloc space for PolyPointLis _P_AUX");
  _V_AUX = (VertexNumList *) malloc(sizeof(VertexNumList));
  if (_V_AUX == NULL) Die("Unable to alloc space for VertexNumList _V_AUX");
  _E_AUX = (EqList *) malloc(sizeof(EqList));
  if (_E_AUX == NULL) Die("Unable to alloc space for EqList _E_AUX");

  if (_P->np >= VERT_Nmax) Die("Need _P->np < VERT_Nmax in PRINT_Fibrations");
  IP_Fiber_Data(_P, _P_AUX, _V->nv, G, dim, &nf, CD); 
  if (nf >= VERT_Nmax) Die("Need  nf < VERT_Nmax in PRINT_Fibrations");
 
  

  if(nf){
    for(j = 0; j < _P->np - 1; j++)
      fprintf(outFILE,(_P->np > 20) ? "----" :"-----");
    fprintf(outFILE," #fibrations=%d\n",nf);
  }
  for(n = 0; n < nf; n++){
    c = 0;
    for(i = 0; i < (_P->np - 1) ; i++){
      j = dim[n]; fib = 1;
      while((j < _P->n) && fib){
	if (G_x_P(G[n][j], _P->x[i], &_P->n))
	  fib = 0;
	j++;
      }
      if(fib){
	for(j = 0; j < dim[n]; j++)
	  _P_AUX->x[c][j] = G_x_P(G[n][j], _P->x[i], &_P->n);
	s[c++] = i;
      }
    }
    _P_AUX->np = c; _P_AUX->n = dim[n];
    assert(Ref_Check(_P_AUX, _V_AUX, _E_AUX));
    for(i = 0; i < (_P->np-1); i++) C[i]='_'; 
    for(i = 0; i < c; i++)
      C[s[i]] = 'p';
    for(i = 0; i < _V_AUX->nv; i++)
      C[s[_V_AUX->v[i]]] = 'v';
    for(i = 0; i < _P->np-1; i++) 
      fprintf(outFILE,"%s%c", (_P->np > 20) ? "   " : "    ",C[i]);
    nv = _V_AUX->nv; np = (_P_AUX->np + 1);
    EL_to_PPL(_E_AUX, _P_AUX, &dim[n]);
    assert(Ref_Check(_P_AUX, _V_AUX, _E_AUX));
    {
      Long X[VERT_Nmax][VERT_Nmax]; 
      Make_VEPM(_P_AUX, _V_AUX, _E_AUX, X);
      Complete_Poly(X, _E_AUX, _V_AUX->nv, _P_AUX);
      /*Dim_Fib_CI(dim, n, _PTL, C);*/
      fprintf(outFILE,"  cd=%d  m:%3d %2d n:%2d %d\n",(_P->n - dim[n]), 
	      _P_AUX->np, _V_AUX->nv, np, nv);
    }

  } 
  free(_P_AUX);free(_V_AUX);free(_E_AUX);
}

/*   ===============	End of FIBRATIONS	        ===================  */

void Die(char *comment){
  printf("\n%s\n",comment); exit(0);
}

void Time_Info(time_t *_Tstart, clock_t *_Cstart, const char *comment){

  fprintf(outFILE, "%s     %dsec  %dcpu\n", comment, (int)
	  /* CLOCKS_PER_SEC::10^6 */
	  difftime(time(NULL), *_Tstart), (int) ((((long long)
	  clock() - *_Cstart) / (long long) CLOCKS_PER_SEC))); fflush(0);
}

void Print_Nefinfo(PartList *_PTL, /* Flags *_F,*/ time_t *_Tstart, 
		   clock_t *_Cstart){

  int i, d = 0, p = 0;
  
  for(i = 0; i < _PTL->n; i++){
    if(_PTL->DirProduct[i]) d++;
  }
  for(i = 0; i < _PTL->n; i++){
    if(_PTL->Proj[i]) p++;
  }
  fprintf(outFILE, "np=%d d:%d p:%d %4dsec  %4dcpu\n", _PTL->n-d-p, d, p,  
	  /* CLOCKS_PER_SEC::10^6 */ (int) difftime(time(NULL), *_Tstart), (int)
	  ((((long long) clock() - *_Cstart) / (long long) CLOCKS_PER_SEC)));
  fflush(0);
}

int N_Part(PartList *_PTL){

  int i;  
  for(i = 0; i < _PTL->n; i++)
    if((!_PTL->DirProduct[i]) && (!_PTL->Proj[i]))
      return 1;
  return 0;
}

void Print_L(LInfo *_L, int p, int v)
{
  int codim, D, i, j, N;

  if(v && p) Die("only -Lp OR -Lv !");

  N = (v ? _L->nv : (_L->nv + 1));
  if(v){
    assert(FIB_POINT_Nmax >= N);
    fprintf(outFILE,"%d %d Vertices in N-lattice:\n",_L->d, N);
    for(j = 0; j < _L->d; j++){
      for(i = 0; i < N; i++)	
				fprintf(outFILE,(_L->nv > 20) ? " %3d" : " %4d", (int) _L->VM[i][j]);
      fprintf(outFILE,"\n");
    }
  }
  for(j = 0; j < N; j++)
    fprintf(outFILE,(_L->nv > 20) ? "----" :"-----");
  fprintf(outFILE,"\n");
  assert(FIB_Nmax >= _L->nw); assert(FIB_POINT_Nmax >= _L->nv);
  for(i = 0; i < _L->nw; i++){
    D = 0; codim = 0;
    for(j = 0; j < _L->nv; j++){
      D +=  _L->W[i][j];
      if(_L->W[i][j] == 0) 
				codim++;
    }
    for(j = 0; j < _L->nv; j++)
      fprintf(outFILE,(N > 20) ? " %3d" : " %4d",(int) _L->W[i][j]);    
    fprintf(outFILE,"  d=%d  ",D);
    fprintf(outFILE,"codim=%d\n",codim+_L->d-_L->nv+1);
  }
}

void Make_L(PolyPointList *_P, VertexNumList *_V, LInfo *_L , int p, int v){
  int i,j;
 
  if(v){
    for(i = 0; i < _V->nv; i++)
      for(j = 0; j < _P->n; j++)
				_L->VM[i][j] = _P->x[_V->v[i]][j];
    _L->nv = _V->nv;
  }
  if(p){
    assert(_P->np <= VERT_Nmax);
    for(i = 0; i < _P->np; i++)
      for(j = 0; j < _P->n; j++)
				_L->VM[i][j] = _P->x[i][j];
    _L->nv = (_P->np-1);
  }
  _L->Wmax = FIB_Nmax; _L->nw = 0; _L->d = _P->n;
  IP_Simplex_Decomp(_L->VM, _L->nv, _P->n, &_L->nw, _L->W, _L->Wmax,0);
}

int IntSqrt(int q)
{				/* sqrt(q) => r=1; r'=(q+r*r)/(2r); */
  if (q == 0)
    return 0;
  if (q < 4)
    return 1;
  else {			/* troubles: e.g. 9408 */
    long long r = (q + 1) / 2, n;
    while (r > (n = (q + r * r) / (2 * r)))
      r = n;
    if (q < r * r)
      r--;
    if ((r * r <= q) && (q < (r + 1) * (r + 1)))
      return (int) r;
    else{
      printf("Error in sqrt(%d)=%d\n", q, (int) n);
      exit(0);
    }
  }
  return 0;

}void INCI_TO(int I[], INCI * _X, int *_n)
{
/*  INCI X -> (0,0,1,0,.....,0,1,0);*/
  int i;

  INCI Y = *_X;
  for (i = 0; i < *_n; i++) {
    if (INCI_M2(Y)) 
      I[i] = 1;
    else 
      I[i] = 0;
    Y = INCI_D2(Y);
  }
}

int Num_Pos(Cone *_C){

  int d, n=2;

  for (d = 1; d < _C->dim; d++)
    n += _C->nface[d];
  return n;
}
     
void Make_PosetList(Cone *_C, Poset_Element_List *_PEL){

  int n = 0, d, i;

  for (d = 0; d <= _C->dim; d++)
    for (i = 0; i < _C->nface[d]; i++) {
      _PEL->L[n].dim = d;
      _PEL->L[n].num = i;
      n++;			
    }
}

int Interval_Check(int rank, Poset_Element * _x, Poset_Element * _y,
		   Cone * _C)
{
  int flag = 0;

  if (_x->dim == (_y->dim - rank))
    if (INCI_LE(_C->edge[_x->dim][_x->num], _C->edge[_y->dim][_y->num]))
      flag = 1;
  return flag;
}

void Make_Intervallist(Interval_List *_IL, Poset_Element_List *_PEL, Cone *_C){

  int d, i, j;

  _IL->n = 0;
  for (d = 0; d <= _C->dim; d++)
    for (i = 0; i < _PEL->n; i++)
      for (j = 0; j <= i; j++)
	if (Interval_Check(d, &_PEL->L[j], &_PEL->L[i], _C) == 1) {
	  _IL->L[_IL->n].min = j;
	  _IL->L[_IL->n].max = i;
	  _IL->n++;		
	}
}

int Make_Mirror(EPoly *_EP, int h[][POLY_Dmax], int D, int dim)
{
  int i, j, k, u, v, chi = 0, H;

  for (u = 0; u < 4*(Pos_Max); u++)
    for (v = 0; v < 3*(Pos_Max); v++){
      if(((-2 * D + u) > dim) || ((-2 * D + u) < 0) ||
	 ((-2 * D + v) > dim) || ((-2 * D + v) < 0)){
	if (_EP->E[u][v] != 0)
	  Die("Something wrong with E poly");
      }
      else {
	h[dim - (u - 2*D)][v - 2*D] = _EP->E[u][v];
	chi += _EP->E[u][v];
	if(((u - 2*D + v - 2*D) %2 ) != 0)
	  h[dim - (u - 2*D)][v - 2*D] *= -1;
      }
    }
  if((dim %2) != 0) chi = chi*(-1); 
  H = - dim*chi;
  for (i = 0; i <= dim; i++)
    for (j = 0; j <= dim; j++){
      if(((i+j) %2) == 0) 
	k = 1;
      else 
	k = -1;
      H += 3*k*(2*i - dim)*(2*i - dim)*h[i][j];
    }
  if (H != 0) Die("paper: CY 4-folds and toric fibrations; equation (8) don't hold");
  if (dim == 4){
    if((-h[2][2] + 44*h[0][0] + 4 * h[1][1] - 2 * h[1][2] + 4 * h[1][3] + 
	20 * h[0][2] - 52 * h[0][1]) != 0)
	Die("paper: CY 4-folds and toric fibrations; equation (9) don't hold");
  }
  return chi;
}

int Max_S(PartList * _PTL, int *_n)
{
  int i, k, m = 0, nv, Nv = 0;
  for (i = 0; i < _PTL->codim; i++) {
    nv = 0;
    for (k = 0; k < _PTL->nv; k++)
      if (_PTL->S[*_n][k] == i)
	nv++;
    if (Nv < nv) {
      Nv = nv;
      m = i;
    }
  }
  return m;
}

void PrintDegrees(	/*Flags *_F,*/ 	LInfo *_L, PartList *_PTL, int m, 
			/*int n,   */ 	int S[VERT_Nmax])
{
  int d[POLY_Dmax], i, j;

  for (i = 0; i < _L->nw; i++){
    fprintf(outFILE," (");
    for (j = 0; j < _PTL->codim; j++)
      d[j] = 0;
    for (j = 0; j < _L->nv; j++)
      d[S[j]] += _L->W[i][j];
    for (j = 0; j < _PTL->codim; j++)
      if(j != m)
	fprintf(outFILE,"%d ",d[j]);
    fprintf(outFILE,"%d",d[m]);
    fprintf(outFILE,")");
  }
}

void PrintWeights(CWS * _W)
{
  int i,j;
  for (i = 0; i < _W->nw; i++) {
    fprintf(outFILE, "%ld", (long) _W->d[i]);
    for (j = 0; j < _W->N; j++)
      fprintf(outFILE, " %ld", (long) _W->W[i][j]);
    if (i != (_W->nw - 1))
      fprintf(outFILE, "  ");
  }
}

void PrintDiamond(int h[][POLY_Dmax], int dim)
{
  int i,j;

  fprintf(outFILE, "\n\n"); fflush(0);
  for (i = 0; i <= dim; i++) {
    fprintf(outFILE, "        ");
    for (j = 0; j <= (dim - i); j++)
      fprintf(outFILE, "     ");
    for (j = 0; j <= i; j++)
      fprintf(outFILE, "   h%2d%2d   ", i - j, j);
    fprintf(outFILE, "\n\n");
    fflush(0);
  }
  for (i = 1; i <= dim; i++) {
    fprintf(outFILE, "        ");
    for (j = 0; j <= i; j++)
      fprintf(outFILE, "     ");
    for (j = i; j <= dim; j++)
      fprintf(outFILE, "   h%2d%2d   ", dim - j + i, j);
    fprintf(outFILE, "\n\n");
    fflush(0);
  }
  fprintf(outFILE, "\n\n");fflush(0);
  for (i = 0; i <= dim; i++) {
    fprintf(outFILE, "     ");
    for (j = 0; j <= (dim - i); j++)
      fprintf(outFILE, "     ");
    for (j = 0; j <= i; j++)
      fprintf(outFILE, "%10d", h[i - j][j]);
    fprintf(outFILE, "\n\n"); 
  }
  for (i = 1; i <= dim; i++) {
    fprintf(outFILE, "     ");
    for (j = 0; j <= i; j++)
      fprintf(outFILE, "     ");
    for (j = i; j <= dim; j++)
      fprintf(outFILE, "%10d", h[dim - j + i][j]);
    fprintf(outFILE, "\n\n"); 
  }
  fflush(0);
}

void Print_Points(PolyPointList *_P, int c, int nv, int S[VERT_Nmax]){

  int i, P=0;
  
  for(i = nv; i < (_P->np - 1); i++)
    if(S[i] == c){
      fprintf(outFILE,"%d ",i);
      P=1;
    }
  if(P)
    fprintf(outFILE," ");
}

void Output(PolyPointList * _P, /* PolyPointList * _DP,*/
	    PolyPointList * _P_D, CWS * _W, EPoly *_EP, /* EqList * _E,*/
	    VertexNumList * _V, int *_n, PartList *_PTL, 
	    int *_codim, FILE *outFILE, Flags * _F, int *_D, LInfo *_L)
{
  int i, j, k, m, chi = 0, D = (_P_D->n + 1), dim = (_P->n - *_codim);
  int h[POLY_Dmax][POLY_Dmax] = {{0}, {0}}, S[VERT_Nmax];
  
  
  m = Max_S(_PTL, _n);
  if (((!_PTL->DirProduct[*_n]) || _F->Dir) && ((!_PTL->Proj[*_n]) || _F->Proj)) {
    if (_F->H == 0){ 
      if (_F->w){ 
	PrintWeights(_W); fprintf(outFILE, " ");
      }
#ifdef	WRITE_CWS
      if (_F->Msum == 1)
	fprintf(outFILE, " d=%d %d", (int) _D[0], (int) _D[1]);
#endif
      if (!_F->p){
	chi =  Make_Mirror(_EP, h, D, dim);
	fprintf(outFILE, "H:");
	for (i = 1; i < dim; i++)
	  fprintf(outFILE, "%d ", h[1][i]);
	fprintf(outFILE, "[%d]", chi);
	for (i = 1; i <=  dim/2; i++)
	  if (h[0][i] != 0)
	    fprintf(outFILE, " h%d=%d", i, h[0][i]);
	if(h[0][0] != 1){
	  if (_PTL->DirProduct[*_n])
	    fprintf(outFILE, " h%d=%d", 0, h[0][0]);
	  else
	    Die("\nh00 not 1 !!!\n");
	}
      }
      if((_P_D->np - *_codim) <= VERT_Nmax)
	for(i = 0; i < (_P_D->np - *_codim); i++){
	  S[i] = 0;
	  for(j = 0; j < (*_codim - 1); j++)
	    if(_P_D->x[i][j])
	      S[i] = (j+1);
	}
      if (*_codim == 2) {
	fprintf(outFILE, " P:%d V:", *_n);
	i = 0;
	if (m == 0)
	  i = 1;
	for (j = 0; j < _PTL->nv; j++)
	  if (_PTL->S[*_n][j] == i)
	    fprintf(outFILE, "%d ", j);
	fprintf(outFILE," ");
	if((_P_D->np - *_codim) <= VERT_Nmax)
	  Print_Points(_P, i, _V->nv, S);
	else
	  fprintf(outFILE, " _P->np > VERT_Nmax! ");
      } 
      else {
	fprintf(outFILE, " P:%d ", *_n);
	j = 0;
	for (i = 0; i < *_codim; i++)
	  if (i != m) {
	    fprintf(outFILE, "V%d:", j);
	    for (k = 0; k < _PTL->nv; k++)
	      if (_PTL->S[*_n][k] == i)
		fprintf(outFILE, "%d ", k);
	    fprintf(outFILE," ");
	    j++;
	    if((_P_D->np - *_codim) <= VERT_Nmax)
	      Print_Points(_P, i, _V->nv, S);
	    else
	      fprintf(outFILE, " _P->np > VERT_Nmax! ");
	  }
      }
      if(_PTL->DProj[*_n]) 
	fprintf(outFILE, "DP ");
      if(_PTL->DirProduct[*_n])
	fprintf(outFILE, " D");
      if(_F->Lv || (_F->Lp && (_P_D->np - *_codim) <= VERT_Nmax)) 
	PrintDegrees(/*_F,*/ _L, _PTL, m, /* *_n,*/ S);
      fflush(0);
    }
    else{
      chi =  Make_Mirror(_EP, h, D, dim);
      PrintDiamond(h, dim);
    }
  }
}


int Min_Dim(int n, int T_flag)
{
  if (!T_flag){
    int i;
    i = n / 2;
    if((n % 2) != 0)
      i += 1;
    return i;
  } else
    return n;
}

void Init_ST(SPoly *_S, SPoly *_T, Poset_Element_List *_PEL)
{
  int i, d;

  for(i = 0; i < _PEL->n; i++){
    for(d = 1; d <= _PEL->L[i].dim; d++){
      _S[i].S[d] = 0;
      _T[i].S[d] = 0;
    }
    _S[i].S[0] = 1; _T[i].S[0] = 0;
  }
}

Long Eval_Eq_on_x(Long *_x, Equation *_E, int dim){

  Long c = _E->c; 
  int d;

  for(d = 0; d < dim; d++)
    c += _x[d] * _E->a[d];
  return c;
}

INCI INCI_to_x(int n, DYN_PPL *_P, EqList *_E){

  INCI X = INCI_0();
  int i;

  for(i = 0; i < _E->ne; i++)
    X = INCI_PN(X,Eval_Eq_on_x(_P->L[n].x, &_E->e[i], _P->n));
  return X;
}

void Poly_To_ST(DYN_PPL *_P, EqList *_E, Cone *_C, SPoly *_S, SPoly *_T,
		Poset_Element_List *_PEL, int l, int T_flag)

{
  int min, i, dim_flag, in_flag, n, dim, num;
  INCI X;
  
  min = Min_Dim(l, T_flag);
  for(i = 0; i < _P->np; i++){
    dim_flag = 1; in_flag = 0;
    n = (_PEL->n - 1); 
    X = INCI_to_x(i, _P, _E);
    while(dim_flag && !in_flag){
      dim = (_C->dim - _PEL->L[n].dim);
      num = (_C->nface[dim] - _PEL->L[n].num -1);
      if (INCI_LE(_C->edge[dim][num], X)){
	_S[n].S[l] += 1;
	if (INCI_LE(X, _C->edge[dim][num])){
	  _T[n].S[l] += 1;
	  in_flag = 1;
	}	      
      }
      n--;
      if (_PEL->L[n].dim < min)
	dim_flag = 0;
    }
  }
}

void New_CPVE(PolyPointList *_P, DYN_PPL *_CP, VertexNumList *_V, 
	    VertexNumList *_CV, EqList *_E, EqList *_CE, int n)
{
  int j, d, nv = 0;
  
  for(j = 0; j < _E->ne; j++){
    _CE->e[j].c = n*(_E->e[j].c);
    for(d = 0; d < _P->n; d++)
      _CE->e[j].a[d] = _E->e[j].a[d];
  }
  for(j = 0; j < _V->nv; j++){
    for(d = 0; d < _P->n; d++)
      _CP->L[nv].x[d] = n*(_P->x[_V->v[j]][d]);
    _CV->v[nv] = nv;
    nv++;
  }
  _CV->nv = _V->nv; _CP->np = _V->nv; _CP->n = _P->n; _CE->ne = _E->ne;
}

void Poly_To_DYNPoly(DYN_PPL *_CP, PolyPointList *_P){
  int i, d;

  assert(_P->np <= _CP->NP_max);
  for(i = 0; i < _P->np; i++)
    for(d = 0; d < _P->n; d++)
      _CP->L[i].x[d] = _P->x[i][d];
  _CP->n = _P->n; _CP->np = _P->np;
}

void Make_S_Poly(Cone *_C, VertexNumList *_V, EqList *_E, PolyPointList *_P,
		 Poset_Element_List *_PEL, SPoly *_S, int SINFO, int CHECK_SERRE)
{
  Long PM[EQUA_Nmax][VERT_Nmax];
  DYN_PPL CP;
  VertexNumList *_CV;
  EqList *_CE;
  SPoly *_T;
  int i=1, j, d, min;

  CP.NP_max = NP_Max;
  min = Min_Dim(_C->dim,CHECK_SERRE);

  CP.L = (Vector *) calloc(CP.NP_max, sizeof(Vector));
  if (CP.L == NULL) Die("Unable to alloc space for PolyPointList _CP.L");
  _CV = (VertexNumList *) malloc(sizeof(VertexNumList));
  if (_CV == NULL) Die("Unable to alloc space for VertexNumList _CV");
  _CE = (EqList *) malloc(sizeof(EqList));
  if (_CE == NULL) Die("Unable to alloc space for EqList _CE");
  _T = (SPoly *) calloc(_PEL->n, sizeof(SPoly ));
  if (_T == NULL) Die("Unable to alloc space for SPoly _T");

  Poly_To_DYNPoly(&CP, _P);
  Init_ST(_S, _T, _PEL);
  Poly_To_ST(&CP, _E, _C, _S, _T, _PEL, i, CHECK_SERRE);
  
  if(SINFO){
    printf("\n\n#points in largest cone:\n");
    printf("layer: %2d #p: %8d #ip: %8d\n", i=1, (int) _S[_PEL->n-1].S[1],
	   (int)_T[_PEL->n-1].S[1]);
  }
  for (i = 2; i <= min; i++) {
    New_CPVE(_P, &CP, _V, _CV, _E, _CE, i);
    DYNMake_VEPM(&CP,_CV,_CE, PM);
    DYNComplete_Poly(PM, _CE, _V->nv, &CP);
    Poly_To_ST(&CP, _CE, _C, _S, _T, _PEL, i, CHECK_SERRE);
    if(SINFO)
      printf("layer: %2d #p: %8d #ip: %8d\n",i, (int)_S[_PEL->n-1].S[i],
						 (int) _T[_PEL->n-1].S[i]);
  }
  for(i = 0; i < _PEL->n; i++){
    for (j = 0; j < _PEL->L[i].dim; j++){
      min = Min_Dim(_PEL->L[i].dim, CHECK_SERRE);
      for (d = min; d > 0; d--) {
				_S[i].S[d] += -_S[i].S[d - 1];
				_T[i].S[d] += -_T[i].S[d - 1];
      }
    }
    if(CHECK_SERRE){
      for (d = 1; d < _PEL->L[i].dim; d++)
				if(_S[i].S[d] != _T[i].S[_PEL->L[i].dim - d]){
					puts("No Serre Duality!");
					for(d = 1; d < _PEL->L[i].dim; d++)
						printf("S[%d]: %3d   T[%d]: %3d\n",d, (int) _S[i].S[d],
									 (_PEL->L[i].dim - d), (int) _T[i].S[_PEL->L[i].dim - d]);
					exit(0);
				}            
    }
    if (!CHECK_SERRE){
      for (d = 1; d < _PEL->L[i].dim / 2; d++)
				_S[i].S[_PEL->L[i].dim - d] = _T[i].S[d];
      if (_PEL->L[i].dim != 0)
				_S[i].S[_PEL->L[i].dim] = 0;
    }
  }
  free(_T);free(_CV);free(_CE);free(CP.L);
}

void SB_To_E(EPoly *_EP, Cone *_C, Poset_Element_List *_PEL, BPoly * _BL, 
		Interval_List * _IL, SPoly *_S_D, SPoly *_S_N, int *_codim){
  /* E[0][0] = E[2*_C->dim][2*_C->dim] */ 

  int i, j, x, y, u, v, s, dim, dmax, dmin;

  for (u = 0; u < 4*(Pos_Max); u++)
    for (v = 0; v < 4*(Pos_Max); v++)
      _EP->E[u][v] = 0;
	
  for (i = 0; i < _IL->n; i++) {
    dmin = _PEL->L[_IL->L[i].min].dim;
    dmax = _PEL->L[_IL->L[i].max].dim;
    dim = dmax - dmin;
    s = 1;
    for (j = 0; j < dmin; j++)
      s *= -1;
    for (y = 0; y < (_C->dim - dmax + 1); y++)
      for (x = 0; x < (dmin + 1); x++)
				for (u = 0; u < (dim + 1); u++)
					for (v = 0;v <= (dim / 2); v++)
						_EP->E[2*_C->dim - x + y + dmax - u - *_codim]
							[2*_C->dim + y + x + v - *_codim] += 
							_S_D[_IL->L[i].min].S[x] * _S_N[_PEL->n - _IL->L[i].max - 1].S[y]
							* _BL[i].B[u][v] * s;
  }
}

void Make_Cone(Cone * _C_D, Cone * _C_N, FaceInfo * _I, PolyPointList * _P)
{
  int d, i;
  
  _C_D->dim = (_P->n + 1); _C_N->dim = (_P->n + 1);
  _C_D->nface[0] = 1; _C_N->nface[0] = 1;
  _C_D->nface[_P->n + 1] = 1; _C_N->nface[_P->n + 1] = 1;
  _C_D->edge[0][0] = INCI_0(); _C_N->edge[0][0] = INCI_0();
  _C_D->edge[_C_D->dim][0] = INCI_0();
  for (i = 0; i < _I->nf[0]; i++)
    _C_D->edge[_C_D->dim][0] = 
      INCI_OR(_C_D->edge[_C_D->dim][0], _I->v[0][i]);
  _C_N->edge[_C_N->dim][0] = INCI_0();
  for (i = 0; i < _I->nf[_P->n - 1]; i++)
    _C_N->edge[_C_N->dim][0] = 
      INCI_OR(_C_N->edge[_C_N->dim][0], _I->f[_P->n - 1][i]);
  for (d = 0; d < _P->n; d++) {
    _C_D->nface[d + 1] = _I->nf[d];
    _C_N->nface[_C_D->dim - d -1] = _I->nf[d];
    for (i = 0; i < _I->nf[d]; i++) {
      _C_D->edge[d + 1][i] = _I->v[d][i];
      _C_N->edge[_C_D->dim - d -1][_I->nf[d] - i -1] = _I->f[d][i];
   }
  }
}

void Make_EN(PolyPointList * _P, VertexNumList * _V, EqList * _EN, int *_codim)
{
  int x, i, j;

  _EN->ne = _V->nv;
  for (i = 0; i < _V->nv; i++) {
    x = 1; j = 0;
    while((x == 1) && (j < (*_codim - 1))){
      if(_P->x[_V->v[i]][j] == 1)
	x = 0; 
      j++;
    }
    _EN->e[i].c = x;
    for (j = 0; j < _P->n; j++){
      if(j < *_codim - 1)
	_EN->e[i].a[j] = _P->x[_V->v[i]][j] - x;
      else
	_EN->e[i].a[j] = _P->x[_V->v[i]][j];  
    }
  }
}

int Remove_Proj(PolyPointList * _P, /* int *_n,*/ int *_codim)
{
  int nv, Nv=0, j, i = 0, proj_flag = 0;

  while (!proj_flag && (i < *_codim - 1)) {
    nv=0;
    for(j = 0; j < _P->np; j++)
      if(_P->x[j][i] == 1)
	nv++;
    Nv += nv;
    if(nv == 2)
      proj_flag = 1;
    i++;
  }
  if((_P->np - Nv) == 2)
    proj_flag = 1;
  return proj_flag;
}

void Make_Gore_Poly(PolyPointList * _P, PolyPointList * _DP,
		    PolyPointList * _P_D, PolyPointList * _P_N,
		    VertexNumList * _V, PartList * _PTL, int *_codim,
		    int *_n)
{
  /*_P_N from _DP (in M-lattice), _P_D from _P (in N-lattice)  */

  int i, l, k, d, sum, ip, c;

  if(*_codim <= 0) 
    Die("Need Codim > 0");
  _P_N->n = _DP->n + *_codim - 1;
  _P_N->np = 0;
  for (l = 0; l < _DP->np; l++) {
    for (i = 0; i < *_codim; i++) {
      ip = 1;
      k = 0;
      while (ip && (k < _PTL->nv)) {
	sum = 0;
	for (d = 0; d < _DP->n; d++)
	  sum += _DP->x[l][d] * _P->x[_V->v[k]][d];
	if (((_PTL->S[*_n][k] == i) && (sum < -1)) ||
	    ((_PTL->S[*_n][k] != i) && (sum < 0)))
	  ip = 0;
	k++;
      }
      if (ip == 1) {
	assert(_P_N->np < POINT_Nmax);
	for (d = 0; d < _DP->n; d++)
	  _P_N->x[_P_N->np][d + *_codim - 1] = _DP->x[l][d];
	for (d = 1; d < *_codim; d++) {
	  if (d == i)
	    _P_N->x[_P_N->np][d - 1] = 1;
	  else
	    _P_N->x[_P_N->np][d - 1] = 0;
	}
	_P_N->np++;
      }
    }
  }
  _P_D->n = _P->n + *_codim - 1;
  _P_D->np = 0;
  for (l = 0; l < _P->np; l++) {
    for (i = 0; i < *_codim; i++) {
      ip = 1;
      k = 0;
      while (ip && (k < _P_N->np)) {
	sum = 0;
	for (d = 0; d < _P->n; d++)
	  sum += _P->x[l][d] * _P_N->x[k][d + *_codim - 1];
	if (i > 0) {
	  if (((_P_N->x[k][i - 1] == 1) && (sum < -1)) ||
	      ((_P_N->x[k][i - 1] == 0) && (sum < 0)))
	    ip = 0;
	} else {
	  c = -1;
	  for (d = 0; d < *_codim - 1; d++)
	    if (_P_N->x[k][d] == 1)
	      c = 0;
	  if (sum < c)
	    ip = 0;
	}
	k++;
      }
      if (ip == 1) {
	assert(_P_D->np < POINT_Nmax);
	for (d = 0; d < _P->n; d++)
	  _P_D->x[_P_D->np][d + *_codim - 1] = _P->x[l][d];
	for (d = 1; d < *_codim; d++) {
	  if (d == i)
	    _P_D->x[_P_D->np][d - 1] = 1;
	  else
	    _P_D->x[_P_D->np][d - 1] = 0;
	}
	_P_D->np++;
      }
    }
  }
}

Poset Int_Pos(int i, Interval_List * _IL, Poset_Element_List * _PEL)
{
  Poset P;

  P.x = _PEL->L[_IL->L[i].min];
  P.y = _PEL->L[_IL->L[i].max];
  return P;
}

void M_To_B(BPoly *_BP,  BPoly *_MP, int d, int rho){

  int M[Pos_Max][Pos_Max] = {{0}, {0}}, i, u, v;

  if(rho == 0)
    M[0][0] = 1;
  else
    for(u = 0; u <= rho; u++)
      for(v = 0; v < rho/2 + (rho % 2); v++)
	M[rho - u][rho - v] = _MP->B[u][v]; 
  for(i = 1; i <= d - rho; i++){
    for(u = rho + i; u > 0; u--)   /* (u degree != 0) && (v degree != 0) */
      for(v = rho + i; v > 0; v--)
	M[u][v] = - M[u-1][v] + M[u][v-1];
    for(v = rho + i; v > 0; v--)   /* (u degree == 0) && (v degree != 0) */ 
      M[0][v] =  M[0][v-1];
    for(u = rho + i; u > 0; u--)   /* (u degree != 0) && (v degree == 0) */ 
      M[u][0] =  - M[u-1][0];
    M[0][0] = 0;                   /* (u degree == 0) && (v degree == 0) */ 
  }
  for(u = 0; u <= d; u++)
    for(v = 0; v < d/2 + (d % 2); v++)
      _BP->B[u][v] += M[u][v];
}

void N_To_B(BPoly *_BP,  BPoly *_NP, int d, int rho){

  int N[Pos_Max][Pos_Max] = {{0}, {0}}, i, u, v;

  if(rho == 0)
    N[0][0] = 1;
  else
    for(u = 0; u <= rho; u++)
      for(v = 0; v < rho/2 + (rho % 2); v++)
	N[u][v] = _NP->B[u][v];
  for(i = 1; i <= d - rho; i++){
    for(u = rho + i; u > 0; u--)      /* (u degree != 0) && (v degree != 0) */
      for(v = rho / 2 + i; v > 0; v--)
	N[u][v] =  N[u-1][v-1] - N[u][v];
    for(v = rho/2 + i; v > 0; v--)    /* (u degree == 0) && (v degree != 0) */
      N[0][v] *= -1;
    for(u = rho + i; u >= 0; u--)     /* (u degree >= 0) && (v degree == 0) */
      N[u][0] *= -1;                       
  }
  for(u = 0; u <= d; u++)
    for(v = 0; v < d/2 + (d % 2); v++)
      _BP->B[u][v] -= N[u][v];  
}

void Make_B_Poly(Cone * _C, Poset_Element_List * _PEL, Interval_List * _IL,
		 BPoly *_BL)
{
  int D, i, j, k, d, n, n_;
  Poset Pos, pos;

  for (i = 0; i < _IL->n; i++) {
    Pos = Int_Pos(i, _IL, _PEL);
    d  = (Pos.y.dim - Pos.x.dim);
    for (j = 0; j <= d; j++)
      for (k = 0; k <= d / 2; k++)
      	  _BL[i].B[j][k] = 0;
    if (d == 0)
      _BL[i].B[0][0] = 1;
  }
  for (D = 1; D <= _C->dim; D++)   
    for (n = 0; n < _IL->n; n++) { 
      Pos = Int_Pos(n, _IL, _PEL);
      if ((Pos.y.dim - Pos.x.dim) == D) {
	for (d = 0; d < D; d++) 
	  for (n_ = 0; n_ < n; n_++) {
	    pos = Int_Pos(n_, _IL, _PEL);
	    if (Interval_Check(0, &pos.x, &Pos.x, _C))
	      if (Interval_Check((D - d), &pos.y, &Pos.y, _C))
		M_To_B(&_BL[n], &_BL[n_], D, d);	
	    if (Interval_Check(0, &pos.y, &Pos.y, _C))
	      if (Interval_Check((D - d), &Pos.x, &pos.x, _C)) 
		N_To_B(&_BL[n], &_BL[n_], D, d);
	  }
      }	
    }	
}
  
void Compute_E_Poly(EPoly *_EP,  
		    PolyPointList * _P_D, VertexNumList * _V_D, EqList * _E_D,
		    PolyPointList * _P_N, VertexNumList * _V_N, EqList * _E_N,
		    int *_codim, Flags * _F, time_t *_Tstart, clock_t *_Cstart){

  /* Should compute _EP from the rest. Probably requires alignment of
     _E_D with _V_N and of _V_D with _E_N                                 */

  Interval_List IL;
  SPoly *_S_D = NULL, *_S_N = NULL;
  BPoly *_BL = NULL;
  Poset_Element_List PEL_D, PEL_N;
  FaceInfo *_I_D;
  Cone *_C_D, *_C_N;

  _I_D = (FaceInfo *) malloc(sizeof(FaceInfo));
  if (_I_D == NULL) Die("Unable to alloc space for FaceInfo _I_D");
  _C_D = (Cone *) malloc(sizeof(Cone));
  if (_C_D == NULL) Die("Unable to alloc space for Cone _C_D");
  _C_N = (Cone *) malloc(sizeof(Cone));
  if (_C_N == NULL) Die("Unable to alloc space for Cone _C_N");
  
      Make_Incidence(_P_D, _V_D, _E_D, _I_D);
      
      Make_Cone(_C_D, _C_N, _I_D, _P_D);
      PEL_D.n = Num_Pos(_C_D); PEL_N.n = Num_Pos(_C_N);
      
      PEL_D.L = (Poset_Element *) calloc(PEL_D.n, sizeof(Poset_Element));
      if (PEL_D.L == NULL) Die("Unable to alloc space for PEL_D.L");
      PEL_N.L = (Poset_Element *) calloc(PEL_D.n, sizeof(Poset_Element));
      if (PEL_N.L == NULL) Die("Unable to alloc space for PEL_N.L");
      
      Make_PosetList(_C_D, &PEL_D); Make_PosetList(_C_N, &PEL_N);
      if(_F->t) Time_Info(_Tstart, _Cstart, "   BEGIN S-Poly");
      
      _S_D = (SPoly *) calloc(PEL_D.n, sizeof(SPoly ));
      if (_S_D == NULL) Die("Unable to alloc space for SPoly _S_D");
      _S_N = (SPoly *) calloc(PEL_D.n, sizeof(SPoly ));
      if (_S_N == NULL) Die("Unable to alloc space for SPoly _S_N");
	
      Make_S_Poly(_C_N, _V_D, _E_D, _P_D, &PEL_D, _S_D, _F->S, _F->T);
      Make_S_Poly(_C_D, _V_N, _E_N, _P_N, &PEL_N, _S_N, _F->S, _F->T);
	
      if(_F->t) Time_Info(_Tstart, _Cstart, "   BEGIN B-Poly");
	
      IL.L = (Interval *) calloc(((1 + PEL_D.n)/2 + 1)*PEL_D.n, 
		sizeof(Interval));
      if (IL.L == NULL) Die("Unable to alloc space for IL.L");

      Make_Intervallist(&IL, &PEL_D, _C_D);

      _BL = (BPoly *) calloc(IL.n, sizeof(BPoly));
      if (_BL == NULL) Die("Unable to alloc space for _BL");
     
      Make_B_Poly(_C_D, &PEL_D, &IL, _BL);
 
      if(_F->t) Time_Info(_Tstart, _Cstart, "   BEGIN E-Poly");
      SB_To_E(_EP, _C_D, &PEL_D, _BL, &IL, _S_D, _S_N, _codim);
      
      free(PEL_D.L); free(PEL_N.L); free(_S_D); free(_S_N); free(IL.L); 
      free(_BL); free(_I_D); free(_C_D); free(_C_N);
}

void Make_E_Poly(FILE * outFILE, CWS * _W, PolyPointList * _CP,
		 VertexNumList * _CV, EqList * _CE, int *_codim,
		 Flags * _F, int *_D)
{
  time_t Tstart;
  clock_t Cstart;
  int n;
  /* Interval_List IL;
  SPoly *_S_D = NULL, *_S_N = NULL;
  BPoly *_BL = NULL; */
  EPoly EP;
  /* Poset_Element_List PEL_D, PEL_N;*/
  PartList *_PTL;
  PolyPointList *_P = NULL, *_DP = NULL, *_P_D, *_P_N;
  VertexNumList *_V = NULL, *_DV = NULL, *_V_D, *_V_N;
  EqList *_E = NULL, *_DE = NULL, *_E_D, *_E_N;
  /*FaceInfo *_I_D;
    Cone *_C_D, *_C_N;*/
  LInfo *_L = NULL;
  
  /*   ===============	Begin of Static Allocation	===================  */
  _PTL = (PartList *) malloc(sizeof(PartList));
  if (_PTL == NULL) Die("Unable to alloc space for PartList _PTL");
  _P_D = (PolyPointList *) malloc(sizeof(PolyPointList));
  if (_P_D == NULL) Die("Unable to alloc space for PolyPointLis _P_D");
  _P_N = (PolyPointList *) malloc(sizeof(PolyPointList));
  if (_P_N == NULL) Die("Unable to alloc space for PolyPointLis _P_N");
  _V_D = (VertexNumList *) malloc(sizeof(VertexNumList));
  if (_V_D == NULL) Die("Unable to alloc space for VertexNumList _V_D");
  _V_N = (VertexNumList *) malloc(sizeof(VertexNumList));
  if (_V_N == NULL) Die("Unable to alloc space for VertexNumList _V_N");
  _E_D = (EqList *) malloc(sizeof(EqList));
  if (_E_D == NULL) Die("Unable to alloc space for EqList _E_D");
  _E_N = (EqList *) malloc(sizeof(EqList));
  if (_E_N == NULL) Die("Unable to alloc space for EqList _E_N");
  if (_F->Lv || _F->Lp){
    _L = (LInfo *) malloc(sizeof(LInfo));
    if (_L == NULL) Die("Unable to alloc space for LInfo _L");
  }
  if(_F->N){
    _DP = (PolyPointList *) malloc(sizeof(PolyPointList));
    if (_DP == NULL) Die("Unable to alloc space for PolyPointLis _DP");
    _DV = (VertexNumList *) malloc(sizeof(VertexNumList));
    if (_DV == NULL) Die("Unable to alloc space for VertexNumList _DV");
    _DE = (EqList *) malloc(sizeof(EqList));
    if (_DE == NULL) Die("Unable to alloc space for EqList _DE");
  }
  else{
    _P = (PolyPointList *) malloc(sizeof(PolyPointList));
    if (_P == NULL) Die("Unable to alloc space for PolyPointLis _P");
    _V = (VertexNumList *) malloc(sizeof(VertexNumList));
    if (_V == NULL) Die("Unable to alloc space for VertexNumList _V");
    _E = (EqList *) malloc(sizeof(EqList));
    if (_E == NULL) Die("Unable to alloc space for EqList _E");
  }
  /*   ===============	End of Static Allocation	===================  */

  if(_F->N){
    _P=_CP; _E=_CE; _V=_CV;
    Make_Dual_Poly(_P, _V, _E, _DP);
    Find_Equations(_DP, _DV, _DE);
    Sort_PPL(_DP, _DV);
  }
  else{
    _DP=_CP; _DE=_CE; _DV=_CV;
    Make_Dual_Poly(_DP, _DV, _DE, _P);
    Find_Equations(_P, _V, _E);
    Sort_PPL(_P, _V);
  }
  Tstart = time(NULL); Cstart = clock(); _F->Test = 0;fflush(0);
       {	NEF_Flags NF; NF.Sym=_F->Sym;
		NF.noconvex=_F->noconvex; NF.Test=0; NF.Sort=_F->Sort;
		part_nef(_P, _V, _E, _PTL, _codim, &NF);
		if(_F->y && (_PTL->n != 0) && (_W->nw == 0))
			fprintf(outFILE, "%d %d Vertices of Poly in M-lattice:  ", _DP->n, _DV->nv);
		if(!(_F->y && (_PTL->n==0))){
			fprintf(outFILE, "M:%d %d N:%d %d ", _DP->np, _E->ne, _P->np, _V->nv);
			fprintf(outFILE, " codim=%d", *_codim);
			fprintf(outFILE, " #part="); fflush(0);
			fprintf(outFILE, "%d\n",_PTL->n);fflush(0);
		}
	}
		if(_F->y && (_W->nw == 0) && (_PTL->n != 0)) 
      PRINT_TMATRIX(&_DP->x[0][0], _DV->nv, _DP->n, POLY_Dmax);
  if(_F->V)
    Print_VL(_P, _V, "Vertices of P:");
  if((!_F->Lv) && (!_F->n) && _F->Lp)
     PRINT_PL(_P, "Points of Poly in N-Lattice:");
  if(_F->Lv || _F->Lp){
    Make_L(_P, _V, _L ,_F->Lp,_F->Lv);
    Print_L(_L,_F->Lp,_F->Lv);  
  }
  if(_F->f){
    if(_F->Lv)
      PRINT_PL(_P, "Points of Poly in N-Lattice:");
    PRINT_Fibrations(_V, _P, _F  /*, _PTL*/ );
  }
  for (n = 0; n < _PTL->n; n++) {
    /*if ((!_PTL->DirProduct[n]) || _F->Dir){*/
    if(_F->Dir==2) if(!_PTL->DirProduct[n]) continue ;
    Tstart = time(NULL);
    Cstart = clock();
    Make_Gore_Poly(_P, _DP, _P_D, _P_N, _V, _PTL, _codim, &n);
    _PTL->Proj[n] = Remove_Proj(_P_D, /*&n,*/ _codim);
    _PTL->DProj[n] = Remove_Proj(_P_N, /*&n,*/ _codim);
    Find_Equations(_P_D, _V_D, _E_D);
    Find_Equations(_P_N, _V_N, _E_N); 
    Sort_PPL(_P_N, _V_N);

    if (((!_PTL->Proj[n]) || _F->Proj) && ((!_PTL->DirProduct[n]) || _F->Dir)
	&& !_F->p && !_F->y) {
      Make_EN(_P_D, _V_D, _E_N, _codim);
      Compute_E_Poly(&EP, _P_D, _V_D, _E_D, _P_N, _V_N, _E_N, _codim, _F,
		     &Tstart, &Cstart);}
    /*for (i=0;i<=10;i++){
      for (j=0;j<=10;j++) fprintf(outFILE,"%d ",EP.E[i][j]);
      fprintf(outFILE,"\n");}*/

    if(!_F->n && !_F->g && !_F->d && !_F->y){
      Output(	_P, /*_DP,*/ _P_D, _W, &EP, /*_E,*/ 
		_V, &n, _PTL, _codim, outFILE, _F, _D, _L);
      if (((!_PTL->Proj[n]) || _F->Proj) && ((!_PTL->DirProduct[n]) || _F->Dir))
	Time_Info(&Tstart, &Cstart, "");
    }
    if (((!_PTL->Proj[n]) || _F->Proj) && ((!_PTL->DirProduct[n]) || _F->Dir)){
      if(_F->g) PRINT_GORE(_P_D, *_codim, _F->gd, "Points of PG:");
      if(_F->d) PRINT_GORE(_P_N, *_codim, _F->dd, "Points of dual PG:");
    }
  }
  if(!_F->n && !_F->g && !_F->d && !_F->y) Print_Nefinfo(_PTL, /* _F,*/ 
							&Tstart, &Cstart);
  if(_F->n && N_Part(_PTL)) PRINT_PL(_P, "Points of Poly in N-Lattice:");
  /*   ===============	Begin of FREE Static Allocation	===================  */
  free(_PTL);free(_P_D); free(_P_N); free(_V_D); free(_V_N);  free(_E_D); 
  free(_E_N); free(_L);
  if(_F->N){free(_DP); free(_DE); free(_DV);  }
  else{free(_P); free(_E); free(_V);  }
  /*   ===============	End of FREE Static Allocation	===================  */
}

void SL2Z_Make_Poly_UTriang(PolyPointList *P);

void AnalyseGorensteinCone(CWS *_CW,  PolyPointList *_P, VertexNumList *_V, 
			   EqList *_E, int *_codim, Flags * _F){
  /* _P should be the Gorenstein-polytope in M - what is called _P_N in
     certain other parts of the program                                */

  time_t Tstart;
  clock_t Cstart;
  EPoly EP;
  int i, j, k, dim = _P->n - (*_codim * 2) + 1, chi, r=1;
  int h[POLY_Dmax][POLY_Dmax] = {{0}, {0}};
  PolyPointList *_P_D = (PolyPointList *) malloc(sizeof(PolyPointList));
  VertexNumList *_V_D = (VertexNumList *) malloc(sizeof(VertexNumList));
  EqList *_E_D =  (EqList *) malloc(sizeof(EqList));
  /* the Gorenstein-polytope in N-lattice and its vertices and equations */
  EqList *_new_E_D =  (EqList *) malloc(sizeof(EqList));
  /* the equations of _P_D in an order corresponding to the vertices of _P */
  PairMat VPM, VPM_D;

  if ((_P_D == NULL)||(_V_D == NULL)||(_E_D == NULL)||(_new_E_D == NULL))
    Die("Unable to allocate space for _P_D in AnalyseGorensteinCone");

  if (POLY_Dmax  < _P->n + 1){/* only relevant if index == 1 */
    printf("Please increase POLY_Dmax to at least %d = %d + 1\n",
	   (_P->n + 1), _P->n);
    printf("(POLY_Dmax >= dim(cone) = dim(support) + 1 required)\n");
    assert(*_codim == 1);
    exit(0);}
  /* Print_PPL(_P, "_P before sorting");fflush(0); */
  Find_Equations(_P,_V,_E);
  Make_VEPM(_P, _V, _E, VPM);
  Complete_Poly(VPM, _E, _V->nv, _P);
  Sort_PPL(_P,_V);
  /* Print_PPL(_P, "_P after sorting"); */
  /* Print_VL(_P,_V, "_V"); */
  /* Print_EL(_E, &_P->n, 0, "_E"); */ 

  /* Create _P_D from the equations of _P:  
     _P_D is the hull of the generators of the cone dual to the one over _P,
     Gorenstein <-> _P_D lies in a plane at distance 1 from the origin,
     index = modulus of the last component of the equation of this plane */
  _P_D->n = _P->n + 1;
  _P_D->np = _E->ne;
  for (i=0; i<_P_D->np; i++){
    for (j=0; j<_P_D->n; j++) 
      _P_D->x[i][j] = _E->e[i].a[j];
    _P_D->x[i][_P->n] = _E->e[i].c;}
  /* Print_PPL(_P_D, "_P_D "); */
  Find_Equations(_P_D,_V_D,_E_D);
  /* Print_EL(_E_D, &_P_D->n, 0, "_E_D"); */
  if (_E_D->ne == 1){
    if ((_E_D->e[0].c != 1)&&(_E_D->e[0].c != -1)) r = 0;
    else if (abs((int) _E_D->e[0].a[_P->n]) != *_codim){
	printf("Warning: Input has index %d, should be %d!   ", 
	       abs(_E_D->e[0].a[_P->n]), *_codim);
	r = 0;}}
  else {assert(_E_D->ne > _P_D->n); r = 0;}
  
  Tstart = time(NULL); Cstart = clock(); _F->Test = 0;
  for (i = 0; i < _CW->nw; i++) {
    fprintf(outFILE, "%d ", (int) _CW->d[i]);
    for (j = 0; j < _CW->N; j++)
      fprintf(outFILE, "%d ", (int) _CW->W[i][j]);
    if (i + 1 < _CW->nw)
    fprintf(outFILE, " "); } 
  fflush(0);

  if (r) {
    /* shift _P_D into a plane through the origin: */
    for (i=0; i<_P_D->np; i++){
      for (j=0; j<_P_D->n; j++) 
	_P_D->x[i][j] -= _P_D->x[_P_D->np-1][j];}
    /* Print_PPL(_P_D, "_P_D before Make_Poly_UTriang");  */
    SL2Z_Make_Poly_UTriang(_P_D);
    /* Print_PPL(_P_D, "_P_D after Make_Poly_UTriang");  */
    _P_D->n--;
    for (i=0; i<_P_D->np; i++) assert(_P_D->x[i][_P->n] == 0);
    /* Print_PPL(_P_D, "_P_D after reduction"); */
    Find_Equations(_P_D,_V_D,_E_D);
    /* Print_EL(_E_D, &_P_D->n, 0, "_E_D"); */
    Sort_VL(_V_D);
    /* Make_VEPM(_P, _V, _E, VPM); */
    /* Print_Matrix(VPM, _E->ne, _V->nv, "VPM"); */
    Make_VEPM(_P_D, _V_D, _E_D, VPM_D);
    /* Print_Matrix(VPM_D, _E_D->ne, _V_D->nv, "VPM_D"); */
    Complete_Poly(VPM_D, _E_D, _V_D->nv, _P_D);
    /* Print_PPL(_P_D, "_P_D after Complete_Poly"); */
    assert(_E_D->ne == _V->nv);
    assert(_E->ne == _V_D->nv);
    /* Compute _new_E_D->ne by comparing VPM and VPM_D:
       if (VPM[j][i] == VPM_D[k][j]) for all j then _new_E_D->e[i] = _E_D[k] */
    _new_E_D->ne = _E_D->ne;
    for (i=0;i<_V->nv;i++){
      for(k=0;k<_V->nv;k++){
	for (j=0;j<_V_D->nv;j++) if (VPM[j][i] != VPM_D[k][j]) break;
	if (j==_V_D->nv) {_new_E_D->e[i] = _E_D->e[k]; break;} }
      if (k >= _V->nv){
	printf("k = %d, _V->nv = %d\n", k, _V->nv);
	Print_Matrix(VPM, _E->ne, _V->nv, "VPM");
	Print_Matrix(VPM_D, _E_D->ne, _V_D->nv, "VPM_D");}
      assert(k<_V->nv);}
    /* Print_EL(_new_E_D, &_P_D->n, 0, "_new_E_D"); 
       Make_VEPM(_P_D, _V_D, _new_E_D, VPM_D);
       Print_Matrix(VPM_D, _E_D->ne, _V_D->nv, "VPM_D"); */
    if (_F->N){ /* swap M and N  */
      PolyPointList *_auxP = _P_D; 
      VertexNumList *_auxV = _V_D; 
      EqList *_auxE = _new_E_D;
      _P_D = _P; _V_D = _V; _new_E_D = _E; 
      _P = _auxP; _V = _auxV; _E = _auxE;}
    fprintf(outFILE,"M:%d %d ",_P->np,_V->nv);
    fprintf(outFILE,"N:%d %d ",_P_D->np,_V_D->nv); 
    if (!_F->g && !_F->d){
      Compute_E_Poly(&EP, _P_D, _V_D, _new_E_D, _P, _V, _E, _codim, _F,
		     &Tstart, &Cstart);
      chi =  Make_Mirror(&EP, h, _P_D->n + 1, dim);
      fprintf(outFILE, "H:");
      for (i = 1; i < dim; i++)
	fprintf(outFILE, "%d ", h[1][i]);
      fprintf(outFILE, "[%d]", chi);
      for (i = 1; i <=  dim/2; i++)
	if (h[0][i] != 0)
	  fprintf(outFILE, " h%d=%d", i, h[0][i]);
      if(h[0][0] != 1) fprintf(outFILE, " h%d=%d", 0, h[0][0]);
      puts("");
      if (_F->H) PrintDiamond(h, dim);}
    else puts("");
    if(_F->V)
      Print_VL(_P_D, _V_D, "Vertices of support in N:");
    if (_F->g) Print_PPL(_P_D, "Points  of support in N:");
    if (_F->d) Print_PPL(_P, "Points  of support in M:");
    if (_F->t) Time_Info(&Tstart, &Cstart, "");  
    if (_F->N){ /* revert M-N-swap (necessary because of alloc/free)  */
      PolyPointList *_auxP = _P_D; 
      VertexNumList *_auxV = _V_D; 
      EqList *_auxE = _new_E_D;
      _P_D = _P; _V_D = _V; _new_E_D = _E; 
      _P = _auxP; _V = _auxV; _E = _auxE;}}
  else{
    if (_F->N) fprintf(outFILE,"N:%d %d ",_P->np,_V->nv); 
    else fprintf(outFILE,"M:%d %d ",_P->np,_V->nv); 
    fprintf(outFILE,"F:%d ",_E->ne);   
    puts("");
    if ((_F->Rv) || ((_F->V)&&(_F->N))) 
      Print_VL(_P, _V, "Vertices of input polytope:");}
  free(_P_D); free(_E_D); free(_V_D); free(_new_E_D);
}

