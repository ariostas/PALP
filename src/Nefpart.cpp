#include <palp/Global.h>
#include <palp/Nef.h>

#include <algorithm>
#include <memory>
#include <vector>

/*   ===============	    Typedefs and Headers	===================  */

int GLZ_Make_Trian_NF(Long X[][VERT_Nmax], int *n, int *nv,
                      GL_Long G[][POLY_Dmax]); /* current=best */

void Poly_Sym(PolyPointList *_P, VertexNumList *_V, EqList *_F, int *sym_num,
              int V_perm[][VERT_Nmax]);

/*   ===============	local Typedefs and Headers	===================  */

typedef struct {
  int nv;           /*   #vertices of face */
  int v[VERT_Nmax]; /*   vertices of face */
} VList;

typedef struct {
  int Nv;    /*   #vertices */
  int nf;    /*   #facets */
  VList *vl; /*   vertices of facets */
} FVList;

typedef struct {
  int f;
  int v;
} Step;

typedef struct {
  int s[VERT_Nmax];
} V_Flag;

typedef struct {
  int m[FACE_Nmax];
} M_Rank;

typedef struct {
  int M[POLY_Dmax][POLY_Dmax];
  int d;
  int codim;
} MMatrix;

typedef struct {
  int d;
  int nv;
  Long X[POLY_Dmax][VERT_Nmax];
} XMatrix;

typedef struct {
  int d;
  GL_Long G[POLY_Dmax][POLY_Dmax];
} GMatrix;

typedef struct {
  int Vp[SYM_Nmax][VERT_Nmax];
  int ns;
} SYM;

typedef struct {
  int A[POLY_Dmax];
  int m;
  int M;
} Subset;

int GLZ_Start_Simplex(PolyPointList *_P, VertexNumList *_V, CEqList *_C);

/*   ===============	End of Typedefs and Headers	===================  */

void Print_FVl(FVList *_FVl, const char *comment, FILE *out) {
  int i, j;

  fprintf(out, "%s\n", comment);
  for (i = 0; i < _FVl->nf; i++) {
    for (j = 0; j < _FVl->vl[i].nv; j++)
      fprintf(out, "  %d  ", _FVl->vl[i].v[j]);
    fprintf(out, "nv: %d\n", _FVl->vl[i].nv);
    fflush(0);
  }
}

void Print_M(MMatrix *_M, int *_nf, const char *comment, FILE *out) {
  int i, j, k;

  fprintf(out, "%s\n", comment);
  for (i = 0; i < *_nf; i++) {
    fprintf(out, "\n\n facet %d:\n", i);
    for (j = 0; j < _M[i].codim; j++) {
      fprintf(out, "m[%d]: ", j);
      for (k = 0; k < _M[i].d; k++)
        fprintf(out, "%d ", _M[i].M[k][j]);
      fprintf(out, "\n");
    }
  }
}

void Dir_Product(PartList *_PTL, VertexNumList *_V, PolyPointList *_P) {

  int i, j, k, d, rank;

  CEqList CEtemp;
  VertexNumList Vtemp;

  auto _PV_up = std::make_unique<PolyPointList>();
  PolyPointList *_PV = _PV_up.get();

  _PV->n = _P->n;
  for (i = 0; i < _PTL->n; i++) {
    rank = 0;
    for (j = 0; j < _PTL->codim; j++) {
      _PV->np = 0;
      for (k = 0; k < _PTL->nv; k++)
        if (_PTL->S[i][k] == j) {
          for (d = 0; d < _P->n; d++)
            _PV->x[_PV->np][d] = _P->x[_V->v[k]][d];
          _PV->np += 1;
        }
      for (d = 0; d < _P->n; d++)
        _PV->x[_PV->np][d] = 0;
      _PV->np += 1;
      rank += GLZ_Start_Simplex(_PV, &Vtemp, &CEtemp);
    }
    if (rank == _P->n)
      _PTL->DirProduct[i] = 1;
    else
      _PTL->DirProduct[i] = 0;
  }
}

void Set_To_Vlist(int *_S, VertexNumList *_V, PolyPointList *_P,
                  PolyPointList *_PV, Subset *_Pset) {
  int i, d;

  _PV->n = _P->n;
  _PV->np = 0;
  for (i = 0; i < _V->nv; i++)
    if (_Pset->A[_S[i]] == 1) {
      for (d = 0; d < _P->n; d++)
        _PV->x[_PV->np][d] = _P->x[_V->v[i]][d];
      _PV->np++;
    }
  for (d = 0; d < _P->n; d++)
    _PV->x[_PV->np][d] = 0;
  _PV->np += 1;
}

int New_Set(int i, Subset *_Pset, Subset *_CPset) {

  if (_Pset->A[i] == 0) {
    _Pset->A[i] = 1;
    _Pset->m++;
    _CPset->A[i] = 0;
    _CPset->m--;
    return 1;
  } else
    return 0;
}

void Old_Set(int i, Subset *_Pset, Subset *_CPset) {
  _Pset->A[i] = 0;
  _Pset->m--;
  _CPset->A[i] = 1;
  _CPset->m++;
}

void Select_Set(int *_S, VertexNumList *_V, PolyPointList *_P,
                PolyPointList *_PV, CEqList *_CEtemp, VertexNumList *_Vtemp,
                int *_DirFlag, Subset *_Pset, Subset *_CPset, int nmax_old) {
  int rank, nmax;

  Set_To_Vlist(_S, _V, _P, _PV, _Pset);
  rank = GLZ_Start_Simplex(_PV, _Vtemp, _CEtemp);
  Set_To_Vlist(_S, _V, _P, _PV, _CPset);
  rank += GLZ_Start_Simplex(_PV, _Vtemp, _CEtemp);

  /*{int i; for(i=0;i<_Pset->M;i++) printf("   %d",_Pset->A[i]);
   * printf("\n");}*/
  if (rank == _P->n)
    *_DirFlag = 1;
  else {
    if (_Pset->m < (_Pset->M / 2 + (_Pset->M % 2))) {
      nmax = (nmax_old + 1);
      while (!*_DirFlag && (nmax < _Pset->M)) {
        /*printf("aaa");*/
        if (New_Set(nmax, _Pset, _CPset)) {
          Select_Set(_S, _V, _P, _PV, _CEtemp, _Vtemp, _DirFlag, _Pset, _CPset,
                     nmax);
          Old_Set(nmax, _Pset, _CPset);
        }
        nmax++;
      }
    }
  }
}

void REC_Dir_Product(PartList *_PTL, VertexNumList *_V, PolyPointList *_P) {

  int i, j, k;
  CEqList *_CEtemp;
  VertexNumList *_Vtemp;
  PolyPointList *_PV;
  Subset Pset, CPset;

  auto _PV_up = std::make_unique<PolyPointList>();
  _PV = _PV_up.get();
  auto _CEtemp_up = std::make_unique<CEqList>();
  _CEtemp = _CEtemp_up.get();
  auto _Vtemp_up = std::make_unique<VertexNumList>();
  _Vtemp = _Vtemp_up.get();

  for (i = 0; i < _PTL->n; i++) {
    _PTL->DirProduct[i] = 0;
    j = 0;
    /*printf("\n*********************************\n");*/
    while ((j < _PTL->codim) && !_PTL->DirProduct[i]) {
      for (k = 0; k < _PTL->codim; k++) {
        if (k == j) {
          Pset.A[k] = 1;
          CPset.A[k] = 0;
        } else {
          Pset.A[k] = 0;
          CPset.A[k] = 1;
        }
      }
      Pset.m = 1;
      Pset.M = _PTL->codim;
      CPset.m = (_PTL->codim - 1);
      CPset.M = _PTL->codim;
      Select_Set(_PTL->S[i], _V, _P, _PV, _CEtemp, _Vtemp, &_PTL->DirProduct[i],
                 &Pset, &CPset, 0);
      j++;
    }
  }
}

int COMP_S(int SA[], int SB[], int *_nv) {

  int eq_flag = 1, i = 0;

  while ((i < *_nv) && eq_flag) {
    if ((SA[i] > SB[i]) || (SA[i] < SB[i]))
      eq_flag = 0;
    i++;
  }
  if (eq_flag)
    return 0;
  else {
    if (SA[i - 1] > SB[i - 1])
      return 1;
    else
      return -1;
  }
}

void NForm_S(int S[], int *_nv) {

  int s[POLY_Dmax], d = 0, n_flag, i, j;

  for (i = 0; i < *_nv; i++) {
    n_flag = 1;
    j = 0;
    while ((j < i) && n_flag) {
      if (S[j] == S[i])
        n_flag = 0;
      j++;
    }
    if (n_flag) {
      s[S[i]] = d;
      d++;
    }
  }
  for (i = 0; i < *_nv; i++)
    S[i] = s[S[i]];
}

int Bisection_PTL(PartList *_PTL, int s[], int S[]) {

  int min_pos = -1, max_pos = _PTL->n, pos = 1, c = 1;

  while ((max_pos - min_pos) > 1) {
    pos = (max_pos + min_pos) / 2;
    c = COMP_S(S, _PTL->S[s[pos]], &_PTL->nv);
    if (c == 1)
      min_pos = pos;
    else if (c == -1)
      max_pos = pos;
    else
      min_pos = max_pos;
  }
  if (c != 0) {
    fputs("Error: Bisection_PTL partition not found\n", stderr);
    exit(1);
  }
  return s[pos];
}

void Bubble_PTL(PartList *_PTL, int s[]) {

  int i, diff = 0;

  for (i = 0; i < _PTL->n; i++)
    s[i] = i;

  std::sort(s, s + _PTL->n, [nv = &_PTL->nv, &PTL = *_PTL](int a, int b) {
    return COMP_S(PTL.S[a], PTL.S[b], nv) < 0;
  });

  for (i = 1; i < _PTL->n; i++) {
    if (COMP_S(_PTL->S[s[i - 1]], _PTL->S[s[i]], &_PTL->nv) == 0)
      diff += 1;
    else
      s[i - diff] = s[i];
  }
  _PTL->n = _PTL->n - diff;
}

void Remove_Sym(SYM *_VP, PartList *_PTL, PartList *_S_PTL) {

  int n = 1, i, j, k;
  std::vector<int> _s(_PTL->n);
  std::vector<int> _s_sym(_VP->ns);
  std::vector<int> _p(_PTL->n);
  auto _SYM_PTL = std::make_unique<PartList>();

  if (Nef_Max < SYM_Nmax)
    Die("\nNeed Nef_Max >= SYM_Nmax!!!\n");
  for (i = 0; i < _PTL->n; i++) {
    _p[i] = 0;
    _s[i] = 0;
  }
  for (i = 0; i < _PTL->n; i++)
    NForm_S(_PTL->S[i], &_PTL->nv);
  Bubble_PTL(_PTL, _s.data());
  for (i = 0; i < _PTL->n; i++) {
    if (_p[_s[i]] == 0) {
      _p[_s[i]] = n;
      for (j = 0; j < _VP->ns; j++) {
        for (k = 0; k < _PTL->nv; k++)
          _SYM_PTL->S[j][k] = _PTL->S[_s[i]][_VP->Vp[j][k]];
      }
      _SYM_PTL->n = _VP->ns;
      _SYM_PTL->nv = _PTL->nv;
      for (j = 0; j < _VP->ns; j++)
        NForm_S(_SYM_PTL->S[j], &_SYM_PTL->nv);
      Bubble_PTL(_SYM_PTL.get(), _s_sym.data());
      for (j = 1; j < _VP->ns; j++)
        if (COMP_S(_SYM_PTL->S[_s_sym[j]], _SYM_PTL->S[_s_sym[j - 1]],
                   &_PTL->nv) != 0)
          _p[Bisection_PTL(_PTL, _s.data(), _SYM_PTL->S[_s_sym[j]])] = n;
      n++;
    }
  }
  n--;
  _S_PTL->nv = _PTL->nv;
  _S_PTL->n = 0;
  _S_PTL->codim = _PTL->codim;
  while (n > 0) {
    i = 0;
    while (_p[_s[i]] != n)
      i++;
    for (j = 0; j < _PTL->nv; j++)
      _S_PTL->S[_S_PTL->n][j] = _PTL->S[_s[i]][j];
    _S_PTL->n += 1;
    n--;
  }
}

void M_TO_MM(MMatrix *_M, MMatrix *_MM, GMatrix *_G, int *_nf) {

  int i, l, c, j;

  for (i = 0; i < *_nf; i++) {
    for (l = 0; l < _M[i].d; l++)
      for (c = 0; c < _M[i].codim; c++) {
        _MM[i].M[l][c] = 0;
        for (j = 0; j < _M[i].d; j++)
          _MM[i].M[l][c] += _G[i].G[j][l] * _M[i].M[j][c];
      }
    _MM[i].d = _M[i].d;
    _MM[i].codim = _M[i].codim;
  }
}

int Convex_Check(MMatrix *_M, GMatrix *_G, XMatrix *_X, int S[], FVList *_FVl,
                 NEF_Flags *_F) {
  int c_flag = 1, i, j, k, l, d, IP;

  if (_F->noconvex)
    return 1;

  std::vector<MMatrix> _MM_vec(_FVl->nf);
  MMatrix *_MM = _MM_vec.data();

  M_TO_MM(_M, _MM, _G, &_FVl->nf);

  i = 0;
  while ((i < _FVl->nf) && c_flag) {
    j = 0;
    while ((j < _FVl->nf) && c_flag) {
      if (i != j) {
        k = 0;
        while ((k < _MM[i].codim) && c_flag) {
          l = 0;
          while ((l < _X[i].nv) && c_flag) {
            if (S[_FVl->vl[i].v[l]] == k)
              IP = 1;
            else
              IP = 0;
            for (d = 0; d < _MM[i].d; d++)
              IP += -_X[i].X[d][l] * _MM[j].M[d][k];
            if (IP < 0)
              c_flag = 0;
            l++;
          }
          k++;
        }
      }
      j++;
    }
    i++;
  }
  if (c_flag && _F->Test)
    Print_M(_MM, &_FVl->nf, "M-Matrix:", outFILE);
  return c_flag;
}

int Codim_Check(int S[], int *_codim, int *_Nv) {

  int gp[POLY_Dmax] = {0}, gp_flag = 1, i, j = 0;

  while ((j < *_codim) && gp_flag) {
    i = 0;
    while ((i < *_Nv) && (gp[j] == 0)) {
      if (S[i] == j)
        gp[j] = 1;
      i++;
    }
    gp_flag = gp[j];
    j++;
  }
  return gp_flag;
}

int Fix_M(Step *_step, XMatrix *_Y, MMatrix *_M, int S[], M_Rank *_MR,
          FVList *_FVl) {

  int f_flag = 1, i = 0, j, m, IP;

  while (f_flag && (i != _M[_step->f].codim)) {
    if (i == S[_FVl->vl[_step->f].v[_step->v]])
      IP = 1;
    else
      IP = 0;
    for (j = 0; j < _MR->m[_step->f]; j++)
      IP += -_Y[_step->f].X[j][_step->v] * _M[_step->f].M[j][i];
    m = (IP / _Y[_step->f].X[_MR->m[_step->f]][_step->v]);
    if ((IP - m * _Y[_step->f].X[_MR->m[_step->f]][_step->v]) != 0)
      f_flag = 0;
    else
      _M[_step->f].M[_MR->m[_step->f]][i] = m;
    i++;
  }
  return f_flag;
}

int Check_Consistence(Step *_step, XMatrix *_Y, MMatrix *_M, int S[],
                      M_Rank *_MR, FVList *_FVl) {

  int c_flag = 1, i = 0, j, IP;

  while (c_flag && (i != _M[_step->f].codim)) {
    IP = 0;
    for (j = 0; j < _MR->m[_step->f]; j++)
      IP += _Y[_step->f].X[j][_step->v] * _M[_step->f].M[j][i];
    if (i == S[_FVl->vl[_step->f].v[_step->v]]) {
      if (IP != 1)
        c_flag = 0;
    } else {
      if (IP != 0)
        c_flag = 0;
    }
    i++;
  }
  return c_flag;
}

int New_V(V_Flag *_VF, int *_i) {

  int new_flag = 0;

  if (!_VF->s[*_i]) {
    _VF->s[*_i] = 1;
    new_flag = 1;
  }
  return new_flag;
}

int Next_Step(FVList *_FVl, Step *_step) {

  int step_flag = 1;

  if (_step->v == _FVl->vl[_step->f].nv - 1) {
    _step->f += 1;
    _step->v = 0;
  } else
    _step->v += 1;
  if (_step->f == _FVl->nf)
    step_flag = 0;
  return step_flag;
}

void New_VFlag(V_Flag *_VF, /* int *_Nv, */ int *_n) { _VF->s[*_n] = 1; }

void Old_VFlag(V_Flag *_VF, /* int *_Nv, */ int *_n) { _VF->s[*_n] = 0; }

void Raise_M_Rank(M_Rank *_MR, int *_facet) { _MR->m[*_facet] += 1; }

void Lower_M_Rank(M_Rank *_MR, int *_facet) { _MR->m[*_facet] -= 1; }

void Zero_M_Rank(M_Rank *_MR, int *_facet) { _MR->m[*_facet] = 0; }

void Initial_Conditions(MMatrix *_M, XMatrix *_Y, M_Rank *_MR, Step *_step,
                        FVList *_FVl, V_Flag *_VF, int S[], int *_codim,
                        int *_dim, PartList *_PTL) {

  int i;
  for (i = 0; i < _FVl->nf; i++) {
    _M[i].codim = *_codim;
    _M[i].d = *_dim;
    _MR->m[i] = 0;
  }
  _step->f = 0;
  _step->v = 0;
  for (i = 0; i < _FVl->Nv; i++)
    _VF->s[i] = 0;
  _VF->s[_FVl->vl[0].v[0]] = 1;
  _MR->m[0] = 1;
  if ((_Y[0].X[0][0] != 1) && (_Y[0].X[0][0] != -1)) {
    fprintf(stderr, "Error: Init_Matrix first coordinate is %ld, expected ±1\n",
            (long)_Y[0].X[0][0]);
    exit(1);
  }
  for (i = 0; i < *_codim; i++)
    _M[0].M[0][i] = 0;
  _M[0].M[0][0] = _Y[0].X[0][0];
  S[_FVl->vl[0].v[0]] = 0;
  _PTL->n = 0;
  _PTL->nv = _FVl->Nv;
  _PTL->codim = *_codim;
}

void INCI_To_VList(INCI *_X, VList *_Vl, int *_N) {
  /* INCI X -> {3, 5, 7, ...} ... List of Vertices in
     VertexNumList of corresponding face */

  int i, n = 0;
  INCI Y = *_X;

  for (i = 0; i < *_N; i++) {
    if (INCI_M2(Y)) {
      _Vl->v[n] = *_N - i - 1;
      n++;
    }
    Y = INCI_D2(Y);
  }
  _Vl->nv = n;
}

void INCI_To_FVList(FaceInfo *_I, PolyPointList *_P, FVList *_FVl) {
  /* FaceInfo I ->  {3, 5, 7, ...}, {1, 4, 8, ...}, ...
     Lists of Vertices in VertexNumList of all facets  */

  int i;

  for (i = 0; i < _I->nf[_P->n - 1]; i++)
    INCI_To_VList(&_I->v[_P->n - 1][i], &_FVl->vl[i], &_I->nf[0]);
  _FVl->nf = _I->nf[_P->n - 1];
  _FVl->Nv = _I->nf[0];
}

void Sort_FVList(FVList *_FVl, FVList *_FVl_new, int f[]) {
  /* gives f[]  = LIST OF SORTED FACETS IN _FVl:
     f[0]: facet with max # of Vertices
     f[k != 0]: facet with max # of Vertices in f[0],...,f[k-1]
     _FVl->vl[k] -> _FVl->vl[f[k]] =
     {v in f[0],...,f[k-1]}U{v not in f[0],...,f[k-1]} =: FVl_new->vl[k] */

  int i, j, k, l, n, Vequal, vequal, newfacet, newvertex, v[VERT_Nmax],
      w[VERT_Nmax], nv, nw;

  f[0] = 0;
  for (i = 1; i < _FVl->nf; i++) /* first facet */
    if (_FVl->vl[i].nv > _FVl->vl[f[0]].nv)
      f[0] = i;
  for (i = 0; i < _FVl->vl[f[0]].nv; i++)
    v[i] = _FVl->vl[f[0]].v[i];
  nv = _FVl->vl[f[0]].nv;
  for (n = 1; n < _FVl->nf; n++) { /* next facet */
    Vequal = 0;
    for (i = 0; i < _FVl->nf; i++) {
      newfacet = 1;
      k = 0;
      while ((k < n) && newfacet) {
        if (f[k] == i)
          newfacet = 0;
        k++;
      }
      if (newfacet) {
        vequal = 0;
        for (j = 0; j < _FVl->vl[i].nv; j++) {
          newvertex = 1;
          k = 0;
          while ((k < nv) && newvertex) {
            if (v[k] == _FVl->vl[i].v[j])
              newvertex = 0;
            k++;
          }
          if (!newvertex)
            vequal++;
        }
        if (vequal >= Vequal) {
          Vequal = vequal;
          f[n] = i;
        }
      }
    }
    nw = 0;
    for (j = 0; j < _FVl->vl[f[n]].nv; j++) {
      newvertex = 1;
      k = 0;
      while ((k < nv) && newvertex) {
        if (v[k] == _FVl->vl[f[n]].v[j])
          newvertex = 0;
        k++;
      }
      if (newvertex) {
        v[nv] = _FVl->vl[f[n]].v[j];
        nv++;
      } else {
        w[nw] = _FVl->vl[f[n]].v[j];
        nw++;
      }
    }
    for (l = 0; l < nw; l++)
      _FVl->vl[f[n]].v[l] = w[l];
    for (l = 0; l < _FVl->vl[f[n]].nv - nw; l++)
      _FVl->vl[f[n]].v[nw + l] = v[nv - l - 1];
  }
  for (i = 0; i < _FVl->nf; i++) {
    for (j = 0; j < _FVl->vl[f[i]].nv; j++)
      _FVl_new->vl[i].v[j] = _FVl->vl[f[i]].v[j];
    _FVl_new->vl[i].nv = _FVl->vl[f[i]].nv;
  }
  _FVl_new->nf = _FVl->nf;
  _FVl_new->Nv = _FVl->Nv;
}

void Make_Matrix(XMatrix *_X, XMatrix *_Y, VList *_Vl, PolyPointList *_P,
                 VertexNumList *_V) {
  /* Vl -> Matrix X */

  int d, j;

  for (d = 0; d < _P->n; d++)
    for (j = 0; j < _Vl->nv; j++) {
      _X->X[d][j] = _P->x[_V->v[_Vl->v[j]]][d];
      _Y->X[d][j] = _P->x[_V->v[_Vl->v[j]]][d];
    }
  _X->d = _P->n;
  _X->nv = _Vl->nv;
  _Y->d = _P->n;
  _Y->nv = _Vl->nv;
}

void Copy_PTL(PartList *_IN_PTL, PartList *_OUT_PTL) {

  int i, j;
  for (i = 0; i < _IN_PTL->n; i++)
    for (j = 0; j < _IN_PTL->nv; j++)
      _OUT_PTL->S[i][j] = _IN_PTL->S[i][j];
  _OUT_PTL->n = _IN_PTL->n;
  _OUT_PTL->nv = _IN_PTL->nv;
  _OUT_PTL->codim = _IN_PTL->codim;
}

void Select_Sv(int S[], V_Flag *_VF, MMatrix *_M, GMatrix *_G, XMatrix *_X,
               XMatrix *_Y, M_Rank *_MR, FVList *_FVl, Step step,
               PartList *_PTL, NEF_Flags *_F, FILE *out) {

  int i;
  if (Next_Step(_FVl, &step)) {
    if (step.v == 0)
      Zero_M_Rank(_MR, &step.f);
    if (!New_V(_VF, &_FVl->vl[step.f].v[step.v])) {
#ifndef NDEBUG
      if (_F->Test) {
        char c;
        fprintf(out, "\nold vertex f:%d v:%d\n", step.f, step.v);
        if (scanf("%c", &c) != 1)
          c = '\n';
      }
#endif
      if ((_MR->m[step.f] == _M[step.f].d) ||
          (_Y[step.f].X[_MR->m[step.f]][step.v] == 0)) {
        if (Check_Consistence(&step, _Y, _M, S, _MR, _FVl)) {
          if (_F->Test)
            fprintf(out, "f:%d v:%d consistent\n", step.f, step.v);
          Select_Sv(S, _VF, _M, _G, _X, _Y, _MR, _FVl, step, _PTL, _F, out);
        }
      } else {
        if (Fix_M(&step, _Y, _M, S, _MR, _FVl)) {
          if (_F->Test)
            fprintf(out, "f:%d v:%d fixed\n", step.f, step.v);
          Raise_M_Rank(_MR, &step.f);
          Select_Sv(S, _VF, _M, _G, _X, _Y, _MR, _FVl, step, _PTL, _F, out);
        }
      }
    } else {
      if (_F->Test) {
        fprintf(out, "\nnew vertex f:%d v:%d", step.f, step.v);
      }
      New_VFlag(_VF, /* &_FVl->Nv,*/ &_FVl->vl[step.f].v[step.v]);
      for (i = 0; i < _M[step.f].codim; i++) {
        S[_FVl->vl[step.f].v[step.v]] = i;
        if (_F->Test)
          fprintf(out, "to partition: %d\n", i);
        if ((_MR->m[step.f] == _M[step.f].d) ||
            (_Y[step.f].X[_MR->m[step.f]][step.v] == 0)) {
          if (Check_Consistence(&step, _Y, _M, S, _MR, _FVl)) {
            if (_F->Test)
              fprintf(out, "f:%d v:%d consistent\n", step.f, step.v);
            Select_Sv(S, _VF, _M, _G, _X, _Y, _MR, _FVl, step, _PTL, _F, out);
          }
        } else {
          if (Fix_M(&step, _Y, _M, S, _MR, _FVl)) {
            if (_F->Test)
              fprintf(out, "f:%d v:%d fixed\n", step.f, step.v);
            Raise_M_Rank(_MR, &step.f);
            Select_Sv(S, _VF, _M, _G, _X, _Y, _MR, _FVl, step, _PTL, _F, out);
            Lower_M_Rank(_MR, &step.f);
          }
        }
      }
      Old_VFlag(_VF, /* &_FVl->Nv, */ &_FVl->vl[step.f].v[step.v]);
    }
  } else {
    if (_F->Test) {
      fprintf(out, "\n********************************************\n");
      for (i = 0; i < _FVl->Nv; i++)
        fprintf(out, " %d ", S[i]);
      fprintf(out, "\n");
      fprintf(out, "**********************************************\n");
      fflush(0);
    }
    if (Codim_Check(S, &_M[step.f - 1].codim, &_FVl->Nv))
      if (Convex_Check(_M, _G, _X, S, _FVl, _F)) {
        if (_PTL->n >= Nef_Max) {
          fputs("Error: Recursive_Nef partition list overflow\n", stderr);
          exit(1);
        }
        for (i = 0; i < _FVl->Nv; i++)
          _PTL->S[_PTL->n][i] = S[i];
        _PTL->n += 1;
      }
  }
}

void part_nef(PolyPointList *_P, VertexNumList *_V, EqList *_E,
              PartList *_OUT_PTL, int *_codim, NEF_Flags *_F) {

  FaceInfo I;
  FVList FVl_temp, FVl;
  XMatrix *_X, *_Y;
  MMatrix *_M;
  GMatrix *_G;
  Step step;
  V_Flag VF;
  M_Rank MR;
  int i, f[FACE_Nmax], S[VERT_Nmax] = {0};

  PartList *_PTL;
  auto _PTL_up = std::make_unique<PartList>();
  _PTL = _PTL_up.get();

  Make_Incidence(_P, _V, _E, &I);

  std::vector<VList> FVl_vl_vec(I.nf[_P->n - 1]);
  FVl.vl = FVl_vl_vec.data();

  std::vector<VList> FVl_temp_vl_vec;
  if (_F->Sort) {
    FVl_temp_vl_vec.resize(I.nf[_P->n - 1]);
    FVl_temp.vl = FVl_temp_vl_vec.data();
    INCI_To_FVList(&I, _P, &FVl_temp);
    Sort_FVList(&FVl_temp, &FVl, f);
  } else
    INCI_To_FVList(&I, _P, &FVl);
  if (_F->Test) {
    Print_VL(_P, _V, "Vertices of P:");
    Print_FVl(&FVl, "Facets/Vertices:", outFILE);
  }
  std::vector<XMatrix> _X_vec(FVl.nf);
  _X = _X_vec.data();
  std::vector<XMatrix> _Y_vec(FVl.nf);
  _Y = _Y_vec.data();
  std::vector<MMatrix> _M_vec(FVl.nf);
  _M = _M_vec.data();
  std::vector<GMatrix> _G_vec(FVl.nf);
  _G = _G_vec.data();

  for (i = 0; i < FVl.nf; i++) {
    Make_Matrix(&_X[i], &_Y[i], &FVl.vl[i], _P, _V);
    GLZ_Make_Trian_NF(_Y[i].X, &_P->n, &FVl.vl[i].nv, _G[i].G);
  }
  Initial_Conditions(_M, _Y, &MR, &step, &FVl, &VF, S, _codim, &_P->n, _PTL);
  Select_Sv(S, &VF, _M, _G, _X, _Y, &MR, &FVl, step, _PTL, _F, outFILE);
  if (_F->Sym) {
    auto _VP = std::make_unique<SYM>();

    Poly_Sym(_P, _V, _E, &_VP->ns, _VP->Vp, outFILE);
    Remove_Sym(_VP.get(), _PTL, _OUT_PTL);
  } else
    Copy_PTL(_PTL, _OUT_PTL);
  /*Dir_Product(_OUT_PTL, _V, _P);*/
  if (*_codim > 1)
    REC_Dir_Product(_OUT_PTL, _V, _P);
}
