#pragma once

#include <memory>
#include <vector>

constexpr int Nef_Max = 500000;
constexpr int NP_Max = 500000;
#define W_Nmax (POLY_Dmax + 1)
constexpr int MAXSTRING = 100;

constexpr int Pos_Max = POLY_Dmax + 2;
constexpr int FIB_POINT_Nmax = VERT_Nmax;

struct LInfo {
  Long W[FIB_Nmax][FIB_POINT_Nmax];
  Long VM[FIB_POINT_Nmax][POLY_Dmax];
  int nw;
  int nv;
  int d;
  int Wmax;
};

struct Poset_Element {
  int num, dim;
};

struct Interval {
  int min, max;
};

struct Interval_List {
  Interval *L;
  std::unique_ptr<Interval[]> L_owner;
  int n;
};

struct Poset {
  Poset_Element x, y;
};

struct Poset_Element_List {
  Poset_Element *L;
  int n;
};

struct Cone {
  int nface[Pos_Max];
  int dim;
  INCI edge[Pos_Max][FACE_Nmax];
};

struct SPoly {
  Long S[2 * Pos_Max];
};

struct BPoly {
  Long B[Pos_Max][Pos_Max];
};

struct EPoly {
  int E[4 * (Pos_Max)][4 * (Pos_Max)];
};

struct AmbiPointList {
  Long x[POINT_Nmax][W_Nmax];
  int N, np;
};

struct PartList {
  int n;
  int nv;
  int codim;
  int S[Nef_Max][VERT_Nmax];
  int DirProduct[Nef_Max];
  int Proj[Nef_Max];
  int DProj[Nef_Max];
};

struct Part {
  int n;
  int nv;
  int S[Nef_Max][VERT_Nmax];
};

struct Flags {
  int n, y, w, p, t, S, Lv, Lp, N, u, d, g, VP, B, T, H, dd, gd, noconvex, Msum,
      Sym, V, Rv, Test, Sort, Dir, Proj, f, G;
};

struct NEF_Flags {
  int noconvex, Sym, Test, Sort;
};

struct Vector {
  Long x[POLY_Dmax];
};

struct DYN_PPL {
  std::vector<Vector> L;
  int n;
  Long np;
};

void part_nef(PolyPointList *, VertexNumList *, EqList *, PartList *, int *,
              NEF_Flags *, FILE *out = stdout);

void Mink_WPCICY(AmbiPointList *_AP_1, AmbiPointList *_AP_2,
                 AmbiPointList *_AP);

void Make_E_Poly(FILE *out, CWS *, PolyPointList *, VertexNumList *, EqList *,
                 int *, Flags *, int *);

int IsDigit(char);

int IntSqrt(int q);

[[noreturn]] void Die(const char *);

void Print_CWS_Zinfo(CWS *);

void AnalyseGorensteinCone(CWS *_CW, PolyPointList *_P, VertexNumList *_V,
                           EqList *_E, int *_codim, Flags *_F,
                           FILE *out = stdout);
