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

#include <palp/Global.h>
#include <palp/LG.h>

#include <memory>
#include <palp/Mori.h>

/*==========================================================*/

PalpContext palpContext;

void PrintUsage(char *c) {
  printf("This is ``%s'':  star triangulations of a polytope P* in N\n", c);
  printf("                     Mori cone of the corresponding toric ambient "
         "spaces\n");
  printf("                     intersection rings of embedded (CY) "
         "hypersurfaces\n");
  printf("Usage:   %s [-<Option-string>] [in-file [out-file]]\n", c);
  printf("Options (concatenate any number of them into <Option-string>):\n");

  printf("    -h      print this information \n");
  printf("    -f      use as filter\n");
  printf(
      "    -g      general output: triangulation and Stanley-Reisner ideal\n");
  printf("    -I      incidence information of the facets (ignoring IPs of "
         "facets)\n");
  printf("    -m      Mori generators of the ambient space\n");
  printf(
      "    -P      IP-simplices among points of P* (ignoring IPs of facets)\n");
  printf("    -K      points of P* in Kreuzer polynomial form\n");
  printf("    -b      arithmetic genera and Euler number\n");
  printf("    -i      intersection ring\n");
  printf("    -c      Chern classes of the (CY) hypersurface\n");
  printf("    -t      triple intersection numbers\n");
  printf("    -d      topological information on toric divisors & del Pezzo "
         "conditions\n");
  printf("    -a      all of the above except h, f, I and K\n");
  printf("    -D      lattice polytope points of P* as input (default CWS)\n");
  printf("    -H      arbitrary (also non-CY) hypersurface `H = c1*D1 + c2*D2 "
         "+ ...'\n");
  printf("            input: coefficients `c1 c2 ...'\n");
  printf("    -M      manual input of triangulation\n");
  puts("Input: 1) standard: degrees and weights `d1 w11 w12 ... d2 w21 w22 "
       "...'");
  puts("       2) alternative (use -D): `d np' or `np d' (d=Dimension, "
       "np=#[points])");
  puts("                                and (after newline) np*d coordinates");
  puts("Output:   as specified by options");
}

int main(int narg, char *fn[]) {

  int n = 0, i, k;

  /* flags */
  MORI_Flags Flag;

  Flag.FilterFlag = 0;            // filter
  Flag.g = 0;                     // -g: general output
  Flag.m = 0;                     // -m: Mori cone
  Flag.P = 0;                     // -P: IP simplices
  Flag.K = 0;                     // -K: Newton polynomial
  Flag.i = 0;                     // -i: intersection ring
  Flag.t = 0;                     // -t: triple intersection number
  Flag.c = 0;                     // -c: Chern classes
  Flag.d = 0;                     // -d; del Pezzo
  Flag.a = 0;                     // -a: all of the above except h,f and K
  Flag.b = 0;                     // -b: Hodge number of toric div
  Flag.D = 0;                     // -D: dual poly as input
  Flag.H = 0;                     // -H: arbitrary hypersurface
  Flag.I = 0;                     // -I: incidence information
  Flag.M = 0;                     // -M: allows to insert a triangulation
  Flag.Read_HyperSurfCounter = 0; // see Mori.h for description
  char c;

  auto CW_up = std::make_unique<CWS>();
  CWS *CW = CW_up.get();

  VertexNumList V;
  auto E_up = std::make_unique<EqList>();
  EqList *E = E_up.get();
  auto DE_up = std::make_unique<EqList>();
  EqList *DE = DE_up.get();

  auto _P_up = std::make_unique<PolyPointList>();
  PolyPointList *_P = _P_up.get();
  auto _DP_up = std::make_unique<PolyPointList>();
  PolyPointList *_DP = _DP_up.get();

  PairMat PM, DPM;

  CW->nw = 0;

  while (narg > ++n) {
    if (fn[n][0] != '-')
      break;
    k = 0;
    while ((c = fn[n][++k]) != '\0') {
      if (c == 'h') {
        PrintUsage(fn[0]);
        exit(1);
      }
      if (c == 'f')
        Flag.FilterFlag = 1;
      if (c == 'g')
        Flag.g = 1;
      if (c == 'm')
        Flag.m = 1;
      if (c == 'P')
        Flag.P = 1;
      if (c == 'K')
        Flag.K = 1;
      if (c == 'i')
        Flag.i = 1;
      if (c == 't')
        Flag.t = 1;
      if (c == 'c')
        Flag.c = 1;
      if (c == 'd')
        Flag.d = 1;
      if (c == 'a')
        Flag.a = 1;
      if (c == 'b')
        Flag.b = 1;
      if (c == 'D')
        Flag.D = 1;
      if (c == 'H')
        Flag.H = 1;
      if (c == 'I')
        Flag.I = 1;
      if (c == 'M')
        Flag.M = 1;
    }
  }
  n--;

  /*if ((Flag.M)&&(!Flag.D)){
    puts("-M works only when combined with -D!");
    exit(1);}*/
  if (Flag.g + Flag.m + Flag.P + Flag.K + Flag.i + Flag.t + Flag.c + Flag.d +
          Flag.a + Flag.b + Flag.H + Flag.I ==
      0)
    Flag.g = 1;

  if (Flag.a) {
    Flag.g = 1;
    Flag.m = 1;
    Flag.P = 1;
    // Flag.K=1;
    Flag.i = 1;
    Flag.t = 1;
    Flag.c = 1;
    Flag.d = 1;
    Flag.b = 1;
  }

  if (Flag.H == 1 && (Flag.g + Flag.m + Flag.P + Flag.K + Flag.i + Flag.t +
                          Flag.c + Flag.d + Flag.a + Flag.b + Flag.I ==
                      0)) {
    Flag.b = 1;
    // Flag.g=1;
  }

  FILE *in, *out;
  if (Flag.FilterFlag) {
    in = NULL;
    out = stdout;
  }

  else {
    if (narg > ++n)
      in = fopen(fn[n], "r");
    else
      in = stdin;

    if (in == NULL) {
      printf("Input file %s not found!\n", fn[n]);
      exit(1);
    }

    if (narg > ++n)
      out = fopen(fn[n], "w");
    else
      out = stdout;
  }

  while ((Flag.D ? Read_PP(_P, in) : Read_CWS(CW, _P, in, out))) {
    if (!Ref_Check(_P, &V, E)) {
      fprintf(out, "Input not reflexive!\n");
      continue;
    }
    if (Flag.D == 0) { /* dualize: _P should become the N-polytope! */
      if (!EL_to_PPL(E, _P, &_P->n)) {
        fputs("Error: mori could not convert equation list to N-polytope\n",
              stderr);
        exit(1);
      }
      if (!Ref_Check(_P, &V, E)) {
        fputs("Error: mori N-polytope is not reflexive\n", stderr);
        exit(1);
      }
    }
    Sort_VL(&V);
    if (!(Flag.D && Flag.M)) {
      Make_VEPM(_P, &V, E, PM);
      Complete_Poly(PM, E, V.nv, _P);
      for (i = V.nv; i < _P->np - 1; i++)
        if (Vec_is_zero(_P->x[i], _P->n)) {
          Swap_Vecs(_P->x[i], _P->x[_P->np - 1], _P->n);
          break;
        }
    } else {
      for (i = 0; i < _P->np; i++)
        if (Vec_is_zero(_P->x[i], _P->n)) {
          Swap_Vecs(_P->x[i], _P->x[_P->np - 1], _P->n);
          break;
        }
      if (i == _P->np) {
        for (k = 0; k < _P->n; k++)
          _P->x[_P->np][k] = 0;
        _P->np++;
      }
    }
    if (Flag.M) {
      if (POLY_Dmax < (_P->np - _P->n)) {
        printf("Please increase POLY_Dmax to at least %d = %d - %d - 1\n",
               (_P->np - _P->n - 1), _P->np, _P->n);
        printf("(%s -M requires POLY_Dmax >= #(points) - dim N -1)\n", fn[0]);
        exit(1);
      }
    }
    HyperSurfDivisorsQ(_P, &V, E, &Flag, out);
    fflush(out);
  }
  return 0;
}