#include <palp/Global.h>
#include <palp/Rat.h>
#include <palp/Subpoly.h>

#include <memory>
#include <string>
#include <vector>
/*  NB mod 2^32, works for #poly<2^32   */

/*  #include <types.h>  ->  defines  _ILP32   (32-bit programs)
         [ -> isa_defs.h ]     or    _LP64    (long and pointer 64 bits)
    #include <stdint.h>	    (on Linux systems?)		... uintptr_t	*/

/*  #include <limits.h>     LONG_MAX = 2147483647 vs. 2^63-1 	  */

/*   L->dbname points at constant string => don't change, use:
 *   dbname, which has "File_Ext_NCmax" extra characters allocated !!!     */

/*   uchar=unsigned char        uint=unsigned int
 *   NP=#Polys  NB=#Bytes  #files=#nv's with NP>0  #lists=#(nv,nuc) with NP>0
 *
 * uchar  rd  k_1 ... k_rd                         // ... not in data base !!
 *
 * uchar  dim  #files  nv_max  nuc_max
 * uint                 #lists  hNF  hSM  hNM  hNB  slNF  slSM  slNM  slNB
 * uchar  v1    #nuc's with v1
 * uchar  nuc1  uint #NF(v1,nuc1)  uchar nuc2  uint #NV(v1,nuc2)  ...
 * uchar  v2    #nuc's with v2
 * uchar  nuc1  uint #NF(v2,nuc1)  uchar nuc2  uint #NV(v2,nuc2)  ...
 * uchar  "all hNF honest nf's"
 * uchar  "all slNF sublattice {nv nuc nf[]}'s"
 */

void Small_Make_Dual(PolyPointList *_P, VertexNumList *_V, EqList *_E) {
  int i;
  EqList AE = *_E;
  VNL_to_DEL(_P, _V, _E);
  if (!EL_to_PPL(&AE, _P, &_P->n)) {
    fputs("Error: Small_Make_Dual could not convert equation list to "
          "polytope\n",
          stderr);
    exit(1);
  }
  _V->nv = _P->np;
  for (i = 0; i < _P->np; i++)
    _V->v[i] = i;
}

void Polyi_2_DBo(char *polyi, char *dbo) {
  std::string dbnames(dbo);
  char *fx;
  FILE *F = fopen(polyi, "rb"), *Finfo, *Fv, *Fsl;
  time_t Tstart = time(nullptr);
  FInfoList L;
  UPint tNF = 0;
  int d, v, nu, i, j, list_num, sl_nNF, sl_SM, sl_NM, sl_NB;
  Along tNB = 0;

  if (!*polyi) {
    puts("With -do you require -pi or -di and -pa");
    exit(1);
  }

  if (F == nullptr) {
    printf("Input file %s not found\n", polyi);
    exit(1);
  }
  dbnames.resize(dbnames.size() + File_Ext_NCmax + 1, '\0');
  dbnames[strlen(dbo)] = '.';
  fx = dbnames.data() + strlen(dbo) + 1;
  fx[0] = '\0';
  dbnames += "info";
  Finfo = fopen(dbnames.c_str(), "w");
  if (Finfo == nullptr) {
    fprintf(stderr, "Error: Polyi_2_DBo cannot create %s.info\n", dbo);
    exit(1);
  }
  printf("Read %s (", polyi);
  fflush(stdout);

  Init_FInfoList(&L); /* start reading the file */
  d = fgetc(F);
  if (d != 0) { /* for(i=0;i<d;i++) fgetc(F); */
    printf("Recursion depth %d forbidden in DB !!\n", d);
    Finfo = stdout;
  }

  d = fgetc(F);
  L.nV = fgetc(F);
  L.nVmax = fgetc(F);
  L.NUCmax = fgetc(F);
  list_num = fgetUI(F);
  L.nNF = fgetUI(F);
  L.nSM = fgetUI(F);
  L.nNM = fgetUI(F);
  L.NB = fgetUI(F);
  sl_nNF = fgetUI(F);
  sl_SM = fgetUI(F);
  sl_NM = fgetUI(F);
  sl_NB = fgetUI(F);

  for (i = 0; i < L.nV; i++) {
    v = fgetc(F);
    L.nNUC[v] = fgetc(F); /* read #nuc's per #Vert */
    for (j = 0; j < L.nNUC[v]; j++) {
      L.NFnum[v][nu = fgetc(F)] = fgetUI(F); /* read nuc and #NF(v,nu)*/
      tNF += L.NFnum[v][nu];
      tNB += nu * L.NFnum[v][nu];
    }
  }
  if (tNF != L.nNF) {
    fprintf(stderr, "Error: Polyi_2_DBo total NF mismatch: tNF=%u L.nNF=%u\n",
            (unsigned)tNF, (unsigned)L.nNF);
    exit(1);
  }
  if ((unsigned int)(tNB - L.NB) != 0) {
    fprintf(stderr, "Error: Polyi_2_DBo byte count mismatch\n");
    exit(1);
  }
  L.NB = tNB;

  printf("%lldpoly +%dsl %lldb)  write %s.* (%d files)  ",
         2 * L.nNF - L.nSM - L.nNM, 2 * sl_nNF - sl_SM - sl_NM, L.NB + sl_NB,
         dbo, L.nV + 1 + (sl_nNF > 0));
  fflush(stdout);

  fprintf(Finfo, "%d  %d %d %d  %d  %lld %d %lld %lld  %d %d %d %d\n\n", d,
          L.nV, L.nVmax, L.NUCmax, list_num, L.nNF, L.nSM, L.nNM, L.NB, sl_nNF,
          sl_SM, sl_NM, sl_NB);

  for (v = d + 1; v <= L.nVmax; v++)
    if (L.nNUC[v]) /* honest info */
    {
      i = 0;
      fprintf(Finfo, "%d %d\n", v, L.nNUC[v]); /*  v  #nuc's  */
      for (nu = 1; nu <= L.NUCmax; nu++)
        if (L.NFnum[v][nu]) {
          fprintf(Finfo, "%d %d%s", nu, L.NFnum[v][nu], /* nuc #NF(v,nuc) */
                  (++i < L.nNUC[v]) ? "  " : "\n");
        }
    }
  if (Finfo == stdout)
    exit(1);
  if (ferror(Finfo)) {
    printf("File error in %s\n", dbnames.c_str());
    exit(1);
  }
  fclose(Finfo);
  fflush(stdout);

  for (v = d + 1; v <= L.nVmax; v++)
    if (L.nNUC[v]) /* write  honest polys */
    {
      char ext[4] = {'v', 0, 0, 0};
      ext[1] = '0' + v / 10;
      ext[2] = '0' + v % 10;
      strcpy(fx, ext);
      Fv = fopen(dbnames.c_str(), "wb");
      if (Fv == nullptr) {
        fprintf(stderr, "Error: Polyi_2_DBo cannot create %s\n",
                dbnames.c_str());
        exit(1);
      }

      for (nu = 1; nu <= L.NUCmax; nu++)
        if (L.NFnum[v][nu]) {
          int vnuNB = nu * L.NFnum[v][nu];
          for (i = 0; i < vnuNB; i++)
            fputc(fgetc(F), Fv);
        }
      if (ferror(Fv)) {
        printf("File error in %s\n", dbnames.c_str());
        exit(1);
      }
      fclose(Fv);
    }

  if (sl_nNF) /* write  sublattice polys */
  {
    strcpy(fx, "sl");
    Fsl = fopen(dbnames.c_str(), "wb");
    if (Fsl == nullptr) {
      fprintf(stderr, "Error: Polyi_2_DBo cannot create %s.sl\n", dbo);
      exit(1);
    }
    for (i = 0; i < sl_NB; i++)
      fputc(fgetc(F), Fsl);
    if (ferror(Fsl)) {
      printf("File error in %s\n", dbnames.c_str());
      exit(1);
    }
    fclose(Fsl);
  }

  printf("done (%ds)\n", (int)difftime(time(nullptr), Tstart));

  if (ferror(F)) {
    printf("File error in %s\n", polyi);
    exit(1);
  }
  fclose(F);
}

void Init_DB(NF_List *_NFL) {
  /* Read the database, create RAM_poly;
     for given, nv, nuc the matching is as follows:
     DB:  |0|1|...|B-1|B|...|2B|...|((n-1)/B)*B|...|n-2|n-1|
     RAM:             |0  |   1|...| (n-1)/B)-1|
     (n...nNF[nv][nuc], B...BLOCK_LENGTH, the offsets are DB->Fv_pos[v][nuc]
     and DB->RAM_pos[v][nuc], respectively;
     each entry |x| corresponds to nuc unsigned characters)  */

  time_t Tstart = time(nullptr);
  std::string dbname(_NFL->dbname);
  char *fx;
  DataBase *DB = &_NFL->DB;
  int d, v, nu, i, j, list_num, sl_nNF, sl_SM, sl_NM, sl_NB, RAM_pos = 0;
  Along RAM_size = 0;

  printf("Reading data-base %s: ", dbname.c_str());
  dbname.resize(dbname.size() + File_Ext_NCmax + 1, '\0');
  dbname[strlen(_NFL->dbname)] = '\0';
  fx = dbname.data() + strlen(_NFL->dbname) + 1;
  dbname[strlen(_NFL->dbname)] = '.';
  fx[0] = '\0';
  dbname += "info";

  /* read the info-file: */
  DB->Finfo = fopen(dbname.c_str(), "r");
  if (DB->Finfo == nullptr) {
    fprintf(stderr, "Error: Open_DB cannot open %s\n", dbname.c_str());
    exit(1);
  }
  if (fscanf(DB->Finfo, "%d  %d %d %d  %d  %lld %d %lld %lld  %d %d %d %d", &d,
             &DB->nV, &DB->nVmax, &DB->NUCmax, &list_num, &DB->nNF, &DB->nSM,
             &DB->nNM, &DB->NB, &sl_nNF, &sl_SM, &sl_NM, &sl_NB) != 13) {
    fputs("Error: Open_DB malformed info file\n", stderr);
    exit(1);
  }
  printf("%lld+%dsl %lldnf %lldb", 2 * (DB->nNF) - DB->nSM - DB->nNM,
         2 * sl_nNF - sl_SM - sl_NM, DB->nNF + sl_nNF, DB->NB + sl_NB);
  /* if( _FILE_OFFSET_BITS < 64 ) assert(DB->NB <= LONG_MAX);	Along */
  if (_NFL->d && (d != _NFL->d)) {
    fprintf(stderr, "Error: Open_DB dimension mismatch %d != %d\n", d, _NFL->d);
    exit(1);
  } else
    _NFL->d = d;

  for (v = 1; v < VERT_Nmax; v++) {
    DB->nNUC[v] = 0;
    for (nu = 0; nu < NUC_Nmax; nu++)
      DB->NFnum[v][nu] = 0;
  }

  for (i = 0; i < DB->nV; i++) {
    if (fscanf(DB->Finfo, "%d", &v) != 1) {
      fputs("Error: Open_DB expected vertex count\n", stderr);
      exit(1);
    }
    if (fscanf(DB->Finfo, "%d", &(DB->nNUC[v])) != 1) {
      fprintf(stderr, "Error: Open_DB expected nuc count for v=%d\n", v);
      exit(1);
    }
    for (j = 0; j < DB->nNUC[v]; j++) {
      if (fscanf(DB->Finfo, "%d", &nu) != 1) {
        fprintf(stderr, "Error: Open_DB expected nu value for v=%d\n", v);
        exit(1);
      }
      if (fscanf(DB->Finfo, "%d", &(DB->NFnum[v][nu])) != 1) {
        fprintf(stderr, "Error: Open_DB expected NF count for v=%d nu=%d\n", v,
                nu);
        exit(1);
      }
      RAM_size += nu * ((DB->NFnum[v][nu] - 1) / BLOCK_LENGTH);
    }
  }

  if (ferror(DB->Finfo)) {
    printf("File error in %s\n", dbname.c_str());
    exit(1);
  }
  fclose(DB->Finfo);
  fflush(stdout);
  if (RAM_size > INT_MAX) {
    fprintf(stderr, "Error: Open_DB RAM_size %lld exceeds INT_MAX\n",
            (long long)RAM_size);
    exit(1);
  }

  DB->RAM_NF_owner =
      std::make_unique<unsigned char[]>(static_cast<size_t>(RAM_size));
  DB->RAM_NF = DB->RAM_NF_owner.get();

  /* read the DB-files and create RAM_NF: */
  for (v = 2; v <= DB->nVmax; v++)
    if (DB->nNUC[v]) {
      char ext[4] = {'v', 0, 0, 0};
      ext[1] = '0' + v / 10;
      ext[2] = '0' + v % 10;
      strcpy(fx, ext);
      DB->Fv[v] = fopen(dbname.c_str(), "rb");
      if (DB->Fv[v] == nullptr) {
        fprintf(stderr, "Error: Open_DB cannot open %s\n", dbname.c_str());
        exit(1);
      }
      FSEEK(DB->Fv[v], 0, SEEK_END);

      DB->Fv_pos[v][0] = 0;
      for (nu = 0; nu <= DB->NUCmax; nu++) {
        DB->RAM_pos[v][nu] = RAM_pos;
        for (i = 0; i < (DB->NFnum[v][nu] - 1) / BLOCK_LENGTH; i++) {
          FSEEK(DB->Fv[v], DB->Fv_pos[v][nu] + (i + 1) * BLOCK_LENGTH * nu,
                SEEK_SET);
          for (j = 0; j < nu; j++)
            DB->RAM_NF[RAM_pos++] = fgetc(DB->Fv[v]);
        }
        DB->Fv_pos[v][nu + 1] = DB->Fv_pos[v][nu] + DB->NFnum[v][nu] * nu;
      }
    }

  printf("  done (%ds)\n", (int)difftime(time(nullptr), Tstart));
  fflush(stdout);
}

char Compare_Poly(int *nuc, unsigned char *uc1, unsigned char *uc2) {
  /* uc1 always encodes a single poly,
     uc2 might encode a single poly or a mirror pair;
     returns: 'm': if uc1 isn't in uc2, but its mirror is;
              'i': if uc1 is in uc2;
              'l': if uc1 is less than uc2;
              'g': if uc1 is greater than uc2; */
  switch (RIGHTminusLEFT(uc1, uc2, nuc)) {
  case -1:
    return 'g';
  case 1:
    return 'l';
  case 0: {
    if ((*uc1 + *uc2) % 4 == 3)
      return 'm';
    else
      return 'i';
  }
  default:
    puts("Sth. wrong in Compare_Poly!!!");
    exit(1);
  }
  return 0;
}

int Is_in_DB(int *nv, int *nuc, unsigned char *uc, NF_List *_NFL) {
  /* Uses the following basic strategy for identifying the position of some
     object x w.r.t. entries of a list l of length n:
     min_pos=-1;
     max_pos=n;
     while (max_pos-min_pos>1){
       int pos=(max_pos+min_pos)/2;
       switch (Compare(x,l[pos]){
         case 'equal': return sth.;
         case 'less': {max_pos=pos; continue:}
         case 'greater': min_pos=pos;}}
     results in a return or l[min_pos] < x < l[max_pos=min_pos+1]
     Applied here to locate x=uc=encoded polyhedron first w.r.t. RAM_poly and
     then w.r.t. the location in the database */

  int pos, min_pos = -1, max_RAM_pos, i, max_Fv_piece;
  Along Fv_pos;
  unsigned char Aux_poly[BLOCK_LENGTH * NUC_Nmax];

  DataBase *DB = &_NFL->DB;
  if (!DB->NFnum[*nv][*nuc])
    return 0;

  /* Analyse the position of uc w.r.t. DB->RAM_NF */
  max_RAM_pos = (DB->NFnum[*nv][*nuc] - 1) / BLOCK_LENGTH;
  while (max_RAM_pos - min_pos > 1) {
    pos = (max_RAM_pos + min_pos) / 2;
    switch (Compare_Poly(
        nuc, uc, &(DB->RAM_NF[DB->RAM_pos[*nv][*nuc] + (*nuc) * pos]))) {
    case 'm':
      return 0;
    case 'i':
      return 1;
    case 'l': {
      max_RAM_pos = pos;
      continue;
    }
    case 'g':
      min_pos = pos;
    }
  }

  /* Look for uc in the database: */ Fv_pos = max_RAM_pos;
  if (Fv_pos == (DB->NFnum[*nv][*nuc] - 1) / BLOCK_LENGTH)
    max_Fv_piece = DB->NFnum[*nv][*nuc] - Fv_pos * BLOCK_LENGTH;
  else
    max_Fv_piece = BLOCK_LENGTH;
  min_pos = -1;
  if (FSEEK(DB->Fv[*nv],
            DB->Fv_pos[*nv][*nuc] + Fv_pos * (Along)((*nuc) * BLOCK_LENGTH),
            SEEK_SET)) {
    printf("Error in fseek in Is_in_DB!");
    exit(1);
  }
  for (i = 0; i < (*nuc) * (max_Fv_piece); i++)
    Aux_poly[i] = fgetc(DB->Fv[*nv]);
  while (max_Fv_piece - min_pos > 1) {
    pos = (max_Fv_piece + min_pos) / 2;
    switch (Compare_Poly(nuc, uc, &(Aux_poly[(*nuc) * pos]))) {
    case 'm': {
      return 0;
    }
    case 'i': {
      return 1;
    }
    case 'l': {
      max_Fv_piece = pos;
      continue;
    }
    case 'g': {
      min_pos = pos;
    }
    }
  }
  return 0;
}

void Add_Polya_2_DBi(char *dbi, char *polya, char *dbo, FILE *out) {
  FInfoList FIi, FIa, FIo;
  Along Apos, HIpos, HApos, Inp, tnb = 0, tNF = 0;
  unsigned char ucI[NUC_Nmax], ucA[NUC_Nmax], *ucSL = nullptr, *uc;
  int SLp[SL_Nmax];
  int d, vI, nuI, IslNF, IslSM, IslNM, nu;
  unsigned Ili, u, Oli = 0, IslNB;
  int v, vA, nuA, AslNF, AslSM, AslNM, i;
  unsigned Ali, a;
  Along AslNB;
  int s, slNF = 0, slSM = 0, slNM = 0, slNB = 0, slNP = 0, j;
  UPint Anp;
  int AmI = 00, ms, newout = strcmp(dbi, dbo) && (*dbo);
  std::string Ifn = dbi;
  std::string Ofn = newout ? dbo : dbi;
  FILE *FI, *FA, *FO;
  if (*polya == 0) {
    puts("-pa file required");
    exit(1);
  }
  Init_FInfoList(&FIi);
  Init_FInfoList(&FIa);
  Ifn += ".info";
  if (*dbo == 0)
    dbo = dbi;

  if (nullptr == (FI = fopen(Ifn.c_str(), "r"))) {
    printf("Cannot open %s", Ifn.c_str());
    exit(1);
  }
  if (nullptr == (FA = fopen(polya, "rb"))) {
    printf("Cannot open %s", polya);
    exit(1);
  }
  if (fscanf(FI, "%d%d%d%d%d%lld%d%lld %lld %d%d%d%d", &d, &i, &j, &nu, &Ili,
             &FIi.nNF, &FIi.nSM, &FIi.nNM, &FIi.NB, &IslNF, &IslSM, &IslNM,
             &IslNB) != 13) {
    fputs("Error: Add_Polya_2_DBi malformed info header\n", stderr);
    exit(1);
  }
  FIi.nV = i;
  FIi.nVmax = j;
  FIi.NUCmax = nu;
  /*	printf("%d %d %d %d %d %d %d %d %lld %d %d %d %d\n",
     d,FIi.nV,FIi.nVmax,FIi.NUCmax,Ili,FIi.nNF,FIi.nSM,
     FIi.nNM,FIi.NB,IslNF,IslSM,IslNM,IslNB); */
  for (i = 0; i < FIi.nV; i++) {
    if (fscanf(FI, "%d", &v) != 1) {
      fputs("Error: Add_Polya_2_DBi expected vertex count\n", stderr);
      exit(1);
    }
    if (fscanf(FI, "%d", &j) != 1) {
      fprintf(stderr, "Error: Add_Polya_2_DBi expected nuc count for v=%d\n",
              v);
      exit(1);
    }
    FIi.nNUC[v] = j;
    for (j = 0; j < FIi.nNUC[v]; j++) {
      if (fscanf(FI, "%d", &nu) != 1) {
        fprintf(stderr, "Error: Add_Polya_2_DBi expected nu value for v=%d\n",
                v);
        exit(1);
      }
      if (fscanf(FI, "%d", &FIi.NFnum[v][nu]) != 1) {
        fprintf(stderr,
                "Error: Add_Polya_2_DBi expected NF count for v=%d nu=%d\n", v,
                nu);
        exit(1);
      }
    }
  }
  if (tNF != FIi.nNF) {
    fprintf(stderr,
            "Error: Add_Aux_to_DB total NF mismatch: tNF=%lld FIi.nNF=%lld\n",
            (long long)tNF, (long long)FIi.nNF);
    exit(1);
  }
  tNF = 0;
  {
    int rd_byte = fgetc(FA);
    if (rd_byte != EOF) {
      ungetc(rd_byte, FA);
      fputs("Error: Add_Aux_to_DB aux-file recursion depth must be 0\n",
            stderr);
      exit(1);
    }
  }
  Read_Bin_Info(FA, &j, &Ali, &AslNF, &AslSM, &AslNM, &AslNB, &FIa);
  if (d != j) {
    fprintf(stderr, "Error: Add_Aux_to_DB dimension mismatch %d != %d\n", d, j);
    exit(1);
  }
  if (IslNF) {
    Ifn.replace(Ifn.size() - 4, 4, ".sl");
    FI = fopen(Ifn.c_str(), "rb");
    if (FI == nullptr) {
      fprintf(stderr, "Error: Add_Aux_to_DB cannot open %s.sl\n", Ifn.c_str());
      exit(1);
    }
  }
  std::vector<unsigned char> ucSL_buffer;
  if ((IslNB + AslNB))
    ucSL_buffer.resize(IslNB + AslNB);
  ucSL = ucSL_buffer.data();
  HApos = FTELL(FA);
  FSEEK(FA, 0, SEEK_END);
  Apos = FTELL(FA);
  Inp = 2 * FIi.nNF - FIi.nSM - FIi.nNM;
  Anp = 2 * FIa.nNF - FIa.nSM - FIa.nNM;
  printf("Data on %s:  %lld+%dsl  %lldb  (%dd)\n", dbi, Inp,
         /* Islp= */ 2 * IslNF - IslNM - IslSM, FIi.NB + slNB, d);
  printf("Data on %s:  %u+%dsl  %lldb  (%dd)\n", polya, Anp,
         /* Aslp= */ 2 * AslNF - AslNM - AslSM, Apos, d);
  if (HApos + FIa.NB + AslNB != Apos) {
    fputs("Error: Add_Aux_to_DB aux-file size mismatch\n", stderr);
    exit(1);
  }
  FSEEK(FA, -AslNB, SEEK_CUR);
  s = 0;
  if (s < AslNF)
    AuxGet_vn_uc(FA, &vA, &nuA, ucA);
  for (i = 0; i < IslNF; i++) {
    AuxGet_vn_uc(FI, &vI, &nuI, ucI);
    uc = &ucSL[2 + (SLp[slNF++] = slNB)];
    while (s < AslNF) {
      if (!(AmI = vA - vI))
        if (!(AmI = nuA - nuI))
          AmI = RIGHTminusLEFT(ucI, ucA, &nuI);
      if (AmI < 0) /* put A */
      {
        uc[-2] = vA;
        uc[-1] = nuA;
        for (j = 0; j < nuA; j++)
          uc[j] = ucA[j];
        slNB += 2 + nuA;
        if ((ms = (*uc % 4))) {
          if (ms < 3)
            slNM++;
        } else
          slSM++;
        if ((++s) < AslNF)
          AuxGet_vn_uc(FA, &vA, &nuA, ucA);
        uc = &ucSL[2 + (SLp[slNF++] = slNB)];
      } else
        break;
    }
    if ((s < AslNF) && (AmI == 0)) /* put I==A */
    {
      uc[-2] = vI;
      uc[-1] = nuI;
      for (j = 0; j < nuI; j++)
        uc[j] = ucI[j];
      if ((*ucI % 4) != (*ucA % 4))
        *uc = 3 + 4 * (*uc / 4);
      if ((ms = (*uc % 4))) {
        if (ms < 3)
          slNM++;
      } else
        slSM++;
      if ((++s) < AslNF)
        AuxGet_vn_uc(FA, &vA, &nuA, ucA);
      tnb += 2 + nuI;
    } else /* put I */
    {
      uc[-2] = vI;
      uc[-1] = nuI;
      for (j = 0; j < nuI; j++)
        uc[j] = ucI[j];
      if ((ms = (*uc % 4))) {
        if (ms < 3)
          slNM++;
      } else
        slSM++;
    }
    slNB += 2 + nuI;
  }
  while (s < AslNF) /* put A */
  {
    uc = &ucSL[2 + (SLp[slNF++] = slNB)];
    slNB += 2 + nuA;
    uc[-2] = vA;
    uc[-1] = nuA;
    for (j = 0; j < nuA; j++)
      uc[j] = ucA[j];
    if ((ms = (*uc % 4))) {
      if (ms < 3)
        slNM++;
    } else
      slSM++;
    if ((++s) < AslNF)
      AuxGet_vn_uc(FA, &vA, &nuA, ucA);
  }
  if (tnb + slNB != IslNB + AslNB) {
    fputs("Error: Add_Aux_to_DB sublattice byte count mismatch\n", stderr);
    exit(1);
  } /* SL done */

  printf("SL: %dnf %dsm %dnm %db -> ", slNF, slSM, slNM, slNB);
  if (IslNF) {
    HIpos = FTELL(FI);
    if (HIpos != IslNB) {
      fprintf(stderr, "Error: Add_Aux_to_DB sublattice read size mismatch\n");
      exit(1);
    }
    if (ferror(FI)) {
      fprintf(stderr, "Error: Add_Aux_to_DB sublattice file read error\n");
      exit(1);
    }
    fclose(FI);
    if (!newout)
      remove(Ifn.c_str());
  } /* SL file done */

  FSEEK(FA, HApos, SEEK_SET);
  Init_FInfoList(&FIo);
  FIo.nVmax = palp::max(FIi.nVmax, FIa.nVmax);
  FIo.NUCmax = palp::max(FIi.NUCmax, FIa.NUCmax); /* Tnb=0; */
  for (v = d + 1; v <= FIo.nVmax; v++)
    for (nu = 1; nu <= FIo.NUCmax; nu++)
      if ((FIo.NFnum[v][nu] = FIi.NFnum[v][nu] + FIa.NFnum[v][nu])) {
        FIo.nNUC[v]++;
        Oli++; /* Tnb+=FIo.NFnum[v][nu]; */
      }
  FIo.nV = 0;
  for (v = d + 1; v <= FIo.nVmax; v++)
    if (FIo.nNUC[v])
      FIo.nV++;

  for (v = d + 1; v <= FIo.nVmax; v++)
    if (FIo.nNUC[v]) {
      std::string vxt = ".v";
      vxt += static_cast<char>(v / 10 + '0');
      vxt += static_cast<char>(v % 10 + '0');
      std::string Ifn_v = Ifn;
      std::string Ofn_v = Ofn;
      Ifn_v.replace(Ifn_v.size() - 4, 4, vxt);
      Ofn_v.replace(Ofn_v.size() - 4, 4, vxt);
      if (FIi.nNUC[v]) {
        if (!newout) {
          std::string backup = Ifn_v;
          backup += SAVE_FILE_EXT;
          if (rename(Ofn_v.c_str(), backup.c_str()) != 0) {
            fprintf(stderr, "Error: Add_Aux_to_DB rename %s -> %s failed\n",
                    Ofn_v.c_str(), backup.c_str());
            exit(1);
          }
        }
        if (nullptr == (FI = fopen(Ifn_v.c_str(), "rb"))) {
          printf("Ifn %s failed", Ifn_v.c_str());
          exit(1);
        }
      }
      if (nullptr == (FO = fopen(Ofn_v.c_str(), "wb"))) {
        printf("Ofn %s failed", Ofn_v.c_str());
        exit(1);
      }

      for (nu = 1; nu <= FIo.NUCmax; nu++)
        if (FIo.NFnum[v][nu]) {
          unsigned int I_NF = FIi.NFnum[v][nu], A_NF = FIa.NFnum[v][nu],
                       O_NF = 0;
          UPint neq = 0, peq = 0, pa = 0;
          Along pi = 0, po = 0;

          a = 0;
          if (0 < A_NF) {
            AuxGet_uc(FA, &nu, ucA);
            pa += 1 + (((*ucA) % 4) == 3);
          }
          for (u = 0; u < I_NF; u++) {
            AuxGet_uc(FI, &nu, ucI);
            pi += 1 + (((*ucI) % 4) == 3);
            while (a < A_NF) {
              AmI = RIGHTminusLEFT(ucI, ucA, &nu);
              if (AmI < 0) /* put A */
              {
                po += 1 + (((*ucA) % 4) == 3);
                AuxPut_hNF(FO, &v, &nu, ucA, &FIo, &slNF, &slSM, &slNM, &slNB,
                           ucSL, SLp);
                O_NF++;
                a++;
                if (a < A_NF) {
                  AuxGet_uc(FA, &nu, ucA);
                  pa += 1 + (((*ucA) % 4) == 3);
                }
              } else
                break;
            }
            if ((a < A_NF) && (AmI == 0)) /* put I==A */
            {
              int mm = 10 * (*ucI % 4) + (*ucA % 4);
              switch (mm) {
              case 00:
                peq++;
                break;
              case 11:
                peq++;
                break;
              case 12:
                *ucI += 2;
                break;
              case 13:
                *ucI += 2;
                peq++;
                break;
              case 21:
                *ucI += 1;
                break;
              case 22:
                peq++;
                break;
              case 23:
                *ucI += 1;
                peq++;
                break;
              case 31:
                peq++;
                break;
              case 32:
                peq++;
                break;
              case 33:
                peq += 2;
                break;
              default:
                puts("inconsistens mirror flags");
                exit(1);
              }
              AuxPut_hNF(FO, &v, &nu, ucI, &FIo, &slNF, &slSM, &slNM, &slNB,
                         ucSL, SLp);
              po += 1 + (((*ucI) % 4) == 3);
              neq++;
              a++;
              if (a < A_NF) {
                AuxGet_uc(FA, &nu, ucA);
                pa += 1 + (((*ucA) % 4) == 3);
              }
            } else {
              AuxPut_hNF(FO, &v, &nu, ucI, &FIo, &slNF, &slSM, &slNM, &slNB,
                         ucSL, SLp);
              po += 1 + (((*ucI) % 4) == 3);
            }
            O_NF++;
          }
          while (a < A_NF) /* put A */
          {
            O_NF++;
            po += 1 + (((*ucA) % 4) == 3);
            AuxPut_hNF(FO, &v, &nu, ucA, &FIo, &slNF, &slSM, &slNM, &slNB, ucSL,
                       SLp);
            ++a;
            if (a < A_NF) {
              AuxGet_uc(FA, &nu, ucA);
              pa += 1 + (((*ucA) % 4) == 3);
            }
          }
          if (pi + pa != peq + po) {
            fprintf(stderr,
                    "Error: Add_Aux_to_DB checksum mismatch for v=%d nu=%d\n",
                    v, nu);
            exit(1);
          } /* checksum(v,nu) */
          FIo.NFnum[v][nu] = O_NF;
          FIo.nNF += O_NF;
          FIo.NB += O_NF * nu;
          if (O_NF + neq != I_NF + A_NF) {
            fprintf(stderr,
                    "Error: Add_Aux_to_DB NF merge mismatch for v=%d nu=%d\n",
                    v, nu);
            exit(1);
          }
          /*
          {static int list;printf("#%d v=%d nu=%d Inf=%d Anf=%d ",++list,v,nu,
          I_NF,A_NF);printf("Onf=%d   pi=%d pa=%d  po=%d\n",O_NF,pi,pa,po);
          }*/
        }
      if (FIi.nNUC[v]) {
        if (ferror(FI)) {
          fprintf(stderr, "Error: Add_Aux_to_DB input file read error\n");
          exit(1);
        }
        fclose(FI);
        if (!newout)
          remove(Ifn_v.c_str());
      }
      if (ferror(FO)) {
        fprintf(stderr, "Error: Add_Aux_to_DB output file write error\n");
        exit(1);
      }
      fclose(FO);
    }
  tnb = 0;

  if (slNF) {
    std::string Ofn_sl = Ofn;
    Ofn_sl.replace(Ofn_sl.size() - 4, 4, ".sl");
    FO = fopen(Ofn_sl.c_str(), "wb");
    if (FO == nullptr) {
      fprintf(stderr, "Error: Add_Aux_to_DB cannot create %s.sl\n", Ofn_sl.c_str());
      exit(1);
    }
    for (i = 0; i < slNF; i++) /* write SL */
    {
      uc = &ucSL[SLp[i] + 2];
      v = uc[-2];
      nu = uc[-1];
      tnb += nu + 2;
      if (uc[-2] >= VERT_Nmax) {
        fprintf(stderr, "Error: Add_Aux_to_DB vertex count %d out of range\n",
                uc[-2]);
        exit(1);
      }
      fputc(uc[-2], FO);
      slNP += 1 + (((*uc) % 4) == 3);
      fputc(nu, FO);
      for (s = 0; s < nu; s++)
        fputc(uc[s], FO);
    }
    if (tnb != slNB) {
      fprintf(stderr, "Error: Add_Aux_to_DB sublattice byte write mismatch\n");
      exit(1);
    }
    if (slNP != 2 * slNF - slNM - slSM) {
      fprintf(stderr, "Error: Add_Aux_to_DB sublattice NP count mismatch\n");
      exit(1);
    }
    if (ferror(FO)) {
      fprintf(stderr, "Error: Add_Aux_to_DB sublattice file write error\n");
      exit(1);
    }
    fclose(FO);
  }
  printf("\nd=%d v%d v<=%d n<=%d vn%d  %lld %d %lld %lld  %d %d %d %d\n", d,
         FIo.nV, FIo.nVmax, FIo.NUCmax, Oli, FIo.nNF, FIo.nSM, FIo.nNM, FIo.NB,
         slNF, slSM, slNM, slNB);

  std::string Ofn_info = Ofn;
  Ofn_info.replace(Ofn_info.size() - 4, 4, ".info");
  FO = fopen(Ofn_info.c_str(), "w");
  if (FO == nullptr) {
    fprintf(stderr, "Error: Add_Aux_to_DB cannot create %s.info\n", Ofn_info.c_str());
    exit(1);
  }
  fprintf(FO, /* write FO.info */
          "%d  %d %d %d  %d  %lld %d %lld %lld  %d %d %d %d\n\n", d, FIo.nV,
          FIo.nVmax, FIo.NUCmax, Oli, FIo.nNF, FIo.nSM, FIo.nNM, FIo.NB, slNF,
          slSM, slNM, slNB);

  for (v = d + 1; v <= FIo.nVmax; v++)
    if (FIo.nNUC[v]) /* honest info */
    {
      i = 0;
      fprintf(FO, "%d %d\n", v, FIo.nNUC[v]); /*  v  #nuc's  */
      for (nu = 1; nu <= FIo.NUCmax; nu++)
        if (FIo.NFnum[v][nu]) {
          fprintf(FO, "%d %d%s", nu, FIo.NFnum[v][nu], /* nuc #NF(v,nuc) */
                  (++i < FIo.nNUC[v]) ? "  " : "\n");
        }
    }
  printf("Writing %s: %lld+%dsl %lldm+%ds %lldb", dbo,
         2 * FIo.nNF - FIo.nNM - FIo.nSM, slNP,
         /*Tnb=*/FIo.nNF - FIo.nNM - FIo.nSM, FIo.nSM, FIo.NB + slNB);
  /*   if(tnb>99)
       {	long long tnp=(2*FIo.nNF-FIo.nNM-FIo.nSM)/10; tnp*=tnp;
          tnp/=20; tnp/=Tnb; printf("   [p^2/2m=%ldk]",tnp);
       }
  */
  Print_Expect(&FIo, out);
  puts("");
  if (ferror(FA)) {
    fprintf(stderr, "Error: Add_Aux_to_DB aux-file read error\n");
    exit(1);
  }
  fclose(FA);
  if (ferror(FO)) {
    fprintf(stderr, "Error: Add_Aux_to_DB info-file write error\n");
    exit(1);
  }
  fclose(FO);
}
namespace {
struct OrderState {
  int initialized = 0;
  int V = 0, NU = 0;
  unsigned char UC[NUC_Nmax] = {};
};
OrderState slOrderState, hnfOrderState;

int checkOrder(OrderState &state, int *v, int *nu, unsigned char *uc) {
  if (state.initialized) {
    if (state.V > (*v))
      return 0;
    if ((state.V == (*v)) && (state.NU > (*nu)))
      return 0;
    if ((state.V == (*v)) && (state.NU == (*nu)) &&
        (RIGHTminusLEFT(state.UC, uc, &state.NU) <= 0))
      return 0;
  }
  state.V = *v;
  state.NU = *nu;
  for (int n = 0; n < state.NU; n++)
    state.UC[n] = uc[n];
  state.initialized = 1;
  return 1;
}
} // namespace

int Check_sl_order(int *v, int *nu, unsigned char *uc) {
  return checkOrder(slOrderState, v, nu, uc);
}
int Check_hnf_order(int *v, int *nu, unsigned char *uc) {
  return checkOrder(hnfOrderState, v, nu, uc);
}
void Print_NF(FILE *F, int *d, int *v, Long NF[POLY_Dmax][VERT_Nmax]) {
  int i, j;
  fprintf(F, "%d %d\n", *d, *v);
  for (i = 0; i < *d; i++)
    for (j = 0; j < *v; j++)
      fprintf(F, "%d%s", (int)NF[i][j], (*v == j + 1) ? "\n" : " ");
}
void Print_Missing_Mirror(int *d, int *v, int *nu, unsigned char *uc,
                          PolyPointList *_P, FILE *out) {
  int I, J, MS;
  Long NF[POLY_Dmax][VERT_Nmax];
  VertexNumList V;
  EqList E;
  UCnf2vNF(d, v, nu, uc, NF, &MS);
  MS %= 4;
  _P->n = *d;
  _P->np = *v;
  for (I = 0; I < *v; I++)
    for (J = 0; J < *d; J++)
      _P->x[I][J] = NF[J][I];
  if (MS == 2)
    Print_NF(out, d, v, NF);
  else if (MS == 1) {
    IP_Check(_P, &V, &E);
    Small_Make_Dual(_P, &V, &E);
    Make_Poly_NF(_P, &V, &E, NF);
    Print_NF(out, d, &(V.nv), NF);
  } else {
    puts("Only use Print_Missing_Mirror for MM!");
    exit(1);
  }
}

/*   cF->{1::c 2::C (extended output)}  cF->{-1::M (missing mirrors)}   */
void Check_NF_Order(char *polyi, char *dbi, int cF, PolyPointList *_P,
                    FILE *out) /* 1=MM */ {
  FILE *F = nullptr;
  FInfoList L; /* time_t Tstart=time(nullptr); */
  unsigned int rd, i, j, list_num, tln = 0, tSM = 0;
  int d, nu, v, si;
  int sl_nNF, sl_SM, sl_NM, sl_NB;
  Along tNF = 0, tNB = 0, tNM = 0, SLpos, Hpos = 0;
  std::string Ifn;
  if ((*polyi) && (*dbi))
    puts("only give one of -pi FILE or -di FILE");
  if ((*polyi == 0) && (*dbi == 0))
    puts("I need one of: -pi FILE or -di FILE");
  if ((*polyi == 0) + (*dbi == 0) != 1)
    exit(1);

  if (*polyi) /*   read INFO PART */
  {
    printf("Checking consistency of Aux/InFile %s\n", polyi);
    F = fopen(polyi, "rb");
    if (F == nullptr) {
      puts("File not found");
      exit(1);
    }
    Init_FInfoList(&L); /* start reading the file */

    rd = fgetc(F);
    if (rd > 127)
      rd = 128 * (rd - 128) + fgetc(F);
    if (rd > MAX_REC_DEPTH) {
      fprintf(stderr, "Error: DB_Check recursion depth %d exceeds maximum\n",
              rd);
      exit(1);
    }
    printf("rd=%d: ", rd);
    for (i = 0; i < rd; i++)
      printf(" %d", fgetc(F));
    if (rd)
      puts("");

    d = fgetc(F);
    L.nV = fgetc(F);
    L.nVmax = fgetc(F);
    L.NUCmax = fgetc(F);
    list_num = fgetUI(F);
    printf("d=%d  nV=%d nVmax=%d NUCmax=%d #lists=%d ", d, L.nV, L.nVmax,
           L.NUCmax, list_num);
    fflush(stdout);
    L.nNF = fgetUI(F);
    L.nSM = fgetUI(F);
    L.nNM = fgetUI(F);
    L.NB = fgetUI(F);
    sl_nNF = fgetUI(F);
    sl_SM = fgetUI(F);
    sl_NM = fgetUI(F);
    sl_NB = fgetUI(F);

    for (i = 0; i < L.nV; i++) {
      v = fgetc(F);
      L.nNUC[v] = fgetc(F);
      if (cF == 2)
        printf("v%dn%d: ", v, L.nNUC[v]);
      fflush(stdout);
      for (j = 0; j < L.nNUC[v]; j++) {
        L.NFnum[v][nu = fgetc(F)] = fgetUI(F); /* read nuc and #NF(v,nu)*/
        tNF += L.NFnum[v][nu];
        if (cF == 2)
          printf("%d:%d ", nu, L.NFnum[v][nu]);
        tln++;
        tNB += L.NFnum[v][nu] * (Along)nu;
      }
      if (cF == 2)
        puts("");
    }
    if ((unsigned int)(tNB - L.NB) != 0) {
      fputs("Error: DB_Check aux-file byte count mismatch\n", stderr);
      exit(1);
    }
    L.NB = tNB;
    if (cF == 2)
      printf("np=%lld+%dsl %lldb  %d files\n", 2 * L.nNF - L.nSM - L.nNM,
             2 * sl_nNF - sl_SM - sl_NM, L.NB + sl_NB, L.nV + 1 + (sl_nNF > 0));
    if (cF == 2)
      printf("%lldnf %dsm %lldnm %lldb   sl: %d %d %d %d\n", L.nNF, L.nSM,
             L.nNM, L.NB, sl_nNF, sl_SM, sl_NM, sl_NB);
    fflush(stdout);
  }
  if (*dbi) {
    printf("Checking consistency of DataBase %s:\n", dbi);
    Ifn = dbi;
    Ifn += ".info";
    F = fopen(Ifn.c_str(), "r");
    if (F == nullptr) {
      puts("Info File not found");
      exit(1);
    }
    Init_FInfoList(&L); /* start reading the file */
    if (fscanf(F, "%d%d%d%d%d%lld%d%lld %lld %d%d%d%d", &d, &i, &j, &nu,
               &list_num, &L.nNF, &L.nSM, &L.nNM, &L.NB, &sl_nNF, &sl_SM,
               &sl_NM, &sl_NB) != 13) {
      fputs("Error: DB_Check malformed info header\n", stderr);
      fclose(F);
      exit(1);
    }
    L.nV = i;
    L.nVmax = j;
    L.NUCmax = nu;
    printf("d=%d  nV=%d nVmax=%d NUCmax=%d #lists=%d  ", d, L.nV, L.nVmax,
           L.NUCmax, list_num);
    printf("np=%lld+%dsl %lldb  %d files\n", 2 * L.nNF - L.nSM - L.nNM,
           2 * sl_nNF - sl_SM - sl_NM, L.NB + sl_NB, L.nV + 1 + (sl_nNF > 0));
    printf("%lldnf %dsm %lldnm %lldb   sl: %d %d %d %d\n", L.nNF, L.nSM, L.nNM,
           L.NB, sl_nNF, sl_SM, sl_NM, sl_NB);
    fflush(stdout);

    for (i = 0; i < L.nV; i++) {
      if (fscanf(F, "%d", &v) != 1) {
        fputs("Error: DB_Check expected vertex count\n", stderr);
        fclose(F);
        exit(1);
      }
      if (fscanf(F, "%d", &j) != 1) {
        fprintf(stderr, "Error: DB_Check expected nuc count for v=%d\n", v);
        fclose(F);
        exit(1);
      }
      L.nNUC[v] = j;
      if (cF == 2)
        printf("v%dn%d: ", v, L.nNUC[v]);
      fflush(stdout);
      for (j = 0; j < L.nNUC[v]; j++) {
        if (fscanf(F, "%d", &nu) != 1) {
          fprintf(stderr, "Error: DB_Check expected nu value for v=%d\n", v);
          fclose(F);
          exit(1);
        }
        if (fscanf(F, "%d", &L.NFnum[v][nu]) != 1) {
          fprintf(stderr, "Error: DB_Check expected NF count for v=%d nu=%d\n",
                  v, nu);
          fclose(F);
          exit(1);
        }
        tln++;
        tNF += L.NFnum[v][nu];
        if (cF == 2)
          printf("%d:%d ", nu, L.NFnum[v][nu]);
      }
      if (cF == 2)
        puts("");
    }
    if (tln != list_num) {
      fprintf(stderr, "Error: DB_Check info list count mismatch\n");
      exit(1);
    }
  }
  printf("#hNF=%lld sum=%lld %s\n", L.nNF, tNF,
         (tNF == L.nNF) ? "o.k." : "Error");
  if (tNF != L.nNF)
    exit(1);
  if (ferror(F)) {
    fprintf(stderr, "Error: DB_Check file read error\n");
    exit(1);
  }
  tNF = 0;
  if (tln != list_num) {
    printf("ERROR: #li=%d != %d\n", list_num, tln);
    exit(1);
  }

  { /* long long np=2*L.nNF-L.nSM-L.nNM, pp2m=L.nNF-L.nSM-L.nNM;
    j=2*sl_nNF-sl_SM-sl_NM; */
    printf("np=%lld+%dsl  %lldnf %dsm %lldnm %lldb  sl: %d %d %d %d\n",
           2 * L.nNF - L.nSM - L.nNM, 2 * sl_nNF - sl_SM - sl_NM, L.nNF, L.nSM,
           L.nNM, L.NB, sl_nNF, sl_SM, sl_NM, sl_NB);
    fflush(stdout); /* if(pp2m) tln=(np/10)*(np/10)/20/pp2m;else tln=0;*/
  }
  printf("sl: "); /* check SL PART */
  if (*polyi) {
    Hpos = FTELL(F);
    FSEEK(F, -sl_NB, SEEK_END);
    SLpos = FTELL(F);
    if (SLpos - Hpos == L.NB)
      printf("NB o.k. ");
  } else if (sl_NB) {
    Ifn.replace(Ifn.size() - 4, 4, ".sl");
    fclose(F);
    if (nullptr == (F = fopen(Ifn.c_str(), "rb"))) {
      printf("Open %s failed", Ifn.c_str());
      exit(1);
    }
  } else
    puts("no .sl file");
  for (si = 0; si < sl_nNF; si++) {
    unsigned char uc[NUC_Nmax];
    int ch = fgetc(F);
    if (ch == EOF) {
      fprintf(stderr, "Error: DB_Check unexpected EOF reading vertex count\n");
      exit(1);
    }
    v = ch;
    if (v > VERT_Nmax) {
      fprintf(stderr, "Error: DB_Check vertex count %d out of range\n", v);
      exit(1);
    }
    ch = fgetc(F);
    if (ch == EOF) {
      fprintf(stderr, "Error: DB_Check unexpected EOF reading nuc count\n");
      exit(1);
    }
    nu = ch; /* assert(nu<=L.NUCmax); */
    AuxGet_uc(F, &nu, uc);
    if (ferror(F)) {
      fprintf(stderr, "Error: DB_Check sublattice read error\n");
      exit(1);
    }
    if ((*uc % 4) == 0)
      tSM++;
    else if ((*uc % 4) < 3)
      tNM++;
    Test_ucNF(&d, &v, &nu, uc, _P);
    if (!Check_sl_order(&v, &nu, uc)) {
      fputs("Error: DB_Check sublattice order violation\n", stderr);
      exit(1);
    }
  }
  if (sl_NB && (*dbi)) {
    SLpos = FTELL(F);
    if (sl_NB != SLpos) {
      fprintf(stderr,
              "Error: DB_Check sublattice byte position mismatch "
              "sl_NB=%lld SLpos=%lld\n",
              (long long)sl_NB, (long long)SLpos);
      exit(1);
    }
    printf("NB o.k. ");
  }
  printf("sm=%d=%d nm=%d=%lld", sl_SM, tSM, sl_NM, tNM);
  fflush(stdout);
  if ((sl_SM != (int)tSM) || (sl_NM != (int)tNM)) {
    puts("ERROR!!");
    exit(1);
  }
  tSM = tNM = 0; /* if(tln>1)printf("  p^2/2m=%ldkCY",tln); */
  Print_Expect(&L, out);

  printf("\nv:");
  if (*polyi)
    FSEEK(F, Hpos, SEEK_SET); /* check h-order */
  if (cF < 0)
    puts("");
  for (v = d + 1; v <= L.nVmax; v++) /* if(cF) */
    if (L.nNUC[v]) {
      Along nbsum = 0;
      if (cF > 0) {
        printf(" %d", v);
        fflush(stdout);
      }
      if (*dbi) {
        std::string Ifn_v = Ifn;
        std::string vxt = ".v";
        vxt += static_cast<char>(v / 10 + '0');
        vxt += static_cast<char>(v % 10 + '0');
        Ifn_v.replace(Ifn_v.size() - 4, 4, vxt);
        fclose(F);
        if (nullptr == (F = fopen(Ifn_v.c_str(), "rb"))) {
          printf("Ifn %s failed", Ifn_v.c_str());
          exit(1);
        }
      }
      for (nu = 1; nu <= L.NUCmax; nu++)
        if (L.NFnum[v][nu]) {
          unsigned char uc[NUC_Nmax];
          for (i = 0; i < L.NFnum[v][nu]; i++) {
            AuxGet_uc(F, &nu, uc);
            Check_hnf_order(&v, &nu, uc);
            switch (*uc % 4) {
            case 0:
              tSM++;
              break;
            case 1:;
            case 2:
              tNM++;
            case 3:;
            }

            if (cF < 0)
              if ((*uc % 4) % 3)
                Print_Missing_Mirror(&d, &v, &nu, uc, _P, out);
          }
          nbsum += nu * L.NFnum[v][nu];
        }
      if (*dbi) {
        Hpos = FTELL(F);
        if (Hpos != nbsum) {
          fprintf(stderr,
                  "Error: DB_Check honest byte position mismatch v=%d\n", v);
          exit(1);
        }
      }
      if (ferror(F)) {
        fprintf(stderr, "Error: DB_Check honest file read error v=%d\n", v);
        exit(1);
      }
    }
  printf("  sm=%d nm=%lld", tSM, tNM);
  printf("  order o.k.");
  if ((L.nSM != (int)tSM) || (L.nNM != tNM))
    puts("    CheckSum ERROR!");
  if (cF < 0)
    return;

  printf("\nv:");
  if (*polyi)
    FSEEK(F, Hpos, SEEK_SET); /* check h-NF */
  for (v = d + 1; v <= L.nVmax; v++)
    if (L.nNUC[v]) {
      int nbsum = 0;
      printf(" %d", v);
      fflush(stdout);
      if (*dbi) {
        std::string Ifn_v = Ifn;
        std::string vxt = ".v";
        vxt += static_cast<char>(v / 10 + '0');
        vxt += static_cast<char>(v % 10 + '0');
        Ifn_v.replace(Ifn_v.size() - 4, 4, vxt);
        fclose(F);
        if (nullptr == (F = fopen(Ifn_v.c_str(), "rb"))) {
          printf("Ifn %s failed", Ifn_v.c_str());
          exit(1);
        }
      }
      for (nu = 1; nu <= L.NUCmax; nu++)
        if (L.NFnum[v][nu]) {
          unsigned char uc[NUC_Nmax];
          for (i = 0; i < L.NFnum[v][nu]; i++) {
            AuxGet_uc(F, &nu, uc);
            Test_ucNF(&d, &v, &nu, uc, _P);
            switch (*uc % 4) {
            case 0:
              tSM++;
              break;
            case 1:;
            case 2:
              tNM++;
            case 3:;
            }
          }
          nbsum += nu * L.NFnum[v][nu];
        }
      if (*dbi) {
        Hpos = FTELL(F);
        if (Hpos != nbsum) {
          fprintf(stderr,
                  "Error: DB_Check honest NF byte position mismatch v=%d\n", v);
          exit(1);
        }
      }
      if (ferror(F)) {
        fprintf(stderr, "Error: DB_Check honest NF read error v=%d\n", v);
        exit(1);
      }
    }
  printf("  NF o.k.\n");
  if (ferror(F)) {
    fprintf(stderr, "Error: DB_Check file error after honest NF pass\n");
    exit(1);
  }
  fclose(F);
}

/*	A=h1&h2 B=h1&s2 C=s1&h2 D=s1&s2  1n=1-2  2n=2-1  1i=1-1n  2i=2-2n
 *      1n=(h1-A, s1-C-D)  2n=(h2-A, s2-B-D)  1i=(A,C+D)  2i=(A,B+D)

ln -s ../zzu.47 zzu.1 ; ln -s ../zzu.58 zzu.2
make && class.x -pi zzu.1 -ps zzu.2 -po zzu.1n
make && class.x -pi zzu.2 -ps zzu.1 -po zzu.2n
        class.x -pi zzu.1 -ps zzu.1n -po zzu.1i
        class.x -pi zzu.2 -ps zzu.2n -po zzu.2i

        class.x -pi zzu.1i -ps zzu.2n -po zzu.AD1
        class.x -pi zzu.2i -ps zzu.1n -po zzu.AD2
        class.x -pi zzu.1i -ps zzu.2i -po zzu.0C
        class.x -pi zzu.2i -ps zzu.1i -po zzu.0B
 *						      1 + 2 = 1n ^ 2n ^ AD  */

void Reduce_Aux_File(char *polyi, char *polys, char *dbsub, char *polyo,
                     FILE *out) {
  FILE *FI = fopen(polyi, "rb"), *FS, *FO = fopen(polyo, "wb");
  FInfoList FIi, FIs, FIo;
  Along Ipos, Spos = 00, HIpos, HSpos = 00, Opos, HOpos, IslNB, SslNB, tnb = 0;
  unsigned char ucI[NUC_Nmax], ucS[NUC_Nmax], *ucSL, *uc;
  int SLp[SL_Nmax];
  int IslNF, IslSM, IslNM, db = 0, vI, nuI, j, i;
  unsigned Ili, u;
  UPint Inp;
  int SslNF, SslSM, SslNM, dv = 0, vS, nuS, s, d;
  unsigned Sli;
  Along Snp;
  int slNF = 0, slSM = 0, slNM = 0, slNB = 0, slNP = 0, SmI = 00, nu, ms, v;
  UPint Oli = 0;
  std::string Sfn;
  Along HIPli[VERT_Nmax][NUC_Nmax], HSPli[VERT_Nmax][NUC_Nmax];

  if ((*polys == 0) != (*dbsub == 0))
    db = (*dbsub != 0);
  else {
    printf("Need ONE of ps=%s and ds=%s\n", polys, dbsub);
    exit(1);
  }
  if (!*polyi || !*polyo) {
    puts("With -ps or -ds you have to specify I/O files via -pi and -po");
    exit(1);
  }

  if (nullptr == FI) {
    printf("Cannot open %s", polyi);
    exit(1);
  }
  if (nullptr == FO) {
    printf("Cannot open %s", polyo);
    exit(1);
  }
  std::vector<unsigned char> ucSL_buffer(SL_Nmax * CperR_MAX * sizeof(char));
  ucSL = ucSL_buffer.data();
  Init_FInfoList(&FIs);

  if (db) {
    unsigned tln = 0;
    Along tNF = 0;
    Sfn = dbsub;
    Sfn += ".info";
    FS = fopen(Sfn.c_str(), "r");
    if (FS == nullptr) {
      puts("Info File not found");
      exit(1);
    }
    {
      Along sNF, sNM; /* start reading the file */
      if (fscanf(FS, "%d%d%d%d%d%lld%d%lld %lld %d%d%d%lld", &d, &i, &j, &nu,
                 &Sli, &sNF, &FIs.nSM, &sNM, &FIs.NB, &SslNF, &SslSM, &SslNM,
                 &SslNB) != 13) {
        fputs("Error: Subtract_Aux_from_DB malformed source info header\n",
              stderr);
        fclose(FS);
        exit(1);
      }
      FIs.nNF = sNF;
      FIs.nNM = sNM;
    }
    FIs.nV = i;
    FIs.nVmax = j;
    FIs.NUCmax = nu;
    for (i = 0; i < FIs.nV; i++) {
      if (fscanf(FS, "%d", &v) != 1) {
        fputs("Error: Subtract_Aux_from_DB expected vertex count\n", stderr);
        fclose(FS);
        exit(1);
      }
      if (fscanf(FS, "%d", &j) != 1) {
        fprintf(stderr,
                "Error: Subtract_Aux_from_DB expected nuc count for v=%d\n", v);
        fclose(FS);
        exit(1);
      }
      FIs.nNUC[v] = j;
      for (j = 0; j < FIs.nNUC[v]; j++) {
        if (fscanf(FS, "%d", &nu) != 1) {
          fprintf(stderr,
                  "Error: Subtract_Aux_from_DB expected nu value for v=%d\n",
                  v);
          fclose(FS);
          exit(1);
        }
        if (fscanf(FS, "%d", &FIs.NFnum[v][nu]) != 1) {
          fprintf(stderr,
                  "Error: Subtract_Aux_from_DB expected NF count for v=%d "
                  "nu=%d\n",
                  v, nu);
          fclose(FS);
          exit(1);
        }
        tln++;
        tNF += FIs.NFnum[v][nu];
      }
    }
    if (tln != Sli) {
      fprintf(stderr,
              "Error: Subtract_Aux_from_DB source list count mismatch\n");
      exit(1);
    }
    s = d;
    if (tNF != FIs.nNF) {
      fprintf(stderr,
              "Error: Subtract_Aux_from_DB source NF total mismatch "
              "tNF=%lld FIs.nNF=%lld\n",
              (long long)tNF, (long long)FIs.nNF);
      exit(1);
    }
    if (ferror(FS)) {
      fprintf(stderr, "Error: Subtract_Aux_from_DB source file read error\n");
      exit(1);
    }
    tNF = 0;
  } else {
    if (nullptr == (FS = fopen(polys, "rb"))) {
      printf("Cannot open %s", polys);
      exit(1);
    }
    if (fgetc(FS)) {
      puts("don't subtract aux files!");
      exit(1);
    }
    Read_Bin_Info(FS, &s, &Sli, &SslNF, &SslSM, &SslNM, &SslNB, &FIs);
    HSpos = FTELL(FS);
    FSEEK(FS, 0, SEEK_END);
    Spos = FTELL(FS);
    if (HSpos + FIs.NB + SslNB != Spos) {
      fprintf(stderr,
              "Error: Subtract_Aux_from_DB source file size mismatch\n");
      exit(1);
    }
    FSEEK(FS, -SslNB, SEEK_CUR);
  }
  d = fgetc(FI);
  if (d > 127)
    d = 128 * (d - 128) + fgetc(FI);
  Init_FInfoList(&FIi);
  if (d)
    printf("RecursionDepth=%d on %s\n", d, polyi);
  if (d < 128)
    fputc(d, FO);
  else {
    fputc(128 + (d) / 128, FO);
    fputc(d % 128, FO);
  }
  for (i = 0; i < d; i++) {
    j = fgetc(FI);
    fputc(j, FO);
  }
  Read_Bin_Info(FI, &d, &Ili, &IslNF, &IslSM, &IslNM, &IslNB, &FIi);
  if (d != s) {
    fprintf(stderr, "Error: Subtract_Aux_from_DB dimension mismatch %d != %d\n",
            d, s);
    exit(1);
  }
  HIpos = FTELL(FI);
  FSEEK(FI, 0, SEEK_END);
  Ipos = FTELL(FI);
  Inp = 2 * FIi.nNF - FIi.nSM - FIi.nNM;
  Snp = 2 * FIs.nNF - FIs.nSM - FIs.nNM;

  printf("Data on %s:  %d+%dsl  %lldb  (%dd)\n", polyi, Inp,
         2 * IslNF - IslNM - IslSM, Ipos, s);
  printf("Data on %s:  %lld+%dsl  %lldb  (%dd)\n", db ? dbsub : polys, Snp,
         2 * SslNF - SslNM - SslSM, db ? FIs.NB + SslNB : Spos, d);
  if (HIpos + FIi.NB + IslNB != Ipos) {
    fprintf(stderr, "Error: Subtract_Aux_from_DB input file size mismatch\n");
    exit(1);
  }
  FSEEK(FI, -IslNB, SEEK_CUR);

  if (db && SslNF) {
    Sfn.replace(Sfn.size() - 4, 4, ".sl");
    fclose(FS);
    if (nullptr == (FS = fopen(Sfn.c_str(), "rb"))) {
      printf("Open %s failed", Sfn.c_str());
      exit(1);
    }
  }
  s = 0;
  if (0 < SslNF)
    AuxGet_vn_uc(FS, &vS, &nuS, ucS);
  for (v = 0; v < IslNF; v++) {
    AuxGet_vn_uc(FI, &vI, &nuI, ucI);
    uc = &ucSL[2 + (SLp[slNF++] = slNB)];
    while (s < SslNF) {
      if (!(SmI = vS - vI))
        if (!(SmI = nuS - nuI))
          SmI = RIGHTminusLEFT(ucI, ucS, &nuI);
      if (SmI < 0) /* next S */
      {
        if ((++s) < SslNF)
          AuxGet_vn_uc(FS, &vS, &nuS, ucS);
      } else
        break;
    }
    if ((s < SslNF) && (SmI == 0)) /* I==S */
    {
      ms = (*ucS % 4);
      if ((ms % 3) && ((*ucI % 4) != ms)) /* put I==S */
      {
        uc[-2] = vI;
        uc[-1] = nuI;
        for (nu = 0; nu < nuI; nu++)
          uc[nu] = ucI[nu];
        *uc = (3 - ms) + 4 * (*uc / 4);
        slNM++;
        slNB += 2 + nuI;
        tnb += 2 + nuI;
      } else
        slNF--;
      if ((++s) < SslNF)
        AuxGet_vn_uc(FS, &vS, &nuS, ucS);
    } else /* put I */
    {
      uc[-2] = vI;
      uc[-1] = nuI;
      for (i = 0; i < nuI; i++)
        uc[i] = ucI[i];
      if ((ms = (*uc % 4))) {
        if (ms < 3)
          slNM++;
      } else
        slSM++;
      slNB += 2 + nuI;
    }
  } /* assert(tnb+slNB==IslNB+AslNB); */ /* SL done */

  printf("SL: %dnf %dsm %dnm %db -> ", slNF, slSM, slNM, slNB);
  tnb = 0;
  fflush(stdout);

  if (!db)
    FSEEK(FS, HSpos, SEEK_SET);
  FSEEK(FI, HIpos, SEEK_SET);
  Opos = HIpos;

  for (v = d + 1; v <= FIi.nVmax; v++)
    if (FIi.nNUC[v]) /* init HIPli */
      for (nu = 1; nu <= FIi.NUCmax; nu++)
        if (FIi.NFnum[v][nu]) {
          HIPli[v][nu] = Opos;
          Opos += nu * FIi.NFnum[v][nu];
        }
  if (Opos != HIpos + FIi.NB) {
    fprintf(stderr,
            "Error: Subtract_Aux_from_DB honest input offset mismatch\n");
    exit(1);
  }
  Opos = HSpos;
  for (v = d + 1; v <= FIs.nVmax; v++)
    if (FIs.nNUC[v]) /* init HSPli */
      for (nu = 1; nu <= FIs.NUCmax; nu++)
        if (FIs.NFnum[v][nu]) {
          if (db)
            if (v > dv) {
              Opos = 0;
              dv = v;
            }
          HSPli[v][nu] = Opos;
          Opos += nu * FIs.NFnum[v][nu];
        }
  if (!db && (Opos != HSpos + FIs.NB)) {
    fprintf(stderr, "Error: Subtract_Aux_from_DB source offset mismatch\n");
    exit(1);
  }
  dv = 0;
  Init_FInfoList(&FIo);

  for (v = d + 1; v <= FIi.nVmax; v++)
    if (FIi.nNUC[v]) /* FIo.NFnum[v][nu]?=0 */
      for (nu = 1; nu <= FIi.NUCmax; nu++)
        if (FIi.NFnum[v][nu]) {
          if (FIs.NFnum[v][nu]) /* check for add polys */
          {
            unsigned int *_nf = &FIo.NFnum[v][nu];
            /* (*_nf)=#NFnum(v,n) ?= 0 */
            unsigned int n, I_NF = FIi.NFnum[v][nu], S_NF = FIs.NFnum[v][nu];
            if (db)
              if (v > dv) {
                std::string vxt = ".v";
                vxt += static_cast<char>(v / 10 + '0');
                vxt += static_cast<char>(v % 10 + '0');
                std::string Sfn_v = Sfn;
                Sfn_v.replace(Sfn_v.size() - 4, 4, vxt);
                if (ferror(FS)) {
                  fprintf(stderr,
                          "Error: Subtract_Aux_from_DB source file read "
                          "error v=%d\n",
                          v);
                  exit(1);
                }
                fclose(FS);
                if (nullptr == (FS = fopen(Sfn_v.c_str(), "rb"))) {
                  printf("%s open failed", Sfn_v.c_str());
                  exit(1);
                }
                dv = v;
              }
            FSEEK(FI, HIPli[v][nu], SEEK_SET);
            FSEEK(FS, HSPli[v][nu], SEEK_SET);
            u = 0;
            if (0 < S_NF)
              AuxGet_uc(FS, &nu, ucS);
            for (n = 0; n < I_NF; n++) {
              AuxGet_uc(FI, &nu, ucI);
              while (u < S_NF) {
                SmI = RIGHTminusLEFT(ucI, ucS, &nu);
                if (SmI < 0) {
                  if ((++u) < S_NF)
                    AuxGet_uc(FS, &nu, ucS); /* get S */
                } else
                  break;
              }
              if ((u < S_NF) && (SmI == 0)) /* found I==S */
              {
                ms = (*ucS) % 4;
                if ((ms % 3) && (((*ucI) % 4) != ms))
                  (*_nf)++;
                if ((++u) < S_NF)
                  AuxGet_uc(FS, &nu, ucS);
              } else
                (*_nf)++;
              if (*_nf)
                break;
            }
          } else
            FIo.NFnum[v][nu] = FIi.NFnum[v][nu];
          if (FIo.NFnum[v][nu]) {
            FIo.nVmax = palp::max(FIo.nVmax, v);
            FIo.NUCmax = palp::max(FIo.NUCmax, nu);
            FIo.nNUC[v]++;
            Oli++;
          }
        }
  dv = 0;
  for (v = d + 1; v <= FIo.nVmax; v++)
    if (FIo.nNUC[v])
      FIo.nV++;

  fputc(d, FO);
  fputc(FIo.nV, FO);
  fputc(FIo.nVmax, FO);
  fputc(FIo.NUCmax, FO);
  fputUI(Oli, FO);
  Opos = FTELL(FO);
  j = InfoSize(0, Oli, &FIo) - 5 - sizeof(int);
  for (i = 0; i < j; i++)
    fputc(0, FO);

  HOpos = FTELL(FO); /* printf("Opos=%d   HOpos=%d\n\n\n",Opos,HOpos); */

  for (v = d + 1; v <= FIo.nVmax; v++)
    if (FIo.nNUC[v])
      for (nu = 1; nu <= FIo.NUCmax; nu++)
        if (FIo.NFnum[v][nu]) {
          unsigned I_NF = FIi.NFnum[v][nu], S_NF = FIs.NFnum[v][nu], O_NF = 0,
                   n;
          UPint neq = 0, peq = 0, pi = 0, po = 0;
          if (v > tnb) {
            printf(" %d", v);
            tnb = v;
            fflush(stdout);
            if (FIs.nNUC[v])
              if (db) {
                std::string vxt = ".v";
                vxt += static_cast<char>(v / 10 + '0');
                vxt += static_cast<char>(v % 10 + '0');
                std::string Sfn_v = Sfn;
                Sfn_v.replace(Sfn_v.size() - 4, 4, vxt);
                if (ferror(FS)) {
                  fprintf(stderr,
                          "Error: Subtract_Aux_from_DB source file read "
                          "error v=%d\n",
                          v);
                  exit(1);
                }
                fclose(FS);
                if (nullptr == (FS = fopen(Sfn_v.c_str(), "rb"))) {
                  printf("%s open failed", Sfn_v.c_str());
                  exit(1);
                }
              }
          }
          FSEEK(FI, HIPli[v][nu], SEEK_SET);
          FSEEK(FS, HSPli[v][nu], SEEK_SET);
          u = 0;
          if (u < S_NF)
            AuxGet_uc(FS, &nu, ucS);
          for (n = 0; n < I_NF; n++) {
            AuxGet_uc(FI, &nu, ucI);
            pi += 1 + (((*ucI) % 4) == 3);
            while (u < S_NF) {
              SmI = RIGHTminusLEFT(ucI, ucS, &nu);
              if (SmI < 0) /* get S */
              {
                u++;
                if (u < S_NF)
                  AuxGet_uc(FS, &nu, ucS);
              } else
                break;
            }
            if ((u < S_NF) && (SmI == 0)) /* put I-S */
            {
              int k, mm = 10 * (*ucI % 4) + (*ucS % 4);
              switch (mm) {
              case 33:
                peq++;
              case 22:;
              case 11:;
              case 00:
                peq++;
                break;
              case 31:;
              case 32:
                *ucI -= (*ucS % 4);
                peq++;
              case 12:;
              case 21:
                neq--;
                po++;
                O_NF++;
                FIo.nNM++;
                for (k = 0; k < nu; k++)
                  fputc(ucI[k], FO);
                break;
              case 13:;
              case 23:
                peq++;
                break;
              default:
                puts("inconsistens mirror flags");
                exit(1);
              }
              neq++;
              u++;
              if (u < S_NF)
                AuxGet_uc(FS, &nu, ucS);
            } else {
              int k;
              for (k = 0; k < nu; k++)
                fputc(ucI[k], FO);
              O_NF++;
              k = (*ucI) % 4;
              po += 1 + (k == 3);
              if (k) {
                if (k < 3)
                  FIo.nNM++;
              } else
                FIo.nSM++;
            }
          }
          if (pi - peq != po) {
            fprintf(stderr,
                    "Error: Subtract_Aux_from_DB checksum mismatch v=%d "
                    "nu=%d\n",
                    v, nu);
            exit(1);
          }
          if (O_NF + neq != I_NF) {
            fprintf(stderr,
                    "Error: Subtract_Aux_from_DB NF mismatch v=%d nu=%d\n", v,
                    nu);
            exit(1);
          } /* checksum(v,nu) */
          FIo.NFnum[v][nu] = O_NF;
          FIo.nNF += O_NF;
          FIo.NB += O_NF * nu;
        }
  tnb = 0;

  uc = &(ucSL[SLp[i = 0]]);
  dv = 0;
  bool sl_done = false;
  for (v = d + 1; v <= FIs.nVmax; v++)
    if (FIs.nNUC[v]) /* subtract S from SL */
      for (nu = 1; nu <= FIs.NUCmax; nu++)
        if (FIs.NFnum[v][nu]) {
          while ((uc[0] < v) || ((uc[0] == v) && (uc[1] < nu))) {
            uc = &(ucSL[SLp[++i]]);
            if (i >= slNF) {
              sl_done = true;
              break;
            }
          } /* go up to (v,nu) in SL-list */
          if (sl_done)
            break;
          if (db)
            if (v > dv) {
              std::string vxt = ".v";
              vxt += static_cast<char>(v / 10 + '0');
              vxt += static_cast<char>(v % 10 + '0');
              std::string Sfn_v = Sfn;
              Sfn_v.replace(Sfn_v.size() - 4, 4, vxt);
              if (ferror(FS)) {
                fprintf(stderr,
                        "Error: Subtract_Aux_from_DB source file read "
                        "error v=%d\n",
                        v);
                exit(1);
              }
              fclose(FS);
              if (nullptr == (FS = fopen(Sfn_v.c_str(), "rb"))) {
                printf("%s open failed", Sfn_v.c_str());
                exit(1);
              }
              dv = v;
            }

          if ((uc[0] == v) && (uc[1] == nu)) /* read from file and sort out */
          {
            FSEEK(FS, HSPli[v][nu], SEEK_SET);
            for (u = 0; u < FIs.NFnum[v][nu]; u++) {
              int HmSL;
              AuxGet_uc(FS, &nu, ucS); /* Test_ucNF(&d,&v,&nu,ucS);*/
              while (0 <
                     (HmSL = RIGHTminusLEFT(&uc[2], ucS, &nu))) /* next SL */
              {
                uc = &(ucSL[SLp[++i]]);
                if (i >= slNF) {
                  sl_done = true;
                  break;
                }
                if ((uc[0] != v) || (uc[1] != nu)) {
                  if ((uc[0] < v) || ((uc[0] == v) && (uc[1] < nu))) {
                    fprintf(stderr,
                            "Error: Subtract_Aux_from_DB SL order violation "
                            "v=%d nu=%d\n",
                            v, nu);
                    exit(1);
                  }
                  break;
                }
              }
              if (sl_done)
                break;
              if (HmSL == 0) /* remove SL */
              {
                int k, sms = uc[2] % 4, hms = (*ucS) % 4;
                switch (10 * hms + sms) {
                case 31:
                case 32:
                case 00:
                case 11:
                case 22:
                case 33:
                  for (k = i + 1; k < slNF; k++)
                    SLp[k - 1] = SLp[k];
                  slNB -= nu + 2;
                  if (sms == 0)
                    --slSM;
                  else if (sms < 3)
                    --slNM;
                  --slNF;
                  if (i >= slNF) {
                    sl_done = true;
                    break;
                  }
                  uc = &(ucSL[SLp[i]]);
                  if ((uc[0] != v) || (uc[1] != nu)) {
                    if ((uc[0] < v) || ((uc[0] == v) && (uc[1] < nu))) {
                      fprintf(stderr,
                              "Error: Subtract_Aux_from_DB SL order violation "
                              "v=%d nu=%d\n",
                              v, nu);
                      exit(1);
                    }
                    break;
                  }
                  break;
                case 12:
                case 21:
                  break;
                case 13:
                case 23:
                  ++slNM;
                  uc[2] -= hms;
                  break;
                default:
                  puts("inconsistent MS flags in SL-H");
                  exit(1);
                }
              } /* else (SL>H): hence next H */
            }
            if (ferror(FS)) {
              fprintf(stderr,
                      "Error: Subtract_Aux_from_DB source file read error\n");
              exit(1);
            }
          } else if ((uc[0] < v) || ((uc[0] == v) && (uc[1] < nu))) {
            fprintf(stderr,
                    "Error: Subtract_Aux_from_DB SL list order violation "
                    "v=%d nu=%d\n",
                    v, nu);
            exit(1);
          }
          if (sl_done)
            break;
        }

  for (i = 0; i < slNF; i++) /* write SL */
  {
    uc = &ucSL[SLp[i] + 2];
    v = uc[-2];
    nu = uc[-1];
    tnb += nu + 2;
    /* printf("#%d:  SLp=%d  v=%d  nu=%d\n",i,SLp[i],v,nu); */
    if (uc[-2] >= VERT_Nmax) {
      fprintf(stderr,
              "Error: Subtract_Aux_from_DB vertex count %d out of range\n",
              uc[-2]);
      exit(1);
    }
    fputc(uc[-2], FO);
    slNP += 1 + (((*uc) % 4) == 3);
    fputc(nu, FO);
    for (s = 0; s < nu; s++)
      fputc(uc[s], FO);
  }
  if (tnb != slNB) {
    fprintf(stderr,
            "Error: Subtract_Aux_from_DB sublattice byte count mismatch\n");
    exit(1);
  }
  if (slNP != 2 * slNF - slNM - slSM) {
    fprintf(stderr,
            "Error: Subtract_Aux_from_DB sublattice NP count mismatch\n");
    exit(1);
  }
  printf("\nd=%d v%d v<=%d n<=%d vn%d  %lld %d %lld %lld  %d %d %d %d\n", d,
         FIo.nV, FIo.nVmax, FIo.NUCmax, Oli, FIo.nNF, FIo.nSM, FIo.nNM, FIo.NB,
         slNF, slSM, slNM, slNB);

  FSEEK(FO, Opos, SEEK_SET);
  fputUI(FIo.nNF, FO); /* write info */
  fputUI(FIo.nSM, FO);
  fputUI(FIo.nNM, FO);
  fputUI(FIo.NB, FO);
  fputUI(slNF, FO);
  fputUI(slSM, FO);
  fputUI(slNM, FO);
  fputUI(slNB, FO);
  for (v = d + 1; v <= FIo.nVmax; v++)
    if (FIo.nNUC[v]) {
      i = 0;
      for (nu = 1; nu <= FIo.NUCmax; nu++)
        if (FIo.NFnum[v][nu])
          i++;
      fputc(v, FO);
      fputc(i, FO); /* v #nuc(v) */
      for (nu = 1; nu <= FIo.NUCmax; nu++)
        if (FIo.NFnum[v][nu]) {
          fputc(nu, FO);
          fputUI(FIo.NFnum[v][nu], FO); /* nuc #NF */
        }
    }

  printf("Writing %s: %lld+%dsl %lldm+%ds %lldb", polyo,
         2 * FIo.nNF - FIo.nNM - FIo.nSM, slNP,
         tnb = FIo.nNF - FIo.nNM - FIo.nSM, FIo.nSM, FIo.NB + slNB);
  /* if(tnb>99)
  {	long long tnp=(2*FIo.nNF-FIo.nNM-FIo.nSM)/1000; tnp*=tnp;
     tnp/=(2*tnb); printf("   [p^2/2m=%ldM]",tnp);
  }	*/
  Print_Expect(&FIo, out);
  puts("");
  if (ferror(FI)) {
    fprintf(stderr, "Error: Reduce_Aux_File input file read error\n");
    exit(1);
  }
  fclose(FI);
  if (ferror(FS)) {
    fprintf(stderr, "Error: Subtract_Aux_from_DB source file read error\n");
    exit(1);
  }
  fclose(FS);
  if (HOpos != FTELL(FO)) {
    fprintf(stderr, "Error: Subtract_Aux_from_DB output position mismatch\n");
    exit(1);
  }
  if (ferror(FO)) {
    fprintf(stderr, "Error: Subtract_Aux_from_DB output file write error\n");
    exit(1);
  }
  fclose(FO);
}
void Bin2a(char *polyi, int max, PolyPointList *_P, FILE *out) {
  FILE *F = fopen(polyi, "rb");
  FInfoList L;
  UPint list_num, tNF = 0;
  Along tNB = 0;
  int d, v, s, sl_nNF, sl_SM, sl_NM, sl_NB, mc = 0, MS, nu;
  unsigned i, j;
  unsigned char uc[POLY_Dmax * VERT_Nmax];
  VertexNumList V;
  EqList E;
  Long NF[POLY_Dmax][VERT_Nmax];
  Init_FInfoList(&L);

  if (F == nullptr) {
    printf("Input file %s not found\n", polyi);
    exit(1);
  }
  d = fgetc(F);
  if (d != 0) { /* for(i=0;i<d;i++) fgetc(F); */
    fprintf(stderr, "Error: Bin2a recursion depth %d must be 0\n", d);
    exit(1);
  }
  d = fgetc(F);
  L.nV = fgetc(F);
  L.nVmax = fgetc(F);
  L.NUCmax = fgetc(F);
  list_num = fgetUI(F);
  L.nNF = fgetUI(F);
  L.nSM = fgetUI(F);
  L.nNM = fgetUI(F);
  L.NB = fgetUI(F);
  sl_nNF = fgetUI(F);
  sl_SM = fgetUI(F);
  sl_NM = fgetUI(F);
  sl_NB = fgetUI(F);

  for (i = 0; i < L.nV; i++) {
    v = fgetc(F);
    L.nNUC[v] = fgetc(F); /* read #nuc's per #Vert */
    for (j = 0; j < L.nNUC[v]; j++) {
      L.NFnum[v][nu = fgetc(F)] = fgetUI(F); /* read nuc and #NF(v,nu)*/
      tNF += L.NFnum[v][nu];
      tNB += L.NFnum[v][nu] * (Along)nu;
    }
  }
  if ((unsigned int)(tNB - L.NB) != 0) {
    fputs("Error: Bin2a byte count mismatch\n", stderr);
    exit(1);
  }
  L.NB = tNB;
  if (tNF != L.nNF) {
    fprintf(stderr, "Error: Bin2a NF total mismatch tNF=%lld L.nNF=%lld\n",
            (long long)tNF, (long long)L.nNF);
    exit(1);
  }

  for (v = d + 1; v <= L.nVmax; v++)
    if (L.nNUC[v]) /* write  honest polys */
    {
      int I, J;
      for (nu = 1; nu <= L.NUCmax; nu++)
        for (j = 0; j < L.NFnum[v][nu]; j++) {
          for (s = 0; s < nu; s++)
            uc[s] = fgetc(F);
          UCnf2vNF(&d, &v, &nu, uc, NF, &MS);
          MS %= 4;
          _P->n = d;
          _P->np = v;
          for (I = 0; I < v; I++)
            for (J = 0; J < d; J++)
              _P->x[I][J] = NF[J][I];
          if (!Ref_Check(_P, &V, &E)) {
            fprintf(stderr,
                    "Error: Bin2a stored polytope not reflexive v=%d nu=%d\n",
                    v, nu);
            exit(1);
          }

          if (MS != 2) /* if(MS!=2) print NF */
            if (!max || Poly_Max_check(_P, &V, &E)) {
              mc++;
              Print_NF(out, &d, &v, NF);
            }
          if (MS > 1) /* if(MS>1); print Mirror */
            if (!max || Poly_Min_check(_P, &V, &E)) {
              mc++;
              Small_Make_Dual(_P, &V, &E);
              Make_Poly_NF(_P, &V, &E, NF);
              Print_NF(out, &d, &(V.nv), NF);
            }
        }
    }
  printf("np=%lld+%dsl  ", 2 * L.nNF - L.nSM - L.nNM,
         2 * sl_nNF - sl_SM - sl_NM);
  printf(/* write Finfo */
         "%dd  %dv<=%d n<=%d  %dnv  %lld %d %lld %lld  %d %d %d %d", d, L.nV,
         L.nVmax, L.NUCmax, list_num, L.nNF, L.nSM, L.nNM, L.NB, sl_nNF, sl_SM,
         sl_NM, sl_NB);
  if (max)
    printf("  r-max=%d", mc);
  puts("");
}
void DB_fromVF_toVT(DataBase *DB, int vf, int vt) {
  int v, n;
  DB->nNF = 0;
  for (v = vf; v <= vt; v++)
    for (n = 0; n < NUC_Nmax; n++)
      DB->nNF += DB->NFnum[v][n];
  if (!DB->nNF) {
    fprintf(stderr, "No NF with %d<=v<=%d\n", vf, vt);
    exit(1);
  }
  while (0 == DB->nNUC[vt]) {
    vt--;
    if (vf > vt) {
      fprintf(stderr, "Error: DB_fromVF_toVT range collapsed %d > %d\n", vf,
              vt);
      exit(1);
    }
  }
  DB->v = vf;
  DB->nVmax = vt;
}
void Bin2aDBsl(char *dbi, int max, int vf, int vt, PolyPointList *_P,
               FILE *out) {
  FILE *F;
  FInfoList L;
  int d, v, nu, i, j, list_num, mc = 0, MS, sl_nNF, sl_SM, sl_NM, sl_NB,
                                tSM = 0, tNM = 0;
  std::string Ifn = dbi;
  Long NF[POLY_Dmax][VERT_Nmax];
  VertexNumList V;
  EqList E;
  Ifn += ".info";
  F = fopen(Ifn.c_str(), "r");
  if (F == nullptr) {
    puts("Info File not found");
    exit(1);
  }
  Init_FInfoList(&L); /* start reading the file */
  if (fscanf(F, "%d%d%d%d%d%lld%d%lld %lld %d%d%d%d", &d, &i, &j, &nu,
             &list_num, &L.nNF, &L.nSM, &L.nNM, &L.NB, &sl_nNF, &sl_SM, &sl_NM,
             &sl_NB) != 13) {
    fputs("Error: Bin2aDBsl malformed info header\n", stderr);
    fclose(F);
    exit(1);
  }
  L.nV = i;
  L.nVmax = j;
  L.NUCmax = nu;
  if (sl_NB && (vf == 2) && (vt == VERT_Nmax - 1)) {
    Ifn.replace(Ifn.size() - 4, 4, ".sl");
    fclose(F);
    if (nullptr == (F = fopen(Ifn.c_str(), "rb"))) {
      printf("Open %s failed", Ifn.c_str());
      exit(1);
    }
  } else /* puts("no .sl file"); */
  {
    std::unique_ptr<DataBase> DB;
    Open_DB(dbi, &DB, 0);
    DataBase *DBRaw = DB.get();
    if ((DBRaw->v < vf) || (DBRaw->nVmax > vt))
      DB_fromVF_toVT(DBRaw, vf, vt); /*  read only vf <= v <= vt :: */
    for (i = 0; Read_H_poly_from_DB(DBRaw, _P); i++)
      Print_PPL(_P, "");
    printf("#poly=%d\n", i);
    Close_DB(DBRaw);
  }
  for (i = 0; i < sl_nNF; i++) {
    int I, J;
    unsigned char uc[NUC_Nmax];
    int ch = fgetc(F);
    if (ch == EOF) {
      fprintf(stderr, "Error: Bin2aDBsl unexpected EOF reading vertex count\n");
      exit(1);
    }
    v = ch;
    if (v > VERT_Nmax) {
      fprintf(stderr, "Error: Bin2aDBsl vertex count %d out of range\n", v);
      exit(1);
    }
    ch = fgetc(F);
    if (ch == EOF) {
      fprintf(stderr, "Error: Bin2aDBsl unexpected EOF reading nuc count\n");
      exit(1);
    }
    nu = ch;
    AuxGet_uc(F, &nu, uc);
    if (ferror(F)) {
      fprintf(stderr, "Error: Bin2aDBsl sublattice read error\n");
      exit(1);
    }
    if ((*uc % 4) == 0)
      tSM++;
    else if ((*uc % 4) < 3)
      tNM++;
    Test_ucNF(&d, &v, &nu, uc, _P);

    UCnf2vNF(&d, &v, &nu, uc, NF, &MS);
    MS %= 4;
    _P->n = d;
    _P->np = v;
    for (I = 0; I < v; I++)
      for (J = 0; J < d; J++)
        _P->x[I][J] = NF[J][I];
    if (!Ref_Check(_P, &V, &E)) {
      fprintf(stderr,
              "Error: Bin2aDBsl sublattice polytope not reflexive v=%d nu=%d\n",
              v, nu);
      exit(1);
    }
    if (MS != 2) /* if(MS!=2) print NF */
      if (!max || Poly_Max_check(_P, &V, &E)) {
        mc++;
        Print_NF(out, &d, &v, NF);
      }
    if (MS > 1) /* if(MS>1); print Mirror */
      if (!max || Poly_Min_check(_P, &V, &E)) {
        mc++;
        Small_Make_Dual(_P, &V, &E);
        Make_Poly_NF(_P, &V, &E, NF);
        Print_NF(out, &d, &(V.nv), NF);
      }
  }
  printf("np=%lld+%dsl  ", 2 * L.nNF - L.nSM - L.nNM,
         2 * sl_nNF - sl_SM - sl_NM);
  printf(                                            /* write Finfo */
         "%dd  %dv<=%d n<=%d  %dnv ... %d %d %d %d", /* %d %d %d %lld */
         d, L.nV, L.nVmax, L.NUCmax, list_num,
         /* L.nNF,L.nSM,L.nNM,L.NB, */ sl_nNF, sl_SM, sl_NM, sl_NB);
  if (max)
    printf("  r-max=%d", mc);
  puts("");
}
void Bin_2_ascii(char *polyi, char *dbin, int max, int vf, int vt,
                 PolyPointList *P, FILE *out) {
  if (*polyi)
    Bin2a(polyi, max, P, out);
  else if (*dbin)
    Bin2aDBsl(dbin, max, vf, vt, P, out);
  else
    puts("With -b[2a] you have to specify input via -pi or -di");
}

/*  ======================================================================  */
/*  ==========                                                  ==========  */
/*  ==========    Hodge-database routines                       ==========  */
/*  ==========                                                  ==========  */
/*  ======================================================================  */

namespace {
constexpr int Hod_Dif_max = 480;
constexpr int Hod_Min_max = 251;
} // namespace

#if (POLY_Dmax < 6)
void DB_to_Hodge(char *dbin, char *dbout, int vfrom, int vto, PolyPointList *_P,
                 FILE *out) {

  /* Read the database, write the Hodge numbers */

  DataBase DB;
  VertexNumList V;
  Long VPM[EQUA_Nmax][VERT_Nmax];
  EqList E;
  time_t Tstart;
  char *fx;
  std::string dbname = dbin;
  std::string dbhname = dbout;
  std::string dbext;
  unsigned char uc_poly[NUC_Nmax];
  int d, v, nu, i, j, list_num, sl_nNF, sl_SM, sl_NM, sl_NB, dh,
      nnf_vd[VERT_Nmax][Hod_Dif_max + 1], nnf_v[VERT_Nmax];
  BaHo BH;
  FILE *Faux[Hod_Dif_max + 1];
  FILE *Fvinfo;
  std::unique_ptr<PolyPointList> _PD = std::make_unique<PolyPointList>();

  if (!*dbin || !*dbout) {
    puts("You have to specify I/O database names via -di and -do");
    exit(1);
  }

  for (i = 0; i <= Hod_Dif_max; i++)
    for (j = 0; j < VERT_Nmax; j++)
      nnf_vd[j][i] = 0;
  dbname += ".info";
  dbhname += ".vinfo";
  Fvinfo = fopen(dbhname.c_str(), "a");

  printf("Reading %s: ", dbname.c_str());
  fflush(0);

  /* read the info-file: */
  DB.Finfo = fopen(dbname.c_str(), "r");
  if (DB.Finfo == nullptr) {
    fprintf(stderr, "Error: CY_hodge_split cannot open %s\n", dbname.c_str());
    exit(1);
  }
  if (fscanf(DB.Finfo, "%d  %d %d %d  %d  %lld %d %lld %lld  %d %d %d %d", &d,
             &DB.nV, &DB.nVmax, &DB.NUCmax, &list_num, &DB.nNF, &DB.nSM,
             &DB.nNM, &DB.NB, &sl_nNF, &sl_SM, &sl_NM, &sl_NB) != 13) {
    fputs("Error: CY_hodge_split malformed info file\n", stderr);
    exit(1);
  }
  printf("%lldp (%dsl) %lldnf %lldb\n",
         2 * (DB.nNF) - DB.nSM - DB.nNM + 2 * sl_nNF - sl_SM - sl_NM,
         2 * sl_nNF - sl_SM - sl_NM, DB.nNF + sl_nNF, DB.NB + sl_NB);

  for (v = 1; v < VERT_Nmax; v++) {
    DB.nNUC[v] = 0;
    for (nu = 0; nu < NUC_Nmax; nu++)
      DB.NFnum[v][nu] = 0;
  }

  for (i = 0; i < DB.nV; i++) {
    if (fscanf(DB.Finfo, "%d", &v) != 1) {
      fputs("Error: CY_hodge_split expected vertex count\n", stderr);
      exit(1);
    }
    if (fscanf(DB.Finfo, "%d", &(DB.nNUC[v])) != 1) {
      fprintf(stderr, "Error: CY_hodge_split expected nuc count for v=%d\n", v);
      exit(1);
    }
    for (j = 0; j < DB.nNUC[v]; j++) {
      if (fscanf(DB.Finfo, "%d", &nu) != 1) {
        fprintf(stderr, "Error: CY_hodge_split expected nu value for v=%d\n",
                v);
        exit(1);
      }
      if (fscanf(DB.Finfo, "%d", &(DB.NFnum[v][nu])) != 1) {
        fprintf(stderr,
                "Error: CY_hodge_split expected NF count for v=%d nu=%d\n", v,
                nu);
        exit(1);
      }
    }
  }

  if (ferror(DB.Finfo)) {
    printf("File error in %s\n", dbname.data());
    exit(1);
  }
  fclose(DB.Finfo);
  fflush(stdout);

  printf("Reading DB-files, calculating Hodge numbers, writing aux-files:\n");

  /* read the DB-files and calculate Hodge numbers */
  for (v = vfrom; v <= vto; v++)
    if (DB.nNUC[v]) {
      int nd = 0;
      std::string dbfile = dbin;
      dbfile += ".v";
      dbfile += static_cast<char>('0' + v / 10);
      dbfile += static_cast<char>('0' + v % 10);
      dbext = dbhname;
      dbext.replace(dbext.size() - 5, 1, std::to_string(0));
      dbext.replace(dbext.size() - 4, 1, std::to_string(0));
      dbext.replace(dbext.size() - 3, 1, std::to_string(0));
      Tstart = time(nullptr);
      printf("v=%d: ", v);
      fflush(0);
      DB.Fv[v] = fopen(dbfile.c_str(), "rb");
      if (DB.Fv[v] == nullptr) {
        fprintf(stderr, "Error: CY_hodge_split cannot open %s\n",
                dbfile.c_str());
        exit(1);
      }
      for (nu = 0; nu <= DB.NUCmax; nu++)
        for (i = 0; i < DB.NFnum[v][nu]; i++) {
          int mirror = 0;
          int MS;
          unsigned char c;
          for (j = 0; j < nu; j++)
            uc_poly[j] = fgetc(DB.Fv[v]);
          uc_nf_to_P(_P, &MS, &d, &v, &nu, uc_poly);
          if (!Ref_Check(_P, &V, &E)) {
            fprintf(stderr,
                    "Error: CY_hodge_split stored polytope not reflexive "
                    "v=%d nu=%d\n",
                    v, nu);
            exit(1);
          }
          if (V.nv != v) {
            fprintf(stderr,
                    "Error: CY_hodge_split vertex count mismatch V.nv=%d "
                    "v=%d\n",
                    V.nv, v);
            exit(1);
          }
          Make_VEPM(_P, &V, &E, VPM);
          Complete_Poly(VPM, &E, V.nv, _P);
          Make_Dual_Poly(_P, &V, &E, _PD.get());
          RC_Calc_BaHo(_P, &V, &E, _PD.get(), &BH);
          if (BH.h1[1] < BH.h1[2])
            mirror = 1;
          if (mirror)
            dh = BH.h1[2] - BH.h1[1];
          else
            dh = BH.h1[1] - BH.h1[2];
          if (!nnf_vd[v][dh] || (dh > 250)) {
            dbext = dbhname;
            dbext.replace(dbext.size() - 5, 1, std::to_string(dh / 100));
            dbext.replace(dbext.size() - 4, 1, std::to_string((dh / 10) % 10));
            dbext.replace(dbext.size() - 3, 1, std::to_string(dh % 10));
            Faux[dh] = fopen(dbext.c_str(), "ab");
          }
          if (!nnf_vd[v][dh])
            nd++;
          c = BH.h1[2 - mirror];
          fputc(c, Faux[dh]);
          c = BH.mp / 256 + 4 * BH.mv;
          fputc(c, Faux[dh]);
          c = BH.mp % 256;
          fputc(c, Faux[dh]);
          c = BH.np / 256 + 4 * BH.nv;
          fputc(c, Faux[dh]);
          c = BH.np % 256;
          fputc(c, Faux[dh]);
          c = nu + 64 * mirror;
          fputc(c, Faux[dh]);
          for (j = 0; j < nu; j++)
            fputc(uc_poly[j], Faux[dh]);
          if (dh > 250)
            fclose(Faux[dh]);
          (nnf_vd[v][dh])++;
          (nnf_v[v])++;
        }

      if (ferror(DB.Fv[v])) {
        printf("File error in %s\n", dbname.c_str());
        exit(1);
      }
      fclose(DB.Fv[v]);
      for (dh = 0; dh <= 250; dh++)
        if (nnf_vd[v][dh]) {
          if (ferror(Faux[dh])) {
            printf("File error at dh=%d\n", dh);
            exit(1);
          }
          fclose(Faux[dh]);
        }

      fprintf(Fvinfo, "%d %d %d\n", v, nd, nnf_v[v]);
      for (dh = 0; dh <= Hod_Dif_max; dh++)
        if (nnf_vd[v][dh])
          fprintf(Fvinfo, "%d %d ", dh, nnf_vd[v][dh]);
      fprintf(Fvinfo, "\n");
      printf(" %d NF (%ds)\n", nnf_v[v], (int)difftime(time(nullptr), Tstart));
      fflush(0);
    }
  fclose(Fvinfo);
}

void Sort_Hodge(char *dbaux, char *dbout) {
  /* Sort from v-chi-format to chi-h12-format */

  time_t Tstart;
  char *fax;
  std::vector<char> dbaname(6 + strlen(dbaux) + File_Ext_NCmax);
  char *fhx;
  std::vector<char> dbhname(6 + strlen(*dbout ? dbout : dbaux) +
                            File_Ext_NCmax);
  int v, i, j, dh, nd, nnf_d[Hod_Dif_max + 1],
      nnf_vd[VERT_Nmax][Hod_Dif_max + 1], nnf_h[Hod_Min_max + 1],
      nnf_v[VERT_Nmax];
  FILE *Fh[Hod_Min_max + 1];
  FILE *Fchia, *Fvinfo, *Fhinfo;

  for (i = 0; i <= Hod_Dif_max; i++) {
    nnf_d[i] = 0;
    for (j = 0; j < VERT_Nmax; j++)
      nnf_vd[j][i] = 0;
  }

  if (!*dbout)
    dbout = dbaux;
  strcpy(dbaname.data(), dbaux);
  strcat(dbaname.data(), ".vinfo");
  strcpy(dbhname.data(), dbout);
  strcat(dbhname.data(), ".hinfo");
  fhx = &dbhname[strlen(dbout) + 1];
  fax = &dbaname[strlen(dbaux) + 1];

  printf("Reading %s\n", dbaname.data());
  fflush(0);

  /* read the info-file: */
  Fvinfo = fopen(dbaname.data(), "r");
  if (Fvinfo == nullptr) {
    fprintf(stderr, "Error: Sort_Hodge_files cannot open %s\n", dbaname.data());
    exit(1);
  }
  while (fscanf(Fvinfo, "%d", &v) == 1) {
    if (fscanf(Fvinfo, "%d  %d", &nd, &(nnf_v[v])) != 2) {
      fprintf(stderr, "Error: Sort_Hodge_files expected nd nnf_v for v=%d\n",
              v);
      fclose(Fvinfo);
      exit(1);
    }
    /* printf("%d %d %d  ", v, nd, nnf_v[v]  ); fflush(0); */
    for (i = 0; i < nd; i++) {
      if (fscanf(Fvinfo, "%d", &dh) != 1) {
        fprintf(stderr,
                "Error: Sort_Hodge_files expected dh entry %d for v=%d\n", i,
                v);
        fclose(Fvinfo);
        exit(1);
      }
      if (fscanf(Fvinfo, "%d", &(nnf_vd[v][dh])) != 1) {
        fprintf(stderr,
                "Error: Sort_Hodge_files expected count for v=%d dh=%d\n", v,
                dh);
        fclose(Fvinfo);
        exit(1);
      }
      nnf_d[dh] += nnf_vd[v][dh];
    }
  }
  if (ferror(Fvinfo)) {
    printf("File error in %s\n", dbaname.data());
    exit(1);
  }
  fclose(Fvinfo);

  printf("Sorting and writing Hodge-files:\n");
  fflush(0);
  Tstart = time(nullptr);

  Fhinfo = fopen(dbhname.data(), "w");
  if (Fhinfo == nullptr) {
    fprintf(stderr, "Error: Sort_Hodge_files cannot create %s\n",
            dbhname.data());
    exit(1);
  }

  /* Sort the Hodge&Poly-Data */
  for (dh = 0; dh <= Hod_Dif_max; dh++)
    if (nnf_d[dh]) {
      char aext[8], hext[9];
      int h12, nh = 0;
      unsigned char c, nuc;
      aext[0] = 'v';
      aext[3] = 'd';
      aext[4] = '0' + dh / 100;
      aext[5] = '0' + (dh / 10) % 10;
      aext[6] = '0' + dh % 10;
      aext[7] = 0;
      hext[0] = 'd';
      hext[1] = '0' + dh / 100;
      hext[2] = '0' + (dh / 10) % 10;
      hext[3] = '0' + dh % 10;
      hext[4] = 'h';
      hext[5] = hext[6] = hext[7] = hext[8] = 0;

      printf("dh=%d: %dNF...", dh, nnf_d[dh]);
      fflush(0);
      for (j = 0; j <= Hod_Min_max; j++)
        nnf_h[j] = 0;

      for (v = 2; v < VERT_Nmax; v++)
        if (nnf_vd[v][dh]) {
          aext[1] = '0' + v / 10;
          aext[2] = '0' + v % 10;
          strcpy(fax, aext);
          Fchia = fopen(dbaname.data(), "rb");
          for (i = 0; i < nnf_vd[v][dh]; i++) {
            h12 = fgetc(Fchia);
            if (!nnf_h[h12]) {
              hext[5] = '0' + h12 / 100;
              hext[6] = '0' + (h12 / 10) % 10;
              hext[7] = '0' + h12 % 10;
              strcpy(fhx, hext);
              Fh[h12] = fopen(dbhname.data(), "wb");
              nh++;
            }
            nnf_h[h12]++;
            c = fgetc(Fchia);
            fputc(c, Fh[h12]);
            if (v != c / 4) {
              printf("v=%d, hp.mv=%d", v, (int)(c / 4));
              exit(1);
            }
            /* for (j=0;j<4;j++) {c=fgetc(Fchia); fputc(c,Fh[h12]);}
            for (j=0;j<c%64;j++) fputc(fgetc(Fchia),Fh[h12]);}*/
            for (j = 0; j < 3; j++)
              fputc(fgetc(Fchia), Fh[h12]);
            nuc = fgetc(Fchia);
            fputc(nuc, Fh[h12]);
            c = fgetc(Fchia);
            {
              int ms = c % 4;
              if (ms)
                c += 3 - ms;
            }
            fputc(c, Fh[h12]);
            for (j = 0; j < nuc % 64 - 1; j++)
              fputc(fgetc(Fchia), Fh[h12]);
          }
          if (ferror(Fchia)) {
            printf("File error in Fchia at dh=%d v=%d\n", dh, v);
            exit(1);
          }
          fclose(Fchia);
        }

      fprintf(Fhinfo, "%d %d %d\n", dh, nh, nnf_d[dh]);
      for (h12 = 0; h12 <= Hod_Min_max; h12++)
        if (nnf_h[h12]) {
          fprintf(Fhinfo, "%d %d ", h12, nnf_h[h12]);
          nnf_d[dh] -= nnf_h[h12];
          if (ferror(Fh[h12])) {
            printf("File error at dh=%d h12=%d\n", dh, h12);
            exit(1);
          }
          fclose(Fh[h12]);
        }
      fprintf(Fhinfo, "\n");
      if (nnf_d[dh]) {
        printf("nnf_d[%d]!=sum nnf_dh[%d][h12]!", dh, dh);
        exit(1);
      }
      printf(" sorted\n");
    }

  if (ferror(Fhinfo)) {
    printf("File error in Fhinfo\n");
    exit(1);
  }
  fclose(Fhinfo);
  printf("  done (%ds)\n", (int)difftime(time(nullptr), Tstart));
  fflush(stdout);
}

void Test_Hodge_db(char *dbname) {

  time_t Tstart;
  char *fhx;
  std::vector<char> filename(6 + strlen(dbname) + File_Ext_NCmax);
  int i, j, dh, h12, nh, nnf_sum, nnf_d[Hod_Dif_max + 1],
      nnf_dh[Hod_Dif_max + 1][Hod_Min_max + 1];
  FILE *Fh;
  FILE *Fhinfo;

  for (i = 0; i <= Hod_Dif_max; i++) {
    nnf_d[i] = 0;
    for (j = 0; j <= Hod_Min_max; j++)
      nnf_dh[i][j] = 0;
  }

  strcpy(filename.data(), dbname);
  strcat(filename.data(), ".hinfo");
  fhx = &filename[strlen(dbname) + 1];

  printf("Reading %s\n", filename.data());
  fflush(0);

  /* read the info-file: */
  Fhinfo = fopen(filename.data(), "r");
  if (Fhinfo == nullptr) {
    fprintf(stderr, "Error: Test_Hodge_db cannot open %s\n", filename.data());
    exit(1);
  }
  while (fscanf(Fhinfo, "%d", &dh) == 1) {
    if (fscanf(Fhinfo, "%d  %d", &nh, &(nnf_d[dh])) != 2) {
      fprintf(stderr, "Error: Test_Hodge_db expected nh nnf_d for dh=%d\n", dh);
      fclose(Fhinfo);
      exit(1);
    }
    nnf_sum = 0;
    for (i = 0; i < nh; i++) {
      if (fscanf(Fhinfo, "%d", &h12) != 1) {
        fprintf(stderr,
                "Error: Test_Hodge_db expected h12 entry %d for dh=%d\n", i,
                dh);
        fclose(Fhinfo);
        exit(1);
      }
      if (fscanf(Fhinfo, "%d", &(nnf_dh[dh][h12])) != 1) {
        fprintf(stderr,
                "Error: Test_Hodge_db expected count for dh=%d h12=%d\n", dh,
                h12);
        fclose(Fhinfo);
        exit(1);
      }
      nnf_sum += nnf_dh[dh][h12];
    }
    if (nnf_sum != nnf_d[dh]) {
      printf("nnf_d[%d]!=sum nnf_dh[%d][h12]!", dh, dh);
      exit(1);
    }
  }
  if (ferror(Fhinfo)) {
    printf("File error in %s\n", filename.data());
    exit(1);
  }
  fclose(Fhinfo);

  printf("Analysing Hodge-files:\n");
  fflush(0);
  Tstart = time(nullptr);

  /* Analyse Hodge&Poly-Data */
  for (dh = 0; dh <= Hod_Dif_max; dh++)
    if (nnf_d[dh]) {
      char hext[9];
      hext[0] = 'd';
      hext[1] = '0' + dh / 100;
      hext[2] = '0' + (dh / 10) % 10;
      hext[3] = '0' + dh % 10;
      hext[4] = 'h';
      hext[5] = hext[6] = hext[7] = hext[8] = 0;

      printf("dh=%d: %dNF...", dh, nnf_d[dh]);
      fflush(0);
      nnf_sum = 0;

      for (h12 = 0; h12 <= Hod_Min_max; h12++)
        if (nnf_dh[dh][h12]) {
          /* unsigned char uc_poly[NUC_Nmax]; */
          int c1, /* c2, */ nuc;
          hext[5] = '0' + h12 / 100;
          hext[6] = '0' + (h12 / 10) % 10;
          hext[7] = '0' + h12 % 10;
          strcpy(fhx, hext);
          Fh = fopen(filename.data(), "rb");
          if (Fh == nullptr) {
            fprintf(stderr, "Error: Test_Hodge_db cannot open %s\n",
                    filename.data());
            exit(1);
          }
          while ((c1 = fgetc(Fh)) != EOF) {
            nnf_sum++;
            /* c2= */ fgetc(Fh);
            /* mv=c1/4; mp=(c1%4)*256+c2; */
            c1 = fgetc(Fh); /* c2= */
            fgetc(Fh);
            /* nv=c1/4; np=(c1%4)*256+c2; */
            c1 = fgetc(Fh);
            /* mirror=c1/64; */
            nuc = c1 % 64;
            for (j = 0; j < nuc; j++) /* uc_poly[j]= */
              fgetc(Fh);
          /* uc_nf_to_P(_P, &MS, &d, &mv, &nuc, uc_poly); */ }
          if (ferror(Fh)) {
            printf("File error in Fh at dh=%d h12=%d\n", dh, h12);
            exit(1);
          }
          fclose(Fh);
        }

      printf("%d\n", nnf_sum);
    }

  printf("  done (%ds)\n", (int)difftime(time(nullptr), Tstart));
  fflush(stdout);
}

void Extract_from_Hodge_db(char *dbname, char *x_string, PolyPointList *_P,
                           FILE *out) {

  time_t Tstart;
  char c = *x_string, hext[9], com[64];
  std::vector<char> filename(6 + strlen(dbname) + File_Ext_NCmax);
  char *fhx;
  unsigned char uc_poly[NUC_Nmax];
  int E = 998, H1 = 0, H2 = 0, M = 0, V = 0, N = 0, F = 0, L = 1000, i = 0, j,
      dh, h12, nh, nnf_sum, nnf_d[Hod_Dif_max + 1],
      nnf_dh[Hod_Dif_max + 1][Hod_Min_max + 1], HD_from = 0,
      HD_to = Hod_Dif_max, HM_from = 0, HM_to = Hod_Min_max, c1, c2, nuc,
      mirror, nv, np, mv, mp, d = 4, MS, is_poly, is_dual, true_H1, true_H2,
      max_mv = 34;
  FILE *Fh;
  FILE *Fhinfo;

  hext[0] = 'd';
  hext[4] = 'h';
  hext[8] = 0;

  while (c) {
    if (c == 'E') {
      int neg = 0;
      c = x_string[++i];
      if (c == '-') {
        neg = 1;
        c = x_string[++i];
      }
      if ((c - '0' >= 0) && (c - '0' <= 9))
        E = 0;
      while ((c - '0' >= 0) && (c - '0' <= 9)) {
        E = 10 * E + c - '0';
        c = x_string[++i];
      }
      if (E % 2) {
        puts("The Euler number E number must be even!");
        return;
      }
      if (neg)
        E *= -1;
    } else if (c == 'H') {
      c = x_string[++i];
      while ((c - '0' >= 0) && (c - '0' <= 9)) {
        H1 = 10 * H1 + c - '0';
        c = x_string[++i];
      }
    } else if (c == ':') {
      c = x_string[++i];
      while ((c - '0' >= 0) && (c - '0' <= 9)) {
        H2 = 10 * H2 + c - '0';
        c = x_string[++i];
      }
    } else if (c == 'M') {
      c = x_string[++i];
      while ((c - '0' >= 0) && (c - '0' <= 9)) {
        M = 10 * M + c - '0';
        c = x_string[++i];
      }
    } else if (c == 'V') {
      c = x_string[++i];
      while ((c - '0' >= 0) && (c - '0' <= 9)) {
        V = 10 * V + c - '0';
        c = x_string[++i];
      }
    } else if (c == 'N') {
      c = x_string[++i];
      while ((c - '0' >= 0) && (c - '0' <= 9)) {
        N = 10 * N + c - '0';
        c = x_string[++i];
      }
    } else if (c == 'F') {
      c = x_string[++i];
      while ((c - '0' >= 0) && (c - '0' <= 9)) {
        F = 10 * F + c - '0';
        c = x_string[++i];
      }
    } else if (c == 'L') {
      c = x_string[++i];
      L = 0;
      while ((c - '0' >= 0) && (c - '0' <= 9)) {
        L = 10 * L + c - '0';
        c = x_string[++i];
      }
    } else {
      printf("`%c' is not valid input", c);
      return;
    }
  }

  if (H1 && H2) {
    E = 2 * (H1 - H2);
  } else if (H1 && (E != 998))
    H2 = H1 - E / 2;
  else if (H2 && (E != 998))
    H1 = H2 + E / 2;

  if (H1 > 491 || H2 > 491 || H1 + H2 > 502 ||
      (abs(E) != 998 && abs(E) > 960)) {
    puts("#NF: 0   Note the range for Hodge numbers:");
    puts("         h11,h12<=491, h11+h12<=502, |E|<=960.");
    return;
  }

  if ((V && (V < 5 || V > 33)) || (F && (F < 5 || F > 33))) {
    puts("#NF: 0   Note the range [5,33] for vertex/facet numbers!");
    return;
  }

  if ((M && (M < 6 || M > 680)) || (N && (N < 6 || N > 680))) {
    puts("#NF: 0   Note the range [6,680] for point numbers!");
    return;
  }

  if (V && (V < max_mv))
    max_mv = V;
  if (F && (F < max_mv))
    max_mv = F;
  if (M && (M - 1 < max_mv))
    max_mv = M - 1;
  if (N && (N - 1 < max_mv))
    max_mv = N - 1;

  for (i = 0; i <= Hod_Dif_max; i++) {
    nnf_d[i] = 0;
    for (j = 0; j <= Hod_Min_max; j++)
      nnf_dh[i][j] = 0;
  }

  strcpy(filename.data(), dbname);
  strcat(filename.data(), ".hinfo");
  fhx = &filename[strlen(dbname) + 1];

  /* printf("Reading %s\n",filename.data()); fflush(0); */

  /* read the info-file: */
  Fhinfo = fopen(filename.data(), "r");
  if (Fhinfo == nullptr) {
    fprintf(stderr, "Error: Extract_from_Hodge_db cannot open %s\n",
            filename.data());
    exit(1);
  }
  while (fscanf(Fhinfo, "%d", &dh) == 1) {
    if (fscanf(Fhinfo, "%d  %d", &nh, &(nnf_d[dh])) != 2) {
      fprintf(stderr,
              "Error: Extract_from_Hodge_db expected nh nnf_d for dh=%d\n", dh);
      fclose(Fhinfo);
      exit(1);
    }
    nnf_sum = 0;
    for (i = 0; i < nh; i++) {
      if (fscanf(Fhinfo, "%d", &h12) != 1) {
        fprintf(stderr,
                "Error: Extract_from_Hodge_db expected h12 entry %d for "
                "dh=%d\n",
                i, dh);
        fclose(Fhinfo);
        exit(1);
      }
      if (fscanf(Fhinfo, "%d", &(nnf_dh[dh][h12])) != 1) {
        fprintf(stderr,
                "Error: Extract_from_Hodge_db expected count for dh=%d "
                "h12=%d\n",
                dh, h12);
        fclose(Fhinfo);
        exit(1);
      }
      nnf_sum += nnf_dh[dh][h12];
    }
    if (nnf_sum != nnf_d[dh]) {
      printf("nnf_d[%d]!=sum nnf_dh[%d][h12]!", dh, dh);
      exit(1);
    }
  }

  if (ferror(Fhinfo)) {
    printf("File error in %s\n", filename.data());
    exit(1);
  }
  fclose(Fhinfo);

  /* printf("Analysing Hodge-files:\n"); fflush(0); */
  Tstart = time(nullptr);

  /* Analyse Hodge&Poly-Data */

  nnf_sum = 0;
  if (H1 && H2) {
    HM_from = (H1 < H2 ? H1 : H2);
    HM_to = HM_from;
  } else if (H1 || H2) {
    HM_to = palp::max(H1, H2);
    HM_from = palp::max(0, HM_to - Hod_Dif_max);
  }

  for (h12 = HM_from; h12 <= HM_to; h12++) {

    if (abs(E) != 998) {
      HD_from = abs(E / 2);
      HD_to = abs(E / 2);
    } else if ((H1 || H2) && (h12 < HM_to)) {
      HD_from = HM_to - h12;
      HD_to = HM_to - h12;
    } else {
      HD_from = 0;
      HD_to = Hod_Dif_max;
    }

    for (dh = HD_from; dh <= HD_to; dh++)
      if (nnf_dh[dh][h12]) {

        hext[1] = '0' + dh / 100;
        hext[2] = '0' + (dh / 10) % 10;
        hext[3] = '0' + dh % 10;
        hext[5] = '0' + h12 / 100;
        hext[6] = '0' + (h12 / 10) % 10;
        hext[7] = '0' + h12 % 10;
        strcpy(fhx, hext);
        Fh = fopen(filename.data(), "rb");
        if (Fh == nullptr) {
          fprintf(stderr, "Error: Extract_from_Hodge_db cannot open %s\n",
                  filename.data());
          exit(1);
        }
        while ((c1 = fgetc(Fh)) != EOF) {
          mv = c1 / 4;
          if (mv > max_mv)
            break;
          c2 = fgetc(Fh);
          mp = (c1 % 4) * 256 + c2;
          c1 = fgetc(Fh);
          c2 = fgetc(Fh);
          nv = c1 / 4;
          np = (c1 % 4) * 256 + c2;
          c1 = fgetc(Fh);
          mirror = c1 / 64;
          nuc = c1 % 64;
          for (j = 0; j < nuc; j++)
            uc_poly[j] = fgetc(Fh);

          if (mirror) {
            true_H1 = h12;
            true_H2 = dh + h12;
          } else {
            true_H2 = h12;
            true_H1 = dh + h12;
          }
          is_poly = 0;
          if ((abs(E) == 998) || (E == 2 * (true_H1 - true_H2)))
            if (!H1 || (true_H1 == H1))
              if (!H2 || (true_H2 == H2))
                if (!M || (M == mp))
                  if (!V || (V == mv))
                    if (!N || (N == np))
                      if (!F || (F == nv))
                        is_poly = 1;
          is_dual = 0;
          if ((abs(E) == 998) || (E == 2 * (true_H2 - true_H1)))
            if (!H1 || (true_H2 == H1))
              if (!H2 || (true_H1 == H2))
                if (!M || (M == np))
                  if (!V || (V == nv))
                    if (!N || (N == mp))
                      if (!F || (F == mv))
                        is_dual = 1;
          if (!is_poly && !is_dual)
            continue;

          uc_nf_to_P(_P, &MS, &d, &mv, &nuc, uc_poly);

          if (is_poly) {
            /* if(!MS) puts("!MS");
               if(!mirror) puts("!mirror"); */
            if (++nnf_sum > L) {
              printf("Exceeded limit of %d\n", L);
              return;
            }
            snprintf(com, sizeof(com), "M:%d %d N:%d %d H:%d,%d [%d]", mp, mv,
                     np, nv, true_H1, true_H2, 2 * (true_H1 - true_H2));
            Print_PPL(_P, com);
          }
          if (is_dual && MS) {
            VertexNumList VNL;
            EqList EL;
            Long NF[POLY_Dmax][VERT_Nmax];
            if (++nnf_sum > L) {
              printf("Exceeded limit of %d\n", L);
              return;
            }
            snprintf(com, sizeof(com), "M:%d %d N:%d %d H:%d,%d [%d]", np, nv,
                     mp, mv, true_H2, true_H1, 2 * (true_H2 - true_H1));
            Find_Equations(_P, &VNL, &EL);
            Small_Make_Dual(_P, &VNL, &EL);
            Make_Poly_NF(_P, &VNL, &EL, NF, out);
            Print_Matrix(NF, _P->n, VNL.nv, com, out);
          }
        }
        if (ferror(Fh)) {
          printf("File error in Fh at dh=%d h12=%d\n", dh, h12);
          exit(1);
        }
        fclose(Fh);
      }
  }
  printf("#NF: %d\n", nnf_sum);
  printf("  done (%ds)\n", (int)difftime(time(nullptr), Tstart));
  fflush(stdout);
}

void Test_Hodge_file(char *filename, PolyPointList *_P) {
  FILE *Ft = fopen(filename, "rb");
  unsigned char uc_poly[NUC_Nmax];

  int nuc, mirror, nv, np, mv, mp, j, c1, c2, d = 4, MS;
  if (Ft == nullptr) {
    fprintf(stderr, "Error: Test_Hodge_file cannot open %s\n", filename);
    exit(1);
  }
  while ((c1 = fgetc(Ft)) != EOF) {
    c2 = fgetc(Ft);
    mv = c1 / 4;
    mp = (c1 % 4) * 256 + c2;
    c1 = fgetc(Ft);
    c2 = fgetc(Ft);
    nv = c1 / 4;
    np = (c1 % 4) * 256 + c2;
    c1 = fgetc(Ft);
    mirror = c1 / 64;
    nuc = c1 % 64;
    for (j = 0; j < nuc; j++)
      uc_poly[j] = fgetc(Ft);
    uc_nf_to_P(_P, &MS, &d, &mv, &nuc, uc_poly);
    Print_PPL(_P, "");
    printf("nuc=%d mirror=%d mv=%d mp=%d nv=%d np=%d MS=%d\n", nuc, mirror, mv,
           mp, nv, np, MS);
  }
  fclose(Ft);
}
#endif

/*  =====================================================================  */
void Open_DB(char *dbin, std::unique_ptr<DataBase> *_DB, int info) {
  int i, j, v, nu;
  char ext[4];
  std::string dbname;
  if (*dbin == 0) {
    *_DB = nullptr;
    return;
  }
  auto DBOwner = std::make_unique<DataBase>();
  DataBase *DB = DBOwner.get();
  DB->readHucNF_TotNF = 0;
  dbname = dbin;
  dbname.resize(dbname.size() + File_Ext_NCmax + 1, '\0');
  dbname[strlen(dbin)] = '.';
  char *fx = dbname.data() + strlen(dbin) + 1;
  fx[0] = '\0';
  dbname += "info";
  if (info) {
    printf("Reading %s: ", dbname.data());
    fflush(0);
  }
  DB->Finfo = fopen(dbname.data(), "r");
  if (DB->Finfo == nullptr) {
    fprintf(stderr, "Error: Open_DB cannot open %s\n", dbname.data());
    exit(1);
  }
  if (fscanf(DB->Finfo, "%d  %d %d %d  %d  %lld %d %lld %lld  %d %d %d %d",
             &DB->d, &DB->nV, &DB->nVmax, &DB->NUCmax, &DB->list_num, &DB->nNF,
             &DB->nSM, &DB->nNM, &DB->NB, &DB->sl_nNF, &DB->sl_SM, &DB->sl_NM,
             &DB->sl_NB) != 13) {
    fputs("Error: Open_DB malformed info file\n", stderr);
    exit(1);
  }
  if (info)
    printf("%lldp (%dsl) %lldnf %lldb\n",
           2 * (DB->nNF) - DB->nSM - DB->nNM + 2 * DB->sl_nNF - DB->sl_SM -
               DB->sl_NM,
           2 * DB->sl_nNF - DB->sl_SM - DB->sl_NM, DB->nNF + DB->sl_nNF,
           DB->NB + DB->sl_NB);
  for (v = 1; v < VERT_Nmax; v++) {
    DB->nNUC[v] = 0;
    for (nu = 0; nu < NUC_Nmax; nu++)
      DB->NFnum[v][nu] = 0;
  }
  for (i = 0; i < DB->nV; i++) {
    if (fscanf(DB->Finfo, "%d", &v) != 1) {
      fputs("Error: Open_DB expected vertex count\n", stderr);
      exit(1);
    }
    if (fscanf(DB->Finfo, "%d", &(DB->nNUC[v])) != 1) {
      fprintf(stderr, "Error: Open_DB expected nuc count for v=%d\n", v);
      exit(1);
    }
    for (j = 0; j < DB->nNUC[v]; j++) {
      if (fscanf(DB->Finfo, "%d", &nu) != 1) {
        fprintf(stderr, "Error: Open_DB expected nu value for v=%d\n", v);
        exit(1);
      }
      if (fscanf(DB->Finfo, "%d", &(DB->NFnum[v][nu])) != 1) {
        fprintf(stderr, "Error: Open_DB expected NF count for v=%d nu=%d\n", v,
                nu);
        exit(1);
      }
    }
  }
  if (ferror(DB->Finfo)) {
    printf("File error in %s\n", dbname.data());
    exit(1);
  }
  fclose(DB->Finfo);
  ext[0] = 'v';
  ext[3] = 0;
  DB->v = DB->p = DB->nu = 0;
  for (v = DB->d + 1; v <= DB->nVmax; v++)
    if (DB->nNUC[v]) {
      ext[1] = '0' + v / 10;
      ext[2] = '0' + v % 10;
      strcpy(fx, ext);
      DB->Fv[v] = fopen(dbname.data(), "rb");
      if (DB->Fv[v] == nullptr) {
        fprintf(stderr, "Error: Open_DB cannot open %s\n", dbname.data());
        exit(1);
      }
      if (0 == DB->v) {
        DB->v = v;
        while (0 == DB->NFnum[v][DB->nu])
          (DB->nu)++;
      }
    }
}
void Close_DB(DataBase *DB) {
  int v;
  if (DB == nullptr)
    return;
  for (v = 1 + DB->d; v <= DB->nVmax; v++)
    if (DB->nNUC[v]) {
      if (ferror(DB->Fv[v])) {
        printf("File error at v=%d\n", v);
        exit(1);
      }
      fclose(DB->Fv[v]);
    }
}
int Read_H_ucNF_from_DB(DataBase *DB, unsigned char *uc) /* p=next read pos */
{
  int rest;
  if (DB == nullptr) {
    fputs("Error: Read_H_ucNF_from_DB called with nullptr database\n", stderr);
    exit(1);
  }
  rest = DB->NFnum[DB->v][DB->nu] - DB->p;
  if (rest < 0) {
    fprintf(stderr, "Error: Read_H_ucNF_from_DB negative rest %d\n", rest);
    exit(1);
  }
  if (rest == 0) /* search next v/nu */
  {
    int v = DB->v, nu = DB->nu; /* search nu(v): */
    while (DB->nu < DB->NUCmax)
      if (DB->NFnum[DB->v][++(DB->nu)])
        break;
    if ((nu == DB->nu) || (DB->NFnum[DB->v][DB->nu] == 0)) /* search v: */
    {
      if (DB->v < DB->nVmax) {
        while (0 == DB->nNUC[++DB->v])
          if (DB->v >= DB->nVmax) {
            fprintf(stderr,
                    "Error: Read_H_ucNF_from_DB exhausted vertex range\n");
            exit(1);
          }
        DB->nu = 0;
        while (DB->nu < DB->NUCmax)
          if (DB->NFnum[DB->v][++(DB->nu)])
            break;
      } else {
        if (DB->readHucNF_TotNF != DB->nNF) {
          fprintf(stderr,
                  "Error: Read_H_ucNF_from_DB total NF mismatch "
                  "totNF=%lld DB->nNF=%lld\n",
                  (long long)DB->readHucNF_TotNF, (long long)DB->nNF);
          exit(1);
        }
        return 0;
      }
    }
    if ((v < DB->v) || (nu < DB->nu)) {
      rest = DB->NFnum[DB->v][DB->nu];
      DB->p = 0;
    }
  }
  if (DB->p > DB->NFnum[DB->v][DB->nu]) {
    fprintf(stderr,
            "Error: Read_H_ucNF_from_DB position exceeds NFnum v=%d nu=%d\n",
            DB->v, DB->nu);
    exit(1);
  }
  if (rest == 0) {
    fputs("Error: Read_H_ucNF_from_DB zero rest\n", stderr);
    exit(1);
  }
  if (DB->readHucNF_TotNF >= DB->nNF) {
    fprintf(
        stderr,
        "Error: Read_H_ucNF_from_DB totNF %lld not less than DB->nNF %lld\n",
        (long long)DB->readHucNF_TotNF, (long long)DB->nNF);
    exit(1);
  }
  AuxGet_uc(DB->Fv[DB->v], &DB->nu, uc);
  ++DB->p;
  DB->readHucNF_TotNF++;
  return 1;
}

int Read_SLucNF_from_DB(void) {
  puts("Read_SLucNF_from_DB: to be implemented");
  exit(1);
  return 0;
}

int Read_SLpoly_from_DB(void) {
  puts("Read_SLpoly_from_DB: to be implemented");
  exit(1);
  return 0;
}

int Read_H_poly_from_DB(DataBase *DB, PolyPointList *P) {
  int i, j, MS;
  Long NF[POLY_Dmax][VERT_Nmax];
  VertexNumList V;
  EqList E;
  if (DB->last_ms3) {
    --DB->last_uc[0];
  } else {
    if (0 == Read_H_ucNF_from_DB(DB, DB->last_uc))
      return 0;
  }
  UCnf2vNF(&DB->d, &DB->v, &DB->nu, DB->last_uc, NF, &MS);
  P->n = DB->d;
  P->np = DB->v;
  for (i = 0; i < P->np; i++)
    for (j = 0; j < P->n; j++)
      P->x[i][j] = NF[j][i];
  MS %= 4;
  DB->last_ms3 = (MS == 3);
  if (MS == 2) {
    if (!Ref_Check(P, &V, &E)) {
      fprintf(stderr,
              "Error: Read_H_poly_from_DB stored polytope not reflexive "
              "v=%d nu=%d\n",
              DB->v, DB->nu);
      exit(1);
    }
    P->np = E.ne;
    for (i = 0; i < P->np; i++)
      for (j = 0; j < P->n; j++)
        P->x[i][j] = E.e[i].a[j];
  }
  return 1;
}

int Read_H_poly_from_DB_or_inFILE(DataBase *DB, PolyPointList *P, FILE *in) {
  if ((DB == nullptr) || (in == nullptr)) {
    CWS cws;
    return Read_CWS_PP(&cws, P, in);
  } else
    return Read_H_poly_from_DB(DB, P);
}
/*  =====================================================================  */

/*  =====================================================================  */
/*      ==============		from Subpoly.c:		=============       */
/* void Make_All_Sublat(NF_List *_L, int n, int v, subl_int diag[POLY_Dmax],
                     subl_int u[][VERT_Nmax], char *mFlag)		*/
/* create all inequivalent decompositions diag=s*t  into upper
   triangular matrices s,t;
   t*(first lines of u) becomes the poly P;
   choose the elements of t (columns rising, line # falling),
   calculate corresponding element of s at the same time;
   finally calculate P and add to list  */
/* void MakePolyOnSublat(NF_List *_L, subl_int x[VERT_Nmax][VERT_Nmax],
                      int v, int f, int *max_order, char *mFlag)	*/
/*   Decompose the VPM x as x=w*diag*u, where w and u are SL(Z);
     the first lines of u are the coordinates on the coarsest lattices */
/*      ==============	     End of "from Subpoly.c"	 =============      */

/*   ===================	Sublattice: phv	      ===================   */
Long AuxGxP(Long *Gi, Long *V, int *d) {
  Long x = 0;
  int j;
  for (j = 0; j < *d; j++)
    x += Gi[j] * V[j];
  return x;
}

/*   Glz x (lincomb(Points)) -> Diag: if(index>1) print vertices of G*P & D */
namespace {
constexpr bool ph_test = true;
}
void Aux_Print_CoverPoly(int *I, int *d, int *N, Long *X[POLY_Dmax],
                         Long G[][POLY_Dmax], Long *D, int *x,
                         Long Z[][VERT_Nmax], Long *M, int r, FILE *out) {
  int i, j, dia = 1, err = 0;
  fprintf(out, "%d %d    index=%d  D=%ld", *d, *N, *I, D[0]);
  for (i = 1; i < *d; i++)
    printf(" %ld", D[i]);
  for (i = 0; i < r; i++) {
    fprintf(out, " /Z%ld:", M[i]);
    for (j = 0; j < *N; j++)
      fprintf(out, " %ld", Z[i][j]);
  }
  printf("  #%d\n", *x);
  for (i = 0; i < *d; i++) {
    for (j = 0; j < *N; j++) {
      Long Xij = AuxGxP(G[i], X[j], d);
      if (0 != (Xij % D[i]))
        err = 1;
      printf("%ld ", Xij);
    }
    puts("");
  }
  if constexpr (ph_test) {
    if ((dia == 0) || err) {
      int i, j;
      printf("D=");
      for (i = 0; i < *d; i++)
        printf("%ld ", D[i]);
      printf("   index[%d]=%d\n", *x, *I);
      for (i = 0; i < *d; i++) {
        printf("G=");
        for (j = 0; j < *d; j++)
          printf("%2ld ", G[i][j]);
        printf("   G.P=");
        for (j = 0; j < *N; j++)
          printf("%2ld ", AuxGxP(G[i], X[j], d));
        puts("");
      }
      exit(1);
    }
  }
}
void Aux_Print_SLpoly(int *I, int *d, int *N, Long *X[POLY_Dmax],
                      Long G[][POLY_Dmax], Long *D, int *x, FILE *out) {
  int i, j, dia = 1, err = 0;
  printf("%d %d    index=%d  D=%ld", *d, *N, *I, D[0]);
  for (i = 1; i < *d; i++)
    printf(" %ld", D[i]);
  printf("  #%d\n", *x);
  for (i = 0; i < *d; i++) { /* int z=1;*/
    for (j = 0; j < *N; j++) {
      Long Xij = AuxGxP(G[i], X[j], d);
      if (0 != (Xij % D[i]))
        err = 1;
      printf("%ld ", Xij / D[i]);
    }
    puts("");
  }
  if constexpr (ph_test) {
    if ((dia == 0) || err) {
      int i, j;
      printf("D=");
      for (i = 0; i < *d; i++)
        printf("%ld ", D[i]);
      printf("   index[%d]=%d\n", *x, *I);
      for (i = 0; i < *d; i++) {
        printf("G=");
        for (j = 0; j < *d; j++)
          printf("%2ld ", G[i][j]);
        printf("   G.P=");
        for (j = 0; j < *N; j++)
          printf("%2ld ", AuxGxP(G[i], X[j], d));
        puts("");
      }
      exit(1);
    }
  }
}
void Aux_Make_Dual(PolyPointList *P, VertexNumList *V, EqList *E) {
  Long VM[VERT_Nmax][POLY_Dmax];
  int i, j, d = P->n, e = E->ne, v = V->nv;
  if (e > VERT_Nmax) {
    fprintf(stderr,
            "Error: Aux_Make_Dual equation count %d exceeds VERT_Nmax\n", e);
    exit(1);
  }
  P->np = V->nv = e;
  E->ne = v;
  for (i = 0; i < v; i++)
    for (j = 0; j < d; j++)
      VM[i][j] = P->x[V->v[i]][j];
  for (i = 0; i < e; i++) {
    for (j = 0; j < d; j++)
      P->x[i][j] = E->e[i].a[j];
    V->v[i] = i;
  }
  for (i = 0; i < v; i++) {
    for (j = 0; j < d; j++)
      E->e[i].a[j] = VM[i][j];
    E->e[i].c = 1;
  }
  if (!Ref_Check(P, V, E)) {
    fputs("Error: Aux_Make_Dual result not reflexive\n", stderr);
    exit(1);
  }
}
void PrintVPHMusage(void);
int Make_Lattice_Basis(int d, int p, Long *P[POLY_Dmax], /* index=det(D) */
                       Long G[][POLY_Dmax],
                       Long *D); /* G x P generates diagonal lattice D */
void PH_Sublat_Polys(char *dbin, int omitFIP, PolyPointList *_P, char sF,
                     FILE *in, FILE *out) {
  EqList E;
  VertexNumList V;
  int x = 0, K, B, I = 1; /* index>I only */
  Long *RelPts[POINT_Nmax];
  std::unique_ptr<DataBase> DB;
  if (*dbin)
    Open_DB(dbin, &DB, 0);
  K = ((sF != 'P') && (sF != 'H') && (sF != 'Q') &&
       (sF != 'B')); /*  K=='CoverPoly'  */
  if (('1' < sF) && (sF <= '9'))
    I = sF - '0';
  B = (omitFIP == 2); /* 'q' for index >I */
  while (Read_H_poly_from_DB_or_inFILE(DB.get(), _P, in)) {
    Long D[POLY_Dmax], G[POLY_Dmax][POLY_Dmax], PM[VERT_Nmax][VERT_Nmax];
    int index, N = 0; /* if(!Ref_Check(_P,&V,&E)) Print_PPL(_P,""); */
    if (B && (_P->n != 4)) {
      fprintf(stderr,
              "Error: PH_Sublat_Polys Brower group requires d=4, got d=%d\n",
              _P->n);
      exit(1);
    }
    if (!Ref_Check(_P, &V, &E)) {
      fprintf(stderr,
              "Error: PH_Sublat_Polys input polytope not reflexive d=%d\n",
              _P->n);
      exit(1);
    }
    /* Aux_Make_Dual(_P,&V,&E); */ /* don't dualize: take M-lattice poly */
    Make_VEPM(_P, &V, &E, PM);
    _P->np = V.nv;
    Complete_Poly(PM, &E, V.nv, _P);
    ++x;
    if (omitFIP == 0)
      for (N = 0; N < _P->np; N++)
        RelPts[N] = _P->x[N];
    else if (omitFIP < 3) /* Omit facet-IPs */
    {
      int p;
      for (p = 0; p < _P->np; p++) {
        int e, z = 0;
        for (e = 0; e < E.ne; e++)
          if (0 == Eval_Eq_on_V(&E.e[e], _P->x[p], _P->n))
            z++;
        if (z > omitFIP)
          RelPts[N++] = _P->x[p];
      }
      if (V.nv > N) {
        fprintf(stderr,
                "Error: PH_Sublat_Polys relevant point count N=%d less than "
                "vertices %d\n",
                N, V.nv);
        exit(1);
      } /* count <E,.>=0; if(n>1) add_to_RelPts; */
    } else if (omitFIP == 3) /* Omit all non-vertices */
    {
      for (N = 0; N < V.nv; N++)
        RelPts[N] = _P->x[N];
      /*   { int tnv=V.nv; Find_Equations(_P,&V,&E);assert(V.nv==tnv);
             for(tnv=0;tnv<V.nv;tnv++)assert(V.v[tnv]<V.nv);} */
    } else {
      puts("something wrong in PH_Sublat_Polys");
      exit(1);
    }
    if (K) {
      Long Z[POLY_Dmax][VERT_Nmax], M[POLY_Dmax];
      int r;
      index = Sublattice_Basis(_P->n, N, RelPts, Z, M, &r, G, D);
      if (B)
        if (D[2] == 1)
          continue;
      if (omitFIP == 3)
        if (D[1] == 1)
          continue;
      if (index <= 0) {
        fprintf(stderr,
                "Error: PH_Sublat_Polys non-positive sublattice index %d\n",
                index);
        exit(1);
      }
      if (index <= I)
        continue;
      if (!Ref_Check(_P, &V, &E)) {
        fprintf(stderr,
                "Error: PH_Sublat_Polys cover polytope not reflexive\n");
        exit(1);
      }
      Aux_Print_CoverPoly(&index, &_P->n, &N, RelPts, G, D, &x, Z, M, r, out);
    } else {
      index = Make_Lattice_Basis(_P->n, N, RelPts, G, D);
      if (1 == index)
        continue;
      if (B)
        if (D[2] == 1)
          continue;
      if (omitFIP == 3)
        if (D[1] == 1)
          continue;
      if (index <= 0) {
        fprintf(stderr,
                "Error: PH_Sublat_Polys non-positive lattice index %d\n",
                index);
        exit(1);
      }
      if (!Ref_Check(_P, &V, &E)) {
        fprintf(stderr,
                "Error: PH_Sublat_Polys sublattice polytope not reflexive\n");
        exit(1);
      }
      Print_VL(_P, &V, "");
      Aux_Print_SLpoly(&index, &_P->n, &N, RelPts, G, D, &x, out);
    }
  }
  if (*dbin)
    Close_DB(DB.get());
}
/*	Lattice generated by vertices; UT-decomp of diag	*/
void V_Sublat_Polys(char mr, char *dbin, char *polyi, char *polyo,
                    PolyPointList *_P, FILE *in, FILE *out) {
  NF_List _L_obj;
  NF_List *_L = &_L_obj;
  int max_order = 1;
  EqList E;
  VertexNumList V;
  int x = 0;
  Long *RelPts[VERT_Nmax];
  std::unique_ptr<DataBase> DB;
  if (*dbin)
    Open_DB(dbin, &DB, 0);
  if (!(*polyo)) {
    puts("You have to specify an output file via -po in -sv-mode.");
    printf("For more help use option `-h'\n");
    exit(1);
  }
  _L->of = 0;
  _L->rf = 0;
  _L->iname = polyi;
  _L->oname = polyo;
  _L->dbname = dbin;
  Init_NF_List(_L);
  _L->SL = 0;
  while (Read_H_poly_from_DB_or_inFILE(DB.get(), _P, in)) {
    Long D[POLY_Dmax], G[POLY_Dmax][POLY_Dmax];
    int index, N;
    if (!Ref_Check(_P, &V, &E)) {
      fprintf(stderr,
              "Error: V_Sublat_Polys input polytope not reflexive d=%d\n",
              _P->n);
      exit(1);
    }
    for (N = 0; N < V.nv; N++)
      RelPts[N] = _P->x[V.v[N]];
    ++x;
    index = Make_Lattice_Basis(_P->n, N, RelPts, G, D);
    if (1 == index)
      continue;
    if (index <= 0) {
      fprintf(stderr, "Error: V_Sublat_Polys non-positive index %d\n", index);
      exit(1);
    }
    if (index > max_order)
      max_order = index;

    {
      int i, j;
      subl_int diag[POLY_Dmax], U[POLY_Dmax][VERT_Nmax];
      for (i = 0; i < _P->n; i++) {
        diag[i] = D[i];
        for (j = 0; j < V.nv; j++) {
          int k;
          U[i][j] = 0;
          for (k = 0; k < _P->n; k++)
            U[i][j] += G[i][k] * RelPts[j][k];
          if ((D[i] != 0) && (U[i][j] % D[i] != 0)) {
            fprintf(stderr,
                    "Error: V_Sublat_Polys lattice vector not divisible by "
                    "D[%d]=%ld\n",
                    i, (long)D[i]);
            exit(1);
          }
          U[i][j] /= D[i];
        }
      }
      Make_All_Sublat(_L, _P->n, V.nv, diag, U, &mr, _P, out);
    }
  }
  if (*dbin)
    Close_DB(DB.get());
  printf("max_order=%d\n", max_order);
  Write_List_2_File(polyo, _L, out);
  _L->TIME = time(nullptr);
  fputs(ctime(&_L->TIME), stdout);
}
void VPHM_Sublat_Polys(char sFlag, char mr, char *dbin, char *polyi,
                       char *polyo, PolyPointList *_P, FILE *in, FILE *out) {
  switch (sFlag) { /* if(dbin=0) read from in; */
  case 'p':
  case 'P':
    PH_Sublat_Polys(dbin, 0, _P, sFlag, in, out);
    break;
  case 'h':
  case 'H':
    PH_Sublat_Polys(dbin, 1, _P, sFlag, in, out);
    break;
  case 'b':
  case 'B':
    PH_Sublat_Polys(dbin, 2, _P, sFlag, in, out);
    break;
  case 'q':
  case 'Q':
    PH_Sublat_Polys(dbin, 3, _P, sFlag, in, out);
    break;
  case 'v':
  case 'V':
    V_Sublat_Polys(mr, dbin, polyi, polyo, _P, in, out);
    break;
  case 'm':
  case 'M':
    Find_Sublat_Polys(mr, dbin, polyi, polyo, _P, in, out);
    break;
  default:
    if (('1' < sFlag) && (sFlag <= '9'))
      PH_Sublat_Polys(dbin, 3, _P, sFlag, in, out);
    else {
      puts("-s# requires that # is in {v,p,h,b,m,q}");
      PrintVPHMusage();
    }
  }
}
void PrintVPHMusage(void) {
  puts("	-sh ... gen by codim>1 points (omit IPs of facets)");
  puts("	-sp ... gen by all points");
  puts("	-sb ... generated by dim<=1 (edges), print if rank=2	");
  puts("	-sq ... generated by vertices,       print if rank=3	");
  puts("	    q,b currently assume that dim=4");
  exit(1);
}

void Bin_2_ANF(char *polyi, int max, PolyPointList *_P, FILE *out) {
  FILE *F = fopen(polyi, "rb");
  FInfoList L;
  UPint list_num, tNF = 0;
  Along tNB = 0;
  int d, v, s, sl_nNF, sl_SM, sl_NM, sl_NB, mc = 0, MS, nu;
  unsigned i, j;
  unsigned char uc[POLY_Dmax * VERT_Nmax];
  VertexNumList V;
  EqList E;
  Long NF[POLY_Dmax][VERT_Nmax];
  Init_FInfoList(&L);

  if (F == nullptr) {
    printf("Input file %s not found\n", polyi);
    exit(1);
  }
  d = fgetc(F);
  if (d != 0) { /* for(i=0;i<d;i++) fgetc(F); */
    fprintf(stderr, "Error: Bin_2_ANF recursion depth %d must be 0\n", d);
    exit(1);
  }
  d = fgetc(F);
  L.nV = fgetc(F);
  L.nVmax = fgetc(F);
  L.NUCmax = fgetc(F);
  list_num = fgetUI(F);
  L.nNF = fgetUI(F);
  L.nSM = fgetUI(F);
  L.nNM = fgetUI(F);
  L.NB = fgetUI(F);
  sl_nNF = fgetUI(F);
  sl_SM = fgetUI(F);
  sl_NM = fgetUI(F);
  sl_NB = fgetUI(F);

  for (i = 0; i < L.nV; i++) {
    v = fgetc(F);
    L.nNUC[v] = fgetc(F); /* read #nuc's per #Vert */
    for (j = 0; j < L.nNUC[v]; j++) {
      L.NFnum[v][nu = fgetc(F)] = fgetUI(F); /* read nuc and #NF(v,nu)*/
      tNF += L.NFnum[v][nu];
      tNB += L.NFnum[v][nu] * (Along)nu;
    }
  }
  if ((unsigned int)(tNB - L.NB) != 0) {
    fputs("Error: Bin_2_ANF byte count mismatch\n", stderr);
    exit(1);
  }
  L.NB = tNB;
  if (tNF != L.nNF) {
    fprintf(stderr, "Error: Bin_2_ANF NF total mismatch tNF=%lld L.nNF=%lld\n",
            (long long)tNF, (long long)L.nNF);
    exit(1);
  }

  for (v = d + 1; v <= L.nVmax; v++)
    if (L.nNUC[v]) /* write  honest polys */
    {
      int I, J;
      for (nu = 1; nu <= L.NUCmax; nu++)
        for (j = 0; j < L.NFnum[v][nu]; j++) {
          for (s = 0; s < nu; s++)
            uc[s] = fgetc(F);
          UCnf_2_ANF(&d, &v, &nu, uc, NF, &MS);
          MS %= 4;
          _P->n = d;
          _P->np = v;
          for (I = 0; I < v; I++)
            for (J = 0; J < d; J++)
              _P->x[I][J] = NF[J][I];
          /* assert(Ref_Check(_P,&V,&E)); */
          if (MS != 1) {
            fprintf(stderr, "Error: Bin_2_ANF MS=%d, expected 1\n", MS);
            exit(1);
          }

          if (MS != 2) /* if(MS!=2) print NF */
            if (!max || Poly_Max_check(_P, &V, &E)) {
              mc++;
              Print_NF(out, &d, &v, NF);
            }
          if (MS > 1) /* if(MS>1); print Mirror */
            if (!max || Poly_Min_check(_P, &V, &E)) {
              mc++;
              Small_Make_Dual(_P, &V, &E);
              Make_Poly_NF(_P, &V, &E, NF);
              Print_NF(out, &d, &(V.nv), NF);
            }
        }
    }
  printf("np=%lld+%dsl  ", 2 * L.nNF - L.nSM - L.nNM,
         2 * sl_nNF - sl_SM - sl_NM);
  printf(/* write Finfo */
         "%dd  %dv<=%d n<=%d  %dnv  %lld %d %lld %lld  %d %d %d %d", d, L.nV,
         L.nVmax, L.NUCmax, list_num, L.nNF, L.nSM, L.nNM, L.NB, sl_nNF, sl_SM,
         sl_NM, sl_NB);
  if (max)
    printf("  r-max=%d", mc);
  puts("");
}

void Bin_2_ANF_DBsl(char *dbi, int max, int vf, int vt, PolyPointList *_P,
                    FILE *out) {
  FILE *F;
  FInfoList L;
  int d, v, nu, i, j, list_num, mc = 0, MS, sl_nNF, sl_SM, sl_NM, sl_NB,
                                tSM = 0, tNM = 0;
  std::vector<char> Ifn(1 + strlen(dbi) + File_Ext_NCmax);
  char *Ifx;
  Long NF[POLY_Dmax][VERT_Nmax];
  VertexNumList V;
  EqList E;
  strcpy(Ifn.data(), dbi);
  Ifx = &Ifn[strlen(dbi)];
  strcpy(Ifx, ".info");
  F = fopen(Ifn.data(), "r");
  if (F == nullptr) {
    puts("Info File not found");
    exit(1);
  }
  Init_FInfoList(&L); /* start reading the file */
  if (fscanf(F, "%d%d%d%d%d%lld%d%lld %lld %d%d%d%d", &d, &i, &j, &nu,
             &list_num, &L.nNF, &L.nSM, &L.nNM, &L.NB, &sl_nNF, &sl_SM, &sl_NM,
             &sl_NB) != 13) {
    fputs("Error: Bin_2_ANF_DBsl malformed info header\n", stderr);
    fclose(F);
    exit(1);
  }
  L.nV = i;
  L.nVmax = j;
  L.NUCmax = nu;
  if (sl_NB && (vf == 2) && (vt == VERT_Nmax - 1)) {
    strcpy(Ifx, ".sl");
    fclose(F);
    if (nullptr == (F = fopen(Ifn.data(), "rb"))) {
      printf("Open %s failed", Ifn.data());
      exit(1);
    }
  } else /* puts("no .sl file"); */
  {
    std::unique_ptr<DataBase> DB;
    Open_DB(dbi, &DB, 0);
    DataBase *DBRaw = DB.get();
    if ((DBRaw->v < vf) || (DBRaw->nVmax > vt))
      DB_fromVF_toVT(DBRaw, vf, vt); /*  read only vf <= v <= vt :: */
    for (i = 0; Read_H_poly_from_DB(DBRaw, _P); i++) {
      int c, l, p = _P->np - 1, off = _P->x[p][0];
      if (off)
        for (l = 0; l < d; l++)
          for (c = l; c <= p; c++)
            _P->x[c][l] -= off;
      /* if(off) Print_PPL(_P,"off"); else */
      Print_PPL(_P, "");
    }
    printf("#poly=%d\n", i);
    Close_DB(DBRaw);
  }
  for (i = 0; i < sl_nNF; i++) {
    int I, J;
    unsigned char uc[NUC_Nmax];
    int ch = fgetc(F);
    if (ch == EOF) {
      fprintf(stderr,
              "Error: Bin_2_ANF_DBsl unexpected EOF reading vertex count\n");
      exit(1);
    }
    v = ch;
    if (v > VERT_Nmax) {
      fprintf(stderr, "Error: Bin_2_ANF_DBsl vertex count %d out of range\n",
              v);
      exit(1);
    }
    ch = fgetc(F);
    if (ch == EOF) {
      fprintf(stderr,
              "Error: Bin_2_ANF_DBsl unexpected EOF reading nuc count\n");
      exit(1);
    }
    nu = ch;
    AuxGet_uc(F, &nu, uc);
    if (ferror(F)) {
      fprintf(stderr, "Error: Bin_2_ANF_DBsl sublattice read error\n");
      exit(1);
    }
    if ((*uc % 4) == 0)
      tSM++;
    else if ((*uc % 4) < 3)
      tNM++;

    UCnf_2_ANF(&d, &v, &nu, uc, NF, &MS);
    MS %= 4;
    _P->n = d;
    _P->np = v;
    for (I = 0; I < v; I++)
      for (J = 0; J < d; J++)
        _P->x[I][J] = NF[J][I];
    if (MS != 2) /* if(MS!=2) print NF */
      if (!max || Poly_Max_check(_P, &V, &E)) {
        mc++;
        Print_NF(out, &d, &v, NF);
      }
    if (MS > 1) /* if(MS>1); print Mirror */
      if (!max || Poly_Min_check(_P, &V, &E)) {
        mc++;
        Small_Make_Dual(_P, &V, &E);
        Make_Poly_NF(_P, &V, &E, NF);
        Print_NF(out, &d, &(V.nv), NF);
      }
  }
  printf("np=%lld+%dsl  ", 2 * L.nNF - L.nSM - L.nNM,
         2 * sl_nNF - sl_SM - sl_NM);
  printf(                                            /* write Finfo */
         "%dd  %dv<=%d n<=%d  %dnv ... %d %d %d %d", /* %d %d %d %lld */
         d, L.nV, L.nVmax, L.NUCmax, list_num,
         /* L.nNF,L.nSM,L.nNM,L.NB, */ sl_nNF, sl_SM, sl_NM, sl_NB);
  if (max)
    printf("  r-max=%d", mc);
  puts("");
}

void Gen_Bin_2_ascii(char *pi, char *dbi, int max, int vf, int vt,
                     PolyPointList *P, FILE *out) {
  if (*pi)
    Bin_2_ANF(pi, max, P, out);
  else if (*dbi)
    Bin_2_ANF_DBsl(dbi, max, vf, vt, P, out);
  else
    puts("With -B[2A] you have to specify input via -pi or -di");
}