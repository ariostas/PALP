/*  ======================================================================  */
/*  ==========                                                  ==========  */
/*  ==========                  C L A S S . C                   ==========  */
/*  ==========                                                  ==========  */
/*  ==========               A p r i l   2 0 0 2                ==========  */
/*  ==========                                                  ==========  */
/*  ======================================================================  */

/*      Coord   get coordinates by reading them or converting weight input
 *      Rat     rational functions
 *      Vertex  computes Vertices and Faces
 *      Subpoly subpolys, sublattices, minimality
 *      Polynf  normal form and symmetries
 *      SubAdd  sorting
 *      Subdb   data base
 */

#include "Global.h"
#include "Subpoly.h"

#if ( POLY_Dmax * POINT_Nmax  > 83400000 )
#error		decrease POLY_Dmax or/and POINT_Nmax for compiling class
#endif

FILE *inFILE, *outFILE;

void PrintExtOptions(void){puts("Extended/experimental options:");
     puts("           -d1 -d2 [-po]      combined mirror info (projected"); 
     exit(0);}
void LocalPrintUsage(char *c, char hc){
printf("This is  `%s', a program for classifying reflexive polytopes\n",c);
while (hc!='e'){
if (hc=='g'){
printf("Usage:     %s  [options] [ascii-input-file [ascii-output-file]]\n",c);
printf("Options:   -h          print this information\n           ");
puts("-f or -     use as filter; otherwise parameters denote I/O files");
puts("           -m*         various types of minimality checks (* ... lvra)");
puts("           -p* NAME    specification of a binary I/O file (* ... ioas)");
puts("           -d* NAME    specification of a binary I/O database (DB) (* ... ios)");
puts("           -r          recover: file=po-file.aux, use same pi-file");
puts("           -o[#]       original lattice [omit up to # points] only");
puts(
"           -s*         subpolytopes on various sublattices (* ... vphmbq)");
puts("           -k          keep some of the vertices");
puts("           -c, -C      check consistency of binary file or DB");
puts("           -M[M]       print missing mirrors to ascii-output");
puts("           -a[2b], -A  create binary file from ascii-input");
puts("           -b[2a], -B  ascii-output from binary file or DB");
puts("           -H*         applications related to Hodge number DBs (* ...cstfe)");
}
else if (hc=='m'){
printf("`%s -m*' reads a list of combined weight systems as ascii-input\n",c);
puts("              and writes the sublist with a particular property, ");
puts("              possibly with extra information:");
printf("`%s -ml' returns only lp-minimal CWS\n",c);
printf("`%s -mv' returns only v-minimal CWS\n",c);
printf("`%s -mr' returns only r-minimal CWS\n",c);
printf(
"`%s -ma' returns the CWS that determine reflexive polytopes, with\n",c);
puts(
"              information on the above properties and the `span property'.");
}
else if (hc=='p'){
puts("With `-p*' you specify a binary I/O file that encodes a sorted list");
puts("           of normal forms of polytopes. In particular, with");
puts("-piNAME or -pi NAME     input and with");
puts("-poNAME or -po NAME     output is specified.");
puts("-paNAME or -pa NAME     specifies a list that should be added to another ");
puts("         list given as a binary file (with -pi) or as a database (with -di).");
puts("-psNAME or -ps NAME     specifies a list that should be subtracted from");
puts("         another list given as a binary file (with -pi) .");
puts("With -pa and -ps you have to specify output via -po or -do.");
}
else if (hc=='d'){
puts("With `-d*' you specify a binary I/O database that encodes a sorted");
puts("           list of normal forms of polytopes. In particular, with");
puts("-diNAME or -di NAME     input and with");
puts("-doNAME or -do NAME     output is specified.");
puts("-dsNAME or -ds NAME     specifies a list in database format that should");
puts(
"         be subtracted from another list given as a binary file (with -pi),");
puts("         with an output file specified via -po.");
}
else if (hc=='r'){
printf("As %s sometimes requires very long running times, intermediate\n",c);
puts("results are regularly written to a file <out-file>.aux. If such a");
printf("file exists, `%s -r -po<out-file>' can be used to recover an\n",c);
printf("unfinished but terminated run of `%s -po<out-file>'.\n",c);
puts("Possible input files or databases should be identical.");
}
else if (hc=='o'){
printf("In normal mode `%s' determines reflexive subpolytopes both on\n",c);
printf("the original lattice and on sublattices. With `%s -o' you can\n",c);
puts("restrict to polytopes on the original lattice only. If you also");
puts("specify a number # via `-o#', then only the polytopes obtained by");
puts("omitting # or less lattice points are determined.");
puts("For `-o0' the recursion breaks at any reflexive polytope.");
printf("For `-oc' complete (including sublattice)");
puts(" by ignoring input polytope in list.");
}
else if (hc=='s'){
printf("`%s -s* [-di<input-db>] [-mr] -po<out-file>' polytopes on ",c);  
puts("sublattices:");  
printf("`%s -sh [-di<input-db>]' finds Calabi-Yau hypersurfaces that are free"
,c);
printf("\n             free quotients (i.e. points on codim>1 faces of ");
printf("the dual\n             polytope do not span the N-lattice). Input ");
printf("can be ascii or DB.\n");
printf(
"`%s -sp [-di<input-db>]' same as `-sh' except that it is checked whether",
  c);
puts("\n             all points of the dual span the N lattice.");
printf("`%s -sv [-di<input-db>] [-mr] -po<out-file>' serves to determine\n",c);
printf("             on which sublattices of the original lattice a given ");
printf("polytope is\n             still a lattice polytope. Input can be "); 
printf("ascii or database. In the\n             former case all sublattice ");
printf("polytopes are determined and in the\n             latter case only ");
printf("those not yet in the database. With the option\n             `-mr' ");
printf("the result, which is written to <out-file>, is\n             ");
puts("restricted to r-maximal polytopes.");
printf("`%s -sm [-di<input-db>] [-mr] -po<out-file>' same as `-sv' but ",c);
printf("now all\n             reflexive polytopes that have the same pairing");
puts(" matrix between\n             vertices and equations are constructed");
puts(
"In addition there are the somewhat experimental options -sb, -sq for dim=4.");
puts("They are similar to -sh, -sp, with the relationship summarized by");
puts("	-sh ... generated by codim>1 points (omit IPs of facets)");
puts("	-sp ... generated by all points");
puts("	-sb ... generated by dim<=1 (edges), print if rank=2	");
puts("	-sq ... generated by vertices,       print if rank=3	");
}
else if (hc=='k'){
printf(
"`%s -k* [-di<input-db>] -po<out-file>' gives you a list of the vertices \n",
c);
puts("   of the input polytope and asks which of them should be kept;");
puts(
"   all reflexive subpolytopes containing the kept vertices are determined.");}
else if (hc=='c'){
printf(
"`%s -c' (or `-C') checks the consistency of a binary file or database\n",c);
puts("specified via `-pi' or `-di'."); 
puts("`-C' results in a more detailed output than `-c'");
}
else if (hc=='M'){
printf("`%s -M[M]' looks for polytopes in a list specified by `-pi' or\n",c);
puts("`-di' whose mirrors are not in the same list. The resulting `missing");
puts("mirrors' are written in ascii format.");
}
else if (hc=='a'){
printf(
"`%s -a[2b] -po<out-file>' converts ascii-input of reflexive polytopes\n",c);
puts("to binary file format.");
puts("If an input file or database is specified via `-pi' or `-di', only");
puts("the polytopes not in one of these lists are written to <out-file>.");
puts("If an ascii output file is explicitly specified, weights corresponding");
puts("to new polytopes are written to that file.");
puts("-A[2B] converts non-reflexive ascii input to binary file format.");
}
else if (hc=='b'){
printf("`%s -b[2a]' converts binary input to a list of normal forms in\n",c);
puts("ascii format. For file input (specified via `-pi') all normal forms of");
puts("polytopes on original lattices are displayed, but for database input");
puts("(`-di') the normal forms of the sublattice polytopes in the database");
puts("are shown. If no sublattice polytopes are left, then all polytopes ");
puts("in the database are displayed.");
puts("For non-reflexive binary input, -B[2A] should be used for conversion");
puts("to ascii format.");
}
else if (hc=='H'){
puts(
"Options of the type `-H*' are used for handling Hodge number databases");
puts("and work only for polytopes of dimension four. In particular,");
printf("`%s -Hc [-vf#] [-vt#] -di<DB> -do<Hodge-DB>' calculates the\n",c);
puts("   Hodge numbers of the polytopes in DB and creates a Hodge number");
puts("   database whose files correspond to fixed vertex numbers and Euler");
puts("   numbers. If -vf and/or -vt are specified, only the Hodge numbers");
puts("   of the polytopes whose vertex numbers are in the corresponding");
puts("   (from/to) range are determined.");
printf(
"`%s -Hs -di<Hodge-DB> [-do<Hodge-DB>]' sorts a Hodge-DB resulting\n",c);
puts("   from -Hc to one consisting of files of fixed Hodge number pairs");
puts("   (with the same name if -do is omitted),");
printf("`%s  -Ht -di<Hodge-DB>'  tests a Hodge-DB for consistency and\n",c);
printf("`%s  -Hf -pi<Hodge-DB-file>' serves for testing a Hodge-DB-file.\n",c);
printf("`%s  -He<search-string> -di<Hodge-DB>' extracts data on \n",c);
puts("   polytopes (in ascii) from a Hodge-DB. A search string may take"); 
puts("   the form `E#H#:#M#V#N#F#L#', where the #'s denote numbers:"); 
puts("      E...Euler characteristic, H#:#...Hodge numbers h11,h12,"); 
puts("      M/V/N/F...numbers of points/vertices/dual points/facets,"); 
puts("      L...Limit on the total number of polytopes displayed."); 
puts("   The ordering is inessential and if a value isn't specified the"); 
puts("   corresponding symbol may be omitted. For example, `-He:1'"); 
puts("   leads to a search for all polytopes with h12=1.");
puts("   Unless at least one of h11, h12, E is specified, the search will");
puts("   take quite long.");
}
else if (hc=='I'){
printf("There are three basic types of I/O structures for %s:\n",c);
puts("ascii files, binary files and binary databases.");
puts("   Binary files and databases always encode ordered lists of normal");
puts("forms of polytopes, and any such structure created by some application");
printf("of %s may be used as input for some other application of %s.\n",c,c);
puts(
"   A database consists of various files NAME.<extension> and is specified");
puts("via `-d* NAME'. It contains one ascii file and several binary files.");
puts("   Ascii input should always correspond to a list of polytopes given");
puts("either by combined weight systems or by lists of lattice points.");
puts("Weight input is specified by a single line of the form");
printf("   d1 w11 w12 ... d2 w21 w22 ... [comments ignored by %s]\n",c);
puts("with sum_j wij=di for every i.");
puts("Lattice point input is specified by a line of the form");
printf("   #colums #lines [comments ignored by %s]\n",c);
puts("followed by (#lines) lines each of which has (#colums) integers such");
puts("that the resulting matrix encodes the coordinates of the polytope");
puts("with lattice points given either as row or as column vectors.");
puts("Sometimes ascii output may also be used as input.");
}
puts("");
puts(
"Type one of [m,p,d,r,o,s,c,M,a,b,H] for help on options,");
printf(
"`g' for general help, `I' for general information on I/O or `e' to exit: ");
scanf("%s",&hc);
puts("");
}
}


int  main (int narg, char* fn[])
{ int n=0, FilterFlag=0, oFlag=0, cFlag=0, rFlag=0, abFlag=0, kFlag=0,
    vf=2, vt=VERT_Nmax-1;
  char Blank=0, *dbin=&Blank, *dbsub=&Blank, *dbout=dbin, *x_string=&Blank,
    *polyi=dbin, *polya=dbin, *polys=dbin, *polyo=dbin, mFlag=0, HFlag=0,
    sFlag=0; static CWS W; PolyPointList *_P;
  if(narg==1) {
    printf("For help type `%s -h'\n", fn[0]);
    exit(0);}
  _P = (PolyPointList *) malloc(sizeof (PolyPointList));
  if(_P==NULL) {puts("Unable to allocate space for _P"); exit(0);}

  while(narg > ++n)
    if(fn[n][0]!='-') break;
    else switch(fn[n][1])
      { case 'h':               LocalPrintUsage(fn[0],'g');     exit(0);
        case 'f': case  0 :     FilterFlag=1;                   break;
        case 'm':               mFlag=fn[n][2];                 break;
        case 's':               sFlag=fn[n][2];                 break;
        case 'c':               cFlag=1;                        break;
        case 'C':               cFlag=2;                        break;
        case 'M':               cFlag=-1;                       break;
        case 'r':               rFlag=1;                        break;
        case 'a':               abFlag=1;                       break;
        case 'b':               abFlag=-1;                      break;
        case 'A':               abFlag=2;                       break;
        case 'B':               abFlag=-2;                      break;
        case 'k':               kFlag=1;                        break;
        case 'H': {
#if (POLY_Dmax != 4)
	  puts("For using Hodge-DB-routines set POLY_Dmax=4!"); exit(0);
#endif
	  HFlag=fn[n][2]; 
	  if(HFlag=='e') x_string = (fn[n][3]) ? &fn[n][3] : fn[++n];  }
                                                                break;
        case 'p':
        {        if(fn[n][2]=='i') polyi = (fn[n][3]) ? &fn[n][3] : fn[++n];
            else if(fn[n][2]=='a') polya = (fn[n][3]) ? &fn[n][3] : fn[++n];
            else if(fn[n][2]=='s') polys = (fn[n][3]) ? &fn[n][3] : fn[++n];
            else if(fn[n][2]=='o') polyo = (fn[n][3]) ? &fn[n][3] : fn[++n];
            else { LocalPrintUsage(fn[0],'g'); exit(0); }
        }                                                       break;
        case 'd':
        {        if(fn[n][2]=='i') dbin  = (fn[n][3]) ? &fn[n][3] : fn[++n];
            else if(fn[n][2]=='s') dbsub = (fn[n][3]) ? &fn[n][3] : fn[++n];
            else if(fn[n][2]=='o') dbout = (fn[n][3]) ? &fn[n][3] : fn[++n];
            else { LocalPrintUsage(fn[0], 'g'); exit(0); }
        }                                                       break;
        case 'o':       oFlag=-1;             /* original lattice only */
                if(fn[n][2])
		  if(fn[n][2]=='c'){oFlag=-3;puts("complete data");}
		  else {assert(('0'<=fn[n][2])&&(fn[n][2]<='9'));
		  oFlag=atoi(&fn[n][2]);
		  if (!oFlag) {oFlag=-2; puts("break recursion at RPs");}
		  else printf("rec-dep<=%d\n",oFlag);}
                else printf("original lattices only\n");
		                                                break;
        case 'v': {
          if(fn[n][2]=='f') vf=atoi((fn[n][3]) ? &fn[n][3] : fn[++n]);
          if(fn[n][2]=='t') vt=atoi((fn[n][3]) ? &fn[n][3] : fn[++n]);}
                                                                break;
	case 'x':	PrintExtOptions();			break;
        default:        printf("Unknown flag %s !!\n",fn[n]); exit(0);
      }
  n--;

  if(FilterFlag)      { inFILE=NULL; outFILE=stdout;     }
  else
    {   if (narg > ++n)  inFILE=fopen(fn[n],"r");  else inFILE=stdin;
        if (inFILE==NULL){printf("Input file %s not found!\n",fn[n]);exit(0);}
        if (narg > ++n) outFILE=fopen(fn[n],"w");  else outFILE=stdout;
    }

       if(sFlag)          VPHM_Sublat_Polys(sFlag,mFlag,dbin,polyi,polyo,_P);
  else if(abFlag==1)	  Ascii_to_Binary(&W,_P,dbin,polyi,polyo);
  else if(abFlag==-1)	  Bin_2_ascii(polyi,dbin,(mFlag=='r'),vf,vt,_P);
  else if(abFlag==2)      Gen_Ascii_to_Binary(&W,_P,dbin,polyi,polyo);
  else if(abFlag==-2)     Gen_Bin_2_ascii(polyi,dbin,(mFlag=='r'),vf,vt,_P);
  else if(cFlag)          Check_NF_Order(polyi,dbin,cFlag,_P);
  else if(mFlag=='a')     while(Read_CWS_PP(&W,_P)) Overall_check(&W,_P);
  else if(mFlag=='r')     while(Read_CWS_PP(&W,_P)) Max_check(&W,_P);
  else if(mFlag=='v')     while(Read_CWS_PP(&W,_P)) DPvircheck(&W,_P);
  else if(mFlag=='l')     while(Read_CWS_PP(&W,_P)) DPircheck(&W,_P);
#if (POLY_Dmax < 6)
  else if(HFlag=='c')     DB_to_Hodge(dbin, dbout, vf, vt,_P);
  else if(HFlag=='s')     Sort_Hodge(dbin, dbout);
  else if(HFlag=='f')     Test_Hodge_file(polyi,_P);
  else if(HFlag=='t')     Test_Hodge_db(dbin);
  else if(HFlag=='e')     Extract_from_Hodge_db(dbin,x_string,_P);
#endif
  else if(*dbin&&!*polyo) Add_Polya_2_DBi(dbin,polya,dbout);
  else if(*dbout)         Polyi_2_DBo(polyi,dbout);
  else if(*polya)         Add_Polya_2_Polyi(polyi,polya,polyo);
  else if(*polys||*dbsub) Reduce_Aux_File(polyi,polys,dbsub,polyo);
  else Do_the_Classification(&W,_P, /* fn[0], */ oFlag,rFlag,kFlag,
				polyi,polyo,dbin);
  return 0;
}
