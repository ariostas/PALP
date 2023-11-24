#include <limits.h>

#ifdef  __DECC				/* use local "/tmp" on clusters: */
#define	USE_TMP_DIR	(1)		/* write aux-files to "/tmp"     */ 
#else
#define	USE_TMP_DIR	(0)
#endif

#define	NUC_Nmax	256

/*	on 32-bit architectures the GNU C compiler requires flags:	     *
 *	-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE			     *
 *	moreover, fseek and ftell have to be replaced by fseeko and ftello   */
#ifdef	_LARGEFILE_SOURCE
#define FTELL ftello
#define FSEEK fseeko
#else
#define FTELL ftell
#define FSEEK fseek
#endif

/*	Target dimensions: 2^32 polytopes (unsigned!) ... 80GB/DB 	*/
/*	Along=64bit for polytopes except on bin-files::unsigned  	*/
/*	SL, self-mirror, NFnum[nv][nu]<2^32,  rec-depth?		*/

#define MAX_REC_DEPTH	16383	 /* present "dirty fix" workes till 128^2-1 */


/* large files with gcc -> -D_FILE_OFFSET_BITS=64 and ftell -> ftello */

#define	Along		long long	/* signed > addresses bytes DBpolys */
#define	UPint		unsigned	/* unsigned > #polys in RAM & aux   */

#define FORCE_SAVE_TIME         (28800) /* real time for Aux-IO in seconds */
#define GOOD_SAVE_TIME          (21600) /* same at multiple of 1000 polys */

#define	WRITE_DIM	4		/* File IO after weight if WD<= dim */
#define MIN_NEW		1		/* write file if MIN_NEW <= #newREF */
#define MIN_W_SAVE_TIME		(7200)  /* min real time for IO after poly */


#if	(POLY_Dmax < 5)

#define	WATCHREF	(100000)	/* print some info after X refs    */
#define	SAVE_INC	(1000000)	/* save PolyNFlist after X refs    */
#define	CperR_MAX	(32)		/* 4*9 is safe for CY (average=10) */
#define	BASE_MAX	(905)		/* 191->317  184->338  218->464 */
#define BLOCK_LENGTH    (64)	     /* fraction of data base stored in RAM */

#else

#define	WATCHREF	(100000)	    /* print some info after X refs */
#define	SAVE_INC	(4000000)	    /* save PolyNFlist after X refs */
#define	CperR_MAX	(40)	    	/*  average is 35 for PolyDim=5  */
#define	BASE_MAX    (1631723)	/* 1631721 1 903 37947 233103 543907 815860 */
#define BLOCK_LENGTH    (128)	     /* fraction of data base stored in RAM */

#endif

#define subl_int	LLong
#define NB_MAX		POLY_Dmax*VERT_Nmax	/* an uneducated guess */

#if(POLY_Dmax<5)
#define SL_Nmax		(65536)
#else
#define SL_Nmax		(300000)
#endif

/* ====	    this should be o.k. (only change if you know what you do)  ==== */

#define File_Ext_NCmax  5		 /* space for following FILE_EXTs */
#define	SAVE_FILE_EXT	".aux"		/* aux. O/I file for big allocation */
#define	TEMP_FILE_EXT	".tmp"	       /* undefine to directly overwrite aux */
#define	MOVE_SAVE_FILE	".bak"	      /* move SAVE-file after read  */
				 /* undef: risk data loss -> save disk space */
#undef	MOVE_SAVE_FILE		 /* undefine to directly overwrite SAVE file */
#undef	TEMP_FILE_EXT		 /* undefine to directly overwrite aux file */


#ifdef	MOVE_SAVE_FILE			/* consistency of these options */
#if	USE_TMP_DIR			/* is not yet implemented:      */
#error		Inconsistent options  USE_TMP_DIR  and  MOVE_SAVE_FILE !!!
#endif
#endif

typedef struct {
  Along nNF, nNM; int nSM,  /* #ref=2*nNF-nSelfMir.-nNoMir. */ 
    d, v, nu, p, sl_nNF, sl_SM, sl_NM, sl_NB, list_num,
    nV, nNUC[VERT_Nmax], nVmax, NUCmax, NFnum[VERT_Nmax][NUC_Nmax];
  Along  Fv_pos[VERT_Nmax][NUC_Nmax]; UPint RAM_pos[VERT_Nmax][NUC_Nmax];
  unsigned char *RAM_NF; long long NB;
  FILE *Finfo, *Fsl, *Fv[VERT_Nmax];   }                           DataBase;

typedef	struct {
  Along nNF, nNM; int nSM;	    /* #ref=2*nNF-nSelfMir.-nNoMir. */ 
  unsigned char nV, nNUC[VERT_Nmax+1], nVmax, NUCmax;
  unsigned int NFnum[VERT_Nmax+1][NUC_Nmax]; long long NB;
  unsigned char *NF[VERT_Nmax+1][NUC_Nmax],*NFli;}                FInfoList;

typedef struct {unsigned int n, c; }  /* below n in NFptr, @NewNF[c] */	PEnt;
typedef struct {int n; PEnt pe; }     /* below #n in PEnt[], PEnt    */	PPEnt;
	
typedef struct {
  /* flags and file names: */
  int SL, of, rf, kf, rd, b[POINT_Nmax];   /* orig/omit-flag, recov-flag, 
                                              keep-flag, rec.depth, branch */
  char *iname, *oname, *dbname;	     /* file names : &(Constant Strings) !! */ 
  /* Statistics: */	
  UPint nNF, nIP; unsigned int nSLNF,  hc, con;         /* #refhit #con.ref */
  int V, F, VN, FN, d, Nmin, Nmax, Xdif, Xnuc;             /* max values */
  time_t TIME, SAVE; clock_t CLOCK; long long IP_Time, NF_Time;
  /* Results */	
  Along NP, savedNP, NC; int nSLP;   		/* NC = # unit-encode */
  /* Lists */	
  DataBase DB; FInfoList Aux, In; PEnt *PE; PPEnt *PPE; 
  unsigned char *NewNF; 
  Along NewNB, ANB; /* allocate*/ 	int RemNB; /* SL-remove */
  int *SLp, SLN, PEN, PPEN, peNM, peSM, slNM, slSM; 
                                               /* no/self mirror's on NewNF */
                                                 		 }    NF_List;

/* 	NP (HP)=#(honest) refpolys         NNF=#nf's, NSM=#self mirror,  */
/*		NP==HRP+SLRP	HRP==2*HNF-HSM-HnoMirror
 *		HNF=PEN+PPEN+Aux.nNF	HSM=peSM+Aux.nSM    HNM=peNM+Aux.nNM
 *		SLRP=2*SLN-slNM-slSM                                      */
/* watch: #R=NP (nSLP) #C=NC #IP=nIP #NF=nNF (nSLNF) Nmin..Nmax vV fF dx=dX */
		
/*  ==========          	checks and warnigs		 ========== */

#if	( FORCE_SAVE_TIME <= MIN_W_SAVE_TIME )
#error	MIN_W_SAVE_TIME should be smaller than AUX-file save_times
#endif


void Make_ANF(PolyPointList *P,VertexNumList *V,       /* affine normal form */
	      EqList *E, Long ANF[POLY_Dmax][VERT_Nmax]);
void Gen_Ascii_to_Binary(CWS *W, PolyPointList *P, 
		     char *dbin, char *polyi, char *polyo);
void Gen_Bin_2_ascii(char *pi,char *dbi,int max,int vf,int vt,PolyPointList *);


/*  ==========          Functions from Subpoly.c                 ========== */
void DPircheck(CWS *_W, PolyPointList *_P);
void DPvircheck(CWS *_W, PolyPointList *_P);
void Max_check(CWS *_W, PolyPointList *_P);
void Overall_check(CWS *_W, PolyPointList *_P);
void Do_the_Classification(CWS *W, PolyPointList *P, /* char *fn, */
       int oFlag, int rFlag, int kFlag, char *polyi, char *polyo, char *dbin);
void Find_Sublat_Polys(char mFlag, char *dbi, char *polyi, char *slout, 
		       PolyPointList *_P);
void Ascii_to_Binary(CWS *W, PolyPointList *P, 
		     char *dbin, char *polyi, char *polyo);
int  Start_Find_Ref_Subpoly(PolyPointList *_P/*, NF_List *_NFL*/);
void uc_nf_to_P(PolyPointList *_P, int *MS, int *d, int *v, int *nuc, 
                unsigned char *uc);
void Make_All_Sublat(NF_List *_L, int n, int v, subl_int diag[POLY_Dmax], 
                     subl_int u[][VERT_Nmax], char *mFlag, PolyPointList *P);
int  Poly_Max_check(PolyPointList *_P, VertexNumList *_V, EqList *_E);
int  Poly_Min_check(PolyPointList *_P, VertexNumList *_V, EqList *_E);

/*  ==========        Functions from Subdb.c                     ========== */
void Init_DB(NF_List *_NFL);
void Check_NF_Order(char *polyi,char *polyo,int cFlag, PolyPointList *P);
void Add_Polya_2_DBi(char *dbi,char *polya,char *dbo);
void Polyi_2_DBo(char *polyi,char *dbo);
void Reduce_Aux_File(char *polyi,char *polys,char *dbsub,char *polyo);
void Bin_2_ascii(char *polyi,char *dbi,int max,int vf,int vt,PolyPointList *P);
int  Is_in_DB(int *nv, int *nuc, unsigned char *uc, NF_List *_L);

#if (POLY_Dmax <6)
  void DB_to_Hodge(char *dbin,char *dbout, int vfrom,int vto, PolyPointList *P);
  void Sort_Hodge(char *dbaux, char *dbout);
  void Test_Hodge_file(char *filename, PolyPointList *_P);
  void Test_Hodge_db(char *dbname);
  void Extract_from_Hodge_db(char *dname, char *x_string, PolyPointList *P);
#endif

void Open_DB(char *dbin, DataBase **DB, int info);
int  Read_H_poly_from_DB(DataBase *DB,PolyPointList *_P);
void Close_DB(DataBase *DB);
void VPHM_Sublat_Polys(char sFlag,char mr,char *dbin,char *polyi,char *polyo, 
		       PolyPointList *P);


/*  ==========        Functions from Subadd.c                    ========== */
void Add_Polya_2_Polyi(char *polyi,char *polya,char *polyo);
void Init_NF_List(NF_List *);
void Init_FInfoList(FInfoList *FI);
void Read_File_2_List(char *polyi,NF_List *_NFL);
void Write_List_2_File(char *polyo,NF_List *_NFL);
void Print_Weight_Info(CWS *_W, NF_List *_L);
void fputUI(unsigned int l,FILE *F);
void UCnf2vNF(int *d, int *v, int *nuc, unsigned char *uc,  	      /* IN */
	      Long NF[POLY_Dmax][VERT_Nmax], int *MS);	             /* OUT */
int  Add_NF_to_List(PolyPointList *_P, VertexNumList *_V, EqList *_F,
		    NF_List *_NFL);				/* 1 if new */
int  RIGHTminusLEFT(unsigned char *ucL, unsigned char *ucR, int *nuc);
unsigned int fgetUI(FILE *F);
void Test_ucNF(int *d, int *tnv, int *tnuc, unsigned char *tuc, 
	       PolyPointList *P);
int  InfoSize(int rd, int lists, FInfoList *FI);
int  Make_Poly_NF(PolyPointList *_P, VertexNumList *_V, EqList *_E,
		Long pNF[POLY_Dmax][VERT_Nmax]);	  /* 1 if reflexive */
void ANF_2_ucNF(PolyPointList *P, VertexNumList *V, EqList *E,	      /* IN */
     		 int *NV, int *nUC, unsigned char *UC);		     /* OUT */
void UCnf_2_ANF(int *d, int *v, int *nuc, unsigned char *uc,          /* IN */
              Long NF[POLY_Dmax][VERT_Nmax], int *MS);                /* OUT */

void Print_Expect(FInfoList *L);

/* ====    headers of aux-Routines for Add_Polya_2_DBi in Subadd.c ==== */

void Read_Bin_Info(FILE *F, int *d, unsigned *li, int *slNF, int *slSM, 
		int *slNM, Along *NewNB, FInfoList *FI);
void AuxGet_vn_uc(FILE *F,int *v, int *nu, unsigned char *uc);
void AuxGet_uc(FILE *F,int *nu, unsigned char *uc);
void AuxPut_hNF(FILE *F,int *v,int *nu,unsigned char *Huc,FInfoList *Io,
        int *slNF,int *slSM,int *slNM,int *slNB,unsigned char *ucSL,int *SLp);

#define min(a,b)	(((a) < (b)) ? (a) : (b))
#define max(a,b)	(((a) > (b)) ? (a) : (b))

