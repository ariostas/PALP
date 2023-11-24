#include "Global.h"
#include "Subpoly.h"


/*   VF_2_ucNF / UCnf2vNF compression/decompression assumes that vNF[0][0]==1
 *   (see UCnf2vNF: "{ off=NF[0][0]-1; ..."
 */

/* Write_List_2_File */

/*   BIN File Format:   		==[FileINFO|data]	MAX_REC_DEPTH
 *   [FileINFO]:: [RecDepth<128] or [128+rd/256,rd%256]		but: .<7 !!!
 *		  L->b[]		rd byte		<256 !!
 * 		  d #Vposs #Vmax NUCmax	4  byte		<256 !!
 *		  #li #NF #SM #NM #NB   5*4byte		<2^32!! but: NB%2^32
 *	   sublattice #N? #SM #NM #NB	4*4byte
 *	   	  #v_ #nu_  #v2 nu2 ...	2  byte		<256 !!
 *		      -> nu nNF		1+4byte	... 	li=sum_{#nu(v)}
 */


#define	TEST_UCnf	    /*  ====  check that decomp(compress)=id  ====  */

#if ( 2 * WATCHREF  > SAVE_INC )
#error		increase SAVE_INC / WATCHREF
#endif

#ifdef	  ADD_LIST_LENGTH		 /*  ... square root of SAVE_INC    */
#if ( 2 * ADD_LIST_LENGTH  > SAVE_INC )	 
#error		increase SAVE_INC / ADD_LIST_LENGTH
#endif
#else
int  IntSqrt(int q)                     /* sqrt(q) => r=1; r'=(q+r*r)/(2r); */
{    assert(q>0); if(q<4) return 1; else {	/* troubles: e.g. 9408 */
     long long r=(q+1)/2,n; while(r > (n=(q+r*r)/(2*r))) r=n; if(q<r*r)r--;
     	if((r*r<=q)&&(q<(r+1)*(r+1))) return (int)r;
	else {printf("Error in sqrt(%d)=%d\n",q,(int)n);exit(0);}} return 0;
}
#define	ADD_LIST_LENGTH		(IntSqrt(SAVE_INC))	
#endif

#undef	More_File_IO_Data

/* === POLY_NF on FILE:  np nv nf complete lines of upper_trian_matrix  === */
/* === POLY_NF on LIST:  np nv nf lines with v00 & v[ij] (i>j) omitted  === */

#define INCREMENTAL_TIME		/* these flags should stay "on" */
#define INCREMENTAL_WRITE		/* write no polys from pi-file */
#define	ACCEL_PEntComp			/* speed optimization flags   */
#define	USE_UNIT_ENCODE

/*  ==========	Auxiliary routines from Polynf.c                ==========  */
void Eval_Poly_NF(int *d,int *v,int *f, Long VM[POLY_Dmax][VERT_Nmax],
                Long VPM[VERT_Nmax][VERT_Nmax],                       /* in */
                Long pNF[POLY_Dmax][VERT_Nmax], int t);              /* out */
int  Init_rVM_VPM(PolyPointList *P, VertexNumList *_V,EqList *_F,     /* in */
                int *d,int *v,int *f, Long VM[POLY_Dmax][VERT_Nmax], /* out */
                Long VPM[VERT_Nmax][VERT_Nmax]);        /* return reflexive */


typedef struct {int base[VERT_Nmax+1][NB_MAX], nuc[VERT_Nmax+1][NB_MAX],
					v[VERT_Nmax+1];}	 Base_List;

void VF_2_ucNF(PolyPointList *P, VertexNumList *V, EqList *E,	      /* IN */
     		 int *NV, int *nUC, unsigned char *UC);		     /* OUT */
void Print_Statistics(NF_List *);
void Init_BaseList(Base_List **BL,int *d);    /* malloc + init.; BL=&(list) */
void Insert_PPent_into_Pent(NF_List *S);
void Read_In_File(NF_List *);
void Read_Aux_File(NF_List *);

/*   =========		End of Headers and TypeDefs		==========  */
/*   =====================================================================  */

void Init_FInfoList(FInfoList *FI)
{    int i, j;
     FI->nNF= FI->nSM= FI->nNM= FI->NB=0; FI->nVmax= FI->NUCmax= FI->nV=0; 
     for(i=0;i<=VERT_Nmax;i++) FI->nNUC[i]=0;
     FI->NFli=NULL;
     for(i=0;i<=VERT_Nmax;i++) for(j=0;j<NUC_Nmax;j++) FI->NFnum[i][j]=0;
}

void Init_New_List(NF_List *S)
{    S->ANB = (SAVE_INC+SL_Nmax) * CperR_MAX;  S->NewNB = S->NP = S->nSLP = 0;
     S->RemNB=S->PEN=S->PPEN=S->SLN=0; S->peNM=S->peSM=S->slNM=S->slSM=0;
     S->PE  = (PEnt  *) malloc( (SAVE_INC+SL_Nmax) * sizeof(PEnt)  );
     S->PPE = (PPEnt *) malloc( ADD_LIST_LENGTH * sizeof(PPEnt) );
     S->SLp = (int   *) malloc( SL_Nmax * sizeof(int)   );
     S->NewNF = (unsigned char *) malloc( S->ANB * sizeof(char) ); S->NC = 0;
     assert((S->PE!=NULL)&&(S->PPE!=NULL)&&(S->SLp!=NULL)&&(S->NewNF!=NULL));
#ifdef  __DECC
        /* printf("NULL=%p S->PE=%p S->NewNF=%p\n",NULL,S->PE,S->NewNF); */
#endif
}

NF_List *AuxNFLptr=NULL;		/* dirty trick for Xmin Xmax Xdif */
void Init_NF_List(NF_List *L)
{    
     L->TIME= L->SAVE= time(NULL); fputs(ctime(&L->TIME),stdout); L->CLOCK= clock();
     L->IP_Time = L->NF_Time = 0;  L->d = L->nNF = L->nIP = L->nSLNF = 0;
     
     Init_FInfoList(&L->In);	if(*L->iname) Read_In_File(L);	
     if(*L->dbname) Init_DB(L);
     Init_FInfoList(&L->Aux);	Init_New_List(L);  L->savedNP= 0; AuxNFLptr=L;
     if(L->rf) Read_Aux_File(L); else L->rd=0;
}

void PrintNumbers(NF_List *S)
{    printf(
	"NP=%lld  nSLP=%d  H=2*%lld - %lldnm - %dsm  SL=2*%d - %dnm - %dsm\n",
	S->NP, S->nSLP, S->Aux.nNF+S->PEN+S->PPEN,
	S->Aux.nNM+S->peNM, S->peSM+S->Aux.nSM, S->SLN,
	S->slNM, S->slSM); fflush(stdout);
}
void Test_SLnbMS(NF_List *S)
{    int i, tNB=0, tNM=0, tSM=0, testNB=!(S->PEN||S->PPEN);
     if(testNB)assert(S->RemNB==0);
     for(i=0;i<S->SLN;i++)
     {	unsigned char *C=&S->NewNF[S->SLp[i]]; int ms=C[2]%4;
	if(C[0]>VERT_Nmax) {printf("v[%d]=%d\n",i,C[0]);exit(0);}
        tNB+=C[1]+2; if(ms) {if(ms<3)tNM++;} else tSM++;
     }	if((S->slNM!=tNM) || (S->slSM!=tSM) || (testNB&&(S->NewNB!=tNB)))
     {	printf("Test_SLnbMS  NM: %d=%d  SM: %d=%d  NM: %lld=%d",S->slNM,tNM,
		S->slSM,tSM,S->NewNB,tNB); exit(0);
     }	assert(S->nSLP==2*S->SLN-tSM-tNM);
}
void Test_ucNF(int *d, int *tnv, int *tnuc, unsigned char *tuc, 
	       PolyPointList *_P);
void Test_NF_List(NF_List *S, PolyPointList *_P)
{    int i; 			if(S->PPEN) Insert_PPent_into_Pent(S);
     printf("test PEN=%d  ",S->PEN);
     for(i=0;i<S->PEN;i++)
     {	unsigned char *C=&S->NewNF[S->PE[i].c], *uc=&C[2];
	int nv=C[0], nuc=C[1]; Test_ucNF(&S->d,&nv,&nuc,uc,_P);
        if(i)
	{   unsigned char *oldC=&S->NewNF[S->PE[i-1].c], *olduc=&C[2];
	    unsigned int oldn=S->PE[i-1].n;
	    if(oldC[0]>nv) { printf("NV failed at i=%d\n",i);exit(0); }
	    if(oldC[0]==nv) 
	    { 	assert(oldC[1]<=nuc);
     		if(oldC[1]==nuc) if(RIGHTminusLEFT(olduc,uc,&nuc)<0)
		{    printf("failed at i=%d\n",i);exit(0); }
		if(oldC[1]==nuc) if(oldn>S->PE[i].n)
		{    printf("oldn>S->PE[%d].n\n",i);exit(0); }
	    }
	}
     }					/* 	printf("test SLN=%d ",S->SLN);
     for(i=0;i<S->SLN;i++)
     {	unsigned char *C=&S->NewNF[S->SLp[i]], *uc=&C[2];
	int nv=C[0], nuc=C[1]; Test_ucNF(&S->d,&nv,&nuc,uc); 
     }							puts("Test o.k."); */
}
void Test_PPEN(NF_List *S)
{    int i; 		/* printf("test PPEN=%d  ",S->PPEN); */
     for(i=0;i<S->PPEN;i++)
     {	unsigned char *C=&S->NewNF[S->PPE[i].pe.c], *uc=&C[2];
	int nv=C[0], nuc=C[1]; /**  / Test_ucNF(&S->d,&nv,&nuc,uc); /  **/
        if(i)
	{   unsigned char *oldC=&S->NewNF[S->PPE[i-1].pe.c], *olduc=&oldC[2];
	    if(oldC[0]>nv) { printf("NV failed at i=%d:  v%d=%d  v%d=%d\n",
		i,i,*C,i-1,*oldC);exit(0); }
	    if(oldC[0]==nv) 
	    { 	assert(oldC[1]<=nuc);
     		if(oldC[1]==nuc) if(RIGHTminusLEFT(olduc,uc,&nuc)<0)
		{    printf("failed at i=%d\n",i);exit(0);
		}
	    }
	}
     }	
}
int  InfoSize(int rd, int lists, FInfoList *FI)
{    return (rd/128+rd+1 + 4+ 2*FI->nV+ lists) + sizeof(int) * (9 + lists); 
}
unsigned int fgetUI(FILE *F)       /* read unsigned int from bin file */
{    unsigned char A,B,C,D;             /* L=D+256*(C+256*(B+256*A)); */ 
     fscanf(F,"%c%c%c%c",&A,&B,&C,&D);  
     return (((unsigned int)A  * 256 + (unsigned int)B) * 256 + 
              (unsigned int)C) * 256 + (unsigned int)D;
}
void Read_Bin_Info(FILE *F, int *d, unsigned *li, int *SLN, int *slSM, 
		int *slNM, Along *NewSLnb, FInfoList *FI)
{    int i,j,v; unsigned tli=0; Along tNF=0,tNB=0; *d=fgetc(F); 
     assert(*d <= POLY_Dmax);
     FI->nV=fgetc(F); FI->nVmax=fgetc(F); FI->NUCmax=fgetc(F);
	*li=fgetUI(F); FI->nNF=fgetUI(F); FI->nSM=fgetUI(F); 
	FI->nNM=fgetUI(F); FI->NB=fgetUI(F); 
	*SLN=fgetUI(F); *slSM=fgetUI(F); *slNM=fgetUI(F); *NewSLnb=fgetUI(F); 
     for(i=0;i<FI->nV;i++)
     {	v=fgetc(F); tli+=(FI->nNUC[v]=fgetc(F)); 
	for(j=0;j<FI->nNUC[v];j++) 
	{   unsigned nu, nnb=FI->NFnum[v][nu=fgetc(F)]=fgetUI(F);
	    tNF+=nnb; tNB+=nnb*nu;
	}
     }	assert(tNF==FI->nNF); assert(tli==*li);
     assert( 0 == (unsigned int) (tNB-FI->NB) ); FI->NB=tNB;
}
void Read_Honest_Poly(FILE *F,FInfoList *FI,NF_List *L)
{    int i, v, nu; unsigned li; Along tNF=0, pos=0;
     Init_FInfoList(FI); L->rd=fgetc(F); if(128<=(L->rd))
     {	assert((L->rd-=128)<7); L->rd = 128*L->rd + fgetc(F); /* DirtyFix rd */
     }
     for(i=0;i<L->rd;i++) L->b[i]=fgetc(F);
     Read_Bin_Info(F,&L->d,&li,&L->SLN,&L->slSM,&L->slNM,&L->NewNB,FI); 
	
     L->nSLP=2*L->SLN-L->slSM-L->slNM;	   L->PEN=L->PPEN=L->peNM=L->peSM=0;
     L->NP=L->savedNP=L->nSLP+2*FI->nNF-FI->nSM-FI->nNM;	L->RemNB=0;
     printf("%lld+%dsl %lldm+%ds %lldb",L->NP-L->nSLP,L->nSLP, 
			FI->nNF-FI->nSM-FI->nNM,FI->nSM, FI->NB+L->NewNB); 
	if(L->rd) printf(" rd=%d",L->rd);
     if(FI->NFli != NULL) printf("WARNing: NFli != NULL");
     fflush(stdout);
     FI->NFli = (unsigned char *) malloc( FI->NB * sizeof(char) );
     if(FI->NFli==NULL) {puts("Aux.NFli allocation failed");exit(0);}
     for(v=L->d+1;v<=FI->nVmax;v++) if(FI->nNUC[v]) 
     for(nu=1;nu<=FI->NUCmax;nu++) if(FI->NFnum[v][nu])
     {	unsigned j; FI->NF[v][nu]=&FI->NFli[pos]; tNF+=FI->NFnum[v][nu];
	for(j=0;j<FI->NFnum[v][nu];j++) for(i=0;i<nu;i++) 
	    FI->NFli[pos++]=fgetc(F); 
     }				assert(pos==FI->NB);
     /* printf("\nRead: tNF=%d  pos=%d\n",tNF,pos); 
	for(v=L->d+1;v<=FI->nVmax;v++) if(FI->nNUC[v]) 
	for(nu=1;nu<=FI->NUCmax;nu++) for(j=0;j<FI->NFnum[v][nu];j++)
	Test_ucNF(&L->d,&v,&nu,&FI->NF[v][nu][nu*j]); puts("o.k.");
      */
}
void Read_SubLat_Poly(FILE *F,NF_List *L)
{    int i, j, pos=0; for(i=0;i<L->SLN;i++)
     {	unsigned char *C=&L->NewNF[pos]; L->SLp[i]=pos;
	L->NewNF[pos++]=fgetc(F); L->NewNF[pos++]=fgetc(F);
        for(j=0;j<C[1];j++) L->NewNF[pos++]=fgetc(F);
     }							assert(pos==L->NewNB);
}
void Read_In_File(NF_List *S)
{    time_t Tstart=time(NULL); FILE *F=fopen(S->iname,"rb");     /* F=fopen */
     printf("Reading In-File %s: ",S->iname);fflush(stdout);
     if(F==NULL) {puts("Cannot open (read)!"); exit(0);}
     Read_Honest_Poly(F,&S->In,S); S->NP=S->nSLP=0; assert(S->rd==0);
     fclose(F); printf("  done (%ds)\n",(int) difftime(time(NULL),Tstart)); 
     fflush(stdout);
}
void Read_File_2_List(char *fn,NF_List *L)	/* ... like Read_Aux_File */
{    time_t Tstart=time(NULL); FILE *F=fopen(fn,"rb");
     printf("Reading %s: ",fn); 
     if(F==NULL) {puts("Cannot open (read)!"); exit(0);}
     Read_Honest_Poly(F,&L->Aux,L);
     Read_SubLat_Poly(F,L);
     printf("  done (%ds)\n",(int) difftime(time(NULL),Tstart)); 
     fclose(F); fflush(stdout);
}
void Read_Aux_File(NF_List *L)
{    time_t Tstart=time(NULL); FILE *F; /* F=fopen */
     int NCalloc=strlen(L->oname)+strlen(SAVE_FILE_EXT)+ (USE_TMP_DIR ? 6 : 1);
     char *auxfn = (char *) malloc(NCalloc);
     if(USE_TMP_DIR) {strcpy(auxfn,"/tmp/"); strcat(auxfn,L->oname);} 
     else strcpy(auxfn,L->oname);
     strcat(auxfn,SAVE_FILE_EXT);
     F=fopen(auxfn,"rb"); printf("Reading %s: ",auxfn); fflush(stdout);
     if(F==NULL) puts("No aux-file found!"); 
     else
     {	Read_Honest_Poly(F,&L->Aux,L); Read_SubLat_Poly(F,L);
     	printf("  done (%ds)\n",(int) difftime(time(NULL),Tstart));
	fclose(F); fflush(stdout); 
#ifdef	MOVE_SAVE_FILE			/* inconsistent with USE_TMP_DIR !! */
     	{char *Mfn=(char *) malloc(1+strlen(L->oname)+strlen(MOVE_SAVE_FILE));
	strcpy(Mfn,L->oname); strcat(Mfn,MOVE_SAVE_FILE);
	assert(!rename(auxfn,Mfn)); free(Mfn);}
#endif
     }	
     fflush(stdout); free(auxfn);
}
void fputUI(unsigned int l,FILE *F) 	/* write unsigned int to bin file */
{    unsigned char A,B,C,D;
     D=l % 256; l/=256; C=l % 256; l/=256; B=l % 256; l/=256;   A=l; 
     fprintf(F,"%c%c%c%c",A,B,C,D);
}
void TestMSbits(NF_List *S, PolyPointList *_P)
{    int i, v, nu, tc, peNF=0,peSM=0,peNM=0, 
	slNF=0, slSM=0,slNM=0,axSM=0;		UPint axNF=0,axNM=0; 
     for(i=0;i<S->PEN;i++) 
     {	unsigned char *C=&S->NewNF[S->PE[i].c]; int ms=C[2]%4; peNF++; 
	{ /* int nv=*C,nuc=C[1]; Test_ucNF(&S->d,&nv,&nuc,&C[2]);*/ }
	if(ms==0) (peSM)++; else if(ms!=3) (peNM)++;
     }
     for(i=0;i<S->PPEN;i++) 
     {	unsigned char *C=&S->NewNF[S->PPE[i].pe.c]; int ms=C[2]%4; peNF++; 
	{ /* int nv=*C,nuc=C[1]; Test_ucNF(&S->d,&nv,&nuc,&C[2]);*/ }
	if(ms==0) (peSM)++; else if(ms!=3) (peNM)++;
     }
     for(v=S->d+1;v<=S->Aux.nVmax;v++) if(S->Aux.nNUC[v])
     for(nu=1;nu<=S->Aux.NUCmax;nu++) if((tc=S->Aux.NFnum[v][nu]))
     for(i=0;i<tc;i++)
     {	unsigned char *C=&S->Aux.NF[v][nu][nu*i]; int ms=C[0]%4; axNF++;
	if(ms==0) (axSM)++; else if(ms!=3) (axNM)++;
     }
     printf("\nTEST: hNP=%lld  hNF=%d+%lld=%d+%d  nNM=%d+%d  nSM=%d+%d\n",
	S->NP-S->nSLP,S->PEN+S->PPEN,S->Aux.nNF,peNF,axNF,peNM,axNM,peSM,
	axSM); fflush(stdout);      
     assert(S->NP-S->nSLP == 2*peNF-peNM-peSM + 2*axNF-axNM-axSM);
     for(i=0;i<S->SLN;i++) 
     {	unsigned char *C=&S->NewNF[S->SLp[i]]; int ms=C[2]%4; slNF++; 
	{int nv=*C,nuc=C[1]; Test_ucNF(&S->d,&nv,&nuc,&C[2],_P); }
	if(ms==0) (slSM)++; else if(ms!=3) (slNM)++;
     }
     printf("  sl: slNP=%d  slNF=%d=%d  slNM=%d  slSM=%d\n\n",
	S->nSLP,S->SLN,slNF,slNM,slSM); fflush(stdout);      
     assert(S->nSLP == 2*slNF-slNM-slSM);
}
void AuxCalcNumbers(NF_List *S,FInfoList *AI)
{    int i,j; Init_FInfoList(AI); if(S->PPEN) Insert_PPent_into_Pent(S);
     for(i=0;i<S->PEN;i++) 
     {	unsigned char *C=&S->NewNF[S->PE[i].c]; int ms=C[2]%4; AI->nNF++; 
	++(AI->NFnum[C[0]][C[1]]); 
	if(ms==0) (AI->nSM)++; else if(ms!=3) (AI->nNM)++;
     }
     for(i=S->d+1;i<=VERT_Nmax;i++) for(j=0;j<NUC_Nmax;j++)if(AI->NFnum[i][j]) 
     {	AI->nNUC[i]++; if(j>AI->NUCmax)AI->NUCmax=j;
	AI->NB += j * AI->NFnum[i][j]; 
     }  						/* TestMSbits(S); */
     for(i=S->d+1;i<=VERT_Nmax;i++) if(AI->nNUC[i]) AI->nVmax=i;
}
/*   uchar=unsigned char	uint=unsigned int   
 *   NP=#Polys  NB=#Bytes  #files=#nv's with NP>0  #lists=#(nv,nuc) with NP>0
 *
 * uchar  rd  k_1 ... k_rd    			   // ... not in data base !! 
 *					   // big rd:   "128+rd/128" "rd%128"
 * uchar  dim  #files  nv_max  nuc_max
 * uint   		#lists  hNF  hSM  hNM  hNB  slNF  slSM  slNM  slNB
 * uchar  v1    #nuc's with v1
 * uchar  nuc1  uint #NF(v1,nuc1)  uchar nuc2  uint #NV(v1,nuc2)  ...
 * uchar  v2    #nuc's with v2
 * uchar  nuc1  uint #NF(v2,nuc1)  uchar nuc2  uint #NV(v2,nuc2)  ...
 * uchar  "all hNF honest nf's"
 * uchar  "all slNF sublattice {nv nuc nf[]}'s"
 */
void Write_Bin_File(FILE *F,NF_List *L)
{    unsigned int i,j,fi=0,li=0,v,nu,tc=0,pen=0, tnb;
     UPint NFnum[VERT_Nmax][NUC_Nmax];
     int n; unsigned char *C, nVmax,NUCmax;
     FInfoList AI, *_AI = (FInfoList *) &AI; 

     printf("%lld+%dsl %lldm+%ds %lldb",L->NP-L->nSLP,L->nSLP,
	L->Aux.nNF - L->Aux.nNM - L->Aux.nSM + L->PEN + L->PPEN - L->peNM 
	- L->peSM,		L->Aux.nSM + L->peSM,
	/* nNF :: L->Aux.nNF+L->PEN+L->PPEN+L->SLN, */
	L->Aux.NB+L->NewNB-L->RemNB-2*(L->PEN+L->PPEN));
     if(L->rd) printf(" rd=%d",L->rd); 
#ifdef	USE_UNIT_ENCODE
	if(L->NC>0) printf("  u%lld",L->NC*100 / (L->Aux.nNF+L->PEN+L->PPEN) );
#endif
     AuxCalcNumbers(L,_AI);
     Print_Expect(_AI);
     fflush(stdout);
     nVmax=max(L->Aux.nVmax,AI.nVmax); NUCmax=max(L->Aux.NUCmax,AI.NUCmax);
     for(i=L->d+1;i<=nVmax;i++)
     { int v_li=0;
       for(j=1;j<=NUCmax;j++)
	 if((NFnum[i][j]=L->Aux.NFnum[i][j]+AI.NFnum[i][j])) v_li++;
       li+=v_li;
       if(v_li) fi++;
     }
     assert(L->rd <= MAX_REC_DEPTH);
     if(L->rd < 128) fputc(L->rd,F);
     else {fputc(128+(L->rd)/128,F); fputc((L->rd)%128,F);} 
     for(n=0;n<L->rd;n++) fputc(L->b[n],F); 
     fputc(L->d,F); fputc(fi,F); fputc(nVmax,F); fputc(NUCmax,F); 
	fputUI(li,F); fputUI(L->Aux.nNF+AI.nNF,F); 
	fputUI(L->Aux.nSM+AI.nSM,F); fputUI(L->Aux.nNM+AI.nNM,F); 
#ifdef TEST_Write
	printf("\nd=%d  nV=%d nVmax=%d NUCmax=%d #lists=%d  ",
	    L->d,fi,nVmax,NUCmax,li); fflush(stdout);
	printf("%dnf %dsm %dnm %lldb   sl: %d %d %d %lld\n",L->Aux.nNF+AI.nNF,
	    L->Aux.nSM+AI.nSM,L->Aux.nNM+AI.nNM,L->Aux.NB+AI.NB,L->SLN,
	    L->slSM,L->slNM,L->NewNB-L->RemNB-AI.NB-2*AI.nNF); fflush(stdout);
#endif
     assert(L->NP-L->nSLP==2*AI.nNF-AI.nNM-AI.nSM
			+2*L->Aux.nNF-L->Aux.nNM-L->Aux.nSM); 
	fputUI(tnb=L->Aux.NB+AI.NB,F); 
	fputUI(L->SLN,F); fputUI(L->slSM,F); fputUI(L->slNM,F);
     assert(L->nSLP==2*L->SLN-L->slNM-L->slSM);
	fputUI(L->NewNB-L->RemNB-AI.NB-2*AI.nNF,F);
     for(v=L->d+1;v<=nVmax;v++) if(AI.nNUC[v]+L->Aux.nNUC[v])
     {  i=0; for(nu=1;nu<=NUCmax;nu++) if(NFnum[v][nu]) i++;
#ifdef TEST_Write
	printf("v%d:n%d ",v,i);
#endif
	fputc(v,F); fputc(i,F);				       /* v #nuc(v) */
	for(nu=1;nu<=NUCmax;nu++) if(NFnum[v][nu])
	{   fputc(nu,F); fputUI(NFnum[v][nu],F);		/* nuc #NF */
	}
     }
     j=0; for(v=L->d+1;v<=nVmax;v++)	if(AI.nNUC[v] || L->Aux.nNUC[v]) 
     for(nu=1;nu<=NUCmax;nu++)	if((tc=NFnum[v][nu]))
     {	unsigned int k; pen+=AI.NFnum[v][nu];	/* n -> i; lipos -> j */
	for(i=0; i < (L->Aux.NFnum[v][nu]); i++)
     	{   if(j<pen) while(L->PE[j].n == i) 
	    {	C=&L->NewNF[L->PE[j++].c+2];for(k=0;k<nu;k++)fputc(C[k],F);
		tc--; tnb-=nu; if(j==pen) break;
            }
	    C=&L->Aux.NF[v][nu][nu*i]; for(k=0;k<nu;k++) fputc(C[k],F); 
	    tnb-=nu; tc--;
     	}
	while(j<pen)
	{   C=&L->NewNF[L->PE[j++].c+2]; for(k=0;k<nu;k++)fputc(C[k],F); 
	    tnb-=nu; tc--;
     	}        				assert(tc==0);
     }						assert(tnb==0);
     for(n=0;n<L->SLN;n++)					/* =S->nSLP */
     {	C=&L->NewNF[L->SLp[n]]; nu=C[1]+2; for(j=0;j<nu;j++) fputc(C[j],F); 
	tc += C[1] + 2;
     }				assert(tc==L->NewNB-L->RemNB-AI.NB-2*AI.nNF);
#ifdef	More_File_IO_Data
     printf("\n    nv<=%d  nuc<=%d  files=%d  lists=%d  nNF=%d  NB=%lld ..",
	nVmax,NUCmax,fi,li,AI.nNF,AI.NB);     
#endif
     if(NULL!=L->Aux.NFli) {free(L->Aux.NFli); L->Aux.NFli=NULL;}
}

void Write_Aux_File(NF_List *S)
{    time_t Tstart=time(NULL); FILE *F; 
     int NCalloc=strlen(S->oname)+strlen(SAVE_FILE_EXT)+ (USE_TMP_DIR ? 6 : 1);
     char *auxfn = (char *) malloc(NCalloc);
#ifdef TEMP_FILE_EXT
     char *tmpfn = (char *) malloc(1+strlen(S->oname)+4);
     strcpy(tmpfn,S->oname); strcat(tmpfn,".tmp"); 
#else
     char *tmpfn = auxfn;
#endif
     if(USE_TMP_DIR) {strcpy(auxfn,"/tmp/"); strcat(auxfn,S->oname);} 
     else strcpy(auxfn,S->oname);
     strcat(auxfn,SAVE_FILE_EXT);
     F=fopen(tmpfn,"wb"); printf("Writing %s: ",auxfn); fflush(stdout); 
     if(F==NULL){puts("Cannot open!");exit(0);}
     Write_Bin_File(F,S); if(ferror(F)) {puts("File ERROR!!");exit(0);}
     fclose(F); printf(" done: %ds\n",(int) difftime(time(NULL),Tstart));
     fflush(stdout); free(auxfn); 
#ifdef  TEMP_FILE_EXT
     rename(tmpfn,auxfn); free(tmpfn);
#endif 
}

void Write_List_2_File(char *fn,NF_List *S)
{    time_t Tstart=time(NULL); FILE *F=fopen(fn,"wb");
     printf("Writing %s: ",fn); fflush(stdout);
     if(F==NULL){puts("Cannot open!");exit(0);}
     Write_Bin_File(F,S); if(ferror(F)) {puts("File ERROR!!");exit(0);}
     fclose(F); printf(" done: %ds\n",(int) difftime(time(NULL),Tstart));
     S->SAVE=time(NULL); 				fflush(stdout);
}
void ReAlloc_SortList(NF_List *L)
{    Write_Aux_File(L); Read_Aux_File(L); /* Test_SLnbMS(L); TestMSbits(L); */
     L->SAVE=time(NULL);
}
void CheckLastSaveTime(NF_List *L, int maxsec)
{    if((int)difftime(time(NULL),L->SAVE) > maxsec) 
     {	Print_Statistics(L);  ReAlloc_SortList(L);
     }
}
int  PEntComp(int *FIpos,int *nv,int *nuc,unsigned char *uc, /* pe/pos - uc */
	PEnt *pe,NF_List *S)
{    unsigned char *C=&S->NewNF[pe->c]; 
     int i = (int) C[0] - *nv; if(i) return i;	
     i = (int) C[1] - *nuc; if(i) return i;
#ifdef ACCEL_PEntComp
     i = pe->n - *FIpos; if(i) return i;
#endif
     return RIGHTminusLEFT(uc,&(C[2]),nuc);
}
int  SearchPEntList(int *nv,int *nuc,unsigned char *uc,	    /* 0 :: new     */
                    int *PEpos, int *FIpos, NF_List *S)	   /* -1 exists(SL) */
{    int i, n0, n1;					   /* +1 ex. honest */
     *PEpos=0; if(!S->PEN) return 0;
     if((i = PEntComp(FIpos,nv,nuc,uc,&S->PE[n0=0],S)))	  /* i = PE[0] - uc */
     {  if(i>0) { return 0;}}
     else return 1;
     if((i = PEntComp(FIpos,nv,nuc,uc,&S->PE[n1=S->PEN-1],S)))
     {  if(i<0) { *PEpos=S->PEN; return 0;}}
     else { *PEpos=n1; return 1; }
     while(n1>n0+1) 
     {  *PEpos=(n0+n1)/2; i=PEntComp(FIpos,nv,nuc,uc,&S->PE[*PEpos],S);
        if(i) if(i>0) n1=*PEpos; else n0=*PEpos; 
        else return 1;
     }
     *PEpos=n1; return 0; /* exists: return 1; else PEpos=where it would go */
}
void Insert_PPent_into_Pent(NF_List *S)
{    int Mpos=S->PEN, Spos=S->PPEN, pos=Mpos+Spos;
     if(S->PPEN<=0)
     {  puts("This should not happen in Insert_PPent_into_Pent");exit(0);
     }						assert(pos<=SAVE_INC+SL_Nmax);
     while(Spos--)
     {  if(S->PPE[Spos].n == Mpos) S->PE[--pos] = S->PPE[Spos].pe;
        else break;
     }
     if(++Spos) while(Spos--)
     {  while(S->PPE[Spos].n < Mpos) S->PE[--pos] = S->PE[--Mpos];
        S->PE[--pos] = S->PPE[Spos].pe;
     }
     S->PEN+=S->PPEN; S->PPEN=0;	/* Test_NF_List(S); */
}
int  PPEntComp(int *FIpos,int *PEpos,int *nv,int *nuc,unsigned char *uc,
	PPEnt *ppe,NF_List *S)				    /* ppe/pos - uc */
{    unsigned char *C=&S->NewNF[ppe->pe.c]; 
     int i = (int) C[0] - *nv; if(i) return i;	
     i = (int) C[1] - *nuc; if(i) return i;
#ifdef ACCEL_PEntComp
     i = ppe->n - *PEpos; if(i) return i; if((i=ppe->pe.n -*FIpos)) return i;
#endif
     return RIGHTminusLEFT(uc,&(C[2]),nuc);
}
void InsertPNFintoPPEntList(int *spos, int *mpos, int *lpos,
                            int *nv, int *nuc, unsigned char *uc, NF_List *S)
{    int l=S->PPEN++; unsigned char *C=&(S->NewNF[S->NewNB]); 
     PPEnt *ppe=&S->PPE[*spos];		static int AddListLength; 
     if(0==AddListLength) AddListLength=ADD_LIST_LENGTH;
#ifdef	USE_UNIT_ENCODE
	if((*uc % 8) > 3) S->NC ++;
#endif
     if(*nuc+2+S->NewNB > S->ANB) 
     {  printf("increase CperR_MAX or write/read %s\n",S->iname); exit(0);
     }
     while(*spos < (l--)) S->PPE[l+1]=S->PPE[l];
     ppe->n=*mpos; ppe->pe.n=*lpos; ppe->pe.c=S->NewNB;
     C[0]=*nv; C[1]=*nuc; C=&C[2]; for(l=0;l<*nuc;l++) C[l]=uc[l];
     S->NewNB += *nuc+2; S->NP++; if(*C % 4) S->peNM++; else S->peSM++;
     if(S->PPEN==AddListLength) Insert_PPent_into_Pent(S);
}
int  SearchPPEntList(int *nv,int *nuc,unsigned char *uc,    	/* new :: 0 */
		int *Ppos, int *Mpos,int *Lpos, NF_List *S)  /* exists :: 1 */
{    int i, n0, n1;	   
     *Ppos=0; if(!S->PPEN) return 0;
     if((i = PPEntComp(Lpos,Mpos,nv,nuc,uc,&S->PPE[n0=0],S)))
     {  if(i>0) { *Ppos=0; return 0;}}
     else return 1;
     if((i = PPEntComp(Lpos,Mpos,nv,nuc,uc,&S->PPE[n1=S->PPEN-1],S)))
     {  if(i<0) { *Ppos=S->PPEN; return 0;}}
     else  {*Ppos=n1; return 1;}
     while(n1>n0+1) 
     {  *Ppos=(n0+n1)/2; i=PPEntComp(Lpos,Mpos,nv,nuc,uc,&S->PPE[*Ppos],S);
        if(i) if(i>0) n1=*Ppos; else n0=*Ppos; 
        else return 1;
     }
     *Ppos=n1; return 0;   /* exists: return 1; else Ppos=where it would go */
}
int  SearchFIList(int *v,int *nu,unsigned char *uc, int *pos,FInfoList *FI)
{    int i, n0, n1; *pos=0; if(!FI->nNF) return 0; 
     if(!FI->NFnum[*v][*nu]) return 0;
     if((i = RIGHTminusLEFT(uc,&FI->NF[*v][*nu][n0=0],nu)))
     {  if(i>0) return 0;}
     else return 1;
     if((i = RIGHTminusLEFT(uc,
			&FI->NF[*v][*nu][(*nu)*(n1=FI->NFnum[*v][*nu]-1)],nu)))
     {  if(i<0) { *pos=n1+1; return 0;}}
     else  {*pos=n1; return 1;}
     while(n1>n0+1) 
     {  *pos=(n0+n1)/2; i=RIGHTminusLEFT(uc,
			&FI->NF[*v][*nu][(*nu)*(*pos)],nu);
        if(i) if(i>0) n1=*pos; else n0=*pos; 
        else return 1;
     }
     *pos=n1; return 0;   /* exists: return 1; else pos=where it would go */
}
int  SL_Comp(int *nv,int *nuc,unsigned char *uc,int *lipos,NF_List *S)
{    unsigned char *C=&S->NewNF[*lipos]; 		     /* lipos - uc */
     int i = (int) C[0] - *nv; if(i) return i;	
     i = (int) C[1] - *nuc; if(i) return i;
     return RIGHTminusLEFT(uc,&(C[2]),nuc);
}
int  SearchSL_List(int *nv,int *nuc,unsigned char *uc, int *SLnum,NF_List *S)
{    int i, n0, n1;		   	/* 1 = exists (mirror or straight) */
     *SLnum=0; if(!S->SLN) return 0;
     if((i = SL_Comp(nv,nuc,uc,&S->SLp[n0=0],S))) 
     {  if(i>0) { return 0;}}
     else return 1;
     if((i = SL_Comp(nv,nuc,uc,&S->SLp[n1=S->SLN-1],S)))
     {  if(i<0) { *SLnum=S->SLN; return 0;}}
     else { *SLnum=n1; return 1;}
     while(n1>n0+1) 
     {  *SLnum=(n0+n1)/2; i=SL_Comp(nv,nuc,uc,&S->SLp[*SLnum],S);
        if(i) if(i>0) n1=*SLnum; else n0=*SLnum; 
        else return 1;
     }
     *SLnum=n1; return 0; /* exists: return 1; else SLnum=where it would go */
}
void SL_List_Insert(int *nv,int *nuc,unsigned char *uc,int *SLnum,NF_List *S)
{    int l=S->SLN++; unsigned char *C=&(S->NewNF[S->NewNB]); 
     if(*nuc+2 + S->NewNB > S->ANB) 
     {  printf("increase CperR_MAX or write/read %s\n",S->iname); exit(0);
     }
     while(*SLnum < (l--)) S->SLp[l+1]=S->SLp[l]; 
     S->SLp[*SLnum]=S->NewNB;
     C[0]=*nv; C[1]=*nuc; C=&C[2]; for(l=0;l<*nuc;l++) C[l]=uc[l];
     S->NewNB+=*nuc+2; S->NP++;S->nSLP++; if(*C % 4)S->slNM++; else S->slSM++;
     if(S->SLN==SL_Nmax) {puts("Increase SL_Nmax");exit(0);}
}
void SL_List_Remove(int *NV,int *nUC, /* unsigned char *UC, */
		    int *SLnum, NF_List *S)
{    int l=*SLnum; unsigned char *C=&S->NewNF[S->SLp[l]];
     S->RemNB += (2+C[1]); if(C[2] % 4) S->slNM--; else S->slSM--; 
     assert((*NV==C[0]) && (*nUC==C[1]));
     /* assert(0==RIGHTminusLEFT(UC,&C[2],nUC)); */
     while ((++l) < S->SLN) S->SLp[l-1]=S->SLp[l]; 
     S->SLN--; S->NP--; S->nSLP--;
}

#define	MirTest(A,B)	((((A)+(B)) % 4) != 3)	   /* 1=same, 0=only mirror */
/*	ucNF_Sort_Add = 1 = continue = V is honest and new ;
 *			0 = SL or Ex(honest)
 *      Search...List = 1 : 0 : -1  iff  honest : new : sublattice
 */
int  ucNF_Sort_Add(int *NV, int *nUC, unsigned char *UC, NF_List *S)
{    int InPos=0, FIpos=0, PEpos=0, NEWpos=0, NewH=0, NewSL=0;
#ifdef  INCREMENTAL_WRITE
     if(SearchFIList(NV,nUC,UC,&InPos, &S->In)) 
     if(MirTest(*UC,S->In.NF[*NV][*nUC][(*nUC)*InPos])) 
     {	if(!S->SL) S->hc++; goto Ret0;}
#endif
     if(*S->dbname) if(Is_in_DB(NV,nUC,UC,S)) {if(!S->SL) S->hc++; goto Ret0;}

     if(SearchFIList(NV,nUC,UC,&FIpos,&S->Aux))
     {	unsigned char *c=&(S->Aux.NF[*NV][*nUC][(*nUC)*FIpos]);
	if(MirTest(*UC,*c)) goto Ret0;   /* if(S->con==S->AP) S->hc++; */
	else if(!S->SL) 
	{   (*c) += ((*UC) % 4); S->Aux.nNM--; S->NP++; NewH=1; 
	}				/* Mir(honest) on Aux => not on New */
     }
     else if(SearchPEntList (NV,nUC,UC,&PEpos,&FIpos,S))
     {	unsigned char *c=&(S->NewNF[S->PE[PEpos].c+2]);
	if(MirTest(*UC,*c)) goto Ret0;
	else if(!S->SL)
	{   (*c) += ((*UC) % 4); S->peNM--; S->NP++; NewH=1; 
	} 			         /* Mir(honest) on PE => not on PPE */
     }
     else if(SearchPPEntList(NV,nUC,UC,&NEWpos,&PEpos,&FIpos,S))
     {	unsigned char *c=&(S->NewNF[S->PPE[NEWpos].pe.c+2]);
	if(MirTest(*UC,*c)) goto Ret0;
	else if(!S->SL)
	{   (*c) += ((*UC) % 4); S->peNM--; S->NP++; NewH=1; 
	}
     }
     else if(!S->SL)
     {	InsertPNFintoPPEntList(&NEWpos,&PEpos,&FIpos,NV,nUC,UC,S); NewH=1;
     }
						assert(NewH + S->SL == 1);
     if(SearchSL_List(NV,nUC,UC,&NEWpos,S))
     {	unsigned char *c=&(S->NewNF[S->SLp[NEWpos]+2]);	
	if(NewH)
	{   if(MirTest(*UC,*c))
	    {	if((*c % 4) != 3) SL_List_Remove(NV,nUC, /*UC,*/ &NEWpos,S);
		else { (*c) -= (*UC % 4); S->slNM++; S->NP--; S->nSLP--;}
		NewSL=-1;
	    }
	}
	else						 /* SubLattice case */
	{   if(MirTest(*UC,*c)) goto Ret0;
	    else { (*c) += ((*UC) % 4); S->NP++; S->slNM--;S->nSLP++;NewSL=1;}
	}
     }
     else if(S->SL) {assert(NewH==0);
			NewSL=1; SL_List_Insert(NV,nUC,UC,&NEWpos,S);}

     if(NewH+NewSL) if(S->NP % WATCHREF == 0) 
     {  Print_Statistics(S); 
	if(S->NP - S->savedNP > SAVE_INC - WATCHREF) 
	/* printf("NP=%d  PEN=%d  peNM=%d  peSM=%d\n", S->NP,S->PEN,S->peNM,
	   S->peSM);fflush(stdout); 
	   if(2 *(S->PEN+S->PPEN) - S->peNM - S->peSM > SAVE_INC - WATCHREF)*/
     	{  ReAlloc_SortList(S); S->con=0; }
     }	if(S->NP%1000 == 0) CheckLastSaveTime(S,(19*GOOD_SAVE_TIME)/20);
     CheckLastSaveTime(S,GOOD_SAVE_TIME);			return NewH;
     Ret0:  CheckLastSaveTime(S,FORCE_SAVE_TIME); return 0;
}
/*      ==============================================================      */


/*      ==============================================================      */
void Print_Statistics(NF_List *_L)
{    clock_t CLOCK=clock(); time_t DATE=time(NULL);
     int CPUsec=(CLOCK-_L->CLOCK)/CLOCKS_PER_SEC;  /* CLOCKS_PER_SEC::10^6 */
     int REALsec= (int) difftime(DATE,_L->TIME), NFsec; /* int IPperNF; */
     char bni[2]; int BNI=(_L->nNF>2000000) ? 1000000 : 1000;
     strcpy(bni,(BNI==1000) ? "k" : "M");
     printf("%dkR-%d %dMB %d%sIP %d%sNF-%dk ", (int)_L->NP/1000,
	(int)_L->nSLP,(int) ((_L->NewNB-_L->RemNB+_L->Aux.NB)/1000000), (int)
	(_L->nIP/BNI),bni,(int)_L->nNF/BNI,bni,(int)_L->nSLNF/1000);  
     printf("%d%s%d v%dr%d f%dr%d %db%d",
	_L->Nmin, _L->hc ? "-" : "_", _L->Nmax, _L->VN,_L->V, _L->FN,_L->F, 
	_L->Xdif, _L->Xnuc); 
#ifdef	INCREMENTAL_TIME
     /* if(_L->NF_Time>0)IPperNF=_L->IP_Time/_L->NF_Time; else IPperNF=-1; */
     NFsec=_L->NF_Time/CLOCKS_PER_SEC; 
     CPUsec=(_L->IP_Time+_L->NF_Time)/CLOCKS_PER_SEC;
     if(CPUsec<1000)printf(" %ds %du %dn",REALsec,CPUsec,NFsec);
     else if((CPUsec)<60000)
	printf(" %dm %du %dn",REALsec/60,CPUsec/60,NFsec/60); 
     else printf(" %dh %du %dn",REALsec/3600,CPUsec/3600,NFsec/3600); 
#else
     CPUsec=(CLOCK-_L->CLOCK)/CLOCKS_PER_SEC;
     if(CPUsec<1000)printf(" %ds %du",REALsec,CPUsec); 
     else if((CPUsec)<60000) printf(" %dm %du",REALsec/60,CPUsec/60); 
     else printf(" %dh %du",REALsec/3600,CPUsec/3600); 
#endif
     _L->Nmin=1000; _L->Nmax=0; puts(""); fflush(stdout);
}
void Print_Expect(FInfoList *L)
{    int s=L->nSM; Along m=L->nNF-L->nNM-L->nSM, p=L->nNF+m; double x=p;
     /* if(d<2) {if(m)printf("  p^2/2m=%d",(p*p)/(2*m));
	      if(s)printf("  pp/ss=%d",(p*p)/(s*s));} */
     /* if(m) {x*=(x+s)/(2*m);printf(" pP/2m=%g",x);} */
		if(m) {x=p+s; x*=x/(2*m+s);printf(" pp/2m=%g",x);}
		if(s) {x=p; x/=s; x*=x; printf(" pp/ss=%g",x);}
}
/*      ==============================================================      */


/*      ==============================================================      */
/*		SearchLongList:= { new(below *pos)::0, found::1 }	    */
/*		SearchPEntLists:={ new::0, found::1 }			    */
/*   Add_NF_to_List()			assumes that _V and _F are assigned
 *
 *	if(SL)	{ if(new) "add to dishonest"; 	  	return 0; }	::stop
 *	if(exists on honest) 				return 0;	::stop
 *	if(exists on dishonest)	{ "move to honest"; 	return 1; }	::cont
 *	"add to honest";				return 1; 	::cont
 *
 *   Make_Poly_NF (Make_VPM_NF; Make_PNF_from_VNF;);	
 *   PNF_Sort_Add (NF[][],*d,*np,*nv,*nf,*S) := { List::+1, new::0, PEnt::-1 }
 *   SearchLongList:= { new(below *pos)::0, found::1 }
 *   SearchPEntLists:={ new::0, found::1 }
 */
int  Add_NF_to_List(PolyPointList *_P, VertexNumList *_V, EqList *_E,
                    NF_List *_L)
{    unsigned char UC[NB_MAX]; int nUC, NV;
     int NewNF; 
#ifdef INCREMENTAL_TIME
     long long cpuT=clock(), incT=cpuT - _L->CLOCK; 
     if(incT>0) _L->IP_Time+=incT;
     _L->CLOCK=cpuT;
#endif
     if((_E->ne)>_L->F) _L->F=_E->ne;		/* check vertex numbers */
     if((_V->nv)>_L->V) _L->V=_V->nv;
     if((_E->ne > VERT_Nmax)||(_V->nv > VERT_Nmax))
     {	printf("Increase VERT_Nmax:  f=%d v=%d\n",_E->ne,_V->nv);exit(0);
     }
     if(_P->np>_L->Nmax) _L->Nmax=_P->np;
     if(_P->np<_L->Nmin) _L->Nmin=_P->np;
     _L->nNF++; if(_L->SL)_L->nSLNF++;  

     VF_2_ucNF(_P, _V, _E,  &NV, &nUC, UC); 

     NewNF=ucNF_Sort_Add(&NV, &nUC, UC, _L);  /* 1::new::cont. */

#ifdef INCREMENTAL_TIME
     cpuT=clock(), incT=cpuT - _L->CLOCK; 
     if(incT>0) _L->NF_Time+=incT;
     _L->CLOCK=cpuT;
#endif
     return NewNF;
}

/* 				     #R sl hit #IP #NF TIME( r > u > AddNF) */
void Print_Weight_Info(CWS *W,NF_List *_L)
{    int i,j;
#if	POLY_Dmax>3
     Print_Statistics(_L);
#endif
     for(i=0;i<W->nw;i++)
     {   fprintf(outFILE,"%ld ",W->d[i]);                 
         for(j=0;j<W->N;j++) fprintf(outFILE,"%ld ", W->W[i][j]);
     }	
     for(i=0;i<W->nz;i++) 
     {	fprintf(outFILE,"/Z%d: ",W->m[i]);
	for(j=0;j<W->N;j++) fprintf(outFILE,"%d ",W->z[i][j]);
     }	
     fprintf(outFILE,"R=%lld +%dsl hit=%d IP=%d NF=%d (%d)\n",
	_L->NP-_L->nSLP, _L->nSLP,_L->hc,_L->nIP,_L->nNF,_L->nSLNF); 
     fflush(outFILE);
}

/*   ===============	    compression package		=================== */

#if	(INT_MAX != 2147483647)		/* UINT_MAX == 4294967295 */	
#error	use other date types	      /* ULLONG_MAX == 18446744073709551615 */
#endif
#define	UCM		256		/* Unsigned Char Modulo  (= max+1)  */
#define USM		65536		/* Unsigned Short Modulo (= max+1)  */
#define Nint_XLong	(NB_MAX + 1) / 2
#define NX(d,v)		((d)*(v)-((d)*(d-1))/2)

#ifdef	USE_UNIT_ENCODE
#define UNIT_NX(d,v)		((d)*((v)-(d))+1)
#define UNIT_OFF		(8)
#else
#define UNIT_OFF		(4)
#endif

#ifndef	LL_BASE		
#define	LL_BASE		32767	  /* limit for 64-bit BaseGetInt: cf. NF of */
#endif			    /* 3198174 49 1723 74375 456882 1066058 1599087 */

typedef struct {int n; unsigned short x[Nint_XLong];}		UXLong;

void XPrint(UXLong *X)
{    int i; printf("X.n=%d:",X->n);for(i=0;i<X->n;i++)printf(" %d",X->x[i]);
}
void VPrint(int *d,int *v,Long V[POLY_Dmax][VERT_Nmax])
{    int i,j; puts("");for(i=0;i<*d;i++)
     {	for(j=0;j<*v;j++)printf("%4ld ",V[i][j]);puts("");}
}
void LLBasePutInt(int *base, int *x, UXLong *X)   /* X = X * base + x */
{    int i=0; long long z, ad=*x; 
     while(i<X->n)
     {	z=ad+(*base)*((long long)X->x[i]); X->x[i++]=z% USM; ad=z / USM; 
     }  
     while(ad) { assert(X->n < Nint_XLong); X->x[X->n++]=ad% USM; ad/=USM; }
}
void LLBaseGetInt(int *base, int *x, UXLong *X)  /* x = X mod b, X /= b */
{    int i; long long a, zq, zr;
     if(!X->n) {*x=0; return;}		a=X->x[(i=X->n-1)];
     while(i) {	zq=a/(*base); zr=a%(*base); X->x[i]=zq; a=zr*USM+X->x[--i]; }
     zq=a/(*base); zr=a%(*base); X->x[i]=zq; *x=zr;
     for(i=X->n-1; 0<=i; i--) if(X->x[i]) return; else X->n--; 
}
void BasePutInt(int *base, int *x, UXLong *X)   /* X = X * base + x */
{    int i=0, z, ad=*x; if(*base > LL_BASE) {LLBasePutInt(base,x,X); return;}
     while(i<X->n)
     {	z=ad+(*base)*X->x[i]; X->x[i++]=z% USM; ad=z / USM; 
     }  
     if(ad) { assert(X->n <= (Nint_XLong-1)); X->x[X->n++]=ad; }
}
void BaseGetInt(int *base, int *x, UXLong *X)  /* x = X mod b, X /= b */
{    int i, a; div_t z;	if(*base > LL_BASE) {LLBaseGetInt(base,x,X); return;}
     if(!X->n) {*x=0; return;}		a=X->x[(i=X->n-1)];
     while(i) {	z=div(a,*base); X->x[i]=z.quot; a=z.rem*USM+X->x[--i]; }
     z=div(a,*base); X->x[i]=z.quot; *x=z.rem; 
     for(i=X->n-1; 0<=i; i--) if(X->x[i]) return; else X->n--; 
}
int  BminOff(Long V[POLY_Dmax][VERT_Nmax], int *d, int *v, int *off, int *bmin)
{    int i, j, moff=0, vmax=0;
     for(i=0;i<*d;i++) for(j=i;j<*v;j++)
     {	Long x=V[i][j]; if(moff>x) moff=x;
       if(vmax<x) vmax=x;
     }
     *off=-moff; *bmin=*off+vmax+1;
#ifdef	USE_UNIT_ENCODE
     for(i=0;i<*d;i++) for(j=i;j<*d;j++) if(V[i][j]!=(i==j)) return 0;
     return 1;
#else
     return 0;
#endif
}
void AuxBase2nUC(int *base, int *nx, int *nuc)      /* b^nx-2=x0+x1*USM+... */
{    int i, m=*base-1; UXLong X; X.n=0; /* m=X->n-1 =>	  ... + xm*USM^xm   */
     for(i=0;i<*nx-1;i++) BasePutInt(base,&m,&X);
     m--; BasePutInt(base,&m,&X);
     i=UNIT_OFF; m=i-1;	BasePutInt(&i,&m,&X); *nuc=2*X.n-(X.x[X.n-1]<256);
}
void AuxNextGoodBase(int *v,int *nx,Base_List *BL)  /* nx= v*d -d*(d-1)/2 */
{    int *bo=&(BL->base[*v][BL->v[*v]-1]), nuct, 
	 *bn=&(BL->base[*v][BL->v[*v]]), *nucn=&(BL->nuc[*v][BL->v[*v]]);
     *bn=*bo+1; AuxBase2nUC(bn,nx,nucn);
     do { (*bn)++; AuxBase2nUC(bn,nx,&nuct); } while (nuct==*nucn);
     (*bn)--; BL->v[*v]++;
}
void Init_BaseList(Base_List **bl, int *d) /* once: alloc and init BaseList */
{    static Base_List *BL=NULL;	int i; 	   /* always: set *bl=BL=&BaseList */ 
     if(BL!=NULL) {*bl=BL; return;}
     *bl = BL = (Base_List *) malloc(sizeof(Base_List)); assert(BL!=NULL);
     for(i=*d; i<=VERT_Nmax; i++)	
     {	int nx= NX(*d,i), bn=3, nuco, nucn; 
        BL->v[i]=1; AuxBase2nUC(&bn,&nx,&nuco); nucn=BL->nuc[i][0]=nuco;
        while(nucn==nuco) { bn++; AuxBase2nUC(&bn,&nx,&nucn); }   
	BL->base[i][0]=bn-1;
     }
}
void Bmin2BaseUCn(int *d, int *v, int *bmin,  int *base, int *nuc)
{    int i, *n, nx = NX(*d,*v);
     Base_List *BL; Init_BaseList(&BL,d); n=&(BL->v[*v]); i=*n-1;

     if(*bmin<=BL->base[*v][i])
     { while(i > 0 && *bmin <= (BL->base[*v][i-1])) --i;
     }
     else
     {	do AuxNextGoodBase(v,&nx,BL); while(*bmin > (BL->base[*v][++i]));
     }
     *base=BL->base[*v][i]; *nuc=BL->nuc[*v][i];
}
void NUCtoBase(int *d, int *v, int *nuc, int *Base)
{    int i, *n, nx = NX(*d,*v);
     Base_List *BL; Init_BaseList(&BL,d); n=&(BL->v[*v]); i=*n-1;
     if(*nuc<=BL->nuc[*v][i])
     { while(i > 0 && *nuc <= (BL->nuc[*v][i-1])) --i;
     }
     else
     {	do AuxNextGoodBase(v,&nx,BL); while(*nuc > (BL->nuc[*v][++i]));
     }
     assert(*nuc==BL->nuc[*v][i]); *Base=BL->base[*v][i]; 
}

#ifdef	USE_UNIT_ENCODE
void UNIT_Init_BaseList(Base_List **bl, int *d) /* once: ainit BaseList */
{    static Base_List *BL=NULL;	int i; 	   /* always: set *bl=BL=&BaseList */ 
     if(BL!=NULL) {*bl=BL; return;}
     *bl = BL = (Base_List *) malloc(sizeof(Base_List)); assert(BL!=NULL);
     for(i=*d; i<=VERT_Nmax; i++)	
     {	int nx= UNIT_NX(*d,i), bn=3, nuco, nucn; 
        BL->v[i]=1; AuxBase2nUC(&bn,&nx,&nuco); nucn=BL->nuc[i][0]=nuco;
        while(nucn==nuco) { bn++; AuxBase2nUC(&bn,&nx,&nucn); }   
	BL->base[i][0]=bn-1;
     }
}
void UNIT_Bmin2BaseUCn(int *d, int *v, int *bmin,  int *base, int *nuc)
{    int i, *n, nx = UNIT_NX(*d,*v);
     Base_List *BL; UNIT_Init_BaseList(&BL,d); n=&(BL->v[*v]); i=*n-1;
     if(*bmin<=BL->base[*v][i])
     { while(i > 0 && *bmin <= (BL->base[*v][i-1])) --i;
     }
     else
     {	do AuxNextGoodBase(v,&nx,BL); while(*bmin > (BL->base[*v][++i]));
     }
     *base=BL->base[*v][i]; *nuc=BL->nuc[*v][i];
}
void UNIT_NUCtoBase(int *d, int *v, int *nuc, int *Base)
{    int i, *n, nx = UNIT_NX(*d,*v);
     Base_List *BL; UNIT_Init_BaseList(&BL,d); n=&(BL->v[*v]); i=*n-1;
     if(*nuc<=BL->nuc[*v][i])
     { while(i > 0 && *nuc <= (BL->nuc[*v][i-1])) --i;
     }
     else
     {	do AuxNextGoodBase(v,&nx,BL); while(*nuc > (BL->nuc[*v][++i]));
     }
     assert(*nuc==BL->nuc[*v][i]); *Base=BL->base[*v][i]; 
}
#endif
							   /* UNIT :: ms+=4 */
void AuxVnf2ucNF(Long nf[POLY_Dmax][VERT_Nmax], int *d, int *v, 
		int *off, int *base, int *nuc, int *ms,     unsigned char *UC) 
{    int i, j, one=(*ms > 3);
     UXLong X; X.n=0;
     for(i=0;i<*d;i++) for(j = (one) ? *d : i ;j<*v;j++)
     {	int x=nf[i][j]+(*off); BasePutInt(base,&x,&X);
     }	
     if(one) BasePutInt(base,off,&X);
     i=UNIT_OFF; BasePutInt(&i,ms,&X);
     for(i=X.n;(2*i)<(*nuc);i++) X.x[i]=0;
     for(i=0;i<*nuc;i++)if(i%2) UC[i]=X.x[i/2] /UCM; else UC[i]=X.x[i/2] %UCM;
}
void UCnf2vNF(int *d, int *v, int *nuc, unsigned char *uc,  	      /* IN */
	      Long NF[POLY_Dmax][VERT_Nmax], int *MS)		     /* OUT */
{    int i, j, off=0, base, one;
     UXLong X; X.n=*nuc/2; if(*nuc % 2) X.x[X.n++]=uc[*nuc-1];
     for(i=0;i<*nuc/2;i++) X.x[i]=uc[2*i]+UCM*uc[2*i+1];
     i=UNIT_OFF; 
     BaseGetInt(&i,MS,&X); one = (*MS > 3);
#ifdef	USE_UNIT_ENCODE
     if(one) {UNIT_NUCtoBase(d,v,nuc,&base); BaseGetInt(&base,&off,&X);} else
#endif
     NUCtoBase(d,v,nuc,&base);
     i=*d; while(i--)
     {	int m=(one) ? *d : i;
	j=*v;while(m<(j--)) {int Iau;BaseGetInt(&base,&Iau,&X);NF[i][j]=Iau;}
	j=i; while(j--) NF[i][j]=0;
     }	
#ifdef	USE_UNIT_ENCODE
     if(one) {	for(i=0;i<*d;i++) for(j=*d;j<*v;j++) NF[i][j]-=off;
		for(i=0;i<*d;i++) for(j=i;j<*d;j++) NF[i][j]=(i==j);
     }	else
#endif
     { off=NF[0][0]-1; for(i=0;i<*d;i++) for(j=i;j<*v;j++) NF[i][j]-=off;}
     for(i=0;i<*d;i++) for(j=0;j<i;j++) NF[i][j]=0;
}

int  RIGHTminusLEFT(unsigned char *ucL, unsigned char *ucR, int *nuc)
{    int i=*nuc; 
     while(--i) if(ucL[i]!=ucR[i]) return (ucL[i] < ucR[i]) ? 1 : -1;
     if(((*ucL)/4)!=((*ucR)/4)) return (((*ucL)/4) < ((*ucR)/4)) ? 1 : -1;
     return 0;
}
void VF_2_ucNF(PolyPointList *P, VertexNumList *V, EqList *E,	      /* IN */
     		 int *NV, int *nUC, unsigned char *UC)		     /* OUT */
{    Long V_NF[POLY_Dmax][VERT_Nmax], F_NF[POLY_Dmax][VERT_Nmax];
     Long VM[POLY_Dmax][VERT_Nmax],  VPM[VERT_Nmax][VERT_Nmax]; 
     unsigned char auxUC[POLY_Dmax*VERT_Nmax];
     int MS=0, vb,vo,vnuc,vbmin, fb,fo,fnuc,fbmin, vone=0,fone=0,MSone;
#ifndef	SORT_NUC_FIRST
     if(V->nv<E->ne) MS=1;
     if(V->nv>E->ne) MS=2;
#endif					/* else: priority to sorting by nuc */
     assert(Init_rVM_VPM(P,V,E,&P->n,&V->nv,&E->ne,VM,VPM));/* make ref VPM */
     if(MS <2)
     {	Eval_Poly_NF(&P->n,&V->nv,&E->ne,VM,VPM,V_NF,0);	    /* V_NF */
        vone=BminOff(V_NF,&P->n,&V->nv,&vo,&vbmin);
#ifdef	USE_UNIT_ENCODE
	if(vone) UNIT_Bmin2BaseUCn(&P->n,&V->nv,&vbmin,&vb,&vnuc); else
#endif
	Bmin2BaseUCn(&P->n,&V->nv,&vbmin,&vb,&vnuc);
     }
     if(MS!=1)
     {	int i, j, mi=(E->ne>V->nv) ? V->nv : E->ne;     	    /* VM=E */
        for(i=0;i<P->n;i++)for(j=0;j<E->ne;j++) VM[i][j]=E->e[j].a[i];
	for(i=0;i<mi;i++)for(j=0;j<i;j++)		   /* transpose VPM */
	{   Long a=VPM[i][j]; VPM[i][j]=VPM[j][i]; VPM[j][i]=a;
	}
	if(E->ne>V->nv) 
	     for(i=mi;i<E->ne;i++)for(j=0;j<V->nv;j++)VPM[j][i]=VPM[i][j];
	else for(i=mi;i<V->nv;i++)for(j=0;j<E->ne;j++)VPM[i][j]=VPM[j][i];
	Eval_Poly_NF(&P->n,&E->ne,&V->nv,VM,VPM,F_NF,0);    /* compute F_NF */
	fone=BminOff(F_NF,&P->n,&E->ne,&fo,&fbmin);
#ifdef	USE_UNIT_ENCODE
	if(fone) UNIT_Bmin2BaseUCn(&P->n,&E->ne,&fbmin,&fb,&fnuc); else
#endif
	Bmin2BaseUCn(&P->n,&E->ne,&fbmin,&fb,&fnuc);
     }
     if(MS==0) { if(vnuc<fnuc)   MS=1; if(fnuc<vnuc)   MS=2; }  /* sort nuc */
     if(MS==0) { if(V->nv<E->ne) MS=1; if(E->ne<V->nv) MS=2; } /* sort vert */
#ifdef	USE_UNIT_ENCODE
     if(MS==0) { if(vone>fone)   MS=1; if(vone<fone)   MS=2; } /* sort UNIT */
#endif

     *nUC = (MS==2) ? fnuc : vnuc;   *NV = (MS==2) ? E->ne : V->nv;  MSone=MS;

     if(AuxNFLptr!=NULL)			    /* base:byte statistics */
     {	/* if(AuxNFLptr->Xmin>-vo) AuxNFLptr->Xmin=-vo;
	   if(AuxNFLptr->Xmax<vbmin-vo-1) AuxNFLptr->Xmax=vbmin-vo-1; */
	if(MS==2) {if(AuxNFLptr->Xdif<fbmin) AuxNFLptr->Xdif=fbmin;}
	else 	  {if(AuxNFLptr->Xdif<vbmin) AuxNFLptr->Xdif=vbmin;}
	if(AuxNFLptr->Xnuc<*nUC) AuxNFLptr->Xnuc=*nUC;
	if(((MS==2)&&(fbmin>BASE_MAX))||((MS<2)&&(vbmin>BASE_MAX))) 
	{   printf("WARNING MS=%d  v=%d vb=%d vnuc=%d  f=%d fb=%d fnuc=%d\n",
		MS,V->nv,vbmin,vnuc,E->ne,fbmin,fnuc);
	    VPrint(&P->n,&V->nv,V_NF);puts("============= V_NF");
	    VPrint(&P->n,&E->ne,F_NF);puts("============= F_NF"); 
	    puts("BASE_MAX exceeded"); exit(0);
	}
     }

#ifdef	USE_UNIT_ENCODE
     if(MS <2) MSone += 4*vone;	else MSone += 4*fone;
#endif
     if(MS <2) AuxVnf2ucNF(V_NF, &P->n, &V->nv, &vo, &vb, &vnuc, &MSone, UC);
     else      AuxVnf2ucNF(F_NF, &P->n, &E->ne, &fo, &fb, &fnuc, &MSone, UC);

     if(MS==0)				/* nuc and nv agree => compare UC[] */
     {	int RmL; 
	AuxVnf2ucNF(F_NF, &P->n, &E->ne, &fo, &fb, &fnuc, &MSone, auxUC);
	RmL=RIGHTminusLEFT(UC,auxUC,nUC);  
	if(RmL>0) {MS=1; (*UC)++;}        	/* o.k. since *UC=4*(...) */
	if(RmL<0) {int i; for(i=0;i<*nUC;i++)UC[i]=auxUC[i]; (*UC) += (MS=2);}
     }
	/* 6+8*(2+17(1+17(0+17(1+17(1+17(1+17(2)))))))  {23  189  20  198}
	*/
#ifdef TEST_UCnf		/* One=K3[3/F]:  2 1 0 1 1 1 2  */
     {	Long tNF[POLY_Dmax][VERT_Nmax]; int tMS, i, j, *d=&P->n, err=0;
	static int pn, Err; 
	pn++;				UCnf2vNF(d, NV, nUC, UC,   tNF, &tMS);
	if(MS<2) for(i=0;i<*d;i++) for(j=i;j<*NV;j++)
				     err+=(tNF[i][j]!=V_NF[i][j]);
	if(err) err*=1000;
	if(MS-1) for(i=0;i<*d;i++) for(j=i;j<*NV;j++)
				     err+=(tNF[i][j]!=F_NF[i][j]);
	if(err) Err=1;
	if(Err) 
	{   int base;
#ifdef	USE_UNIT_ENCODE
	    if(MSone>3) UNIT_NUCtoBase(&P->n,NV,nUC,&base); else
#endif
	    NUCtoBase(&P->n,NV,nUC,&base);
	    printf("\n#%d  MS=%d  100*e(V)-e(F) = %d\n",pn,MS,err);
	    for(i=0;i<*nUC;i++)printf("%d ",UC[i]);
	    printf("= UC[%d]  base=%d  vo=%d  fo=%d\n",*nUC,base,vo,fo);
	    VPrint(&P->n,&V->nv,V_NF);puts("============= V_NF");
	    VPrint(&P->n,&E->ne,F_NF);puts("============= F_NF");
	    VPrint(&P->n,NV,tNF);puts("============= t_NF");	exit(0); 
	}
     }
#endif
     assert((*UC % 4)==MS);
}
void Test_ucNF(int *d, int *v, int *nuc, unsigned char *uc, PolyPointList *_P)
{    Long tNF[POLY_Dmax][VERT_Nmax]; int i, j, tMS, NV, NUC; 
     unsigned char UC[POLY_Dmax*VERT_Nmax];
     VertexNumList V; EqList E; 
     UCnf2vNF(d, v, nuc, uc, tNF, &tMS); _P->n=*d; _P->np=*v; tMS %= 4;
     for(i=0;i<*v;i++)for(j=0;j<*d;j++)_P->x[i][j]=tNF[j][i];
     assert(Ref_Check(_P,&V,&E)); VF_2_ucNF(_P,&V,&E,&NV,&NUC,UC);
     assert(*v==NV); assert(*nuc==NUC);
     assert(RIGHTminusLEFT(uc,UC,nuc)==0);
     assert((tMS==0)==((*uc%4)==0)); 
}

void AuxGet_vn_uc(FILE *F,int *v, int *nu, unsigned char *uc)
{    int i; *v=fgetc(F); *nu=fgetc(F); for(i=0;i<*nu;i++) uc[i]=fgetc(F);
}
void AuxGet_uc(FILE *F,int *nu, unsigned char *uc)
{    int i; for(i=0;i<*nu;i++) uc[i]=fgetc(F);
}

void AuxPut_hNF(FILE *F,int *v,int *nu,unsigned char *Huc,FInfoList *Io,
	int *slNF,int *slSM,int *slNM,int *slNB,unsigned char *ucSL,int *SLp)
{    int i, Hms=(*Huc % 4); static int pos; unsigned char *Suc=NULL;
     for(i=0;i<*nu;i++) fputc(Huc[i],F);
     if(Hms) {if(Hms<3) Io->nNM++;} else Io->nSM++;
     while(pos<*slNF)
     {	int HmSL; Suc=&ucSL[SLp[pos]]; HmSL=*v-*Suc; 
	if(HmSL>0) {pos++; continue;} else if(HmSL<0) return;
	else HmSL=*nu-Suc[1];
	if(HmSL>0) {pos++; continue;} else if(HmSL<0) return; 
	else HmSL=RIGHTminusLEFT(&Suc[2],Huc,nu);
	if(HmSL>0) {pos++; continue;} else if(HmSL<0) return; else break;
     }
     if(pos<*slNF) 
     {	int Sms=Suc[2]%4; pos++; switch(10*Hms+Sms) {
	case 00: case 11: case 22: case 31: case 32: case 33: break;
	case 12: case 21: return;
	case 13: Suc[2]-=1; (*slNM)++; return;		/* remove 1st */
	case 23: Suc[2]-=2; (*slNM)++; return;		/* remove 2nd */
	default: puts("inconsistent MS flags in AuxPut_hNF"); exit(0);}

	if(Sms) {if(Sms<3) (*slNM)--;} else (*slSM)--;	/* remove SL-entry */
	for(i=pos--;i<*slNF;i++) SLp[i-1]=SLp[i];
	(*slNF)--; (*slNB)-= *nu+2;
     }
}
/*	Alloc & read SL; go thru pi / pa; write SL + statistics
 */
void Add_Polya_2_Polyi(char *polyi,char *polya,char *polyo)
{    FILE *FI=fopen(polyi,"rb"), *FA=fopen(polya,"rb"), *FO; 
     FInfoList FIi, FIa, FIo;	  Along Ipos, Apos, HIpos, HApos;
     unsigned char ucI[NUC_Nmax], ucA[NUC_Nmax], *ucSL, *uc; int SLp[SL_Nmax];
     int IslNF, IslSM, IslNM, vI, nuI, i, d; unsigned Ili; Along IslNB; 
     int AslNF, AslSM, AslNM, vA=0, nuA=0, a, j; unsigned Ali; Along AslNB; 
     int slNF=0,slSM=0,slNM=0, slNB=0, slNP=0, Oli=0, v, nu, AmI=0, ms, tnb=0; 
     Init_FInfoList(&FIi); Init_FInfoList(&FIa);        if(!*polyi||!*polyo) {
       puts("With -pa you require -pi and -po or -di and -do");exit(0);}
     if(NULL==FI) {printf("Cannot open %s\n",polyi);exit(0);}
     if(NULL==FA) {printf("Cannot open %s\n",polya);exit(0);}
     if(NULL==(FO=fopen(polyo,"wb"))){printf("Cannot open %s",polyo);exit(0);}
     ucSL = (unsigned char *) malloc( SL_Nmax * CperR_MAX * sizeof(char) );
     assert(ucSL!=NULL); assert(!fgetc(FI)); assert(!fgetc(FA)); 
     Read_Bin_Info(FI,&i,&Ili,&IslNF,&IslSM,&IslNM,&IslNB,&FIi); 
     Read_Bin_Info(FA,&d,&Ali,&AslNF,&AslSM,&AslNM,&AslNB,&FIa); assert(d==i);
     HIpos=FTELL(FI); FSEEK(FI,0,SEEK_END); Ipos=FTELL(FI);
     HApos=FTELL(FA); FSEEK(FA,0,SEEK_END); Apos=FTELL(FA);
     printf("Data on %s:  %lld+%dsl  %lldb  (%dd)\n", polyi,
	(FIi.nNF-FIi.nSM)+(FIi.nNF-FIi.nNM),
	/* Islp= */ 2*IslNF-IslNM-IslSM,Ipos,i);
     printf("Data on %s:  %lld+%dsl  %lldb  (%dd)\n", polya,
	(FIa.nNF-FIa.nSM)+(FIa.nNF-FIa.nNM),
	/* Aslp= */ 2*AslNF-AslNM-AslSM,Apos,d);
/* printf("HIpos=%lld NB=%lld sl=%lld pos=%lld\n",HIpos,FIi.NB,IslNB,Ipos); */
     assert(HIpos+FIi.NB+IslNB==Ipos); FSEEK(FI,-IslNB,SEEK_CUR);
/* printf("HApos=%lld NB=%lld sl=%lld pos=%lld\n",HApos,FIa.NB,AslNB,Apos); */
     assert(HApos+FIa.NB+AslNB==Apos); FSEEK(FA,-AslNB,SEEK_CUR);
     a=0; if(a<AslNF) AuxGet_vn_uc(FA,&vA,&nuA,ucA);
     for(i=0;i<IslNF;i++)
     {	AuxGet_vn_uc(FI,&vI,&nuI,ucI); uc = &ucSL[ 2 + (SLp[slNF++]=slNB) ];
	while(a<AslNF)
	{   if(!(AmI=vA-vI)) if(!(AmI=nuA-nuI))
			       AmI=RIGHTminusLEFT(ucI,ucA,&nuI);
	    if(AmI<0)						   /* put A */
	    {	uc[-2]=vA; uc[-1]=nuA; for(j=0;j<nuA;j++)uc[j]=ucA[j];
		slNB+=2+nuA; if((ms=(*uc%4))) {if(ms<3) slNM++;} else slSM++;
		if((++a)<AslNF) AuxGet_vn_uc(FA,&vA,&nuA,ucA);
		uc = &ucSL[ 2 + (SLp[slNF++]=slNB) ];
	    }
	    else break;
	}
	if((a<AslNF)&&(AmI==0))					/* put I==A */
	{   uc[-2]=vI; uc[-1]=nuI; for(j=0;j<nuI;j++) uc[j]=ucI[j];
	    if((*ucI%4)!=(*ucA%4)) *uc = 3+4*(*uc/4);
	    if((ms=(*uc%4))) {if(ms<3) slNM++;} else slSM++;
	    if((++a)<AslNF) AuxGet_vn_uc(FA,&vA,&nuA,ucA);
	    tnb+=2+nuI;
	}
	else							   /* put I */
	{   uc[-2]=vI; uc[-1]=nuI; for(j=0;j<nuI;j++) uc[j]=ucI[j];
	    if((ms=(*uc%4))) {if(ms<3) slNM++;} else slSM++;
	}   slNB+=2+nuI;
     }     
     while(a<AslNF)						   /* put A */
     {	uc = &ucSL[ 2 + (SLp[slNF++]=slNB) ]; slNB+=2+nuA;
	uc[-2]=vA; uc[-1]=nuA; for(j=0;j<nuA;j++)uc[j]=ucA[j];
	if((ms=(*uc%4))) {if(ms<3) slNM++;} else slSM++;
	if((++a)<AslNF) AuxGet_vn_uc(FA,&vA,&nuA,ucA);
     }	assert(tnb+slNB==IslNB+AslNB);				 /* SL done */

     printf("SL: %dnf %dsm %dnm %db -> ",slNF,slSM,slNM,slNB);

     FSEEK(FI,HIpos,SEEK_SET); FSEEK(FA,HApos,SEEK_SET); Init_FInfoList(&FIo);
     FIo.nVmax=max(FIi.nVmax,FIa.nVmax);
     FIo.NUCmax=max(FIi.NUCmax,FIa.NUCmax); tnb=0;
     for(v=d+1;v<=FIo.nVmax;v++) for(nu=1;nu<=FIo.NUCmax;nu++) 
     if((FIo.NFnum[v][nu]=FIi.NFnum[v][nu]+FIa.NFnum[v][nu]))
     {	FIo.nNUC[v]++; Oli++; tnb+=FIo.NFnum[v][nu];
     }	FIo.nV=0; for(v=d+1;v<=FIo.nVmax;v++) if(FIo.nNUC[v]) FIo.nV++;
     fputc(0,FO); fputc(d,FO); fputc(FIo.nV,FO); fputc(FIo.nVmax,FO); 
     fputc(FIo.NUCmax,FO); fputUI(Oli,FO);
     j=InfoSize(0,Oli,&FIo)-5-sizeof(int); for(i=0;i<j;i++) fputc(0,FO);

     for(v=d+1;v<=FIo.nVmax;v++) 	if(FIo.nNUC[v])
     for(nu=1;nu<=FIo.NUCmax;nu++) 	if(FIo.NFnum[v][nu])
     {	int I_NF=FIi.NFnum[v][nu], A_NF=FIa.NFnum[v][nu], O_NF=0;
	int neq=0, peq=0, pi=0, pa=0, po=0;
	a=0; if(a<A_NF) {AuxGet_uc(FA,&nu,ucA); pa += 1+(((*ucA)%4)==3);}
     	for(i=0;i<I_NF;i++)
     	{   AuxGet_uc(FI,&nu,ucI); pi += 1+(((*ucI)%4)==3);
	    while(a<A_NF)
	    {   AmI=RIGHTminusLEFT(ucI,ucA,&nu);
	        if(AmI<0)					   /* put A */
	        {   po += 1+(((*ucA)%4)==3); 
		    AuxPut_hNF(FO,&v,&nu,ucA,&FIo,
			       &slNF,&slSM,&slNM,&slNB,ucSL,SLp); O_NF++; a++;
		    if(a<A_NF){AuxGet_uc(FA,&nu,ucA);pa += 1+(((*ucA)%4)==3);}
	        }
	        else break;
	    }
	    if((a<A_NF)&&(AmI==0))				/* put I==A */
	    {	int mm=10*(*ucI%4) + (*ucA%4); switch(mm) {
		case 00: peq++;		break;	case 11: peq++;		break;
		case 12: *ucI+=2; 	break;	case 13: *ucI+=2;peq++; break;
		case 21: *ucI+=1;	break;	case 22: peq++;		break;
		case 23: *ucI+=1;peq++; break;	case 31: peq++;		break;
		case 32: peq++;		break;	case 33: peq+=2;	break; 
		default: puts("inconsistens mirror flags");exit(0);} 
		AuxPut_hNF(FO,&v,&nu,ucI,&FIo,&slNF,&slSM,&slNM,&slNB,ucSL,
			   SLp); po += 1+(((*ucI)%4)==3); neq++; a++;
	    	if(a<A_NF) {AuxGet_uc(FA,&nu,ucA); pa += 1+(((*ucA)%4)==3);}
	    }
	    else
	    {	AuxPut_hNF(FO,&v,&nu,ucI,&FIo,&slNF,
			&slSM,&slNM,&slNB,ucSL,SLp); po += 1+(((*ucI)%4)==3);
	    }   O_NF++;
	}     
     	while(a<A_NF)						   /* put A */
     	{   O_NF++; po += 1+(((*ucA)%4)==3);
	    AuxPut_hNF(FO,&v,&nu,ucA,&FIo,&slNF,&slSM,&slNM,&slNB,ucSL,SLp);
	    ++a; if(a<A_NF) {AuxGet_uc(FA,&nu,ucA); pa += 1+(((*ucA)%4)==3);}
     	}   assert(pi+pa==peq+po);			  /* checksum(v,nu) */
	FIo.NFnum[v][nu]=O_NF; FIo.nNF+=O_NF; FIo.NB+=O_NF*nu;
	assert(O_NF+neq==I_NF+A_NF); 
	/*
	{static int list;printf("#%d v=%d nu=%d Inf=%d Anf=%d ",++list,v,nu,
	I_NF,A_NF);printf("Onf=%d   pi=%d pa=%d  po=%d\n",O_NF,pi,pa,po); 
	 */
     }							tnb=0;
     for(i=0;i<slNF;i++)					/* write SL */
     {	uc=&ucSL[SLp[i]+2]; v=uc[-2]; nu=uc[-1]; tnb+=nu+2;
	/* printf("#%d:  SLp=%d  v=%d  nu=%d\n",i,SLp[i],v,nu); */
	assert(uc[-2]<VERT_Nmax); fputc(uc[-2],FO); slNP+=1+(((*uc)%4)==3);
	fputc(nu,FO); for(j=0;j<nu;j++) fputc(uc[j],FO);
     }	assert(tnb==slNB);
	assert(slNP==2*slNF-slNM-slSM);
     printf("\nd=%d v%d v<=%d n<=%d vn%d  %lld %d %lld %lld  %d %d %d %d\n",
	d, FIo.nV, FIo.nVmax, FIo.NUCmax,Oli,
	FIo.nNF,FIo.nSM,FIo.nNM,FIo.NB,	 slNF, slSM, slNM, slNB);
/*     for(v=d+1;v<=FIo.nVmax;v++) if(FIo.nNUC[v])	     
     {	i=0; printf("%d %d\n",v,FIo.nNUC[v]);	     
	for(nu=1;nu<=FIo.NUCmax;nu++) if(FIo.NFnum[v][nu])
	{   printf("%d %d%s",nu,FIo.NFnum[v][nu],  
		(++i<FIo.nNUC[v]) ? "  " : "\n");
	}
     }	*/					

     FSEEK(FO,5+sizeof(int),SEEK_SET); 
				fputUI(FIo.nNF,FO);    /* write info */
     fputUI(FIo.nSM,FO); fputUI(FIo.nNM,FO); fputUI(FIo.NB,FO);
     fputUI(slNF,FO); fputUI(slSM,FO); fputUI(slNM,FO); fputUI(slNB,FO); 
     for(v=d+1;v<=FIo.nVmax;v++) if(FIo.nNUC[v])
     {  i=0; for(nu=1;nu<=FIo.NUCmax;nu++) if(FIo.NFnum[v][nu]) i++;
        fputc(v,FO); fputc(i,FO);                              /* v #nuc(v) */
        for(nu=1;nu<=FIo.NUCmax;nu++) if(FIo.NFnum[v][nu])
        {   fputc(nu,FO); fputUI(FIo.NFnum[v][nu],FO);           /* nuc #NF */
        }
     }

     printf("Writing %s: %lld+%dsl %lldm+%ds %lldb",polyo,2*FIo.nNF-FIo.nNM
	-FIo.nSM,slNP,/*tnb=*/ FIo.nNF-FIo.nNM-FIo.nSM,FIo.nSM,FIo.NB+slNB);
     /* if(tnb>99)
     {	long long tnp=(2*FIo.nNF-FIo.nNM-FIo.nSM)/1000; tnp*=tnp; 
	tnp/=(2*tnb); printf("   [p^2/2m=%ldM]",tnp);
     }*/
     Print_Expect(&FIo); puts("");
     assert(ferror(FI)==0); fclose(FI); 
     assert(ferror(FA)==0); fclose(FA);
     assert(ferror(FO)==0); fclose(FO);	
}


/*   ==================	    AffineNF modifications	==================  */

void UCnf_2_ANF(int *d, int *v, int *nuc, unsigned char *uc,          /* IN */
              Long NF[POLY_Dmax][VERT_Nmax], int *MS)                /* OUT */
{    int i,off;UCnf2vNF(d,v,nuc,uc,NF,MS);off=NF[0][*v-1];  /* offset not 1 */
     for(i=1;i<*d;i++)assert(NF[i][*v-1]==off);        /* but last column ! */
     
     for(i=0;i<*d;i++) {int c; for(c=0;c<i;c++) assert(NF[i][c]==0);
	for(c=i;c<*v;c++) NF[i][c]-=off;}	assert(*MS%4==1);
}
void ANF_2_ucNF(PolyPointList *P, VertexNumList *V, EqList *E,	      /* IN */
     		 int *NV, int *nUC, unsigned char *UC)		     /* OUT */
{    int vb,vo,vnuc,vbmin, vone=0, MSone; Long V_NF[POLY_Dmax][VERT_Nmax];
     Make_ANF(P,V,E,V_NF);		    /* Affine V_NF */
        vone=BminOff(V_NF,&P->n,&V->nv,&vo,&vbmin);
#ifdef	USE_UNIT_ENCODE
	if(vone) UNIT_Bmin2BaseUCn(&P->n,&V->nv,&vbmin,&vb,&vnuc); else
#endif
	Bmin2BaseUCn(&P->n,&V->nv,&vbmin,&vb,&vnuc);
     
     *nUC = vnuc;   *NV = V->nv;  MSone=1;

     if(AuxNFLptr!=NULL)			    /* base:byte statistics */
     {	if(AuxNFLptr->Xdif<vbmin) AuxNFLptr->Xdif=vbmin;
	if(AuxNFLptr->Xnuc<*nUC) AuxNFLptr->Xnuc=*nUC;
	if(vbmin>BASE_MAX) 
	{   printf("WARNING MS=%d  v=%d vb=%d vnuc=%d\n",
		1,V->nv,vbmin,vnuc);
	    VPrint(&P->n,&V->nv,V_NF);puts("============= V_NF");
	    puts("BASE_MAX exceeded"); exit(0);
	}
     }
#ifdef	USE_UNIT_ENCODE
     MSone += 4*vone;
#endif
     AuxVnf2ucNF(V_NF, &P->n, &V->nv, &vo, &vb, &vnuc, &MSone, UC);

#ifdef TEST_UCnf		/* One=K3[3/F]:  2 1 0 1 1 1 2  */
     {	Long tNF[POLY_Dmax][VERT_Nmax]; int tMS, i, j, *d=&P->n, err=0;
	static int pn; 
	pn++;		UCnf_2_ANF(d, NV, nUC, UC,   tNF, &tMS);
	for(i=0;i<*d;i++) for(j=0;j<*NV;j++) 
	    err+=(tNF[i][j]!=V_NF[i][j]);
	if(err) err*=1000;
	if(err) 
	{   int base;
#ifdef	USE_UNIT_ENCODE
	    if(MSone>3) UNIT_NUCtoBase(&P->n,NV,nUC,&base); else
#endif
	    NUCtoBase(&P->n,NV,nUC,&base);
	    printf("\n#%d  MS=%d  100*e(V)-e(F) = %d\n",pn,1,err);
	    for(i=0;i<*nUC;i++)printf("%d ",UC[i]);
	    printf("= UC[%d]  base=%d  vo=%d\n",*nUC,base,vo);
	    VPrint(&P->n,&V->nv,V_NF);puts("============= V_NF");
	    VPrint(&P->n,NV,tNF);puts("============= t_NF");	exit(0); 
	}
     }
#endif
     assert((*UC % 4)==1);
}

int  Add_ANF_to_List(PolyPointList *_P, VertexNumList *_V, EqList *_E,
                    NF_List *_L)
{    unsigned char UC[NB_MAX]; int nUC, NV, NewNF; 
#ifdef INCREMENTAL_TIME
     long long cpuT=clock(), incT=cpuT - _L->CLOCK; 
     if(incT>0) _L->IP_Time+=incT;
     _L->CLOCK=cpuT;
#endif
     if((_E->ne)>_L->F) _L->F=_E->ne;           /* check vertex numbers */
     if((_V->nv)>_L->V) _L->V=_V->nv;
     if((_E->ne > VERT_Nmax)||(_V->nv > VERT_Nmax))
     {  printf("Increase VERT_Nmax:  f=%d v=%d\n",_E->ne,_V->nv);exit(0);
     }
     if(_P->np>_L->Nmax) _L->Nmax=_P->np;
     if(_P->np<_L->Nmin) _L->Nmin=_P->np;
     _L->nNF++; if(_L->SL)_L->nSLNF++;  

     ANF_2_ucNF(_P, _V, _E,  &NV, &nUC, UC); 

     NewNF=ucNF_Sort_Add(&NV, &nUC, UC, _L);  /* 1::new::cont. */

#ifdef INCREMENTAL_TIME
     cpuT=clock(), incT=cpuT - _L->CLOCK; 
     if(incT>0) _L->NF_Time+=incT;
     _L->CLOCK=cpuT;
#endif
     return NewNF;
}

void Gen_Ascii_to_Binary(CWS *W, PolyPointList *P, 
		     char *dbin, char *polyi, char *polyo)
{ NF_List *_NFL=(NF_List *) malloc(sizeof(NF_List)); 
  VertexNumList V; EqList F;     assert(_NFL!=NULL);
  if(!(*polyo)) {
    puts("You have to specify an output file via -po in -a-mode!\n"); 
    printf("For more help use option '-h'\n");
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
    
    Find_Equations(P,&V,&F);
    if (Add_ANF_to_List(P,&V,&F,_NFL)) if (outFILE!=stdout){
      int i, j;
      for(i=0;i<W->nw;i++)     {        
        fprintf(outFILE,"%ld ",W->d[i]);
        for(j=0;j<W->N;j++) fprintf(outFILE,"%ld ",W->W[i][j]);
        if(i+1<W->nw) fprintf(outFILE," "); else fprintf(outFILE,"\n");   } 
      fflush(0);}  }
  Write_List_2_File(polyo,_NFL); 
  free(_NFL);
}
