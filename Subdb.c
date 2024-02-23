#include "Global.h"
#include "Subpoly.h"
#include "Rat.h"
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

void Small_Make_Dual(PolyPointList *_P, VertexNumList *_V, EqList *_E){
  int i;
  EqList AE=*_E;
  VNL_to_DEL(_P, _V, _E);
  assert(EL_to_PPL(&AE, _P, &_P->n));
  _V->nv=_P->np;
  for (i=0;i<_P->np;i++) _V->v[i]=i;
}

void Polyi_2_DBo(char *polyi,char *dbo)
{    char *dbnames = (char *) malloc(1+strlen(dbo)+File_Ext_NCmax), *fx;
     FILE *F=fopen(polyi,"rb"), *Finfo, *Fv, *Fsl;
     time_t Tstart=time(NULL);          FInfoList L; 		UPint tNF=0;
     int d, v, nu, i, j, list_num, sl_nNF, sl_SM, sl_NM, sl_NB; Along tNB=0;

     if (!*polyi) {puts("With -do you require -pi or -di and -pa"); exit(0);}

     if(F==NULL) {printf("Input file %s not found\n",polyi); exit(0);}
     strcpy(dbnames,dbo); fx=&dbnames[strlen(dbo)+1]; 
     strcat(dbnames,".info"); Finfo=fopen(dbnames,"w"); assert(Finfo!=NULL);
     printf("Read %s (",polyi); fflush(stdout);

     Init_FInfoList(&L);                        /* start reading the file */
     d=fgetc(F); assert(d==0);	/* for(i=0;i<d;i++) fgetc(F); */
     if(d){printf("Recursion depth %d forbidden in DB !!\n",d);Finfo=stdout;}
     
     d=fgetc(F); L.nV=fgetc(F); L.nVmax=fgetc(F); L.NUCmax=fgetc(F);
     list_num=fgetUI(F); 
     L.nNF=fgetUI(F); L.nSM=fgetUI(F); L.nNM=fgetUI(F);   L.NB=fgetUI(F); 
     sl_nNF=fgetUI(F); sl_SM=fgetUI(F); sl_NM=fgetUI(F);     sl_NB=fgetUI(F); 

     for(i=0;i<L.nV;i++)
     {  v=fgetc(F); L.nNUC[v]=fgetc(F);   /* read #nuc's per #Vert */
        for(j=0;j<L.nNUC[v];j++)                  
        {   L.NFnum[v][nu=fgetc(F)]=fgetUI(F);  /* read nuc and #NF(v,nu)*/
            tNF+=L.NFnum[v][nu]; tNB+=nu*L.NFnum[v][nu];
        }
     }  assert(tNF==L.nNF); assert(0== (unsigned int)(tNB-L.NB)); L.NB=tNB;

     printf("%lldpoly +%dsl %lldb)  write %s.* (%d files)  ",
        2*L.nNF-L.nSM-L.nNM, 2*sl_nNF-sl_SM-sl_NM, L.NB+sl_NB,
	dbo,L.nV+1+(sl_nNF>0)); fflush(stdout);

      fprintf(Finfo,                                          /* write Finfo */
        "%d  %d %d %d  %d  %lld %d %lld %lld  %d %d %d %d\n\n",
        d, L.nV, L.nVmax, L.NUCmax,list_num,
        L.nNF,L.nSM,L.nNM,L.NB,  sl_nNF, sl_SM, sl_NM, sl_NB);

     for(v=d+1;v<=L.nVmax;v++) if(L.nNUC[v])                 /* honest info */
     {  i=0; fprintf(Finfo,"%d %d\n",v,L.nNUC[v]);           /*  v  #nuc's  */
        for(nu=1;nu<=L.NUCmax;nu++) if(L.NFnum[v][nu])
        {   fprintf(Finfo,"%d %d%s",nu,L.NFnum[v][nu],    /* nuc #NF(v,nuc) */
                (++i<L.nNUC[v]) ? "  " : "\n");
        }
     }                                          if(Finfo==stdout) exit(0);
     if(ferror(Finfo)) {printf("File error in %s\n",dbnames);exit(0);}
     fclose(Finfo); fflush(stdout);

     for(v=d+1;v<=L.nVmax;v++) if(L.nNUC[v])         /* write  honest polys */
     {  char ext[4]={'v',0,0,0};
        ext[1]='0' + v / 10; ext[2]='0' + v % 10;  
        strcpy(fx,ext); Fv=fopen(dbnames,"wb"); assert(Fv!=NULL); 

        for(nu=1;nu<=L.NUCmax;nu++) if(L.NFnum[v][nu])
        {   int vnuNB=nu*L.NFnum[v][nu]; 
            for(i=0;i<vnuNB;i++) fputc(fgetc(F),Fv);
        }
        if(ferror(Fv)) {printf("File error in %s\n",dbnames);exit(0);} 
	fclose(Fv);
     }

     if(sl_nNF)                                  /* write  sublattice polys */
     {  strcpy(fx,"sl"); Fsl=fopen(dbnames,"wb");
        assert(Fsl!=NULL); for(i=0;i<sl_NB;i++) fputc(fgetc(F),Fsl);
        if(ferror(Fsl)) {printf("File error in %s\n",dbnames);exit(0);}
	fclose(Fsl);
     }

     printf("done (%ds)\n",(int) difftime(time(NULL),Tstart));

     if(ferror(F)) {printf("File error in %s\n",polyi);exit(0);} fclose(F);
}

void Init_DB(NF_List *_NFL){    
  /* Read the database, create RAM_poly;
     for given, nv, nuc the matching is as follows:
     DB:  |0|1|...|B-1|B|...|2B|...|((n-1)/B)*B|...|n-2|n-1|
     RAM:             |0  |   1|...| (n-1)/B)-1|
     (n...nNF[nv][nuc], B...BLOCK_LENGTH, the offsets are DB->Fv_pos[v][nuc]
     and DB->RAM_pos[v][nuc], respectively; 
     each entry |x| corresponds to nuc unsigned characters)  */

  time_t Tstart=time(NULL); 
  char *dbname = (char *) malloc(1+strlen(_NFL->dbname)+File_Ext_NCmax), *fx;
  DataBase *DB=&_NFL->DB;
  int d, v, nu, i, j, list_num, sl_nNF, sl_SM, sl_NM, sl_NB, 
    RAM_pos=0; Along RAM_size=0;

  strcpy(dbname,_NFL->dbname);
     printf("Reading data-base %s: ",dbname);
  strcat(dbname,".info");
  fx=&dbname[strlen(_NFL->dbname)+1];  

  /* read the info-file: */
  DB->Finfo=fopen(dbname,"r");
  assert(DB->Finfo!=NULL);
  fscanf(DB->Finfo, "%d  %d %d %d  %d  %lld %d %lld %lld  %d %d %d %d",
          &d,   &DB->nV, &DB->nVmax, &DB->NUCmax,   &list_num,
      &DB->nNF, &DB->nSM, &DB->nNM, &DB->NB, &sl_nNF, &sl_SM, &sl_NM, &sl_NB);
  printf("%lld+%dsl %lldnf %lldb", 
         2*(DB->nNF)-DB->nSM-DB->nNM,
         2*sl_nNF-sl_SM-sl_NM, DB->nNF+sl_nNF, DB->NB+sl_NB);
  /* if( _FILE_OFFSET_BITS < 64 ) assert(DB->NB <= LONG_MAX);	Along */
  if(_NFL->d) assert(d==_NFL->d); else _NFL->d=d;		

  for (v=1;v<VERT_Nmax;v++){ 
    DB->nNUC[v]=0;
    for(nu=0; nu<NUC_Nmax; nu++) DB->NFnum[v][nu]=0;}

  for(i=0; i<DB->nV; i++){
    fscanf(DB->Finfo,"%d",&v);
    fscanf(DB->Finfo,"%d",&(DB->nNUC[v]));
    for(j=0;j<DB->nNUC[v];j++){
      fscanf(DB->Finfo,"%d", &nu);
      fscanf(DB->Finfo,"%d", &(DB->NFnum[v][nu]));
      RAM_size+=nu*((DB->NFnum[v][nu]-1)/BLOCK_LENGTH);   } } 

  if(ferror(DB->Finfo)) {printf("File error in %s\n",dbname); exit(0);}
  fclose(DB->Finfo); 
  fflush(stdout);				assert(RAM_size<=INT_MAX);

  DB->RAM_NF=(unsigned char *) malloc(RAM_size);
  assert(DB->RAM_NF!=NULL); 

  /* read the DB-files and create RAM_NF: */
  for (v=2;v<=DB->nVmax;v++) if(DB->nNUC[v]){
    char ext[4]={'v',0,0,0};
    ext[1]='0' + v / 10; ext[2]='0' + v % 10;  
    strcpy(fx,ext); 
    DB->Fv[v]=fopen(dbname,"rb"); 
    assert(DB->Fv[v]!=NULL); 
    FSEEK(DB->Fv[v],0,SEEK_END);

    DB->Fv_pos[v][0]=0;
    for (nu=0;nu<=DB->NUCmax;nu++){
      DB->RAM_pos[v][nu]=RAM_pos;
      for (i=0; i<(DB->NFnum[v][nu]-1)/BLOCK_LENGTH; i++){
        FSEEK(DB->Fv[v], DB->Fv_pos[v][nu]+(i+1)*BLOCK_LENGTH*nu, SEEK_SET);
        for (j=0; j<nu; j++) DB->RAM_NF[RAM_pos++]=fgetc(DB->Fv[v]); }
      DB->Fv_pos[v][nu+1]=DB->Fv_pos[v][nu]+DB->NFnum[v][nu]*nu;    }    }

  printf("  done (%ds)\n",(int) difftime(time(NULL),Tstart));
  fflush(stdout);
}

char Compare_Poly(int *nuc, unsigned char *uc1, unsigned char *uc2){
  /* uc1 always encodes a single poly, 
     uc2 might encode a single poly or a mirror pair;
     returns: 'm': if uc1 isn't in uc2, but its mirror is;
              'i': if uc1 is in uc2;
              'l': if uc1 is less than uc2;
              'g': if uc1 is greater than uc2; */
  switch (RIGHTminusLEFT(uc1,uc2,nuc)){
  case -1: return 'g';
  case  1: return 'l';
  case  0: {if ((*uc1 + *uc2) % 4 == 3) return 'm'; else return 'i';}
  default: puts("Sth. wrong in Compare_Poly!!!"); exit(0);} return 0;
}

int  Is_in_DB(int *nv, int *nuc, unsigned char *uc, NF_List *_NFL){
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

  int pos, min_pos=-1, max_RAM_pos, i, max_Fv_piece; Along Fv_pos;
  unsigned char Aux_poly[BLOCK_LENGTH*NUC_Nmax];

  DataBase *DB=&_NFL->DB;	if (!DB->NFnum[*nv][*nuc]) return 0;

  /* Analyse the position of uc w.r.t. DB->RAM_NF */
  max_RAM_pos=(DB->NFnum[*nv][*nuc]-1)/BLOCK_LENGTH;
  while (max_RAM_pos-min_pos>1){
    pos=(max_RAM_pos+min_pos)/2;
    switch
      (Compare_Poly(nuc, uc, 
                    &(DB->RAM_NF[DB->RAM_pos[*nv][*nuc]+(*nuc)*pos]))){
      case 'm': return 0;
      case 'i': return 1;
      case 'l': {max_RAM_pos=pos; continue;}
      case 'g': min_pos=pos;}}

  /* Look for uc in the database: */			Fv_pos=max_RAM_pos;
  if (Fv_pos==(DB->NFnum[*nv][*nuc]-1)/BLOCK_LENGTH) 
    max_Fv_piece=DB->NFnum[*nv][*nuc]-Fv_pos*BLOCK_LENGTH;
  else max_Fv_piece=BLOCK_LENGTH;
  min_pos=-1;
  if (FSEEK(DB->Fv[*nv], 
            DB->Fv_pos[*nv][*nuc]+Fv_pos*
		(Along) ( (*nuc)*BLOCK_LENGTH) , SEEK_SET)){
    printf("Error in fseek in Is_in_DB!"); exit(0);}
  for (i=0;i<(*nuc)*(max_Fv_piece);i++) Aux_poly[i]=fgetc(DB->Fv[*nv]); 
  while (max_Fv_piece-min_pos>1){
    pos=(max_Fv_piece+min_pos)/2;
    switch(Compare_Poly(nuc, uc, &(Aux_poly[(*nuc)*pos]))){
    case 'm': {return 0;}
    case 'i': {return 1;}
    case 'l': {max_Fv_piece=pos; continue;}
    case 'g': {min_pos=pos;}  }  }
  return 0;
}

void Add_Polya_2_DBi(char *dbi,char *polya,char *dbo)
{    FInfoList FIi, FIa, FIo;  	Along Apos, HIpos, HApos, Inp, tnb=0, tNF=0;
     unsigned char ucI[NUC_Nmax],ucA[NUC_Nmax],*ucSL=NULL,*uc;int SLp[SL_Nmax];
     int d, vI, nuI, IslNF, IslSM, IslNM, nu; unsigned Ili, u, Oli=0, IslNB;
     int v, vA, nuA, AslNF, AslSM, AslNM, i;  unsigned Ali, a;	Along AslNB;
     int s, slNF=0, slSM=0, slNM=0, slNB=0, slNP=0; 		UPint Anp;
     int AmI=00,ms, newout=strcmp(dbi,dbo) && (*dbo),j=1+strlen(SAVE_FILE_EXT);
     char *Ifx, *Ifn = (char *) malloc(j+strlen(dbi)+File_Ext_NCmax), *Ofx,
     	*Ofn = (char *) malloc(j+strlen(newout ? dbo : dbi)+File_Ext_NCmax);
     FILE *FI, *FA, *FO;  if(*polya==0) {puts("-pa file required"); exit(0);}
	Init_FInfoList(&FIi); Init_FInfoList(&FIa);
     strcpy(Ifn,dbi); Ifx=&Ifn[strlen(dbi)]; strcpy(Ifx,".info");
     if(*dbo==0)dbo=dbi;
     strcpy(Ofn,dbo); Ofx=&Ofn[strlen(dbo)];

     if(NULL==(FI=fopen(Ifn,"r")))   {printf("Cannot open %s",Ifn);exit(0);}
     if(NULL==(FA=fopen(polya,"rb"))){printf("Cannot open %s",polya);exit(0);}
     fscanf(FI,"%d%d%d%d%d%lld%d%lld %lld %d%d%d%d",&d,&i,&j,&nu,&Ili,
	&FIi.nNF,&FIi.nSM,&FIi.nNM,&FIi.NB,&IslNF,&IslSM,&IslNM,&IslNB);
     FIi.nV=i; FIi.nVmax=j; FIi.NUCmax=nu;
     /*	printf("%d %d %d %d %d %d %d %d %lld %d %d %d %d\n",
	d,FIi.nV,FIi.nVmax,FIi.NUCmax,Ili,FIi.nNF,FIi.nSM,
	FIi.nNM,FIi.NB,IslNF,IslSM,IslNM,IslNB); */
     for(i=0;i<FIi.nV;i++)
     {	fscanf(FI,"%d",&v); fscanf(FI,"%d",&j); FIi.nNUC[v]=j;
        for(j=0;j<FIi.nNUC[v];j++) 
	{   fscanf(FI,"%d",&nu); fscanf(FI,"%d",&FIi.NFnum[v][nu]);
	    tNF+=FIi.NFnum[v][nu];
	}
     }                                     assert(tNF==FIi.nNF); tNF=0;
     assert(!fgetc(FA));	/*  rd==0  (recursion depth::no aux-file)  */
     Read_Bin_Info(FA,&j,&Ali,&AslNF,&AslSM,&AslNM,&AslNB,&FIa);assert(d==j);
     strcpy(Ifx,".sl"); if(IslNF) assert(NULL != (FI=fopen(Ifn,"rb")));
     if((IslNB+AslNB))
     assert(NULL != (ucSL = (unsigned char *) malloc( (IslNB+AslNB) )));
     HApos=FTELL(FA); FSEEK(FA,0,SEEK_END); Apos=FTELL(FA);
     Inp=2*FIi.nNF-FIi.nSM-FIi.nNM; Anp=2*FIa.nNF-FIa.nSM-FIa.nNM;
     printf("Data on %s:  %lld+%dsl  %lldb  (%dd)\n", dbi,Inp,
	/* Islp= */ 2*IslNF-IslNM-IslSM,FIi.NB+slNB,d);
     printf("Data on %s:  %u+%dsl  %lldb  (%dd)\n", polya,Anp,
	/* Aslp= */ 2*AslNF-AslNM-AslSM, Apos,d);
     assert(HApos+FIa.NB+AslNB==Apos); FSEEK(FA,-AslNB,SEEK_CUR);
     s=0; if(s<AslNF) AuxGet_vn_uc(FA,&vA,&nuA,ucA);
     for(i=0;i<IslNF;i++)
     {	AuxGet_vn_uc(FI,&vI,&nuI,ucI); uc = &ucSL[ 2 + (SLp[slNF++]=slNB) ];
	while(s<AslNF)
	{   if(!(AmI=vA-vI)) if(!(AmI=nuA-nuI))
			       AmI=RIGHTminusLEFT(ucI,ucA,&nuI);
	    if(AmI<0)						   /* put A */
	    {	uc[-2]=vA; uc[-1]=nuA; for(j=0;j<nuA;j++)uc[j]=ucA[j];
		slNB+=2+nuA; if( (ms=(*uc%4)) ) {if(ms<3) slNM++;} else slSM++;
		if((++s)<AslNF) AuxGet_vn_uc(FA,&vA,&nuA,ucA);
		uc = &ucSL[ 2 + (SLp[slNF++]=slNB) ];
	    }
	    else break;
	}
	if((s<AslNF)&&(AmI==0))					/* put I==A */
	{   uc[-2]=vI; uc[-1]=nuI; for(j=0;j<nuI;j++) uc[j]=ucI[j];
	    if((*ucI%4)!=(*ucA%4)) *uc = 3+4*(*uc/4);
	    if( (ms=(*uc%4)) ) {if(ms<3) slNM++;} else slSM++;
	    if((++s)<AslNF) AuxGet_vn_uc(FA,&vA,&nuA,ucA);
	    tnb+=2+nuI;
	}
	else							   /* put I */
	{   uc[-2]=vI; uc[-1]=nuI; for(j=0;j<nuI;j++) uc[j]=ucI[j];
	    if( (ms=(*uc%4)) ) {if(ms<3) slNM++;} else slSM++;
	}   slNB+=2+nuI;
     }     
     while(s<AslNF)						   /* put A */
     {	uc = &ucSL[ 2 + (SLp[slNF++]=slNB) ]; slNB+=2+nuA;
	uc[-2]=vA; uc[-1]=nuA; for(j=0;j<nuA;j++)uc[j]=ucA[j];
	if( (ms=(*uc%4)) ) {if(ms<3) slNM++;} else slSM++;
	if((++s)<AslNF) AuxGet_vn_uc(FA,&vA,&nuA,ucA);
     }	assert(tnb+slNB==IslNB+AslNB);				 /* SL done */

     printf("SL: %dnf %dsm %dnm %db -> ",slNF,slSM,slNM,slNB);
     if(IslNF) 
     {	HIpos=FTELL(FI); assert(HIpos==IslNB); assert(!ferror(FI)); 
	fclose(FI); if(!newout) remove(Ifn);
     }							    /* SL file done */

     FSEEK(FA,HApos,SEEK_SET); Init_FInfoList(&FIo);
     FIo.nVmax=max(FIi.nVmax,FIa.nVmax);
     FIo.NUCmax=max(FIi.NUCmax,FIa.NUCmax); /* Tnb=0; */
     for(v=d+1;v<=FIo.nVmax;v++) for(nu=1;nu<=FIo.NUCmax;nu++) 
     if( (FIo.NFnum[v][nu]=FIi.NFnum[v][nu]+FIa.NFnum[v][nu]) )
     {	FIo.nNUC[v]++; Oli++; /* Tnb+=FIo.NFnum[v][nu]; */
     }	FIo.nV=0; for(v=d+1;v<=FIo.nVmax;v++) if(FIo.nNUC[v]) FIo.nV++;

     for(v=d+1;v<=FIo.nVmax;v++) 	if(FIo.nNUC[v])
     { 	char vxt[5]; strcpy(vxt,".v"); vxt[2]=v/10+'0'; vxt[3]=v%10+'0'; 
	vxt[4]=0; strcpy(Ofx,vxt); strcpy(Ifx,vxt);
	if(FIi.nNUC[v])
	{   if(!newout) {strcat(Ifx,SAVE_FILE_EXT); assert(!rename(Ofn,Ifn));}
	  if(NULL==(FI=fopen(Ifn,"rb"))){printf("Ifn %s failed",Ifn);exit(0);}
	} if(NULL==(FO=fopen(Ofn,"wb"))){printf("Ofn %s failed",Ofn);exit(0);}

     for(nu=1;nu<=FIo.NUCmax;nu++) 	if(FIo.NFnum[v][nu])
     {	unsigned int I_NF=FIi.NFnum[v][nu], A_NF=FIa.NFnum[v][nu], O_NF=0;
	UPint neq=0, peq=0, pa=0; Along pi=0, po=0;

	a=0; if(0<A_NF) {AuxGet_uc(FA,&nu,ucA); pa += 1+(((*ucA)%4)==3);}
     	for(u=0;u<I_NF;u++)
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
	}*/
     }
	if(FIi.nNUC[v])
	{   assert(!ferror(FI)); fclose(FI); if(!newout) remove(Ifn); 
     	}   assert(!ferror(FO)); fclose(FO);
     }							tnb=0;

     if(slNF)
     {	strcpy(Ofx,".sl"); assert(NULL != (FO=fopen(Ofn,"wb")));
     for(i=0;i<slNF;i++)					/* write SL */
     {	uc=&ucSL[SLp[i]+2]; v=uc[-2]; nu=uc[-1]; tnb+=nu+2;
	assert(uc[-2]<VERT_Nmax); fputc(uc[-2],FO); slNP+=1+(((*uc)%4)==3);
	fputc(nu,FO); for(s=0;s<nu;s++) fputc(uc[s],FO);
     }	assert(tnb==slNB);
	assert(slNP==2*slNF-slNM-slSM);
	assert(!ferror(FO)); fclose(FO);					
     }
     printf("\nd=%d v%d v<=%d n<=%d vn%d  %lld %d %lld %lld  %d %d %d %d\n",
	d, FIo.nV, FIo.nVmax, FIo.NUCmax,Oli,
	FIo.nNF,FIo.nSM,FIo.nNM,FIo.NB,	 slNF, slSM, slNM, slNB);

     strcpy(Ofx,".info"); assert(NULL != (FO=fopen(Ofn,"w")));
     fprintf(FO,                                           /* write FO.info */
        "%d  %d %d %d  %d  %lld %d %lld %lld  %d %d %d %d\n\n",
        d, FIo.nV, FIo.nVmax, FIo.NUCmax,Oli,
        FIo.nNF,FIo.nSM,FIo.nNM,FIo.NB,  slNF, slSM, slNM, slNB);

     for(v=d+1;v<=FIo.nVmax;v++) if(FIo.nNUC[v])             /* honest info */
     {  i=0; fprintf(FO,"%d %d\n",v,FIo.nNUC[v]);            /*  v  #nuc's  */
        for(nu=1;nu<=FIo.NUCmax;nu++) if(FIo.NFnum[v][nu])
        {   fprintf(FO,"%d %d%s",nu,FIo.NFnum[v][nu],  	  /* nuc #NF(v,nuc) */
                (++i<FIo.nNUC[v]) ? "  " : "\n");
        }
     }
     printf("Writing %s: %lld+%dsl %lldm+%ds %lldb",dbo,
	2*FIo.nNF-FIo.nNM-FIo.nSM,
	slNP,/*Tnb=*/ FIo.nNF-FIo.nNM-FIo.nSM,FIo.nSM,FIo.NB+slNB);
/*   if(tnb>99)
     {	long long tnp=(2*FIo.nNF-FIo.nNM-FIo.nSM)/10; tnp*=tnp; 
	tnp/=20; tnp/=Tnb; printf("   [p^2/2m=%ldk]",tnp);
     }
*/   Print_Expect(&FIo);
     puts("");
     assert(ferror(FA)==0); fclose(FA);
     assert(ferror(FO)==0); fclose(FO);	
}
int  Check_sl_order(int *v, int *nu, unsigned char *uc)
{    static int n, V, NU; static unsigned char UC[NUC_Nmax]; if(n)
     {	if(V>(*v)) return 0; if((V==(*v))&&(NU>(*nu))) return 0;
     	if((V==(*v))&&(NU==(*nu))&&(RIGHTminusLEFT(UC,uc,&NU)<=0)) return 0;
     }	V=*v; NU=*nu; for(n=0;n<NU;n++) UC[n]=uc[n]; return 1;
}
int  Check_hnf_order(int *v, int *nu, unsigned char *uc)
{    static int n, V, NU; static unsigned char UC[NUC_Nmax]; if(n)
     {	if(V>(*v)) return 0; if((V==(*v))&&(NU>(*nu))) return 0;
     	if((V==(*v))&&(NU==(*nu))&&(RIGHTminusLEFT(UC,uc,&NU)<=0)) return 0;
     }	V=*v; NU=*nu; for(n=0;n<NU;n++) UC[n]=uc[n]; return 1;
}
void Print_NF(FILE *F,int *d, int *v,Long NF[POLY_Dmax][VERT_Nmax])
{    int i,j; fprintf(F,"%d %d\n",*d,*v); for(i=0;i<*d;i++)
     for(j=0;j<*v;j++)fprintf(F,"%d%s",(int)NF[i][j],(*v==j+1) ? "\n" : " ");
}
void Print_Missing_Mirror(int *d, int *v, int *nu, unsigned char *uc, 
			  PolyPointList *_P)
{    int I,J,MS; Long NF[POLY_Dmax][VERT_Nmax];
     VertexNumList V; EqList E;
     UCnf2vNF(d,v,nu,uc,NF,&MS); MS %= 4; _P->n=*d; _P->np=*v;
     for(I=0;I<*v;I++) for(J=0;J<*d;J++) _P->x[I][J]=NF[J][I];
     if(MS==2)	Print_NF(outFILE,d,v,NF);
     else if(MS==1)
       {	IP_Check(_P,&V,&E); Small_Make_Dual(_P,&V,&E);
                Make_Poly_NF(_P,&V,&E,NF); Print_NF(outFILE,d,&(V.nv),NF);
       }
     else {puts("Only use Print_Missing_Mirror for MM!");exit(0);}
}

/*   cF->{1::c 2::C (extended output)}  cF->{-1::M (missing mirrors)}   */
void Check_NF_Order(char *polyi,char *dbi, int cF, PolyPointList *_P)/* 1=MM */
{    FILE *F=NULL; FInfoList L; 		/* time_t Tstart=time(NULL); */
     unsigned int rd, i, j, list_num, tln=0, tSM=0; int d, nu, v, si;
     int sl_nNF, sl_SM, sl_NM, sl_NB; Along tNF=0, tNB=0, tNM=0, SLpos, Hpos=0;
     char *Ifx=NULL, *Ifn = (char *) malloc(1+strlen(dbi)+File_Ext_NCmax);
     if((*polyi)&&(*dbi)) puts("only give one of -pi FILE or -di FILE");
     if((*polyi==0)&&(*dbi==0)) puts("I need one of: -pi FILE or -di FILE");
     if((*polyi==0)+(*dbi==0)!=1) exit(0);	

     if(*polyi)						/*   read INFO PART */
     {	printf("Checking consistency of Aux/InFile %s\n",polyi); 
	F=fopen(polyi,"rb"); if(F==NULL) {puts("File not found");exit(0);}
     	Init_FInfoList(&L);                       /* start reading the file */

    	rd=fgetc(F); if(rd>127) rd=128*(rd-128)+fgetc(F);
	assert(rd <= MAX_REC_DEPTH);    printf("rd=%d: ",rd);
	for(i=0;i<rd;i++) printf(" %d",fgetc(F));
	if(rd) puts("");
     
     	d=fgetc(F); L.nV=fgetc(F); L.nVmax=fgetc(F); L.NUCmax=fgetc(F);
     	list_num=fgetUI(F); 
	printf("d=%d  nV=%d nVmax=%d NUCmax=%d #lists=%d ",
	    d,L.nV,L.nVmax,L.NUCmax,list_num); fflush(stdout);	
     	L.nNF=fgetUI(F); L.nSM=fgetUI(F); L.nNM=fgetUI(F); L.NB=fgetUI(F);
     	sl_nNF=fgetUI(F); sl_SM=fgetUI(F); sl_NM=fgetUI(F);  sl_NB=fgetUI(F); 

     	for(i=0;i<L.nV;i++)
     	{   v=fgetc(F); L.nNUC[v]=fgetc(F);   
	    if(cF==2) printf("v%dn%d: ",v,L.nNUC[v]);
	    fflush(stdout);
            for(j=0;j<L.nNUC[v];j++)                  
            {   L.NFnum[v][nu=fgetc(F)]=fgetUI(F); /* read nuc and #NF(v,nu)*/
                tNF+=L.NFnum[v][nu]; 
		if(cF==2) printf("%d:%d ",nu,L.NFnum[v][nu]); 
		tln++; tNB+=L.NFnum[v][nu]*(Along)nu;
            }	if(cF==2) puts("");
        }   assert( 0 == (unsigned int) (tNB-L.NB)  );  L.NB=tNB;
     	if(cF==2) printf("np=%lld+%dsl %lldb  %d files\n",2*L.nNF-L.nSM-L.nNM, 
	    2*sl_nNF-sl_SM-sl_NM,L.NB+sl_NB,L.nV+1+(sl_nNF>0)); 
	if(cF==2) printf("%lldnf %dsm %lldnm %lldb   sl: %d %d %d %d\n",
	    L.nNF,L.nSM,L.nNM,L.NB,sl_nNF,sl_SM,sl_NM,sl_NB); fflush(stdout);
     }	
     if(*dbi)
     {	printf("Checking consistency of DataBase %s:\n",dbi);
     	strcpy(Ifn,dbi); Ifx=&Ifn[strlen(dbi)]; strcpy(Ifx,".info");
	F=fopen(Ifn,"r"); if(F==NULL) {puts("Info File not found");exit(0);}
     	Init_FInfoList(&L);                       /* start reading the file */
	fscanf(F,"%d%d%d%d%d%lld%d%lld %lld %d%d%d%d",
        &d,&i,&j,&nu,&list_num,&L.nNF,&L.nSM,&L.nNM,&L.NB,
        &sl_nNF,&sl_SM,&sl_NM,&sl_NB);   L.nV=i; L.nVmax=j; L.NUCmax=nu;
	printf("d=%d  nV=%d nVmax=%d NUCmax=%d #lists=%d  ",
	    d,L.nV,L.nVmax,L.NUCmax,list_num);
     	printf("np=%lld+%dsl %lldb  %d files\n",2*L.nNF-L.nSM-L.nNM,
	    2*sl_nNF-sl_SM-sl_NM,L.NB+sl_NB,L.nV+1+(sl_nNF>0));
	printf("%lldnf %dsm %lldnm %lldb   sl: %d %d %d %d\n",
	    L.nNF,L.nSM,L.nNM,L.NB,sl_nNF,sl_SM,sl_NM,sl_NB); fflush(stdout);

     	for(i=0;i<L.nV;i++)
     	{   fscanf(F,"%d",&v); fscanf(F,"%d",&j); L.nNUC[v]=j;
	    if(cF==2) printf("v%dn%d: ",v,L.nNUC[v]);
	    fflush(stdout);
            for(j=0;j<L.nNUC[v];j++) 
            {	fscanf(F,"%d",&nu); fscanf(F,"%d",&L.NFnum[v][nu]); tln++;
            	tNF+=L.NFnum[v][nu];  
		if(cF==2) printf("%d:%d ",nu,L.NFnum[v][nu]);
            }	if(cF==2) puts("");
     	}			assert(tln==list_num);
     } printf("#hNF=%lld sum=%lld %s\n",L.nNF,tNF,
	(tNF==L.nNF) ? "o.k." : "Error");
     if(tNF!=L.nNF) exit(0); assert(!ferror(F));  tNF=0;
     if(tln!=list_num){printf("ERROR: #li=%d != %d\n",list_num,tln);exit(0);}

     {	/* long long np=2*L.nNF-L.nSM-L.nNM, pp2m=L.nNF-L.nSM-L.nNM; 
	j=2*sl_nNF-sl_SM-sl_NM; */
     	printf("np=%lld+%dsl  %lldnf %dsm %lldnm %lldb  sl: %d %d %d %d\n",
	    2*L.nNF-L.nSM-L.nNM,2*sl_nNF-sl_SM-sl_NM,L.nNF,L.nSM,L.nNM,L.NB, 
	    sl_nNF,sl_SM,sl_NM,sl_NB);
	fflush(stdout); /* if(pp2m) tln=(np/10)*(np/10)/20/pp2m;else tln=0;*/ 
     }
     printf("sl: ");					/* check SL PART */
     if(*polyi)
     {	Hpos=FTELL(F); FSEEK(F,-sl_NB,SEEK_END); SLpos=FTELL(F);
	if(SLpos-Hpos==L.NB) printf("NB o.k. ");
     }
     else if(sl_NB)
     {	strcpy(Ifx,".sl"); fclose(F);
        if(NULL==(F=fopen(Ifn,"rb"))){printf("Open %s failed",Ifn);exit(0);}
     }	else puts("no .sl file");
     for(si=0;si<sl_nNF;si++)
     {	unsigned char uc[NUC_Nmax]; v=fgetc(F); assert(v<=VERT_Nmax); 
	nu=fgetc(F); /* assert(nu<=L.NUCmax); */
	AuxGet_uc(F,&nu,uc); assert(!ferror(F));
        if((*uc%4)==0)tSM++; 
	else if((*uc%4)<3)tNM++;
	Test_ucNF(&d,&v,&nu,uc,_P); 
	assert(Check_sl_order(&v,&nu,uc));
     }  if(sl_NB&&(*dbi))
     {  SLpos=FTELL(F); assert(sl_NB==SLpos); printf("NB o.k. ");
     }	printf("sm=%d=%d nm=%d=%lld",sl_SM,tSM,sl_NM,tNM); fflush(stdout);
     if((sl_SM!=(int)tSM)||(sl_NM!=(int)tNM)){puts("ERROR!!");exit(0);} 
     tSM=tNM=0; /* if(tln>1)printf("  p^2/2m=%ldkCY",tln); */
     Print_Expect(&L);
    
     printf("\nv:"); if(*polyi) FSEEK(F,Hpos,SEEK_SET);	   /* check h-order */
     if(cF<0) puts("");
     for(v=d+1;v<=L.nVmax;v++)    /* if(cF) */ if(L.nNUC[v])
     {  Along nbsum=0; if(cF>0){printf(" %d",v);fflush(stdout);}  if(*dbi) 
     	{   char vxt[5]; strcpy(vxt,".v"); vxt[2]=v/10+'0'; vxt[3]=v%10+'0'; 
	    vxt[4]=0; strcpy(Ifx,vxt); fclose(F);
            if(NULL==(F=fopen(Ifn,"rb"))){printf("Ifn %s failed",Ifn);exit(0);}
	}
        for(nu=1;nu<=L.NUCmax;nu++)      if(L.NFnum[v][nu])
	{   unsigned char uc[NUC_Nmax];  for(i=0;i<L.NFnum[v][nu];i++)
            {   AuxGet_uc(F,&nu,uc); 	 Check_hnf_order(&v,&nu,uc);
		switch(*uc%4)
		{ case 0: tSM++; break; case 1:; case 2:tNM++; case 3:;}

		if(cF<0)if((*uc%4)%3)Print_Missing_Mirror(&d,&v,&nu,uc,_P);

	    }	nbsum+=nu*L.NFnum[v][nu];
	}
        if(*dbi){Hpos=FTELL(F); assert(Hpos==nbsum);} assert(!ferror(F));
     }	
     printf("  sm=%d nm=%lld",tSM,tNM); printf("  order o.k.");
     if((L.nSM!=(int)tSM)||(L.nNM!=tNM)) puts("    CheckSum ERROR!");
     if(cF<0) return;

     printf("\nv:"); if(*polyi) FSEEK(F,Hpos,SEEK_SET);	      /* check h-NF */
     for(v=d+1;v<=L.nVmax;v++)        	if(L.nNUC[v])
     {  int nbsum=0; printf(" %d",v);fflush(stdout);	if(*dbi) 
     	{   char vxt[5]; strcpy(vxt,".v"); vxt[2]=v/10+'0'; vxt[3]=v%10+'0'; 
	    vxt[4]=0; strcpy(Ifx,vxt); fclose(F);
            if(NULL==(F=fopen(Ifn,"rb"))){printf("Ifn %s failed",Ifn);exit(0);}
	}
        for(nu=1;nu<=L.NUCmax;nu++)      if(L.NFnum[v][nu])
	{   unsigned char uc[NUC_Nmax];  for(i=0;i<L.NFnum[v][nu];i++)
            {   AuxGet_uc(F,&nu,uc); 	 Test_ucNF(&d,&v,&nu,uc,_P);
		switch(*uc%4)
		{ case 0: tSM++; break; case 1:; case 2:tNM++; case 3:;}
	    }	nbsum+=nu*L.NFnum[v][nu];
	}
        if(*dbi){Hpos=FTELL(F); assert(Hpos==nbsum);} assert(!ferror(F));
     }	printf("  NF o.k.\n");	   assert(!ferror(F)); fclose(F);
}

#define SUBTRACT_H_FROM_SL

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

void Reduce_Aux_File(char *polyi,char *polys,char *dbsub,char *polyo)
{    FILE *FI=fopen(polyi,"rb"), *FS, *FO=fopen(polyo,"wb"); 
     FInfoList FIi, FIs, FIo; 
     Along Ipos,Spos=00, HIpos,HSpos=00, Opos,HOpos, IslNB, SslNB, tnb=0;
     unsigned char ucI[NUC_Nmax], ucS[NUC_Nmax], *ucSL, *uc; int SLp[SL_Nmax];
     int IslNF, IslSM, IslNM, db=0, vI, nuI, j, i; unsigned Ili, u;  UPint Inp;
     int SslNF, SslSM, SslNM, dv=0, vS, nuS, s, d; unsigned Sli;    Along Snp;
     int slNF=0, slSM=0, slNM=0, slNB=0, slNP=0, SmI=00, nu,ms, v; UPint Oli=0;
     char *Sfx=NULL, *Sfn = (*polys) ? (char *) NULL :
	(char *) malloc(1+strlen(dbsub)+File_Ext_NCmax);   
     Along HIPli[VERT_Nmax][NUC_Nmax], HSPli[VERT_Nmax][NUC_Nmax];

     if((*polys==0) != (*dbsub==0)) db=(*dbsub!=0); else
     {	printf("Need ONE of ps=%s and ds=%s\n",polys,dbsub); exit(0); }
     if(!*polyi||!*polyo) {
       puts("With -ps or -ds you have to specify I/O files via -pi and -po");
       exit(0);}

     if(NULL==FI) {printf("Cannot open %s",polyi);exit(0);}
     if(NULL==FO) {printf("Cannot open %s",polyo);exit(0);}
     ucSL = (unsigned char *) malloc( SL_Nmax * CperR_MAX * sizeof(char) );
     assert(ucSL!=NULL);  Init_FInfoList(&FIs); 

     if(db)
     {	unsigned tln=0; Along tNF=0; strcpy(Sfn,dbsub);
	Sfx=&Sfn[strlen(dbsub)]; strcpy(Sfx,".info");  FS=fopen(Sfn,"r"); 
	if(FS==NULL) {puts("Info File not found");exit(0);}
     	{Along sNF,sNM;  /* start reading the file */
	fscanf(FS,"%d%d%d%d%d%lld%d%lld %lld %d%d%d%lld",&d,&i,&j,&nu,&Sli,
            &sNF,&FIs.nSM,&sNM,&FIs.NB,&SslNF,&SslSM,&SslNM,&SslNB); 
            FIs.nNF=sNF; FIs.nNM=sNM; } FIs.nV=i; FIs.nVmax=j; FIs.NUCmax=nu;
     	for(i=0;i<FIs.nV;i++)
     	{   fscanf(FS,"%d",&v); fscanf(FS,"%d",&j); FIs.nNUC[v]=j;
            for(j=0;j<FIs.nNUC[v];j++) 
            {  fscanf(FS,"%d",&nu);fscanf(FS,"%d",&FIs.NFnum[v][nu]);
	       tln++;  tNF+=FIs.NFnum[v][nu];
            }
	}	assert(tln==Sli); s=d;
     	assert(tNF==FIs.nNF); assert(!ferror(FS)); tNF=0;
     	if(tln!=Sli){printf("ERROR: #li=%d != %d\n",Sli,tln);exit(0);}
     }	
     else
     {	if(NULL==(FS=fopen(polys,"rb")))
	{printf("Cannot open %s",polys);exit(0);} 
     	if(fgetc(FS)) {puts("don't subtract aux files!");exit(0);} 
     	Read_Bin_Info(FS,&s,&Sli,&SslNF,&SslSM,&SslNM,&SslNB,&FIs);
     	HSpos=FTELL(FS); FSEEK(FS,0,SEEK_END); Spos=FTELL(FS);
     	assert(HSpos+FIs.NB+SslNB==Spos); FSEEK(FS,-SslNB,SEEK_CUR);
     }
     d=fgetc(FI); if(d>127) d=128*(d-128)+fgetc(FI); Init_FInfoList(&FIi); 
     if(d) printf("RecursionDepth=%d on %s\n",d,polyi);
     if(d < 128) fputc(d,FO); else {fputc(128+(d)/128,FO); fputc(d%128,FO);} 
     for(i=0;i<d;i++) {j=fgetc(FI); fputc(j,FO);} 
     Read_Bin_Info(FI,&d,&Ili,&IslNF,&IslSM,&IslNM,&IslNB,&FIi); assert(d==s);
     HIpos=FTELL(FI); FSEEK(FI,0,SEEK_END); Ipos=FTELL(FI);
     Inp=2*FIi.nNF-FIi.nSM-FIi.nNM; Snp=2*FIs.nNF-FIs.nSM-FIs.nNM;

     printf("Data on %s:  %d+%dsl  %lldb  (%dd)\n", polyi,Inp,
	2*IslNF-IslNM-IslSM,Ipos,s);
     printf("Data on %s:  %lld+%dsl  %lldb  (%dd)\n", db ? dbsub : polys,Snp,
	2*SslNF-SslNM-SslSM,db ? FIs.NB+SslNB : Spos,d);
     assert(HIpos+FIi.NB+IslNB==Ipos); FSEEK(FI,-IslNB,SEEK_CUR);

     if(db && SslNF)
     {	strcpy(Sfx,".sl"); fclose(FS);
        if(NULL==(FS=fopen(Sfn,"rb"))){printf("Open %s failed",Sfn);exit(0);}
     }
     s=0; if(0<SslNF) AuxGet_vn_uc(FS,&vS,&nuS,ucS);
     for(v=0;v<IslNF;v++)
     {	AuxGet_vn_uc(FI,&vI,&nuI,ucI); uc = &ucSL[ 2 + (SLp[slNF++]=slNB) ];
	while(s<SslNF)
	{   if(!(SmI=vS-vI)) if(!(SmI=nuS-nuI))
			       SmI=RIGHTminusLEFT(ucI,ucS,&nuI);
	    if(SmI<0)						  /* next S */
	    {	if((++s)<SslNF) AuxGet_vn_uc(FS,&vS,&nuS,ucS);
	    }
	    else break;
	}
	if((s<SslNF)&&(SmI==0))					    /* I==S */
	{   ms=(*ucS%4); if( (ms%3) && ((*ucI%4)!=ms) )		/* put I==S */
	    {	uc[-2]=vI; uc[-1]=nuI; for(nu=0;nu<nuI;nu++) uc[nu]=ucI[nu];
		*uc=(3-ms)+4*(*uc/4); slNM++; slNB+=2+nuI; tnb+=2+nuI; 
	    }
	    else slNF--;
	    if((++s)<SslNF) AuxGet_vn_uc(FS,&vS,&nuS,ucS);
	}
	else							   /* put I */
	{   uc[-2]=vI; uc[-1]=nuI; for(i=0;i<nuI;i++) uc[i]=ucI[i];
	    if( (ms=(*uc%4)) ) {if(ms<3) slNM++;} else slSM++; slNB+=2+nuI;
	}   
     }			/* assert(tnb+slNB==IslNB+AslNB); */	 /* SL done */

     printf("SL: %dnf %dsm %dnm %db -> ",slNF,slSM,slNM,slNB);	tnb=0;
							fflush(stdout);

     if(!db) FSEEK(FS,HSpos,SEEK_SET);   FSEEK(FI,HIpos,SEEK_SET); Opos=HIpos;
     
     for(v=d+1;v<=FIi.nVmax;v++) if(FIi.nNUC[v])	      /* init HIPli */
     for(nu=1;nu<=FIi.NUCmax;nu++)	if(FIi.NFnum[v][nu])
     {	HIPli[v][nu]=Opos; Opos+=nu*FIi.NFnum[v][nu]; 
     }	assert(Opos==HIpos+FIi.NB);				Opos=HSpos;
     for(v=d+1;v<=FIs.nVmax;v++) if(FIs.nNUC[v])  	      /* init HSPli */
     for(nu=1;nu<=FIs.NUCmax;nu++)	if(FIs.NFnum[v][nu])
     {	if(db) if(v>dv) {Opos=0; dv=v;}
	HSPli[v][nu]=Opos; Opos+=nu*FIs.NFnum[v][nu]; 
     }  if(!db) assert(Opos==HSpos+FIs.NB); dv=0;	 Init_FInfoList(&FIo);

     for(v=d+1;v<=FIi.nVmax;v++)   if(FIi.nNUC[v])   /* FIo.NFnum[v][nu]?=0 */
     for(nu=1;nu<=FIi.NUCmax;nu++) if(FIi.NFnum[v][nu])
     {	if(FIs.NFnum[v][nu])			     /* check for add polys */
     	{   volatile unsigned int *_nf=&FIo.NFnum[v][nu]; 
						 /* (*_nf)=#NFnum(v,n) ?= 0 */
     	    unsigned int n, I_NF=FIi.NFnum[v][nu], S_NF=FIs.NFnum[v][nu];
	    if(db) if(v>dv) 
	    {	char vxt[5]; strcpy(vxt,".v"); vxt[2]=v/10+'0';
		vxt[3]=v%10+'0'; vxt[4]=0; strcpy(Sfx,vxt); 
		assert(!ferror(FS));fclose(FS); if(NULL==(FS=fopen(Sfn,"rb")))
		{printf("%s open failed",Sfn);exit(0);} dv=v;
	    }
	    FSEEK(FI,HIPli[v][nu],SEEK_SET); FSEEK(FS,HSPli[v][nu],SEEK_SET);
	    u=0; if(0<S_NF) AuxGet_uc(FS,&nu,ucS);
     	    for(n=0;n<I_NF;n++)
     	    {   AuxGet_uc(FI,&nu,ucI);
	        while(u<S_NF)
	        {   SmI=RIGHTminusLEFT(ucI,ucS,&nu);
	            if(SmI<0) 
	            {	if((++u)<S_NF)AuxGet_uc(FS,&nu,ucS);	   /* get S */
		    }
	    	    else break;
	    	}
	        if((u<S_NF)&&(SmI==0))			      /* found I==S */
		{   ms=(*ucS)%4; if((ms%3) && (((*ucI)%4)!=ms)) (*_nf)++;
		    if((++u)<S_NF)AuxGet_uc(FS,&nu,ucS);
		}
		else (*_nf)++;
		if(*_nf) break;
	    }     
     	}
     	else FIo.NFnum[v][nu]=FIi.NFnum[v][nu];
	if(FIo.NFnum[v][nu]) 
	{   FIo.nVmax=max(FIo.nVmax,v); FIo.NUCmax=max(FIo.NUCmax,nu); 
	    FIo.nNUC[v]++; Oli++;
	}
     }	dv=0;
     for(v=d+1;v<=FIo.nVmax;v++) if(FIo.nNUC[v]) FIo.nV++;

     fputc(d,FO); fputc(FIo.nV,FO); fputc(FIo.nVmax,FO); 
     fputc(FIo.NUCmax,FO); fputUI(Oli,FO); Opos=FTELL(FO); 
     j=InfoSize(0,Oli,&FIo)-5-sizeof(int); for(i=0;i<j;i++) fputc(0,FO);

     HOpos=FTELL(FO); /* printf("Opos=%d   HOpos=%d\n\n\n",Opos,HOpos); */

     for(v=d+1;v<=FIo.nVmax;v++) 	if(FIo.nNUC[v])
     for(nu=1;nu<=FIo.NUCmax;nu++) 	if(FIo.NFnum[v][nu])
     {	unsigned I_NF=FIi.NFnum[v][nu], S_NF=FIs.NFnum[v][nu], O_NF=0, n;
	UPint neq=0, peq=0, pi=0, po=0;
        if(v>tnb) 
	{   printf(" %d",v); tnb=v; fflush(stdout);	if(FIs.nNUC[v])	if(db) 
	    {	char vxt[5]; strcpy(vxt,".v"); vxt[2]=v/10+'0';
		vxt[3]=v%10+'0'; vxt[4]=0; strcpy(Sfx,vxt); 
		assert(!ferror(FS));fclose(FS); if(NULL==(FS=fopen(Sfn,"rb")))
		{printf("%s open failed",Sfn);exit(0);} 
	    }
	}
        FSEEK(FI,HIPli[v][nu],SEEK_SET); FSEEK(FS,HSPli[v][nu],SEEK_SET);
	u=0; if(u<S_NF) AuxGet_uc(FS,&nu,ucS);
     	for(n=0;n<I_NF;n++)
     	{   AuxGet_uc(FI,&nu,ucI); pi += 1+(((*ucI)%4)==3);
	    while(u<S_NF)
	    {   SmI=RIGHTminusLEFT(ucI,ucS,&nu);
	        if(SmI<0)					   /* get S */
	        {   u++; if(u<S_NF) AuxGet_uc(FS,&nu,ucS);
	        }
	        else break;
	    }
	    if((u<S_NF)&&(SmI==0))				 /* put I-S */
	    {	int k, mm=10*(*ucI%4) + (*ucS%4); switch(mm) {
		case 33: peq++; case 22:; case 11:; case 00: peq++;    	break;
		case 31:; case 32: *ucI-=(*ucS%4); peq++; 
		case 12:; case 21: 		neq--;	  po++; O_NF++; 
			FIo.nNM++; for(k=0;k<nu;k++) fputc(ucI[k],FO); 	break;
		case 13:; case 23: peq++; 				break;
		default: puts("inconsistens mirror flags");exit(0);} 
		neq++; u++; if(u<S_NF) AuxGet_uc(FS,&nu,ucS);
	    }
	    else
	    {	int k;for(k=0;k<nu;k++)fputc(ucI[k],FO);O_NF++;k=(*ucI)%4;
		po += 1+(k==3); if(k) {if(k<3)FIo.nNM++;} else FIo.nSM++;
	    }
	}   assert(pi-peq==po);	assert(O_NF+neq==I_NF);	  /* checksum(v,nu) */
	FIo.NFnum[v][nu]=O_NF; FIo.nNF+=O_NF; FIo.NB+=O_NF*nu;
     }							tnb=0;

#ifdef	SUBTRACT_H_FROM_SL
     uc = & (ucSL[SLp[i=0]]); dv=0;
     for(v=d+1;v<=FIs.nVmax;v++)   if(FIs.nNUC[v])    /* subtract S from SL */
     for(nu=1;nu<=FIs.NUCmax;nu++) if(FIs.NFnum[v][nu])
     {	while( (uc[0]<v) || ((uc[0]==v)&&(uc[1]<nu))) 
	{   uc = & (ucSL[SLp[++i]]); if(i>=slNF) goto END_SL;
	}				      /* go up to (v,nu) in SL-list */
	    if(db) if(v>dv) 
	    {	char vxt[5]; strcpy(vxt,".v"); vxt[2]=v/10+'0';
		vxt[3]=v%10+'0'; vxt[4]=0; strcpy(Sfx,vxt); 
		assert(!ferror(FS));fclose(FS); if(NULL==(FS=fopen(Sfn,"rb")))
		{printf("%s open failed",Sfn);exit(0);} dv=v;
	    }

	if((uc[0]==v)&&(uc[1]==nu))	     /* read from file and sort out */
	{   FSEEK(FS,HSPli[v][nu],SEEK_SET); 
	    for(u=0;u<FIs.NFnum[v][nu];u++)
	    {	int HmSL; AuxGet_uc(FS,&nu,ucS);/* Test_ucNF(&d,&v,&nu,ucS);*/
		while(0 < (HmSL=RIGHTminusLEFT(&uc[2],ucS,&nu))) /* next SL */
		{   uc = & (ucSL[SLp[++i]]); if(i>=slNF) goto END_SL;
		    if((uc[0]!=v)||(uc[1]!=nu)) goto END_VN;
		}
		if(HmSL==0)				       /* remove SL */
		{   int k, sms=uc[2]%4, hms=(*ucS)%4;
		    switch(10*hms + sms){ 
		    case 31: case 32: case 00: case 11: case 22: case 33: 
		    	for(k=i+1;k<slNF;k++) SLp[k-1]=SLp[k];
			slNB -= nu+2;
			if(sms==0) --slSM; else if(sms<3) --slNM;
			--slNF;
			if(i>=slNF) goto END_SL;
			uc = & (ucSL[SLp[i]]);
			if((uc[0]!=v)||(uc[1]!=nu)) goto END_VN;
			break;
		    case 12: case 21: 					break;
		    case 13: case 23: ++slNM; uc[2]-=hms; 		break;
		    default: puts("inconsistent MS flags in SL-H"); exit(0);}
		}			       /* else (SL>H): hence next H */
	    }	assert(!ferror(FS));
	}
	else 	END_VN: assert( (uc[0]>v) || ((uc[0]==v)&&(uc[1]>nu)) );
     }		END_SL: ;
#endif

     for(i=0;i<slNF;i++)					/* write SL */
     {	uc=&ucSL[SLp[i]+2]; v=uc[-2]; nu=uc[-1]; tnb+=nu+2;
	/* printf("#%d:  SLp=%d  v=%d  nu=%d\n",i,SLp[i],v,nu); */
	assert(uc[-2]<VERT_Nmax); fputc(uc[-2],FO); slNP+=1+(((*uc)%4)==3);
	fputc(nu,FO); for(s=0;s<nu;s++) fputc(uc[s],FO);
     }	assert(tnb==slNB);
	assert(slNP==2*slNF-slNM-slSM);
     printf("\nd=%d v%d v<=%d n<=%d vn%d  %lld %d %lld %lld  %d %d %d %d\n",
	d, FIo.nV, FIo.nVmax, FIo.NUCmax,Oli,
	FIo.nNF,FIo.nSM,FIo.nNM,FIo.NB,	 slNF, slSM, slNM, slNB);

     FSEEK(FO,Opos,SEEK_SET); 
				fputUI(FIo.nNF,FO);    	      /* write info */
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
	-FIo.nSM,slNP,tnb=FIo.nNF-FIo.nNM-FIo.nSM,FIo.nSM,FIo.NB+slNB);
     /* if(tnb>99)
     {	long long tnp=(2*FIo.nNF-FIo.nNM-FIo.nSM)/1000; tnp*=tnp; 
	tnp/=(2*tnb); printf("   [p^2/2m=%ldM]",tnp);
     }	*/
     Print_Expect(&FIo);
     puts("");
     assert(ferror(FI)==0); fclose(FI); 
     assert(ferror(FS)==0); fclose(FS);
     assert(HOpos==FTELL(FO)); assert(ferror(FO)==0); fclose(FO);	
}
void Bin2a(char *polyi, int max, PolyPointList *_P)
{    FILE *F=fopen(polyi,"rb"); FInfoList L; UPint list_num,tNF=0; Along tNB=0;
     int d, v, s, sl_nNF, sl_SM, sl_NM, sl_NB, mc=0, MS, nu; unsigned i, j; 
     unsigned char uc[POLY_Dmax*VERT_Nmax]; VertexNumList V; EqList E;	
     Long NF[POLY_Dmax][VERT_Nmax];      		Init_FInfoList(&L);

     if(F==NULL) {printf("Input file %s not found\n",polyi); exit(0);}
     d=fgetc(F); assert(d==0);  /* for(i=0;i<d;i++) fgetc(F); */
     d=fgetc(F); L.nV=fgetc(F); L.nVmax=fgetc(F); L.NUCmax=fgetc(F);
     list_num=fgetUI(F); 
     L.nNF=fgetUI(F); L.nSM=fgetUI(F); L.nNM=fgetUI(F);   L.NB=fgetUI(F); 
     sl_nNF=fgetUI(F); sl_SM=fgetUI(F); sl_NM=fgetUI(F);     sl_NB=fgetUI(F); 

     for(i=0;i<L.nV;i++)
     {  v=fgetc(F); L.nNUC[v]=fgetc(F);   /* read #nuc's per #Vert */
        for(j=0;j<L.nNUC[v];j++)                  
        {   L.NFnum[v][nu=fgetc(F)]=fgetUI(F);  /* read nuc and #NF(v,nu)*/
            tNF+=L.NFnum[v][nu]; tNB+=L.NFnum[v][nu]*(Along)nu;
        }
     }	assert( 0 == (unsigned int)(tNB-L.NB) );  L.NB=tNB; assert(tNF==L.nNF);

     for(v=d+1;v<=L.nVmax;v++) if(L.nNUC[v])         /* write  honest polys */
     {	int I,J; for(nu=1;nu<=L.NUCmax;nu++) for(j=0;j<L.NFnum[v][nu];j++)
        {   for(s=0;s<nu;s++) uc[s]=fgetc(F);
	    UCnf2vNF(&d,&v,&nu,uc,NF,&MS); MS %= 4; _P->n=d; _P->np=v;
	    for(I=0;I<v;I++) for(J=0;J<d;J++) _P->x[I][J]=NF[J][I];
	    assert(Ref_Check(_P,&V,&E));

	    if(MS!=2)	    			      /* if(MS!=2) print NF */
	       if(!max || Poly_Max_check(_P,&V,&E)) 
		{   mc++; Print_NF(outFILE,&d,&v,NF);}
	    if(MS>1) 			  /* if(MS>1); print Mirror */
	      if(!max || Poly_Min_check(_P,&V,&E))
	        {   mc++; Small_Make_Dual(_P,&V,&E); 
	            Make_Poly_NF(_P,&V,&E,NF); Print_NF(outFILE,&d,&(V.nv),NF);
		}
        }
     }
     printf("np=%lld+%dsl  ",2*L.nNF-L.nSM-L.nNM,2*sl_nNF-sl_SM-sl_NM);
     printf(                                          /* write Finfo */
        "%dd  %dv<=%d n<=%d  %dnv  %lld %d %lld %lld  %d %d %d %d",
        d, L.nV, L.nVmax, L.NUCmax,list_num,
        L.nNF,L.nSM,L.nNM,L.NB,  sl_nNF, sl_SM, sl_NM, sl_NB);
     if(max) printf("  r-max=%d",mc);
     puts("");
}
void DB_fromVF_toVT(DataBase *DB,int vf,int vt)
{    int v,n; DB->nNF=0;
     for(v=vf;v<=vt;v++) for(n=0;n<NUC_Nmax;n++) DB->nNF+=DB->NFnum[v][n]; 
     if(!DB->nNF){fprintf(stderr,"No NF with %d<=v<=%d\n",vf,vt);exit(0);}
     while(0==DB->nNUC[vt]) {vt--;assert(vf<=vt);} DB->v=vf; DB->nVmax=vt;
}
void Bin2aDBsl(char *dbi, int max, int vf, int vt, PolyPointList *_P)
{    FILE *F; FInfoList L;
     int d, v, nu, i, j, list_num, mc=0, MS,
        sl_nNF, sl_SM, sl_NM, sl_NB, tSM=0, tNM=0;
     char *Ifx, *Ifn = (char *) malloc(1+strlen(dbi)+File_Ext_NCmax);
     Long NF[POLY_Dmax][VERT_Nmax]; 
     VertexNumList V;EqList E;
     strcpy(Ifn,dbi); Ifx=&Ifn[strlen(dbi)]; strcpy(Ifx,".info");
     F=fopen(Ifn,"r"); if(F==NULL) {puts("Info File not found");exit(0);}
     Init_FInfoList(&L);                       /* start reading the file */
     fscanf(F,"%d%d%d%d%d%lld%d%lld %lld %d%d%d%d",
	    &d,&i,&j,&nu,&list_num,&L.nNF,&L.nSM,&L.nNM,&L.NB,
	    &sl_nNF,&sl_SM,&sl_NM,&sl_NB);   L.nV=i; L.nVmax=j; L.NUCmax=nu;
     if(sl_NB && (vf==2)&&(vt==VERT_Nmax-1))
     {  strcpy(Ifx,".sl"); fclose(F);
        if(NULL==(F=fopen(Ifn,"rb"))){printf("Open %s failed",Ifn);exit(0);}
     }  
     else /* puts("no .sl file"); */ 
     {	DataBase *DB; Open_DB(dbi, &DB, 0); if((DB->v<vf)||(DB->nVmax>vt))
	    DB_fromVF_toVT(DB,vf,vt); /*  read only vf <= v <= vt :: */
	for(i=0; Read_H_poly_from_DB(DB,_P); i++) Print_PPL(_P,"");
	printf("#poly=%d\n",i); Close_DB(DB);
     }
     for(i=0;i<sl_nNF;i++)
     {  int I,J; unsigned char uc[NUC_Nmax]; v=fgetc(F); assert(v<=VERT_Nmax); 
        nu=fgetc(F); AuxGet_uc(F,&nu,uc); assert(!ferror(F));
        if((*uc%4)==0)tSM++; else if((*uc%4)<3)tNM++; 
	Test_ucNF(&d,&v,&nu,uc,_P); 
	
	UCnf2vNF(&d,&v,&nu,uc,NF,&MS); MS %= 4; _P->n=d; _P->np=v;
	for(I=0;I<v;I++) for(J=0;J<d;J++) _P->x[I][J]=NF[J][I];
	assert(Ref_Check(_P,&V,&E));
	if(MS!=2)                                 /* if(MS!=2) print NF */
	  if(!max || Poly_Max_check(_P,&V,&E)) 
	    {   mc++; Print_NF(outFILE,&d,&v,NF);}
	if(MS>1)                              /* if(MS>1); print Mirror */
	  if(!max || Poly_Min_check(_P,&V,&E))
	    {   mc++; Small_Make_Dual(_P,&V,&E); 
	        Make_Poly_NF(_P,&V,&E,NF); Print_NF(outFILE,&d,&(V.nv),NF);
	    }
     }  
     printf("np=%lld+%dsl  ",2*L.nNF-L.nSM-L.nNM,2*sl_nNF-sl_SM-sl_NM);
     printf(                                          /* write Finfo */
        "%dd  %dv<=%d n<=%d  %dnv ... %d %d %d %d",  /* %d %d %d %lld */ 
        d, L.nV, L.nVmax, L.NUCmax,list_num,
        /* L.nNF,L.nSM,L.nNM,L.NB, */ sl_nNF, sl_SM, sl_NM, sl_NB);
     if(max) printf("  r-max=%d",mc);
     puts("");
}
void Bin_2_ascii(char *polyi,char *dbin,int max,int vf,int vt,PolyPointList *P)
{    if(*polyi)     Bin2a(polyi, max, P);
     else if(*dbin) Bin2aDBsl(dbin, max, vf, vt, P);
     else puts("With -b[2a] you have to specify input via -pi or -di"); 
}     

/*  ======================================================================  */
/*  ==========                                                  ==========  */
/*  ==========    Hodge-database routines                       ==========  */
/*  ==========                                                  ==========  */
/*  ======================================================================  */

#define Hod_Dif_max 480
#define Hod_Min_max 251

#if (POLY_Dmax < 6)
void DB_to_Hodge(char *dbin, char *dbout, int vfrom, int vto, 
		 PolyPointList *_P){    

  /* Read the database, write the Hodge numbers */

  DataBase DB;
  VertexNumList V;
  Long VPM[EQUA_Nmax][VERT_Nmax];
  static EqList E;
  time_t Tstart;
  char *dbname = (char *) malloc(1+strlen(dbin)+File_Ext_NCmax), *fx;
  char *dbhname = (char *) malloc(6+strlen(dbout)+File_Ext_NCmax), *fhx;
  unsigned char uc_poly[NUC_Nmax];
  int d, v, nu, i, j, list_num, sl_nNF, sl_SM, sl_NM, sl_NB, dh,
    nnf_vd[VERT_Nmax][Hod_Dif_max+1], 
    nnf_v[VERT_Nmax];
  BaHo BH;
  FILE *Faux[Hod_Dif_max+1];
  FILE *Fvinfo;
  PolyPointList *_PD = (PolyPointList *) malloc(sizeof(PolyPointList));
  assert(_PD!=NULL);

  if (!*dbin||!*dbout){
    puts("You have to specify I/O database names via -di and -do");
    exit(0);}

  for(i=0;i<=Hod_Dif_max;i++) for(j=0;j<VERT_Nmax;j++) nnf_vd[j][i]=0;
  strcpy(dbname,dbin);
  strcpy(dbhname,dbout);
  strcat(dbname,".info");
  strcat(dbhname,".vinfo");
  fx=&dbname[strlen(dbin)+1]; 
  fhx=&dbhname[strlen(dbout)+1];
  Fvinfo=fopen(dbhname,"a");

  printf("Reading %s: ",dbname); fflush(0);

  /* read the info-file: */
  DB.Finfo=fopen(dbname,"r");
  assert(DB.Finfo!=NULL);
  fscanf(DB.Finfo, "%d  %d %d %d  %d  %lld %d %lld %lld  %d %d %d %d",
          &d,   &DB.nV, &DB.nVmax, &DB.NUCmax,   &list_num,
      &DB.nNF, &DB.nSM, &DB.nNM, &DB.NB, &sl_nNF, &sl_SM, &sl_NM, &sl_NB);
  printf("%lldp (%dsl) %lldnf %lldb\n", 
         2*(DB.nNF)-DB.nSM-DB.nNM+2*sl_nNF-sl_SM-sl_NM,
         2*sl_nNF-sl_SM-sl_NM, DB.nNF+sl_nNF, DB.NB+sl_NB);

  for (v=1;v<VERT_Nmax;v++){ 
    DB.nNUC[v]=0;
    for(nu=0; nu<NUC_Nmax; nu++) DB.NFnum[v][nu]=0;}

  for(i=0; i<DB.nV; i++){
    fscanf(DB.Finfo,"%d",&v);
    fscanf(DB.Finfo,"%d",&(DB.nNUC[v]));
    for(j=0;j<DB.nNUC[v];j++){
      fscanf(DB.Finfo,"%d", &nu);
      fscanf(DB.Finfo,"%d", &(DB.NFnum[v][nu]));   } }

  if(ferror(DB.Finfo)) {printf("File error in %s\n",dbname); exit(0);}
  fclose(DB.Finfo); 
  fflush(stdout);

  printf("Reading DB-files, calculating Hodge numbers, writing aux-files:\n");

  /* read the DB-files and calculate Hodge numbers */
  for (v=vfrom;v<=vto;v++) if(DB.nNUC[v]){
    int nd=0;
    char ext[4], aext[8];
    ext[0]='v'; ext[1]='0' + v / 10; ext[2]='0' + v % 10; ext[3]=0;
    aext[0]='v'; aext[1]='0' + v / 10; aext[2]='0' + v % 10; aext[3]='d'; 
    aext[4]=aext[5]=aext[6]=aext[7]=0;
    Tstart=time(NULL); 
    printf("v=%d: ", v); fflush(0);
    strcpy(fx,ext); 
    DB.Fv[v]=fopen(dbname,"rb"); 
    assert(DB.Fv[v]!=NULL); 
    for (nu=0;nu<=DB.NUCmax;nu++) for (i=0; i<DB.NFnum[v][nu]; i++){
      int mirror=0;
      int MS;
      unsigned char c;
      for (j=0; j<nu; j++) uc_poly[j]=fgetc(DB.Fv[v]);
      uc_nf_to_P(_P, &MS, &d, &v, &nu, uc_poly);  
      assert(Ref_Check(_P, &V, &E));
      assert(V.nv==v);
      Make_VEPM(_P,&V,&E,VPM);
      Complete_Poly(VPM,&E,V.nv,_P);
      Make_Dual_Poly(_P,&V,&E,_PD);
      RC_Calc_BaHo(_P,&V,&E,_PD,&BH);
      if (BH.h1[1]<BH.h1[2]) mirror=1;
      if (mirror) dh=BH.h1[2]-BH.h1[1];
      else dh=BH.h1[1]-BH.h1[2];
      if (!nnf_vd[v][dh]||(dh>250)){
	aext[4]='0'+dh/100; aext[5]='0'+(dh/10)%10;  aext[6]='0'+dh%10;  
	strcpy(fhx,aext); 
	Faux[dh]=fopen(dbhname,"ab"); }
      if (!nnf_vd[v][dh]) nd++; 
      c=BH.h1[2-mirror]; fputc(c,Faux[dh]);
      c=BH.mp/256+4*BH.mv; fputc(c,Faux[dh]);
      c=BH.mp%256; fputc(c,Faux[dh]);
      c=BH.np/256+4*BH.nv; fputc(c,Faux[dh]);
      c=BH.np%256; fputc(c,Faux[dh]);
      c=nu+64*mirror; fputc(c,Faux[dh]);
      for (j=0;j<nu;j++) fputc(uc_poly[j],Faux[dh]); 
      if (dh>250) fclose(Faux[dh]);
      (nnf_vd[v][dh])++;
      (nnf_v[v])++;    }

    if(ferror(DB.Fv[v])) {printf("File error in %s\n",dbname); exit(0);}
    fclose(DB.Fv[v]);
    for (dh=0;dh<=250;dh++) if (nnf_vd[v][dh]) {
      if(ferror(Faux[dh])) {printf("File error at dh=%d\n",dh); exit(0);}
      fclose(Faux[dh]);}

    fprintf(Fvinfo, "%d %d %d\n", v, nd, nnf_v[v]);
    for (dh=0;dh<=Hod_Dif_max;dh++) if (nnf_vd[v][dh]) 
      fprintf(Fvinfo, "%d %d ", dh, nnf_vd[v][dh]);
    fprintf(Fvinfo, "\n");
    printf(" %d NF (%ds)\n", nnf_v[v], (int) difftime(time(NULL),Tstart)); 
    fflush(0);  }
  free(_PD);
  fclose(Fvinfo);
}

void Sort_Hodge(char *dbaux, char *dbout){    
  /* Sort from v-chi-format to chi-h12-format */

  time_t Tstart;
  char *dbaname = (char *) malloc(6+strlen(dbaux)+File_Ext_NCmax), *fax;
  char *dbhname = 
    (char *) malloc(6+strlen(*dbout ? dbout : dbaux)+File_Ext_NCmax), *fhx;
  int v, i, j, dh, nd,
    nnf_d[Hod_Dif_max+1], nnf_vd[VERT_Nmax][Hod_Dif_max+1], 
    nnf_h[Hod_Min_max+1], nnf_v[VERT_Nmax];
  FILE *Fh[Hod_Min_max+1];
  FILE *Fchia, *Fvinfo, *Fhinfo;

  for(i=0;i<=Hod_Dif_max;i++){
    nnf_d[i]=0;
    for(j=0;j<VERT_Nmax;j++) nnf_vd[j][i]=0;}

  if(!*dbout) dbout=dbaux;
  strcpy(dbaname,dbaux);
  strcat(dbaname,".vinfo");
  strcpy(dbhname,dbout);
  strcat(dbhname,".hinfo");
  fhx=&dbhname[strlen(dbout)+1];
  fax=&dbaname[strlen(dbaux)+1];

  printf("Reading %s\n",dbaname); fflush(0);

  /* read the info-file: */
  Fvinfo=fopen(dbaname,"r");
  assert(Fvinfo!=NULL);
  while ((fscanf(Fvinfo, "%d",&v))!=EOF){
    fscanf(Fvinfo, "%d  %d", &nd, &(nnf_v[v]));
    /* printf("%d %d %d  ", v, nd, nnf_v[v]  ); fflush(0); */
    for (i=0;i<nd;i++) {
      fscanf(Fvinfo, "%d", &dh);
      fscanf(Fvinfo, "%d", &(nnf_vd[v][dh]));
      nnf_d[dh]+=nnf_vd[v][dh];}}
  if(ferror(Fvinfo)) {printf("File error in %s\n",dbaname); exit(0);}
  fclose(Fvinfo);

  printf("Sorting and writing Hodge-files:\n"); fflush(0);
  Tstart=time(NULL); 

  Fhinfo=fopen(dbhname,"w");
  assert(Fhinfo!=NULL);

  /* Sort the Hodge&Poly-Data */
  for (dh=0;dh<=Hod_Dif_max;dh++) if (nnf_d[dh]){
    char aext[8], hext[9];
    int h12, nh=0;
    unsigned char c, nuc;
	aext[0]='v'; aext[3]='d'; aext[4]='0'+dh/100;
	aext[5]='0'+(dh/10)%10; aext[6]='0'+dh%10; aext[7]=0;
	hext[0]='d'; hext[1]='0'+dh/100; hext[2]='0'+(dh/10)%10; 
	hext[3]='0'+dh%10; hext[4]='h'; hext[5]=hext[6]=hext[7]=hext[8]=0;

    printf("dh=%d: %dNF...", dh, nnf_d[dh]); fflush(0);
    for(j=0;j<=Hod_Min_max;j++) nnf_h[j]=0;

    for (v=2;v<VERT_Nmax;v++) if(nnf_vd[v][dh]){
      aext[1]='0'+v/10; aext[2]='0'+v%10;
      strcpy(fax,aext); 
      Fchia=fopen(dbaname,"rb");
      for(i=0;i<nnf_vd[v][dh];i++){
	h12=fgetc(Fchia);
	if (!nnf_h[h12]){
	  hext[5]='0'+h12/100; hext[6]='0'+(h12/10)%10; hext[7]='0'+h12%10;  
	  strcpy(fhx,hext); 
	  Fh[h12]=fopen(dbhname,"wb");
	  nh++;}
	nnf_h[h12]++; 
	c=fgetc(Fchia); fputc(c,Fh[h12]);
	if (v!=c/4) {printf("v=%d, hp.mv=%d", v, (int) (c/4)); exit(0);}
	/* for (j=0;j<4;j++) {c=fgetc(Fchia); fputc(c,Fh[h12]);}
	for (j=0;j<c%64;j++) fputc(fgetc(Fchia),Fh[h12]);}*/
        for (j=0;j<3;j++) fputc(fgetc(Fchia),Fh[h12]);
        nuc=fgetc(Fchia); fputc(nuc,Fh[h12]);
        c=fgetc(Fchia); 
        { int ms=c%4; if(ms) c+=3-ms; }
        fputc(c,Fh[h12]);
        for (j=0;j<nuc%64-1;j++) fputc(fgetc(Fchia),Fh[h12]);}
      if(ferror(Fchia)) {
	printf("File error in Fchia at dh=%d v=%d\n",dh,v); exit(0);}
      fclose(Fchia); }

    fprintf(Fhinfo, "%d %d %d\n", dh, nh, nnf_d[dh]);
    for (h12=0;h12<=Hod_Min_max;h12++) if (nnf_h[h12]) {
      fprintf(Fhinfo, "%d %d ", h12, nnf_h[h12]);
      nnf_d[dh]-=nnf_h[h12];
      if(ferror(Fh[h12])) {
	printf("File error at dh=%d h12=%d\n",dh, h12); exit(0);}
      fclose(Fh[h12]);} 
    fprintf(Fhinfo, "\n");
    if (nnf_d[dh]) 
      {printf("nnf_d[%d]!=sum nnf_dh[%d][h12]!", dh, dh); exit(0);}
    printf(" sorted\n");  }

  if(ferror(Fhinfo)) {printf("File error in Fhinfo\n"); exit(0);}
  fclose(Fhinfo);
  printf("  done (%ds)\n",(int) difftime(time(NULL),Tstart)); 
  fflush(stdout);
}

void Test_Hodge_db(char *dbname){

  time_t Tstart;
  char *filename = (char *) malloc(6+strlen(dbname)+File_Ext_NCmax), *fhx;
  int i, j, dh, h12, nh, nnf_sum,
    nnf_d[Hod_Dif_max+1], nnf_dh[Hod_Dif_max+1][Hod_Min_max+1];
  FILE *Fh;
  FILE *Fhinfo;

  for(i=0;i<=Hod_Dif_max;i++){
    nnf_d[i]=0;
    for(j=0;j<=Hod_Min_max;j++) nnf_dh[i][j]=0;}

  strcpy(filename,dbname);
  strcat(filename,".hinfo");
  fhx=&filename[strlen(dbname)+1];

  printf("Reading %s\n",filename); fflush(0);

  /* read the info-file: */
  Fhinfo=fopen(filename,"r");
  assert(Fhinfo!=NULL);
  while ((fscanf(Fhinfo, "%d",&dh))!=EOF){
    fscanf(Fhinfo, "%d  %d", &nh, &(nnf_d[dh]));
    nnf_sum=0;
    for (i=0;i<nh;i++) {
      fscanf(Fhinfo, "%d", &h12);
      fscanf(Fhinfo, "%d", &(nnf_dh[dh][h12]));
      nnf_sum+=nnf_dh[dh][h12];}
    if (nnf_sum!=nnf_d[dh]) {
      printf("nnf_d[%d]!=sum nnf_dh[%d][h12]!", dh, dh); exit(0);}}
  if(ferror(Fhinfo)) {printf("File error in %s\n",filename); exit(0);}
  fclose(Fhinfo);

  printf("Analysing Hodge-files:\n"); fflush(0);
  Tstart=time(NULL); 

  /* Analyse Hodge&Poly-Data */
  for (dh=0;dh<=Hod_Dif_max;dh++) if (nnf_d[dh]){
    char hext[9];
    hext[0]='d'; hext[1]='0'+dh/100; hext[2]='0'+(dh/10)%10; 
    hext[3]='0'+dh%10; hext[4]='h'; hext[5]=hext[6]=hext[7]=hext[8]=0;
    
    printf("dh=%d: %dNF...", dh, nnf_d[dh]); fflush(0);
    nnf_sum=0;
    
    for (h12=0;h12<=Hod_Min_max;h12++) if(nnf_dh[dh][h12]){
      /* unsigned char uc_poly[NUC_Nmax]; */
      int c1, /* c2, */ nuc;
      hext[5]='0'+h12/100; hext[6]='0'+(h12/10)%10; hext[7]='0'+h12%10;  
      strcpy(fhx,hext); 
      Fh=fopen(filename,"rb");
      assert (Fh!=0);
      while((c1=fgetc(Fh))!=EOF){
	nnf_sum++;
	/* c2= */ fgetc(Fh); 
	/* mv=c1/4; mp=(c1%4)*256+c2; */
	c1=fgetc(Fh); /* c2= */ fgetc(Fh); 
	/* nv=c1/4; np=(c1%4)*256+c2; */
	c1=fgetc(Fh);
	/* mirror=c1/64; */
	nuc=c1%64;
	for (j=0;j<nuc;j++) /* uc_poly[j]= */ fgetc(Fh);
	/* uc_nf_to_P(_P, &MS, &d, &mv, &nuc, uc_poly); */} 
      if(ferror(Fh)) {
	printf("File error in Fh at dh=%d h12=%d\n",dh,h12); exit(0);}
      fclose(Fh); }
    
    printf("%d\n",nnf_sum); }
    
    printf("  done (%ds)\n",(int) difftime(time(NULL),Tstart)); 
  fflush(stdout);
}

void Extract_from_Hodge_db(char *dbname, char *x_string, PolyPointList *_P){

  time_t Tstart;
  char c=*x_string, hext[9], com[35],
    *filename = (char *) malloc(6+strlen(dbname)+File_Ext_NCmax), *fhx;
  unsigned char uc_poly[NUC_Nmax];
  int E=998, H1=0, H2=0, M=0, V=0, N=0, F=0, L=1000, i=0, j, dh, h12, nh, 
    nnf_sum, nnf_d[Hod_Dif_max+1], nnf_dh[Hod_Dif_max+1][Hod_Min_max+1], 
    HD_from=0, HD_to=Hod_Dif_max, HM_from=0, HM_to=Hod_Min_max, c1, c2, nuc,
    mirror, nv, np, mv, mp, d=4, MS, is_poly, is_dual, true_H1, true_H2, 
    max_mv=34;
  FILE *Fh;
  FILE *Fhinfo;

  hext[0]='d'; hext[4]='h'; hext[8]=0;

  while(c){
    if(c=='E'){
      int neg=0;
      c=x_string[++i];
      if(c=='-') {neg=1; c=x_string[++i];}
      if((c-'0'>=0)&&(c-'0'<=9)) E=0;
      while((c-'0'>=0)&&(c-'0'<=9)) {E=10*E+c-'0'; c=x_string[++i];}
      if(E%2) {puts("The Euler number E number must be even!"); return;}
      if(neg) E*=-1;}
    else if(c=='H'){ 
      c=x_string[++i]; 
      while((c-'0'>=0)&&(c-'0'<=9)) {H1=10*H1+c-'0'; c=x_string[++i];}}
    else if(c==':'){
      c=x_string[++i];
      while((c-'0'>=0)&&(c-'0'<=9)) {H2=10*H2+c-'0'; c=x_string[++i];}}
    else if(c=='M'){
      c=x_string[++i];
      while((c-'0'>=0)&&(c-'0'<=9)) {M=10*M+c-'0'; c=x_string[++i];}}
    else if(c=='V'){
      c=x_string[++i];
      while((c-'0'>=0)&&(c-'0'<=9)) {V=10*V+c-'0'; c=x_string[++i];}}
    else if(c=='N'){
      c=x_string[++i];
      while((c-'0'>=0)&&(c-'0'<=9)) {N=10*N+c-'0'; c=x_string[++i];}}
    else if(c=='F'){
      c=x_string[++i];
      while((c-'0'>=0)&&(c-'0'<=9)) {F=10*F+c-'0'; c=x_string[++i];}}
    else if(c=='L'){
      c=x_string[++i];
      L=0;
      while((c-'0'>=0)&&(c-'0'<=9)) {L=10*L+c-'0'; c=x_string[++i];}}
    else {printf("`%c' is not valid input", c); return;}}

  if(H1&&H2){E=2*(H1-H2);}
  else if(H1&&(E!=998)) H2=H1-E/2;
  else if(H2&&(E!=998)) H1=H2+E/2;

  if( H1>491 || H2>491 || H1+H2>502 || (abs(E)!=998 && abs(E)>960) ){
    puts("#NF: 0   Note the range for Hodge numbers:");
    puts("         h11,h12<=491, h11+h12<=502, |E|<=960.");
    return;}

 if( (V&&(V<5||V>33)) || (F&&(F<5||F>33))){
    puts("#NF: 0   Note the range [5,33] for vertex/facet numbers!");
    return;}

  if( (M&&(M<6||M>680)) || (N&&(N<6||N>680)) ){
    puts("#NF: 0   Note the range [6,680] for point numbers!");
    return;}

  if(V&&(V<max_mv)) max_mv=V; 
  if(F&&(F<max_mv)) max_mv=F;
  if(M&&(M-1<max_mv)) max_mv=M-1;
  if(N&&(N-1<max_mv)) max_mv=N-1; 
 
  for(i=0;i<=Hod_Dif_max;i++){
    nnf_d[i]=0;
    for(j=0;j<=Hod_Min_max;j++) nnf_dh[i][j]=0;}

  strcpy(filename,dbname);
  strcat(filename,".hinfo");
  fhx=&filename[strlen(dbname)+1];

  /* printf("Reading %s\n",filename); fflush(0); */

  /* read the info-file: */
  Fhinfo=fopen(filename,"r");
  assert(Fhinfo!=NULL);
  while ((fscanf(Fhinfo, "%d",&dh))!=EOF){
    fscanf(Fhinfo, "%d  %d", &nh, &(nnf_d[dh]));
    nnf_sum=0;
    for (i=0;i<nh;i++) {
      fscanf(Fhinfo, "%d", &h12);
      fscanf(Fhinfo, "%d", &(nnf_dh[dh][h12]));
      nnf_sum+=nnf_dh[dh][h12];}
    if (nnf_sum!=nnf_d[dh]) {
      printf("nnf_d[%d]!=sum nnf_dh[%d][h12]!", dh, dh); exit(0);}}
  if(ferror(Fhinfo)) {printf("File error in %s\n",filename); exit(0);}
  fclose(Fhinfo);

  /* printf("Analysing Hodge-files:\n"); fflush(0); */
  Tstart=time(NULL); 

  /* Analyse Hodge&Poly-Data */

  nnf_sum=0;
  if (H1&&H2) {HM_from=(H1<H2 ? H1 : H2); HM_to=HM_from;}
  else if (H1||H2) {HM_to=max(H1,H2); HM_from=max(0,HM_to-Hod_Dif_max);}

  for (h12=HM_from; h12<=HM_to; h12++) {

    if(abs(E)!=998) {HD_from=abs(E/2); HD_to=abs(E/2);}
    else if ((H1||H2)&&(h12<HM_to)) {HD_from=HM_to-h12; HD_to=HM_to-h12;}
    else {HD_from=0; HD_to=Hod_Dif_max;}

    for (dh=HD_from; dh<=HD_to; dh++) if (nnf_dh[dh][h12]){
    
      hext[1]='0'+dh/100; hext[2]='0'+(dh/10)%10; hext[3]='0'+dh%10;
      hext[5]='0'+h12/100; hext[6]='0'+(h12/10)%10; hext[7]='0'+h12%10;      
      strcpy(fhx,hext); 
      Fh=fopen(filename,"rb");
      assert (Fh!=0);
      while((c1=fgetc(Fh))!=EOF){
	mv=c1/4; 
	if(mv>max_mv) break;
	c2=fgetc(Fh); 
	mp=(c1%4)*256+c2;
	c1=fgetc(Fh); c2=fgetc(Fh); 
	nv=c1/4; np=(c1%4)*256+c2;
	c1=fgetc(Fh);
	mirror=c1/64;
	nuc=c1%64;
	for (j=0;j<nuc;j++) uc_poly[j]=fgetc(Fh);

	if (mirror) {true_H1=h12; true_H2=dh+h12;}
	else {true_H2=h12; true_H1=dh+h12;}
	is_poly=0;
	if ((abs(E)==998)||(E==2*(true_H1-true_H2))) if (!H1||(true_H1==H1)) 
	  if (!H2||(true_H2==H2)) if(!M||(M==mp)) if (!V||(V==mv)) 
	    if (!N||(N==np)) if (!F||(F==nv)) is_poly=1;
	is_dual=0;
	if ((abs(E)==998)||(E==2*(true_H2-true_H1))) if (!H1||(true_H2==H1)) 
	  if (!H2||(true_H1==H2)) if(!M||(M==np)) if (!V||(V==nv)) 
	    if (!N||(N==mp)) if (!F||(F==mv)) is_dual=1;
	if (!is_poly&&!is_dual) continue;

	uc_nf_to_P(_P, &MS, &d, &mv, &nuc, uc_poly);

	if (is_poly){
	  /* if(!MS) puts("!MS");
	     if(!mirror) puts("!mirror"); */
	  if (++nnf_sum>L) {printf("Exceeded limit of %d\n",L); return;}
	  sprintf(com,"M:%d %d N:%d %d H:%d,%d [%d]",
		  mp,mv,np,nv,true_H1,true_H2,2*(true_H1-true_H2));
	  Print_PPL(_P,com);  }
	if (is_dual&&MS){
	  VertexNumList VNL; EqList EL;
	  Long NF[POLY_Dmax][VERT_Nmax];
	  if (++nnf_sum>L) {printf("Exceeded limit of %d\n",L); return;}
	  sprintf(com,"M:%d %d N:%d %d H:%d,%d [%d]",
		  np,nv,mp,mv,true_H2,true_H1,2*(true_H2-true_H1));
	  Find_Equations(_P,&VNL,&EL);
	  Small_Make_Dual(_P, &VNL, &EL);
	  Make_Poly_NF(_P, &VNL, &EL, NF);
	  Print_Matrix(NF, _P->n, VNL.nv, com);	}
      } 
      if(ferror(Fh)) {
	printf("File error in Fh at dh=%d h12=%d\n",dh,h12); exit(0);}
      fclose(Fh); }    
    }
  printf("#NF: %d\n",nnf_sum); 
  printf("  done (%ds)\n",(int) difftime(time(NULL),Tstart)); 
  fflush(stdout);
}

void Test_Hodge_file(char *filename, PolyPointList *_P){
  FILE *Ft=fopen(filename,"rb");
  unsigned char uc_poly[NUC_Nmax];
  
  int nuc, mirror, nv, np, mv, mp, j, c1, c2, d=4, MS;
  assert (Ft!=0);
  while((c1=fgetc(Ft))!=EOF){
    c2=fgetc(Ft); 
    mv=c1/4; mp=(c1%4)*256+c2;
    c1=fgetc(Ft); c2=fgetc(Ft); 
    nv=c1/4; np=(c1%4)*256+c2;
    c1=fgetc(Ft);
    mirror=c1/64; nuc=c1%64;
    for (j=0;j<nuc;j++) uc_poly[j]=fgetc(Ft);
    uc_nf_to_P(_P, &MS, &d, &mv, &nuc, uc_poly);  
    Print_PPL(_P,"");
    printf("nuc=%d mirror=%d mv=%d mp=%d nv=%d np=%d MS=%d\n", 
	   nuc,mirror,mv,mp,nv,np,MS); }
  fclose(Ft);
}
#endif

/*  =====================================================================  */
void Open_DB(char *dbin, DataBase **_DB, int info)
{    int i,j,v,nu; DataBase *DB; char *dbname,*fx,ext[4]; 
     if(*dbin==0){*_DB=NULL;return;}
     dbname = (char *) malloc(1+strlen(dbin)+File_Ext_NCmax);
     DB = (DataBase *) malloc(sizeof(DataBase)); assert(DB != NULL); *_DB=DB;
     strcpy(dbname,dbin); strcat(dbname,".info"); fx=&dbname[strlen(dbin)+1]; 
     if(info) {printf("Reading %s: ",dbname); fflush(0);}
     DB->Finfo=fopen(dbname,"r"); assert(DB->Finfo!=NULL);
     fscanf(DB->Finfo, "%d  %d %d %d  %d  %lld %d %lld %lld  %d %d %d %d", 
	&DB->d,&DB->nV,&DB->nVmax,&DB->NUCmax,&DB->list_num, &DB->nNF,&DB->nSM,
	&DB->nNM, &DB->NB, &DB->sl_nNF, &DB->sl_SM, &DB->sl_NM, &DB->sl_NB);
     if(info) printf("%lldp (%dsl) %lldnf %lldb\n", 
         2*(DB->nNF)-DB->nSM-DB->nNM+2*DB->sl_nNF-DB->sl_SM-DB->sl_NM, 2*
         DB->sl_nNF-DB->sl_SM-DB->sl_NM, DB->nNF+DB->sl_nNF,DB->NB+DB->sl_NB);
     for (v=1;v<VERT_Nmax;v++){ DB->nNUC[v]=0;
	for(nu=0; nu<NUC_Nmax; nu++) DB->NFnum[v][nu]=0;}
     for(i=0; i<DB->nV; i++){
	fscanf(DB->Finfo,"%d",&v); fscanf(DB->Finfo,"%d",&(DB->nNUC[v]));
    	for(j=0;j<DB->nNUC[v];j++){
      	    fscanf(DB->Finfo,"%d", &nu);
      	    fscanf(DB->Finfo,"%d", &(DB->NFnum[v][nu])); } }
     if(ferror(DB->Finfo)) {printf("File error in %s\n",dbname); exit(0);}
     fclose(DB->Finfo); ext[0]='v'; ext[3]=0; 		DB->v=DB->p=DB->nu=0;
     for(v=DB->d+1; v<=DB->nVmax; v++) if(DB->nNUC[v])
     {	ext[1]='0' + v / 10; ext[2]='0' + v % 10; strcpy(fx,ext);
	DB->Fv[v]=fopen(dbname,"rb"); assert(DB->Fv[v]!=NULL); 
	if(0==DB->v) {DB->v=v; while(0==DB->NFnum[v][DB->nu]) (DB->nu)++;}
     }
}
void Close_DB(DataBase *DB)
{    int v; if(DB==NULL) return; 
     for(v=1+DB->d; v<=DB->nVmax; v++) if(DB->nNUC[v]) {
	if(ferror(DB->Fv[v])) {printf("File error at v=%d\n",v); exit(0);}
	fclose(DB->Fv[v]);} free(DB);
}
#undef TEST
int  Read_H_ucNF_from_DB(DataBase *DB, unsigned char *uc)/* p=next read pos */
{    static Along totNF; int rest; assert(DB!=NULL);
     rest = DB->NFnum[DB->v][DB->nu] - DB->p; 
#ifdef	TEST
     {static int vNF, xvNF, nuNF=DB->NFnum[DB->v][DB->nu]; 
     if(totNF==0){int i,j;for(i=1;i<=DB->NUCmax;i++)vNF+=DB->NFnum[DB->v][i];
	for(i=DB->v;i<=DB->nVmax;i++) for(j=1;j<=DB->NUCmax;j++)
	    xvNF+=DB->NFnum[i][j]; 
	if(xvNF!=DB->nNF) {printf("xvNF=%d  nNF=%d",xvNF,DB->nNF);exit(0);}
     	xvNF=vNF;
     }}
#endif
     assert(0<=rest); if(rest==0)			/* search next v/nu */
     {	int v=DB->v, nu=DB->nu; 			   /* search nu(v): */
	while(DB->nu < DB->NUCmax) if(DB->NFnum[DB->v][++(DB->nu)]) break;
	if((nu==DB->nu)||(DB->NFnum[DB->v][DB->nu]==0))	       /* search v: */
	{ if(DB->v < DB->nVmax)   
	  { while(0==DB->nNUC[++DB->v]) assert(DB->v < DB->nVmax); DB->nu=0;
	    while(DB->nu < DB->NUCmax) if(DB->NFnum[DB->v][++(DB->nu)]) break;
          }   else { assert(totNF==DB->nNF); return 0; } }
	if((v<DB->v)||(nu<DB->nu)) {rest=DB->NFnum[DB->v][DB->nu]; DB->p=0;}
#ifdef	TEST
if((v<DB->v)||(nu<DB->nu)){assert(totNF==nuNF);nuNF+=DB->NFnum[DB->v][DB->nu];}
	if(v<DB->v) {
if(vNF){printf("vNF[%d]=%d was=%d v_new=%d\n",v,xvNF,vNF,DB->v);exit(0);}
		for(i=1;i<=DB->NUCmax;i++) vNF+=DB->NFnum[DB->v][i];}
#endif

     }	assert(DB->p<=DB->NFnum[DB->v][DB->nu]); assert(rest);
     assert(totNF < DB->nNF);
#ifdef	TEST
	printf("return 1 at  totNF=%d  v=%d  nu=%d  p=%d\n",
	    totNF,DB->v,DB->nu,DB->p); assert(0<vNF--);
#endif
     AuxGet_uc(DB->Fv[DB->v],&DB->nu,uc); ++DB->p; totNF++; return 1;
}

int  Read_SLucNF_from_DB(void)
{    puts("Read_SLucNF_from_DB: to be implemented"); exit(0); return 0;
}

int  Read_SLpoly_from_DB(void)
{    puts("Read_SLpoly_from_DB: to be implemented"); exit(0); return 0;
}

int  Read_H_poly_from_DB(DataBase *DB,PolyPointList *P)
{    static unsigned char uc[NUC_Nmax]; static int ms3; 
     int i,j, MS; Long NF[POLY_Dmax][VERT_Nmax]; VertexNumList V; EqList E; 
     if(0==ms3) { if(0==Read_H_ucNF_from_DB(DB,uc)) return 0;} else --uc[0];
     UCnf2vNF(&DB->d,&DB->v,&DB->nu,uc,NF,&MS); P->n=DB->d; P->np=DB->v;
     for(i=0;i<P->np;i++) for(j=0;j<P->n;j++) P->x[i][j]=NF[j][i];
     MS%=4;
     ms3=(MS==3); if(MS==2)
     {	assert(Ref_Check(P,&V,&E)); P->np=E.ne;
	for(i=0;i<P->np;i++)for(j=0;j<P->n;j++) P->x[i][j]=E.e[i].a[j];
     }	return 1;
}

int  Read_H_poly_from_DB_or_inFILE(DataBase *DB,PolyPointList *P)
{    if((DB==NULL)||(inFILE==NULL)) {CWS cws; return Read_CWS_PP(&cws,P);}
     else return Read_H_poly_from_DB(DB,P);
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
Long AuxGxP(Long *Gi,Long *V,int *d)
{    Long x=0; int j; for(j=0;j<*d;j++) x+=Gi[j]*V[j]; return x;
}

/*   Glz x (lincomb(Points)) -> Diag: if(index>1) print vertices of G*P & D */
#define TEST
#define TEST_OUT
void Aux_Print_CoverPoly(int *I,int *d, int *N,Long *X[POLY_Dmax],
	Long G[][POLY_Dmax],Long *D,int *x,Long Z[][VERT_Nmax],Long *M,int r)
{    int i,j,dia=1, err=0; fprintf(outFILE,"%d %d    index=%d  D=%ld",*d,*N,
	*I,D[0]); for(i=1;i<*d;i++) printf(" %ld",D[i]); 
     for(i=0;i<r;i++)
     {  fprintf(outFILE," /Z%ld:",M[i]);
        for(j=0;j<*N;j++) fprintf(outFILE," %ld",Z[i][j]);
     }	printf("  #%d\n",*x);
     for(i=0;i<*d;i++){
       for(j=0;j<*N;j++){
	 Long Xij=AuxGxP(G[i],X[j],d);
	 if(0!=(Xij % D[i])) err=1;
	 printf("%ld ",Xij);}
       puts("");
     }
#ifdef	TEST
	if((dia==0)||err){
	  int i,j;printf("D=");for(i=0;i<*d;i++)printf("%ld ",D[i]);
	  printf("   index[%d]=%d\n",*x,*I);
	  for(i=0;i<*d;i++)	    {
	    printf("G=");
	    for(j=0;j<*d;j++) printf("%2ld ",G[i][j]);
	    printf("   G.P=");
	    for(j=0;j<*N;j++) printf("%2ld ",AuxGxP(G[i],X[j],d));
	    puts("");}
	  exit(0);}
#endif
}
void Aux_Print_SLpoly(int *I,int *d, int *N,Long *X[POLY_Dmax],
	Long G[][POLY_Dmax],Long *D,int *x)
{    int i,j,dia=1, err=0; printf("%d %d    index=%d  D=%ld",*d,*N,*I,D[0]);
     for(i=1;i<*d;i++) printf(" %ld",D[i]);
     printf("  #%d\n",*x);
     for(i=0;i<*d;i++)
     {	/* int z=1;*/ for(j=0;j<*N;j++)
	{    Long Xij=AuxGxP(G[i],X[j],d);
#ifdef	TEST_STUFF
	if(Xij){if(z&&dia)dia=(labs(Xij)==D[i]);z=0;} /* SUFFICIENT only !!!*/
        if(dia==0)if(i==(*d-1)){Long g=Xij,T; int J; for(J=j+1;J<*N;J++){T=
	    labs(AuxGxP(G[i],X[J],d)); if(T) g=Fgcd(g,T);} if(g==D[i]) dia=1;}
#endif
	    if(0!=(Xij % D[i])) err=1;
	    printf("%ld ",Xij/D[i]);}
       puts("");
     }
#ifdef	TEST
	if((dia==0)||err){
	    int i,j;printf("D=");for(i=0;i<*d;i++)printf("%ld ",D[i]);
	    printf("   index[%d]=%d\n",*x,*I);for(i=0;i<*d;i++)
	    {printf("G=");for(j=0;j<*d;j++)printf("%2ld ",G[i][j]);
	    printf("   G.P=");
	    for(j=0;j<*N;j++) printf("%2ld ",AuxGxP(G[i],X[j],d));
	    puts("");}
	    exit(0);}
#endif
}
void Aux_Make_Dual(PolyPointList *P, VertexNumList *V, EqList *E)
{    Long VM[VERT_Nmax][POLY_Dmax]; int i,j, d=P->n, e=E->ne, v=V->nv; 
     assert(e<=VERT_Nmax); P->np=V->nv=e; E->ne=v;
     for(i=0;i<v;i++)for(j=0;j<d;j++) VM[i][j]=P->x[V->v[i]][j]; 
     for(i=0;i<e;i++){for(j=0;j<d;j++)P->x[i][j]=E->e[i].a[j]; V->v[i]=i;}
     for(i=0;i<v;i++){for(j=0;j<d;j++)E->e[i].a[j]=VM[i][j]; E->e[i].c=1;}
     assert(Ref_Check(P,V,E));
}
void PrintVPHMusage(void);
int  Make_Lattice_Basis(int d, int p, Long *P[POLY_Dmax],  /* index=det(D) */
        Long G[][POLY_Dmax], Long *D);/* G x P generates diagonal lattice D */
void PH_Sublat_Polys(char *dbin, int omitFIP, PolyPointList *_P, char sF) 
{     static EqList E; VertexNumList V; int x=0, K, B, I=1; /* index>I only */
     Long *RelPts[POINT_Nmax];DataBase *DB=NULL; if(*dbin)Open_DB(dbin,&DB,0);
     K=((sF!='P')&&(sF!='H')&&(sF!='Q')&&(sF!='B'));    /*  K=='CoverPoly'  */
     if(('1'<sF)&&(sF<='9')) I=sF-'0';
     B=(omitFIP==2);	/* 'q' for index >I */
     while(Read_H_poly_from_DB_or_inFILE(DB,_P))
     {	Long D[POLY_Dmax],G[POLY_Dmax][POLY_Dmax],PM[VERT_Nmax][VERT_Nmax];
 	int index, N=0; /* if(!Ref_Check(_P,&V,&E)) Print_PPL(_P,""); */
	if(B)assert(_P->n==4);/* b:Brower group only for CY hypersurface d=4 */
	assert(Ref_Check(_P,&V,&E));  
        /* Aux_Make_Dual(_P,&V,&E); */ /* don't dualize: take M-lattice poly */
	Make_VEPM(_P,&V,&E,PM); _P->np=V.nv; Complete_Poly(PM,&E,V.nv,_P); ++x;
	if(omitFIP==0) for(N=0;N<_P->np;N++) RelPts[N]=_P->x[N];
	else if(omitFIP<3 )		  		  /* Omit facet-IPs */
	{   int p; for(p=0;p<_P->np;p++)	        
	    {	int e, z=0;
	      for(e=0;e<E.ne;e++) 
	        if(0==Eval_Eq_on_V(&E.e[e],_P->x[p],_P->n)) z++;
	      if(z>omitFIP) RelPts[N++]=_P->x[p];
	    }	assert(V.nv<=N);   /* count <E,.>=0; if(n>1) add_to_RelPts; */
	}
	else if(omitFIP==3)			   /* Omit all non-vertices */
        {    for(N=0;N<V.nv;N++) RelPts[N]=_P->x[N];
	/*   { int tnv=V.nv; Find_Equations(_P,&V,&E);assert(V.nv==tnv);
	       for(tnv=0;tnv<V.nv;tnv++)assert(V.v[tnv]<V.nv);} */
 	}
        else {puts("something wrong in PH_Sublat_Polys"); exit(0);}
	if(K)
	{   Long Z[POLY_Dmax][VERT_Nmax], M[POLY_Dmax]; int r;
	    index=Sublattice_Basis(_P->n,N,RelPts,Z,M,&r,G,D);
	    if(B)if(D[2]==1)continue;
	    if(omitFIP==3)if(D[1]==1)continue;
	    assert(index>0); if(index<=I)continue; assert(Ref_Check(_P,&V,&E));
	    Aux_Print_CoverPoly(&index,&_P->n,&N,RelPts,G,D,&x,Z,M,r);
	}
	else 
	{   index=Make_Lattice_Basis(_P->n,N,RelPts,G,D); if(1==index)continue;
	    if(B)if(D[2]==1)continue;
	    if(omitFIP==3)if(D[1]==1)continue;
	    assert(index>0); assert(Ref_Check(_P,&V,&E)); Print_VL(_P,&V,"");
            Aux_Print_SLpoly(&index,&_P->n,&N,RelPts,G,D,&x);
	}
     }	if(*dbin) Close_DB(DB);	
}
#undef	TEST
/*	Lattice generated by vertices; UT-decomp of diag	*/
void V_Sublat_Polys(char mr,char *dbin,char *polyi,char *polyo, 
	PolyPointList *_P)
{    NF_List *_L=(NF_List *) malloc(sizeof(NF_List)); 
     int max_order=1;
     static EqList E; VertexNumList V; int x=0;
     Long *RelPts[VERT_Nmax];DataBase *DB=NULL; if(*dbin)Open_DB(dbin,&DB,0);
     assert(_L!=NULL); 
     if(!(*polyo)) {
	puts("You have to specify an output file via -po in -sv-mode."); 
    	printf("For more help use option `-h'\n");
	exit(0);}
     _L->of=0; _L->rf=0; _L->iname=polyi; _L->oname=polyo; _L->dbname=dbin;
     Init_NF_List(_L); _L->SL=0;
     while(Read_H_poly_from_DB_or_inFILE(DB,_P))
     {	Long D[POLY_Dmax],G[POLY_Dmax][POLY_Dmax];
 	int index, N; assert(Ref_Check(_P,&V,&E)); 
        for(N=0;N<V.nv;N++) RelPts[N]=_P->x[V.v[N]];
	++x;
	index=Make_Lattice_Basis(_P->n,N,RelPts,G,D); if(1==index) continue;
	assert(index>0); if(index>max_order) max_order=index;
#ifdef	TEST
	Aux_Print_SLpoly(&index,&_P->n, &N,RelPts,G,D,&x);
#endif
	{   int i,j; subl_int diag[POLY_Dmax], U[POLY_Dmax][VERT_Nmax];
	    for(i=0;i<_P->n;i++){diag[i]=D[i]; for(j=0;j<V.nv;j++) {int k;
		U[i][j]=0; for(k=0;k<_P->n;k++) U[i][j]+=G[i][k]*RelPts[j][k];
		assert(0==(U[i][j]%D[i])); U[i][j]/=D[i];}}
	    Make_All_Sublat(_L, _P->n, V.nv, diag, U, &mr, _P);
	}
     }	if(*dbin) Close_DB(DB);
     printf("max_order=%d\n", max_order); Write_List_2_File(polyo,_L); 
     _L->TIME=time(NULL); fputs(ctime(&_L->TIME),stdout);
     free(_L);
}
void VPHM_Sublat_Polys(char sFlag,char mr,char *dbin,char *polyi,char *polyo, 
		       PolyPointList *_P)
{    switch(sFlag){	/* if(dbin=0) read from inFILE; */
	case 'p': case 'P': PH_Sublat_Polys(dbin,0,_P,sFlag); 	       break ;
	case 'h': case 'H': PH_Sublat_Polys(dbin,1,_P,sFlag); 	       break ;
	case 'b': case 'B': PH_Sublat_Polys(dbin,2,_P,sFlag); 	       break ;
	case 'q': case 'Q': PH_Sublat_Polys(dbin,3,_P,sFlag); 	       break ;
	case 'v': case 'V': V_Sublat_Polys(mr,dbin,polyi,polyo,_P);    break ;
	case 'm': case 'M': Find_Sublat_Polys(mr,dbin,polyi,polyo,_P); break ;
        default:if(('1'<sFlag)&&(sFlag<='9'))PH_Sublat_Polys(dbin,3,_P,sFlag);
        else{puts("-s# requires that # is in {v,p,h,b,m,q}");PrintVPHMusage();
     } }
}
void PrintVPHMusage(void)
{puts("	-sh ... gen by codim>1 points (omit IPs of facets)");
 puts("	-sp ... gen by all points");
 puts("	-sb ... generated by dim<=1 (edges), print if rank=2	");
 puts("	-sq ... generated by vertices,       print if rank=3	");
 puts("	    q,b currently assume that dim=4");
 exit(0);
}



void Bin_2_ANF(char *polyi, int max, PolyPointList *_P)
{    FILE *F=fopen(polyi,"rb"); FInfoList L; UPint list_num,tNF=0; Along tNB=0;
     int d, v, s, sl_nNF, sl_SM, sl_NM, sl_NB, mc=0, MS, nu; unsigned i, j; 
     unsigned char uc[POLY_Dmax*VERT_Nmax]; VertexNumList V; EqList E;  
     Long NF[POLY_Dmax][VERT_Nmax];                     Init_FInfoList(&L);

     if(F==NULL) {printf("Input file %s not found\n",polyi); exit(0);}
     d=fgetc(F); assert(d==0);  /* for(i=0;i<d;i++) fgetc(F); */
     d=fgetc(F); L.nV=fgetc(F); L.nVmax=fgetc(F); L.NUCmax=fgetc(F);
     list_num=fgetUI(F); 
     L.nNF=fgetUI(F); L.nSM=fgetUI(F); L.nNM=fgetUI(F);   L.NB=fgetUI(F); 
     sl_nNF=fgetUI(F); sl_SM=fgetUI(F); sl_NM=fgetUI(F);     sl_NB=fgetUI(F); 

     for(i=0;i<L.nV;i++)
     {  v=fgetc(F); L.nNUC[v]=fgetc(F);   /* read #nuc's per #Vert */
        for(j=0;j<L.nNUC[v];j++)                  
        {   L.NFnum[v][nu=fgetc(F)]=fgetUI(F);  /* read nuc and #NF(v,nu)*/
            tNF+=L.NFnum[v][nu]; tNB+=L.NFnum[v][nu]*(Along)nu;
        }
     }  assert( 0 == (unsigned int)(tNB-L.NB) );  L.NB=tNB; assert(tNF==L.nNF);

     for(v=d+1;v<=L.nVmax;v++) if(L.nNUC[v])         /* write  honest polys */
     {  int I,J; for(nu=1;nu<=L.NUCmax;nu++) for(j=0;j<L.NFnum[v][nu];j++)
        {   for(s=0;s<nu;s++) uc[s]=fgetc(F);
            UCnf_2_ANF(&d,&v,&nu,uc,NF,&MS); MS %= 4; _P->n=d; _P->np=v;
            for(I=0;I<v;I++) for(J=0;J<d;J++) _P->x[I][J]=NF[J][I];
            /* assert(Ref_Check(_P,&V,&E)); */ assert(MS==1);

            if(MS!=2)                                 /* if(MS!=2) print NF */
               if(!max || Poly_Max_check(_P,&V,&E)) 
                {   mc++; Print_NF(outFILE,&d,&v,NF);}
            if(MS>1)                      /* if(MS>1); print Mirror */
              if(!max || Poly_Min_check(_P,&V,&E))
                {   mc++; Small_Make_Dual(_P,&V,&E); 
                    Make_Poly_NF(_P,&V,&E,NF); Print_NF(outFILE,&d,&(V.nv),NF);
                }
        }
     }
     printf("np=%lld+%dsl  ",2*L.nNF-L.nSM-L.nNM,2*sl_nNF-sl_SM-sl_NM);
     printf(                                          /* write Finfo */
        "%dd  %dv<=%d n<=%d  %dnv  %lld %d %lld %lld  %d %d %d %d",
        d, L.nV, L.nVmax, L.NUCmax,list_num,
        L.nNF,L.nSM,L.nNM,L.NB,  sl_nNF, sl_SM, sl_NM, sl_NB);
     if(max) printf("  r-max=%d",mc);
     puts("");
}

void Bin_2_ANF_DBsl(char *dbi, int max, int vf, int vt, PolyPointList *_P)
{    FILE *F; FInfoList L;
     int d, v, nu, i, j, list_num, mc=0, MS,
        sl_nNF, sl_SM, sl_NM, sl_NB, tSM=0, tNM=0;
     char *Ifx, *Ifn = (char *) malloc(1+strlen(dbi)+File_Ext_NCmax);
     Long NF[POLY_Dmax][VERT_Nmax]; 
     VertexNumList V;EqList E;
     strcpy(Ifn,dbi); Ifx=&Ifn[strlen(dbi)]; strcpy(Ifx,".info");
     F=fopen(Ifn,"r"); if(F==NULL) {puts("Info File not found");exit(0);}
     Init_FInfoList(&L);                       /* start reading the file */
     fscanf(F,"%d%d%d%d%d%lld%d%lld %lld %d%d%d%d",
            &d,&i,&j,&nu,&list_num,&L.nNF,&L.nSM,&L.nNM,&L.NB,
            &sl_nNF,&sl_SM,&sl_NM,&sl_NB);   L.nV=i; L.nVmax=j; L.NUCmax=nu;
     if(sl_NB && (vf==2)&&(vt==VERT_Nmax-1))
     {  strcpy(Ifx,".sl"); fclose(F);
        if(NULL==(F=fopen(Ifn,"rb"))){printf("Open %s failed",Ifn);exit(0);}
     }  
     else /* puts("no .sl file"); */ 
     {  DataBase *DB; Open_DB(dbi, &DB, 0); if((DB->v<vf)||(DB->nVmax>vt))
            DB_fromVF_toVT(DB,vf,vt); /*  read only vf <= v <= vt :: */
        for(i=0; Read_H_poly_from_DB(DB,_P); i++) 
	{   int c, l, p=_P->np-1, off=_P->x[p][0]; 
	    if(off) for(l=0;l<d;l++) for(c=l;c<=p;c++) _P->x[c][l]-=off;
					/* if(off) Print_PPL(_P,"off"); else */
	    Print_PPL(_P,"");
	}   printf("#poly=%d\n",i); Close_DB(DB);
     }
     for(i=0;i<sl_nNF;i++)
     {  int I,J; unsigned char uc[NUC_Nmax]; v=fgetc(F); assert(v<=VERT_Nmax); 
        nu=fgetc(F); AuxGet_uc(F,&nu,uc); assert(!ferror(F));
        if((*uc%4)==0)tSM++; else if((*uc%4)<3)tNM++; 
        
        UCnf_2_ANF(&d,&v,&nu,uc,NF,&MS); MS %= 4; _P->n=d; _P->np=v;
        for(I=0;I<v;I++) for(J=0;J<d;J++) _P->x[I][J]=NF[J][I];
        if(MS!=2)                                 /* if(MS!=2) print NF */
          if(!max || Poly_Max_check(_P,&V,&E)) 
            {   mc++; Print_NF(outFILE,&d,&v,NF);}
        if(MS>1)                              /* if(MS>1); print Mirror */
          if(!max || Poly_Min_check(_P,&V,&E))
            {   mc++; Small_Make_Dual(_P,&V,&E); 
                Make_Poly_NF(_P,&V,&E,NF); Print_NF(outFILE,&d,&(V.nv),NF);
            }
     }  
     printf("np=%lld+%dsl  ",2*L.nNF-L.nSM-L.nNM,2*sl_nNF-sl_SM-sl_NM);
     printf(                                          /* write Finfo */
        "%dd  %dv<=%d n<=%d  %dnv ... %d %d %d %d",  /* %d %d %d %lld */ 
        d, L.nV, L.nVmax, L.NUCmax,list_num,
        /* L.nNF,L.nSM,L.nNM,L.NB, */ sl_nNF, sl_SM, sl_NM, sl_NB);
     if(max) printf("  r-max=%d",mc);
     puts("");
}

void Gen_Bin_2_ascii(char *pi,char *dbi,int max,int vf,int vt,PolyPointList *P)
{    if(*pi)     Bin_2_ANF(pi, max, P);
     else if(*dbi) Bin_2_ANF_DBsl(dbi, max, vf, vt, P);
     else puts("With -B[2A] you have to specify input via -pi or -di"); 
}

