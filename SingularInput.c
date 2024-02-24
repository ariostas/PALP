/* ======================================================== */
/* ===                                                  === */
/* ===            S i n g u l a r I n p u t . c         === */
/* ===                                                  === */
/* ===	Authors: Maximilian Kreuzer, Nils-Ole Walliser	=== */
/* ===	Emanuel Scheidegger                             === */
/* ===	Last update: 05/04/12                           === */
/* ===                                                  === */
/* ======================================================== */


/* ======================================================== */
/* =========            H E A D E R s             ========= */

#include <unistd.h> /* close */

#include "Global.h"
#include "Mori.h"

/* ======================================================== */
/* =========            D E F I N I T I O N s     ========= */

/***  local for Singularinput.c ***/
#define DijkEQ                  "->"	/* Useful for Mathematica rules: "->" */
#define	T_DIV                   "d"	    /* Toric divisors */
#define DIVclassBase            "J"	    /* Basis of the divisor classes */

/* diagnostic stuff */
#define TEST_PRINT_SINGULAR_IO  (0)     /* IO of SINGULAR */

#define NORM_SIMP_NUM           (0)

/*=========================================================*/

#if (TEST_PRINT_SINGULAR_IO)
void CatFile(char *fn){
  char* CAT = (char*)malloc(30 + strlen(fn));
  strcpy(CAT,"cat ");
  strcat(CAT,fn);
  printf("======= FILE content of %s:\n",fn);
  fflush(0);
  assert(0==system(CAT));
  printf("====== End of FILE content of %s\n\n",fn);
  fflush(0);
  free(CAT);
}
#endif

int Read_HyperSurf(int *he, int divclassnr, int maxline, char filename[20], MORI_Flags *_Flag){

    FILE *stream;
    
    int i;
    char string[maxline];
    char delims[] = " ";
    char *result = NULL;
	
	if(_Flag->Read_HyperSurfCounter==0){	
    		if( (stream = fopen(filename, "w")) == NULL) {
        		printf("Error: cannot open file!\n");
            		exit(1);
        	}     
 
		fgets(string, sizeof string, stdin);
		fprintf(stream, "%s\n", string);
	}
	
	if(_Flag->Read_HyperSurfCounter != 0){
		if( (stream = fopen(filename,"r")) == NULL) {
			printf("Error: cannot read file!\n");
			exit(1);	
		}
		fgets(string, maxline ,stream);
	}	

        i=0;
        result = strtok( string, delims );
        while( result != NULL ) { 
            he[i] = atoi(result); 
            i++;
            if(i == divclassnr)
                break;
            result = strtok( NULL, delims );
        }   		
        fclose(stream);

        return i;
}

void HyperSurfSingular(PolyPointList *P,triang *T, triang *SR ,MORI_Flags *_Flag , FibW *F,int *cp){
  int p=SR->v,d=SR->d,i,j,r, N=0;/* N=not(0th)=offset(IP) */
  int divclassnr = *cp;
  int TORDIM=P->n;
  int CODIM=1;
  int DIM=TORDIM-CODIM;
  char *D=T_DIV,*B=DIVclassBase;

  /* Put temporary files in $TMPDIR if it is set */
  char* tmpdir = getenv("TMPDIR");
  if (tmpdir == NULL) {
    tmpdir = "/tmp";
  }

  /* Add one to ensure room for the null byte at the end */
  size_t template_size = strlen(tmpdir) + strlen("/SFnameXXXXXX") + 1;
  char* SFname = (char*)malloc(template_size);
  snprintf(SFname, template_size, "%s/SFnameXXXXXX", tmpdir);

  int SF = mkstemp(SFname);
  assert(-1 != SF);

  /* Construct the singular command. In this case, template_size is
     already padded by 1 byte (see above). */
  size_t singular_cmd_size = strlen("Singular -q < ") + template_size;
  char* SingularCall = (char*)malloc(singular_cmd_size);
  snprintf(SingularCall, singular_cmd_size, "Singular -q < %s", SFname);

  dprintf(SF,"LIB \"general.lib\";\n");
  dprintf(SF,"option(noredefine);\n");
  dprintf(SF,"ring r=0,(t,%s%d",D,N+1);	// ring r=0,(t,D1,...,Dp,B1,...,Bp-d),ds;
  for(i=N+2;i<=p;i++)
	  dprintf(SF,",%s%d",D,i);
  for(i=1;i<=p-d;i++)
	  dprintf(SF,",%s%d",B,i);
  dprintf(SF,"),(ds,dp(%d));\n",2*p-d);

  /*LINEAR EQUIVALENCES*/
  dprintf(SF,"ideal lin=");
  for(r=0;r<d;r++){
    if(r)dprintf(SF,",");
    dprintf(SF,"0");
    for(j=0;j<p;j++){
    	int X=P->x[j][r];
    	if((X>0))
    		dprintf(SF,"+");
    	if(X)
    		dprintf(SF,"%d*%s%d",X,D,j+1);
    }
  }
  dprintf(SF,";");
  DivClassBasis(SF,P,p,D,B);

  /* STANLEY REISNER */
  dprintf(SF,"ideal sr=");
  for(r=0;r<SR->n;r++){
	  int f=1;
	  if(r)
		  dprintf(SF,",");
	  for(j=0;j<p;j++)
		  if(getN(j,SR->I[r])){
			  if(f)
				  f=0;
			  else
				  dprintf(SF,"*");
			  dprintf(SF,"%s%d",D,j+1);
		  }
  }

  /* CHOW RING */
  dprintf(SF,";\nideal chow=std(lin+sr+DCbase);\n");

    /*reads input for the hyper surface
     * creates HEInput.txt
     * reads HEInput.txt
     * and defines poly cy as H*/
  if(_Flag->H==1){

	  if(!(_Flag->FilterFlag) && _Flag->Read_HyperSurfCounter == 0)
		  printf("Type the %d (integer) entries for the hypersurface class:\n", divclassnr);


 	int *he, i;
        he=malloc(divclassnr*sizeof(int));         
	
	for(i=0;i<divclassnr;i++)        
		he[i]=0;

	int control = Read_HyperSurf(he, divclassnr, 5*divclassnr ,"HEInput.txt", _Flag);
	
if(_Flag->Read_HyperSurfCounter==0){
         if(control != divclassnr)
                 printf("Warning: %d entries for the hypersurface class were expected\n"
                 "         (i.e. as many entries as toric divisor classes) \n",divclassnr );

dprintf(SF,"\"Hypersurface class: ");
	  for(i=0; i<divclassnr; i++)
		  if(he[i]!=0)
			  dprintf(SF,"%d*d%d ",he[i],i+1);
	  dprintf(SF,"\";\n");
}
	  dprintf(SF,"poly HySurf=(0");
	  for(i=0;i<divclassnr;i++)		 // Hypersurface equation
		  if(he[i]!=0)
			  dprintf(SF,"+%d*d%d",he[i],i+1);
	  dprintf(SF,"); \n");

	  int dim = F->nw;
	  int DegreeVec[dim];
	  for(i=0;i<F->nw;i++){
		  DegreeVec[i]=0;
	  }

	  for(i=0;i<F->nw;i++){
    		  for(j=0;j<divclassnr;j++){
    			  DegreeVec[i] += he[j]*(F->W[i][j]);
    		  }
	  }
if(_Flag->Read_HyperSurfCounter==0){
	  fprintf(outFILE,"Hypersurface degrees: (");
	  for(i=0;i<F->nw;i++){
		  fprintf(outFILE," %d ",DegreeVec[i]);
	  }
	  fprintf(outFILE,")\n");
}

	
//	_Flag->Read_HyperSurfCounter++;
	
	  free(he); fflush(0);

  }

  else {
	  dprintf(SF,"poly HySurf=(%s%d",D,N+1);
	  for(i=N+1;i<p;i++)		 // CALABI-YAU
		  dprintf(SF,"+%s%d",D,i+1);
	  dprintf(SF,");\n");
  }

  dprintf(SF,"number vol=");{
	  int NN=NORM_SIMP_NUM;          // NORMALIZATION
	  Long *X[POLY_Dmax], Vijkl;
	  int x[POLY_Dmax], n=0; i=j=0;
	  if((NN<0)||(NN > T->n)) NN=0;
//  for(j=1;j<T->n;j++)if(T->I[j]<T->I[NN])NN=j; //avoid D_I.DHySurf=0 for P3xP3/Z3
	  while(i<p){
    	if(getN(i++,T->I[ NN ])){
    		x[n]=i-1;
    		X[n++]=P->x[i-1];
    	}
	  }
    assert(n==P->n);
    Vijkl=SimplexVolume(X,P->n);
    dprintf(SF,"%ld;\n",Vijkl);
    dprintf(SF,"poly norm=vol*reduce(%s%d",D,x[0]+1);
	for(j=1;j<TORDIM;j++) dprintf(SF,"*%s%d",D,x[j]+1);
	dprintf(SF,",chow);\n");
  }

  if(_Flag->i || _Flag->c){
  dprintf(SF,"\"SINGULAR -> divisor classes (integral basis %s1",B);
  if(p-d>1)
	  dprintf(SF," ... %s%d",B,p-d);
  dprintf(SF,"):\";\n");
  dprintf(SF,"string LR=string(%s1)+\"=\"+string(reduce(%s1,chow));\n",D,D);
  for(j=2;j<=p;j++){
	  dprintf(SF,"LR=LR+\", \"+string(%s%d)+\"=\"+",D,j);
	  dprintf(SF,"string(reduce(%s%d,chow));\n ",D,j);
  }
  dprintf(SF,"LR;\n");
  }

	dprintf(SF, "list dli=d1");
	for(j=2;j<=p;j++)
		dprintf(SF,",d%d",j);
	dprintf(SF,";\n");
	dprintf(SF, "list Jli=J1");
	for(j=2;j<=p-d;j++)
		dprintf(SF,",J%d",j);
	dprintf(SF,";\n");
	
  dprintf(SF,"int i,j,k;\n");

  /*linear basis ideal*/
  dprintf(SF,"ideal lb=std(lin+DCbase);\n");

	/* Test for nonintersecting divisors */
	dprintf(SF,"ideal hypideal = quotient(chow,HySurf);\n");
	if(_Flag->t || _Flag->d){
		dprintf(SF,"list nonintersectingd = list();\n");
		dprintf(SF,"for(i=1;i<=%d;i++){\n",p);
		dprintf(SF,"  if(reduce(dli[i],std(hypideal))==0){\n");
		dprintf(SF,"    nonintersectingd=nonintersectingd+list(dli[i]);\n");
	    dprintf(SF,"  }\n");
		dprintf(SF,"}\n");
		dprintf(SF,"if(size(nonintersectingd)>0){\n");
		dprintf(SF,"  \"SINGULAR -> nonintersecting divisor classes : \";\n");
		dprintf(SF,"  string nd=string(nonintersectingd[1]);\n");
		dprintf(SF,"  for(i=2;i<=size(nonintersectingd);i++){\n");
		dprintf(SF,"    nd=nd+\", \"+string(nonintersectingd[i]);\n");
		dprintf(SF,"  }\n");
		dprintf(SF,"nd;}\n");
	}
	if(_Flag->i || _Flag->c){
		dprintf(SF,"list nonintersectingJ = list();\n");
		dprintf(SF,"for(i=1;i<=%d;i++){\n",p-d);
		dprintf(SF,"if(reduce(Jli[i],std(hypideal))==0){nonintersectingJ=nonintersectingJ+list(Jli[i]);}\n");
		dprintf(SF,"}\n");
		dprintf(SF,"if(size(nonintersectingJ)>0){\n");
		dprintf(SF,"  \"SINGULAR -> nonintersecting divisor classes (integral basis): \";\n");
		dprintf(SF,"  string nJ=string(nonintersectingJ[1]);\n");
		dprintf(SF,"  for(i=2;i<=size(nonintersectingJ);i++){\n");
		dprintf(SF,"    nJ=nJ+\", \"+string(nonintersectingJ[i]);\n");
		dprintf(SF,"  }\n");
		dprintf(SF,"nJ;}\n");
	}

	
/* Singular procedure to compute the Todd class of a vector bundle given by its Chern character */		
	dprintf(SF,"proc Todd(poly ch, int Dim) {\n");
	dprintf(SF,"  poly ToddD = 1;\n");
	dprintf(SF,"  for(i=2;i<=Dim+1;i++){ToddD=ToddD+(-dli[1]*t)^(i-1)/factorial(i);}\n");
	dprintf(SF,"  ToddD = jet(1,ToddD,2*Dim);\n");
	dprintf(SF,"  matrix c=coeffs(ToddD,t);\n");
	dprintf(SF,"  matrix p[Dim][1] = c[2..Dim+1,1];\n");
	dprintf(SF,"  for(i=1;i<=Dim;i++){\n");
	dprintf(SF,"    p[i,1]=-i*p[i,1];\n");
	dprintf(SF,"    for(j=1;j<=i-1;j++){\n");
	dprintf(SF,"      p[i,1]=p[i,1]-p[j,1]*c[i-j+1,1];\n");
	dprintf(SF,"    };\n");
	dprintf(SF,"  };\n");
	dprintf(SF,"  poly tmp=0;\n");
	dprintf(SF,"  for(i=1;i<=Dim;i++){\n");
	dprintf(SF,"    tmp = tmp+p[i,1]/factorial(i)*(-t)^i;\n");
	dprintf(SF,"  };\n");
	dprintf(SF,"  matrix ctmp=coeffs(tmp,dli[1]);\n");
	dprintf(SF,"  matrix cch=coeffs(ch,t);\n");
	dprintf(SF,"  int m=nrows(ctmp);\n");
	dprintf(SF,"  if(nrows(cch) < m){m=nrows(cch)};\n");
	dprintf(SF,"  tmp=0;\n");
	dprintf(SF,"  for(i=0;i<=m-1;i++){\n");
	dprintf(SF,"    tmp=tmp+ctmp[i+1,1]*cch[i+1,1]*factorial(i);\n");
	dprintf(SF,"  }\n");
	dprintf(SF,"  p=coeffs(tmp,t);\n");
	dprintf(SF,"  for(i=1;i<=nrows(p)-1;i++){\n");
	dprintf(SF,"    p[i+1,1]=p[i+1,1]*(-1)^(i)*factorial(i);\n");
	dprintf(SF,"  }\n");
	dprintf(SF,"  c[1,1]=1;\n");
	dprintf(SF,"  for(i=1;i<=Dim;i++){\n");
	dprintf(SF,"    c[i+1,1]=0;\n");
	dprintf(SF,"    m=i;\n");
	dprintf(SF,"    if(nrows(p) < i+1){m=nrows(p)-1};\n");
	dprintf(SF,"    for(j=1;j<=m;j++){\n");
	dprintf(SF,"      c[i+1,1]=c[i+1,1]-p[j+1,1]*c[i+1-j,1];\n");
	dprintf(SF,"    };\n");
	dprintf(SF,"    c[i+1,1]=c[i+1,1]/(i);\n");
	dprintf(SF,"  };\n");
	dprintf(SF,"  tmp=0;\n");
	dprintf(SF,"  for(i=0;i<=Dim;i++){\n");
	dprintf(SF,"    tmp=tmp+c[i+1,1]*t^i;\n");
	dprintf(SF,"  }\n");	  
	dprintf(SF,"  return(tmp);\n");
	dprintf(SF,"};\n");
		

/*####### Chern classes of CY/H ################*/
	  dprintf(SF,"poly ChernUp=1;\n");
	  dprintf(SF,"for(i=1;i<=%d;i++){ChernUp=ChernUp*(1+dli[i]);}\n",p);
  	  dprintf(SF,"poly ChernBottomHypersurface=1;\n");
	  dprintf(SF,"for(i=1;i<=%d;i++){ChernBottomHypersurface=ChernBottomHypersurface+(-1)^i*HySurf^i;}\n",DIM);
	  dprintf(SF,"poly ChernHySurf=ChernUp*ChernBottomHypersurface;\n");
	  for (j=0;j<DIM; j++) dprintf(SF,"poly c%dHySurf=jet(ChernHySurf,%d)-jet(ChernHySurf,%d);\n",j+1,j+1,j);

	
	  // PRINT HODGE and EULER
	  if(_Flag->b){
		  if(_Flag->H !=1){dprintf(SF,"\"SINGULAR  -> Arithmetic genera and Euler number of the CY:\";\n");}
		  else if (_Flag->H ==1){dprintf(SF,"\"SINGULAR  -> Arithmetic genera and Euler number of H:\";\n");}
		  /* Chern character and the Todd class of the tangent bundle */
		  dprintf(SF,"matrix ch[%d][1] = 1",DIM+1);
		  for(i=2;i<=DIM+1;i++){
			  dprintf(SF,",c%dHySurf",i-1);
		  };
		  dprintf(SF,";\n");
		  dprintf(SF,"matrix pp[%d][1] = -ch[2,1]",DIM);
		  for(i=2;i<=DIM;i++){
			  dprintf(SF,",-%d*ch[%d,1]",i,i+1);
		  };
		  dprintf(SF,";\n");
		  dprintf(SF,"for(i=1;i<=%d;i++){\n",DIM);
		  dprintf(SF,"  for(j=1;j<=i-1;j++){\n");
		  dprintf(SF,"    pp[i,1]=pp[i,1]-pp[j,1]*ch[i-j+1,1];\n");
		  dprintf(SF,"  };\n");
		  dprintf(SF,"};\n");
		  dprintf(SF,"poly tmp=%d;\n",DIM);
		  dprintf(SF,"for(i=1;i<=%d;i++){\n",DIM);
		  dprintf(SF,"  tmp = tmp+pp[i,1]/factorial(i)*(-t)^i;\n");
		  dprintf(SF,"};\n");
		  dprintf(SF,"poly ToddHySurf=subst(Todd(tmp,%d),t,1);\n",DIM);
		  dprintf(SF,"tmp=%d;\n",DIM);
		  dprintf(SF,"for(i=1;i<=%d;i++){\n",DIM);
		  dprintf(SF,"  tmp = tmp+pp[i,1]/factorial(i)*(t)^i;\n");
		  dprintf(SF,"};\n");
		  dprintf(SF,"tmp=subst(tmp,t,1);\n");
		  
		  
		  dprintf(SF,"\"chi_0: \",");
		  dprintf(SF,"reduce(ToddHySurf*HySurf,chow)/norm,");
		  dprintf(SF,"\",\", ""\"chi_1:\",");
		  dprintf(SF,"reduce(tmp*ToddHySurf*HySurf,chow)/norm,");
		  dprintf(SF,"\" [\",");
		  dprintf(SF,"reduce(c%dHySurf*HySurf,chow)/norm,",DIM);
		  dprintf(SF,"\"]\";\n");	  
	}

/*#############################################*/

/*####### Intersection Polynomial #############*/
  if(_Flag->i){	
	  dprintf(SF,"\"SINGULAR -> intersection polynomial:\";\n");
	  dprintf(SF,"poly f=1-t*(reduce(%s%d,std(hypideal))",B,1);
	  for (j=2; j<=p-d; j++) {
		  dprintf(SF,"+reduce(%s%d,std(hypideal))",B,j);
	  }
	  dprintf(SF,");\n");
	  dprintf(SF,"poly m=jet(1,f,%d)-jet(1,f,%d);\n",2*DIM,2*DIM-2);
	  dprintf(SF,"matrix M=coef(m,%s%d",B,1);
	  for (j=2; j<=p-d; j++) {
		  dprintf(SF,"*%s%d",B,j);
	  }
	  dprintf(SF,");\n");
	  dprintf(SF,"for(i=1;i<=ncols(M);i++){\n");
	  dprintf(SF,"  M[2,i]=reduce(M[1,i]*HySurf,chow)/norm;\n");
	  dprintf(SF,"}\n");
	  dprintf(SF,"poly p=M[2,1]*M[1,1];\n");
	  dprintf(SF,"for(i=2;i<=ncols(M);i++){\n");
	  dprintf(SF,"  p=p+M[2,i]*M[1,i];\n");
	  dprintf(SF,"}\n");
	  dprintf(SF,"p;\n");
  }
	

  //PRINT Chern classes of CY/H
  if(_Flag->H !=1 && _Flag->c){
	  dprintf(SF,"\"SINGULAR  -> Chern classes of the CY-hypersurface:\";\n");
	  for(i=1;i<=DIM-1;i++){
		  dprintf(SF,"\"c%d(CY)= \",reduce(c%dHySurf,chow);\n",i,i);
	  }
	  dprintf(SF,"\"c%d(CY)= \",reduce(c%dHySurf*HySurf,chow)/norm,\"*[pt]\";\n",DIM,DIM);
  }
  if(_Flag->H==1 && _Flag->c){
	  dprintf(SF,"\"SINGULAR  -> Chern classes of the hypersurface H:\";\n");
	  for(i=1;i<=DIM-1;i++){
		  dprintf(SF,"\"c%d(H)= \",reduce(c%dHySurf,chow);\n",i,i);
	  }
	  dprintf(SF,"\"c%d(H)= \",reduce(c%dHySurf*HySurf,chow)/norm,\"*[pt]\";\n",DIM,DIM);
  }

	
/*####### triple intersection numbers #########*/
  if(_Flag->t){
	  dprintf(SF,"\"SINGULAR -> triple intersection numbers:\";\n");
	  dprintf(SF,"poly f=1;\n");
	  dprintf(SF,"for(i=1;i<=%d;i++){\n",p);
	  dprintf(SF,"  k=1;\n");
	  dprintf(SF,"  for(j=1;j<=size(nonintersectingd);j++){\n");
	  dprintf(SF,"    if(dli[i]==nonintersectingd[j]){k=0;}\n");
	  dprintf(SF,"  }\n");     
	  dprintf(SF,"  f=f-t*k*dli[i];\n");
	  dprintf(SF,"}\n");
	  dprintf(SF,"poly m=jet(1,f,%d)-jet(1,f,%d);\n",2*DIM,2*DIM-2);
	  dprintf(SF,"matrix M=coef(m,dli[%d]",1);
	  for (j=2; j<=p; j++) {
		  dprintf(SF,"*dli[%d]",j);
	  }
	  dprintf(SF,");\n");
	  dprintf(SF,"for(i=1;i<=ncols(M);i++){\n");
	  dprintf(SF,"  M[2,i]=reduce(M[1,i]*HySurf,chow)/norm;\n");
	  dprintf(SF,"}\n");
	  dprintf(SF,"string colon=\",\";\n");
	  dprintf(SF,"string arrow=\"->\";\n");
	  dprintf(SF,"for(i=1;i<=ncols(M);i++){\n");
	  dprintf(SF,"  string(M[1,i],arrow,M[2,i],colon);}\n");
  }
	

	
	
/*####### Top. quant'ies of the divisors '#######*/
  if(_Flag->d && DIM>1){
	dprintf(SF,"\"SINGULAR -> topological quantities of the toric divisors:\";\n");
	for(i=1; i<=p; i++){
		dprintf(SF,"k=1;\n");
		dprintf(SF,"for(j=1;j<=size(nonintersectingd);j++){\n");
		dprintf(SF,"    if(dli[%d]==nonintersectingd[j]){k=0;}\n",i);
		dprintf(SF,"}\n");     
		dprintf(SF,"if(k==1){\n");
		dprintf(SF,"poly ChernBottomD%d=1;\n",i);
		dprintf(SF,"for(i=1;i<=%d;i++){ChernBottomD%d=ChernBottomD%d+(-1)^i*(dli[%d])^i;}\n",DIM-1,i,i,i);
		dprintf(SF,"poly ChernD%d=ChernUp*ChernBottomD%d*ChernBottomHypersurface;\n",i,i);
		/*Chern classes of divisors*/
		for (j=0;j<DIM-1; j++) dprintf(SF,"poly c%dD%d=jet(ChernD%d,%d)-jet(ChernD%d,%d);\n",j+1,i,i,j+1,i,j);
		dprintf(SF,"matrix ch[%d][1] = 1",DIM);
		for(j=2;j<=DIM;j++){
			dprintf(SF,",c%dD%d",j-1,i);
		};
		dprintf(SF,";\n");
		dprintf(SF,"matrix pp[%d][1] = -ch[2,1]",DIM-1);
		for(j=2;j<=DIM-1;j++){
			dprintf(SF,",-%d*ch[%d,1]",j,j+1);
		};
		dprintf(SF,";\n");
		dprintf(SF,"for(i=1;i<=%d;i++){\n",DIM-1);
		dprintf(SF,"  for(j=1;j<=i-1;j++){\n");
		dprintf(SF,"    pp[i,1]=pp[i,1]-pp[j,1]*ch[i-j+1,1];\n");
		dprintf(SF,"  };\n");
		dprintf(SF,"};\n");
		dprintf(SF,"poly tmp=%d;\n",DIM-1);
		dprintf(SF,"for(i=1;i<=%d;i++){\n",DIM-1);
		dprintf(SF,"  tmp = tmp+pp[i,1]/factorial(i)*(-t)^i;\n");
		dprintf(SF,"};\n");
		dprintf(SF,"poly ToddD%d=subst(Todd(tmp,%d),t,1);\n",i,DIM-1);
		  
		/*Euler chrachteristics and arithmetic genera*/
		dprintf(SF,"poly EulerD%d=reduce(dli[%d]*c%dD%d*HySurf,chow)/norm;\n",i,i,DIM-1,i);
		dprintf(SF,"poly EulerHD%d=reduce(ToddD%d*dli[%d]*HySurf,chow)/norm;\n",i,i,i);
		dprintf(SF,"};\n");
	}
	dprintf(SF,"string EulerD;");
	dprintf(SF,"string EulerHD;");
    dprintf(SF,"list c1D; list eulerHD; list eulerD;");
	for(i=1;i<=p;i++){
		dprintf(SF,"if(defined(EulerD%d)){\n",i);
  		dprintf(SF,"  EulerD=EulerD+string(EulerD%d)+\" \";",i);
		dprintf(SF,"  EulerHD=EulerHD+string(EulerHD%d)+\" \";",i);
		dprintf(SF,"  c1D[%d]=c1D%d;",i,i);
		dprintf(SF,"  eulerHD[%d]=EulerHD%d;",i,i);
		dprintf(SF,"  eulerD[%d]=EulerD%d;",i,i);
        dprintf(SF,"};\n");		
	}
	dprintf(SF,"\"Euler characteristics:\", EulerD;\n ");
	dprintf(SF,"\"Arithmetic genera:\", EulerHD;\n");
/*#############################################*/

/*####### dP conditions #######################*/
	  if(DIM == 3){
		  dprintf(SF,"kill k;\n");
		  dprintf(SF,"if(defined(f)){kill f;}\n");
		  dprintf(SF,"if(defined(p)){kill p;}\n");
		  dprintf(SF,"list dP; list N; int l=0; int p=0; int t=0; int s=0; int k=0; int u=0; int v=0;");
			dprintf(SF,""
					"for(i=1;i<=size(dli);i++){"
					"if(reduce(dli[i],std(hypideal))!=0){"
					"if(eulerHD[i]==1 && eulerD[i]>=3 && eulerD[i]<=11 && reduce(c1D[i]^2*dli[i]*HySurf,chow)/norm >=1 && reduce(c1D[i]^2*dli[i]*HySurf,chow)/norm <= 9){"
	 					  "for(j=1;j<=size(dli);j++){"
	 							  "if(dli[j] != dli[i]){"
									  "for(k=1;k<=size(dli);k++){"
										  "v++;"
										  "if(reduce(dli[i]*dli[j]*dli[k]*HySurf,chow)/norm == 0){u++;}"
									  "}"
									  "if(u < v){"
	 								  "if(reduce(dli[j]*dli[i]*c1D[i]*HySurf,chow)/norm <= 0){"
	 									  "p++;"
	 								  "}"
									  "} u=0; v=0;"
	 							  "}"
	 						  "}"
	 					  "if(p==0){"
	 							  "s++;"
	 						  "dP[s]=dli[i];"
	 						  "N[s]=eulerD[i]-3;"
	 					  "}"
	 				  "}"
					  "}"
	 			  "p=0;"
	 			  "t=0;"	/* index for no repetitions of div. classes */
	 			  "}	");

			dprintf(SF,"string DelPezzo;for(i=1;i<=size(dP);i++){DelPezzo=DelPezzo+string(dP[i])+\"(\"+string(N[i])+\")\"+\" \";};");
			dprintf(SF,"list nonintDP; int g=0; int f=0; if(size(dP)>0){if(size(dP)>1){for(i=1;i<=size(dP);i++){for(j=1;j<=size(dP);j++){if(dP[j]!=dP[i]){for(k=1;k<=size(dli);k++){if(reduce(dP[i]*dP[j]*dli[k]*HySurf,chow)/norm!=0){g++;} ;};};}; if(g==0){f++; nonintDP[f]=dP[i];}; g=0;};} else {nonintDP[1]=dP[1];};};");
			dprintf(SF,"string niDP;for(i=1;i<=size(nonintDP);i++){niDP=niDP+string(nonintDP[i])+\" \";};");
			dprintf(SF, "\"dPs:\",size(dP),\";\",DelPezzo,\"nonint:\",size(nonintDP),\";\",niDP;");
	  }
/*#############################################*/
}

  dprintf(SF,"quit;\n");
  close(SF);

#if (TEST_PRINT_SINGULAR_IO)
				CatFile(SFname);
#endif

  if( system(SingularCall) ) {puts("Check Singular installation");exit(1);}
  remove(SFname);
  free(SFname);
  free(SingularCall);
}



