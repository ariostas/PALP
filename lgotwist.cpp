/*  ======================================================================  */
/*  ==========            lgotwist.c     Sep 29 / 02	 =================  */
/*  ======================================================================  */
/*  aao2.6
002244 4 3 4 3 4 3
-21 21 0 1
-19 19 0 2
-9 9 0 3
*/


#include <stdio.h>
#include <string.h>
#ifdef    __MSDOS__
#define   NM   9 /*13*/          /*  maximum number of fields               */
#define   WM       1024   /*8192*/  /*  maximum number of words             */
		      /*  (a word is binary code for a subset of X_i) */
#define   NS   10 /*14*/         /*  maximum number of generators           */
#define   HODDIM   500
#else
#define   NM   13          /*  maximum number of fields               */
#define   NS   14          /*  maximum number of generators           */
#define   WM   8192        /*  maximum number of words             */
#define   HODDIM   4000
#endif
#define   NPN  1000
#define   NP   10                /*  maximum number of prime factors        */
#define   LONGOUT  (0)           /*  (1) for long output, otherwise (0)     */
#define   ONLYINV  (0)           /*  (1) iff only invertibles are computed  */
#define   MAXPN    512           /* maximum number of pointers at one point */
#define SAFER  (1)     /* SAFER=1 ==> b01 with rat, less danger of overflow */
#define ADDZD (0)             /* ADDZD = 1 (0) iff Z_d is (not) to be added */
int mask[]={1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384};
int  interb01();                /* ...  interface to "mhodge" conventions   */

FILE *infi, *outfi;
int  stdi=0, bugcount=0, invertible;
            /* stdi=1 if called without "infile"  (see: main, readline) */


/*  ======================================================================  */
/*  ==========             integer and rational stuff           ==========  */
/*  ======================================================================  */
/*  ==========     rat.h   (header -> #include "rat.h")         ==========  */
/*  ======================================================================  */
#define  lint      long
#define  sint      int
#define  mod(a,b)  ((a)%(b))                         /* ((a)-(b)*((a)/(b))) */
#define  MOD(a,b)  ((b) ? mod(a,b) : (a))
#define  min(a,b)  (((a)<(b)) ? (a) : (b))
#define  max(a,b)  (((a)>(b)) ? (a) : (b))
#define  lcm(a,b)  ((a)*(b)/gcd((a),(b)))
#define  plcm(a,b) (a)*((b)/pgcd(a,b))
#define  abs(a)    (((a)<0) ? (-(a)) : (a))

lint  gcd(lint a, lint b); /* modulus of greatest common div.; gcd(0,n)=|n| */
typedef   struct {lint num; lint den;} rat;
rat   rI(lint a);          /*  conversion  lint -> rat  */
rat   rR(lint a, lint b);  /*  conversion  a/b  -> rat  */
rat   rS(rat a, rat b);    /*  a + b  */
rat   rD(rat a, rat b);    /*  a - b  */
rat   rP(rat a, rat b);    /*  a * b  */
rat   rQ(rat a, rat b);    /*  a / b  */
void  fS(rat *a, rat *b);  /*  fast sum: add rat  *b  to rat  *a  */
void  iS(rat *a, int *b);  /*  fast sum: add int  *b  to rat  *a  */
/*  ======================================================================  */
/*  ========== rat.c (source code: cc -c rat.c generates rat.o) ==========  */
/*  ========== cc -o prog prog.c rat.o   incl. obj. code rat.o  ==========  */
/*  #include "rat.h"      ================================================  */
rat   rI(register lint a) {rat c; c.num=a; c.den=1; return c; }    /*  a/1  */
rat   rR(register lint a, register lint b)                         /*  a/b  */
{     register lint g=gcd(a,b); rat c; g=(b<0) ? -g:g; c.num=a/g;c.den=b/g;
      return c;
}
rat   rS(rat a, rat b)                                           /*  a + b  */
{     register rat c; register lint g=gcd(a.den,b.den);
      g=gcd(c.den=a.den*(b.den/g),c.num=a.num*(b.den/g)+b.num*(a.den/g));
      c.num/=g; c.den/=g; return c;
}
rat   rD(rat a, rat b)                                           /*  a - b  */
{     register rat c; register lint g=gcd(a.den,b.den);
      g=gcd(c.den=a.den*(b.den/g),c.num=a.num*(b.den/g)-b.num*(a.den/g));
      c.num/=g; c.den/=g; return c;
/**      register rat c; register lint g=gcd(c.num=a.num*b.den-b.num*a.den,
      c.den=a.den*b.den);      c.num/=g; c.den/=g; return c; */
}
rat   rP(rat a, rat b)                                           /*  a * b  */
{     register lint g=gcd(a.num,b.den); register lint h=gcd(b.num,a.den);
      register rat c; c.num=(a.num/g)*(b.num/h); c.den=(a.den/h)*(b.den/g);
      return c;
}
rat   rQ(rat a, rat b)                                           /*  a / b  */
{     register lint g=gcd(a.num,b.num); register lint h=gcd(b.den,a.den);
      register rat c; c.num=(a.num/g)*(b.den/h); c.den=(a.den/h)*(b.num/g);
#ifdef TEST
    if (!c.den)
    fprintf(outfi,"warning: vanishing denominator in rQ in %s!\n",infun);
#endif
      if (c.den<0) {c.num=-c.num; c.den=-c.den;} return c;
}
lint  gcd(register lint a, register lint b)
{     a = (a<0) ? -a:a; b = (b<0) ? -b:b; if (!a) return b; if (!b) return a;
       {register lint c;  while(c=mod(a,b))  {a=b;b=c;} return  b;}
}
lint  argint(char* s)                                  /* string to integer */
{     lint d=(*s)-'0'; sint i=1; while(s[i]) d=10*d+s[i++]-'0'; return d;
}
void  prat(rat a) {fprintf(outfi,"%ld/%ld ",a.num,a.den);}
/*  ======================================================================  */
/*  ==========          end of integer and rational stuff       ==========  */
/*  ======================================================================  */

/****************************************************************************/
/*     Construction of the maximal abelian symmetry                         */
/****************************************************************************/

typedef sint spectrum[3];               /*  b01 change  */
typedef sint pointlist[NM+1];
typedef lint symlist[NS];
typedef sint primelist[NP];
typedef sint prsymlist[NS][NP];
typedef symlist symsymlist[NS];
typedef lint sympointlist[NM+1][NS];
typedef sint prsympointlist[NM+1][NS][NP];
typedef sint prsymsymlist[NS][NS][NP];
typedef lint worte[WM+1][3];

lint lcmd, specnum, symnum, totsymnum=0, totspecnum=0, maxsymnum=0,
      maxspecnum=0, modelnum=0;
    /* lcmd = least common multiple of the orders of the generators of the
	  maximal group                                                 */
sint n, norig, ns, npr, addsyms, evenn, oddd, over=1, D;  /* D=sum(1-2q_i) */
	/* n = number of points X_i;                                    */
	/* ns = number of phase symmetries, ns<NS                       */
	/* npr = number of prime numbers involved, npr<NP               */
	/* redundant is the number of det=1 symmetries acting only in the
	   trivial sector                                               */
symlist d, det;
	     /* d[0]=d, d[i]=order of i'th phase symmetry           */
	     /* det[i]/d[i] = determinant of i'th symmetry          */
/* the letters 'pr' in the name of a variable/field indicate that it is
   the prime decomposed version of some other v./f., x=prod_i{prx[i]};
   for explanations of spr.., npr.., aux.. and norm see reccon and fillup   */
rat prdet[NS][NP];
primelist pr, prns, auxns;
    /* pr is the list of prime numbers occurring in a specific model        */
prsymlist prd, sprd, norm, auxd;
sympointlist wei;
       /*    n_i=wei[i-1][0]; phase of X_i under k'th symmetry =wei[i-1][k] */
prsympointlist prwei, nprwei, auxwei;
prsymsymlist sprwei;
spectrum hodlist[HODDIM];

/*  ======================================================================  */
/*  ==========                abelmax                           ==========  */
/*  ======================================================================  */
typedef struct {int p[NM], a[NM], N;} skelet;    /* p: pointer; a: exponent */
typedef struct {int p[NPN], m[NPN];} prili;/* p[0]=#(primes); m[0]=max mult.*/

int readline(skelet *);      /* reads: "string[#+1] exp_0 ... exp_# ... \n" */
void analy(skelet); /* calculate: order of evaluation; loops; pointed at by;*/
prili pmax;  /* global var. for prime decomposistion of lcm of group orders */
void printpri(prili);

int readline(skelet *s)      /* reads: "string[#+1] exp_0 ... exp_# ... \n" */
{   int i; s->N=0;       if(stdi) printf("skeleton? ");
    while(' '!=(i=fgetc(infi))) if(i==EOF) return 0; else s->p[s->N++]=i-'0';
    for(i=0;i<(s->N);i++) fscanf(infi,"%d",&s->a[i]);
    n=(*s).N;
    while(fgetc(infi)-'\n'); return 1;
}

long pgcd(register long a, register long b)
{    register long c;   while(c=mod(a,b)) {a=b;b=c;} return  b;
}

int prime[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,
89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,
307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,
421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,
547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,
659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,
797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,
929,937,941,947,953,967,971,977,983,991,997};

void printpri(prili li)
{    int i; fprintf(outfi,"\nlist of primes:\n");
     for(i=1;i<=*li.p;i++) {fprintf(outfi,"%d,",li.p[i]);
      if(!mod(i,20)) fprintf(outfi,"\n");}
     fprintf(outfi,"\nmaximal multiplicities (in one generator):\n");
     for(i=1;i<=*li.p;i++) {fprintf(outfi,"(%d,%d) ",li.p[i],li.m[i]);
      if(!mod(i,10)) fprintf(outfi,"\n");}
     fprintf(outfi,"\n");
}

prili prideco(long x)
{     int n=0,p; long q; prili l; *l.p=*l.m=1; l.m[1]=0;
      while(x>=((p=prime[n++])^2)) {
       while(!(x-p*(q=x/p))) {x=q; l.m[*l.p]++;}
       if(l.m[*l.p]){l.p[*l.p]=p;*l.m=max(*l.m,l.m[*l.p]);l.m[++(*l.p)]=0;}}
      if(x==1) (*l.p)--; else {l.m[*l.p]=1; l.p[*l.p]=x;}   return l;
}

void maxpri(prili pn) /* pmax = prime decomposistion of lcm of group orders */
{    int i=1,j=1,n=1; prili po=pmax; *pmax.m=max(*po.m,*pn.m);
     po.p[*po.p+1]=pn.p[*pn.p+1]=po.p[*po.p]+pn.p[*pn.p]+1;
     if(*po.p){while((i<=*po.p)||(j<=*pn.p)) {
       if(po.p[i]<pn.p[j]) {pmax.p[n]=po.p[i];pmax.m[n]=po.m[i++];}
       else if(po.p[i]>pn.p[j]) {pmax.p[n]=pn.p[j];pmax.m[n]=pn.m[j++];}
       else{pmax.m[n]=max(po.m[i],pn.m[j]);pmax.p[n]=po.p[i++];j++;}n++;}
       *pmax.p=n-1;}
     else pmax=pn; /* printpri(pmax); */
}

/*   analy first finds all minimal member of loops (or fermats): value = i  *
 *   for any such i its connected component of the graph is reconstructed:  *
 *   i=lo[1] -> lo[2] -> ... -> lo[*lo] -> i;                               *
 *   the remaining variables are on "ord[1],...,ord[*ord]" such that all    *
 *   pointers go from right to left (or from ord to lo, of course).         *
 *   inv[k][1],...,inv[k][*inv[k]] are all variables pointing at k.         */
/*   Now go thru ord from right to left and calculate lcm=prod(t) with t_k  *
 *   dividing the corresponding exponent s.a[], divide the orders of the    *
 *   groups corresponding to inv[] by t_k and evaluate the new generator    */
/*   finally repeat this for the loop of i and evaluate its symm.-generator */
void analy(skelet s)
{    long ph[NM][NM+1],      /* p[i][0] = order of sym. i<s.N */
      d,num[NM],den[NM];/* (n_i)/d == (num_i)/(den_i) */
     int i,j,n,ord[NM+1],  /* order of evaluation (right to left) */
       lo[NM+1],lopo,   /* loop of length *lo; lopo=#(pointers at loop) */
       inv[NM][NM];    /* j is pointed at by inv[j][1],...,inv[j][*inv[j] */
     for(i=0;i<s.N;i++) {*ph[i]=1;for(j=0;j<s.N; ph[i][++j]=0);} /* init. ph*/
     for(i=0;i<s.N;i++) *inv[i]=0;    /* add pointer at s.p[i]: */
     for(i=0;i<s.N;i++) if(s.p[i]!=i) {inv[s.p[i]][++(*inv[s.p[i]])]=i;}
     invertible=1; for(i=0;(i<s.N)&&invertible;i++)if(*inv[i]>1) invertible=0;
     for(i=0;i<s.N;i++) { /* find j in comp. of i: j==i iff min.loop memb.: */
     n=s.N; j=s.p[i]; while((--n)&&(j>i)) j=s.p[j];
     if(j==i){                   /* check for i == minimal loop member '{{' */
      int m=0; *lo=s.N-n;                  /* comlete connected component */
      if(*lo==1) {*lo=0; n=2;ord[m=1]=i;}      /* if i fermat -> no loop! */
      else {lo[1]=i; n=1;       /* initialize lo[] and evaluation of ord: */
          for(j=2;j<=*lo;j++) lo[j]=s.p[lo[j-1]];
          while(m<*inv[i]) {if(lo[*lo]!=(ord[n]=inv[i][++m])) n++;}
          for(j=2;j<=*lo;j++) {m=0; while(m<*inv[lo[j]])
                  if(lo[j-1]!=(ord[n]=inv[lo[j]][++m])) n++;}}
      lopo=n-1;   /* now ord[1],...,ord[lopo] are the pointers at loop(i) */
      for(m=1;m<n;m++) for(j=0;j<*inv[ord[m]];) ord[n++]=inv[ord[m]][++j];
                             /* end of: comlete connected component */
/*   still with i == minimal loop member:                                   */
/*   calculate symmetries: first for vars. on ord[], then for the loop      */
     for(*ord=--n;n;n--) {   /* go backwards thru trees according to "ord"  */
      int *nt=inv[ord[n]];  /* nt points at the vector of ptrs. at ord[n] */
      if(*nt) {prili pli; long tli[NM], lcm=*ph[nt[1]]; tli[1]=1;
         for(j=2;j<=*nt;j++){tli[j]=1;lcm*=*ph[nt[j]]/pgcd(lcm,*ph[nt[j]]);}
       pli=prideco(lcm);
       for(j=1;j<=*pli.p;j++) {int po=pli.p[j]; m=0;
             while(++m<pli.m[j]) po*=pli.p[j];
           for(m=1;mod(*ph[nt[m]],po);m++); tli[m]*=po;}  /* choose t_m */
       *ph[ord[n]]=s.a[ord[n]]*lcm;    /* assign order to new generator */
       ph[ord[n]][ord[n]+1]=lcm;
       for(m=1;m<=*nt;m++){long fac=lcm/(*ph[nt[m]]);for(j=n+1;j<=*ord;j++)
                   ph[ord[n]][ord[j]+1]-=fac*ph[nt[m]][ord[j]+1];}
       for(m=1;m<=(*nt);m++) *ph[nt[m]]/=tli[m];/* divide orders by t_m */
       }
      else {*ph[ord[n]]=s.a[ord[n]]; ph[ord[n]][ord[n]+1]=1;}}
/*   now (if *lo>0) calculate symmetry for the loop of the component of i:  */
/*   remember that ord[1],...,ord[lopo] point at the loop                   */
     if(*lo>1) {prili pli; long P=s.a[lo[*lo]],b[NM+1], tli[NM], lcm;
      b[1]=1;     for(j=1;j<*lo;++j) b[j+1]=-b[j]*s.a[lo[j]];
      P *= b[j]; if((++P)<0) P=-P; den[i]=P;    /* P = prod(exp(j)) +/- 1 */
      lcm=1;                                    /*  nt[*nt] -> ord[lopo]  */
      for(j=0;(j++)<lopo;tli[j]=1) lcm*=*ph[ord[j]]/pgcd(lcm,*ph[ord[j]]);
      pli=prideco(lcm);
      for(j=1;j <= *pli.p;j++) {int po=pli.p[j]; m=0;
        while(++m<pli.m[j]) po*=pli.p[j];
        for(m=1;mod(*ph[ord[m]],po);m++); tli[m]*=po;}    /* choose t_m */
      *ph[i]=P*lcm;                      /* assign order to new generator */
      for(n=1;n<=*lo;n++) ph[i][lo[n]+1]=b[n]*lcm;
      for(n=1;n<=lopo;n++) {long fac=1; while(lo[fac]!=s.p[ord[n]]) {fac++;
/* if(fac>NM) {printf("infinit loop(1): fac=%d ",fac);return;} ... debug */}
         fac=b[fac]; while(pgcd(abs(fac),*ph[ord[n]])>1) {fac-=P;
/* if(abs(fac/P)>NM) {printf("infinit loop(2): fac=%d ",fac);return;} debug*/}
       fac*=lcm/(*ph[ord[n]]);
       for(j=1;j<=*ord;j++) ph[i][ord[j]+1]-=fac*ph[ord[n]][ord[j]+1];}
      for(j=1;j<=lopo;j++) *ph[ord[j]]/=tli[j];   /* divide orders by t_j */
     }  maxpri(prideco(*ph[i]));
     if(*lo<2) {num[i]=1;den[i]=s.a[i];} /* calculate weights num[j]/den[j] */
     else {num[i]=1-s.a[lo[2]]; for(j=2;j<*lo;num[i]=1-s.a[lo[++j]]*num[i]);
         d=pgcd(den[i],num[i]=abs(num[i])); num[i]/=d; den[i]/=d;}
     for(j=1;j<*lo;j++) {den[n=lo[j+1]]=den[lo[j]]; d=pgcd(den[n],
       num[n]=den[n]-s.a[lo[j]]*num[lo[j]]); num[n]/=d;den[n]/=d;}
     for(j=1+(*lo==0);j<=*ord;j++){den[ord[j]]=s.a[ord[j]]*den[n=s.p[ord[j]]];
       d=pgcd(den[ord[j]],num[ord[j]]=den[n]-num[n]);
       num[ord[j]]/=d;den[ord[j]]/=d;}   /* end of calculation of weights */
     }}                                                             /* '}}' */
     d=*den; for(j=1;j<s.N;j++) d*=den[j]/pgcd(d,den[j]);/* calculate d and */
     D=0; for(j=0;j<s.N;j++) D+=(num[j]=num[j]*(d/den[j])); /* n_j==num[j]  */
     if((D*=2)% d)fprintf(outfi,"D not int!\n"); D=s.N-D/d; /* D=sum(1-2q_i)*/
for(j=0;j<s.N;j++)if( (num[j]*(s.a[j]-(s.p[j]==j))) != (d-num[s.p[j]]) )
	{printf("Error in calculation of weights!\n"); /* check weights */
	fprintf(outfi,"Error in calculation of weights!: ");bugcount++;}
/*   throw away order 1; write result to outfi (generators, skeleton, #gen.)*/
     {int k,l; long G=1; n=0;
     if(LONGOUT){fprintf(outfi,"N=%d ", s.N); fprintf(outfi,"d=%ld, n_i=",d);}
     /* else {fprintf(outfi,"%d ", s.N); fprintf(outfi,"%ld",d);} */
     wei[s.N][0]=d;
     for(j=0;j<s.N;j++){if(LONGOUT)fprintf(outfi," %ld",num[j]);
      wei[j][0]=num[j];   }
     for(j=0;j<s.N;j++) if((*ph[j])>1) {
      n++;
      if(LONGOUT)fprintf(outfi,"\n%d->Z[%ld]: ",j,*ph[j]);
      wei[s.N][n]=ph[j][0];
      for(k=1;k<=s.N;k++){
       wei[k-1][n]=ph[j][k];
       if(LONGOUT)fprintf(outfi," %ld",ph[j][k]);}
      G*=*ph[j];}
     ns=n;
     if(LONGOUT) {fprintf(outfi,"\nskeleton =");  fprintf(outfi," ");}
     if(invertible || !ONLYINV) {
	for(l=0;l<s.N;l++) fprintf(outfi,"%d",s.p[l]);
	for(l=0;l<s.N;l++) fprintf(outfi," %d", s.a[l]);
	fprintf(outfi," inv=%d ", invertible);  }
     if(LONGOUT)fprintf(outfi,"  #generators=%d, Order=%ld\n",n,G);  }
}

/************  End of maxsym                                        *********/


/****************************************************************************/
/*    prime decomposition and reduction to allowed determinants             */
/****************************************************************************/

void primedecomp()
/* decomposition of wei, d, det into their components corresponding to
   the primenumbers pr[j], namely prwei[j], prd[j], prdet[j]                */
{     sint i, j=0, k, l, p;
      lint dk;
      for (i=0;(i<NPN)&&(lcmd!=1);i++){
       if (!mod(lcmd,p=prime[i])) {
        pr[j]=p;
        for (k=0;k<=ns;k++){
	   prd[k][j]=1; dk=d[k];
           while (!mod(dk,p)){dk/=p; prd[k][j]*=p;};
           for (l=0;l<n;l++) prwei[l][k][j]=mod(wei[l][k],prd[k][j]);
           prdet[k][j]=rR(mod(det[k],prd[k][j]),prd[k][j]);       };
        while (!mod(lcmd,p)) lcmd/=p;
        j++;   }   }
      npr=j;
      if (lcmd!=1) {fprintf(outfi,"caution: Big prime number!!!"); bugcount++;}
}

void gooddets()
/* recombines the generators in such a way that they have det=1 (det^d=1 for
   torsion), orders them (decreasing orders), counts them                   */
{     sint i,j,k,l,ig,is,imax,maxdetden,w;
       /* maxdetden is the maximal denominator of the determinants        */
      lint auxlong;
      for (j=0;j<=npr;j++){
       maxdetden=1;
       for (i=1;i<=ns;i++) maxdetden=max(maxdetden, prdet[i][j].den);
       while (maxdetden>1 /* prd[0][j] for torsion  */){
	imax=0;
        for (i=1;i<=ns;i++) if (prdet[i][j].den==maxdetden){
           if (!imax) imax=i;
           else {
            if (prd[i][j]>prd[imax][j]) {ig=i; is=imax;}
              else {ig=imax; is=i;};
            k=0;
            while (mod(prdet[ig][j].num+k*prdet[is][j].num,maxdetden))
                   k++;
            for (l=0;l<n;l++){
             auxlong=k;
             auxlong*=prwei[l][is][j];
             auxlong*=prd[ig][j];
             auxlong/=prd[is][j];
             auxlong+=prwei[l][ig][j];
             prwei[l][ig][j]=mod(auxlong,prd[ig][j]);}
            prdet[ig][j]=rI(0);
            imax=is;    }   }
        prd[imax][j]/=pr[j];
        for (l=0;l<n;l++)
	    prwei[l][imax][j]=mod(prwei[l][imax][j], prd[imax][j]);
          prdet[imax][j]=rP(prdet[imax][j],rI(pr[j]));
        maxdetden/=pr[j];     }
       /* ordering of the generators:                                     */
       for (i=1;i<ns;i++){
        imax=i;
        for (k=i+1;k<=ns;k++) if (prd[k][j]>prd[imax][j]) imax=k;
        if (imax>i) {
             for (k=0;k<n;k++){
              w=prwei[k][imax][j];
            prwei[k][imax][j]=prwei[k][i][j];
            prwei[k][i][j]=w;    }
           w=prd[imax][j];
             prd[imax][j]=prd[i][j];
             prd[i][j]=w;
           w=prdet[imax][j].den;
           prdet[imax][j].den=prdet[i][j].den;
             prdet[i][j].den=w;
             w=prdet[imax][j].num;
             prdet[imax][j].num=prdet[i][j].num;
	   prdet[i][j].num=w;    }   }
       /* calculation of number prns[j] of non-trivial generators         */
       i=1;
       while ((prd[i][j]>1)&&(i<ns+1)) i++;
       prns[j]=i-1;            }
}

/***** End of primedecomp&gooddets ******************************************/


/****************************************************************************/
/*     Routines for checking the link criterion                             */
/****************************************************************************/

typedef int smon[1];

pointlist pointernum;
int linklist[NM][MAXPN], targlist[NM][MAXPN], monlist[NM][MAXPN];
       /*    #pointers at X_i = pointernum[i-1]                             */
       /*    linklist[i-1] contains list of all links pointing at X_i       */
       /*    targlist[i][j] indicates subtargets of linklist[i][j]          */
       /*    monlist[i][j] is the actual monomial realising linklist[i][j]  */

sint symcheck(symlist sum, int link, smon mon){
/* symcheck checks whether there is a monomial mon in
 * the variables indicated by link whose total weight is sum[0] and which
 * transforms under the k'th symmetry with a phase sum[k];
 * if X_i occurs in the monomial then mon[i] is set to 1.                   */
   sint i, j, k, check, expo;
   int newlink;
   symlist newsum;
   for (i=0;i<n;i++) if (link&mask[i]){
      if (!mod(sum[0], wei[i][0])){
       expo=sum[0]/wei[i][0];
       check=1;
       for (j=1;(j<=ns)&&check;j++)
          if (mod(expo*wei[i][j]-sum[j],d[j])) check=0;
       if (check) {mon[0]=mon[0]|mask[i]; return 1;};   };
       newlink=link-mask[i];
       for (j=0;j*wei[i][0]<sum[0];j++){
        for (k=0;k<=ns;k++) newsum[k]=sum[k]-j*wei[i][k];
        if (symcheck(newsum,newlink,mon)) {
           if (j) mon[0]=mon[0]|mask[i];
           return 1; };   };   };
   return 0;
}

sint checklink(int link, int targets){
/* checklink checks recursively whether a specific link with subtargets
 * indicated by targets exists and can be added to the graph.
 * If a pointer is required, the existence of all further links required
 * by our theorem is checked by recursive calls of checklink.               */
   sint i, j, k, pn, check;
   int newtarg, newlink;
   smon mon;
   symlist dw;
   if (symcheck(d,link,mon)) return 0;
   for (i=0;i<n;i++) if(!((targets|link)&mask[i])){
      mon[0]=0;
      for (j=0;j<=ns;j++) dw[j]=d[j]-wei[i][j];
      if (symcheck(dw,link,mon)){
       for (j=0;j<pointernum[i];j++){
           newlink=(link|linklist[i][j]);
           newtarg=(targets|targlist[i][j]);
           if (mon[0]!=monlist[i][j])
            if (checklink(newlink, newtarg)) return 1; };
       pn=pointernum[i]++;
       if (pn>=MAXPN) {
	  printf("pointernum too large\n");
          return 1;}
       monlist[i][pn]=mon[0];
       linklist[i][pn]=link;
       targlist[i][pn]=targets|mask[i];
       return 0;   };   };
   return 1;
}

sint checkweight(){ /* checks whether our weight system allows a
		   non-degenerate symmetry-respecting potential         */
   sint i, j;
   for (i=0;i<n;i++) pointernum[i]=0;
   for (i=0;i<n;i++) if (checklink(mask[i], 0)) return 1;
   return 0;
}

/****************** End of linkcheck part ***********************************/

/****************************************************************************/
/*      Calculation of Hodge numbers                                        */
/****************************************************************************/

void addhod(spectrum hodge);

void proced()
{     int i, j, k, fac, wort, chi, nvar, expo;
      spectrum spec;
      lint zsum1=0, zsum2=0, mo=1, ng, ngb, kgV=1;
      rat prod;
      symlist omega, ele;   /* omega_i = kgV/O_i, ele encodes group element */
      worte wo;        /*  wo[word][0(2)] = contribution of word to ng(ngb) */
		   /* wo[word][1] = number of group elements whose sets
			       of survivors are indicated by word */
      for (i=0;i<=WM;i++) for (k=0;k<=2;k++)  wo[i][k] = 0;
      for (k=1;k<=ns;k++) kgV = lcm(kgV,wei[n][k]);
      for (k=1;k<=ns;k++) { omega[k] = kgV/wei[n][k]; mo*=wei[n][k]; }
      ele[j=ns]=0;
      while(j<=ns){
      if(ele[j]==wei[n][j])ele[++j]++;
      else if(j-1) ele[--j]=0;
      else {
	 nvar=0;
	 for (i=0;i<n;i++) {
	    expo=0;
	    for (k=1;k<=ns;k++) expo+=ele[k]*omega[k]*wei[i][k];
	    if (!mod(expo,kgV))  nvar+=mask[i];     }
	 wo[nvar][1]++;
	 ele[1]++;}}
      for (i=0;i<mask[n];i++) {
       fac = 1;
       for (j=0;j<n;j++) if ((mask[j] & i) == mask[j])  fac*=-1;
       if (wo[i][1]) for (k=0;k<mask[n];k++) if (wo[k][1]) {
	    wort = i & k;
	    wo[wort][fac+1]+=wo[i][1]*wo[k][1]; }  }
      for (i=mask[n]-1;i>=0;i--) if (wo[i][0]+wo[i][2]) {
       prod=rI(1);
       for (j=0;j<n;j++) if ((mask[j]&i)==mask[j])
	  prod = rP(prod,rR(wei[j][0]-wei[n][0],wei[j][0]));
       if (prod.den != 1){
	  fprintf(outfi,"\ncaution: prod.den != 1\n");bugcount++;}
       zsum1+= wo[i][0]*prod.num; zsum2+= wo[i][2]*prod.num;        }
      spec[2]=interb01();                      /*  b01 change  */
/*      over = wo[mask[n]-1][1];                                            */
      spec[2]/=over; /* bisher b01 nicht durch over dividiert; rueckrechnen:*/
		     /* n n 0 2 -> n+2 n+2 0 1 und n n 0 6 -> n+6 n+6 0 3.  */
      if ((spec[2]!=0)&&(spec[2]!=1)&&(spec[2]!=3)&&(D==3)) {
	  fprintf(outfi,"caution: b01=%d\n", spec[2]); bugcount++;}
      ng = -(zsum1/mo/over+2)/2; ngb = (zsum2/mo/over-2)/2;
      if (zsum1 != -2*(ng+1)*mo*over)
      {fprintf(outfi,"ng is not integer!\n");bugcount++;}
      if (zsum2 != 2*(ngb+1)*mo*over)
      {fprintf(outfi,"ngb is not integer!\n");bugcount++;}
      chi = 2*(ngb-ng);
      if(LONGOUT) fprintf(outfi,"ngb: %ld ng: %ld chi: %d\n",ngb,ng,chi);
      spec[0]=ngb; spec[1]=chi;
      addhod(spec);
}

/*************   End of calculation of Hodge numbers               **********/


/****************************************************************************/
/*    Construction of all symmetries                                        */
/****************************************************************************/

sint forbidden(sint i, sint j, sint sprns)
{     sint l;
      for (l=0;l<sprns;l++) if (i==norm[l][j]) return (l+1);
      return 0;
}

void reccon(sint j, sint sprns);

void processym()
{     sint i, j;
      ns=auxns[npr-1];
      if(LONGOUT){
       fprintf(outfi,"\nwei:\n");
       for(j=0;j<n;j++) fprintf(outfi," %d", wei[j][0]);
       fprintf(outfi,"  %d\n", d[0]);   }
      for (i=1;i<=ns;i++){
       d[i]=auxd[i-1][npr-1];
       wei[n][i]=d[i];
       /* if torsion also dets    */
       for (j=0;j<n;j++){
	  wei[j][i]=auxwei[j][i-1][npr-1];
	if(LONGOUT)fprintf(outfi," %d",wei[j][i]); }
       if(LONGOUT)fprintf(outfi,"  %d\n",d[i]);   }
      if (checkweight()){if(LONGOUT)fprintf(outfi,"degenerate!!\n"); }
      else{if(LONGOUT)fprintf(outfi,"Hodge numbers: "); proced(); symnum++;};
}

sint rectest(sint j, sint k, pointlist auxel)
/* the recursive part of zdtest; checks whether auxel and the generators of
   nprwei with index <= k can combine to the generator of the Z_d           */
{     sint i, m;
      lint auxlong;
      pointlist newauxel;
      if (k<0) {
       for (i=0;i<n;i++)
        if (mod(auxel[i]-prwei[i][0][j],prd[0][j])) return 0;
       return 1;}
      for (m=0;m<prd[0][j];m+=max(1,prd[0][j]/sprd[k][j])){
       for (i=0;i<n;i++){
        auxlong=m; auxlong*=nprwei[i][k][j]; auxlong+=auxel[i];
        newauxel[i]=mod(auxlong,prd[0][j]);}
       if (rectest(j,k-1,newauxel)) return 1;   }
      return 0;
}

sint zdtest(sint j, sint sprns)
/* checks whether the group generated by nprwei contains the generator of
   the Z_d by calling rectest                                               */
{     sint i;
      pointlist auxel;
      if (prd[0][j]<=1) return 1;
      for (i=0; i<n; i++) auxel[i]=0;
      return rectest(j, sprns-1, auxel);
}

void fillup(sint j, sint sprns, sint k, sint l)
/* j: 0..npr-1, sprns: 0..prns[j], k: 0..sprns-1, l: 1..prns[j]             */
/* given norm and ord, fillup recursively constructs all possibilities for
   sprwei, calculates nprwei and auxwei from sprwei and calls (depending on
   j) reccon(j+1,0) or hodge;
   sprwei encodes the pr[j]-components of the symmetry in terms of the
   generators of prwei, nprwei is the explicit form, auxwei is the product of
   the nprwei's for the prime factors up to pr[j];
   the new symmetry is given by auxwei[npr-1]                               */
{     sint imax, m, i;
      lint auxlong;
      if (j>=npr) processym();
      else if (k>=sprns){ if (zdtest(j,sprns)){
       /* if nprwei contains projection of Z_d calculate auxwei:  */
       if (j==0) {imax=0; auxns[0]=sprns;}
       else {imax=min(auxns[j-1],sprns); auxns[j]=max(auxns[j-1],sprns);}
       for (i=0;i<imax;i++){
        for (m=0;m<n;m++)
           auxwei[m][i][j]=sprd[i][j]*auxwei[m][i][j-1]+
            auxd[i][j-1]*nprwei[m][i][j];
	auxd[i][j]=auxd[i][j-1]*sprd[i][j];
           /* if Torsion also auxdet   */    }
       for (i=imax;i<sprns;i++){
        for (m=0;m<n;m++)
           auxwei[m][i][j]=nprwei[m][i][j];
        auxd[i][j]=sprd[i][j];
           /* if Torsion also auxdet   */    }
       if (j>0) for (i=imax;i<auxns[j-1];i++){
	for (m=0;m<n;m++)
           auxwei[m][i][j]=auxwei[m][i][j-1];
	auxd[i][j]=auxd[i][j-1];
           /* if Torsion also auxdet   */    }
      if(LONGOUT){
      fprintf(outfi,"\nsprwei[%d]:\n",j);
      for (i=0;i<sprns;i++){
        for (m=1;m<=prns[j];m++) fprintf(outfi," %d",sprwei[m][i][j]);
        fprintf(outfi,"  %d\n", sprd[i][j]);        }
      fprintf(outfi,"nprwei[%d]:\n",j);
      for (i=0;i<sprns;i++){
        for (m=0;m<n;m++) fprintf(outfi," %d",nprwei[m][i][j]);
	fprintf(outfi,"  %d\n", sprd[i][j]);        }
      fprintf(outfi,"auxwei[%d]:\n",j);
      for (i=0;i<auxns[j];i++){
        for (m=0;m<n;m++) fprintf(outfi," %d",auxwei[m][i][j]);
        fprintf(outfi,"  %d\n", auxd[i][j]);        }     }
       reccon(j+1, 0);     }          }
      else if (l>prns[j]){
       /* calculate nprwei[.][k][j] from sprwei[l][k][j]             */
       for (m=0;m<n;m++){
       nprwei[m][k][j]=0;
       for (i=1;i<=prns[j];i++){
          auxlong=prwei[m][i][j];
          auxlong*=sprwei[i][k][j];
          auxlong*=sprd[k][j];
          auxlong/=prd[i][j];
          nprwei[m][k][j]+=mod(auxlong,sprd[k][j]);}
       nprwei[m][k][j]=mod(nprwei[m][k][j],sprd[k][j]);
           /* if Torsion also nprdet     */      }
       fillup(j, sprns, k+1, 1);   }
      else{
       /* for all possible values of sprwei[l][k][j] call
        fillup(j, sprns, k, l+1)                                   */
       m=forbidden(l, j, sprns);
       if (!m) {
        sprwei[l][k][j]=0;
        if (l<norm[k][j]) while (sprwei[l][k][j]<prd[l][j]){
           fillup(j, sprns, k, l+1);
           sprwei[l][k][j]+=max(1,pr[j]*prd[l][j]/sprd[k][j]);  }
	if (l>norm[k][j]) while (sprwei[l][k][j]<prd[l][j]){
           fillup(j, sprns, k, l+1);
	   sprwei[l][k][j]+=max(1,prd[l][j]/sprd[k][j]);  }   }
       m--;
       if (m==k){
        sprwei[l][k][j]=prd[l][j]/sprd[k][j];
        fillup(j, sprns, k, l+1);    }
       else if (m>=0) {
        sprwei[l][k][j]=0;
        if (l<norm[k][j]) while (sprwei[l][k][j]<prd[l][j]/sprd[m][j]){
	   fillup(j, sprns, k, l+1);
           sprwei[l][k][j]+=max(1,pr[j]*prd[l][j]/sprd[k][j]);  }
	if (l>norm[k][j]) while (sprwei[l][k][j]<prd[l][j]/sprd[m][j]){
           fillup(j, sprns, k, l+1);
           sprwei[l][k][j]+=max(1,prd[l][j]/sprd[k][j]);  }   }    } ;
}

void reccon(sint j, sint sprns)
/* reccon is the starting point for the recursive construction of all
   subgroups of the maximal symmetry group;
   reccon assigns orders sprd[i][j] to the i'th generators of the new group
   (in the sector corresponding to the prime number pr[j]) and indicates by
   norm[i][j] which generator of prwei is normalized in the i'th generator;
   reccon calls fillup to assign values to the other components             */
{     sint i;
      fillup(j,sprns,0,1);
      if (j<npr) for (i=1;i<=prns[j];i++)
      if (!forbidden(i,j,sprns)) {
       norm[sprns][j]=i;
       if (sprns==0) sprd[sprns][j]=prd[i][j];
       else if (i>norm[sprns-1][j])
        sprd[sprns][j]=min(prd[i][j],sprd[sprns-1][j]);
       else sprd[sprns][j]=min(prd[i][j],sprd[sprns-1][j]/pr[j]);
       while (sprd[sprns][j]>1) {
        reccon(j,sprns+1);
        sprd[sprns][j]/=pr[j]; };   }
}

/****** End of construction of symmetries  *********************************/


/****************************************************************************/
/*     Input/output/statistics part                                         */
/****************************************************************************/

void datain()      /*   asks for input, reads input, calculates det         */
{     int i, k, evenn=mod(n-D,2), addtriv=2*((addsyms+1)/2)+evenn;
      norig=n;                           /* n-1 -> n-D   with D=sum(1-2q_i) */
      d[0]=wei[n][0];
      oddd=mod(d[0],2);
      lcmd=d[0];
      det[0]=0;
      for (k=1;k<=ns;k++){
       det[k]=0;
       for (i=0;i<=n;i++) det[k]+=wei[i][k];
       d[k]=wei[n][k];
       det[k]=mod(det[k], d[k]);
       lcmd=lcm(lcmd,d[k]);      }
      if (addtriv){             /* add correct trivial factor + triv. symm. */
       if(oddd){
	for (i=0;i<=n;i++) wei[i][0]*=2;
	d[0]=wei[n][0];  }
       wei[n+addtriv][0] = d[0];
       for (i=0;i<addtriv;i++) wei[n+i][0] = d[0]/2;
       for (k=1;k<=ns;k++){
	wei[n+addtriv][k] = wei[n][k];
	for (i=0;i<addtriv;i++) wei[n+i][k] = 0;      }
       n+=addtriv;
       for (i=1;i<=addsyms+evenn;i++){
	for (k=0;k<n;k++) wei[k][ns+i]=0;
	wei[norig+i-1][ns+i]=1;
	wei[norig+i-1][ns+1]=1;
	wei[n][ns+i]=2;
	d[ns+i]=2;
	det[ns+i]=1;}
       for (k=norig;k<n;k++) wei[k][ns+1]=1;
       det[ns+1]=evenn;
       ns+=addsyms+evenn; }
       specnum=0;
       symnum=0;
}

void longoutput0()
{     sint j, k;
      fprintf(outfi,"\nwei:\n");
       for (j=0;j<=ns;j++){
        for (k=0;k<n;k++) fprintf(outfi," %ld",wei[k][j]);
        fprintf(outfi,"  %ld\n", d[j]);}
}

void longoutput1()
{     sint i, j, k;
      fprintf(outfi,"\nDecomposition into prime numbers:\n");
      for (i=0;i<npr;i++){
       fprintf(outfi,"p=%d:\n",pr[i]);
       for (j=0;j<=ns;j++){
          for (k=0;k<n;k++) fprintf(outfi," %d",prwei[k][j][i]);
          fprintf(outfi,"  %d", prd[j][i]);
          fprintf(outfi,"  det: %ld/%ld\n",prdet[j][i].num,prdet[j][i].den); 
}  }  }

void longoutput2()
{     sint i, j, k;
      fprintf(outfi,"\nSubgroup with det=1, ordered:\n");
      for (i=0;i<npr;i++){
       fprintf(outfi,"p=%d:\n",pr[i]);
       for (j=0;j<=prns[i];j++){
        for (k=0;k<n;k++) fprintf(outfi," %d",prwei[k][j][i]);
        fprintf(outfi,"  %d", prd[j][i]);
        fprintf(outfi,"  det: %ld/%ld\n",prdet[j][i].num,prdet[j][i].den);
};  }  }

sint hodcomp(spectrum hodge1, spectrum hodge2){
      if (hodge1[2]<hodge2[2]) return -1; /*  b01 change  */
      if (hodge1[2]>hodge2[2]) return 1;
      if (hodge1[1]<hodge2[1]) return -1;
      if (hodge1[1]>hodge2[1]) return 1;
      if (hodge1[0]<hodge2[0]) return -1;
      if (hodge1[0]>hodge2[0]) return 1;
      return 0;
}

spectrum search={0,0,0};
void searchspec(spectrum);
void addhod(spectrum hodge){
      sint i, j, k;
      searchspec(hodge);
      for(j=0; ((hodcomp(hodge, hodlist[j])>0)&&(j<specnum)); ++j);
      if ((j==specnum)||hodcomp(hodge, hodlist[j])){
       for(k=specnum;k>j;k--)
	for (i=0;i<3;i++) hodlist[k][i]=hodlist[k-1][i];
       for (i=0;i<3;i++) hodlist[j][i]=hodge[i];
       specnum++;   }
}

void finishmodel()
{     sint j;
      if(LONGOUT)fprintf(outfi,"\nspecnum, symnum: ");
      fprintf(outfi,"sp=%ld sy=%ld\n", specnum, symnum);
      for (j=0;j<specnum;j++)
       {int B01=hodlist[j][2], NGB=hodlist[j][0]-2*B01,CHI=hodlist[j][1];
         if(B01) fprintf(outfi,"-%d %d %d %d\n",NGB-CHI/2,NGB,CHI,B01);
	   else    fprintf(outfi,"%d %d %d\n",NGB-CHI/2,NGB,CHI);}
      /*fprintf(outfi,"%d %d\n",hodlist[j][0],hodlist[j][1]);*//* b01 change */
      totsymnum+=symnum;
      totspecnum+=specnum;
      maxspecnum=max(maxspecnum,specnum);
      maxsymnum=max(maxsymnum,symnum);
      modelnum++;
      if (!specnum) bugcount++;
}
void ErrEx(char *c){puts(c);exit(0);}
void ReadEOL(){char c;
  while('\n'!=(c=fgetc(infi)))if(c==EOF){puts("End of File");exit(0);}}

void ReadSpec(){int g,a,c,b=0;       /* g=h[0]-(h[1]=chi)/2; a=h[0]; b=h[3]; */
  if(stdi)printf("Type 'g a c' or '-g a c h01' with g=h12 and a=h11: ");
  fscanf(infi,"%d%d%d",&g,&a,&c); if(g<0){ fscanf(infi,"%d",&b); g=-g;}
  if(c!=2*(a-g))ErrEx("inconsistent a g c [need c=2(a-g)]"); 
  if(b&&c)ErrEx("inconsistent: b>0 => c=0");        /* search={a,c=2(a-g),b} */
  ReadEOL(); search[0]=a; search[1]=c; search[2]=b;}
void PrintUse(char *s){puts(s);
  puts("Either '-s' or '-g # -a #' is required, the rest is optional");
  puts(" -s      ask for spectrum");
  puts(" -g h12  #generations");    puts(" -a h11  #anti-generations");
  puts(" -b h01  #(01)-forms");     puts(" -t #    # of trivial pairs [1]");
  puts(" -i InFile "); puts(" -o OutFile"); exit(0);}
void LgoTwistInit(int narg, char* fn[]){ int n=1,t=0, /* t-> # trivial pairs */
  g=0,a=0,b=0; char *c;  infi=stdin; stdi=1; outfi=stdout;
  if(narg<2)PrintUse(""); while(n<narg) if(fn[n][0]!='-') PrintUse(""); else 
  switch(fn[n][1]){
    case 's': ReadSpec(); n++; break;
    case 't': c=(fn[n][2]) ? &fn[n][2] : fn[++n]; t=atoi(c); n++; break;
    case 'g': c=(fn[n][2]) ? &fn[n][2] : fn[++n]; if((*c)=='-')PrintUse("");
	g=atoi(c); n++; break;
    case 'a': c=(fn[n][2]) ? &fn[n][2] : fn[++n]; 
	a=atoi(c); n++; break;
    case 'b': c=(fn[n][2]) ? &fn[n][2] : fn[++n]; 
	b=atoi(c); n++; break;
    case 'i': c=(fn[n][2]) ? &fn[n][2] : fn[++n]; infi=fopen(c,"r"); 
	if(infi==NULL) PrintUse("Open infile failed"); n++; stdi=0; break; 
    case 'o': c=(fn[n][2]) ? &fn[n][2] : fn[++n]; outfi=fopen(c,"w"); 
	if(outfi==NULL) PrintUse("Open outfile failed"); n++; break; 
    default: PrintUse("Wrong option");}
  if(g+a) {search[0]=a; search[1]=2*(a-g); search[2]=b;}
  else {a=search[0]; g=search[0]-search[1]/2; b=search[2]; }  addsyms=2*t; 
  printf("searching for '%d %d %d", b ? -g : g,a,2*(a-g));
  if(b)printf(" %d",b);
  printf("' with %d trivial pair(s); stop skeleton input with EOF\n",t); 
  }

int  main(int narg, char* fn[])
{    skelet s;
     LgoTwistInit(narg,fn);
/*   if (narg>1) infi=fopen(fn[1],"r");		aao2.6	  002244 4 3 4 3 4 3
     else {
       infi=stdin; stdi=1;
       printf("usage: arg1=input file [stdin]; arg2=output file [stdout];\n");
       printf("       arg3=#pairs of trivial fields to be added;\n");
       printf("skeleton = string[#] exp_1 ... exp_#\n");}
     if (narg>2) outfi=fopen(fn[2],"w"); else outfi=stdout;
     if (narg>3) addsyms= *fn[3]-'0'; else addsyms=0;
     if (narg>4) search[1]=atoi(fn[4]); if (narg>5) search[0]=atoi(fn[5]);
     if (narg>6) search[2]=atoi(fn[6]); search[1]=2*(search[0]-search[1]);
 */  while(readline(&s)&&!bugcount){               /* search={a,c=2(a-g),b} */
       analy(s);   if(invertible || !ONLYINV) {
	  datain();
	  if(LONGOUT)longoutput0();
	  primedecomp();
	  if(LONGOUT)longoutput1();
	  gooddets();
	  if(LONGOUT)longoutput2();
	  reccon(0,0);
	  finishmodel();}      }
/*  fprintf(outfi,"modelnum, totsymnum, maxsymnum, totspecnum, maxspecnum: ");
    fprintf(outfi,"%ld %ld %ld %ld %ld", modelnum, totsymnum, maxsymnum,
	     totspecnum, maxspecnum); */
    fprintf(outfi,
       "#skel=%ld, #sym=%ld, sym/skel<=%ld, #spectra=%ld, spec/skel<=%ld",
      modelnum,  totsymnum, maxsymnum,    totspecnum,   maxspecnum);
    printpri(pmax);
}

/*  =========================         b01         ========================= */
int  interb01(/*struct ein */); /* ...  interface to "mhodge" conventions   */
/*  calculate b01 for NG generators: Z_{*gen[j]}: (gen[j][1],...,gen[j][N]) *
 *  0<=j<NG on N fields X_i. zd[i]/(*zd) is the weight of X_i,0<i<=N.       *
 *          !!! note that gen is changed (in version without rat) !!!       */
/*  how to avoid the recursion: pow[j] determines the current group element,*
 *  i.e. we need all possible vectors with  0<=pow[j]<O_j  for  0<=j<J:     *
 *  logic: if pow_j==O_j ==> "exit" (i.e. pow[j++]++); else                 *
 *  if j==0 ==> "do it"; else "call" (i.e. {pow[--j]=0; initializations;})  *
 *  code:   pow[j=NG-1]=0;while(j<NG){if(pow[j]==O_j)pow[++j]++; else       *
 *          if(j)pow[--j]=0; else {"do it for pow[]!";  *pow++;}}         *
 *  where O_j=*gen[j], J=NG;      !!! note that pow[J+1] is set to 0 !!!    */
#if SAFER
void mod1(rat *a) {a->num=mod(a->num,a->den);}
int b01(sint N, sint NG, sint zd[NM+1], lint gen[NS][NM+1])
{   int b=0, i,j, eq[NM+1],pow[NS+1]; /* eq=th_i==q_i; pow(er) of generator */
over=0;
    pow[j=NG-1]=0; while(j<NG) {      /*     -- begin generate group --     */
    if(pow[j]==*gen[j]) pow[++j]++; else if(j) pow[--j]=0; else { /* do it: */

int can=1, sum=0, k;        /* can(didate) for b01-contribution */
       for(i=1;can&&(i<=N);i++) {rat th=rI(0);    /* check th=0/q -> set eq */
        for(k=0;k<NG;k++){th=rS(th,rR(pow[k]*gen[k][i],*gen[k]));mod1(&th);}
        if(th.num) {if(th.num<0) th.num+=th.den;
                     th=rP(th,rR(*zd,zd[i])); if(th.num-th.den) can=0;
                 else {eq[i]=1;sum+=*zd-2*zd[i];} }    else eq[i]=0; }
if(sum&&(sum!= *zd)) can=0;                                /* i.e. charge=1? */
       for(k=0;can&&(k<NG);k++){lint gph=0;/*check inv. under G_k: gph(ase)*/
         for(i=1;i<=N;i++) {if(eq[i]) gph+=gen[k][i];gph=mod(gph,*gen[k]);}
	 if(gph) can=0; }
if(can) {if(sum) b++; else over++;}

    (*pow)++;}}                                    /* end of generate group */
    return b;
}
#else
int b01(sint N, sint NG, sint zd[NM+1], lint gen[NS][NM+1])
{   int b=0, i,j, eq[NM+1],pow[NS+1]; /* eq=th_i==q_i; pow(er) of generator */
    lint G=**gen; for(j=1;j<NG;j++) G*=(*gen[j]/gcd(G,*gen[j]));/* lcm(O_j) */
    {lint inv; for(j=0;j<NG;j++) {inv=G/(*gen[j]);             /* inv=G/O_j */
	for(i=1;i<=N;i++) gen[j][i]*=inv;}} /* bring phases to common den.*/
    pow[j=NG-1]=0; while(j<NG) {      /* begin generate group (+ next line) */
    if(pow[j]==*gen[j]) pow[++j]++; else if(j) pow[--j]=0; else { /* do it: */

       int can=1, sum=-(*zd), k;        /* can(didate) for b01-contribution */
       for(i=1;can&&(i<=N);i++)                   /* check th=0/q -> set eq */
         {lint th=0; for(k=0;k<NG;k++) th+=pow[k]*gen[k][i]; th=mod(th,G);
	  if(th<0){printf("warning: th<0!\n"); th+=G;}
          if(th) if(th-zd[i]*(G/ *zd)) can=0;           /* no candidate */
	       else {eq[i]=1;sum+=*zd-2*zd[i];} else eq[i]=0; }
       if(sum) can=0;                                /* i.e. charge=1? */
       for(k=0;can&&(k<NG);k++){lint gph=0;/*check inv. under G_k: gph(ase)*/
         for(i=1;i<=N;i++) {if(eq[i]) gph+=gen[k][i];gph=mod(gph,G);}
         if(gph) can=0; }
       if(can) b++;
    (*pow)++;}}                                    /* end of generate group */
    return b;
}
#endif
int interb01()    /* aao.c:  *e.np -> ns;  *e.n -> n;  *e.w -> wei;  PM->NS */
{   sint i,j,zd[NM+1]; lint gen[NS][NM+1]; zd[0]=wei[n][0];
    for(i=1-ADDZD;i<=ns;i++) gen[i][0]=wei[n][i];
    for(i=1;i<n+1;i++) {zd[i]=wei[i-1][0];
       for(j=1-ADDZD;j<=ns;j++) gen[j][i]=wei[i-1][j]; }
    return b01(n, ns+ADDZD, zd, &(gen[1-ADDZD]));
}
/*  =========================     end of b01      ========================= */

/*  ========================    search for g,a,b  ========================= */
#define   fpri(list,num)            fprintf(outfi,list,(long) num)
void searchspec(spectrum h)         /* g=h[0]-(h[1]=chi)/2; a=h[0]; b=h[3]; */
{   /* printf("h=%d %d %d  search=%d %d %d\n",h[0],h[1],h[2],
	search[0],search[1],search[2]);*/
    if(h[0]!=search[0]+2*search[2]) return;  /* h[0]=trace=s[0]+2*s[2] !!! */
     if(h[1]!=search[1]) return;
     if(h[2]!=search[2]) return;
#ifdef	OLD_FORMAT
     {int i,j; fprintf(outfi,"\n"); fpri("C_{(%ld",wei[0][0]);
	 for(i=1;i<n;i++) fpri(",%ld",wei[i][0]); fpri(")[%ld]\n",wei[n][0]);
	 for(j=1;j<=ns;j++) {fpri("Z_{%ld}(",wei[n][j]);
	    for(i=0;i<n-1;i++) fpri("%ld,",wei[i][j]);
	    fpri("%ld)\n",wei[i][j]); }}
#else
     {int i,j; fprintf(outfi,"\n"); fpri("%ld",wei[n][0]);
	 for(i=0;i<n;i++) fpri(" %ld",wei[i][0]);
	if(n%2==0) fpri("%ld",wei[n][0]/2);
	 for(j=1;j<=ns;j++) {fpri(" /Z%ld:",wei[n][j]);
	    for(i=0;i<n;i++) fpri(" %ld",wei[i][j]);}}
#endif
/*     exit(); */
}
