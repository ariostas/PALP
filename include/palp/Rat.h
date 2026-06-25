#ifndef __Rat__
#define __Rat__

struct Rat { Long N; Long D; Rat(Long n=0, Long d=1) : N(n), D(d) {} };
struct LRat { LLong N; LLong D; LRat(LLong n=0, LLong d=1) : N(n), D(d) {} };

Long Fgcd(Long a, Long b);		   /* Fast greatest common divisor  */
Long NNgcd(Long a, Long b); 		   /* NonNegative gcd handling zero */

Long Egcd(Long, Long, Long*,Long*);		     /*  extended gcd(a,b)  */
Long REgcd(Long *vec_in, int *d, Long *vec_out);/* extended gcd(a_1,...a_d) */

LLong LFgcd(LLong a, LLong b);		   /* Fast greatest common divisor  */
LLong LNNgcd(LLong a, LLong b); 	   /* NonNegative gcd handling zero */

LLong LEgcd(LLong, LLong, LLong*,LLong*);	     /*  extended gcd(a,b)  */
LLong LREgcd(LLong *vec_in, int *d, LLong *vec_out);
                                                /* extended gcd(a_1,...a_d) */

Rat  rI(Long a);		/*  conversion  Long -> Rat  */
Rat  rR(Long a, Long b);	/*  conversion  a/b  -> Rat  */
Rat  irP(Long a, Rat b);	/*  a b		integer * Rat */
Rat  rS(Rat a, Rat b);		/*  a + b	rational Sum */
Rat  rD(Rat a, Rat b);		/*  a - b	rational Difference */
Rat  rP(Rat a, Rat b);		/*  a * b	rational Product   */
Rat  rQ(Rat a, Rat b);    	/*  a / b	rational Quotient */
int  rC(Rat a, Rat b);          /* Compare = [1 / 0 / -1] if a [gt/eq/lt] b */
void Rpr(Rat c);		/*  write  "c.N/c.D"  to outFN */

LRat  LrI(LLong a);		/*  conversion  LLong -> LRat  	*/
LRat  LrR(LLong a, LLong b);	/*  conversion  a/b   -> LRat  	*/
LRat  LirP(LLong a, LRat b);	/*  a b		integer * LRat 		*/
LRat  LrS(LRat a, LRat b);	/*  a + b	LRational Sum 		*/
LRat  LrD(LRat a, LRat b);	/*  a - b	LRational Difference	*/
LRat  LrP(LRat a, LRat b);	/*  a * b	LRational Product   	*/
LRat  LrQ(LRat a, LRat b);    	/*  a / b	LRational Quotient 	*/
int   LrC(LRat a, LRat b);      /* Compare = [1 / 0 / -1] if a [gt/eq/lt] b */
void  LRpr(LRat c);		/*  write  "c.N/c.D"  to outFN */

/*   Map Permutations: Do "ArgFun" for all permutations pi of *d elements */
#define ARG_FUN		void (*ArgFun)(int *d,int *pi,int *pinv,void *info)
void  Map_Permut(int *d,int *pi,int *pinv,ARG_FUN,void *AuxPtr);

Long  W_to_GLZ(Long *W, int *d, Long **GLZ);	/* "triangluar" form of GLZ */
Long  PW_to_GLZ(Long *W, int *d, Long **GLZ);	/* improved by permutations */

#endif
