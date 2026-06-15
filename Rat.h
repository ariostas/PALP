#ifndef PALP_RAT_H
#define PALP_RAT_H

#ifdef __cplusplus
extern "C" {
#endif

typedef	struct {Long N; Long D;} 		             Rat;  /* = N/D */

Long Fgcd(Long a, Long b);		   /* Fast greatest common divisor  */
Long NNgcd(Long a, Long b); 		   /* NonNegative gcd handling zero */

Long Egcd(Long, Long, Long*,Long*);		     /*  extended gcd(a,b)  */
Long REgcd(Long *vec_in, int *d, Long *vec_out);/* extended gcd(a_1,...a_d) */

Rat  rI(Long a);		/*  conversion  Long -> Rat  */
Rat  rR(Long a, Long b);	/*  conversion  a/b  -> Rat  */
Rat  irP(Long a, Rat b);	/*  a b		integer * Rat */
Rat  rS(Rat a, Rat b);		/*  a + b	rational Sum */
Rat  rD(Rat a, Rat b);		/*  a - b	rational Difference */
Rat  rP(Rat a, Rat b);		/*  a * b	rational Product   */
Rat  rQ(Rat a, Rat b);    	/*  a / b	rational Quotient */
int  rC(Rat a, Rat b);          /* Compare = [1 / 0 / -1] if a [gt/eq/lt] b */
void Rpr(Rat c);		/*  write  "c.N/c.D"  to outFN */

typedef	struct {LLong N; LLong D;} 		             LRat;  /* = N/D */

LLong LFgcd(LLong a, LLong b);		   /* Fast greatest common divisor  */
LLong LNNgcd(LLong a, LLong b); 	   /* NonNegative gcd handling zero */

LLong LEgcd(LLong, LLong, LLong*,LLong*);	     /*  extended gcd(a,b)  */
LLong LREgcd(LLong *vec_in, int *d, LLong *vec_out);
                                                /* extended gcd(a_1,...a_d) */

LRat  LrI(LLong a);		/*  conversion  LLong -> LRat  		*/
LRat  LrR(LLong a, LLong b);	/*  conversion  a/b   -> LRat  		*/
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

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

#include <stdexcept>

namespace palp {

class Rational {
public:
  Rational() : value_(rI(0)) {}
  explicit Rational(Long value) : value_(rI(value)) {}

  Rational(Long numerator, Long denominator)
  {
    if (denominator == 0) {
      throw std::domain_error("Rational denominator must be nonzero");
    }
    value_ = rR(numerator, denominator);
  }

  explicit Rational(Rat value) : Rational(value.N, value.D) {}

  Long numerator() const noexcept { return value_.N; }
  Long denominator() const noexcept { return value_.D; }
  Rat legacy() const noexcept { return value_; }

  friend Rational operator+(Rational lhs, Rational rhs)
  {
    return from_legacy(rS(lhs.value_, rhs.value_));
  }

  friend Rational operator-(Rational lhs, Rational rhs)
  {
    return from_legacy(rD(lhs.value_, rhs.value_));
  }

  friend Rational operator*(Rational lhs, Rational rhs)
  {
    return from_legacy(rP(lhs.value_, rhs.value_));
  }

  friend Rational operator/(Rational lhs, Rational rhs)
  {
    if (rhs.value_.N == 0) {
      throw std::domain_error("Rational division by zero");
    }
    return from_legacy(rQ(lhs.value_, rhs.value_));
  }

  friend bool operator==(Rational lhs, Rational rhs)
  {
    return rC(lhs.value_, rhs.value_) == 0;
  }

  friend bool operator!=(Rational lhs, Rational rhs) { return !(lhs == rhs); }
  friend bool operator<(Rational lhs, Rational rhs)
  {
    return rC(lhs.value_, rhs.value_) < 0;
  }
  friend bool operator>(Rational lhs, Rational rhs) { return rhs < lhs; }
  friend bool operator<=(Rational lhs, Rational rhs) { return !(rhs < lhs); }
  friend bool operator>=(Rational lhs, Rational rhs) { return !(lhs < rhs); }

private:
  static Rational from_legacy(Rat value)
  {
    Rational result;
    result.value_ = value;
    return result;
  }

  Rat value_;
};

class LongRational {
public:
  LongRational() : value_(LrI(0)) {}
  explicit LongRational(LLong value) : value_(LrI(value)) {}

  LongRational(LLong numerator, LLong denominator)
  {
    if (denominator == 0) {
      throw std::domain_error("LongRational denominator must be nonzero");
    }
    value_ = LrR(numerator, denominator);
  }

  explicit LongRational(LRat value) : LongRational(value.N, value.D) {}

  LLong numerator() const noexcept { return value_.N; }
  LLong denominator() const noexcept { return value_.D; }
  LRat legacy() const noexcept { return value_; }

  friend LongRational operator+(LongRational lhs, LongRational rhs)
  {
    return from_legacy(LrS(lhs.value_, rhs.value_));
  }

  friend LongRational operator-(LongRational lhs, LongRational rhs)
  {
    return from_legacy(LrD(lhs.value_, rhs.value_));
  }

  friend LongRational operator*(LongRational lhs, LongRational rhs)
  {
    return from_legacy(LrP(lhs.value_, rhs.value_));
  }

  friend LongRational operator/(LongRational lhs, LongRational rhs)
  {
    if (rhs.value_.N == 0) {
      throw std::domain_error("LongRational division by zero");
    }
    return from_legacy(LrQ(lhs.value_, rhs.value_));
  }

  friend bool operator==(LongRational lhs, LongRational rhs)
  {
    return LrC(lhs.value_, rhs.value_) == 0;
  }

  friend bool operator!=(LongRational lhs, LongRational rhs)
  {
    return !(lhs == rhs);
  }
  friend bool operator<(LongRational lhs, LongRational rhs)
  {
    return LrC(lhs.value_, rhs.value_) < 0;
  }
  friend bool operator>(LongRational lhs, LongRational rhs)
  {
    return rhs < lhs;
  }
  friend bool operator<=(LongRational lhs, LongRational rhs)
  {
    return !(rhs < lhs);
  }
  friend bool operator>=(LongRational lhs, LongRational rhs)
  {
    return !(lhs < rhs);
  }

private:
  static LongRational from_legacy(LRat value)
  {
    LongRational result;
    result.value_ = value;
    return result;
  }

  LRat value_;
};

}  // namespace palp

#endif

#endif
