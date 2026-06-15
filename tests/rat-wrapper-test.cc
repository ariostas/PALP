#include "Global.h"
#include "Rat.h"

#include <cassert>
#include <stdexcept>

extern "C" {
FILE *inFILE = stdin;
FILE *outFILE = stdout;
}

int main()
{
  using palp::LongRational;
  using palp::Rational;

  const Rational half(1, 2);
  const Rational third(1, 3);
  const Rational minus_half(-2, 4);

  assert(half.numerator() == 1);
  assert(half.denominator() == 2);
  assert(minus_half.numerator() == -1);
  assert(minus_half.denominator() == 2);
  assert((half + third) == Rational(5, 6));
  assert((half - third) == Rational(1, 6));
  assert((half * Rational(6, 5)) == Rational(3, 5));
  assert((half / third) == Rational(3, 2));
  assert(third < half);
  assert(Rational(2) == Rational(2, 1));
  const Rat legacy_half = half.legacy();
  assert(legacy_half.N == 1);
  assert(legacy_half.D == 2);
  assert(Rational(rR(2, 4)) == half);

  bool threw = false;
  try {
    (void)Rational(1, 0);
  } catch (const std::domain_error &) {
    threw = true;
  }
  assert(threw);

  threw = false;
  try {
    (void)(half / Rational(0));
  } catch (const std::domain_error &) {
    threw = true;
  }
  assert(threw);

  const LongRational large(6000000000LL, 4000000000LL);
  assert(large == LongRational(3, 2));
  assert((large + LongRational(1, 2)) == LongRational(2));
  const LRat legacy_large = large.legacy();
  assert(legacy_large.N == 3);
  assert(legacy_large.D == 2);
  assert(LongRational(LrR(6, 4)) == large);

  threw = false;
  try {
    (void)LongRational(5, 0);
  } catch (const std::domain_error &) {
    threw = true;
  }
  assert(threw);

  return 0;
}
