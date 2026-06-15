#include "Global.h"
#include "Rat.h"
#include "LG.h"
#include "Subpoly.h"
#include "Nef.h"
#include "Mori.h"

#include "Global.h"
#include "Rat.h"
#include "LG.h"
#include "Subpoly.h"
#include "Nef.h"
#include "Mori.h"

static_assert(sizeof(Long) * CHAR_BIT >= 32,
              "Long must provide at least 32 bits");
static_assert(sizeof(LLong) * CHAR_BIT >= 64,
              "LLong must provide at least 64 bits");

int main()
{
  return sizeof(PolyPointList) == 0;
}
