#include "Global.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

FILE *inFILE;
FILE *outFILE;

/* The tested Read_PP/Print_PPL path never reaches quotient reduction. */
extern "C" void QuotZ_2_SublatG(Long[][VERT_Nmax], int *, Long *, int *,
				Long[POLY_Dmax][POLY_Dmax])
{
  std::fprintf(stderr, "unexpected QuotZ_2_SublatG call\n");
  std::exit(1);
}

static void require(bool condition, const char *message)
{
  if(!condition) {
    std::fprintf(stderr, "%s\n", message);
    std::exit(1);
  }
}

int main()
{
  inFILE = stdin;
  outFILE = stdout;
  FILE *input = std::tmpfile();
  FILE *output = std::tmpfile();
  require(input != nullptr, "tmpfile failed for input");
  require(output != nullptr, "tmpfile failed for output");

  std::fputs("2 3\n1 0 0\n0 1 0\n", input);
  std::rewind(input);

  PALP_RuntimeContext ctx;
  ctx.in = input;
  ctx.out = output;

  static PolyPointList poly;
  require(PALP_Read_PP(&ctx, &poly) == 1, "PALP_Read_PP failed");
  require(poly.n == 2, "unexpected polytope dimension");
  require(poly.np == 3, "unexpected point count");
  require(inFILE == stdin, "PALP_Read_PP did not restore inFILE");
  require(outFILE == stdout, "PALP_Read_PP did not restore outFILE");

  FILE *cws_input = std::tmpfile();
  require(cws_input != nullptr, "tmpfile failed for CWS input");
  std::fputs("4 1 1 1 1\n", cws_input);
  std::rewind(cws_input);
  ctx.in = cws_input;
  static CWS cws;
  static PolyPointList cws_poly;
  require(PALP_ReadCwsPp(&ctx, &cws, &cws_poly, 1, 1) == 1,
	  "PALP_ReadCwsPp failed");
  require(cws.nw == 1, "unexpected CWS weight-system count");
  require(cws_poly.n == 3, "unexpected CWS polytope dimension");
  require(inFILE == stdin, "PALP_ReadCwsPp did not restore inFILE");
  require(outFILE == stdout, "PALP_ReadCwsPp did not restore outFILE");

  ctx.in = input;
  PALP_Print_PPL(&ctx, &poly, "adapter");
  require(inFILE == stdin, "PALP_Print_PPL did not restore inFILE");
  require(outFILE == stdout, "PALP_Print_PPL did not restore outFILE");

  cws.nw = 1;
  cws.nz = 1;
  cws.N = 3;
  cws.m[0] = 5;
  cws.z[0][0] = 4;
  cws.z[0][1] = 1;
  cws.z[0][2] = 0;
  PALP_Print_CWS_Zinfo(&ctx, &cws);
  require(inFILE == stdin, "PALP_Print_CWS_Zinfo did not restore inFILE");
  require(outFILE == stdout, "PALP_Print_CWS_Zinfo did not restore outFILE");

  std::rewind(output);
  char buffer[256];
  const std::size_t nread = std::fread(buffer, 1, sizeof(buffer) - 1, output);
  buffer[nread] = '\0';
  require(std::strstr(buffer, "2 3  adapter\n") != nullptr,
	  "PALP_Print_PPL wrote unexpected output");
  require(std::strstr(buffer, "/Z5: 4 1 0 ") != nullptr,
	  "PALP_Print_CWS_Zinfo wrote unexpected output");

  std::fclose(input);
  std::fclose(cws_input);
  std::fclose(output);
  return 0;
}
