#pragma once

/*
 * C++-friendly type aliases for the integer types used throughout PALP.
 * These mirror the traditional C macros:
 *     #define Long  long
 *     #define LLong long long
 * but are safe in C++ where the macros could otherwise interfere with
 * standard-library names or template parameters.
 */

#include <cstdint>

using Long = long;
using LLong = long long;

static_assert(sizeof(Long) == 8, "Long must be 64-bit");
static_assert(sizeof(LLong) == 8, "LLong must be 64-bit");
