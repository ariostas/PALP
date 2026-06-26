#pragma once

/*
 * C++-friendly type aliases for the integer types used throughout PALP.
 * These mirror the traditional C macros:
 *     #define Long  long
 *     #define LLong long long
 * but are safe in C++ where the macros could otherwise interfere with
 * standard-library names or template parameters.
 */

using Long = long;
using LLong = long long;
