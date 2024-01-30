
#pragma once

#pragma message "Using double as type for real numbers"
typedef double real;
#pragma message "Using unsigned long long as type for integer numbers"
typedef unsigned long long integer;

#define SHIFT_SIGN_RIGHT(x) \
  static_cast<integer>(static_cast<integer>(1) << (static_cast<integer>(64) - static_cast<integer>(x)))

/*
 * Set control macros
 *
 * Declare the vector length
 * Additionally declare a configuration specific macro of the form
 * 	VECTOR_extension_precision
 *
 * Moreover, declare an integer, representing the number,
 * which is returned by the instruction MOVEMASK, if the instruction is called on a vector with ALL SIGN BITS set
 * Use is as
 *
 * 	const real_vector vector = CONDITION(operand_a, operand_b);
 *
 * 	if (MOVEMASK(vector) == VECTOR_FULL_MASK) {
 *		// all components fullfill the condition
 *	} else {
 *		// some components do not fullfill the condition
 *	}
 */
#pragma message "Using AVX for vectorization"
#include <immintrin.h>

typedef __m256i integer_vector;

#define VECTOR_LENGTH 4
#define VECTOR_AVX_FLOAT64
#define VECTOR_FULL_MASK 0x0000000F

typedef __m256d real_vector;
#pragma message "Using vectors of 4 doubles"

/*
 * Control macros are set
 *
 * Include the function macros
 */
#include "SIMDDefinitionsDouble.hpp"


