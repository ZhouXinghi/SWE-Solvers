/***********************************************************************************/ /**
                                                                                       *
                                                                                       * \file SIMD_DEFINITIONS.hpp
                                                                                       *
                                                                                       * \brief Contains macro
                                                                                       *definitions for the intrinsics
                                                                                       *used for the vectorization
                                                                                       *
                                                                                       * \author Wolfgang Hölzl
                                                                                       *(hoelzlw), hoelzlw AT in.tum.de
                                                                                       *
                                                                                       **************************************************************************************/

#pragma once

#include <cmath>
#include <limits>
#include <immintrin.h>


/*
 * Map single precision AVX intrinsics
 */
#define ADDV _mm256_add_pd
#define SUBV _mm256_sub_pd
#define MULV _mm256_mul_pd
#define DIVV _mm256_div_pd
#define SQRTV _mm256_sqrt_pd

#define LOADU _mm256_loadu_pd
#define STOREU _mm256_storeu_pd
#define SETV_R _mm256_set1_pd
#define SETV_I _mm256_set1_epi64x
#define ZEROV_R _mm256_setzero_pd
#define ZEROV_I _mm256_setzero_si256

#define MAXV _mm256_max_pd
#define MINV _mm256_min_pd

#define CMP_LT(x, y) _mm256_cmp_pd((x), (y), _CMP_LT_OS)
#define CMP_LE(x, y) _mm256_cmp_pd((x), (y), _CMP_LE_OS)
#define CMP_GT(x, y) _mm256_cmp_pd((x), (y), _CMP_GT_OS)
#define CMP_GE(x, y) _mm256_cmp_pd((x), (y), _CMP_GE_OS)

#define CMP_EQ_R(x, y) _mm256_cmp_pd((x), (y), _CMP_EQ_OS)

#define CMP_EQ_I(x, y) _mm256_cmpeq_epi64((x), (y))

#define ANDV_R _mm256_and_pd
#define ORV_R _mm256_or_pd
#define XORV_R _mm256_xor_pd
#define ORV_I(x, y) CAST_REAL_TO_INT_V(_mm256_or_pd(CAST_INT_TO_REAL_V(x), CAST_INT_TO_REAL_V(y)))

#define NOTV_R not_pd
#define NOTV_I not_si256
#define ANDNOTV_R _mm256_andnot_pd

#define BLENDV _mm256_blendv_pd
#define BLENDV_I(else_part, if_part, mask) \
  CAST_REAL_TO_INT_V(_mm256_blendv_pd(CAST_INT_TO_REAL_V(else_part), CAST_INT_TO_REAL_V(if_part), mask))
#define MOVEMASK _mm256_movemask_pd

#define SHIFT_LEFT(x, y) _mm256_slli_epi32((x), (y))

#define CAST_INT_TO_REAL_V _mm256_castsi256_pd
#define CAST_REAL_TO_INT_V _mm256_castpd_si256
#define FABS fabs_pd

/*
 * Compute the absolute value of a vector by forcing the sign bit to be zero
 */
inline __m256d fabs_pd(const __m256d x) {
//   static const __m256d sign_mask = CAST_INT_TO_REAL_V(_mm256_set1_epi64x(1 << 63));
  const __m256d sign_mask = _mm256_set1_pd(-0.0);
  return _mm256_andnot_pd(sign_mask, x);
}

/*
 * Bitwise NOT operation for integers
 */
inline __m256i not_si256(const __m256i x) {
  static const __m256i mask = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFF);
  return CAST_REAL_TO_INT_V(_mm256_xor_pd(CAST_INT_TO_REAL_V(mask), CAST_INT_TO_REAL_V(x)));
}

/*
 * Bitwise NOT operation for reals
 */
inline __m256d not_pd(const __m256d x) {
  static const __m256i mask = _mm256_set1_epi64x(0xFFFFFFFFFFFFFFFF);
  return _mm256_xor_pd(CAST_INT_TO_REAL_V(mask), x);
}

/*
 * Check, whether a real_vector contains infinity or NaN
 */
inline bool checkVector(const __m256d x) {
  static const real_vector infinity = SETV_R(std ::numeric_limits<float>::infinity());
  return MOVEMASK(ANDV_R(CMP_EQ_R(x, x), NOTV_R(CMP_EQ_R(x, infinity)))) == VECTOR_FULL_MASK;
}

