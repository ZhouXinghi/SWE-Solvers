
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "SIMDTypesDouble.hpp"
#include "WavePropagationSolver.hpp"

namespace Solvers {
  class AugRie_SIMD;
} // namespace solver
/**
 * Approximate Augmented Riemann Solver for the Shallow Water Equations.
 *
 * T should be double or float.
 */
class Solvers::AugRie_SIMD: public Solvers::WavePropagationSolver<real> {
private:
  // explicit for unit tests
  // use nondependent names (template base class)
  using WavePropagationSolver<real>::dryTol_;
  using WavePropagationSolver<real>::gravity_;
  using WavePropagationSolver<real>::zeroTol_;

  using WavePropagationSolver<real>::hLeft_;
  using WavePropagationSolver<real>::hRight_;
  using WavePropagationSolver<real>::huLeft_;
  using WavePropagationSolver<real>::huRight_;
  using WavePropagationSolver<real>::bLeft_;
  using WavePropagationSolver<real>::bRight_;
  using WavePropagationSolver<real>::uLeft_;
  using WavePropagationSolver<real>::uRight_;

  using Solvers::WavePropagationSolver<real>::storeParameters;

  integer              wetDryState;
  static const integer DryDry               = 0;
  static const integer WetWet               = SHIFT_SIGN_RIGHT(1);
  static const integer WetDryInundation     = SHIFT_SIGN_RIGHT(2);
  static const integer WetDryWall           = SHIFT_SIGN_RIGHT(3);
  static const integer WetDryWallInundation = SHIFT_SIGN_RIGHT(4);
  static const integer DryWetInundation     = SHIFT_SIGN_RIGHT(5);
  static const integer DryWetWall           = SHIFT_SIGN_RIGHT(6);
  static const integer DryWetWallInundation = SHIFT_SIGN_RIGHT(7);

#ifndef VECTOR_NOVEC
  // Vector components
  const integer_vector DryDry_V;
  const integer_vector WetWet_V;
  const integer_vector DryWetInundation_V;
  const integer_vector WetDryInundation_V;
  const integer_vector DryWetWallInundation_V;
  const integer_vector DryWetWall_V;
  const integer_vector WetDryWallInundation_V;
  const integer_vector WetDryWall_V;
  integer_vector       wetDryState_v;
#endif /* VECTOR_NOVEC */

  //! Newton-tolerance (exit Newton-Raphson-method, if we are close enough to the root)
  const real newtonTol;
  //! maximum number of Newton-Raphson-Iterations
  const int maxNumberOfNewtonIterations;

#ifndef VECTOR_NOVEC
  // Vector members
  real_vector hLeft_v;
  real_vector hRight_v;
  real_vector huLeft_v;
  real_vector huRight_v;
  real_vector bLeft_v;
  real_vector bRight_v;
  real_vector uLeft_v;
  real_vector uRight_v;
#endif /* VECTOR_NOVEC */

  //! height of our homogeneous Riemann-problem at middle state (computed by determineMiddleState)
  real hMiddle;
#ifndef VECTOR_NOVEC
  real_vector hMiddle_v;
#endif /* VECTOR_NOVEC */

  //! shock or inner rarefaction speeds of our homogeneous Riemann-problem (computed by computeMiddleState)
  real middleStateSpeeds[2];
#ifndef VECTOR_NOVEC
  real_vector middleStateSpeeds_v[2];
#endif /* VECTOR_NOVEC */

  /**
   * the Riemann-struture of the homogeneous Riemann-problem.
   */
  static const integer ShockShock             = SHIFT_SIGN_RIGHT(1);
  static const integer RarefactionRarefaction = SHIFT_SIGN_RIGHT(2);
  static const integer ShockRarefaction       = SHIFT_SIGN_RIGHT(3);
  static const integer RarefactionShock       = SHIFT_SIGN_RIGHT(4);
  static const integer DrySingleRarefaction   = SHIFT_SIGN_RIGHT(5);
  static const integer SingleRarefactionDry   = SHIFT_SIGN_RIGHT(6);

#ifndef VECTOR_NOVEC
  integer_vector       riemannStructure_v;
  const integer_vector ShockShock_V;
  const integer_vector RarefactionRarefaction_V;
  const integer_vector ShockRarefaction_V;
  const integer_vector RarefactionShock_V;
  const integer_vector DrySingleRarefaction_V;
  const integer_vector SingleRarefactionDry_V;
#endif /* VECTOR_NOVEC */

  //! Riemann-structure of our homogeneous Riemann-problem (determined by determineRiemannStructure)
  integer riemannStructure;

public:
  /**
   * Constructor of the Augmented Riemann solver with optional parameters.
   *
   * @param i_dryTolerance numerical definition of "dry".
   * @param i_gravity gravity constant.
   * @param i_newtonTolerance numerical definition of "convergence" (used in the AugRie solver only).
   * @param i_maxNumberOfNewtonIterations maximum steps for the Newton-Raphson method (used in the AugRie solver only).
   * @param i_zeroTolerance numerical definition of zero.
   */
  AugRie_SIMD(
    real i_dryTolerance                = static_cast<real>(0.01),
    real i_gravity                     = static_cast<real>(9.81),
    real i_newtonTolerance             = static_cast<real>(0.000001),
    int  i_maxNumberOfNewtonIterations = 10,
    real i_zeroTolerance               = static_cast<real>(0.00001)
  ):
    WavePropagationSolver<real>(i_dryTolerance, i_gravity, i_zeroTolerance),

    DryDry_V(SETV_I(DryDry)),
    WetWet_V(SETV_I(WetWet)),
    DryWetInundation_V(SETV_I(DryWetInundation)),
    WetDryInundation_V(SETV_I(WetDryInundation)),
    DryWetWallInundation_V(SETV_I(DryWetWallInundation)),
    DryWetWall_V(SETV_I(DryWetWall)),
    WetDryWallInundation_V(SETV_I(WetDryWallInundation)),
    WetDryWall_V(SETV_I(WetDryWall)),

    newtonTol(i_newtonTolerance),
    maxNumberOfNewtonIterations(i_maxNumberOfNewtonIterations),

    ShockShock_V(SETV_I(ShockShock)),
    RarefactionRarefaction_V(SETV_I(RarefactionRarefaction)),
    ShockRarefaction_V(SETV_I(ShockRarefaction)),
    RarefactionShock_V(SETV_I(RarefactionShock)),
    DrySingleRarefaction_V(SETV_I(DrySingleRarefaction)),
    SingleRarefactionDry_V(SETV_I(SingleRarefactionDry)),

    sqrt_g(std ::sqrt(gravity_)),

    sqrt_g_v(SETV_R(sqrt_g)),
    zeroTol_v(SETV_R(zeroTol_)),
    zeroTol_v_neg(SETV_R(-zeroTol_)),
    dryTol_v(SETV_R(dryTol_)),
    dryTol_v_neg(SETV_R(-dryTol_)),
    g_v(SETV_R(gravity_))
  {
  }

  ~AugRie_SIMD() {}

  // declare variables which are used over and over again
  real sqrt_g_hLeft;
  real sqrt_g_hRight;

  real sqrt_hLeft;
  real sqrt_hRight;

  real_vector sqrt_g_hLeft_v;
  real_vector sqrt_g_hRight_v;

  real_vector sqrt_hLeft_v;
  real_vector sqrt_hRight_v;

  const real sqrt_g;
  const real_vector sqrt_g_v;

  const real_vector zeroTol_v;
  const real_vector zeroTol_v_neg;
  const real_vector dryTol_v;
  const real_vector dryTol_v_neg;
  const real_vector g_v;

public:
  void computeNetUpdates
  (
      const real& hLeft,
      const real& hRight,
      const real& huLeft,
      const real& huRight,
      const real& bLeft,
      const real& bRight,
      real&       o_hUpdateLeft,
      real&       o_hUpdateRight,
      real&       o_huUpdateLeft,
      real&       o_huUpdateRight,
      real&       o_maxWaveSpeed
  )
  {
    int i = 0;
  }

  void computeNetUpdates_SIMD(
    const real* const i_hLeft,
    const real* const i_hRight,
    const real* const i_huLeft,
    const real* const i_huRight,
    const real* const i_bLeft,
    const real* const i_bRight,

    real* const o_hUpdateLeft,
    real* const o_hUpdateRight,
    real* const o_huUpdateLeft,
    real* const o_huUpdateRight,
    real&       o_maxWaveSpeed
  )
  {
    hLeft_v   = LOADU(i_hLeft);
    hRight_v  = LOADU(i_hRight);
    huLeft_v  = LOADU(i_huLeft);
    huRight_v = LOADU(i_huRight);
    bLeft_v   = LOADU(i_bLeft);
    bRight_v  = LOADU(i_bRight);
    uLeft_v   = ZEROV_R();
    uRight_v  = ZEROV_R();

    real_vector hUpdateLeft_v   = ZEROV_R();
    real_vector hUpdateRight_v  = ZEROV_R();
    real_vector huUpdateLeft_v  = ZEROV_R();
    real_vector huUpdateRight_v = ZEROV_R();

    determineWetDryState();
    // all 4 is not drydry
    const real_vector wetDryMask = CAST_INT_TO_REAL_V(NOTV_I(CMP_EQ_I(wetDryState_v, DryDry_V)));

    // if all 4 is drydry
    if (MOVEMASK(wetDryMask) == 0)  {
      STOREU(o_hUpdateLeft, ZEROV_R());
      STOREU(o_huUpdateLeft, ZEROV_R());
      STOREU(o_hUpdateRight, ZEROV_R());
      STOREU(o_huUpdateRight, ZEROV_R());
      o_maxWaveSpeed = static_cast<real>(0.0);
      return;
    }

    // there is at least one that is not dry dry
    sqrt_hLeft_v = SQRTV(hLeft_v);
    sqrt_hRight_v = SQRTV(hRight_v);
    sqrt_g_hLeft_v = MULV(sqrt_g_v, sqrt_hLeft_v);
    sqrt_g_hRight_v = MULV(sqrt_g_v, sqrt_hRight_v);

    real_vector fWaves[3][2];
    real_vector waveSpeeds[3];

    computeWaveDecomposition_SIMD(fWaves, waveSpeeds, wetDryMask);

    for (int waveNumber = 0; waveNumber < 3; ++waveNumber) {
      const real_vector leftMask = ANDV_R(CMP_LT(waveSpeeds[waveNumber], zeroTol_v_neg), wetDryMask);
      hUpdateLeft_v = BLENDV(hUpdateLeft_v, ADDV(hUpdateLeft_v, fWaves[waveNumber][0]), leftMask);
      huUpdateLeft_v = BLENDV(huUpdateLeft_v, ADDV(huUpdateLeft_v, fWaves[waveNumber][1]), leftMask);

      const real_vector rightMask = ANDNOTV_R(leftMask,(CMP_GT(waveSpeeds[waveNumber], zeroTol_v), wetDryMask));
      hUpdateRight_v  = BLENDV(hUpdateRight_v, ADDV(hUpdateRight_v, fWaves[waveNumber][0]), rightMask);
      huUpdateRight_v = BLENDV(huUpdateRight_v, ADDV(huUpdateRight_v, fWaves[waveNumber][1]), rightMask);

      const real_vector middleMask = ANDV_R(ANDV_R(NOTV_R(leftMask), NOTV_R(rightMask)), wetDryMask);
      if (MOVEMASK(middleMask) != 0) {
        //should not happen
        static const real_vector half          = SETV_R(static_cast<real>(0.5));
        const real_vector        half_hUpdate  = MULV(half, fWaves[waveNumber][0]);
        const real_vector        half_huUpdate = MULV(half, fWaves[waveNumber][1]);

        hUpdateLeft_v  = BLENDV(hUpdateLeft_v, ADDV(hUpdateLeft_v, half_hUpdate), middleMask);
        huUpdateLeft_v = BLENDV(huUpdateLeft_v, ADDV(huUpdateLeft_v, half_huUpdate), middleMask);

        hUpdateRight_v  = BLENDV(hUpdateRight_v, ADDV(hUpdateRight_v, half_hUpdate), middleMask);
        huUpdateRight_v = BLENDV(huUpdateRight_v, ADDV(huUpdateRight_v, half_huUpdate), middleMask);
      }
      assert(MOVEMASK(ORV_R(NOTV_R(wetDryMask), ORV_R(leftMask, ORV_R(rightMask, middleMask)))) == VECTOR_FULL_MASK);
    }


    STOREU(o_hUpdateLeft, BLENDV(ZEROV_R(), hUpdateLeft_v, wetDryMask));
    STOREU(o_huUpdateLeft, BLENDV(ZEROV_R(), huUpdateLeft_v, wetDryMask));
    STOREU(o_hUpdateRight, BLENDV(ZEROV_R(), hUpdateRight_v, wetDryMask));
    STOREU(o_huUpdateRight, BLENDV(ZEROV_R(), huUpdateRight_v, wetDryMask));

    const real_vector maxEdgeSpeed_v = MAXV(
      MAXV(FABS(ANDV_R(waveSpeeds[0], wetDryMask)), FABS(ANDV_R(waveSpeeds[1], wetDryMask))),
      FABS(ANDV_R(waveSpeeds[2], wetDryMask))
    );
    const real* const    p_speed  = reinterpret_cast<const real*>(&maxEdgeSpeed_v);
    const integer* const p_wetDry = reinterpret_cast<const integer*>(&wetDryState_v);

    o_maxWaveSpeed = static_cast<real>(0);
    for (size_t i = 0; i < VECTOR_LENGTH; ++i) {
      if (p_wetDry[i] != 0) {
        o_maxWaveSpeed = std ::max(o_maxWaveSpeed, p_speed[i]);
      }
    }
  }

private:

  void computeWaveDecomposition_SIMD(real_vector o_fWaves[3][2], real_vector o_waveSpeeds[3], const real_vector wetDryMask)
  {
    static const real_vector half = SETV_R(static_cast<real>(0.5));

    // compute eigenvalues of the jacobian matrices in states Q_{i-1} and Q_{i} (char. speeds)
    const real_vector characteristicSpeeds[2] = {SUBV(uLeft_v, sqrt_g_hLeft_v), ADDV(uRight_v, sqrt_g_hRight_v)};

    // compute "Roe speeds"
    const real_vector hRoe = MULV(ADDV(hRight_v, hLeft_v), half);
    const real_vector uRoe = BLENDV(
      ZEROV_R(),
      DIVV(ADDV(MULV(uLeft_v, sqrt_hLeft_v), MULV(uRight_v, sqrt_hRight_v)), ADDV(sqrt_hLeft_v, sqrt_hRight_v)),
      wetDryMask
    );

    const real_vector sqrt_g_hRoe = SQRTV(MULV(g_v, hRoe));

    const real_vector roeSpeeds[2] = {SUBV(uRoe, sqrt_g_hRoe), ADDV(uRoe, sqrt_g_hRoe)};

    const real_vector middleState_blendMask = ANDV_R(
      NOTV_R(ORV_R(
        CAST_INT_TO_REAL_V(CMP_EQ_I(wetDryState_v, WetDryWall_V)),
        CAST_INT_TO_REAL_V(CMP_EQ_I(wetDryState_v, DryWetWall_V))
      )),
      wetDryMask
    );
    // compute the middle state of the homogeneous Riemann-Problem
    computeMiddleState_SIMD(hLeft_v, hRight_v, uLeft_v, uRight_v, middleState_blendMask);

    // compute extended eindfeldt speeds (einfeldt speeds + middle state speeds)
    //   \cite[ch. 5.2]{george2008augmented}, \cite[ch. 6.8]{george2006finite}
    real_vector extEinfeldtSpeeds[2] = {ZEROV_R(), ZEROV_R()};

    extEinfeldtSpeeds[0] = MINV(characteristicSpeeds[0], roeSpeeds[0]);
    extEinfeldtSpeeds[1] = MAXV(roeSpeeds[1], middleStateSpeeds_v[0]);

    const real_vector speeds_t0[2] = {
      MINV(extEinfeldtSpeeds[0], middleStateSpeeds_v[1]), MAXV(extEinfeldtSpeeds[1], characteristicSpeeds[1])};
    const real_vector speeds_t1[2] = {
      MINV(roeSpeeds[0], middleStateSpeeds_v[1]), MAXV(characteristicSpeeds[1], roeSpeeds[1])};

    const real_vector wetMask = ANDV_R(
      CAST_INT_TO_REAL_V(ORV_I(
        ORV_I(CMP_EQ_I(wetDryState_v, WetWet_V), CMP_EQ_I(wetDryState_v, WetDryWall_V)),
        CMP_EQ_I(wetDryState_v, DryWetWall_V)
      )),
      wetDryMask
    );

    real_vector mask_accumulator = wetMask;

    extEinfeldtSpeeds[0] = BLENDV(extEinfeldtSpeeds[0], speeds_t0[0], wetMask);
    extEinfeldtSpeeds[1] = BLENDV(extEinfeldtSpeeds[1], speeds_t0[1], wetMask);

    const real_vector leftDryMask = ANDV_R(ANDNOTV_R(mask_accumulator, CMP_LT(hLeft_v, dryTol_v)), wetDryMask);

    extEinfeldtSpeeds[0] = BLENDV(extEinfeldtSpeeds[0], speeds_t1[0], leftDryMask);
    extEinfeldtSpeeds[1] = BLENDV(extEinfeldtSpeeds[1], speeds_t1[1], leftDryMask);

    assert(MOVEMASK(ANDV_R(leftDryMask, CMP_GE(middleStateSpeeds_v[0], extEinfeldtSpeeds[1]))) == 0);

    // HLL middle state
    //   \cite[theorem 3.1]{george2006finite}, \cite[ch. 4.1]{george2008augmented}
    const real_vector hLLMiddleHeight = MAXV(
      DIVV(
        SUBV(ADDV(SUBV(huLeft_v, huRight_v), MULV(extEinfeldtSpeeds[1], hRight_v)), MULV(extEinfeldtSpeeds[0], hLeft_v)),
        SUBV(extEinfeldtSpeeds[1], extEinfeldtSpeeds[0])
      ),
      ZEROV_R()
    );

#pragma message \
  "strongRarefaction set to false; further investigation needed, about senseless initializations of eigenValues[1]"
    //		static const bool strongRarefaction = false;
    static const real_vector one            = SETV_R(static_cast<real>(1));
    const real_vector        eigenValues[3] = {
             extEinfeldtSpeeds[0], MULV(half, ADDV(extEinfeldtSpeeds[0], extEinfeldtSpeeds[1])), extEinfeldtSpeeds[1]};
    assert(checkVector(eigenValues[0]));
    assert(checkVector(eigenValues[1]));
    assert(checkVector(eigenValues[2]));

    const real_vector eigenVectors[3][3] = {
      {one, ZEROV_R(), one},
      {eigenValues[0], ZEROV_R(), eigenValues[2]},
      {MULV(eigenValues[0], eigenValues[0]), one, MULV(eigenValues[2], eigenValues[2])}};
    assert(checkVector(eigenVectors[0][0]));
    assert(checkVector(eigenVectors[0][1]));
    assert(checkVector(eigenVectors[0][2]));
    assert(checkVector(eigenVectors[1][0]));
    assert(checkVector(eigenVectors[1][1]));
    assert(checkVector(eigenVectors[1][2]));
    assert(checkVector(eigenVectors[2][0]));
    assert(checkVector(eigenVectors[2][1]));
    assert(checkVector(eigenVectors[2][2]));

    // compute steady state wave
    //   \cite[ch. 4.2.1 \& app. A]{george2008augmented}, \cite[ch. 6.2 \& ch. 4.4]{george2006finite}
    const real_vector hBar = MULV(ADDV(hLeft_v, hRight_v), half);

    const real_vector difB = SUBV(bLeft_v, bRight_v);

    real_vector steadyStateWave[2] = {difB, MULV(difB, MULV(hBar, g_v))};

    // preserve depth-positivity
    //   \cite[ch. 6.5.2]{george2006finite}, \cite[ch. 4.2.3]{george2008augmented}
    mask_accumulator                = ZEROV_R();
    const real_vector subsonic_mask = ANDV_R(
      ANDV_R(CMP_LT(eigenValues[0], zeroTol_v_neg), CMP_GT(eigenValues[2], zeroTol_v)), wetDryMask
    );
    mask_accumulator                        = ORV_R(mask_accumulator, subsonic_mask);
    const real_vector supersonic_right_mask = ANDV_R(
      ANDNOTV_R(mask_accumulator, CMP_GT(eigenValues[0], zeroTol_v)), wetDryMask
    );
    mask_accumulator                       = ORV_R(mask_accumulator, supersonic_right_mask);
    const real_vector supersonic_left_mask = ANDV_R(
      ANDNOTV_R(mask_accumulator, CMP_LT(eigenValues[2], zeroTol_v_neg)), wetDryMask
    );
    mask_accumulator = ORV_R(mask_accumulator, supersonic_left_mask);

    if (MOVEMASK(mask_accumulator)) {
      static const real_vector minus_one           = SETV_R(static_cast<real>(-1));
      const real_vector        hLeft_neg           = MULV(hLeft_v, minus_one);
      const real_vector        depth_positivity_t0 = MULV(hLLMiddleHeight, SUBV(eigenValues[2], eigenValues[0]));

      const real_vector depth_positivity_t1 = DIVV(depth_positivity_t0, eigenValues[0]);
      const real_vector depth_positivity_t2 = DIVV(depth_positivity_t0, eigenValues[2]);

      const real_vector depth_positivity_max1 = MAXV(steadyStateWave[0], depth_positivity_t1);
      const real_vector depth_positivity_max2 = MAXV(steadyStateWave[0], hLeft_neg);
      const real_vector depth_positivity_max3 = MAXV(steadyStateWave[0], depth_positivity_t2);

      const real_vector depth_positivity_min1 = MINV(depth_positivity_max1, depth_positivity_t2);
      const real_vector depth_positivity_min2 = MINV(depth_positivity_max2, depth_positivity_t1);
      const real_vector depth_positivity_min3 = MINV(depth_positivity_max3, hRight_v);

      steadyStateWave[0] = BLENDV(steadyStateWave[0], depth_positivity_min1, subsonic_mask);
      steadyStateWave[0] = BLENDV(steadyStateWave[0], depth_positivity_min2, supersonic_right_mask);
      steadyStateWave[0] = BLENDV(steadyStateWave[0], depth_positivity_min3, supersonic_left_mask);
    }

    const real_vector difB_left  = MULV(hLeft_v, difB);
    const real_vector difB_right = MULV(hRight_v, difB);

    const real_vector difB_max = MULV(g_v, MAXV(difB_left, difB_right));
    const real_vector difB_min = MULV(g_v, MINV(difB_left, difB_right));

    // Limit the effect of the source term
    //   \cite[ch. 6.4.2]{george2006finite}
    steadyStateWave[1] = MINV(steadyStateWave[1], difB_max);
    steadyStateWave[1] = MAXV(steadyStateWave[1], difB_min);

    // TODO: Optimize
    const real_vector rightHandSide[3] = {
      SUBV(SUBV(hRight_v, hLeft_v), steadyStateWave[0]),
      SUBV(huRight_v, huLeft_v),
      SUBV(
        SUBV(
          ADDV(MULV(huRight_v, uRight_v), MULV(MULV(MULV(half, g_v), hRight_v), hRight_v)),
          ADDV(MULV(huLeft_v, uLeft_v), MULV(MULV(MULV(half, g_v), hLeft_v), hLeft_v))
        ),
        steadyStateWave[1]
      )};

    assert(checkVector(rightHandSide[0]));
    assert(checkVector(rightHandSide[1]));
    assert(checkVector(rightHandSide[2]));

    real_vector beta[3];

    // everything is ready, solve the equations!
    solveLinearEquation_SIMD(eigenVectors, rightHandSide, beta);
    beta[0] = BLENDV(ZEROV_R(), beta[0], wetDryMask);
    beta[1] = BLENDV(ZEROV_R(), beta[1], wetDryMask);
    beta[2] = BLENDV(ZEROV_R(), beta[2], wetDryMask);
    assert(checkVector(beta[0]));
    assert(checkVector(beta[1]));
    assert(checkVector(beta[2]));

    // compute f-waves and wave-speeds
    const real_vector DryWetWall_Mask  = ANDV_R(CAST_INT_TO_REAL_V(CMP_EQ_I(wetDryState_v, DryWetWall_V)), wetDryMask);
    const integer     DryWetWall_State = MOVEMASK(DryWetWall_Mask);
    const real_vector WetDryWall_Mask  = ANDV_R(
      ANDNOTV_R(DryWetWall_Mask, CAST_INT_TO_REAL_V(CMP_EQ_I(wetDryState_v, WetDryWall_V))), wetDryMask
    );
    const integer WetDryWall_State = MOVEMASK(WetDryWall_Mask);

    const real_vector noWall_Mask  = ANDNOTV_R(ORV_R(DryWetWall_Mask, WetDryWall_Mask), wetDryMask);
    const integer     noWall_State = MOVEMASK(noWall_Mask);

    assert(
      MOVEMASK(ORV_R(NOTV_R(wetDryMask), ORV_R(DryWetWall_Mask, ORV_R(WetDryWall_Mask, noWall_Mask))))
      == VECTOR_FULL_MASK
    );


    if (noWall_State) {
      // compute f-waves (default)
      for (int waveNumber = 0; waveNumber < 3; ++waveNumber) {
        o_fWaves[waveNumber][0] = BLENDV(
          o_fWaves[waveNumber][0], MULV(beta[waveNumber], eigenVectors[1][waveNumber]), noWall_Mask
        );
        o_fWaves[waveNumber][1] = BLENDV(
          o_fWaves[waveNumber][1], MULV(beta[waveNumber], eigenVectors[2][waveNumber]), noWall_Mask
        );
      }

      o_waveSpeeds[0] = BLENDV(o_waveSpeeds[0], eigenValues[0], noWall_Mask);
      o_waveSpeeds[1] = BLENDV(o_waveSpeeds[1], eigenValues[1], noWall_Mask);
      o_waveSpeeds[2] = BLENDV(o_waveSpeeds[2], eigenValues[2], noWall_Mask);
    }

    if (WetDryWall_State) {
      // zero ghost updates (wall boundary)
      // care about the left going wave (0) only
      o_fWaves[0][0] = BLENDV(o_fWaves[0][0], MULV(beta[0], eigenVectors[1][0]), WetDryWall_Mask);
      o_fWaves[0][1] = BLENDV(o_fWaves[0][1], MULV(beta[0], eigenVectors[2][0]), WetDryWall_Mask);

      o_waveSpeeds[0] = BLENDV(o_waveSpeeds[0], eigenValues[0], WetDryWall_Mask);

      // set the rest to zero
      o_fWaves[1][0] = BLENDV(o_fWaves[1][0], ZEROV_R(), WetDryWall_Mask);
      o_fWaves[1][1] = BLENDV(o_fWaves[1][1], ZEROV_R(), WetDryWall_Mask);
      o_fWaves[2][0] = BLENDV(o_fWaves[2][0], ZEROV_R(), WetDryWall_Mask);
      o_fWaves[2][1] = BLENDV(o_fWaves[2][1], ZEROV_R(), WetDryWall_Mask);

      o_waveSpeeds[1] = BLENDV(o_waveSpeeds[1], ZEROV_R(), WetDryWall_Mask);
      o_waveSpeeds[2] = BLENDV(o_waveSpeeds[2], ZEROV_R(), WetDryWall_Mask);
      assert(MOVEMASK(ANDV_R(WetDryWall_Mask, CMP_GE(eigenValues[0], zeroTol_v))) == 0);
    }
    if (DryWetWall_State) {
      // zero ghost updates (wall boundary)
      // care about the right going wave (2) only
      o_fWaves[2][0] = BLENDV(o_fWaves[2][0], MULV(beta[2], eigenVectors[1][2]), DryWetWall_Mask);
      o_fWaves[2][1] = BLENDV(o_fWaves[2][1], MULV(beta[2], eigenVectors[2][2]), DryWetWall_Mask);

      o_waveSpeeds[2] = BLENDV(o_waveSpeeds[2], eigenValues[2], DryWetWall_Mask);

      // set the rest to zero
      o_fWaves[0][0] = BLENDV(o_fWaves[0][0], ZEROV_R(), DryWetWall_Mask);
      o_fWaves[0][1] = BLENDV(o_fWaves[0][1], ZEROV_R(), DryWetWall_Mask);
      o_fWaves[1][0] = BLENDV(o_fWaves[1][0], ZEROV_R(), DryWetWall_Mask);
      o_fWaves[1][1] = BLENDV(o_fWaves[1][1], ZEROV_R(), DryWetWall_Mask);

      o_waveSpeeds[0] = BLENDV(o_waveSpeeds[0], ZEROV_R(), DryWetWall_Mask);
      o_waveSpeeds[1] = BLENDV(o_waveSpeeds[1], ZEROV_R(), DryWetWall_Mask);
      assert(MOVEMASK(ANDV_R(DryWetWall_Mask, CMP_LE(eigenValues[2], zeroTol_v_neg))) == 0);
    }
  }
  
  void determineWetDryState() 
  {
    const real_vector mask_left = CMP_GT(hLeft_v, dryTol_v);

    uLeft_v  = BLENDV(ZEROV_R(), DIVV(huLeft_v, hLeft_v), mask_left);
    bLeft_v  = BLENDV(ADDV(bLeft_v, hLeft_v), bLeft_v, mask_left);
    hLeft_v  = BLENDV(ZEROV_R(), hLeft_v, mask_left);
    huLeft_v = BLENDV(ZEROV_R(), huLeft_v, mask_left);

    const real_vector mask_right = CMP_GT(hRight_v, dryTol_v);

    uRight_v  = BLENDV(ZEROV_R(), DIVV(huRight_v, hRight_v), mask_right);
    bRight_v  = BLENDV(ADDV(bRight_v, hRight_v), bRight_v, mask_right);
    hRight_v  = BLENDV(ZEROV_R(), hRight_v, mask_right);
    huRight_v = BLENDV(ZEROV_R(), huRight_v, mask_right);

    const real_vector wetWetMask = ANDV_R(CMP_GE(hLeft_v, dryTol_v), CMP_GE(hRight_v, dryTol_v));
    wetDryState_v                = BLENDV_I(wetDryState_v, WetWet_V, wetWetMask);

    if (MOVEMASK(wetWetMask) == VECTOR_FULL_MASK) {
      return;
    }

    real_vector mask_accumulator = wetWetMask;

    {
      const real_vector dryDryMask = ANDNOTV_R(
        mask_accumulator, ANDV_R(CMP_LT(hLeft_v, dryTol_v), CMP_LT(hRight_v, dryTol_v))
      );
      wetDryState_v    = BLENDV_I(wetDryState_v, DryDry_V, dryDryMask);
      mask_accumulator = ORV_R(mask_accumulator, dryDryMask);
    }
    {
      const real_vector dryWetInMask = ANDNOTV_R(
        mask_accumulator, ANDV_R(CMP_LT(hLeft_v, dryTol_v), CMP_GT(ADDV(hRight_v, bRight_v), bLeft_v))
      );
      wetDryState_v    = BLENDV_I(wetDryState_v, DryWetInundation_V, dryWetInMask);
      mask_accumulator = ORV_R(mask_accumulator, dryWetInMask);
    }
    {
      const real_vector wetDryInMask = ANDNOTV_R(
        mask_accumulator, ANDV_R(CMP_LT(hRight_v, dryTol_v), CMP_GT(ADDV(hLeft_v, bLeft_v), bRight_v))
      );
      wetDryState_v    = BLENDV_I(wetDryState_v, WetDryInundation_V, wetDryInMask);
      mask_accumulator = ORV_R(mask_accumulator, wetDryInMask);
    }
    const real_vector        wallMask_left = ANDNOTV_R(mask_accumulator, CMP_LT(hLeft_v, dryTol_v));
    static const real_vector minus_one     = SETV_R(static_cast<real>(-1));
    if (MOVEMASK(wallMask_left)) {
      const real_vector uRight_v_neg = MULV(uRight_v, minus_one);

      computeMiddleState_SIMD(hRight_v, hRight_v, uRight_v_neg, uRight_v, wallMask_left, maxNumberOfNewtonIterations);

      const real_vector middleMask = CMP_GT(ADDV(hMiddle_v, bRight_v), bLeft_v);
      {
        const real_vector dryWetWallInMask = ANDV_R(wallMask_left, middleMask);
        wetDryState_v                      = BLENDV_I(wetDryState_v, DryWetWallInundation_V, dryWetWallInMask);
      }
      {
        const real_vector dryWetWallMask = ANDV_R(wallMask_left, NOTV_R(middleMask));
        wetDryState_v                    = BLENDV_I(wetDryState_v, DryWetWall_V, dryWetWallMask);

        hLeft_v  = BLENDV(hLeft_v, hRight_v, dryWetWallMask);
        uLeft_v  = BLENDV(uLeft_v, uRight_v_neg, dryWetWallMask);
        huLeft_v = BLENDV(huLeft_v, MULV(huRight_v, minus_one), dryWetWallMask);
        bLeft_v  = BLENDV(bLeft_v, ZEROV_R(), dryWetWallMask);
        bRight_v = BLENDV(bRight_v, ZEROV_R(), dryWetWallMask);
      }
    }

    mask_accumulator = ORV_R(mask_accumulator, wallMask_left);

    const real_vector wallMask_right = ANDNOTV_R(mask_accumulator, CMP_LT(hRight_v, dryTol_v));
    if (MOVEMASK(wallMask_right)) {
      const real_vector uLeft_v_neg = MULV(uLeft_v, minus_one);

      computeMiddleState_SIMD(hLeft_v, hLeft_v, uLeft_v, uLeft_v_neg, wallMask_right, maxNumberOfNewtonIterations);

      const real_vector middleMask = CMP_GT(ADDV(hMiddle_v, bLeft_v), bRight_v);
      {
        const real_vector wetDryWallInMask = ANDV_R(wallMask_right, middleMask);
        wetDryState_v                      = BLENDV_I(wetDryState_v, WetDryWallInundation_V, wetDryWallInMask);
      }
      {
        const real_vector wetDryWallMask = ANDV_R(wallMask_right, NOTV_R(middleMask));
        wetDryState_v                    = BLENDV_I(wetDryState_v, WetDryWall_V, wetDryWallMask);

        hRight_v  = BLENDV(hRight_v, hLeft_v, wetDryWallMask);
        uRight_v  = BLENDV(uRight_v, uLeft_v_neg, wetDryWallMask);
        huRight_v = BLENDV(huRight_v, MULV(huLeft_v, minus_one), wetDryWallMask);
        bRight_v  = BLENDV(bRight_v, ZEROV_R(), wetDryWallMask);
        bLeft_v   = BLENDV(bLeft_v, ZEROV_R(), wetDryWallMask);
      }
    }

    mask_accumulator = ORV_R(mask_accumulator, wallMask_right);
    assert(MOVEMASK(mask_accumulator) == VECTOR_FULL_MASK);

    const real_vector dryWetWallInMask = CAST_INT_TO_REAL_V(CMP_EQ_I(wetDryState_v, DryWetWallInundation_V));
    const real_vector wetDryWallInMask = CAST_INT_TO_REAL_V(CMP_EQ_I(wetDryState_v, WetDryWallInundation_V));

    bLeft_v  = BLENDV(bLeft_v, ADDV(hRight_v, bRight_v), dryWetWallInMask);
    bRight_v = BLENDV(bRight_v, ADDV(hLeft_v, bLeft_v), wetDryWallInMask);
  }

private:
  inline void computeMiddleState_SIMD(
    const real_vector i_hLeft,
    const real_vector i_hRight,
    const real_vector i_uLeft,
    const real_vector i_uRight,
    const real_vector blendMask,
    const int         i_maxNumberOfNewtonIterations = 1
  ) {
    static const real_vector two       = SETV_R(static_cast<real>(2));
    static const real_vector three     = SETV_R(static_cast<real>(3));
    static const real_vector half_g    = MULV(SETV_R(static_cast<real>(0.5)), g_v);
    static const real_vector quarter_g = MULV(SETV_R(static_cast<real>(0.25)), g_v);

    const real_vector newtonTol_v = SETV_R(newtonTol);

    // set everything to zero
    real_vector l_hMiddle              = ZEROV_R();
    real_vector l_middleStateSpeeds[2] = {ZEROV_R(), ZEROV_R()};

    // compute local square roots
    //(not necessarily the same ones as in computeNetUpdates!)
    const real_vector l_sqrt_g_hRight = SQRTV(MULV(g_v, i_hRight));
    const real_vector l_sqrt_g_hLeft  = SQRTV(MULV(g_v, i_hLeft));

    // single rarefaction in the case of a wet/dry interface
    const real_vector DrySingleMask = BLENDV(ZEROV_R(), CMP_LT(i_hLeft, dryTol_v), blendMask);
    const real_vector SingleDryMask = BLENDV(
      ZEROV_R(), ANDNOTV_R(DrySingleMask, CMP_LT(i_hRight, dryTol_v)), blendMask
    );
    const real_vector Riemann_blendMask = BLENDV(ZEROV_R(), NOTV_R(ORV_R(DrySingleMask, SingleDryMask)), blendMask);

    const real_vector speeds_dry_single = SUBV(i_uRight, MULV(two, l_sqrt_g_hRight));
    riemannStructure_v                  = BLENDV_I(riemannStructure_v, DrySingleRarefaction_V, DrySingleMask);
    const real_vector speeds_single_dry = ADDV(i_uLeft, MULV(two, l_sqrt_g_hLeft));
    riemannStructure_v                  = BLENDV_I(riemannStructure_v, SingleRarefactionDry_V, SingleDryMask);

    if (MOVEMASK(Riemann_blendMask) != 0) {
      // determine the wave structure of the Riemann-problem
      riemannStructure_v = determineRiemannStructure_SIMD(i_hLeft, i_hRight, i_uLeft, i_uRight, Riemann_blendMask);

      // will be computed later
      real_vector sqrt_g_hMiddle = ZEROV_R();

      const real_vector ShockShockMask = CAST_INT_TO_REAL_V(riemannStructure_v);
      if (MOVEMASK(ShockShockMask)) {
        /* Compute All-Shock Riemann Solution
         * \cite[ch. 13.7]{leveque2002finite}
         */
        real_vector iteration_variable = MINV(i_hLeft, i_hRight);

        for (int i = 0; i < i_maxNumberOfNewtonIterations; ++i) {
          const real_vector l_sqrtTermH[2] = {
            SQRTV(MULV(half_g, DIVV(ADDV(iteration_variable, i_hLeft), MULV(iteration_variable, i_hLeft)))),
            SQRTV(MULV(half_g, DIVV(ADDV(iteration_variable, i_hRight), MULV(iteration_variable, i_hRight))))};

          const real_vector delta_h_left  = SUBV(iteration_variable, i_hLeft);
          const real_vector delta_h_right = SUBV(iteration_variable, i_hRight);

          const real_vector phi = ADDV(
            SUBV(i_uRight, i_uLeft), ADDV(MULV(delta_h_left, l_sqrtTermH[0]), MULV(delta_h_right, l_sqrtTermH[1]))
          );

          if (MOVEMASK(ANDV_R(ShockShockMask, CMP_GT(FABS(phi), newtonTol_v))) == 0) {
            break;
          }

          const real_vector quarter_g_div_hMiddle_squared = DIVV(
            quarter_g, MULV(iteration_variable, iteration_variable)
          );

          const real_vector derivativePhi = SUBV(
            ADDV(l_sqrtTermH[0], l_sqrtTermH[1]),
            MULV(
              quarter_g_div_hMiddle_squared,
              ADDV(DIVV(delta_h_left, l_sqrtTermH[0]), DIVV(delta_h_right, l_sqrtTermH[1]))
            )
          );

          iteration_variable = SUBV(iteration_variable, DIVV(phi, derivativePhi));

          assert(MOVEMASK(ORV_R(NOTV_R(ShockShockMask), CMP_GE(iteration_variable, dryTol_v))) == VECTOR_FULL_MASK);
        }

        l_hMiddle = BLENDV(l_hMiddle, iteration_variable, ShockShockMask);
      }

      const real_vector RarefactionRarefactionMask = CAST_INT_TO_REAL_V(SHIFT_LEFT(riemannStructure_v, 1));
      if (MOVEMASK(RarefactionRarefactionMask)) {
        /* Compute All-Rarefaction Riemann Solution
         * \cite[ch. 13.8.6]{leveque2002finite}
         */
        static const real_vector one                = SETV_R(static_cast<real>(1));
        static const real_vector sixteen            = SETV_R(static_cast<real>(16));
        real_vector              iteration_variable = MAXV(
          ZEROV_R(), ADDV(SUBV(i_uLeft, i_uRight), MULV(two, ADDV(l_sqrt_g_hLeft, l_sqrt_g_hRight)))
        );
        iteration_variable = MULV(DIVV(one, MULV(sixteen, g_v)), MULV(iteration_variable, iteration_variable));

        l_hMiddle = BLENDV(l_hMiddle, iteration_variable, RarefactionRarefactionMask);
      }

      sqrt_g_hMiddle = BLENDV(
        sqrt_g_hMiddle,
        SQRTV(MULV(g_v, l_hMiddle)),
        BLENDV(ZEROV_R(), ORV_R(ShockShockMask, RarefactionRarefactionMask), blendMask)
      );

      const real_vector ShockRarefactionMask = BLENDV(
        ZEROV_R(), CAST_INT_TO_REAL_V(SHIFT_LEFT(riemannStructure_v, 2)), blendMask
      );
      const real_vector RarefactionShockMask = BLENDV(
        ZEROV_R(), CAST_INT_TO_REAL_V(SHIFT_LEFT(riemannStructure_v, 3)), blendMask
      );
      const real_vector MixShockRarefactionMask = ORV_R(ShockRarefactionMask, RarefactionShockMask);
      if (MOVEMASK(MixShockRarefactionMask)) {
        /* Compute dam-break Riemann-solution
         * TODO: reference
         */
        const real_vector hMin = BLENDV(i_hRight, i_hLeft, ShockRarefactionMask);
        const real_vector hMax = BLENDV(i_hLeft, i_hRight, ShockRarefactionMask);

        real_vector       iteration_variable      = hMin;
        real_vector       iteration_variable_sqrt = SQRTV(MULV(g_v, iteration_variable));
        const real_vector sqrt_g_hMax             = SQRTV(MULV(g_v, hMax));

        for (int i = 0; i < i_maxNumberOfNewtonIterations; ++i) {
          const real_vector sqrtTermHMin = SQRTV(
            MULV(half_g, DIVV(ADDV(iteration_variable, hMin), MULV(iteration_variable, hMin)))
          );

          const real_vector delta_h_min = SUBV(iteration_variable, hMin);

          const real_vector phi = ADDV(
            SUBV(i_uRight, i_uLeft),
            ADDV(MULV(delta_h_min, sqrtTermHMin), MULV(two, SUBV(iteration_variable_sqrt, sqrt_g_hMax)))
          );

          if (MOVEMASK(ANDV_R(MixShockRarefactionMask, CMP_GT(FABS(phi), newtonTol_v))) == 0) {
            break;
          }

          const real_vector derivativePhi = ADDV(
            SUBV(
              sqrtTermHMin,
              MULV(quarter_g, DIVV(delta_h_min, MULV(MULV(iteration_variable, iteration_variable), sqrtTermHMin)))
            ),
            DIVV(sqrt_g_v, iteration_variable_sqrt)
          );

          iteration_variable = SUBV(iteration_variable, DIVV(phi, derivativePhi));

          iteration_variable_sqrt = SQRTV(MULV(g_v, iteration_variable));
        }

        l_hMiddle      = BLENDV(l_hMiddle, iteration_variable, MixShockRarefactionMask);
        sqrt_g_hMiddle = BLENDV(sqrt_g_hMiddle, iteration_variable_sqrt, MixShockRarefactionMask);
        assert(checkVector(BLENDV(ZEROV_R(), sqrt_g_hMiddle, blendMask)));
      }

      const real_vector three_sqrt_g_hMiddle = MULV(three, sqrt_g_hMiddle);
      l_middleStateSpeeds[0]                 = ADDV(i_uLeft, SUBV(MULV(two, l_sqrt_g_hLeft), three_sqrt_g_hMiddle));
      l_middleStateSpeeds[1]                 = ADDV(SUBV(i_uRight, MULV(two, l_sqrt_g_hRight)), three_sqrt_g_hMiddle);
      assert(
        MOVEMASK(ORV_R(
          NOTV_R(blendMask),
          ORV_R(
            DrySingleMask,
            ORV_R(SingleDryMask, ORV_R(ShockShockMask, ORV_R(RarefactionRarefactionMask, MixShockRarefactionMask)))
          )
        ))
        == VECTOR_FULL_MASK
      );
    }
    l_middleStateSpeeds[0] = BLENDV(l_middleStateSpeeds[0], speeds_dry_single, DrySingleMask);
    l_middleStateSpeeds[1] = BLENDV(l_middleStateSpeeds[1], speeds_dry_single, DrySingleMask);

    l_middleStateSpeeds[0] = BLENDV(l_middleStateSpeeds[0], speeds_single_dry, SingleDryMask);
    l_middleStateSpeeds[1] = BLENDV(l_middleStateSpeeds[1], speeds_single_dry, SingleDryMask);

    assert(checkVector(BLENDV(ZEROV_R(), l_middleStateSpeeds[0], blendMask)));
    assert(checkVector(BLENDV(ZEROV_R(), l_middleStateSpeeds[1], blendMask)));

    hMiddle_v              = BLENDV(hMiddle_v, l_hMiddle, blendMask);
    middleStateSpeeds_v[0] = BLENDV(middleStateSpeeds_v[0], l_middleStateSpeeds[0], blendMask);
    middleStateSpeeds_v[1] = BLENDV(middleStateSpeeds_v[1], l_middleStateSpeeds[1], blendMask);

    assert(MOVEMASK(BLENDV(ZEROV_R(), CMP_LT(hMiddle_v, ZEROV_R()), blendMask)) == 0);
    assert(checkVector(hMiddle_v));
    assert(checkVector(middleStateSpeeds_v[0]));
    assert(checkVector(middleStateSpeeds_v[1]));
  }

  inline integer_vector determineRiemannStructure_SIMD(
    const real_vector i_hLeft,
    const real_vector i_hRight,
    const real_vector i_uLeft,
    const real_vector i_uRight,
    const real_vector blendMask
  ) const {
    static const real_vector two = SETV_R(static_cast<real>(2));

    integer_vector l_riemannStructure = RarefactionShock_V;

    const real_vector hMin = MINV(i_hLeft, i_hRight);
    const real_vector hMax = MAXV(i_hLeft, i_hRight);
    const real_vector uDif = SUBV(i_uLeft, i_uRight);

    const real_vector condition_value_first      = MULV(two, SUBV(SQRTV(MULV(g_v, hMin)), SQRTV(MULV(g_v, hMax))));
    const real_vector RarefactionRarefactionMask = BLENDV(ZEROV_R(), CMP_GE(condition_value_first, uDif), blendMask);
    l_riemannStructure = BLENDV_I(l_riemannStructure, RarefactionRarefaction_V, RarefactionRarefactionMask);

    if (MOVEMASK(ORV_R(NOTV_R(blendMask), RarefactionRarefactionMask)) != VECTOR_FULL_MASK) {
      real_vector mask_accumulator = RarefactionRarefactionMask;

      static const real_vector one = SETV_R(static_cast<real>(1));
      // Note that the following line will cause errors, if the gravitational constant g changes during the simulation
      static const real_vector half_g = MULV(g_v, SETV_R(static_cast<real>(0.5)));

      const real_vector condition_value_second = MULV(
        SUBV(hMax, hMin), SQRTV(MULV(half_g, ADDV(DIVV(one, hMax), DIVV(one, hMin))))
      );
      const real_vector ShockShockMask = BLENDV(
        ZEROV_R(), ANDNOTV_R(mask_accumulator, CMP_LE(condition_value_second, uDif)), blendMask
      );
      l_riemannStructure = BLENDV_I(l_riemannStructure, ShockShock_V, ShockShockMask);

      mask_accumulator = ORV_R(mask_accumulator, ShockShockMask);

      const real_vector ShockRarefactionMask = BLENDV(
        ZEROV_R(), ANDNOTV_R(mask_accumulator, CMP_LT(i_hLeft, i_hRight)), blendMask
      );
      l_riemannStructure = BLENDV_I(l_riemannStructure, ShockRarefaction_V, ShockRarefactionMask);

    }

    return BLENDV_I(riemannStructure_v, l_riemannStructure, blendMask);
  }

  void solveLinearEquation_SIMD(const real_vector i_matrix[3][3], const real_vector i_b[3], real_vector o_x[3]) const
  {
    // Vectorization works, but results are different from file 5 up to the end
    const real_vector m[3][3] = {
      {SUBV(MULV(i_matrix[1][1], i_matrix[2][2]), MULV(i_matrix[1][2], i_matrix[2][1])),
       SUBV(MULV(i_matrix[0][2], i_matrix[2][1]), MULV(i_matrix[0][1], i_matrix[2][2])),
       SUBV(MULV(i_matrix[0][1], i_matrix[1][2]), MULV(i_matrix[0][2], i_matrix[1][1]))},
      {SUBV(MULV(i_matrix[1][2], i_matrix[2][0]), MULV(i_matrix[1][0], i_matrix[2][2])),
       SUBV(MULV(i_matrix[0][0], i_matrix[2][2]), MULV(i_matrix[0][2], i_matrix[2][0])),
       SUBV(MULV(i_matrix[0][2], i_matrix[1][0]), MULV(i_matrix[0][0], i_matrix[1][2]))},
      {SUBV(MULV(i_matrix[1][0], i_matrix[2][1]), MULV(i_matrix[1][1], i_matrix[2][0])),
       SUBV(MULV(i_matrix[0][1], i_matrix[2][0]), MULV(i_matrix[0][0], i_matrix[2][1])),
       SUBV(MULV(i_matrix[0][0], i_matrix[1][1]), MULV(i_matrix[0][1], i_matrix[1][0]))}};

    /*
     * this version is tested, even though nobody knows whether it is the right one
     */
    const real_vector d = ADDV(
      MULV(i_matrix[0][0], m[0][0]), ADDV(MULV(i_matrix[0][1], m[1][0]), MULV(i_matrix[0][2], m[2][0]))
    );

    static const real_vector one = SETV_R(static_cast<real>(1));
    const real_vector        s   = DIVV(one, d);

    o_x[0] = MULV(ADDV(MULV(m[0][0], i_b[0]), ADDV(MULV(m[0][1], i_b[1]), MULV(m[0][2], i_b[2]))), s);
    o_x[1] = MULV(ADDV(MULV(m[1][0], i_b[0]), ADDV(MULV(m[1][1], i_b[1]), MULV(m[1][2], i_b[2]))), s);
    o_x[2] = MULV(ADDV(MULV(m[2][0], i_b[0]), ADDV(MULV(m[2][1], i_b[1]), MULV(m[2][2], i_b[2]))), s);
  }
};