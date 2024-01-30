#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include "WavePropagationSolver.hpp"
#include "SIMDTypesDouble.hpp"

namespace Solvers
{
    class FWaveSolver_SIMD;
}

class Solvers::FWaveSolver_SIMD {
private:
    // real hLeft_;
    // real hRight_;
    // real huLeft_;
    // real huRight_;
    // real bLeft_;
    // real bRight;
    // real uLeft;
    // real uRight;

    real_vector hLeft_v;
    real_vector hRight_v;
    real_vector huLeft_v;
    real_vector huRight_v;
    real_vector bLeft_v;
    real_vector bRight_v;
    real_vector uLeft_v;
    real_vector uRight_v;

    real_vector notDryDryMask;



    // constants
    const real_vector g_v;
    const real_vector sqrt_g_v;
    const real_vector half_v;
    const real_vector zeroTol_v;
    const real_vector neg_zeroTol_v;
    const real_vector dryTol_v;

    static const integer DryDry               = 0;
    static const integer WetWet               = SHIFT_SIGN_RIGHT(1);
    static const integer WetDryInundation     = SHIFT_SIGN_RIGHT(2);
    static const integer WetDryWall           = SHIFT_SIGN_RIGHT(3);
    static const integer WetDryWallInundation = SHIFT_SIGN_RIGHT(4);
    static const integer DryWetInundation     = SHIFT_SIGN_RIGHT(5);
    static const integer DryWetWall           = SHIFT_SIGN_RIGHT(6);
    static const integer DryWetWallInundation = SHIFT_SIGN_RIGHT(7);
    // static const integer WetWet               = 1;
    // static const integer WetDryInundation     = 2;
    // static const integer WetDryWall           = 3;
    // static const integer WetDryWallInundation = 4;
    // static const integer DryWetInundation     = 5;
    // static const integer DryWetWall           = 6;
    // static const integer DryWetWallInundation = 7;

    const integer_vector DryDry_v;
    const integer_vector WetWet_v;
    const integer_vector WetDryInundation_v;
    const integer_vector WetDryWall_v;
    const integer_vector WetDryWallInundation_v;
    const integer_vector DryWetInundation_v;
    const integer_vector DryWetWall_v;
    const integer_vector DryWetWallInundation_v;
    integer_vector       wetDryState_v;
public:
    FWaveSolver_SIMD(
        real i_dryTolerance                = static_cast<real>(0.01),
        real i_gravity                     = static_cast<real>(9.81),
        real i_newtonTolerance             = static_cast<real>(0.000001),
        int  i_maxNumberOfNewtonIterations = 10,
        real i_zeroTolerance               = static_cast<real>(0.00001)
    ): 
    DryDry_v(SETV_I(DryDry)),
    WetWet_v(SETV_I(WetWet)),
    DryWetInundation_v(SETV_I(DryWetInundation)),
    WetDryInundation_v(SETV_I(WetDryInundation)),
    DryWetWallInundation_v(SETV_I(DryWetWallInundation)),
    DryWetWall_v(SETV_I(DryWetWall)),
    WetDryWallInundation_v(SETV_I(WetDryWallInundation)),
    WetDryWall_v(SETV_I(WetDryWall)),
    g_v(SETV_R(i_gravity)),
    sqrt_g_v(SETV_R(std::sqrt(i_gravity))),
    half_v(SETV_R(static_cast<real>(0.5))),
    zeroTol_v(SETV_R(i_zeroTolerance)),
    neg_zeroTol_v(SETV_R(-i_zeroTolerance)),
    dryTol_v(SETV_R(i_dryTolerance))
    {}

    ~FWaveSolver_SIMD()
    {}

    void computeNetUpdates(
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
        real_vector hUpdateLeft_v   = ZEROV_R();
        real_vector hUpdateRight_v  = ZEROV_R();
        real_vector huUpdateLeft_v  = ZEROV_R();
        real_vector huUpdateRight_v = ZEROV_R();
        real_vector maxWaveSpeed_v = ZEROV_R();

        uLeft_v = ZEROV_R();
        uRight_v = ZEROV_R();

        // o_maxWaveSpeed = 0;

        real_vector waveSpeeds[2];

        hLeft_v = LOADU(i_hLeft);
        hRight_v  = LOADU(i_hRight);
        huLeft_v  = LOADU(i_huLeft);
        huRight_v = LOADU(i_huRight);
        bLeft_v   = LOADU(i_bLeft);
        bRight_v  = LOADU(i_bRight);

        determineWetDryState();

        notDryDryMask = CAST_INT_TO_REAL_V(NOTV_I(CMP_EQ_I(wetDryState_v, DryDry_v)));
        if (MOVEMASK(notDryDryMask) == 0) {
            STOREU(o_hUpdateLeft, ZEROV_R());
            STOREU(o_huUpdateLeft, ZEROV_R());
            STOREU(o_hUpdateRight, ZEROV_R());
            STOREU(o_huUpdateRight, ZEROV_R());
            o_maxWaveSpeed = static_cast<real>(0.0);
            return;
        }
        // 接下来每一步操作都要排除掉vector中DryDry的部分，即与notDryDryMask and 以后为1的部分再继续操作

        computeWaveSpeeds(waveSpeeds);

        computeNetUpdatesWithWaveSpeeds(waveSpeeds, hUpdateLeft_v, hUpdateRight_v, huUpdateLeft_v, huUpdateRight_v, maxWaveSpeed_v);

        // WetDryWall与DryDry相矛盾，因此已经排除了DryDry
        const real_vector rightMask = CMP_EQ_R(wetDryState_v, WetDryWall_v);
        hUpdateRight_v = BLENDV(hUpdateRight_v, ZEROV_R(), rightMask);
        huUpdateRight_v = BLENDV(huUpdateRight_v, ZEROV_R(), rightMask);

        const real_vector leftMask = CMP_EQ_R(wetDryState_v, DryWetWall_v);
        hUpdateLeft_v = BLENDV(hUpdateLeft_v, ZEROV_R(), leftMask);
        huUpdateLeft_v = BLENDV(huUpdateLeft_v, ZEROV_R(), leftMask);

        STOREU(o_hUpdateLeft, hUpdateLeft_v);
        STOREU(o_hUpdateRight, hUpdateRight_v);
        STOREU(o_huUpdateLeft, huUpdateLeft_v);
        STOREU(o_huUpdateRight, huUpdateRight_v);
        o_maxWaveSpeed = static_cast<real>(0.0);
        const real* const pSpeed = reinterpret_cast<const real*>(&maxWaveSpeed_v);
        const integer* const pMask = reinterpret_cast<const integer*>(&wetDryState_v);
        for (size_t i = 0; i < VECTOR_LENGTH; ++i) {
            if (pMask[i] != 0) {
                o_maxWaveSpeed = std::max(o_maxWaveSpeed, pSpeed[i]);
            }
        }
    }

private:
    void determineWetDryState()
    {
        const real_vector mask_left_lt = CMP_LT(hLeft_v, dryTol_v);
        const real_vector mask_right_lt = CMP_LT(hRight_v, dryTol_v);
        const real_vector mask_left_ge = CMP_GE(hLeft_v, dryTol_v);
        const real_vector mask_right_ge = CMP_GE(hRight_v, dryTol_v);
        const real_vector maskDryDry = ANDV_R(mask_left_lt, mask_right_lt);
        const real_vector maskDryWet = ANDV_R(mask_left_lt, mask_right_ge);
        const real_vector maskWetDry = ANDV_R(mask_left_ge, mask_right_lt);
        const real_vector maskWetWet = ANDV_R(mask_left_ge, mask_right_ge);

        hLeft_v = BLENDV(hLeft_v, ADDV(hLeft_v, dryTol_v), mask_left_lt);
        hRight_v = BLENDV(hRight_v, ADDV(hRight_v, dryTol_v), mask_right_lt);

        // DryDry
        wetDryState_v = BLENDV_I(wetDryState_v, DryDry_v, maskDryDry);
        wetDryState_v = BLENDV_I(wetDryState_v, DryWetWall_v, maskDryWet);
        wetDryState_v = BLENDV_I(wetDryState_v, WetDryWall_v, maskWetDry);
        wetDryState_v = BLENDV_I(wetDryState_v, WetWet_v, maskWetWet);
        hLeft_v = BLENDV(hLeft_v, hRight_v, maskDryWet);
        hRight_v = BLENDV(hRight_v, hLeft_v, maskWetDry);
        // DryWet
        uRight_v = BLENDV(uRight_v, DIVV(huRight_v, hRight_v), maskDryWet);
        bLeft_v = BLENDV(bLeft_v, bRight_v, maskDryWet);
        huLeft_v = BLENDV(huLeft_v, MULV(SETV_R(-1.0), huRight_v) , maskDryWet);
        uLeft_v = BLENDV(uLeft_v, MULV(SETV_R(-1.0), uRight_v) , maskDryWet);
        // WetDry
        uLeft_v = BLENDV(uLeft_v, DIVV(huLeft_v, hLeft_v), maskWetDry);
        bRight_v = BLENDV(bRight_v, bLeft_v, maskWetDry);
        huRight_v = BLENDV(huRight_v, MULV(SETV_R(-1.0), huLeft_v) , maskWetDry);
        uRight_v = BLENDV(uRight_v, MULV(SETV_R(-1.0), uLeft_v) , maskWetDry);
        //WetWet
        uLeft_v = BLENDV(uLeft_v, DIVV(huLeft_v, hLeft_v), maskWetWet);
        uRight_v = BLENDV(uRight_v, DIVV(huRight_v, hRight_v), maskWetWet);
        // real_vector a = CAST_INT_TO_REAL_V(WetWet_v);
    }

    void computeWaveSpeeds(real_vector o_waveSpeeds[2]) const
    {
        real_vector sqrt_hLeft_v    = SQRTV(hLeft_v);
        real_vector sqrt_hRight_v   = SQRTV(hRight_v);
        real_vector sqrt_g_hLeft_v  = MULV(sqrt_g_v, sqrt_hLeft_v);
        real_vector sqrt_g_hRight_v = MULV(sqrt_g_v, sqrt_hRight_v);

        real_vector characteristicSpeeds[2];
        characteristicSpeeds[0] = SUBV(uLeft_v, sqrt_g_hLeft_v);
        characteristicSpeeds[1] = ADDV(uRight_v, sqrt_g_hRight_v);

        const real_vector hRoe = MULV(half_v, ADDV(hRight_v, hLeft_v));

        const real_vector uRoe = BLENDV(
            ZEROV_R(),
            DIVV(ADDV(MULV(uLeft_v, sqrt_hLeft_v), MULV(uRight_v, sqrt_hRight_v)), ADDV(sqrt_hLeft_v, sqrt_hRight_v)),
            notDryDryMask
        );

        real_vector roeSpeeds[2];
        const real_vector sqrt_hRoe_v = SQRTV(hRoe);
        roeSpeeds[0] = SUBV(uRoe, MULV(sqrt_g_v, sqrt_hRoe_v));
        roeSpeeds[1] = ADDV(uRoe, MULV(sqrt_g_v, sqrt_hRoe_v));
        
        o_waveSpeeds[0] = MINV(characteristicSpeeds[0], roeSpeeds[0]);
        o_waveSpeeds[1] = MAXV(characteristicSpeeds[1], roeSpeeds[1]);
    }

    void computeNetUpdatesWithWaveSpeeds(
        const real_vector waveSpeeds[2],
        real_vector&      o_hUpdateLeft,
        real_vector&      o_hUpdateRight,
        real_vector&      o_huUpdateLeft,
        real_vector&      o_huUpdateRight,
        real_vector&      o_maxWaveSpeed
    )
    {
        real_vector fWaves[2][2];

        computeWaveDecomposition(waveSpeeds, fWaves);


        for (int waveNumber = 0; waveNumber < 2; ++waveNumber) {
            const real_vector leftMask = CMP_LT(waveSpeeds[waveNumber], neg_zeroTol_v);
            o_hUpdateLeft = BLENDV(o_hUpdateLeft, ADDV(o_hUpdateLeft, fWaves[waveNumber][0]), leftMask);
            o_huUpdateLeft = BLENDV(o_huUpdateLeft, ADDV(o_huUpdateLeft, fWaves[waveNumber][1]), leftMask);

            const real_vector rightMask = CMP_GT(waveSpeeds[waveNumber], zeroTol_v);
            o_hUpdateRight = BLENDV(o_hUpdateRight, ADDV(o_hUpdateRight, fWaves[waveNumber][0]), rightMask);
            o_huUpdateRight = BLENDV(o_huUpdateRight, ADDV(o_huUpdateRight, fWaves[waveNumber][1]), rightMask);

            const real_vector middleMask = ANDV_R(NOTV_R(leftMask), NOTV_R(rightMask));
            if (MOVEMASK(middleMask) != 0) {
                o_hUpdateLeft = BLENDV(o_hUpdateLeft, ADDV(o_hUpdateLeft, MULV(half_v, fWaves[waveNumber][0])), middleMask);
                o_huUpdateLeft = BLENDV(o_huUpdateLeft, ADDV(o_huUpdateLeft, MULV(half_v, fWaves[waveNumber][1])), middleMask);
                o_hUpdateRight = BLENDV(o_hUpdateRight, ADDV(o_hUpdateRight, MULV(half_v, fWaves[waveNumber][0])), middleMask);
                o_huUpdateRight = BLENDV(o_huUpdateRight, ADDV(o_huUpdateRight, MULV(half_v, fWaves[waveNumber][1])), middleMask);
            }
        }
        real_vector a = FABS(waveSpeeds[0]);
        real_vector b = FABS(waveSpeeds[1]);
        o_maxWaveSpeed = MAXV(FABS(ANDV_R(waveSpeeds[0], notDryDryMask)), FABS(ANDV_R(waveSpeeds[1], notDryDryMask)));
    }

    void computeWaveDecomposition(const real_vector waveSpeeds[2], real_vector o_fWaves[2][2]) const
    {
        real_vector lambdaDif = SUBV(waveSpeeds[1], waveSpeeds[0]);

        real_vector Rinv[2][2];

        // may happen: divide by zero
        real_vector oneDivLambdaDif = DIVV(SETV_R(1.0), lambdaDif);
        Rinv[0][0] = MULV(oneDivLambdaDif, waveSpeeds[1]);
        Rinv[0][1] = MULV(SETV_R(-1.0), oneDivLambdaDif);
        Rinv[1][0] = MULV(oneDivLambdaDif, MULV(SETV_R(-1.0), waveSpeeds[0]));
        Rinv[1][1] = oneDivLambdaDif;

        real_vector fDif[2];
        fDif[0] = SUBV(huRight_v, huLeft_v);
        fDif[1] = SUBV(
            ADDV(MULV(huRight_v, uRight_v), MULV(half_v, MULV(g_v, MULV(hRight_v, hRight_v)))), 
            ADDV(MULV(huLeft_v, uLeft_v), MULV(half_v, MULV(g_v, MULV(hLeft_v, hLeft_v))))
        );

        real_vector psi = MULV(SETV_R(-1.0), MULV(g_v, MULV(half_v, MULV(ADDV(hRight_v, hLeft_v), SUBV(bRight_v, bLeft_v)))));
        fDif[1] = SUBV(fDif[1], psi);

        real_vector beta[2];
        beta[0] = ADDV(MULV(Rinv[0][0], fDif[0]), MULV(Rinv[0][1], fDif[1]));
        beta[1] = ADDV(MULV(Rinv[1][0], fDif[0]), MULV(Rinv[1][1], fDif[1]));

        o_fWaves[0][0] = beta[0];
        o_fWaves[0][1] = MULV(beta[0], waveSpeeds[0]);

        o_fWaves[1][0] = beta[1];
        o_fWaves[1][1] = MULV(beta[1], waveSpeeds[1]);
    }
    
};
    
 
 