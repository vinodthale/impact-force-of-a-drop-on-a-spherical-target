#ifndef HUANG_AXI_ALGORITHM3_CONTACT_NORMAL_H
#define HUANG_AXI_ALGORITHM3_CONTACT_NORMAL_H

#include <float.h>

#include "huang_axi_common_geometry.h"

typedef enum
{
    HUANG_AXI_ALG3_OK = 0,
    HUANG_AXI_ALG3_INVALID_SOLID_NORMAL,
    HUANG_AXI_ALG3_INVALID_GRADIENT,
    HUANG_AXI_ALG3_INVALID_ANGLE,
    HUANG_AXI_ALG3_DEGENERATE_ORIENTATION
} HuangAxiAlgorithm3Status;

typedef enum
{
    HUANG_AXI_ALG3_DEGENERATE = 0,
    HUANG_AXI_ALG3_POSITIVE = 1,
    HUANG_AXI_ALG3_NEGATIVE = -1
} HuangAxiAlgorithm3Orientation;

typedef struct
{
    HuangAxiAlgorithm3Status status;
    HuangAxiAlgorithm3Orientation orientation;
    HuangAxiVec2 solid_normal;
    HuangAxiVec2 liquid_normal;
    double cross_value;
    double contact_angle_residual;
} HuangAxiAlgorithm3Result;

/*
 * Huang Algorithm 3 in the meridional x-y plane.
 *
 * solid_normal is the oriented embedded-solid normal and liquid_gradient
 * supplies only the orientation branch.  The returned normal satisfies
 *
 *     solid_normal . liquid_normal = -cos(contact_angle).
 *
 * The two formula branches are transcribed from Chongsen/TPR2D.h,
 * normal_contact().  A zero cross product is reported as a degeneracy;
 * the routine never guesses an orientation.
 */
static inline HuangAxiAlgorithm3Result
huang_axi_algorithm3_contact_normal(const HuangAxiVec2 solid_normal,
                                    const HuangAxiVec2 liquid_gradient,
                                    const double contact_angle)
{
    const double tolerance = 64.*DBL_EPSILON;
    HuangAxiAlgorithm3Result result = {
        .status = HUANG_AXI_ALG3_INVALID_SOLID_NORMAL,
        .orientation = HUANG_AXI_ALG3_DEGENERATE,
        .solid_normal = {NAN, NAN},
        .liquid_normal = {NAN, NAN},
        .cross_value = NAN,
        .contact_angle_residual = NAN
    };

    if (!isfinite(tolerance) || tolerance < 0. ||
        !huang_axi_vec2_normalize(solid_normal, tolerance,
                                  &result.solid_normal))
        return result;

    HuangAxiVec2 gradient_unit;
    if (!huang_axi_vec2_normalize(liquid_gradient, tolerance, &gradient_unit)) {
        result.status = HUANG_AXI_ALG3_INVALID_GRADIENT;
        return result;
    }

    if (!isfinite(contact_angle) || contact_angle < 0. || contact_angle > M_PI) {
        result.status = HUANG_AXI_ALG3_INVALID_ANGLE;
        return result;
    }

    result.cross_value = huang_axi_vec2_cross(gradient_unit,
                                               result.solid_normal);
    const double cosine = cos(contact_angle);
    const double sine = sin(contact_angle);

    if (result.cross_value > tolerance) {
        result.orientation = HUANG_AXI_ALG3_POSITIVE;
        result.liquid_normal.x = -result.solid_normal.x*cosine +
                                  result.solid_normal.y*sine;
        result.liquid_normal.y = -result.solid_normal.x*sine -
                                  result.solid_normal.y*cosine;
    }
    else if (result.cross_value < -tolerance) {
        result.orientation = HUANG_AXI_ALG3_NEGATIVE;
        result.liquid_normal.x = -result.solid_normal.x*cosine -
                                  result.solid_normal.y*sine;
        result.liquid_normal.y =  result.solid_normal.x*sine -
                                  result.solid_normal.y*cosine;
    }
    else {
        result.status = HUANG_AXI_ALG3_DEGENERATE_ORIENTATION;
        return result;
    }

    HuangAxiVec2 normalized_liquid_normal;
    if (!huang_axi_vec2_normalize(result.liquid_normal, tolerance,
                                  &normalized_liquid_normal)) {
        result.status = HUANG_AXI_ALG3_INVALID_GRADIENT;
        return result;
    }
    result.liquid_normal = normalized_liquid_normal;
    result.contact_angle_residual =
        fabs(huang_axi_vec2_dot(result.solid_normal, result.liquid_normal) +
             cosine);
    result.status = HUANG_AXI_ALG3_OK;
    return result;
}

#endif
