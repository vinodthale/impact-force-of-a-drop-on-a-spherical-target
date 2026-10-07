#ifndef HUANG_AXI_COMMON_GEOMETRY_H
#define HUANG_AXI_COMMON_GEOMETRY_H

#include <math.h>
#include <stdbool.h>

typedef struct
{
    double x;
    double y;
} HuangAxiVec2;

static inline bool huang_axi_vec2_is_finite(const HuangAxiVec2 v)
{
    return isfinite(v.x) && isfinite(v.y);
}

static inline double huang_axi_vec2_dot(const HuangAxiVec2 a, const HuangAxiVec2 b)
{
    return a.x*b.x + a.y*b.y;
}

static inline double huang_axi_vec2_cross(const HuangAxiVec2 a, const HuangAxiVec2 b)
{
    return a.x*b.y - a.y*b.x;
}

static inline double huang_axi_vec2_norm(const HuangAxiVec2 v)
{
    return hypot(v.x, v.y);
}

static inline bool huang_axi_vec2_normalize(const HuangAxiVec2 input,
                                             const double tolerance,
                                             HuangAxiVec2 * output)
{
    if (!output || !huang_axi_vec2_is_finite(input) ||
        !isfinite(tolerance) || tolerance < 0.)
        return false;

    const double magnitude = huang_axi_vec2_norm(input);
    if (!isfinite(magnitude) || magnitude <= tolerance)
        return false;

    output->x = input.x/magnitude;
    output->y = input.y/magnitude;
    return huang_axi_vec2_is_finite(*output);
}

#endif
