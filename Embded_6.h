#ifndef HUANG_AXI_ALGORITHM4_CONTACT_CELL_DETECTION_H
#define HUANG_AXI_ALGORITHM4_CONTACT_CELL_DETECTION_H

#include <stdbool.h>

#include "huang_axi_algorithm3_contact_normal.h"

typedef struct
{
    bool available;
    double cs;
    double c;
} HuangAxiAlgorithm4Cell;

typedef struct
{
    HuangAxiAlgorithm4Cell cells[3][3];
} HuangAxiAlgorithm4Stencil;

typedef enum
{
    HUANG_AXI_ALG4_NOT_MIXED_CELL = 0,
    HUANG_AXI_ALG4_NOT_INTERFACIAL_CANDIDATE,
    HUANG_AXI_ALG4_NEIGHBOR_CONTAINS_INTERFACE,
    HUANG_AXI_ALG4_DIAGONAL_CONTAINS_INTERFACE,
    HUANG_AXI_ALG4_INVALID_TANGENT,
    HUANG_AXI_ALG4_DOMAIN_BOUNDARY,
    HUANG_AXI_ALG4_CONFIRMED_CONTACT_LINE
} HuangAxiAlgorithm4Reason;

typedef struct
{
    HuangAxiAlgorithm4Reason reason;
    bool is_contact_line;
    HuangAxiVec2 tangent;
    int tangent_offset_x;
    int tangent_offset_y;
    int diagonal_offset_x;
    int diagonal_offset_y;
} HuangAxiAlgorithm4Result;

static inline bool huang_axi_algorithm4_contains_interface(
    const HuangAxiAlgorithm4Cell stencil_cell, const double tolerance)
{
    return stencil_cell.available && isfinite(stencil_cell.cs) &&
           isfinite(stencil_cell.c) && stencil_cell.cs > tolerance &&
           stencil_cell.c > tolerance &&
           stencil_cell.c < stencil_cell.cs - tolerance;
}

/*
 * Huang Algorithm 4 using the complete-cell fraction contract:
 * c is physical liquid area divided by the complete Cartesian cell area,
 * while cs is open-fluid area divided by that same complete cell area.
 *
 * The tangent orientation and the two axis-neighbor plus diagonal checks
 * follow Chongsen/tmp_fraction_field.h, sort_cell().  Domain-neighbor
 * availability is explicit so a boundary case cannot be silently accepted.
 */
static inline HuangAxiAlgorithm4Result
huang_axi_algorithm4_detect(const HuangAxiAlgorithm4Stencil * stencil,
                            const HuangAxiAlgorithm3Result normal,
                            const double tolerance)
{
    HuangAxiAlgorithm4Result result = {
        .reason = HUANG_AXI_ALG4_INVALID_TANGENT,
        .is_contact_line = false,
        .tangent = {NAN, NAN},
        .tangent_offset_x = 0,
        .tangent_offset_y = 0,
        .diagonal_offset_x = 0,
        .diagonal_offset_y = 0
    };

    if (!stencil || !isfinite(tolerance) || tolerance < 0. ||
        normal.status != HUANG_AXI_ALG3_OK ||
        !huang_axi_vec2_is_finite(normal.solid_normal) ||
        !huang_axi_vec2_is_finite(normal.liquid_normal))
        return result;

    const HuangAxiAlgorithm4Cell center = stencil->cells[1][1];
    if (!center.available || !isfinite(center.cs) || !isfinite(center.c) ||
        center.cs <= tolerance || center.cs >= 1. - tolerance) {
        result.reason = HUANG_AXI_ALG4_NOT_MIXED_CELL;
        return result;
    }
    if (center.c <= tolerance || center.c >= center.cs - tolerance) {
        result.reason = HUANG_AXI_ALG4_NOT_INTERFACIAL_CANDIDATE;
        return result;
    }

    const double cross_ns_nl =
        huang_axi_vec2_cross(normal.solid_normal, normal.liquid_normal);
    if (cross_ns_nl > tolerance) {
        result.tangent.x = normal.liquid_normal.y;
        result.tangent.y = -normal.liquid_normal.x;
    }
    else if (cross_ns_nl < -tolerance) {
        result.tangent.x = -normal.liquid_normal.y;
        result.tangent.y = normal.liquid_normal.x;
    }
    else
        return result;

    HuangAxiVec2 tangent_unit;
    if (!huang_axi_vec2_normalize(result.tangent, tolerance, &tangent_unit))
        return result;
    result.tangent = tangent_unit;
    /* Huang's released sort_cell() uses sign2(), including the zero case.
     * A zero component must probe the center index, not be coerced to a
     * diagonal neighbor. */
    result.tangent_offset_x = result.tangent.x > 0. ? 1 :
                              (result.tangent.x < 0. ? -1 : 0);
    result.tangent_offset_y = result.tangent.y > 0. ? 1 :
                              (result.tangent.y < 0. ? -1 : 0);
    result.diagonal_offset_x = result.tangent_offset_x;
    result.diagonal_offset_y = result.tangent_offset_y;

    const HuangAxiAlgorithm4Cell x_neighbor =
        stencil->cells[1 + result.tangent_offset_x][1];
    const HuangAxiAlgorithm4Cell y_neighbor =
        stencil->cells[1][1 + result.tangent_offset_y];
    const HuangAxiAlgorithm4Cell diagonal =
        stencil->cells[1 + result.diagonal_offset_x]
                      [1 + result.diagonal_offset_y];

    if (!x_neighbor.available || !y_neighbor.available || !diagonal.available) {
        result.reason = HUANG_AXI_ALG4_DOMAIN_BOUNDARY;
        return result;
    }
    if (huang_axi_algorithm4_contains_interface(x_neighbor, tolerance) ||
        huang_axi_algorithm4_contains_interface(y_neighbor, tolerance)) {
        result.reason = HUANG_AXI_ALG4_NEIGHBOR_CONTAINS_INTERFACE;
        return result;
    }
    if (huang_axi_algorithm4_contains_interface(diagonal, tolerance)) {
        result.reason = HUANG_AXI_ALG4_DIAGONAL_CONTAINS_INTERFACE;
        return result;
    }

    result.reason = HUANG_AXI_ALG4_CONFIRMED_CONTACT_LINE;
    result.is_contact_line = true;
    return result;
}

#endif
