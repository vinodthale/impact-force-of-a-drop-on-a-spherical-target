#ifndef HUANG_AXI_ALGORITHM5_EMBEDDED_HEIGHT_FUNCTION_H
#define HUANG_AXI_ALGORITHM5_EMBEDDED_HEIGHT_FUNCTION_H

#include <float.h>
#include <math.h>
#include <stdbool.h>

#include "huang_axi_common_geometry.h"

typedef struct
{
    double p;
    double q;
} HuangAxiHeightPoint2;

typedef enum
{
    HUANG_AXI_ALG5_INVALID = 0,
    HUANG_AXI_ALG5_PRIMARY_HF,
    HUANG_AXI_ALG5_SECONDARY_HF,
    HUANG_AXI_ALG5_RECONSTRUCTED_INTERFACE_HF,
    HUANG_AXI_ALG5_EXTENDED_FRACTION_FALLBACK,
    HUANG_AXI_ALG5_PATH_COUNT
} HuangAxiAlgorithm5Path;

typedef enum
{
    HUANG_AXI_ALG5_OK = 0,
    HUANG_AXI_ALG5_INVALID_INPUT,
    HUANG_AXI_ALG5_NO_VALID_PATH,
    HUANG_AXI_ALG5_SINGULAR_FIT,
    HUANG_AXI_ALG5_NONFINITE_RESULT
} HuangAxiAlgorithm5Status;

typedef struct
{
    bool accumulation_valid;
    HuangAxiHeightPoint2 height_point;
    HuangAxiHeightPoint2 neighbor_height_point;
    HuangAxiHeightPoint2 interface_segment_a;
    HuangAxiHeightPoint2 interface_segment_b;
    HuangAxiHeightPoint2 solid_segment_a;
    HuangAxiHeightPoint2 solid_segment_b;
    double contact_slope;
    double center_q;
    double ghost_q;
    double curvature_q;

    /* Test/reference metadata; not used by the implementation. */
    double expected_coefficients[3];
} HuangAxiAlgorithm5Candidate;

typedef struct
{
    HuangAxiAlgorithm5Candidate primary;
    HuangAxiAlgorithm5Candidate secondary;
    HuangAxiAlgorithm5Candidate reconstructed;
    bool extended_fraction_valid;
    double extended_height_minus;
    double extended_height_center;
    double extended_height_plus;
    double physical_liquid_measure;
} HuangAxiAlgorithm5Input;

typedef struct
{
    HuangAxiAlgorithm5Status status;
    HuangAxiAlgorithm5Path path;
    HuangAxiHeightPoint2 contact_point;
    double coefficients[3];
    double matrix_condition;
    double fit_residual;
    double height_minus;
    double height_center;
    double height_plus;
    double height_gradient;
    double meridional_curvature;
    double axisymmetric_curvature;
    double physical_liquid_measure;
} HuangAxiAlgorithm5Result;

static inline const char *
huang_axi_algorithm5_path_name(const HuangAxiAlgorithm5Path path)
{
    static const char * names[HUANG_AXI_ALG5_PATH_COUNT] = {
        "INVALID",
        "PRIMARY_HF",
        "SECONDARY_HF",
        "RECONSTRUCTED_INTERFACE_HF",
        "EXTENDED_FRACTION_FALLBACK"
    };
    return path >= 0 && path < HUANG_AXI_ALG5_PATH_COUNT ?
           names[path] : "INVALID";
}

static inline bool huang_axi_algorithm5_point_is_finite(
    const HuangAxiHeightPoint2 point)
{
    return isfinite(point.p) && isfinite(point.q);
}

static inline double huang_axi_algorithm5_cross(
    const HuangAxiHeightPoint2 a, const HuangAxiHeightPoint2 b)
{
    return a.p*b.q - a.q*b.p;
}

/* Intersection of two closed 2D segments in local (height, transverse)
 * coordinates.  Collinear overlap is deliberately rejected because it does
 * not define Huang's unique contact point M. */
static inline bool huang_axi_algorithm5_segment_intersection(
    const HuangAxiHeightPoint2 a, const HuangAxiHeightPoint2 b,
    const HuangAxiHeightPoint2 c, const HuangAxiHeightPoint2 d,
    const double tolerance, HuangAxiHeightPoint2 * intersection)
{
    if (!intersection || !huang_axi_algorithm5_point_is_finite(a) ||
        !huang_axi_algorithm5_point_is_finite(b) ||
        !huang_axi_algorithm5_point_is_finite(c) ||
        !huang_axi_algorithm5_point_is_finite(d))
        return false;
    const HuangAxiHeightPoint2 r = {b.p - a.p, b.q - a.q};
    const HuangAxiHeightPoint2 s = {d.p - c.p, d.q - c.q};
    const HuangAxiHeightPoint2 ca = {c.p - a.p, c.q - a.q};
    const double denominator = huang_axi_algorithm5_cross(r, s);
    const double scale = fmax(1., hypot(r.p, r.q)*hypot(s.p, s.q));
    if (!isfinite(denominator) || fabs(denominator) <= tolerance*scale)
        return false;
    const double t = huang_axi_algorithm5_cross(ca, s)/denominator;
    const double u = huang_axi_algorithm5_cross(ca, r)/denominator;
    if (t < -tolerance || t > 1. + tolerance ||
        u < -tolerance || u > 1. + tolerance)
        return false;
    intersection->p = a.p + t*r.p;
    intersection->q = a.q + t*r.q;
    return huang_axi_algorithm5_point_is_finite(*intersection);
}

/* Scaled partial-pivoting solve for a 3 by 3 system. */
static inline bool huang_axi_algorithm5_solve3(
    const double matrix[3][3], const double right_hand_side[3],
    const double tolerance, double solution[3])
{
    double augmented[3][4];
    double matrix_norm = 0.;
    for (int row = 0; row < 3; row++) {
        double row_norm = 0.;
        for (int column = 0; column < 3; column++) {
            augmented[row][column] = matrix[row][column];
            row_norm += fabs(matrix[row][column]);
        }
        augmented[row][3] = right_hand_side[row];
        matrix_norm = fmax(matrix_norm, row_norm);
    }
    if (!isfinite(matrix_norm) || matrix_norm == 0.)
        return false;

    for (int pivot_column = 0; pivot_column < 3; pivot_column++) {
        int pivot_row = pivot_column;
        double pivot_magnitude = fabs(augmented[pivot_row][pivot_column]);
        for (int row = pivot_column + 1; row < 3; row++)
            if (fabs(augmented[row][pivot_column]) > pivot_magnitude) {
                pivot_row = row;
                pivot_magnitude = fabs(augmented[row][pivot_column]);
            }
        if (!isfinite(pivot_magnitude) ||
            pivot_magnitude <= tolerance*fmax(1., matrix_norm))
            return false;
        if (pivot_row != pivot_column)
            for (int column = pivot_column; column < 4; column++) {
                const double temporary = augmented[pivot_column][column];
                augmented[pivot_column][column] = augmented[pivot_row][column];
                augmented[pivot_row][column] = temporary;
            }
        for (int row = pivot_column + 1; row < 3; row++) {
            const double factor = augmented[row][pivot_column]/
                                  augmented[pivot_column][pivot_column];
            augmented[row][pivot_column] = 0.;
            for (int column = pivot_column + 1; column < 4; column++)
                augmented[row][column] -= factor*augmented[pivot_column][column];
        }
    }

    for (int row = 2; row >= 0; row--) {
        double value = augmented[row][3];
        for (int column = row + 1; column < 3; column++)
            value -= augmented[row][column]*solution[column];
        solution[row] = value/augmented[row][row];
        if (!isfinite(solution[row]))
            return false;
    }
    return true;
}

static inline double huang_axi_algorithm5_condition_inf(
    const double matrix[3][3], const double tolerance)
{
    double matrix_norm = 0.;
    for (int row = 0; row < 3; row++) {
        double sum = 0.;
        for (int column = 0; column < 3; column++)
            sum += fabs(matrix[row][column]);
        matrix_norm = fmax(matrix_norm, sum);
    }
    double inverse[3][3];
    for (int column = 0; column < 3; column++) {
        const double unit[3] = {column == 0, column == 1, column == 2};
        double solution[3] = {0., 0., 0.};
        if (!huang_axi_algorithm5_solve3(matrix, unit, tolerance, solution))
            return INFINITY;
        for (int row = 0; row < 3; row++)
            inverse[row][column] = solution[row];
    }
    double inverse_norm = 0.;
    for (int row = 0; row < 3; row++) {
        double sum = 0.;
        for (int column = 0; column < 3; column++)
            sum += fabs(inverse[row][column]);
        inverse_norm = fmax(inverse_norm, sum);
    }
    return matrix_norm*inverse_norm;
}

static inline double huang_axi_algorithm5_polynomial(
    const double coefficients[3], const double q)
{
    return (coefficients[0]*q + coefficients[1])*q + coefficients[2];
}

static inline bool huang_axi_algorithm5_fit_candidate(
    const HuangAxiAlgorithm5Candidate * candidate,
    const HuangAxiAlgorithm5Path path, const double tolerance,
    HuangAxiAlgorithm5Result * result)
{
    if (!candidate || !result || !candidate->accumulation_valid ||
        !isfinite(candidate->contact_slope) ||
        !isfinite(candidate->center_q) || !isfinite(candidate->ghost_q) ||
        !isfinite(candidate->curvature_q) ||
        !huang_axi_algorithm5_point_is_finite(candidate->height_point) ||
        !huang_axi_algorithm5_point_is_finite(
             candidate->neighbor_height_point))
        return false;

    HuangAxiHeightPoint2 contact_point;
    if (!huang_axi_algorithm5_segment_intersection(
            candidate->interface_segment_a, candidate->interface_segment_b,
            candidate->solid_segment_a, candidate->solid_segment_b,
            tolerance, &contact_point))
        return false;

    const double qf = candidate->height_point.q;
    const double qs = contact_point.q;
    const double q1 = candidate->neighbor_height_point.q;
    const double matrix[3][3] = {
        {qf*qf, qf, 1.},
        {2.*qs, 1., 0.},
        {q1*q1, q1, 1.}
    };
    const double right_hand_side[3] = {
        candidate->height_point.p,
        candidate->contact_slope,
        candidate->neighbor_height_point.p
    };
    double coefficients[3] = {0., 0., 0.};
    if (!huang_axi_algorithm5_solve3(matrix, right_hand_side,
                                     tolerance, coefficients))
        return false;

    double residual = 0.;
    for (int row = 0; row < 3; row++) {
        double value = 0.;
        for (int column = 0; column < 3; column++)
            value += matrix[row][column]*coefficients[column];
        residual = fmax(residual, fabs(value - right_hand_side[row]));
    }
    const double condition =
        huang_axi_algorithm5_condition_inf(matrix, tolerance);
    if (!isfinite(condition) || !isfinite(residual))
        return false;

    result->path = path;
    result->contact_point = contact_point;
    for (int i = 0; i < 3; i++)
        result->coefficients[i] = coefficients[i];
    result->matrix_condition = condition;
    result->fit_residual = residual;
    result->height_minus =
        huang_axi_algorithm5_polynomial(coefficients, candidate->ghost_q);
    result->height_center =
        huang_axi_algorithm5_polynomial(coefficients, candidate->center_q);
    result->height_plus = candidate->neighbor_height_point.p;
    result->height_gradient = 2.*coefficients[0]*candidate->center_q +
                              coefficients[1];
    const double denominator = pow(1. + result->height_gradient*
                                         result->height_gradient, 1.5);
    result->meridional_curvature = 2.*coefficients[0]/denominator;

    /* For x = h(r), the second principal curvature is
     * h'(r)/(r*sqrt(1+h'(r)^2)).  At r = 0 regularity requires h'(0)=0,
     * and l'Hopital's rule gives h''(0).  No arbitrary 2*pi or 1/r patch
     * is introduced.  Native Basilisk curvature remains authoritative when
     * these corrected heights are installed in a Basilisk height field. */
    const double radial_denominator =
        sqrt(1. + result->height_gradient*result->height_gradient);
    double azimuthal_curvature;
    if (fabs(candidate->curvature_q) > tolerance)
        azimuthal_curvature = result->height_gradient/
                              (candidate->curvature_q*radial_denominator);
    else if (fabs(result->height_gradient) <=
             sqrt(tolerance)*fmax(1., fabs(coefficients[0])))
        azimuthal_curvature = 2.*coefficients[0];
    else
        return false;
    result->axisymmetric_curvature = result->meridional_curvature +
                                     azimuthal_curvature;
    if (!isfinite(result->height_minus) ||
        !isfinite(result->height_center) ||
        !isfinite(result->height_plus) ||
        !isfinite(result->height_gradient) ||
        !isfinite(result->meridional_curvature) ||
        !isfinite(result->axisymmetric_curvature))
        return false;
    result->status = HUANG_AXI_ALG5_OK;
    return true;
}

/*
 * Huang Algorithm 5 branch order.  primary and secondary represent the two
 * ordinary height directions.  reconstructed is the Algorithm-3-normal
 * interface path used when both ordinary 1x9 accumulations/intersections
 * fail.  Extended fractions are diagnostic stencil data only and are used
 * last; physical_liquid_measure is copied unchanged in every branch.
 */
static inline HuangAxiAlgorithm5Result huang_axi_algorithm5_compute(
    const HuangAxiAlgorithm5Input * input, const double tolerance)
{
    HuangAxiAlgorithm5Result result = {
        .status = HUANG_AXI_ALG5_INVALID_INPUT,
        .path = HUANG_AXI_ALG5_INVALID,
        .contact_point = {NAN, NAN},
        .coefficients = {NAN, NAN, NAN},
        .matrix_condition = NAN,
        .fit_residual = NAN,
        .height_minus = NAN,
        .height_center = NAN,
        .height_plus = NAN,
        .height_gradient = NAN,
        .meridional_curvature = NAN,
        .axisymmetric_curvature = NAN,
        .physical_liquid_measure = NAN
    };
    if (!input || !isfinite(tolerance) || tolerance <= 0. ||
        !isfinite(input->physical_liquid_measure) ||
        input->physical_liquid_measure < 0.)
        return result;
    result.physical_liquid_measure = input->physical_liquid_measure;

    if (huang_axi_algorithm5_fit_candidate(
            &input->primary, HUANG_AXI_ALG5_PRIMARY_HF,
            tolerance, &result) ||
        huang_axi_algorithm5_fit_candidate(
            &input->secondary, HUANG_AXI_ALG5_SECONDARY_HF,
            tolerance, &result) ||
        huang_axi_algorithm5_fit_candidate(
            &input->reconstructed,
            HUANG_AXI_ALG5_RECONSTRUCTED_INTERFACE_HF,
            tolerance, &result))
        return result;

    if (input->extended_fraction_valid &&
        isfinite(input->extended_height_minus) &&
        isfinite(input->extended_height_center) &&
        isfinite(input->extended_height_plus)) {
        result.status = HUANG_AXI_ALG5_OK;
        result.path = HUANG_AXI_ALG5_EXTENDED_FRACTION_FALLBACK;
        result.height_minus = input->extended_height_minus;
        result.height_center = input->extended_height_center;
        result.height_plus = input->extended_height_plus;
        result.height_gradient =
            (input->extended_height_plus - input->extended_height_minus)/2.;
        result.meridional_curvature =
            (input->extended_height_plus + input->extended_height_minus -
             2.*input->extended_height_center)/
            pow(1. + result.height_gradient*result.height_gradient, 1.5);
        result.axisymmetric_curvature = result.meridional_curvature;
        return result;
    }

    result.status = HUANG_AXI_ALG5_NO_VALID_PATH;
    return result;
}

#ifdef HSHIFT
static inline double huang_axi_algorithm5_encode_height(
    const double value, const int orientation_value)
{
    return value + (orientation_value ? HSHIFT : 0.);
}

/* Install one corrected value at the current cell of a native Basilisk
 * height vector.  Neighbor writes are deliberately forbidden: every cell
 * owns its own value, which keeps foreach() updates race-free under MPI and
 * OpenMP.  The caller selects result.height_minus/center/plus for the cell
 * being visited.  These two functions are explicit because the target is
 * strictly 2D. */
static inline bool huang_axi_algorithm5_install_x_current(
    Point point, vector height_field, const double corrected_value,
    const int orientation_value)
{
    if (!isfinite(corrected_value) ||
        (orientation_value != 0 && orientation_value != 1))
        return false;
    height_field.x[] = huang_axi_algorithm5_encode_height(
        corrected_value, orientation_value);
    return true;
}

static inline bool huang_axi_algorithm5_install_y_current(
    Point point, vector height_field, const double corrected_value,
    const int orientation_value)
{
    if (!isfinite(corrected_value) ||
        (orientation_value != 0 && orientation_value != 1))
        return false;
    height_field.y[] = huang_axi_algorithm5_encode_height(
        corrected_value, orientation_value);
    return true;
}
#endif

#endif
