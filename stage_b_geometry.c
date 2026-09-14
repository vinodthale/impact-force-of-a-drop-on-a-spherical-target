#include "stage_b_geometry.h"

#include <math.h>
#include <stddef.h>

/*
 * This is a 2D, fixed-stencil geometry kernel only.  It has no transported
 * field, no grid mutation, no metric operation and no Navier--Stokes event.
 *
 * Source map:
 *   complete fraction: Huang c=f*cs bridge for 0<=c<=cs
 *   contact normal:    TPR2D.h normal_contact(), lines 899--913
 *   contact marking:   tmp_fraction_field.h sort_cell(), lines 441--512
 *   line fraction:     Basilisk geometry.h line_area(), lines 163--202
 *   extended height:   Huang paper Eq. (13), nine-cell maximum stencil
 *   parabola:          Huang paper Eqs. (25)--(26), Algorithm 5
 *   curvature:         Huang paper Eq. (14)
 *   fallback order:    Huang paper Fig. 8 cases 1--3
 */

static int stage_b_sign2 (double value)
{
  return value > 0. ? 1 : value < 0. ? -1 : 0;
}

static bool stage_b_unit_interval (double value)
{
  return isfinite (value) && value >= 0. && value <= 1.;
}

bool stage_b_complete_fraction (double f, double cs, double * c_complete)
{
  if (c_complete == NULL)
    return false;
  if (!stage_b_unit_interval (f) || !stage_b_unit_interval (cs)) {
    *c_complete = NAN;
    return false;
  }
  *c_complete = f*cs;
  return true;
}

StageBVec2 stage_b_contact_normal_2d (StageBVec2 ns, StageBVec2 nf,
                                      double angle_radians)
{
  StageBVec2 normal = {NAN, NAN};
  if (!isfinite (ns.x) || !isfinite (ns.y) ||
      !isfinite (nf.x) || !isfinite (nf.y) ||
      !isfinite (angle_radians))
    return normal;

  if (-ns.x*nf.y + ns.y*nf.x > 0.) {
    normal.x = -ns.x*cos (angle_radians) + ns.y*sin (angle_radians);
    normal.y = -ns.x*sin (angle_radians) - ns.y*cos (angle_radians);
  }
  else {
    normal.x = -ns.x*cos (angle_radians) - ns.y*sin (angle_radians);
    normal.y =  ns.x*sin (angle_radians) - ns.y*cos (angle_radians);
  }
  return normal;
}

StageBIndex2 stage_b_algorithm4_direction (StageBVec2 ns, StageBVec2 nf,
                                            double angle_radians)
{
  StageBVec2 nc = stage_b_contact_normal_2d (ns, nf, angle_radians);
  if (!isfinite (nc.x) || !isfinite (nc.y))
    return (StageBIndex2){0, 0};

  StageBVec2 tangent;
  if (ns.x*nc.y - ns.y*nc.x > 0.) {
    tangent.x =  nc.y;
    tangent.y = -nc.x;
  }
  else {
    tangent.x = -nc.y;
    tangent.y =  nc.x;
  }
  return (StageBIndex2){stage_b_sign2 (tangent.x),
                        stage_b_sign2 (tangent.y)};
}

bool stage_b_algorithm4_is_contact (
  double c_center, double cs_center,
  const StageBFractionCell neighbors[3], double fraction_tolerance)
{
  if (neighbors == NULL || !isfinite (fraction_tolerance) ||
      fraction_tolerance < 0. ||
      !stage_b_unit_interval (cs_center) || !isfinite (c_center) ||
      c_center < 0. || c_center > cs_center)
    return false;

  if (!(cs_center > 0. && cs_center < 1. &&
        c_center > fraction_tolerance &&
        c_center < cs_center - fraction_tolerance))
    return false;

  for (int k = 0; k < 3; k++) {
    double c = neighbors[k].c;
    double cs = neighbors[k].cs;
    if (!stage_b_unit_interval (cs) || !isfinite (c) ||
        c < 0. || c > cs)
      return false;
    if (cs > 0. && cs < 1. &&
        c >= fraction_tolerance && c <= cs - fraction_tolerance)
      return false;
  }
  return true;
}

double stage_b_line_area_2d (StageBVec2 normal, double alpha)
{
  double nx = normal.x, ny = normal.y;
  if (!isfinite (nx) || !isfinite (ny) || !isfinite (alpha) ||
      fabs (nx) + fabs (ny) == 0.)
    return NAN;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = -nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = -ny;
  }

  if (alpha <= 0.)
    return 0.;
  if (alpha >= nx + ny)
    return 1.;

  double area;
  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    double volume = alpha*alpha;
    double offset = alpha - nx;
    if (offset > 0.)
      volume -= offset*offset;
    offset = alpha - ny;
    if (offset > 0.)
      volume -= offset*offset;
    area = volume/(2.*nx*ny);
  }

  /* This bound is geometric roundoff protection for auxiliary ce only. */
  return area < 0. ? 0. : area > 1. ? 1. : area;
}

bool stage_b_extend_line_3x3 (StageBVec2 normal, double alpha,
                              double ce[3][3])
{
  if (ce == NULL || !isfinite (normal.x) || !isfinite (normal.y) ||
      !isfinite (alpha) || fabs (normal.x) + fabs (normal.y) == 0.)
    return false;

  for (int j = -1; j <= 1; j++)
    for (int i = -1; i <= 1; i++) {
      double local_alpha = alpha - normal.x*i - normal.y*j;
      double fraction = stage_b_line_area_2d (normal, local_alpha);
      if (!stage_b_unit_interval (fraction))
        return false;
      ce[j + 1][i + 1] = fraction;
    }
  return true;
}

StageBHeightStatus stage_b_extended_height_9 (
  const double ce[9], int direction, double delta,
  double fraction_tolerance, double * height_value)
{
  if (height_value == NULL)
    return STAGE_B_HEIGHT_INVALID_INPUT;
  *height_value = NAN;
  if (ce == NULL || (direction != -1 && direction != 1) ||
      !isfinite (delta) || delta <= 0. ||
      !isfinite (fraction_tolerance) || fraction_tolerance < 0.)
    return STAGE_B_HEIGHT_INVALID_INPUT;

  double sum = 0.;
  for (int i = 0; i < 9; i++) {
    if (!stage_b_unit_interval (ce[i]))
      return STAGE_B_HEIGHT_INVALID_INPUT;
    sum += ce[i];
  }

  int start = direction > 0 ? 0 : 8;
  int stop = direction > 0 ? 9 : -1;
  int step = direction;
  bool full_start = false, empty_end = false;
  for (int i = start; i != stop; i += step) {
    if (!full_start) {
      if (ce[i] >= 1. - fraction_tolerance)
        full_start = true;
    }
    else if (ce[i] <= fraction_tolerance) {
      empty_end = true;
      break;
    }
  }
  if (!full_start)
    return STAGE_B_HEIGHT_NO_FULL_LIQUID_START;
  if (!empty_end)
    return STAGE_B_HEIGHT_NO_EMPTY_END;

  *height_value = sum*delta;
  return STAGE_B_HEIGHT_OK;
}

StageBHeightStatus stage_b_contact_height_9 (
  const double c_complete[9], const double cs[9], int direction,
  double delta, double fraction_tolerance, double * height_value)
{
  if (height_value == NULL)
    return STAGE_B_HEIGHT_INVALID_INPUT;
  *height_value = NAN;
  if (c_complete == NULL || cs == NULL ||
      (direction != -1 && direction != 1) ||
      !isfinite (delta) || delta <= 0. ||
      !isfinite (fraction_tolerance) || fraction_tolerance < 0.)
    return STAGE_B_HEIGHT_INVALID_INPUT;

  double sum = 0.;
  for (int i = 0; i < 9; i++) {
    if (!stage_b_unit_interval (cs[i]) || !isfinite (c_complete[i]) ||
        c_complete[i] < 0. || c_complete[i] > cs[i])
      return STAGE_B_HEIGHT_INVALID_INPUT;
    sum += c_complete[i];
  }

  int start = direction > 0 ? 0 : 8;
  int stop = direction > 0 ? 9 : -1;
  int step = direction;
  bool full_start = false, empty_end = false;
  for (int i = start; i != stop; i += step) {
    if (!full_start) {
      if (cs[i] > fraction_tolerance &&
          c_complete[i] >= cs[i] - fraction_tolerance)
        full_start = true;
    }
    else if (cs[i] > fraction_tolerance &&
             c_complete[i] <= fraction_tolerance) {
      empty_end = true;
      break;
    }
  }
  if (!full_start)
    return STAGE_B_HEIGHT_NO_FULL_LIQUID_START;
  if (!empty_end)
    return STAGE_B_HEIGHT_NO_EMPTY_END;

  *height_value = sum*delta;
  return STAGE_B_HEIGHT_OK;
}

StageBIntersectionStatus stage_b_height_solid_intersection (
  StageBHeightAxis height_axis, double height_coordinate,
  StageBVec2 solid_normal, double solid_alpha,
  double parallel_tolerance, double * transverse_coordinate)
{
  if (transverse_coordinate == NULL)
    return STAGE_B_INTERSECTION_INVALID_INPUT;
  *transverse_coordinate = NAN;
  if ((height_axis != STAGE_B_HEIGHT_AXIS_X &&
       height_axis != STAGE_B_HEIGHT_AXIS_Y) ||
      !isfinite (height_coordinate) ||
      !isfinite (solid_normal.x) || !isfinite (solid_normal.y) ||
      !isfinite (solid_alpha) || !isfinite (parallel_tolerance) ||
      parallel_tolerance < 0. ||
      fabs (solid_normal.x) + fabs (solid_normal.y) == 0.)
    return STAGE_B_INTERSECTION_INVALID_INPUT;

  double denominator, numerator;
  if (height_axis == STAGE_B_HEIGHT_AXIS_X) {
    denominator = solid_normal.y;
    numerator = solid_alpha - solid_normal.x*height_coordinate;
  }
  else {
    denominator = solid_normal.x;
    numerator = solid_alpha - solid_normal.y*height_coordinate;
  }

  if (fabs (denominator) <= parallel_tolerance)
    return STAGE_B_INTERSECTION_PARALLEL;

  double coordinate = numerator/denominator;
  *transverse_coordinate = coordinate;
  if (coordinate < -0.5 - parallel_tolerance ||
      coordinate > 0.5 + parallel_tolerance)
    return STAGE_B_INTERSECTION_OUTSIDE_CELL;
  return STAGE_B_INTERSECTION_OK;
}

StageBParabolaStatus stage_b_algorithm5_parabola (
  double hf, double yf, double h1, double y1,
  double contact_slope, double ys, double singular_tolerance,
  StageBParabola * fit)
{
  if (fit == NULL)
    return STAGE_B_PARABOLA_INVALID_INPUT;
  fit->a0 = fit->a1 = fit->a2 = NAN;
  if (!isfinite (hf) || !isfinite (yf) ||
      !isfinite (h1) || !isfinite (y1) ||
      !isfinite (contact_slope) || !isfinite (ys) ||
      !isfinite (singular_tolerance) || singular_tolerance < 0.)
    return STAGE_B_PARABOLA_INVALID_INPUT;

  double dy = y1 - yf;
  double denominator = dy*(y1 + yf - 2.*ys);
  if (fabs (denominator) <= singular_tolerance)
    return STAGE_B_PARABOLA_SINGULAR;

  fit->a0 = (h1 - hf - contact_slope*dy)/denominator;
  fit->a1 = contact_slope - 2.*fit->a0*ys;
  fit->a2 = hf - fit->a0*yf*yf - fit->a1*yf;
  return STAGE_B_PARABOLA_OK;
}

double stage_b_parabola_value (StageBParabola fit, double coordinate)
{
  if (!isfinite (fit.a0) || !isfinite (fit.a1) ||
      !isfinite (fit.a2) || !isfinite (coordinate))
    return NAN;
  return fit.a0*coordinate*coordinate + fit.a1*coordinate + fit.a2;
}

double stage_b_height_curvature (double h_minus, double h_center,
                                 double h_plus, double delta)
{
  if (!isfinite (h_minus) || !isfinite (h_center) ||
      !isfinite (h_plus) || !isfinite (delta) || delta <= 0.)
    return NAN;
  double first = (h_plus - h_minus)/(2.*delta);
  double second = (h_plus - 2.*h_center + h_minus)/(delta*delta);
  return second/pow (1. + first*first, 1.5);
}
