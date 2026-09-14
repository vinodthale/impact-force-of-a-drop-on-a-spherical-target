#ifndef STAGE_B_GEOMETRY_H
#define STAGE_B_GEOMETRY_H

#include <stdbool.h>

typedef struct {
  double x;
  double y;
} StageBVec2;

typedef struct {
  int x;
  int y;
} StageBIndex2;

typedef struct {
  double c;
  double cs;
} StageBFractionCell;

typedef struct {
  double a0;
  double a1;
  double a2;
} StageBParabola;

typedef enum {
  STAGE_B_HEIGHT_OK = 0,
  STAGE_B_HEIGHT_INVALID_INPUT,
  STAGE_B_HEIGHT_NO_FULL_LIQUID_START,
  STAGE_B_HEIGHT_NO_EMPTY_END
} StageBHeightStatus;

typedef enum {
  STAGE_B_PARABOLA_OK = 0,
  STAGE_B_PARABOLA_INVALID_INPUT,
  STAGE_B_PARABOLA_SINGULAR
} StageBParabolaStatus;

typedef enum {
  STAGE_B_HEIGHT_AXIS_X = 0,
  STAGE_B_HEIGHT_AXIS_Y
} StageBHeightAxis;

typedef enum {
  STAGE_B_INTERSECTION_OK = 0,
  STAGE_B_INTERSECTION_INVALID_INPUT,
  STAGE_B_INTERSECTION_PARALLEL,
  STAGE_B_INTERSECTION_OUTSIDE_CELL
} StageBIntersectionStatus;

bool stage_b_complete_fraction (double f, double cs, double * c_complete);

StageBVec2 stage_b_contact_normal_2d (StageBVec2 ns, StageBVec2 nf,
                                      double angle_radians);

StageBIndex2 stage_b_algorithm4_direction (StageBVec2 ns, StageBVec2 nf,
                                            double angle_radians);

bool stage_b_algorithm4_is_contact (
  double c_center, double cs_center,
  const StageBFractionCell neighbors[3], double fraction_tolerance);

double stage_b_line_area_2d (StageBVec2 normal, double alpha);

bool stage_b_extend_line_3x3 (StageBVec2 normal, double alpha,
                              double ce[3][3]);

StageBHeightStatus stage_b_extended_height_9 (
  const double ce[9], int direction, double delta,
  double fraction_tolerance, double * height_value);

StageBHeightStatus stage_b_contact_height_9 (
  const double c_complete[9], const double cs[9], int direction,
  double delta, double fraction_tolerance, double * height_value);

StageBIntersectionStatus stage_b_height_solid_intersection (
  StageBHeightAxis height_axis, double height_coordinate,
  StageBVec2 solid_normal, double solid_alpha,
  double parallel_tolerance, double * transverse_coordinate);

StageBParabolaStatus stage_b_algorithm5_parabola (
  double hf, double yf, double h1, double y1,
  double contact_slope, double ys, double singular_tolerance,
  StageBParabola * fit);

double stage_b_parabola_value (StageBParabola fit, double coordinate);

double stage_b_height_curvature (double h_minus, double h_center,
                                 double h_plus, double delta);

#endif
