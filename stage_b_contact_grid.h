#ifndef STAGE_B_CONTACT_GRID_H
#define STAGE_B_CONTACT_GRID_H

#if dimension != 2
# error "Stage B contact geometry is strictly two-dimensional"
#endif

#include "curvature.h"
#include "stage_b_geometry.h"
#include "stage_b_geometry.c"
#include "huang_axi_locked_include/huang_axi_algorithm4_contact_cell_detection.h"
#include "huang_axi_locked_include/huang_axi_algorithm5_embedded_height_function.h"
#include "stage_b_huang_embed_curvature.h"
#ifndef VFTL
/* Required by Huang's released tmp_fraction_field.h.  This is the same
 * tolerance default provided by the released embed_height_normal.h. */
# define VFTL 1e-10
#endif
#include "huang_original/embed_correct_height.h"

/* Huang's marked height constructor is source-compatible with the released
 * EBM_VOF path but has the same public helper names as Basilisk's native
 * heights.h.  Prefix the helpers so both owners remain explicit in this
 * isolated geometry stage. */
#define half_column huang_half_column
#define column_propagation huang_column_propagation
#define refine_h_x huang_refine_h_x
#define heights huang_heights
#define height huang_height
#define orientation huang_orientation
#include "huang_original/embed_heights.h"
#undef orientation
#undef height
#undef heights
#undef refine_h_x
#undef column_propagation
#undef half_column

scalar stage_b_c_complete[], stage_b_ce[];
scalar stage_b_mark[];
scalar stage_b_angle[], stage_b_alphacs[];
scalar stage_b_contact[], stage_b_support[], stage_b_force_support[];
scalar stage_b_nx[], stage_b_ny[], stage_b_plic_alpha[];
scalar stage_b_contact_x[], stage_b_contact_y[];
scalar stage_b_contact_path_valid[];
vector stage_b_ce_height[];
vector stage_b_mc[], stage_b_ms[], stage_b_oxyi[];
vector stage_b_height_before[];
scalar stage_b_algorithm5_active[];
scalar stage_b_algorithm5_minus[];
scalar stage_b_algorithm5_center[];
scalar stage_b_algorithm5_plus[];
scalar stage_b_algorithm5_orientation[];

typedef struct {
  long mixed_cut_cells;
  long algorithm4_contact_cells;
  long open_polygon_valid;
  long liquid_polygon_valid;
  long contact_intersections;
  long contact_cells;
  long extended_cells;
  long algorithm5_marked_contacts;
  long algorithm5_height_valid;
  long algorithm5_height_installs;
  long invalid_cells;
  double ce_min;
  double ce_max;
  double open_cell_ce_error;
} StageBGridStats;

static StageBGridStats stage_b_grid_stats;

static inline bool stage_b_fraction_valid (double value)
{
  return isfinite (value) && value >= -1e-12 && value <= 1. + 1e-12;
}

static inline bool stage_b_interfacial_w (Point point, scalar solid_fraction)
{
  if (solid_fraction[] >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
        if (solid_fraction[i] <= 0.)
          return true;
  }
  else if (solid_fraction[] <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
        if (solid_fraction[i] >= 1.)
          return true;
  }
  else
    return true;
  return false;
}

/* Direct transcription of Huang's contact_embed_cell() predicate. */
static inline bool stage_b_contact_embed_cell (
  Point point, scalar complete_fraction, scalar solid_fraction,
  scalar contact_mark, double angle_radians)
{
  const double vftl = 1e-10;
  if (contact_mark[] != 10.)
    return false;

  if (complete_fraction[] > vftl &&
      complete_fraction[] < solid_fraction[] - vftl) {
    if (angle_radians <= pi/2.) {
      for (int i = -1; i <= 1; i += 2) {
        foreach_dimension()
          if (contact_mark[i] == 10. && complete_fraction[i] <= vftl)
            return true;
        if ((contact_mark[i,i] == 10. &&
             complete_fraction[i,i] <= vftl) ||
            (contact_mark[i,-i] == 10. &&
             complete_fraction[i,-i] <= vftl))
          return true;
      }
    }
    else {
      for (int i = -1; i <= 1; i += 2) {
        foreach_dimension()
          if (contact_mark[i] == 10. &&
              complete_fraction[i] >= solid_fraction[i] - vftl)
            return true;
        if ((contact_mark[i,i] == 10. &&
             complete_fraction[i,i] >= solid_fraction[i,i] - vftl) ||
            (contact_mark[i,-i] == 10. &&
             complete_fraction[i,-i] >= solid_fraction[i,-i] - vftl))
          return true;
      }
    }
  }
  else if (complete_fraction[] >= solid_fraction[] - vftl) {
    for (int i = -1; i <= 1; i += 2) {
      foreach_dimension()
        if (contact_mark[i] == 10. && complete_fraction[i] <= vftl)
          return true;
    }
  }
  else if (complete_fraction[] <= vftl) {
    for (int i = -1; i <= 1; i += 2) {
      foreach_dimension()
        if (contact_mark[i] == 10. &&
            complete_fraction[i] >= solid_fraction[i] - vftl)
          return true;
    }
  }
  return false;
}

static inline bool stage_b_contact_candidate (
  Point point, double tolerance, double angle_radians, scalar contact_mark)
{
  return cs[] > tolerance &&
         stage_b_interfacial_w (point, cs) &&
         stage_b_contact_embed_cell (point, stage_b_c_complete, cs,
                                     contact_mark, angle_radians);
}

static inline bool stage_b_interfacial_i (Point point, scalar liquid_fraction)
{
  if (liquid_fraction[] >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
        if (liquid_fraction[i] <= 0.)
          return true;
  }
  else if (liquid_fraction[] <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      foreach_dimension()
        if (liquid_fraction[i] >= 1.)
          return true;
  }
  else
    return true;
  return false;
}

/* Huang sort_cell() geometry only.  It produces marks for the auxiliary
 * reconstruction; it never writes the physical VOF field. */
static void stage_b_sort_contact_marks (scalar f, scalar c_complete,
                                        scalar mark,
                                        scalar contact_mark,
                                        double angle_radians)
{
  const double vftl = 1e-10;
  foreach() {
    mark[] = 0.;
    if (cs[] <= 0.)
      mark[] = 1.;
    else if (stage_b_interfacial_w (point, cs)) {
      /* Released sort_cell(): pure cut cells are classified first.  The
       * contact-line stencil is evaluated only for mixed cells. */
      if (stage_b_c_complete[] <= vftl)
        mark[] = 2.;
      else if (stage_b_c_complete[] >= cs[] - vftl)
        mark[] = 3.;
      else {
        mark[] = 5.;
        coord ns = facet_normal (point, cs, fs);
        coord nf = mycs (point, f);
        coord nc = normal_contact (ns, nf, angle_radians);
        coord mnc = {0, 0};
        if (ns.x*nc.y - ns.y*nc.x > 0.) {
          mnc.x = nc.y;
          mnc.y = -nc.x;
        }
        else {
          mnc.x = -nc.y;
          mnc.y = nc.x;
        }
        int a = sign2 (mnc.x), b = sign2 (mnc.y);
        int markk4 = 1;
        if (cs[a,0] > 0. && cs[a,0] < 1. &&
            stage_b_c_complete[a,0] >= vftl &&
            stage_b_c_complete[a,0] <= cs[a,0] - vftl)
          markk4 = 0;
        if (cs[0,b] > 0. && cs[0,b] < 1. &&
            stage_b_c_complete[0,b] >= vftl &&
            stage_b_c_complete[0,b] <= cs[0,b] - vftl)
          markk4 = 0;
        if (cs[a,b] > 0. && cs[a,b] < 1. &&
            stage_b_c_complete[a,b] >= vftl &&
            stage_b_c_complete[a,b] <= cs[a,b] - vftl)
          markk4 = 0;
        if ((cs[a,0] <= 0. || cs[a,0] >= 1.) &&
            (cs[0,b] <= 0. || cs[0,b] >= 1.) &&
            (cs[a,b] <= 0. || cs[a,b] >= 1.))
          markk4 = 1;
        if (markk4)
          mark[] = 4.;
      }
    }
    else if (cs[] >= 1.) {
      if (stage_b_interfacial_i (point, c_complete))
        mark[] = 6.;
      else if (c_complete[] <= 0.)
        mark[] = 7.;
      else if (c_complete[] >= 1.)
        mark[] = 8.;
    }
  }
  boundary ({mark});
}

static inline bool stage_b_build_contact_line (Point point, scalar f,
                                                scalar c_complete,
                                                double tolerance,
                                                double angle_radians,
                                                coord * liquid_normal,
                                                double * liquid_alpha,
                                                coord * contact_point,
                                                int * completed_stage)
{
  *completed_stage = 0;
  coord ns = facet_normal (point, cs, fs);
  coord nf = mycs (point, f);
  coord nc = normal_contact (ns, nf, angle_radians);
  double normal_length = hypot (nc.x, nc.y);
#ifdef STAGE_B_DEBUG_INTERSECTION
  fprintf (stderr, "STAGE_B_CONTACT_ENTRY x=%.17g y=%.17g cs=%.17g "
           "c=%.17g ns=(%.17g,%.17g) nf=(%.17g,%.17g) "
           "nc=(%.17g,%.17g)\n", x, y, cs[], stage_b_c_complete[],
           ns.x, ns.y, nf.x, nf.y, nc.x, nc.y);
#endif
  if (!isfinite (normal_length) || normal_length <= tolerance)
    return false;

  HuangAxiAlgorithm4Stencil stencil = {0};
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++) {
      HuangAxiAlgorithm4Cell * stencil_cell =
        &stencil.cells[i + 1][j + 1];
      stencil_cell->available = true;
      stencil_cell->cs = cs[i,j];
      stencil_cell->c = stage_b_c_complete[i,j];
    }
  HuangAxiAlgorithm3Result normal = {
    .status = HUANG_AXI_ALG3_OK,
    .solid_normal = {ns.x, ns.y},
    .liquid_normal = {nc.x, nc.y}
  };
  HuangAxiAlgorithm4Result classification =
    huang_axi_algorithm4_detect (&stencil, normal, tolerance);
  const bool mixed_interface =
    cs[] > tolerance && cs[] < 1. - tolerance &&
    stage_b_c_complete[] > tolerance &&
    stage_b_c_complete[] < cs[] - tolerance;
#ifdef STAGE_B_DEBUG_INTERSECTION
  fprintf (stderr, "STAGE_B_CONTACT_CLASS x=%.17g y=%.17g mixed=%d "
           "is_contact=%d reason=%d\n", x, y, mixed_interface,
           classification.is_contact_line, classification.reason);
#endif
  if (mixed_interface && !classification.is_contact_line)
    return false;
  *completed_stage = 1;

  /*
   * This is the released Huang contact-cell construction: polygon_alpha()
   * reconstructs the liquid interface in the open EMBED polygon using the
   * complete-cell fraction c_complete, the solid facet normal and plane
   * alpha.  In the native Basilisk harness c_complete is the read-only
   * derived field f*cs; physical f is open-fluid-relative.
   * It is deliberately not replaced by a second square-cell PLIC solve.
   */
  coord p_mof[2] = {{nodata, nodata}, {nodata, nodata}};
  coord polygon[5] = {{nodata, nodata}, {nodata, nodata},
                      {nodata, nodata}, {nodata, nodata},
                      {nodata, nodata}};
  double solid_alpha = plane_alpha (cs[], ns);
  double plic_alpha = polygon_alpha (stage_b_c_complete[], nc, ns,
                                     solid_alpha, p_mof, polygon);
#ifdef STAGE_B_DEBUG_INTERSECTION
  fprintf (stderr, "STAGE_B_CONTACT_PLIC x=%.17g y=%.17g alpha=%.17g "
           "p0=(%.17g,%.17g) p1=(%.17g,%.17g)\n", x, y, plic_alpha,
           p_mof[0].x, p_mof[0].y, p_mof[1].x, p_mof[1].y);
#endif
  if (!isfinite (plic_alpha))
    return false;

  liquid_normal->x = nc.x/normal_length;
  liquid_normal->y = nc.y/normal_length;
  *liquid_alpha = plic_alpha/normal_length;

  coord solid_segment[2] = {{nodata, nodata}, {nodata, nodata}};
  int nsolid = myfacets (ns, solid_alpha, solid_segment);
  HuangAxiHeightPoint2 intersection;
  bool intersects = nsolid == 2 &&
    huang_axi_algorithm5_segment_intersection (
      (HuangAxiHeightPoint2){solid_segment[0].x, solid_segment[0].y},
      (HuangAxiHeightPoint2){solid_segment[1].x, solid_segment[1].y},
      (HuangAxiHeightPoint2){p_mof[0].x, p_mof[0].y},
      (HuangAxiHeightPoint2){p_mof[1].x, p_mof[1].y},
      tolerance, &intersection);
#ifdef STAGE_B_DEBUG_INTERSECTION
  fprintf (stderr,
           "STAGE_B_INTERSECTION x=%.17g y=%.17g cs=%.17g c=%.17g "
           "ns=(%.17g,%.17g) nc=(%.17g,%.17g) solid_alpha=%.17g "
           "nsolid=%d solid0=(%.17g,%.17g) solid1=(%.17g,%.17g) "
           "plic0=(%.17g,%.17g) plic1=(%.17g,%.17g) intersects=%d\n",
           x, y, cs[], stage_b_c_complete[], ns.x, ns.y, nc.x, nc.y,
           solid_alpha, nsolid,
           solid_segment[0].x, solid_segment[0].y,
           solid_segment[1].x, solid_segment[1].y,
           p_mof[0].x, p_mof[0].y, p_mof[1].x, p_mof[1].y,
           intersects);
#endif
  contact_point->x = intersects ? intersection.p : nodata;
  contact_point->y = intersects ? intersection.q : nodata;
  *completed_stage = intersects ? 4 : 3;
  return true;
}

static void stage_b_build_auxiliary_fraction (scalar f,
                                               double angle_radians)
{
  const double tolerance = 1e-12;
  long mixed_cut_cells = 0, contact_cells = 0, extended_cells = 0;
  long algorithm4_contact_cells = 0, open_polygon_valid = 0;
  long liquid_polygon_valid = 0, contact_intersections = 0;
  scalar contact_mark[];
  memset (&stage_b_grid_stats, 0, sizeof (stage_b_grid_stats));
  stage_b_grid_stats.ce_min = HUGE;
  stage_b_grid_stats.ce_max = -HUGE;

  boundary ({f, cs});
  foreach() {
    /* Native Basilisk f is relative to the open-fluid volume.  Huang's
     * released embedded reconstruction requires the complete-cell fraction
     * c in [0,cs].  This is a rebuilt read-only bridge, never a transport
     * field and never a replacement for physical f. */
    stage_b_c_complete[] = f[]*cs[];
    contact_mark[] = stage_b_interfacial_w (point, cs) && cs[] > 0. ? 10. : 0.;
    stage_b_mark[] = 0.;
    stage_b_ce[] = cs[] >= 1. - tolerance ? f[] : 0.;
    stage_b_contact[] = stage_b_support[] = stage_b_force_support[] = 0.;
    stage_b_nx[] = stage_b_ny[] = stage_b_plic_alpha[] = nodata;
    stage_b_contact_x[] = stage_b_contact_y[] = nodata;
    stage_b_contact_path_valid[] = 0.;
  }
  boundary ({stage_b_c_complete, contact_mark, stage_b_mark});

  reconstruction_cs (cs, fs, stage_b_ms, stage_b_alphacs);
  reconstruction_mc_myc (stage_b_c_complete, stage_b_mc);
  stage_b_sort_contact_marks (f, stage_b_c_complete,
                              stage_b_mark, contact_mark,
                              angle_radians);
  foreach()
    stage_b_angle[] = angle_radians*180./pi;
  boundary ({stage_b_angle});
  reconstruction_tmp_embed_fraction_field (
    stage_b_c_complete, cs, stage_b_mc, stage_b_ms, stage_b_alphacs, stage_b_angle,
    stage_b_ce, stage_b_mark);

  foreach (reduction(+:mixed_cut_cells) reduction(+:contact_cells)
           reduction(+:algorithm4_contact_cells)
           reduction(+:open_polygon_valid)
           reduction(+:liquid_polygon_valid)
           reduction(+:contact_intersections)) {
    if (stage_b_mark[] != 4.)
      continue;
    mixed_cut_cells += cs[] > tolerance && cs[] < 1. - tolerance;
    coord normal, contact;
    double plic_alpha;
    int completed_stage = 0;
    if (stage_b_build_contact_line (point, f, stage_b_c_complete,
                                    tolerance, angle_radians,
                                    &normal, &plic_alpha, &contact,
                                    &completed_stage)) {
      stage_b_contact[] = 1.;
      stage_b_nx[] = normal.x;
      stage_b_ny[] = normal.y;
      stage_b_plic_alpha[] = plic_alpha;
      stage_b_contact_x[] = contact.x;
      stage_b_contact_y[] = contact.y;
      contact_cells++;
    }
    algorithm4_contact_cells += completed_stage >= 1;
    open_polygon_valid += completed_stage >= 2;
    liquid_polygon_valid += completed_stage >= 3;
    contact_intersections += completed_stage >= 4;
  }
  boundary ({stage_b_contact, stage_b_nx, stage_b_ny,
             stage_b_plic_alpha, stage_b_contact_x, stage_b_contact_y});

  foreach (reduction(+:extended_cells))
    if (stage_b_mark[] == 4.) {
      stage_b_support[] = 1.;
      extended_cells += cs[] < 1. - tolerance;
    }
  boundary ({stage_b_ce, stage_b_support});
  stage_b_grid_stats.mixed_cut_cells = mixed_cut_cells;
  stage_b_grid_stats.algorithm4_contact_cells = algorithm4_contact_cells;
  stage_b_grid_stats.open_polygon_valid = open_polygon_valid;
  stage_b_grid_stats.liquid_polygon_valid = liquid_polygon_valid;
  stage_b_grid_stats.contact_intersections = contact_intersections;
  stage_b_grid_stats.contact_cells = contact_cells;
  stage_b_grid_stats.extended_cells = extended_cells;
}

static void stage_b_build_auxiliary_heights (scalar f,
                                             double angle_radians)
{
  long algorithm5_marked_contacts = 0;
  long algorithm5_height_valid = 0;
  long algorithm5_height_installs = 0;

  stage_b_ce.height = stage_b_ce_height;
  /* Use Huang's released mark-aware height operator.  Its mark stencil is
   * part of the curvature contract: replacing it with native heights would
   * change which columns are admissible and would no longer be the released
   * Huang operator. */
  stage_b_huang_heights (stage_b_ce, stage_b_mark, stage_b_ce_height);

  foreach() {
    stage_b_height_before.x[] = stage_b_ce_height.x[];
    stage_b_height_before.y[] = stage_b_ce_height.y[];
    stage_b_oxyi.x[] = nodata;
    stage_b_oxyi.y[] = nodata;
    stage_b_angle[] = angle_radians*180./pi;
    stage_b_algorithm5_active[] = 0.;
    stage_b_algorithm5_minus[] = nodata;
    stage_b_algorithm5_center[] = nodata;
    stage_b_algorithm5_plus[] = nodata;
    stage_b_algorithm5_orientation[] = 0.;
  }
  boundary ({stage_b_height_before, stage_b_oxyi, stage_b_angle,
             stage_b_algorithm5_active, stage_b_algorithm5_minus,
             stage_b_algorithm5_center, stage_b_algorithm5_plus,
             stage_b_algorithm5_orientation});

  /* Exact released Huang recompute_h()/embed_h() path.  Its first fraction
   * argument is the complete-cell field used by embed_h(); the native physical
   * f remains separate and is never modified. */
  recompute_h (stage_b_c_complete, stage_b_ce, cs, stage_b_mark, stage_b_angle,
               stage_b_mc, stage_b_ce_height, stage_b_oxyi, 1);
  boundary ((scalar *){stage_b_ce_height, stage_b_oxyi});

  /* The released Huang/Working_case_1 implementation has a named
   * reconstructed-interface Algorithm-5 path for the case where the
   * contact cell has no ordinary height but its adjacent open column does.
   * This is required by a plane normal to the AXI axis: the contact column
   * remains cut by the plane, so the ordinary height is legitimately
   * nodata.  Reuse that source-grounded path verbatim in the isolated
   * auxiliary field; do not manufacture a height or alter physical f. */
  double physical_liquid_measure = 0.;
  foreach (reduction(+:physical_liquid_measure))
    physical_liquid_measure += f[]*dv();

  foreach() {
    if (stage_b_mark[] != 4.)
      continue;

    const bool center_height_valid =
      stage_b_ce_height.x[] != nodata &&
      isfinite (stage_b_ce_height.x[]);
    const bool plus_height_valid =
      stage_b_ce_height.x[0,1] != nodata &&
      isfinite (stage_b_ce_height.x[0,1]);
    const bool neighbor_height_valid = plus_height_valid &&
      (!center_height_valid ||
       orientation(stage_b_ce_height.x[0,1]) ==
       orientation(stage_b_ce_height.x[]));

    coord ns = facet_normal (point, cs, fs);
    coord nf = mycs (point, f);
    coord nc = normal_contact (ns, nf, angle_radians);
    coord p_mof[2] = {{nodata, nodata}, {nodata, nodata}};
    coord polygon[5] = {{nodata, nodata}, {nodata, nodata},
                        {nodata, nodata}, {nodata, nodata},
                        {nodata, nodata}};
    coord solid_segment[2] = {{nodata, nodata}, {nodata, nodata}};
    const double solid_alpha = plane_alpha (cs[], ns);
    polygon_alpha (stage_b_c_complete[], nc, ns, solid_alpha,
                   p_mof, polygon);
    const int nsolid = myfacets (ns, solid_alpha, solid_segment);
    HuangAxiHeightPoint2 intersection = {NAN, NAN};
    const bool intersects = nsolid == 2 &&
      huang_axi_algorithm5_segment_intersection (
        (HuangAxiHeightPoint2){solid_segment[0].x, solid_segment[0].y},
        (HuangAxiHeightPoint2){solid_segment[1].x, solid_segment[1].y},
        (HuangAxiHeightPoint2){p_mof[0].x, p_mof[0].y},
        (HuangAxiHeightPoint2){p_mof[1].x, p_mof[1].y},
        1.e-12, &intersection);

    HuangAxiAlgorithm5Input input = {0};
    input.physical_liquid_measure = physical_liquid_measure;
    HuangAxiAlgorithm5Candidate * candidate = &input.primary;
    candidate->accumulation_valid = intersects && center_height_valid &&
                                    neighbor_height_valid;
    candidate->height_point = (HuangAxiHeightPoint2){
      center_height_valid ? height(stage_b_ce_height.x[]) : NAN, 0.};
    candidate->neighbor_height_point = (HuangAxiHeightPoint2){
      neighbor_height_valid ? height(stage_b_ce_height.x[0,1]) : NAN, 1.};
    if (nsolid == 2) {
      candidate->solid_segment_a =
        (HuangAxiHeightPoint2){solid_segment[0].y, solid_segment[0].x};
      candidate->solid_segment_b =
        (HuangAxiHeightPoint2){solid_segment[1].y, solid_segment[1].x};
    }
    candidate->interface_segment_a =
      (HuangAxiHeightPoint2){p_mof[0].y, p_mof[0].x};
    candidate->interface_segment_b =
      (HuangAxiHeightPoint2){p_mof[1].y, p_mof[1].x};
    candidate->contact_slope = fabs (nc.x) > 1.e-14 ?
      -nc.y/nc.x : NAN;
    candidate->center_q = 0.;
    candidate->ghost_q = -1.;
    candidate->curvature_q = intersects ?
      (y + intersection.q*Delta)/Delta : NAN;

    HuangAxiAlgorithm5Candidate * reconstructed = &input.reconstructed;
    *reconstructed = *candidate;
    reconstructed->accumulation_valid = false;
    if (intersects && !center_height_valid && neighbor_height_valid &&
        huang_axi_algorithm5_point_is_finite(
          reconstructed->interface_segment_a) &&
        huang_axi_algorithm5_point_is_finite(
          reconstructed->interface_segment_b)) {
      reconstructed->height_point = (HuangAxiHeightPoint2){
        0.5*(reconstructed->interface_segment_a.p +
             reconstructed->interface_segment_b.p),
        0.5*(reconstructed->interface_segment_a.q +
             reconstructed->interface_segment_b.q)};
      reconstructed->neighbor_height_point = (HuangAxiHeightPoint2){
        height(stage_b_ce_height.x[0,1]), 1.};
      reconstructed->accumulation_valid = true;
    }

    const HuangAxiAlgorithm5Result correction =
      huang_axi_algorithm5_compute (&input, 1.e-12);
    if (correction.status == HUANG_AXI_ALG5_OK) {
      stage_b_algorithm5_active[] = 1.;
      stage_b_algorithm5_minus[] = correction.height_minus;
      stage_b_algorithm5_center[] = correction.height_center;
      stage_b_algorithm5_plus[] = correction.height_plus;
      stage_b_algorithm5_orientation[] = center_height_valid ?
        orientation(stage_b_ce_height.x[]) :
        orientation(stage_b_ce_height.x[0,1]);
    }
  }
  boundary ({stage_b_algorithm5_active, stage_b_algorithm5_minus,
             stage_b_algorithm5_center, stage_b_algorithm5_plus,
             stage_b_algorithm5_orientation});

  /* Each cell owns its own x-height entry.  The two adjacent writes mirror
   * the already validated Working_case_1 staging rule and stay on this fixed
   * level; no cross-level neighbor write is introduced. */
  foreach() {
    if (stage_b_algorithm5_active[] > 0.5)
      stage_b_ce_height.x[] = huang_axi_algorithm5_encode_height (
        stage_b_algorithm5_center[],
        (int)stage_b_algorithm5_orientation[]);
    else if (stage_b_algorithm5_active[0,1] > 0.5)
      stage_b_ce_height.x[] = huang_axi_algorithm5_encode_height (
        stage_b_algorithm5_minus[0,1],
        (int)stage_b_algorithm5_orientation[0,1]);
    else if (stage_b_algorithm5_active[0,-1] > 0.5)
      stage_b_ce_height.x[] = huang_axi_algorithm5_encode_height (
        stage_b_algorithm5_plus[0,-1],
        (int)stage_b_algorithm5_orientation[0,-1]);
  }
  boundary ({stage_b_ce_height.x});

  foreach (reduction(+:algorithm5_marked_contacts)
           reduction(+:algorithm5_height_valid)
           reduction(+:algorithm5_height_installs)) {
    bool marked = stage_b_mark[] == 4.;
    /* AXI contact curvature uses x=h(r); a valid y-height alone is not an
     * acceptable substitute for this axisymmetric contact path. */
    bool has_height = stage_b_ce_height.x[] != nodata &&
                      isfinite (stage_b_ce_height.x[]);
    bool changed = stage_b_ce_height.x[] != stage_b_height_before.x[] ||
                   stage_b_ce_height.y[] != stage_b_height_before.y[];
    stage_b_contact_path_valid[] = marked && has_height;
    algorithm5_marked_contacts += marked;
    algorithm5_height_valid += marked && has_height;
    algorithm5_height_installs += marked && changed;
  }
  boundary ({stage_b_contact_path_valid});

  foreach() {
    bool force_owner = false;
    /* A face gradient has only two cell endpoints.  The capillary owner is
     * therefore the Algorithm-5 contact cell plus its one-cell face stencil;
     * the two-cell Huang height reach is not a force-support reach. */
    foreach_neighbor (1)
      if (stage_b_contact_path_valid[] > 0.5)
        force_owner = true;
    stage_b_force_support[] = force_owner;
  }
  boundary ({stage_b_force_support});

  stage_b_grid_stats.algorithm5_marked_contacts = algorithm5_marked_contacts;
  stage_b_grid_stats.algorithm5_height_valid = algorithm5_height_valid;
  stage_b_grid_stats.algorithm5_height_installs = algorithm5_height_installs;
}

static void stage_b_rebuild_contact_geometry (scalar f,
                                              double angle_radians)
{
  double ce_min = HUGE, ce_max = -HUGE, open_cell_ce_error = 0.;
  long invalid_cells = 0;
  stage_b_build_auxiliary_fraction (f, angle_radians);
  stage_b_build_auxiliary_heights (f, angle_radians);

  foreach (reduction(min:ce_min) reduction(max:ce_max)
           reduction(max:open_cell_ce_error) reduction(+:invalid_cells)) {
    ce_min = min (ce_min, stage_b_ce[]);
    ce_max = max (ce_max, stage_b_ce[]);
    if (cs[] >= 1. - 1e-12)
      open_cell_ce_error = max (open_cell_ce_error,
                                fabs (stage_b_ce[] - f[]));
    if (!stage_b_fraction_valid (stage_b_ce[]) ||
        !stage_b_fraction_valid (f[]) || !stage_b_fraction_valid (cs[]))
      invalid_cells++;
  }
  stage_b_grid_stats.ce_min = ce_min;
  stage_b_grid_stats.ce_max = ce_max;
  stage_b_grid_stats.open_cell_ce_error = open_cell_ce_error;
  stage_b_grid_stats.invalid_cells = invalid_cells;
}

#endif
