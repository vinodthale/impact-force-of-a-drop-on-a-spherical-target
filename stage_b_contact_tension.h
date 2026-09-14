#ifndef STAGE_B_CONTACT_TENSION_H
#define STAGE_B_CONTACT_TENSION_H

#if dimension != 2
# error "Stage B contact tension is strictly two-dimensional"
#endif

#ifndef STAGE_B_CONTACT_ANGLE_RADIANS
# error "Define STAGE_B_CONTACT_ANGLE_RADIANS before including this header"
#endif

/*
 * Include this header before tension.h in the isolated B2 driver.  Basilisk's
 * native tension/iforce event is then registered after this header and, by
 * Basilisk's event-overload ordering, executes before this correction event.
 * The native event therefore retains ownership of physical-f clamping. Stage
 * B replaces the native capillary contribution only on the deterministic
 * Algorithm-5 contact support.  An auxiliary-field jump outside that support
 * remains native Basilisk-owned.
 */
#include "stage_b_contact_grid.h"

#ifndef STAGE_B_SURFACE_TENSION
# error "Define STAGE_B_SURFACE_TENSION before including this header"
#endif

scalar stage_b_native_kappa[], stage_b_contact_kappa[];
scalar stage_b_combined_kappa[];

typedef struct {
  long noncontact_faces;
  long contact_faces;
  long one_sided_contact_faces;
  long invalid_contact_faces;
  long nonfinite_acceleration_faces;
  double physical_f_change;
  double acceleration_max;
} StageBForceStats;

static StageBForceStats stage_b_force_stats;

/* Native iforce.h also owns the acceleration allocation and native tension.h
 * owns the capillary timestep.  These fallbacks are compiled only for a
 * direct geometry smoke driver that omits native tension.h. */
#ifndef STAGE_B_NATIVE_TENSION_AVAILABLE
event defaults (i = 0)
{
  if (is_constant (a.x)) {
    a = new face vector;
    foreach_face() {
      a.x[] = 0.;
      dimensional (a.x[] == Delta/sq(DT));
    }
  }
}

event stability (i++)
{
  /* The native tension.h stability event is authoritative in B2.  This
   * fallback keeps the isolated geometry-only header self-contained when no
   * native tension header is included. */
  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax)
                reduction(min:dmin))
    if (fm.x[] > 0.) {
      amin = min (amin, alpha.x[]/fm.x[]);
      amax = max (amax, alpha.x[]/fm.x[]);
      dmin = min (dmin, Delta);
    }
  if (STAGE_B_SURFACE_TENSION > 0. && amin < HUGE && amax > -HUGE &&
      dmin < HUGE) {
    double rhom = (1./amin + 1./amax)/2.;
    double capillary_dt = sqrt (rhom*cube(dmin)/
                                (pi*STAGE_B_SURFACE_TENSION));
    if (capillary_dt < dtmax)
      dtmax = capillary_dt;
  }
}
#endif

static void stage_b_contact_geometry_post_vof (void)
{
  double before = 0., after = 0.;
  foreach (reduction(+:before))
    before += f[]*dv();
  stage_b_rebuild_contact_geometry (f, STAGE_B_CONTACT_ANGLE_RADIANS);
  foreach (reduction(+:after))
    after += f[]*dv();
  /* Stage-B is forbidden from changing physical f.  The native iforce.h
   * clamp, when present, is outside this measured adapter action. */
  stage_b_force_stats.physical_f_change = after - before;
}

static inline double stage_b_face_potential (double left, double right)
{
  bool left_valid = left < nodata && isfinite (left);
  bool right_valid = right < nodata && isfinite (right);
  return left_valid && right_valid ? (left + right)/2. :
         left_valid ? left : right_valid ? right : nodata;
}

/*
 * The event is registered before tension.h and is placed after native
 * tension/iforce by Basilisk's overload chain. It therefore performs the
 * source-equivalent replacement on every contact-supported face with an
 * auxiliary-field jump:
 *
 *   a <- a_native - sigma*kappa_f*grad(f)/rho
 *                    + sigma*kappa_ce*grad(ce)/rho.
 *
 * The two terms use the exact native alpha/fm/Delta placement.  No second
 * radial or azimuthal metric factor is introduced.
 */
event acceleration (i++)
{
  stage_b_contact_geometry_post_vof ();
  curvature (f, stage_b_native_kappa);
  stage_b_huang_curvature (f, stage_b_ce, stage_b_mark,
                           stage_b_contact_kappa, 1., false);

  foreach()
    stage_b_combined_kappa[] = stage_b_contact_kappa[] < nodata &&
      isfinite (stage_b_contact_kappa[]) ?
      stage_b_contact_kappa[] : stage_b_native_kappa[];
  boundary ({stage_b_combined_kappa});

  long noncontact_faces = 0, contact_faces = 0;
  long one_sided_contact_faces = 0;
  long invalid_contact_faces = 0, nonfinite_acceleration_faces = 0;
  double acceleration_max = 0.;
  face vector ia = a;

  foreach_face (reduction(+:noncontact_faces) reduction(+:contact_faces)
                reduction(+:one_sided_contact_faces)
                reduction(+:invalid_contact_faces)
                reduction(+:nonfinite_acceleration_faces)
                reduction(max:acceleration_max))
    if (fm.x[] > 0.) {
      const double native_jump = f[] - f[-1];
      const bool contact_owned = stage_b_force_support[] > 0.5 ||
                                 stage_b_force_support[-1] > 0.5;
      if (!contact_owned) {
        if (native_jump != 0.)
          noncontact_faces++;
        continue;
      }

      const double contact_jump = stage_b_ce[] - stage_b_ce[-1];
      if (contact_jump == 0.) {
        if (native_jump != 0.)
          noncontact_faces++;
        continue;
      }

      const double native_potential = stage_b_face_potential (
        stage_b_native_kappa[], stage_b_native_kappa[-1]);
      /* Preserve the released Huang embed_iforce.h contract here.  A
       * missing contact curvature is nodata, not native curvature.  If only
       * one endpoint has a valid Huang curvature, the released operator uses
       * that endpoint alone.  Substituting native curvature before this
       * average would halve the correction on one-sided contact faces. */
      const bool left_contact_kappa_valid =
        stage_b_contact_kappa[] < nodata &&
        isfinite (stage_b_contact_kappa[]);
      const bool right_contact_kappa_valid =
        stage_b_contact_kappa[-1] < nodata &&
        isfinite (stage_b_contact_kappa[-1]);
      const long one_sided_contact =
        left_contact_kappa_valid != right_contact_kappa_valid;
      const double contact_potential = stage_b_face_potential (
        stage_b_contact_kappa[], stage_b_contact_kappa[-1]);

#ifdef STAGE_B_DEBUG_FORCE
      fprintf (stderr,
               "STAGE_B_FORCE_FACE x=%.17g y=%.17g f=(%.17g,%.17g) "
               "ce=(%.17g,%.17g) kf=(%.17g,%.17g) kc=(%.17g,%.17g) "
               "jf=%.17g jc=%.17g\n", x - Delta/2., y, f[], f[-1],
               stage_b_ce[], stage_b_ce[-1], stage_b_native_kappa[],
               stage_b_native_kappa[-1], stage_b_contact_kappa[],
               stage_b_contact_kappa[-1], native_jump, contact_jump);
#endif

      if (contact_jump != 0. &&
          (contact_potential >= nodata || !isfinite (contact_potential))) {
        /* Leave the native force untouched when the contact reconstruction is
         * invalid; this is logged and is never converted into a guessed force. */
        invalid_contact_faces++;
        continue;
      }

      double native_increment = 0., contact_increment = 0.;
      if (native_jump != 0. && native_potential < nodata &&
          isfinite (native_potential))
        native_increment = alpha.x[]/(fm.x[] + SEPS)*
                           STAGE_B_SURFACE_TENSION*
                           native_potential*native_jump/Delta;
      if (contact_jump != 0.)
        contact_increment = alpha.x[]/(fm.x[] + SEPS)*
                            STAGE_B_SURFACE_TENSION*
                            contact_potential*contact_jump/Delta;
      const double correction = contact_increment - native_increment;
      if (!isfinite (correction)) {
        nonfinite_acceleration_faces++;
        continue;
      }
      ia.x[] += correction;
      acceleration_max = max (acceleration_max, fabs (correction));
      contact_faces++;
      one_sided_contact_faces += one_sided_contact;
    }

  stage_b_force_stats.noncontact_faces = noncontact_faces;
  stage_b_force_stats.contact_faces = contact_faces;
  stage_b_force_stats.one_sided_contact_faces = one_sided_contact_faces;
  stage_b_force_stats.invalid_contact_faces = invalid_contact_faces;
  stage_b_force_stats.nonfinite_acceleration_faces =
    nonfinite_acceleration_faces;
  stage_b_force_stats.acceleration_max = acceleration_max;
}

#endif
