#ifndef STAGE_B_HUANG_EMBED_CURVATURE_H
#define STAGE_B_HUANG_EMBED_CURVATURE_H

/*
 * Namespace wrapper for the released Chongsen EBM_VOF curvature operator.
 *
 * The native Basilisk curvature.h is already required by the physical f
 * path.  Chongsen's operator has the additional (tmp_c, mark) contract, so
 * it is included under private names rather than replacing or modifying the
 * native operator.  The two source files included below are byte-copied from
 * Basiliks_sanbox/sandbox/Chongsen/src/EBM_VOF and are not a new curvature
 * algorithm.
 */

#define cstats                    stage_b_huang_cstats
#define curvature                 stage_b_huang_curvature
#define curvature_prolongation    stage_b_huang_curvature_prolongation
#define curvature_restriction     stage_b_huang_curvature_restriction
#define independents               stage_b_huang_independents
#define height_curvature            stage_b_huang_height_curvature
#define height_curvature_fit        stage_b_huang_height_curvature_fit
#define centroids_curvature_fit    stage_b_huang_centroids_curvature_fit
#define height_normal               stage_b_huang_height_normal
#define height_normal_z             stage_b_huang_height_normal_z
#define height_position             stage_b_huang_height_position
#define pos_x                       stage_b_huang_pos_x

#define kappa_x                     stage_b_huang_kappa_x
#define kappa_y                     stage_b_huang_kappa_y
#define kappa_z                     stage_b_huang_kappa_z
#define kappa1_x                    stage_b_huang_kappa1_x
#define kappa1_y                    stage_b_huang_kappa1_y
#define kappa1_z                    stage_b_huang_kappa1_z
#define normal_x                    stage_b_huang_normal_x
#define normal_y                    stage_b_huang_normal_y
#define normal_z                    stage_b_huang_normal_z
#define normal2_x                   stage_b_huang_normal2_x
#define normal2_y                   stage_b_huang_normal2_y
#define normal2_z                   stage_b_huang_normal2_z

#define height                      stage_b_huang_height
#define orientation                 stage_b_huang_orientation
#define heights                     stage_b_huang_heights
#define half_column                 stage_b_huang_half_column
#define column_propagation          stage_b_huang_column_propagation
#define refine_h_x                  stage_b_huang_refine_h_x
#define STAGE_B_HUANG_SKIP_INTERFACE_NORMAL 1

#include "huang_original/embed_curvature.h"

#undef refine_h_x
#undef column_propagation
#undef half_column
#undef heights
#undef STAGE_B_HUANG_SKIP_INTERFACE_NORMAL
#undef orientation
#undef height
#undef normal2_z
#undef normal2_y
#undef normal2_x
#undef normal_z
#undef normal_y
#undef normal_x
#undef kappa1_z
#undef kappa1_y
#undef kappa1_x
#undef kappa_z
#undef kappa_y
#undef kappa_x
#undef pos_x
#undef height_position
#undef height_normal_z
#undef height_normal
#undef centroids_curvature_fit
#undef height_curvature_fit
#undef height_curvature
#undef independents
#undef curvature_restriction
#undef curvature_prolongation
#undef curvature
#undef cstats

#endif
