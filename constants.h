/*
 * UNDERGRADUATE USER CONTRACT -- READ BEFORE EDITING
 * ==================================================
 * This file is the only supported parameter-control file for this case.
 * Students may change values only in sections explicitly marked USER INPUT.
 * Do not rename macros, add forces, alter equations, or edit the files below:
 *
 *   stage_b_contact_grid.h
 *   stage_b_contact_tension.h
 *   stage_b_geometry.c
 *   stage_b_geometry.h
 *   stage_b_huang_embed_curvature.h
 *   huang_axi_locked_include/ (all files)
 *   huang_original/ (all files)
 *
 * Those files implement the reviewed Huang/Basilisk contact-line machinery;
 * their numerical tolerances and reconstruction logic are not case inputs.
 * A parameter absent from the USER INPUT sections is not student-tunable.
 *
 * MODEL CONTRACT
 * --------------
 * - Dimension: 2D two-dimensional axisymmetric (AXI), No  three-dimensional solver .
 * - Coordinates: x is axial, y is radial, and y = 0 is the symmetry axis.
 * - Solver reference scales: D = 1, U0 = 1, and rho_l = 1.
 * - Runtime variables are nondimensional. SI values, when selected below,
 *   are used only to derive Re, We, Fr, density ratio and viscosity ratio.
 * - Basilisk owns physical-f VOF transport, AMR and AXI+EMBED metrics.
 * - Stage B reconstructs Huang contact geometry/curvature and replaces the
 *   native capillary contribution only on its certified contact support.
 * - No clipping, smoothing, damping, mass repair or ad-hoc force is permitted.
 */

/* USER INPUT 1: static contact angle measured through the liquid [degrees].
 * This target-curvature study keeps liquid wetting at 90 degrees.  The
 * zero-degree case label refers to solid curvature, not liquid wetting. */
#define CONTACT_ANGLE_DEGREES 90.

/* DERIVED -- radians required by the trigonometric contact reconstruction.
 * Do not edit this expression; change CONTACT_ANGLE_DEGREES above. */
#define CONTACT_ANGLE_RADIANS (CONTACT_ANGLE_DEGREES*pi/180.)

/* USER INPUT 2: nominal axial position of the flat embedded plane divided by
 * drop diameter D.  The solid fills x < EMBED_PLANE_X*D and the fluid
 * occupies the region to its right; Bdropimpact.c adds a documented subcell
 * alignment displacement. */
#define EMBED_PLANE_X                  0.25

/* CASE IDENTITY -- synchronized with the two enclosing folder names.
 * These strings label outputs and diagnostics only; they do not alter the
 * governing equations, geometry, material properties, or numerical methods. */
#define CASE_VELOCITY_FOLDER            "Vel5.00"
#define CASE_GEOMETRY_FOLDER            "0_degree"
#define CASE_ID                         "Vel5.00_0_degree"
#define SNAPSHOT_DIRECTORY              "intermediate"
#define SNAPSHOT_PREFIX                 SNAPSHOT_DIRECTORY "/snapshot-"
#define FORCE_POST_DIRECTORY            "Post_" CASE_ID
#define FORCE_HISTORY_BASENAME          "impact_force_history_" CASE_ID ".txt"
#define BVIEW_DIRECTORY                 "bviewfiles_" CASE_ID
#define BVIEW_MOVIE_BASENAME            "Bview_" CASE_ID ".mp4"

/* LOCKED SOLVER STACK -- include order is part of the AXI+EMBED formulation.
 * embed.h must precede axi.h so cm/fm combine cut-cell and radial metrics. */
#include "grid/quadtree.h"            // 2D adaptive quadtree and MPI backend
#include "embed.h"                    // cs/fs embedded-solid geometry
#include "axi.h"                      // AXI metrics; x=axial and y=radial
#include "navier-stokes/centered.h"   // incompressible centered momentum solver
#define FILTERED  1                   // native filtered density/viscosity jump
#include "two-phase.h"                // physical VOF f and two-fluid properties

/* LOCKED STAGE-B BRIDGE -- these are aliases, not additional user parameters. */
#define STAGE_B_CONTACT_ANGLE_RADIANS CONTACT_ANGLE_RADIANS // one angle owner
/* Stage B reads Basilisk's physical-interface surface tension. f.sigma is
 * assigned from cfdbv.Sigma in main(); no duplicate sigma is introduced. */
#define STAGE_B_SURFACE_TENSION (f.sigma) // one capillary-coefficient owner
#define STAGE_B_NATIVE_TENSION_AVAILABLE 1 // native capillary dt/force enabled
#include "stage_b_contact_tension.h"  // reviewed Huang contact reconstruction
#include "tension.h"                  // native Basilisk CSF surface tension
#include "reduced.h"                  // reduced-gravity acceleration field G
#include "tag.h"                      // connected-component identification

/* USER INPUT 3: parameterization mode.
 * 'd' = derive nondimensional groups from the SI reference inputs below.
 * 'n' = prescribe nondimensional groups and property ratios directly.
 * 'e' = experimental reconstruction mode using prescribed Re and We.
 * The solver remains nondimensional in every mode. Undergraduate production
 * runs must retain 'd'; modes 'n' and 'e' require supervisor code review. */
#define DIM_NONDIM_EXP			'd'

#if DIM_NONDIM_EXP == 'd' || DIM_NONDIM_EXP == 'D'

/* USER INPUT 3A: dimensional reference data used only for conversion.
 * Re = rho_l U0 D/mu_l; We = rho_l U0^2 D/sigma;
 * Fr = U0/sqrt(gD). The solver subsequently advances D=U0=rho_l=1. */
#define VELOCITY			5.00       // U0 [m/s], Drop impact velocity reference
#define DROP_DIAMETER		2.050e-03  // D [m], initial liquid-drop diameter
#define RHO_L				998.0      // rho_l [kg/m^3], liquid density
#define RHO_G				1.21       // rho_g [kg/m^3], gas density
#define MU_L				0.001      // mu_l [Pa s], liquid dynamic viscosity
#define MU_G				1.81e-5    // mu_g [Pa s], gas dynamic viscosity
#define SIGMA				0.073      // sigma [N/m], liquid-gas surface tension
#define GRAVITY				9.81       // g [m/s^2], gravitational magnitude

/* INACTIVE PLACEHOLDERS IN 'd' MODE -- do not tune these values. The code
 * derives the corresponding groups and ratios from the SI inputs above. */
#define RHO_GL				0.0        // unused rho_g/rho_l placeholder
#define MU_GL				0.0        // unused mu_g/mu_l placeholder
#define REYNOLDS			0.0        // derived from SI inputs in this mode
#define WEBER				0.0        // derived from SI inputs in this mode
#define FROUDE              0.0        // derived from SI inputs in this mode

#elif DIM_NONDIM_EXP == 'n' || DIM_NONDIM_EXP == 'N'

/* SUPERVISOR-ONLY MODE: prescribe nondimensional groups directly.
 * These values determine mu_l*=1/Re, sigma*=1/We and |G*|=1/Fr^2. */
#define WEBER				100.0      // We = rho_l U0^2 D/sigma
#define REYNOLDS			100.0      // Re = rho_l U0 D/mu_l
#define FROUDE              70.0       // Fr = U0/sqrt(gD)
#define RHO_GL				(0.0012)   // density ratio rho_g/rho_l
#define MU_GL				(0.0210)   // viscosity ratio mu_g/mu_l

/* INACTIVE PLACEHOLDERS IN 'n' MODE -- physical SI inputs are not used. */
#define VELOCITY			0.0        // unused SI-speed placeholder
#define DROP_DIAMETER		0.0        // unused SI-diameter placeholder
#define RHO_L				0.0        // unused SI-density placeholder
#define MU_L				0.0        // unused SI-viscosity placeholder
#define SIGMA				0.0        // unused SI-surface-tension placeholder
#define RHO_G				0.0        // unused SI-density placeholder
#define MU_G				0.0        // unused SI-viscosity placeholder
#define GRAVITY				0.0        // unused SI-gravity placeholder

#elif DIM_NONDIM_EXP == 'e' || DIM_NONDIM_EXP == 'E'

/* SUPERVISOR-ONLY LEGACY MODE: reconstruct U0 and mu_l from prescribed
 * Re/We plus selected SI properties. Do not use for student production runs
 * without reviewing numericalmainvalues() and all derived diagnostics. */
#define WEBER				300.0      // prescribed Weber number
#define REYNOLDS			1000.0     // prescribed Reynolds number
#define FROUDE              70.0       // legacy Fr entry; verify before use
#define DROP_DIAMETER		2.0e-3     // D [m], reconstruction reference
#define SIGMA				17.6e-3    // sigma [N/m], reconstruction reference
#define RHO_L				816.0      // rho_l [kg/m^3], reconstruction reference
#define RHO_G				1.2041     // rho_g [kg/m^3], gas property
#define MU_G				1.94e-5    // mu_g [Pa s], gas property

/* INACTIVE/DERIVED PLACEHOLDERS IN 'e' MODE -- do not tune directly. */
#define VELOCITY			0.0        // U0 is reconstructed from We
#define RHO_GL				0.0        // density ratio is derived from SI inputs
#define MU_GL				0.0        // viscosity ratio is derived from SI inputs
#define MU_L				0.0        // mu_l is reconstructed from Re
#define GRAVITY				0.0        // legacy placeholder; review Fr handling

#endif

/* USER INPUT 4: grid, domain, geometry spacing and output cadence.
 * All lengths are normalized by D; all times are normalized by D/U0. */
#define INITAL_GRID_LEVEL		9     // base level; initial grid has 2^9 cells/axis
#define MAX_GRID_LEVEL			12   // finest AMR level; Delta_min=L0/2^12
#define DOMAIN_WIDTH			4.00  // L0/D, square meridional domain width
#define INITIAL_DISTANCE		0.04  // initial drop-surface/solid-surface gap g/D
#define BUBBLE_DIAMETER		    0.00  // energy-diagnostic normalization only; no bubble created
#define DBDELTA       		    0.00  // legacy drop-bubble spacing; inactive here
#define REFINE_GAP				0.02  // retained geometry-band scale; active AMR is controlled by REFINE_VAR and its tolerances
#define MAX_TIME				2.00  // simulated time after nominal contact, in D/U0
#define SAVE_FILE_EVERY			0.01  // snapshot/diagnostic interval, in D/U0

/* USER INPUT 5: AMR error tolerances.
 * adapt_wavelet() monitors variables in this exact order. Each integer q
 * means an absolute wavelet tolerance 10^q for the corresponding field.
 * REFINE_VAR and REFINE_VAR_TEXT must remain synchronized; edit tolerances
 * only after a documented grid-sensitivity study. */
#define REFINE_VAR				{cs, f, u.x, u.y} // monitored EMBED, VOF and velocity fields
#define REFINE_VAR_TEXT			"cs, f, u.x, u.y" // same order, for metadata only
#define REFINE_VALUE_0			-6 // cs tolerance: 10^-6
#define REFINE_VALUE_1			-6 // f tolerance: 10^-6
#define REFINE_VALUE_2			-3 // u.x tolerance: 10^-3
#define REFINE_VALUE_3			-3 // u.y tolerance: 10^-3

/* LEGACY REPORTING CONTROLS -- the current Bdropimpact.c has no component-
 * deletion event. Keep both switches at 'n'. These entries are printed for
 * compatibility and must not be presented as active numerical operations. */
#define REMOVE_DROP_YESNO		'n' // inactive; no drop deletion in active solver
#define REMOVE_DROP_SIZE		4.0 // inactive legacy cell-size threshold
#define REMOVE_DROP_PERIOD		4   // inactive legacy iteration interval
#define REMOVE_BUBBLE_YESNO		'n' // inactive; no bubble deletion in active solver
#define REMOVE_BUBBLE_SIZE		4.0 // inactive legacy cell-size threshold
#define REMOVE_BUBBLE_PERIOD	4   // inactive legacy iteration interval

/* LOCKED OUTPUT BASENAMES -- change only if a downstream reader is updated. */
#define FILENAME_DATA			"data_" CASE_ID           // unused legacy basename; retained for compatibility
#define FILENAME_DURATION		"duration_" CASE_ID       // per-rank timing basename
#define FILENAME_PARAMETERS		"parameters_" CASE_ID ".txt" // resolved case-parameter record
#define FILENAME_ENDOFRUN		"endofrun_" CASE_ID       // per-rank completion record
#define FILENAME_LASTFILE		"lastfile"       // latest restart dump

/* LEGACY COMPATIBILITY CONSTANTS -- not active student controls. */
#define R_VOFLIMIT				1.0e-9 // retained VOF tolerance symbol; unused here
#define R_PI					3.1415926535897932384626433832795 // legacy pi; use Basilisk pi

/* INTERNAL RUNTIME STATE -- do not edit below this boundary. */
int LEVELmin = INITAL_GRID_LEVEL, LEVELmax = MAX_GRID_LEVEL; // CLI-overridable levels
double maxruntime = HUGE;       // unused legacy wall-clock-limit placeholder
scalar fdrop[], pressure[];     // visualization field and dumped pressure copy
scalar fb[];                    // legacy bubble identity field; not in interfaces
vector huang_contact_height[];  // legacy field declaration; inactive in Stage B

struct CFDValues {
	/* rhoL/rhoG: nondimensional densities; muL/muG: nondimensional dynamic
	 * viscosities; Sigma: sigma*=1/We. */
	double rhoL, rhoG, muL, muG, Sigma;
	/* vel: U0*=1; Reynolds/Weber: governing groups; Froude: gravity group;
	 * Bond/Oh/GXnormlised: 'd'-mode diagnostics, with GXnormlised=1/Fr^2. */
	double vel, Reynolds, Weber, Froude,Bond, Oh,GXnormlised;
	/* diameter: D*=1; domainsize/refinegap/embed_plane_x/initialdis: solver-
	 * coordinate lengths normalized by D. */
	double diameter, domainsize, refinegap, embed_plane_x, initialdis;
	/* timecontact/timeend: nondimensional times; timestep is a legacy stored
	 * cadence while the active output event uses SAVE_FILE_EVERY. */
	double timecontact, timeend, timestep;
	double bubblediameter;    // legacy energy-diagnostic normalization
	double dbdelta;           // legacy drop-bubble separation
};

/* INTERNAL CLI PARSER -- coursework should normally change USER INPUT values
 * above and launch without overrides. The implementation is not tunable. */
void readfromarg(char **argv, int argc, struct CFDValues *bvalues);

/* INTERNAL PARAMETER PROJECTION
 * Converts the selected input mode into the nondimensional solver state.
 * The invariant reference values D*=1, U0*=1 and rho_l*=1 are assigned first;
 * physical properties enter the solver only through dimensionless groups.
 * Do not edit equations in this function as part of a parameter study. */
int numericalmainvalues(char **argv, int argc, struct CFDValues *bvalues)
{
	double velocity = VELOCITY, mu_l = MU_L; // conversion-layer references
	bvalues->rhoL = 1.0;                    // enforce rho_l*=1
	bvalues->vel = 1.0;                     // enforce U0*=1
	bvalues->diameter = 1.0;                // enforce D*=1
	;
	bvalues->Reynolds = -1.0; // negative sentinel: use selected-mode default
	bvalues->Weber = -1.0;    // negative sentinel: use selected-mode default
	bvalues->timeend = -1.0;  // negative sentinel: use MAX_TIME
	bvalues->timestep = -1.0; // negative sentinel: use SAVE_FILE_EVERY
	;
	readfromarg(argv, argc, bvalues); // apply optional expert CLI overrides
	switch(DIM_NONDIM_EXP)
	{
	case 'd':
	case 'D':
	{
		bvalues->Reynolds = (RHO_L * VELOCITY * DROP_DIAMETER / mu_l); // Re
		bvalues->Weber = (RHO_L * VELOCITY * VELOCITY * DROP_DIAMETER / SIGMA); // We
		bvalues->Froude = (VELOCITY / sqrt (GRAVITY * DROP_DIAMETER)); // Fr
		;
		bvalues->Bond = (RHO_L * GRAVITY * DROP_DIAMETER  * DROP_DIAMETER / SIGMA); // Bo
		bvalues->Oh = (mu_l/ sqrt(RHO_L * SIGMA * DROP_DIAMETER)); // Oh
		;
		bvalues->muL = (bvalues->rhoL * bvalues->vel * bvalues->diameter / bvalues->Reynolds); // mu_l*=1/Re
		bvalues->Sigma = (bvalues->rhoL * bvalues->vel * bvalues->vel * bvalues->diameter / bvalues->Weber); // sigma*=1/We
		bvalues->rhoG = (RHO_G / RHO_L) * bvalues->rhoL; // rho_g*=rho_g/rho_l
		bvalues->muG = (MU_G / mu_l)  * bvalues->muL; // mu_g*=(mu_g/mu_l)mu_l*
		bvalues->GXnormlised = (GRAVITY * DROP_DIAMETER / (VELOCITY*VELOCITY)); // g*=1/Fr^2
		;
	    bvalues->bubblediameter = BUBBLE_DIAMETER * bvalues->diameter; // legacy Db*
	    bvalues->dbdelta = DBDELTA * bvalues->diameter; // legacy separation*
		break;
	}
	case 'n':
	case 'N':
	{
		if (bvalues->Reynolds < 0.0)
			bvalues->Reynolds = REYNOLDS; // use configured Re unless CLI supplied it
		if (bvalues->Weber < 0.0)
			bvalues->Weber = WEBER; // use configured We unless CLI supplied it
		if (bvalues->Froude < 0.0)
			bvalues->Froude = FROUDE; // use configured Fr unless CLI supplied it
		bvalues->muL = (bvalues->rhoL * bvalues->vel * bvalues->diameter / bvalues->Reynolds); // mu_l*=1/Re
		bvalues->Sigma = (bvalues->rhoL * bvalues->vel * bvalues->vel * bvalues->diameter / bvalues->Weber); // sigma*=1/We
		bvalues->rhoG = RHO_GL * bvalues->rhoL; // apply configured density ratio
		bvalues->muG = MU_GL * bvalues->muL; // apply configured viscosity ratio
		break;
	}
	case 'e':
	case 'E':
	{
		if (bvalues->Reynolds < 0.0)
			bvalues->Reynolds = REYNOLDS; // prescribed reconstruction Re
		if (bvalues->Weber < 0.0)
			bvalues->Weber = WEBER; // prescribed reconstruction We
		bvalues->muL = (bvalues->rhoL * bvalues->vel * bvalues->diameter / bvalues->Reynolds); // solver mu_l*
		bvalues->Sigma = (bvalues->rhoL * bvalues->vel * bvalues->vel * bvalues->diameter / bvalues->Weber); // solver sigma*
		velocity = sqrt (bvalues->Weber * SIGMA / (DROP_DIAMETER * RHO_L)); // reconstructed U0 [m/s]
		mu_l = (RHO_L * velocity * DROP_DIAMETER / bvalues->Reynolds); // reconstructed mu_l [Pa s]
		bvalues->muG = (MU_G / mu_l) * bvalues->muL; // solver gas viscosity
		bvalues->rhoG = (RHO_G / RHO_L) * bvalues->rhoL; // solver gas density
		break;
	}
	}
	bvalues->domainsize = DOMAIN_WIDTH * bvalues->diameter; // L0*=DOMAIN_WIDTH D*
	bvalues->embed_plane_x = EMBED_PLANE_X * bvalues->diameter; // nominal plane x*=EMBED_PLANE_X D*
	if (!(bvalues->embed_plane_x > 0.)) {
		fprintf(ferr, "EMBED_PLANE_X must be positive for the planar surface position.\n");
		exit(1);
	}
	bvalues->initialdis = INITIAL_DISTANCE * bvalues->diameter; // physical surface gap g*
	bvalues->refinegap = REFINE_GAP * bvalues->diameter; // retained geometry-band scale for provenance
	;
	bvalues->timecontact = bvalues->initialdis / bvalues->vel; // nominal g*/U0*
	if (bvalues->timeend < 0.0)
		bvalues->timeend = MAX_TIME; // default post-contact duration
	if (bvalues->timestep < 0.0)
		bvalues->timestep = SAVE_FILE_EVERY; // default output cadence
	;
	/* Rank zero writes the resolved inputs and derived nondimensional values.
	 * This report is provenance only; it does not alter solver fields. */
	switch (pid())
	{
	case 0:
	{
		printf("CASE: %s --- R: %f --- W: %f --- EMBED: plane_normal_to_axis\r\n",
		       CASE_ID, bvalues->Reynolds, bvalues->Weber);
		FILE *fp;
		fp = fopen (FILENAME_PARAMETERS, "w");
		fprintf (fp, "Case ID: %s\r\n", CASE_ID);
		fprintf (fp, "Velocity Folder: %s\r\n", CASE_VELOCITY_FOLDER);
		fprintf (fp, "Geometry Folder: %s\r\n", CASE_GEOMETRY_FOLDER);
	    fprintf (fp, "Name of Liquid : Water drop  \r\n");
	    fprintf (fp, "Experimental parameters  / Numerical simulation parameters in Basilisk unit\r\n");
		fprintf (fp, "Diameter_Experimental: %.3e / Normalized diameter  of drop (D): %.3e\r\n", DROP_DIAMETER, bvalues->diameter);
		fprintf (fp, "Velocity_Experimental: %.3e / Normalized Velocity  of liquid (V) : %.3e\r\n", velocity, bvalues->vel);
		fprintf (fp, "Rho(L)_Experimental: %.3e /  Normalized density of liquid (rho1): %.3e\r\n", RHO_L, bvalues->rhoL); 
		fprintf (fp, "Rho(G)_Experimental: %.3e / Normalized density  of Air (rho2) : %.3e\r\n", RHO_G, bvalues->rhoG);
		fprintf (fp, "Mu(L)_Experimental: %.3e / Normalized viscosity of liquid (mu1): %.3e\r\n", mu_l, bvalues->muL);
		fprintf (fp, "Mu(G)_Experimental: %.3e / Normalized viscosity of Air (mu2): %.3e\r\n", MU_G, bvalues->muG);
		fprintf (fp, "Sigma(L-G)_Experimental: %.3e / Normalized surface tension (f.sigma): %.3e\r\n", SIGMA, bvalues->Sigma);
		fprintf (fp, "Acceleration due to gravity: %.3e / Normalized gravity (G.x): %.3e\r\n", GRAVITY, bvalues->GXnormlised);
		fprintf (fp, "\r\n");
		fprintf (fp, "Reynolds: %.10f\r\n", bvalues->Reynolds);
		fprintf (fp, "Weber: %.10f\r\n", bvalues->Weber);
		fprintf (fp, "Froude: %.10f\r\n", bvalues->Froude);
		fprintf (fp, "Bond: %.10f\r\n", bvalues->Bond);
		fprintf (fp, "Ohsorge: %.10f\r\n", bvalues->Oh);
		fprintf (fp, "\r\n");
		fprintf (fp, "Level Max: %d\r\n", LEVELmax);
		fprintf (fp, "Level Min: %d\r\n", LEVELmin);
		fprintf (fp, "Domain Size: %.2f\r\n", bvalues->domainsize);
		const double embed_plane_eps =
			bvalues->domainsize/(1 << LEVELmax)/1000.;
		fprintf (fp, "EMBED Geometry: plane_normal_to_axis\r\n");
		fprintf (fp, "EMBED Nominal Plane x: %.17e\r\n", bvalues->embed_plane_x);
		fprintf (fp, "EMBED Alignment Epsilon: %.17e\r\n",
		         embed_plane_eps);
		fprintf (fp, "EMBED Actual Plane x: %.17e\r\n",
		         bvalues->embed_plane_x + embed_plane_eps);
		fprintf (fp, "Initial Surface Gap: %.17e\r\n", bvalues->initialdis);
		fprintf (fp, "Configured Refine Gap Scale: %.2f\r\n",
		         bvalues->refinegap);
		fprintf (fp, "Contact Time: %.2f\r\n", bvalues->timecontact);
		fprintf (fp, "Domain Size: %.2f\r\n", bvalues->domainsize);
		fprintf (fp, "\r\n");
		fprintf (fp, "Bubble Diameter: %.2f\r\n", bvalues->bubblediameter);     
		fprintf (fp, "Dbdelta (disatance btw drop and bubble) : %.6f\r\n", bvalues->dbdelta); 
		fprintf (fp, "\r\n");  
		fprintf (fp, "Refine Variables: %s\r\n", REFINE_VAR_TEXT);
		fprintf (fp, "Refine Variables powers: %d, %d, %d, %d\r\n",
		         REFINE_VALUE_0, REFINE_VALUE_1, REFINE_VALUE_2,
		         REFINE_VALUE_3);
		fprintf (fp, "\r\n");
		fprintf (fp, "Remove Drop YesNo: %c\r\n", REMOVE_DROP_YESNO);
		fprintf (fp, "Remove Drop Size: %f\r\n", REMOVE_DROP_SIZE);
		fprintf (fp, "Remove Drop Period: %d\r\n", REMOVE_DROP_PERIOD);
		fprintf (fp, "\r\n");
		fprintf (fp, "Remove Bubble YesNo: %c\r\n", REMOVE_BUBBLE_YESNO);
		fprintf (fp, "Remove Bubble Size: %f\r\n", REMOVE_BUBBLE_SIZE);
		fprintf (fp, "Remove Bubble Period: %d\r\n", REMOVE_BUBBLE_PERIOD);
		fclose (fp);
		break;
	}
	}
	return 1;
}

/* INTERNAL OPTIONAL COMMAND-LINE OVERRIDES
 * R/W/F are ignored or only partly honored outside their compatible legacy
 * modes; in the supported 'd' mode Re, We and Fr are recomputed from SI data.
 * R<value>  Reynolds number; W<value> Weber number; F<value> Froude number.
 * X<level>  maximum AMR level; N<level> minimum/base AMR level.
 * TE<time>  active post-contact end time override.
 * TS<time>  legacy stored value; active output event uses SAVE_FILE_EVERY.
 * Students should use the documented USER INPUT macros instead. Do not add
 * new parser keys without supervisor review and matching provenance output. */
void readfromarg(char **argv, int argc, struct CFDValues *bvalues)
{
	int i, j;
	char tmp[100];
	if (argc < 2)
		return;
	for (i = 1; i < argc; i++)
	{
		switch(argv[i][0])
		{
		case 'r':
		case 'R':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			bvalues->Reynolds = atof(tmp);
			break;
		}
		case 'w':
		case 'W':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			bvalues->Weber = atof(tmp);
			break;
		}
		case 'F':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			bvalues->Froude = atof(tmp);
			break;
		}
		case 'x':
		case 'X':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			LEVELmax = atoi(tmp);
			break;
		}
		case 'n':
		case 'N':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			LEVELmin = atoi(tmp);
			break;
		}
		case 't':
		case 'T':
		{
			switch(argv[i][1])
			{
			case 'e':
			case 'E':
			{
				for(j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				bvalues->timeend = atof(tmp);
				break;
			}
			case 's':
			case 'S':
			{
				for(j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				bvalues->timestep = atof(tmp);
				break;
			}
			}
			break;
		}
		}
	}
}

/* INTERNAL DIAGNOSTIC FORMATTER
 * Converts elapsed wall-clock seconds into d:hh:mm:ss text. This function has
 * no influence on physical time, adaptive dt, fields or numerical results. */
int timecalculation(double t, char *chartime)
{
	int d, h, m, s;
	if(t < 60.0)
	{
		d = 0;
		h = 0;
		m = 0;
		s = (int) t;
	}
	else if(t < 3600.0)
	{
		d = 0;
		h = 0;
		m = (int) (t / 60.0);
		s = (int) (t - m*60.0);
	}
	else if(t < 3600.0*24.0)
	{
		d = 0;
		h = (int) (t / 3600.0);
		m = (int) ((t - h*3600.0) / 60.0);
		s = (int) (t - h*3600.0 - m*60.0);
	}
	else
	{
		d = (int) (t / 3600.0 / 24.0);
		h = (int) ((t - d*3600.0*24.0) / 3600.0);
		m = (int) ((t - d*3600.0*24.0 - h*3600.0) / 60.0);
		s = (int) (t - d*3600.0*24.0 - h*3600.0 - m*60.0);
	}
	sprintf(chartime, "%d:%02d:%02d:%02d", d, h, m, s);
	return 1;
}
