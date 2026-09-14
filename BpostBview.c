#include "constants.h"
#define BVIEWFOLDER	BVIEW_DIRECTORY
double bview_time_1file, bview_time_total;
double SCALE = 1.5, YPOSITION = 0.50, NOPIXELS = 2000; // it was 1000 before 
double j;
#include "view.h"

// install bview: http://basilisk.fr/src/gl/INSTALL
// qcc -Wall -O2 BpostBview.c -o BpostBview -L$BASILISK/gl -lglutils -lfb_osmesa -lGLU -lOSMesa -lm
//ffmpeg -framerate 5 -pattern_type glob -i '*.png' -c:v libx264 -r 5 -pix_fmt yuv420p output.mp4

#if dimension == 3
#include "lambda2.h"
#endif
struct CFDValues cfdbv;

double total_time_1file, total_time_total;
clock_t simulation_str_time, simulation_end_time;

double Time_BGN =0.00, Time_STP = 0.01, Time_END = 1.50; // Time_STP = 0.01 it was this value before this okay i  have changed it okay 

void readfromargPOST(int argc, char **argv);

int main(int argc, char **argv)
{
	simulation_str_time = clock();
	bview_time_1file = 0.0;
	bview_time_total = 0.0;
	total_time_1file = 0.0;
	total_time_total = 0.0;
	cfdbv.diameter = 1.;
	cfdbv.embed_plane_x = EMBED_PLANE_X;
	readfromargPOST(argc, argv);
	if (argc < 4)
	{
		printf("Run like this:\r\n./BpostBview tb0.00 ts0.01 te1.50\r\nStart Time, Step Size, End Time\r\n");
		return 1;
	}
	printf("Case: %s\r\n", CASE_ID);
	printf("Time_BGN-%f__Time_STP-%f__Time_END-%f\r\n", Time_BGN, Time_STP, Time_END);
	printf("Bview geometry: flat EMBED plane at nominal x=%.17g, liquid wetting %.0f deg\r\n",
	       cfdbv.embed_plane_x, CONTACT_ANGLE_DEGREES);
	size(DOMAIN_WIDTH);
	origin(0., 0.);
	run();
	return 0;
}

void readfromargPOST(int argc, char **argv)
{
	int i, j;
	char tmp[100];
	if (argc < 2)
		return;
	for (i = 1; i < argc; i++)
	{
		switch (argv[i][0])
		{
		case 's':
		case 'S':
		{
			for (j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			SCALE = atof(tmp);
			break;
		}
		case 'y':
		case 'Y':
		{
			for (j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			YPOSITION = atof(tmp);
			break;
		}
		case 'p':
		case 'P':
		{
			for (j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			NOPIXELS = (int)atof(tmp);
			break;
		}
		case 't':
		case 'T':
		{
			switch (argv[i][1])
			{
			case 'b':
			case 'B':
			{
				for (j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				Time_BGN = atof(tmp);
				break;
			}
			case 's':
			case 'S':
			{
				for (j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				Time_STP = atof(tmp);
				break;
			}
			case 'e':
			case 'E':
			{
				for (j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				Time_END = atof(tmp);
				break;
			}
			}
			break;
		}
		}
	}
}

event defaults(i = 0)
{
	interfaces = list_add(NULL, f);
	interfaces = list_add(interfaces, fdrop);
}

event loadfiles(i = 0)
{
	int iloop;
	const double iloopmax = (Time_END - Time_BGN) / Time_STP + 1.0001;
	char nameloadfile[500];
	char ETL[500], TMThis1[500], TMTotal[500];
	double tc, estimatetimeleft;
	int cellnumber;
	scalar varLeft[], varRight[];
	scalar alpha[], kappa[], omega[];
	vector vn[];
	clock_t timestr, timeend, timestrtotal, timeendtotal;
	char BVThis1[500], BVTotal[500];
    
	strcpy(BVThis1, "mkdir -p ");
	strcat(BVThis1, BVIEWFOLDER);
	system(BVThis1);

	timestrtotal = clock();
	for (tc = Time_BGN, iloop = 0; tc <= Time_END + Time_STP * 0.01; tc += Time_STP, iloop++)
	{
		snprintf(nameloadfile, sizeof(nameloadfile), SNAPSHOT_PREFIX "%5.4f", tc);
		if (!restore(file = nameloadfile)) {
			fprintf(stderr, "snapshot not found; stopping cleanly: %s\n",
			        nameloadfile);
			break;
		}

			/* Basilisk does not dump embedded face fractions. Reconstruct exactly
			 * the fixed planar surface used by Bdropimpact.c. Native Basilisk
			 * convention: phi < 0 is solid and phi > 0 is fluid. */
			const double eps = L0/(1 << LEVELmax)/1000.;
			const double xplane = cfdbv.embed_plane_x + eps;
			solid(cs, fs, x - xplane);
			fractions_cleanup(cs, fs);
		cm_update(cm, cs, fs);
		fm_update(fm, cs, fs);
		restriction({cs, fs, cm, fm});
		boundary({cs, fs, cm, fm});

		vorticity(u, omega);
		curvature(f, kappa);
		reconstruction(f, vn, alpha);
		cellnumber = 0;
		foreach ()
			cellnumber++;
		printf("==========----------==========----------==========\r\n");
		printf("time: %.4f ****** cell number: %d\r\n", tc, cellnumber);
		printf("==========----------==========----------==========\r\n");

		timestr = clock();
		printf("==========----------==========----------==========\r\n");
		printf("time: %.4f ****** bview image producing.\r\n", tc);
		char nameBview[500], textBview[500];

		// Vorticity calculation
		foreach () {
		    varLeft[] = f[];// top panel
		    varRight[] = f[]; // bottom  grid panel
		}
		boundary({varLeft});
		boundary({varRight});
		sprintf (nameBview, "%s/out-bview-" CASE_ID "-VOF-%09d.png", BVIEWFOLDER, (int)(round (fabs (tc) * 1000000)));
        //view (width = 1.5 * NOPIXELS, height = NOPIXELS, sx = SCALE, sy = SCALE, fov = 20, tx = -0.48, ty = 0.0); // woking 26 april 2025
		view (width = 1398, height = 1020, quat = { 0, 0, -0.707, 0.707 }, sx = SCALE, sy = SCALE, fov = 17.9799, ty = -0.45);
        // --- CAMERA ANIMATION AND CONDITIONAL VIEW ---
		/*view( width = 1.5 * NOPIXELS,
			height = NOPIXELS,
			sx = SCALE,
			sy = SCALE,
			fov = 6,
			tx = -0.41,
			ty = 0.0
		);*/

        // ============================================================
        // AXIS / DOMAIN ORIENTATION INDICATOR FOR BVIEW IMAGE
        // ============================================================
        box(lw = 2);
		clear();
		draw_vof("f");  // Yellow color, thick lines   lc = {1, 1, 0},
		squares("varLeft", linear = false, min = 0.0, max = 1.0);
         /* EMBED wall fill colour (RGB values range from 0 to 1).
          * Active grey:  fc = {0.3, 0.3, 0.3}
          * Red:          fc = {1.0, 0.0, 0.0}
          * Green:        fc = {0.0, 0.8, 0.0}
          * Blue:         fc = {0.0, 0.0, 1.0}
          * Orange:       fc = {1.0, 0.4, 0.0}
          * Yellow:       fc = {1.0, 1.0, 0.0}
          */
         draw_vof ("cs", "fs", filled = -1, fc = {0.3, 0.3, 0.3});
		mirror({0, 1})
		{
			draw_vof("varRight"); //lc = {1, 1, 0},
			cells(lw = 1);
			//squares("varLeft", linear = true, min = 0.0, max = 1.0); //
            //box(lw = 2);
             /* Keep the mirrored EMBED wall colour identical to the main panel.
              * Use any RGB alternative listed above for both draw_vof() calls. */
             draw_vof ("cs", "fs", filled = -1, fc = {0.3, 0.3, 0.3});
		}

		sprintf(textBview, "%s, time t*: %.4f", CASE_ID, tc);
		draw_string(textBview, pos = 4, lw = 3, size = 80);
		save(nameBview);
        save(BVIEW_MOVIE_BASENAME);

		printf("done!\r\n");
		printf("==========----------==========----------==========\r\n");
		timeend = clock();
		bview_time_1file += (double)(timeend - timestr) / CLOCKS_PER_SEC;

		bview_time_total += bview_time_1file;
		timecalculation(bview_time_1file, BVThis1);
		timecalculation(bview_time_total, BVTotal);

		timeendtotal = clock();
		total_time_1file += (double)(timeendtotal - timestrtotal) / CLOCKS_PER_SEC;
		total_time_total += total_time_1file;

		estimatetimeleft = total_time_total * iloopmax / (iloop + 1.0) - total_time_total;
		timecalculation(estimatetimeleft, ETL);
		timecalculation(total_time_1file, TMThis1);
		timecalculation(total_time_total, TMTotal);

		printf("==========----------==========----------==========\r\n");
		printf("LAST FILE DURATIONS:\r\n");
		printf("bview duration: %s\r\n", BVThis1);
		printf("total duration: %s\r\n", TMThis1);
		printf("ALL FILES DURATIONS UNTIL NOW:\r\n");
		printf("bview total duration: %s\r\n", BVThis1);
		bview_time_1file = 0.0;
		printf("total duration: %s\r\n", TMTotal);
		total_time_1file = 0.0;
		printf("==========----------==========----------==========\r\n");
		printf("Estimated time left: %s\r\n", ETL);
		printf("done with t = %.4f\r\n", tc);
		printf("==========----------==========----------==========\r\n");
	}
}

event end(i = 0)
{
}
