/* geofac
 *  USAGE:
 *  	geofac [-r]
 *  		-r    use a known seed value (SEED)
 */

/*
*	gf_MC.c
*	jkl 3/10/93
*	
*	This program calculates the differential geometry factor for
*	a detector consisting of any number of planes of either rectangular
*	or circular geometry.  The only requirement is that the detectors
*	be centered relative to each other.  
*
*	outline:
*       input:  number of detectors;
*		for each:
*		if rectangle x-side y-side z-coord.
*		else         radius   0	   z-coord.
*
*	The number of iterations specified will result in that number of
*	"hits",  the number of tries in the Monte Carlo run will be
*	somewhat greater.
*
*	The output will be numeric only so that it can be used
*	directly in another routine:
*	for each detector:  x_side y_side z_coordinate,
*		(y_side=0 indicates circle of radius x_side.)
*		number_of_hits  number_of_tries
*
*		table of zenith angle bins (5 degree separation) and
*			differential geometry factors.
*
*	Total geometry factor.
*
*	usage: gf_MC.
*
*	geofac - fixed a few bugs, re-indented code a little.
*		Marty Olevitch, 9/23/03
*/




//


#include <stdio.h>
#include <math.h>
#include <stdlib.h>	// getopt
#include <unistd.h>
#include <string.h>
#include <iostream>

#define		PI	3.14159265
#define		M	2147483647.0

#define SEED		1227	// non-random seed
#define MAXLAYERS	10
#define NHITS		90
#define N_SECBINS	50
#define NBINS		18

FILE *infile;
FILE *outfile;

void
usage(void)
{
    fprintf(stderr, "USAGE: geofac [-r]\n");
}

double
rand_no(void)
{
    return rand()/M;
}

int
main(int argc, char *argv[])
{
    // for getopt
    int c;
    extern char *optarg;
    extern int optind;

    int theta_bin;
	int phi_bin;
    int circle[MAXLAYERS];
    double cos_phi;
    long hits;
    int i;
    int j;
    int k;
    int m;
    int missed_detector;
    double min_area;
    int min_index;
    long n_good;
    long n_hit[NHITS][4*NHITS];
    int n_layer;
    long n_try;
    double phi;
    double rx;
    double ry;
    int sec_angle_bin;
    long sec_n_hit[N_SECBINS];
    double sec_theta;
    int	seed;
    double sin_phi;
    double tan_theta;
    double theta;
    int use_non_random_seed = 0;
    double x[MAXLAYERS];
    double x_hit;
    double y[MAXLAYERS];
    double y_hit;
    double z[MAXLAYERS];
    double z_rel[MAXLAYERS];
	
		// Begin variables added by Brian
	double gamma_r;
	double gamma_d;
	double delta_d;
	double delta_r;
	long n_hit_gamma[2*NHITS][2*NHITS];
	double z_arb=10.0;
	int gamma_bin;
	int delta_bin;
	long sum_gamma[2*NBINS][2*NBINS];
    long n_hit_theta_gamma[2*NHITS][NHITS];
    long sum_theta_gamma[2*NHITS][NHITS];

    double theta_dgf_sum;
    double theta_dgf_sum_total=0.0;
    double ang_off=0.0;
    double cos_phi_off;
    double sin_phi_off;
	
	int j_2;
	int k_2;
	int m_2;
	
    
	char extension[20];
	char outfile_name[40];

        char filename[100];
    
		// End variable added by Brian
	
	
	
    while ((c = getopt(argc, argv, "r")) != EOF) {
	switch (c) {
	    case 'r':
		use_non_random_seed = 1;
		break;
	    case '?':
	    default:
		usage();
		exit(1);
	}
    }

    if (use_non_random_seed) {
	seed = SEED;
    } else {
	seed = getpid();
    }
    srand(seed);

    /* dimensions of detector */

    strcpy(filename,"angle_in.txt\0");
    if((infile = fopen(filename,"r")) == NULL)
      {    fprintf(stderr,"Error opening %s for input\n",filename);
      exit(0);    }
    fprintf(stderr,"Reading from %s\n",filename);


       fscanf(infile,"%d", &n_layer);
/*
    fprintf(stderr,"Enter number of layers: ");
    fflush(stderr);
    if(scanf("%d", &n_layer) == EOF) {
	exit(1);
    }
    if (n_layer > MAXLAYERS) {
	fprintf(stderr, "Sorry maximum number of layers is %d, not %d\n",
	    MAXLAYERS, n_layer);
	exit(1);
    }

    fprintf(stderr,"For each of the %d detectors enter either\n",n_layer);
    fprintf(stderr,"  x-side y-side z-coordinate --if rectangular\n");
    fprintf(stderr,"  radius    0   z-coordinate --if circular\n");
*/

    for(i=0; i<n_layer; i++) {
        /*
	fprintf(stderr, "Layer %d. ", i+1);
	fflush(stderr);
        */
	if(fscanf(infile,"%lf %lf %lf",&x[i],&y[i],&z[i]) != 3)	{
	    exit(1);
	}
	if (y[i] == 0) {
	    circle[i] = 1;
	} else {
	    circle[i] = 0;
	    x[i] /= 2.0;
	    y[i] /= 2.0;
	}
    }

    /*
    fprintf(stderr,"Number of iterations (hits)? ");
    fflush(stderr);
    */
    fscanf(infile,"%ld",&hits);


    fscanf(infile,"%lf",&ang_off);
    /*    
    fprintf(stderr,"output file extension (10 characters of less)? ");
	fflush(stderr);
    */
	fscanf(infile,"%s",extension);
	
	fprintf(stderr,"output file extension: %s",extension);
	fflush(stderr);

  fclose(infile);
    
    /* creating log file for Monte Carlo inputs */
    
    sprintf(outfile_name,"run_log%s.txt",extension);
    
    if((outfile = fopen(outfile_name,"w")) == NULL)
    {    fprintf(stderr,"Error opening %s for output\n",outfile_name);
        exit(0);
    }
    
    fprintf(outfile,"number of layers: %2d\n", n_layer);
    for(i=0; i<n_layer; i++) {
        if(circle[i] == 1) {
            fprintf(outfile,"%6lf\t%6lf\t%6lf\n", x[i], y[i], z[i]);
        } else {
            fprintf(outfile,"%6lf\t%6lf\t%6lf\n", 2*x[i], 2*y[i], z[i]);
        }
    }
    fprintf(outfile,"number of hits: %ld\n", hits);
    fprintf(outfile,"output file extension: %s",extension);

    
    /* Find smallest detector */
    if (circle[0] == 1) {
	min_area = PI * x[0] * x[0];
    } else {
	min_area = x[0] * y[0] * 4;
    }

    /* determine smallest detector */
    min_index = 0;
    for (i=1; i<n_layer; i++) {
	double area;
	if (circle[i] == 1) {
	    area = PI * x[i] * x[i];
	} else {
	    area = x[i] * y[i] * 4;
	}
	if (area < min_area) {
	    min_area = area;
	    min_index = i;
	}
    }

    /* calculate z-positions from minimal detector */
    for (i=0; i<n_layer; i++) {
	z_rel[i] = z[i] - z[min_index];
    }

    printf("\nlayer\tx\ty\tz\n");
    for(i=0; i<n_layer; i++) {
	printf("%2d\t%.2lf\t%.2lf\t%.2lf\n", i+1, 2 * x[i], 2 * y[i], z[i]);
    }

    for(j=0; j<NHITS; j++) {
		for (j_2=0; j_2 < 4*NHITS; j_2++) {

			n_hit[j][j_2] = 0.0;
		}
    }
	
    for(j=0; j<N_SECBINS; j++) {
	sec_n_hit[j] = 0.0;
    }
	
    for(j=0; j<NHITS; j++) {
		for (j_2=0; j_2<2*NHITS; j_2++) {
			n_hit_gamma[j][j_2] = 0.0;
		}
    }	

    for(j=0; j<2*NHITS; j++) {
		for (j_2=0; j_2<NHITS; j_2++) {
			sum_theta_gamma[j][j_2] = 0.0;
		}
    }
    
    n_good = 0;
    n_try = 0;
    while(n_good < hits) {
	n_try++;
        
        if (n_good%(int)(hits/100)==0) {
            fprintf(stderr,"n_good: %ld n_try: %ld\n",n_good,n_try);
        }
	/* for greater efficiency choose points on smallest detector */
	if (circle[min_index]) {
	    double phi_0;
	    double r;
	    r = x[min_index] * sqrt(rand_no());
	    phi_0 = 2 * PI * rand_no();
	    rx = r * cos(phi_0);
	    ry = r * sin(phi_0);
	} else {
	    ry = (y[min_index] * 2 * (rand_no())) - y[min_index];
	    rx = (x[min_index] * 2 * (rand_no())) - x[min_index];
	}

	/* choose angles from isotropic distribution */
	theta = acos( sqrt(rand_no()) );
	theta_bin = (int)( (theta*180) / PI );

	tan_theta = tan(theta);
	sec_theta = 1.0 / cos(theta);
	sec_angle_bin = (sec_theta - 1) * 10;

			//	fprintf(stderr, "sec theta %f\n", sec_theta);
			//	fprintf(stderr, "theta %2d\tcos %f\tsec %f\n",
			//	    angle_bin, cos(theta), sec_theta);
	
	phi = 2 * PI * (rand_no());
	cos_phi = cos(phi);
        cos_phi_off = cos(phi + ang_off*PI/180.0);
	sin_phi = sin(phi);
	sin_phi_off = sin(phi + ang_off*PI/180.0);

		phi_bin = (int)( (phi*180) /PI );
		
	/*  Check if trajectory hits lowest detectors;
	 *  if so, continue checking;
	 *  if not, make new choice.
	 */

	missed_detector = 0;
	for(i=0; i<min_index; i++) {
	    x_hit = rx + (z_rel[i] * tan_theta * cos_phi);
	    y_hit = ry + (z_rel[i] * tan_theta * sin_phi);
	    if(circle[i]) {
		if ( ((x_hit * x_hit) + (y_hit * y_hit)) > x[i] * x[i]) {
		    missed_detector = 1;
		    break;
		}
	    } else {
		if ((x_hit < -x[i]) || (x_hit > x[i]))  {
		    missed_detector = 1;
		    break;
		}
		if ((y_hit < -y[i]) || (y_hit > y[i]))  {
		    missed_detector = 1;
		    break;
		}
	    }
	}
	if (missed_detector) {
	    continue;
	}

		
	/* same for upper detectors */

	    for(i=min_index+1; i<n_layer; i++) {
		x_hit = rx + (z_rel[i] * tan_theta * cos_phi);
		y_hit = ry + (z_rel[i] * tan_theta * sin_phi);

		if(circle[i]) {
		    if( ((x_hit * x_hit) + (y_hit * y_hit)) > x[i] * x[i]) {
			    missed_detector = 1;
			    break;
		    }
		} else {
		    if ((x_hit < -x[i]) || (x_hit > x[i]))  {
			missed_detector = 1;
			break;
		    }
		    if ((y_hit < -y[i]) || (y_hit > y[i]))  {
			missed_detector = 1;
			break;
		    }
		}
	    }
	    if (missed_detector) {
		continue;
	    }

	    n_hit[theta_bin][phi_bin]++;   /* Got a hit! Record it. */
	    if (sec_angle_bin < N_SECBINS) {
		// ignores angles greater than 80 degrees for sec theta hits
		sec_n_hit[sec_angle_bin]++;
	    }
		
		
		gamma_r=atan2(z_arb,z_arb * tan_theta * cos_phi_off);
		gamma_d=gamma_r*180/PI;
			// fprintf(stderr,"gamma %.4lf %.2lf\n",gamma_r,gamma_d);
		
		gamma_bin = (int)( (gamma_r*180) / PI );
		
		delta_r=atan2(sqrt(pow(z_arb * tan_theta * cos_phi_off, 2.0)+pow(z_arb, 2.0)),z_arb * tan_theta * sin_phi_off);
		delta_d=delta_r*180/PI;
		
		delta_bin = (int)( (delta_r*180) / PI );
		
		n_hit_gamma[gamma_bin][delta_bin]++;
		
        n_hit_theta_gamma[gamma_bin][theta_bin]++;
		
	    n_good++;
    }

    fprintf(outfile,"\nhits\ttries\n");
    fprintf(outfile,"%ld\t%ld\n\n", hits, n_try);


    {
	/* Calculate geometry factors */
	long sum[NBINS][NBINS*4];
	double area;
	double diff_geo_fac;
	double tot;
		int check_count=0;


			// summing up for 5 degree angle bins
	for(j=0; j<NBINS; j++) {
		for (j_2=0; j_2<4*NBINS; j_2++) {

	    sum[j][j_2] = 0;
			for(k=0; k<5; k++) {
				for (k_2=0; k_2<5; k_2++) {
					m = k + (j*5);
					m_2 = k_2 + (j_2*5);
					sum[j][j_2] += n_hit[m][m_2];
					check_count += n_hit[m][m_2];
				}
			}
			//fprintf(stderr,"sum[%d][%d] %ld\n",j,j_2,sum[j][j_2]);
		}

	}
	fprintf(stderr,"check_count %d\n",check_count);
		
		
			// summing up for 5 degree angle bins
		for(j=0; j<2*NBINS; j++) {
			for (j_2=0; j_2<2*NBINS; j_2++) {
				
				sum_gamma[j][j_2] = 0;
				for(k=0; k<5; k++) {
					for (k_2=0; k_2<5; k_2++) {
						m = k + (j*5);
						m_2 = k_2 + (j_2*5);
						sum_gamma[j][j_2] += n_hit_gamma[m][m_2];
						check_count += n_hit_gamma[m][m_2];
					}
				}
				//fprintf(stderr,"sum[%d][%d] %ld\n",j,j_2,sum_gamma[j][j_2]);
			}
			
		}
		fprintf(stderr,"check_count %d\n",check_count);

		
	/* Total geometry factor */
	area = (PI * (double)hits * min_area) / (double)n_try;
	printf("total geometry factor: ");
	printf("%.2lf\n", area);

    // last bit for the log file
    fprintf(outfile,"total geometry factor: ");
    fprintf(outfile,"%.2lf\n", area);
    fclose(outfile);
			// Handling output file
		
        // Handling output file
		sprintf(outfile_name,"dgf_2d%s.txt",extension);
        
		if((outfile = fopen(outfile_name,"w")) == NULL)
		{    fprintf(stderr,"Error opening %s for output\n","dgf_2d.txt");
			exit(0);    }

	/* differential geometry factors */
	//printf("\nangle\tgeo_factor\n");
	tot = 0;
	for(j=0; j<NBINS; j++) {
		for (j_2=0; j_2<4*NBINS; j_2++) {
			diff_geo_fac = ( (double)sum[j][j_2] / hits) * area;
			tot += diff_geo_fac;
			//printf("%4.1lf\t%4.1lf\t%10.4lf\n", (double)(2.5 + (j*5)), (double)(2.5 + (j_2*5)), diff_geo_fac);
			fprintf(outfile,"%4.1lf\t%4.1lf\t%10.4lf\n", (double)(2.5 + (j*5)), (double)(2.5 + (j_2*5)), diff_geo_fac);

		}
	}
	printf("total\t%10.4f\n", tot);
		fclose(outfile);
		
		
			// Generate 1 degree by 1 degree binned differential geometry factors
			// Handling output file
		
        sprintf(outfile_name,"dgf_2d_bin1%s.txt",extension);

        
		if((outfile = fopen(outfile_name,"w")) == NULL)
		{    fprintf(stderr,"Error opening %s for output\n","dgf_2d_bin1.txt");
			exit(0);    }
		
		/* differential geometry factors */
		printf("\nangle\tgeo_factor\n");
		tot = 0;
		for(j=0; j<NHITS; j++) {
			for (j_2=0; j_2<4*NHITS; j_2++) {
				diff_geo_fac = ( (double)n_hit[j][j_2] / hits) * area;
				tot += diff_geo_fac;
				//printf("%4.1lf\t%4.1lf\t%10.4lf\n", (double)(0.5 + j), (double)(0.5 + j_2), diff_geo_fac);
				fprintf(outfile,"%4.1lf\t%4.1lf\t%10.4lf\n", (double)(0.5 + j), (double)(0.5 + j_2), diff_geo_fac);
				
			}
		}
		printf("total\t%10.4f\n", tot);
		fclose(outfile);
		

			// Generate 1 dimensional 1 degree binned differential geometry factors
			// Handling output file
		
        sprintf(outfile_name,"dgf_1d_bin1%s.txt",extension);

        
		if((outfile = fopen(outfile_name,"w")) == NULL)
		{    fprintf(stderr,"Error opening %s for output\n","dgf_1d_bin1.txt");
			exit(0);    }
		
		/* differential geometry factors */
		printf("\nangle\tgeo_factor\n");
		tot = 0;
		for(j=0; j<NHITS; j++) {
                        theta_dgf_sum=0;
			for (j_2=0; j_2<4*NHITS; j_2++) {
				diff_geo_fac = ( (double)n_hit[j][j_2] / hits) * area;
                                theta_dgf_sum += diff_geo_fac;
                                theta_dgf_sum_total += diff_geo_fac;
				tot += diff_geo_fac;
				//printf("%4.1lf\t%4.1lf\t%10.4lf\n", (double)(0.5 + j), (double)(0.5 + j_2), diff_geo_fac);
				//fprintf(outfile,"%4.1lf\t%4.1lf\t%10.4lf\n", (double)(0.5 + j), (double)(0.5 + j_2), diff_geo_fac);
				
			}
			fprintf(outfile,"%4.1lf\t%10.4lf\n", (double)(0.5 + j), theta_dgf_sum);
		}
		printf("total\t%10.4f\t theta_dgf_sum_total\t%10.4f\n", tot, theta_dgf_sum_total);
		fclose(outfile);
		
		
		
		/* differential geometry factors */
		double diff_geo_fac_gamma;

			// Handling output files

        sprintf(outfile_name,"dgf_2d_g%s.txt",extension);

		if((outfile = fopen(outfile_name,"w")) == NULL)
		{    fprintf(stderr,"Error opening %s for output\n",outfile_name);
			exit(0);    }
		
		printf("\nangle\tgeo_factor\n");
		tot = 0;
		for(j=0; j<2*NBINS; j++) {
			for (j_2=0; j_2<2*NBINS; j_2++) {
				diff_geo_fac_gamma = ( (double)sum_gamma[j][j_2] / hits) * area;
				tot += diff_geo_fac_gamma;
				//printf("%4.1lf\t%4.1lf\t%10.4lf\n", (double)(2.5 + (j*5)),(double)(2.5 + (j_2*5)), diff_geo_fac_gamma);
				fprintf(outfile,"%4.1lf\t%4.1lf\t%10.4lf\n", (double)(2.5 + (j*5)), (double)(2.5 + (j_2*5)), diff_geo_fac_gamma);

			}
		}
		printf("total\t%10.4f\n", tot);
		fclose(outfile);

			// 1 degree by 1 degree binned differential geometry factors
			// Handling output file
		
        sprintf(outfile_name,"dgf_2d_g_bin1%s.txt",extension);

        
		if((outfile = fopen(outfile_name,"w")) == NULL)
		{    fprintf(stderr,"Error opening %s for output\n",outfile_name);
			exit(0);    }
		
		printf("\nangle\tgeo_factor\n");
		tot = 0;
		for(j=0; j<2*NHITS; j++) {
			for (j_2=0; j_2<2*NHITS; j_2++) {
				diff_geo_fac_gamma = ( (double)n_hit_gamma[j][j_2] / hits) * area;
				tot += diff_geo_fac_gamma;
				//printf("%4.1lf\t%4.1lf\t%10.4lf\n", (double)(0.5 + j),(double)(0.5 + j_2), diff_geo_fac_gamma);
				fprintf(outfile,"%4.1lf\t%4.1lf\t%10.4lf\n", (double)(0.5 + j), (double)(0.5 + j_2), diff_geo_fac_gamma);
				
			}
		}
		printf("total\t%10.4f\n", tot);
		fclose(outfile);

        // 1 degree by 1 degree binned differential geometry factors for gamma by theta
        // Handling output file
		
        sprintf(outfile_name,"dgf_2d_g_theta_bin1%s.txt",extension);
        
        
		if((outfile = fopen(outfile_name,"w")) == NULL)
		{    fprintf(stderr,"Error opening %s for output\n",outfile_name);
			exit(0);    }
		
		printf("\nangle\tgeo_factor\n");
		tot = 0;
		for(j=0; j<2*NHITS; j++) {
			for (j_2=0; j_2<NHITS; j_2++) {
				diff_geo_fac_gamma = ( (double)n_hit_theta_gamma[j][j_2] / hits) * area;
				tot += diff_geo_fac_gamma;
				//printf("%4.1lf\t%4.1lf\t%10.4lf\n", (double)(0.5 + j),(double)(0.5 + j_2), diff_geo_fac_gamma);
				fprintf(outfile,"%4.1lf\t%4.1lf\t%10.4lf\n", (double)(0.5 + j), (double)(0.5 + j_2), diff_geo_fac_gamma);
				
			}
		}
		printf("total\t%10.4f\n", tot);
		fclose(outfile);
        
        
        
	/* sec theta geometry factors */
	printf("\nsecant\tgeo_factor\n");
	tot = 0;
	for (j=0; j<N_SECBINS; j++) {
	    diff_geo_fac = ( (double)sec_n_hit[j] / hits) * area;
	    tot += diff_geo_fac;
	    printf("%4.1lf\t%10.4lf\n", (j/10.0) + 1, diff_geo_fac);
	}
	printf("total\t%10.2f\n", tot);
    }

    return 0;
}

