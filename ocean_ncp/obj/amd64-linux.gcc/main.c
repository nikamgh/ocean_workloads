#line 233 "/mnt/panzer/mnika/Workload_Pool/parsec-3.0/bin/../pkgs/libs/parmacs/inst/amd64-linux.gcc/m4/parmacs.pthreads.c.m4"

#line 1 "main.C"
/*************************************************************************/
/*                                                                       */
/*  SPLASH Ocean Code                                                    */
/*                                                                       */
/*  This application studies the role of eddy and boundary currents in   */
/*  influencing large-scale ocean movements.  This implementation uses   */
/*  statically allocated two-dimensional arrays for grid data storage.   */
/*                                                                       */
/*  Command line options:                                                */
/*                                                                       */
/*     -nN : Simulate NxN ocean.  N must be (power of 2)+2.              */
/*     -pP : P = number of processors.  P must be power of 2.            */
/*     -eE : E = error tolerance for iterative relaxation.               */
/*     -rR : R = distance between grid points in meters.                 */
/*     -tT : T = timestep in seconds.                                    */
/*     -s  : Print timing statistics.                                    */
/*     -o  : Print out relaxation residual values.                       */
/*     -h  : Print out command line options.                             */
/*                                                                       */
/*  Default: OCEAN -n130 -p1 -e1e-7 -r20000.0 -t28800.0                  */
/*                                                                       */
/*  NOTE: This code works under both the FORK and SPROC models.          */
/*                                                                       */
/*************************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "decs.h"
//#if defined(ZSIM_TRACE_1) || defined(ZSIM_TRACE_2) || defined(ZSIM_TRACE_3) || defined(ZSIM_TRACE_4)
    #include "zsim_hooks.h"
//#endif
#ifdef ENABLE_PARSEC_HOOKS
#include <hooks.h>
#endif

#line 37
#include <pthread.h>
#line 37
#include <sys/time.h>
#line 37
#include <unistd.h>
#line 37
#include <stdlib.h>
#line 37
#define MAX_THREADS 1024
#line 37
pthread_t PThreadTable[MAX_THREADS];
#line 37


#define DEFAULT_N      258
#define DEFAULT_P        1
#define DEFAULT_E        1e-7
#define DEFAULT_T    28800.0
#define DEFAULT_R    20000.0
#define PAGE_SIZE     4096

struct global_struct *global;
struct fields_struct *fields;
struct fields2_struct *fields2;
struct wrk1_struct *wrk1;
struct wrk3_struct *wrk3;
struct wrk2_struct *wrk2;
struct wrk4_struct *wrk4;
struct wrk6_struct *wrk6;
struct wrk5_struct *wrk5;
struct frcng_struct *frcng;
struct iter_struct *iter;
struct guess_struct *guess;
struct multi_struct *multi;
struct locks_struct *locks;
struct bars_struct *bars;

long startcol[2][INPROCS];
long nprocs = DEFAULT_P;
long startrow[2][INPROCS];
double h1 = 1000.0;
double h3 = 4000.0;
double h = 5000.0;
double lf = -5.12e11;
double eps = 0;
double res = DEFAULT_R;
double dtau = DEFAULT_T;
double f0 = 8.3e-5;
double beta = 2.0e-11;
double gpr = 0.02;
long im = DEFAULT_N;
long jm;
double tolerance = DEFAULT_E;
double eig2;
double ysca;
long jmm1;
double pi;
double t0 = 0.5e-4 ;
double outday0 = 1.0;
double outday1 = 2.0;
double outday2 = 2.0;
double outday3 = 2.0;
double factjacob;
double factlap;
long numlev;
long minlev;
long imx[MAX_LEVELS];
long jmx[MAX_LEVELS];
double lev_res[MAX_LEVELS];
double lev_tol[MAX_LEVELS];
double maxwork = 10000.0;

struct Global_Private *gp;

double i_int_coeff[MAX_LEVELS];
double j_int_coeff[MAX_LEVELS];
long xprocs;
long yprocs;
long do_stats = 0;
long do_output = 0;

int main(int argc, char *argv[])
{
#ifdef ENABLE_PARSEC_HOOKS
	__parsec_bench_begin (__splash2_ocean_ncp);
#endif
   long i;
   long j;
   long xextra;
   long xportion;
   long yextra;
   long yportion;
   long lower;
   double procsqrt;
   long k;
   long logtest;
   long my_num;
   unsigned long computeend;
   double min_total;
   double max_total;
   double avg_total;
   double min_multi;
   double max_multi;
   double avg_multi;
   double min_frac;
   double max_frac;
   double avg_frac;
   extern char *optarg;
   long ch;
   unsigned long start;
   double *tmp_ptr;

   {
#line 137
	struct timeval	FullTime;
#line 137

#line 137
	gettimeofday(&FullTime, NULL);
#line 137
	(start) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 137
}

   while ((ch = getopt(argc, argv, "n:p:e:r:t:soh")) != -1) {
     switch(ch) {
     case 'n': im = atoi(optarg);
               if (im > IMAX) {
                 printerr("Max grid size exceeded\n");
                 exit(-1);
               }
               if (log_2(im-2) == -1) {
                 printerr("Grid must be ((power of 2)+2) in each dimension\n");
                 exit(-1);
               }
               break;
     case 'p': nprocs = atoi(optarg);
               if (nprocs < 1) {
                 printerr("P must be >= 1\n");
                 exit(-1);
               }
               if (log_2(nprocs) == -1) {
                 printerr("P must be a power of 2\n");
                 exit(-1);
               }
               break;
     case 'e': tolerance = atof(optarg); break;
     case 'r': res = atof(optarg); break;
     case 't': dtau = atof(optarg); break;
     case 's': do_stats = !do_stats; break;
     case 'o': do_output = !do_output; break;
     case 'h': printf("Usage: OCEAN <options>\n\n");
               printf("options:\n");
               printf("  -nN : Simulate NxN ocean.  N must be (power of 2)+2.\n");
               printf("  -pP : P = number of processors.  P must be power of 2.\n");
               printf("  -eE : E = error tolerance for iterative relaxation.\n");
               printf("  -rR : R = distance between grid points in meters.\n");
               printf("  -tT : T = timestep in seconds.\n");
               printf("  -s  : Print timing statistics.\n");
               printf("  -o  : Print out relaxation residual values.\n");
               printf("  -h  : Print out command line options.\n\n");
               printf("Default: OCEAN -n%1d -p%1d -e%1g -r%1g -t%1g\n",
                       DEFAULT_N,DEFAULT_P,DEFAULT_E,DEFAULT_R,DEFAULT_T);
               exit(0);
               break;
     }
   }

   {;}

   

   logtest = im-2;
   numlev = 1;
   while (logtest != 1) {
     if (logtest%2 != 0) {
       printerr("Cannot determine number of multigrid levels\n");
       exit(-1);
     }
     logtest = logtest / 2;
     numlev++;
   }

   if (numlev > MAX_LEVELS) {
     printerr("Max grid levels exceeded for multigrid\n");
     exit(-1);
   }

   jm = im;
   printf("\n");
   printf("Ocean simulation with W-cycle multigrid solver\n");
   printf("    Processors                         : %1ld\n",nprocs);
   printf("    Grid size                          : %1ld x %1ld\n",im,jm);
   printf("    Grid resolution (meters)           : %0.2f\n",res);
   printf("    Time between relaxations (seconds) : %0.0f\n",dtau);
   printf("    Error tolerance                    : %0.7g\n",tolerance);
   printf("\n");
//#if defined(ZSIM_TRACE_1) || defined(ZSIM_TRACE_2) || defined(ZSIM_TRACE_3) || defined(ZSIM_TRACE_4)
 //   zsim_roi_begin();
//#endif

   gp = (struct Global_Private *) malloc((nprocs+1)*sizeof(struct Global_Private));;
   for (i=0;i<nprocs;i++) {
     gp[i].multi_time = 0;
     gp[i].total_time = 0;
   }
   global = (struct global_struct *) malloc(sizeof(struct global_struct));;
   fields = (struct fields_struct *) malloc(sizeof(struct fields_struct));;
   fields2 = (struct fields2_struct *) malloc(sizeof(struct fields2_struct));;
   wrk1 = (struct wrk1_struct *) malloc(sizeof(struct wrk1_struct));;
   wrk3 = (struct wrk3_struct *) malloc(sizeof(struct wrk3_struct));;
   wrk2 = (struct wrk2_struct *) malloc(sizeof(struct wrk2_struct));;
   wrk4 = (struct wrk4_struct *) malloc(sizeof(struct wrk4_struct));;
   wrk6 = (struct wrk6_struct *) malloc(sizeof(struct wrk6_struct));;
   wrk5 = (struct wrk5_struct *) malloc(sizeof(struct wrk5_struct));;
   frcng = (struct frcng_struct *) malloc(sizeof(struct frcng_struct));;
   iter = (struct iter_struct *) malloc(sizeof(struct iter_struct));;
   guess = (struct guess_struct *) malloc(sizeof(struct guess_struct));;
   multi = (struct multi_struct *) malloc(sizeof(struct multi_struct));;
   locks = (struct locks_struct *) malloc(sizeof(struct locks_struct));;
   bars = (struct bars_struct *) malloc(sizeof(struct bars_struct));;

   {pthread_mutex_init(&(locks->idlock), NULL);}
   {pthread_mutex_init(&(locks->psiailock), NULL);}
   {pthread_mutex_init(&(locks->psibilock), NULL);}
   {pthread_mutex_init(&(locks->donelock), NULL);}
   {pthread_mutex_init(&(locks->error_lock), NULL);}
   {pthread_mutex_init(&(locks->bar_lock), NULL);}

#if defined(MULTIPLE_BARRIERS)
   {
#line 245
	unsigned long	Error;
#line 245

#line 245
	Error = pthread_mutex_init(&(bars->iteration).mutex, NULL);
#line 245
	if (Error != 0) {
#line 245
		printf("Error while initializing barrier.\n");
#line 245
		exit(-1);
#line 245
	}
#line 245

#line 245
	Error = pthread_cond_init(&(bars->iteration).cv, NULL);
#line 245
	if (Error != 0) {
#line 245
		printf("Error while initializing barrier.\n");
#line 245
		pthread_mutex_destroy(&(bars->iteration).mutex);
#line 245
		exit(-1);
#line 245
	}
#line 245

#line 245
	(bars->iteration).counter = 0;
#line 245
	(bars->iteration).cycle = 0;
#line 245
}
   {
#line 246
	unsigned long	Error;
#line 246

#line 246
	Error = pthread_mutex_init(&(bars->gsudn).mutex, NULL);
#line 246
	if (Error != 0) {
#line 246
		printf("Error while initializing barrier.\n");
#line 246
		exit(-1);
#line 246
	}
#line 246

#line 246
	Error = pthread_cond_init(&(bars->gsudn).cv, NULL);
#line 246
	if (Error != 0) {
#line 246
		printf("Error while initializing barrier.\n");
#line 246
		pthread_mutex_destroy(&(bars->gsudn).mutex);
#line 246
		exit(-1);
#line 246
	}
#line 246

#line 246
	(bars->gsudn).counter = 0;
#line 246
	(bars->gsudn).cycle = 0;
#line 246
}
   {
#line 247
	unsigned long	Error;
#line 247

#line 247
	Error = pthread_mutex_init(&(bars->p_setup).mutex, NULL);
#line 247
	if (Error != 0) {
#line 247
		printf("Error while initializing barrier.\n");
#line 247
		exit(-1);
#line 247
	}
#line 247

#line 247
	Error = pthread_cond_init(&(bars->p_setup).cv, NULL);
#line 247
	if (Error != 0) {
#line 247
		printf("Error while initializing barrier.\n");
#line 247
		pthread_mutex_destroy(&(bars->p_setup).mutex);
#line 247
		exit(-1);
#line 247
	}
#line 247

#line 247
	(bars->p_setup).counter = 0;
#line 247
	(bars->p_setup).cycle = 0;
#line 247
}
   {
#line 248
	unsigned long	Error;
#line 248

#line 248
	Error = pthread_mutex_init(&(bars->p_redph).mutex, NULL);
#line 248
	if (Error != 0) {
#line 248
		printf("Error while initializing barrier.\n");
#line 248
		exit(-1);
#line 248
	}
#line 248

#line 248
	Error = pthread_cond_init(&(bars->p_redph).cv, NULL);
#line 248
	if (Error != 0) {
#line 248
		printf("Error while initializing barrier.\n");
#line 248
		pthread_mutex_destroy(&(bars->p_redph).mutex);
#line 248
		exit(-1);
#line 248
	}
#line 248

#line 248
	(bars->p_redph).counter = 0;
#line 248
	(bars->p_redph).cycle = 0;
#line 248
}
   {
#line 249
	unsigned long	Error;
#line 249

#line 249
	Error = pthread_mutex_init(&(bars->p_soln).mutex, NULL);
#line 249
	if (Error != 0) {
#line 249
		printf("Error while initializing barrier.\n");
#line 249
		exit(-1);
#line 249
	}
#line 249

#line 249
	Error = pthread_cond_init(&(bars->p_soln).cv, NULL);
#line 249
	if (Error != 0) {
#line 249
		printf("Error while initializing barrier.\n");
#line 249
		pthread_mutex_destroy(&(bars->p_soln).mutex);
#line 249
		exit(-1);
#line 249
	}
#line 249

#line 249
	(bars->p_soln).counter = 0;
#line 249
	(bars->p_soln).cycle = 0;
#line 249
}
   {
#line 250
	unsigned long	Error;
#line 250

#line 250
	Error = pthread_mutex_init(&(bars->p_subph).mutex, NULL);
#line 250
	if (Error != 0) {
#line 250
		printf("Error while initializing barrier.\n");
#line 250
		exit(-1);
#line 250
	}
#line 250

#line 250
	Error = pthread_cond_init(&(bars->p_subph).cv, NULL);
#line 250
	if (Error != 0) {
#line 250
		printf("Error while initializing barrier.\n");
#line 250
		pthread_mutex_destroy(&(bars->p_subph).mutex);
#line 250
		exit(-1);
#line 250
	}
#line 250

#line 250
	(bars->p_subph).counter = 0;
#line 250
	(bars->p_subph).cycle = 0;
#line 250
}
   {
#line 251
	unsigned long	Error;
#line 251

#line 251
	Error = pthread_mutex_init(&(bars->sl_prini).mutex, NULL);
#line 251
	if (Error != 0) {
#line 251
		printf("Error while initializing barrier.\n");
#line 251
		exit(-1);
#line 251
	}
#line 251

#line 251
	Error = pthread_cond_init(&(bars->sl_prini).cv, NULL);
#line 251
	if (Error != 0) {
#line 251
		printf("Error while initializing barrier.\n");
#line 251
		pthread_mutex_destroy(&(bars->sl_prini).mutex);
#line 251
		exit(-1);
#line 251
	}
#line 251

#line 251
	(bars->sl_prini).counter = 0;
#line 251
	(bars->sl_prini).cycle = 0;
#line 251
}
   {
#line 252
	unsigned long	Error;
#line 252

#line 252
	Error = pthread_mutex_init(&(bars->sl_psini).mutex, NULL);
#line 252
	if (Error != 0) {
#line 252
		printf("Error while initializing barrier.\n");
#line 252
		exit(-1);
#line 252
	}
#line 252

#line 252
	Error = pthread_cond_init(&(bars->sl_psini).cv, NULL);
#line 252
	if (Error != 0) {
#line 252
		printf("Error while initializing barrier.\n");
#line 252
		pthread_mutex_destroy(&(bars->sl_psini).mutex);
#line 252
		exit(-1);
#line 252
	}
#line 252

#line 252
	(bars->sl_psini).counter = 0;
#line 252
	(bars->sl_psini).cycle = 0;
#line 252
}
   {
#line 253
	unsigned long	Error;
#line 253

#line 253
	Error = pthread_mutex_init(&(bars->sl_onetime).mutex, NULL);
#line 253
	if (Error != 0) {
#line 253
		printf("Error while initializing barrier.\n");
#line 253
		exit(-1);
#line 253
	}
#line 253

#line 253
	Error = pthread_cond_init(&(bars->sl_onetime).cv, NULL);
#line 253
	if (Error != 0) {
#line 253
		printf("Error while initializing barrier.\n");
#line 253
		pthread_mutex_destroy(&(bars->sl_onetime).mutex);
#line 253
		exit(-1);
#line 253
	}
#line 253

#line 253
	(bars->sl_onetime).counter = 0;
#line 253
	(bars->sl_onetime).cycle = 0;
#line 253
}
   {
#line 254
	unsigned long	Error;
#line 254

#line 254
	Error = pthread_mutex_init(&(bars->sl_phase_1).mutex, NULL);
#line 254
	if (Error != 0) {
#line 254
		printf("Error while initializing barrier.\n");
#line 254
		exit(-1);
#line 254
	}
#line 254

#line 254
	Error = pthread_cond_init(&(bars->sl_phase_1).cv, NULL);
#line 254
	if (Error != 0) {
#line 254
		printf("Error while initializing barrier.\n");
#line 254
		pthread_mutex_destroy(&(bars->sl_phase_1).mutex);
#line 254
		exit(-1);
#line 254
	}
#line 254

#line 254
	(bars->sl_phase_1).counter = 0;
#line 254
	(bars->sl_phase_1).cycle = 0;
#line 254
}
   {
#line 255
	unsigned long	Error;
#line 255

#line 255
	Error = pthread_mutex_init(&(bars->sl_phase_2).mutex, NULL);
#line 255
	if (Error != 0) {
#line 255
		printf("Error while initializing barrier.\n");
#line 255
		exit(-1);
#line 255
	}
#line 255

#line 255
	Error = pthread_cond_init(&(bars->sl_phase_2).cv, NULL);
#line 255
	if (Error != 0) {
#line 255
		printf("Error while initializing barrier.\n");
#line 255
		pthread_mutex_destroy(&(bars->sl_phase_2).mutex);
#line 255
		exit(-1);
#line 255
	}
#line 255

#line 255
	(bars->sl_phase_2).counter = 0;
#line 255
	(bars->sl_phase_2).cycle = 0;
#line 255
}
   {
#line 256
	unsigned long	Error;
#line 256

#line 256
	Error = pthread_mutex_init(&(bars->sl_phase_3).mutex, NULL);
#line 256
	if (Error != 0) {
#line 256
		printf("Error while initializing barrier.\n");
#line 256
		exit(-1);
#line 256
	}
#line 256

#line 256
	Error = pthread_cond_init(&(bars->sl_phase_3).cv, NULL);
#line 256
	if (Error != 0) {
#line 256
		printf("Error while initializing barrier.\n");
#line 256
		pthread_mutex_destroy(&(bars->sl_phase_3).mutex);
#line 256
		exit(-1);
#line 256
	}
#line 256

#line 256
	(bars->sl_phase_3).counter = 0;
#line 256
	(bars->sl_phase_3).cycle = 0;
#line 256
}
   {
#line 257
	unsigned long	Error;
#line 257

#line 257
	Error = pthread_mutex_init(&(bars->sl_phase_4).mutex, NULL);
#line 257
	if (Error != 0) {
#line 257
		printf("Error while initializing barrier.\n");
#line 257
		exit(-1);
#line 257
	}
#line 257

#line 257
	Error = pthread_cond_init(&(bars->sl_phase_4).cv, NULL);
#line 257
	if (Error != 0) {
#line 257
		printf("Error while initializing barrier.\n");
#line 257
		pthread_mutex_destroy(&(bars->sl_phase_4).mutex);
#line 257
		exit(-1);
#line 257
	}
#line 257

#line 257
	(bars->sl_phase_4).counter = 0;
#line 257
	(bars->sl_phase_4).cycle = 0;
#line 257
}
   {
#line 258
	unsigned long	Error;
#line 258

#line 258
	Error = pthread_mutex_init(&(bars->sl_phase_5).mutex, NULL);
#line 258
	if (Error != 0) {
#line 258
		printf("Error while initializing barrier.\n");
#line 258
		exit(-1);
#line 258
	}
#line 258

#line 258
	Error = pthread_cond_init(&(bars->sl_phase_5).cv, NULL);
#line 258
	if (Error != 0) {
#line 258
		printf("Error while initializing barrier.\n");
#line 258
		pthread_mutex_destroy(&(bars->sl_phase_5).mutex);
#line 258
		exit(-1);
#line 258
	}
#line 258

#line 258
	(bars->sl_phase_5).counter = 0;
#line 258
	(bars->sl_phase_5).cycle = 0;
#line 258
}
   {
#line 259
	unsigned long	Error;
#line 259

#line 259
	Error = pthread_mutex_init(&(bars->sl_phase_6).mutex, NULL);
#line 259
	if (Error != 0) {
#line 259
		printf("Error while initializing barrier.\n");
#line 259
		exit(-1);
#line 259
	}
#line 259

#line 259
	Error = pthread_cond_init(&(bars->sl_phase_6).cv, NULL);
#line 259
	if (Error != 0) {
#line 259
		printf("Error while initializing barrier.\n");
#line 259
		pthread_mutex_destroy(&(bars->sl_phase_6).mutex);
#line 259
		exit(-1);
#line 259
	}
#line 259

#line 259
	(bars->sl_phase_6).counter = 0;
#line 259
	(bars->sl_phase_6).cycle = 0;
#line 259
}
   {
#line 260
	unsigned long	Error;
#line 260

#line 260
	Error = pthread_mutex_init(&(bars->sl_phase_7).mutex, NULL);
#line 260
	if (Error != 0) {
#line 260
		printf("Error while initializing barrier.\n");
#line 260
		exit(-1);
#line 260
	}
#line 260

#line 260
	Error = pthread_cond_init(&(bars->sl_phase_7).cv, NULL);
#line 260
	if (Error != 0) {
#line 260
		printf("Error while initializing barrier.\n");
#line 260
		pthread_mutex_destroy(&(bars->sl_phase_7).mutex);
#line 260
		exit(-1);
#line 260
	}
#line 260

#line 260
	(bars->sl_phase_7).counter = 0;
#line 260
	(bars->sl_phase_7).cycle = 0;
#line 260
}
   {
#line 261
	unsigned long	Error;
#line 261

#line 261
	Error = pthread_mutex_init(&(bars->sl_phase_8).mutex, NULL);
#line 261
	if (Error != 0) {
#line 261
		printf("Error while initializing barrier.\n");
#line 261
		exit(-1);
#line 261
	}
#line 261

#line 261
	Error = pthread_cond_init(&(bars->sl_phase_8).cv, NULL);
#line 261
	if (Error != 0) {
#line 261
		printf("Error while initializing barrier.\n");
#line 261
		pthread_mutex_destroy(&(bars->sl_phase_8).mutex);
#line 261
		exit(-1);
#line 261
	}
#line 261

#line 261
	(bars->sl_phase_8).counter = 0;
#line 261
	(bars->sl_phase_8).cycle = 0;
#line 261
}
   {
#line 262
	unsigned long	Error;
#line 262

#line 262
	Error = pthread_mutex_init(&(bars->sl_phase_9).mutex, NULL);
#line 262
	if (Error != 0) {
#line 262
		printf("Error while initializing barrier.\n");
#line 262
		exit(-1);
#line 262
	}
#line 262

#line 262
	Error = pthread_cond_init(&(bars->sl_phase_9).cv, NULL);
#line 262
	if (Error != 0) {
#line 262
		printf("Error while initializing barrier.\n");
#line 262
		pthread_mutex_destroy(&(bars->sl_phase_9).mutex);
#line 262
		exit(-1);
#line 262
	}
#line 262

#line 262
	(bars->sl_phase_9).counter = 0;
#line 262
	(bars->sl_phase_9).cycle = 0;
#line 262
}
   {
#line 263
	unsigned long	Error;
#line 263

#line 263
	Error = pthread_mutex_init(&(bars->sl_phase_10).mutex, NULL);
#line 263
	if (Error != 0) {
#line 263
		printf("Error while initializing barrier.\n");
#line 263
		exit(-1);
#line 263
	}
#line 263

#line 263
	Error = pthread_cond_init(&(bars->sl_phase_10).cv, NULL);
#line 263
	if (Error != 0) {
#line 263
		printf("Error while initializing barrier.\n");
#line 263
		pthread_mutex_destroy(&(bars->sl_phase_10).mutex);
#line 263
		exit(-1);
#line 263
	}
#line 263

#line 263
	(bars->sl_phase_10).counter = 0;
#line 263
	(bars->sl_phase_10).cycle = 0;
#line 263
}
   {
#line 264
	unsigned long	Error;
#line 264

#line 264
	Error = pthread_mutex_init(&(bars->error_barrier).mutex, NULL);
#line 264
	if (Error != 0) {
#line 264
		printf("Error while initializing barrier.\n");
#line 264
		exit(-1);
#line 264
	}
#line 264

#line 264
	Error = pthread_cond_init(&(bars->error_barrier).cv, NULL);
#line 264
	if (Error != 0) {
#line 264
		printf("Error while initializing barrier.\n");
#line 264
		pthread_mutex_destroy(&(bars->error_barrier).mutex);
#line 264
		exit(-1);
#line 264
	}
#line 264

#line 264
	(bars->error_barrier).counter = 0;
#line 264
	(bars->error_barrier).cycle = 0;
#line 264
}
#else
   {
#line 266
	unsigned long	Error;
#line 266

#line 266
	Error = pthread_mutex_init(&(bars->barrier).mutex, NULL);
#line 266
	if (Error != 0) {
#line 266
		printf("Error while initializing barrier.\n");
#line 266
		exit(-1);
#line 266
	}
#line 266

#line 266
	Error = pthread_cond_init(&(bars->barrier).cv, NULL);
#line 266
	if (Error != 0) {
#line 266
		printf("Error while initializing barrier.\n");
#line 266
		pthread_mutex_destroy(&(bars->barrier).mutex);
#line 266
		exit(-1);
#line 266
	}
#line 266

#line 266
	(bars->barrier).counter = 0;
#line 266
	(bars->barrier).cycle = 0;
#line 266
}
#endif

   imx[numlev-1] = im;
   jmx[numlev-1] = jm;
   lev_res[numlev-1] = res;
   lev_tol[numlev-1] = tolerance;
   multi->err_multi = 0.0;
   multi->numspin = 0;
   for (i=0;i<nprocs;i++) {
     multi->spinflag[i] = 0;
   }

   for (i=numlev-2;i>=0;i--) {
     imx[i] = ((imx[i+1] - 2) / 2) + 2;
     jmx[i] = ((jmx[i+1] - 2) / 2) + 2;
     lev_res[i] = lev_res[i+1] * 2;
   }

   xprocs = 0;
   yprocs = 0;
   procsqrt = sqrt((double) nprocs);
   j = (long) procsqrt;
   while ((xprocs == 0) && (j > 0)) {
     k = nprocs / j;
     if (k * j == nprocs) {
       if (k > j) {
         xprocs = j;
         yprocs = k;
       } else {
         xprocs = k;
         yprocs = j;
       }
     }
     j--;
   }
   if (xprocs == 0) {
     printerr("Could not find factors for subblocking\n");
     exit(-1);
   }

/* Determine starting coord and number of points to process in     */
/* each direction                                                  */

   for (i=0;i<numlev;i++) {
     xportion = (jmx[i] - 2) / xprocs;
     xextra = (jmx[i] - 2) % xprocs;
     for (j=0;j<xprocs;j++) {
       if (xextra == 0) {
         for (k=0;k<yprocs;k++) {
           gp[k*xprocs+j].rel_start_x[i] = j * xportion + 1;
           gp[k*xprocs+j].rel_num_x[i] = xportion;
         }
       } else {
         if (j + 1 > xextra) {
           for (k=0;k<yprocs;k++) {
             lower = xextra * (xportion + 1);
             gp[k*xprocs+j].rel_start_x[i] = lower + (j - xextra) * xportion + 1;
             gp[k*xprocs+j].rel_num_x[i] = xportion;
           }
         } else {
           for (k=0;k<yprocs;k++) {
             gp[k*xprocs+j].rel_start_x[i] = j * (xportion + 1) + 1;
             gp[k*xprocs+j].rel_num_x[i] = xportion + 1;
           }
         }
       }
     }
     yportion = (imx[i] - 2) / yprocs;
     yextra = (imx[i] - 2) % yprocs;
     for (j=0;j<yprocs;j++) {
       if (yextra == 0) {
         for (k=0;k<xprocs;k++) {
           gp[j*xprocs+k].rel_start_y[i] = j * yportion + 1;
           gp[j*xprocs+k].rel_num_y[i] = yportion;
         }
       } else {
         if (j + 1 > yextra) {
           for (k=0;k<xprocs;k++) {
             lower = yextra * (yportion + 1);
             gp[j*xprocs+k].rel_start_y[i] = lower + (j - yextra) * yportion + 1;
             gp[j*xprocs+k].rel_num_y[i] = yportion;
           }
         } else {
           for (k=0;k<xprocs;k++) {
             gp[j*xprocs+k].rel_start_y[i] = j * (yportion + 1) + 1;
             gp[j*xprocs+k].rel_num_y[i] = yportion + 1;
           }
         }
       }
     }
   }

   i_int_coeff[0] = 0.0;
   j_int_coeff[0] = 0.0;
   for (i=0;i<numlev;i++) {
     i_int_coeff[i] = 1.0/(imx[i]-1);
     j_int_coeff[i] = 1.0/(jmx[i]-1);
   }

   for (my_num=0;my_num<nprocs;my_num++) {
     for (i=0;i<numlev;i++) {
       gp[my_num].rlist[i] = gp[my_num].rel_start_y[i];
       gp[my_num].rljst[i] = gp[my_num].rel_start_x[i];
       gp[my_num].rlien[i] = gp[my_num].rlist[i] + gp[my_num].rel_num_y[i] - 1;
       gp[my_num].rljen[i] = gp[my_num].rljst[i] + gp[my_num].rel_num_x[i] - 1;
       gp[my_num].iist[i] = gp[my_num].rel_start_y[i];
       gp[my_num].ijst[i] = gp[my_num].rel_start_x[i];
       gp[my_num].iien[i] = gp[my_num].iist[i] + gp[my_num].rel_num_y[i] - 1;
       gp[my_num].ijen[i] = gp[my_num].ijst[i] + gp[my_num].rel_num_x[i] - 1;
       gp[my_num].pist[i] = gp[my_num].rel_start_y[i];
       gp[my_num].pjst[i] = gp[my_num].rel_start_x[i];
       gp[my_num].pien[i] = gp[my_num].pist[i] + gp[my_num].rel_num_y[i] - 1;
       gp[my_num].pjen[i] = gp[my_num].pjst[i] + gp[my_num].rel_num_x[i] - 1;

       if (gp[my_num].pist[i] == 1) {
         gp[my_num].pist[i] = 0;
       }
       if (gp[my_num].pjst[i] == 1) {
         gp[my_num].pjst[i] = 0;
       }
       if (gp[my_num].pien[i] == imx[i] - 2) {
         gp[my_num].pien[i] = imx[i]-1;
       }
       if (gp[my_num].pjen[i] == jmx[i] - 2) {
         gp[my_num].pjen[i] = jmx[i]-1;
       }

       if (gp[my_num].rlist[i] % 2 == 0) {
         gp[my_num].eist[i] = gp[my_num].rlist[i];
         gp[my_num].oist[i] = gp[my_num].rlist[i] + 1;
       } else {
         gp[my_num].eist[i] = gp[my_num].rlist[i] + 1;
         gp[my_num].oist[i] = gp[my_num].rlist[i];
       }
       if (gp[my_num].rljst[i] % 2 == 0) {
         gp[my_num].ejst[i] = gp[my_num].rljst[i];
         gp[my_num].ojst[i] = gp[my_num].rljst[i] + 1;
       } else {
         gp[my_num].ejst[i] = gp[my_num].rljst[i] + 1;
         gp[my_num].ojst[i] = gp[my_num].rljst[i];
       }
       if (gp[my_num].rlien[i] == imx[i]-2) {
         gp[my_num].rlien[i] = gp[my_num].rlien[i] - 1;
         if (gp[my_num].rlien[i] % 2 == 0) {
           gp[my_num].ojest[i] = gp[my_num].ojst[i];
           gp[my_num].ejest[i] = gp[my_num].ejst[i];
         } else {
           gp[my_num].ojest[i] = gp[my_num].ejst[i];
           gp[my_num].ejest[i] = gp[my_num].ojst[i];
         }
       }
       if (gp[my_num].rljen[i] == jmx[i]-2) {
         gp[my_num].rljen[i] = gp[my_num].rljen[i] - 1;
         if (gp[my_num].rljen[i] % 2 == 0) {
           gp[my_num].oiest[i] = gp[my_num].oist[i];
           gp[my_num].eiest[i] = gp[my_num].eist[i];
         } else {
           gp[my_num].oiest[i] = gp[my_num].eist[i];
           gp[my_num].eiest[i] = gp[my_num].oist[i];
         }
       }
     }
   }

/* initialize constants and variables

   id is a global shared variable that has fetch-and-add operations
   performed on it by processes to obtain their pids.   */

   global->id = 0;
   global->psibi = 0.0;
   pi = atan(1.0);
   pi = 4.*pi;

   factjacob = -1./(12.*res*res);
   factlap = 1./(res*res);
   eig2 = -h*f0*f0/(h1*h3*gpr);
   jmm1 = jm-1 ;
   ysca = ((double) jmm1)*res ;
   for (i=0;i<im;i++) {
     for (j=0;j<jm;j++) {
       guess->oldga[i][j] = 0.0;
       guess->oldgb[i][j] = 0.0;
     }
   }

   if (do_output) {
     printf("                       MULTIGRID OUTPUTS\n");
   }

#ifdef ENABLE_PARSEC_HOOKS
	__parsec_roi_begin();
#endif
   {
#line 460
	long	i, Error;
#line 460

#line 460
	for (i = 0; i < (nprocs) - 1; i++) {
#line 460
		Error = pthread_create(&PThreadTable[i], NULL, (void * (*)(void *))(slave), NULL);
#line 460
		if (Error != 0) {
#line 460
			printf("Error in pthread_create().\n");
#line 460
			exit(-1);
#line 460
		}
#line 460
	}
#line 460

#line 460
	slave();
#line 460
};
   {
#line 461
	unsigned long	i, Error;
#line 461
	for (i = 0; i < (nprocs) - 1; i++) {
#line 461
		Error = pthread_join(PThreadTable[i], NULL);
#line 461
		if (Error != 0) {
#line 461
			printf("Error in pthread_join().\n");
#line 461
			exit(-1);
#line 461
		}
#line 461
	}
#line 461
};
#ifdef ENABLE_PARSEC_HOOKS
	__parsec_roi_end();
#endif
   {
#line 465
	struct timeval	FullTime;
#line 465

#line 465
	gettimeofday(&FullTime, NULL);
#line 465
	(computeend) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 465
}

   printf("\n");
   printf("                       PROCESS STATISTICS\n");
   printf("                  Total          Multigrid         Multigrid\n");
   printf(" Proc             Time             Time            Fraction\n");
   printf("    0   %15.0f    %15.0f        %10.3f\n", gp[0].total_time,gp[0].multi_time, gp[0].multi_time/gp[0].total_time);

   if (do_stats) {
     min_total = max_total = avg_total = gp[0].total_time;
     min_multi = max_multi = avg_multi = gp[0].multi_time;
     min_frac = max_frac = avg_frac = gp[0].multi_time/gp[0].total_time;
     for (i=1;i<nprocs;i++) {
       if (gp[i].total_time > max_total) {
         max_total = gp[i].total_time;
       }
       if (gp[i].total_time < min_total) {
         min_total = gp[i].total_time;
       }
       if (gp[i].multi_time > max_multi) {
         max_multi = gp[i].multi_time;
       }
       if (gp[i].multi_time < min_multi) {
         min_multi = gp[i].multi_time;
       }
       if (gp[i].multi_time/gp[i].total_time > max_frac) {
         max_frac = gp[i].multi_time/gp[i].total_time;
       }
       if (gp[i].multi_time/gp[i].total_time < min_frac) {
         min_frac = gp[i].multi_time/gp[i].total_time;
       }
       avg_total += gp[i].total_time;
       avg_multi += gp[i].multi_time;
       avg_frac += gp[i].multi_time/gp[i].total_time;
     }
     avg_total = avg_total / nprocs;
     avg_multi = avg_multi / nprocs;
     avg_frac = avg_frac / nprocs;
     for (i=1;i<nprocs;i++) {
       printf("  %3ld   %15.0f    %15.0f        %10.3f\n", i, gp[i].total_time, gp[i].multi_time, gp[i].multi_time/gp[i].total_time);
     }
     printf("  Avg   %15.0f    %15.0f        %10.3f\n", avg_total,avg_multi,avg_frac);
     printf("  Min   %15.0f    %15.0f        %10.3f\n", min_total,min_multi,min_frac);
     printf("  Max   %15.0f    %15.0f        %10.3f\n", max_total,max_multi,max_frac);
   }
   printf("\n");

   global->starttime = start;
   printf("                       TIMING INFORMATION\n");
   printf("Start time                        : %16lu\n", global->starttime);
   printf("Initialization finish time        : %16lu\n", global->trackstart);
   printf("Overall finish time               : %16lu\n", computeend);
   printf("Total time with initialization    : %16lu\n", computeend-global->starttime);
   printf("Total time without initialization : %16lu\n", computeend-global->trackstart);
   printf("    (excludes first timestep)\n");
   printf("\n");
   {exit(0);}
#ifdef ENABLE_PARSEC_HOOKS
	__parsec_bench_end();
#endif
}

long log_2(long number)
{
  long cumulative = 1;
  long out = 0;
  long done = 0;

  while ((cumulative < number) && (!done) && (out < 50)) {
    if (cumulative == number) {
      done = 1;
    } else {
      cumulative = cumulative * 2;
      out ++;
    }
  }

  if (cumulative == number) {
    return(out);
  } else {
    return(-1);
  }
}

void printerr(char *s)
{
  fprintf(stderr,"ERROR: %s\n",s);
}
