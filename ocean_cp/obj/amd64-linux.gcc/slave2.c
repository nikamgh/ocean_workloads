#line 233 "/mnt/panzer/mnika/Workload_Pool/parsec-3.0/bin/../pkgs/libs/parmacs/inst/amd64-linux.gcc/m4/parmacs.pthreads.c.m4"

#line 1 "slave2.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/

/*    ****************
      subroutine slave2
      ****************  */


#line 21
#include <pthread.h>
#line 21
#include <sys/time.h>
#line 21
#include <unistd.h>
#line 21
#include <stdlib.h>
#line 21
extern pthread_t PThreadTable[];
#line 21


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "decs.h"
#include "/mnt/panzer/mnika/Workload_Pool/zsim_hooks.h"

void slave2(long procid, long firstrow, long lastrow, long numrows, long firstcol, long lastcol, long numcols)
{

//zsim_PIM_function_begin();
   long i;
   long j;
   long iindex;
   double hh1;
   double hh3;
   double hinv;
   double h1inv;
   long istart;
   long iend;
   long jstart;
   long jend;
   long ist;
   long ien;
   long jst;
   long jen;
   double fac;
   double ressqr;
   double psiaipriv;
   double f4;
   double timst;
   long psiindex;
   long i_off;
   long j_off;
   long multi_start;
   long multi_end;
   double **t2a;
   double **t2b;
   double **t2c;
   double **t2d;
   double **t2e;
   double **t2f;
   double **t2g;
   double **t2h;
   double *t1a;
   double *t1b;
   double *t1c;
   double *t1d;
   double *t1e;
   double *t1f;
   double *t1g;
   double *t1h;

   ressqr = lev_res[numlev-1] * lev_res[numlev-1];
   i_off = gp[procid].rownum*numrows;
   j_off = gp[procid].colnum*numcols;
   zsim_stamp();
/*   ***************************************************************

          f i r s t     p h a s e   (of timestep calculation)

     ***************************************************************/

   t2a = (double **) ga[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0]=0.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0]=0.0;
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1]=0.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1]=0.0;
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = 0.0;
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = 0.0;
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = 0.0;
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = 0.0;
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       t1a[iindex] = 0.0;
     }
   }

   t2a = (double **) gb[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0]=0.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0]=0.0;
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1]=0.0;
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1]=0.0;
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = 0.0;
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = 0.0;
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = 0.0;
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = 0.0;
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       t1a[iindex] = 0.0;
     }
   }

/* put the laplacian of psi{1,3} in work1{1,2}
   note that psi(i,j,2) represents the psi3 array in
   the original equations  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     t2a = (double **) work1[procid][psiindex];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = 0;
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = 0;
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = 0;
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1] = 0;
     }
     laplacalc(procid,psi,work1,psiindex,
	       firstrow,lastrow,firstcol,lastcol);
   }

/* set values of work2 array to psi1 - psi3   */

   t2a = (double **) work2[procid];
   t2b = (double **) psi[procid][0];
   t2c = (double **) psi[procid][1];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0] = t2b[0][0]-t2c[0][0];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0] = t2b[im-1][0]-t2c[im-1][0];
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1] = t2b[0][jm-1]-t2c[0][jm-1];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1] = t2b[im-1][jm-1] -
				 t2c[im-1][jm-1];
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) t2b[0];
     t1c = (double *) t2c[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1b[j]-t1c[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) t2b[im-1];
     t1c = (double *) t2c[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1b[j]-t1c[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = t2b[j][0]-t2c[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = t2b[j][jm-1]-t2c[j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     t1c = (double *) t2c[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = t1b[iindex] - t1c[iindex];
     }
   }

/* set values of work3 array to h3/h * psi1 + h1/h * psi3  */

   t2a = (double **) work3[procid];
   hh3 = h3/h;
   hh1 = h1/h;
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0] = hh3*t2a[0][0]+hh1*t2c[0][0];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0] = hh3*t2a[im-1][0] +
			      hh1*t2c[im-1][0];
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1] = hh3*t2a[0][jm-1] +
			      hh1*t2c[0][jm-1];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1] = hh3*t2a[im-1][jm-1] +
				 hh1*t2c[im-1][jm-1];
   }
   if (gp[procid].neighbors[UP] == -1) {
     for(j=firstcol;j<=lastcol;j++) {
       t2a[0][j] = hh3*t2a[0][j]+hh1*t2c[0][j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     for(j=firstcol;j<=lastcol;j++) {
       t2a[im-1][j] = hh3*t2a[im-1][j] +
				hh1*t2c[im-1][j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = hh3*t2a[j][0]+hh1*t2c[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = hh3*t2a[j][jm-1] +
				hh1*t2c[j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     t1c = (double *) t2c[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
        t1a[iindex] = hh3*t1a[iindex] + hh1*t1c[iindex];
     }
   }

/* set values of temparray{1,3} to psim{1,3}  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     t2a = (double **) temparray[procid][psiindex];
     t2b = (double **) psi[procid][psiindex];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = t2b[0][0];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = t2b[im-1][0];
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = t2b[0][jm-1];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1] = t2b[im-1][jm-1];
     }
     if (gp[procid].neighbors[UP] == -1) {
       for(j=firstcol;j<=lastcol;j++) {
         t2a[0][j] = t2b[0][j];
       }
     }
     if (gp[procid].neighbors[DOWN] == -1) {
       for(j=firstcol;j<=lastcol;j++) {
         t2a[im-1][j] = t2b[im-1][j];
       }
     }
     if (gp[procid].neighbors[LEFT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][0] = t2b[j][0];
       }
     }
     if (gp[procid].neighbors[RIGHT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][jm-1] = t2b[j][jm-1];
       }
     }

     for(i=firstrow;i<=lastrow;i++) {
       t1a = (double *) t2a[i];
       t1b = (double *) t2b[i];
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = t1b[iindex];
       }
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 342
	unsigned long	Error, Cycle;
#line 342
	int		Cancel, Temp;
#line 342

#line 342
	Error = pthread_mutex_lock(&(bars->sl_phase_1).mutex);
#line 342
	if (Error != 0) {
#line 342
		printf("Error while trying to get lock in barrier.\n");
#line 342
		exit(-1);
#line 342
	}
#line 342

#line 342
	Cycle = (bars->sl_phase_1).cycle;
#line 342
	if (++(bars->sl_phase_1).counter != (nprocs)) {
#line 342
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 342
		while (Cycle == (bars->sl_phase_1).cycle) {
#line 342
			Error = pthread_cond_wait(&(bars->sl_phase_1).cv, &(bars->sl_phase_1).mutex);
#line 342
			if (Error != 0) {
#line 342
				break;
#line 342
			}
#line 342
		}
#line 342
		pthread_setcancelstate(Cancel, &Temp);
#line 342
	} else {
#line 342
		(bars->sl_phase_1).cycle = !(bars->sl_phase_1).cycle;
#line 342
		(bars->sl_phase_1).counter = 0;
#line 342
		Error = pthread_cond_broadcast(&(bars->sl_phase_1).cv);
#line 342
	}
#line 342
	pthread_mutex_unlock(&(bars->sl_phase_1).mutex);
#line 342
}
#else
   {
#line 344
	unsigned long	Error, Cycle;
#line 344
	int		Cancel, Temp;
#line 344

#line 344
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 344
	if (Error != 0) {
#line 344
		printf("Error while trying to get lock in barrier.\n");
#line 344
		exit(-1);
#line 344
	}
#line 344

#line 344
	Cycle = (bars->barrier).cycle;
#line 344
	if (++(bars->barrier).counter != (nprocs)) {
#line 344
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 344
		while (Cycle == (bars->barrier).cycle) {
#line 344
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 344
			if (Error != 0) {
#line 344
				break;
#line 344
			}
#line 344
		}
#line 344
		pthread_setcancelstate(Cancel, &Temp);
#line 344
	} else {
#line 344
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 344
		(bars->barrier).counter = 0;
#line 344
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 344
	}
#line 344
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 344
}
#endif
   zsim_stamp();
/*     *******************************************************

              s e c o n d   p h a s e

       *******************************************************

   set values of psi{1,3} to psim{1,3}   */

   for(psiindex=0;psiindex<=1;psiindex++) {
     t2a = (double **) psi[procid][psiindex];
     t2b = (double **) psim[procid][psiindex];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = t2b[0][0];
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = t2b[0][jm-1];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = t2b[im-1][0];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1] = t2b[im-1][jm-1];
     }
     if (gp[procid].neighbors[UP] == -1) {
       for(j=firstcol;j<=lastcol;j++) {
         t2a[0][j] = t2b[0][j];
       }
     }
     if (gp[procid].neighbors[DOWN] == -1) {
       for(j=firstcol;j<=lastcol;j++) {
         t2a[im-1][j] = t2b[im-1][j];
       }
     }
     if (gp[procid].neighbors[LEFT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][0] = t2b[j][0];
       }
     }
     if (gp[procid].neighbors[RIGHT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][jm-1] = t2b[j][jm-1];
       }
     }

     for(i=firstrow;i<=lastrow;i++) {
       t1a = (double *) t2a[i];
       t1b = (double *) t2b[i];
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = t1b[iindex];
       }
     }
   }

/* put the laplacian of the psim array
   into the work7 array; first part of a three-laplacian
   calculation to compute the friction terms  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     t2a = (double **) work7[procid][psiindex];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = 0;
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = 0;
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = 0;
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1] = 0;
     }
     laplacalc(procid,psim,work7,psiindex,
	       firstrow,lastrow,firstcol,lastcol);
   }

/* to the values of the work1{1,2} arrays obtained from the
   laplacians of psi{1,2} in the previous phase, add to the
   elements of every column the corresponding value in the
   one-dimenional f array  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     t2a = (double **) work1[procid][psiindex];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = t2a[0][0] + f[0];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = t2a[im-1][0] + f[0];
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = t2a[0][jm-1] + f[jmx[numlev-1]-1];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1]=t2a[im-1][jm-1] + f[jmx[numlev-1]-1];
     }
     if (gp[procid].neighbors[UP] == -1) {
       for(j=firstcol;j<=lastcol;j++) {
         t2a[0][j] = t2a[0][j] + f[j+j_off];
       }
     }
     if (gp[procid].neighbors[DOWN] == -1) {
       for(j=firstcol;j<=lastcol;j++) {
         t2a[im-1][j] = t2a[im-1][j] + f[j+j_off];
       }
     }
     if (gp[procid].neighbors[LEFT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][0] = t2a[j][0] + f[j+i_off];
       }
     }
     if (gp[procid].neighbors[RIGHT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][jm-1] = t2a[j][jm-1] + f[j+i_off];
       }
     }
     for(i=firstrow;i<=lastrow;i++) {
       t1a = (double *) t2a[i];
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex]=t1a[iindex] + f[iindex+j_off];
       }
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 469
	unsigned long	Error, Cycle;
#line 469
	int		Cancel, Temp;
#line 469

#line 469
	Error = pthread_mutex_lock(&(bars->sl_phase_2).mutex);
#line 469
	if (Error != 0) {
#line 469
		printf("Error while trying to get lock in barrier.\n");
#line 469
		exit(-1);
#line 469
	}
#line 469

#line 469
	Cycle = (bars->sl_phase_2).cycle;
#line 469
	if (++(bars->sl_phase_2).counter != (nprocs)) {
#line 469
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 469
		while (Cycle == (bars->sl_phase_2).cycle) {
#line 469
			Error = pthread_cond_wait(&(bars->sl_phase_2).cv, &(bars->sl_phase_2).mutex);
#line 469
			if (Error != 0) {
#line 469
				break;
#line 469
			}
#line 469
		}
#line 469
		pthread_setcancelstate(Cancel, &Temp);
#line 469
	} else {
#line 469
		(bars->sl_phase_2).cycle = !(bars->sl_phase_2).cycle;
#line 469
		(bars->sl_phase_2).counter = 0;
#line 469
		Error = pthread_cond_broadcast(&(bars->sl_phase_2).cv);
#line 469
	}
#line 469
	pthread_mutex_unlock(&(bars->sl_phase_2).mutex);
#line 469
}
#else
   {
#line 471
	unsigned long	Error, Cycle;
#line 471
	int		Cancel, Temp;
#line 471

#line 471
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 471
	if (Error != 0) {
#line 471
		printf("Error while trying to get lock in barrier.\n");
#line 471
		exit(-1);
#line 471
	}
#line 471

#line 471
	Cycle = (bars->barrier).cycle;
#line 471
	if (++(bars->barrier).counter != (nprocs)) {
#line 471
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 471
		while (Cycle == (bars->barrier).cycle) {
#line 471
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 471
			if (Error != 0) {
#line 471
				break;
#line 471
			}
#line 471
		}
#line 471
		pthread_setcancelstate(Cancel, &Temp);
#line 471
	} else {
#line 471
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 471
		(bars->barrier).counter = 0;
#line 471
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 471
	}
#line 471
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 471
}
#endif

   zsim_stamp();
/* 	*******************************************************

                 t h i r d   p h a s e

 	*******************************************************

   put the jacobian of the work1{1,2} and psi{1,3} arrays
   (the latter currently in temparray) in the work5{1,2} arrays  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     jacobcalc2(work1,temparray,work5,psiindex,procid,firstrow,lastrow,
	       firstcol,lastcol);
   }

/* set values of psim{1,3} to temparray{1,3}  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     t2a = (double **) psim[procid][psiindex];
     t2b = (double **) temparray[procid][psiindex];
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[0][0] = t2b[0][0];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
       t2a[im-1][0] = t2b[im-1][0];
     }
     if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[0][jm-1] = t2b[0][jm-1];
     }
     if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
       t2a[im-1][jm-1] = t2b[im-1][jm-1];
     }
     if (gp[procid].neighbors[UP] == -1) {
       t1a = (double *) t2a[0];
       t1b = (double *) t2b[0];
       for(j=firstcol;j<=lastcol;j++) {
         t1a[j] = t1b[j];
       }
     }
     if (gp[procid].neighbors[DOWN] == -1) {
       t1a = (double *) t2a[im-1];
       t1b = (double *) t2b[im-1];
       for(j=firstcol;j<=lastcol;j++) {
         t1a[j] = t1b[j];
       }
     }
     if (gp[procid].neighbors[LEFT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][0] = t2b[j][0];
       }
     }
     if (gp[procid].neighbors[RIGHT] == -1) {
       for(j=firstrow;j<=lastrow;j++) {
         t2a[j][jm-1] = t2b[j][jm-1];
       }
     }
     for(i=firstrow;i<=lastrow;i++) {
       t1a = (double *) t2a[i];
       t1b = (double *) t2b[i];
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = t1b[iindex];
       }
     }
   }

/* put the laplacian of the work7{1,2} arrays in the work4{1,2}
   arrays; second step in the three-laplacian friction calculation  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     laplacalc(procid,work7,work4,psiindex,
	       firstrow,lastrow,firstcol,lastcol);
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 547
	unsigned long	Error, Cycle;
#line 547
	int		Cancel, Temp;
#line 547

#line 547
	Error = pthread_mutex_lock(&(bars->sl_phase_3).mutex);
#line 547
	if (Error != 0) {
#line 547
		printf("Error while trying to get lock in barrier.\n");
#line 547
		exit(-1);
#line 547
	}
#line 547

#line 547
	Cycle = (bars->sl_phase_3).cycle;
#line 547
	if (++(bars->sl_phase_3).counter != (nprocs)) {
#line 547
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 547
		while (Cycle == (bars->sl_phase_3).cycle) {
#line 547
			Error = pthread_cond_wait(&(bars->sl_phase_3).cv, &(bars->sl_phase_3).mutex);
#line 547
			if (Error != 0) {
#line 547
				break;
#line 547
			}
#line 547
		}
#line 547
		pthread_setcancelstate(Cancel, &Temp);
#line 547
	} else {
#line 547
		(bars->sl_phase_3).cycle = !(bars->sl_phase_3).cycle;
#line 547
		(bars->sl_phase_3).counter = 0;
#line 547
		Error = pthread_cond_broadcast(&(bars->sl_phase_3).cv);
#line 547
	}
#line 547
	pthread_mutex_unlock(&(bars->sl_phase_3).mutex);
#line 547
}
#else
   {
#line 549
	unsigned long	Error, Cycle;
#line 549
	int		Cancel, Temp;
#line 549

#line 549
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 549
	if (Error != 0) {
#line 549
		printf("Error while trying to get lock in barrier.\n");
#line 549
		exit(-1);
#line 549
	}
#line 549

#line 549
	Cycle = (bars->barrier).cycle;
#line 549
	if (++(bars->barrier).counter != (nprocs)) {
#line 549
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 549
		while (Cycle == (bars->barrier).cycle) {
#line 549
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 549
			if (Error != 0) {
#line 549
				break;
#line 549
			}
#line 549
		}
#line 549
		pthread_setcancelstate(Cancel, &Temp);
#line 549
	} else {
#line 549
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 549
		(bars->barrier).counter = 0;
#line 549
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 549
	}
#line 549
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 549
}
#endif

   zsim_stamp();
/*     *******************************************************

                f o u r t h   p h a s e

       *******************************************************

   put the jacobian of the work2 and work3 arrays in the work6
   array  */

   jacobcalc(work2,work3,work6,procid,firstrow,lastrow,firstcol,lastcol);

/* put the laplacian of the work4{1,2} arrays in the work7{1,2}
   arrays; third step in the three-laplacian friction calculation  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     laplacalc(procid,work4,work7,psiindex,
	       firstrow,lastrow,firstcol,lastcol);
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 572
	unsigned long	Error, Cycle;
#line 572
	int		Cancel, Temp;
#line 572

#line 572
	Error = pthread_mutex_lock(&(bars->sl_phase_4).mutex);
#line 572
	if (Error != 0) {
#line 572
		printf("Error while trying to get lock in barrier.\n");
#line 572
		exit(-1);
#line 572
	}
#line 572

#line 572
	Cycle = (bars->sl_phase_4).cycle;
#line 572
	if (++(bars->sl_phase_4).counter != (nprocs)) {
#line 572
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 572
		while (Cycle == (bars->sl_phase_4).cycle) {
#line 572
			Error = pthread_cond_wait(&(bars->sl_phase_4).cv, &(bars->sl_phase_4).mutex);
#line 572
			if (Error != 0) {
#line 572
				break;
#line 572
			}
#line 572
		}
#line 572
		pthread_setcancelstate(Cancel, &Temp);
#line 572
	} else {
#line 572
		(bars->sl_phase_4).cycle = !(bars->sl_phase_4).cycle;
#line 572
		(bars->sl_phase_4).counter = 0;
#line 572
		Error = pthread_cond_broadcast(&(bars->sl_phase_4).cv);
#line 572
	}
#line 572
	pthread_mutex_unlock(&(bars->sl_phase_4).mutex);
#line 572
}
#else
   {
#line 574
	unsigned long	Error, Cycle;
#line 574
	int		Cancel, Temp;
#line 574

#line 574
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 574
	if (Error != 0) {
#line 574
		printf("Error while trying to get lock in barrier.\n");
#line 574
		exit(-1);
#line 574
	}
#line 574

#line 574
	Cycle = (bars->barrier).cycle;
#line 574
	if (++(bars->barrier).counter != (nprocs)) {
#line 574
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 574
		while (Cycle == (bars->barrier).cycle) {
#line 574
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 574
			if (Error != 0) {
#line 574
				break;
#line 574
			}
#line 574
		}
#line 574
		pthread_setcancelstate(Cancel, &Temp);
#line 574
	} else {
#line 574
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 574
		(bars->barrier).counter = 0;
#line 574
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 574
	}
#line 574
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 574
}
#endif
   zsim_stamp();
/*     *******************************************************

                f i f t h   p h a s e

       *******************************************************

   use the values of the work5, work6 and work7 arrays
   computed in the previous time-steps to compute the
   ga and gb arrays   */

   hinv = 1.0/h;
   h1inv = 1.0/h1;

   t2a = (double **) ga[procid];
   t2b = (double **) gb[procid];
   t2c = (double **) work5[procid][0];
   t2d = (double **) work5[procid][1];
   t2e = (double **) work7[procid][0];
   t2f = (double **) work7[procid][1];
   t2g = (double **) work6[procid];
   t2h = (double **) tauz[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0] = t2c[0][0]-t2d[0][0] +
			eig2*t2g[0][0]+h1inv*t2h[0][0] +
			lf*t2e[0][0]-lf*t2f[0][0];
     t2b[0][0] = hh1*t2c[0][0]+hh3*t2d[0][0] +
			hinv*t2h[0][0]+lf*hh1*t2e[0][0] +
		        lf*hh3*t2f[0][0];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0] = t2c[im-1][0]-t2d[im-1][0] +
	   eig2*t2g[im-1][0] + h1inv*t2h[im-1][0] +
	   lf*t2e[im-1][0] - lf*t2f[im-1][0];
     t2b[im-1][0] = hh1*t2c[im-1][0] +
	   hh3*t2d[im-1][0] + hinv*t2h[im-1][0] +
	   lf*hh1*t2e[im-1][0] + lf*hh3*t2f[im-1][0];
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1] = t2c[0][jm-1]-t2d[0][jm-1]+
	   eig2*t2g[0][jm-1]+h1inv*t2h[0][jm-1] +
	   lf*t2e[0][jm-1]-lf*t2f[0][jm-1];
     t2b[0][jm-1] = hh1*t2c[0][jm-1] +
	   hh3*t2d[0][jm-1]+hinv*t2h[0][jm-1] +
	   lf*hh1*t2e[0][jm-1]+lf*hh3*t2f[0][jm-1];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1] = t2c[im-1][jm-1] -
	   t2d[im-1][jm-1]+eig2*t2g[im-1][jm-1] +
	   h1inv*t2h[im-1][jm-1]+lf*t2e[im-1][jm-1] -
	   lf*t2f[im-1][jm-1];
     t2b[im-1][jm-1] = hh1*t2c[im-1][jm-1] +
	   hh3*t2d[im-1][jm-1]+hinv*t2h[im-1][jm-1] +
	   lf*hh1*t2e[im-1][jm-1] +
	   lf*hh3*t2f[im-1][jm-1];
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) t2b[0];
     t1c = (double *) t2c[0];
     t1d = (double *) t2d[0];
     t1e = (double *) t2e[0];
     t1f = (double *) t2f[0];
     t1g = (double *) t2g[0];
     t1h = (double *) t2h[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1c[j]-t1d[j] +
	   eig2*t1g[j]+h1inv*t1h[j] +
	   lf*t1e[j]-lf*t1f[j];
       t1b[j] = hh1*t1c[j] +
	   hh3*t1d[j]+hinv*t1h[j] +
	   lf*hh1*t1e[j]+lf*hh3*t1f[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) t2b[im-1];
     t1c = (double *) t2c[im-1];
     t1d = (double *) t2d[im-1];
     t1e = (double *) t2e[im-1];
     t1f = (double *) t2f[im-1];
     t1g = (double *) t2g[im-1];
     t1h = (double *) t2h[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1c[j] -
	   t1d[j]+eig2*t1g[j] +
	   h1inv*t1h[j]+lf*t1e[j] -
	   lf*t1f[j];
       t1b[j] = hh1*t1c[j] +
	   hh3*t1d[j]+hinv*t1h[j] +
	   lf*hh1*t1e[j]+lf*hh3*t1f[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = t2c[j][0]-t2d[j][0] +
	   eig2*t2g[j][0]+h1inv*t2h[j][0] +
	   lf*t2e[j][0]-lf*t2f[j][0];
       t2b[j][0] = hh1*t2c[j][0] +
	   hh3*t2d[j][0]+hinv*t2h[j][0] +
	   lf*hh1*t2e[j][0]+lf*hh3*t2f[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = t2c[j][jm-1] -
	   t2d[j][jm-1]+eig2*t2g[j][jm-1] +
	   h1inv*t2h[j][jm-1]+lf*t2e[j][jm-1] -
	   lf*t2f[j][jm-1];
       t2b[j][jm-1] = hh1*t2c[j][jm-1] +
	   hh3*t2d[j][jm-1]+hinv*t2h[j][jm-1] +
	   lf*hh1*t2e[j][jm-1]+lf*hh3*t2f[j][jm-1];
     }
   }

   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     t1c = (double *) t2c[i];
     t1d = (double *) t2d[i];
     t1e = (double *) t2e[i];
     t1f = (double *) t2f[i];
     t1g = (double *) t2g[i];
     t1h = (double *) t2h[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       t1a[iindex] = t1c[iindex] -
	   t1d[iindex]+eig2*t1g[iindex] +
	   h1inv*t1h[iindex]+lf*t1e[iindex] -
	   lf*t1f[iindex];
       t1b[iindex] = hh1*t1c[iindex] +
	   hh3*t1d[iindex]+hinv*t1h[iindex] +
	   lf*hh1*t1e[iindex] +
	   lf*hh3*t1f[iindex];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 712
	unsigned long	Error, Cycle;
#line 712
	int		Cancel, Temp;
#line 712

#line 712
	Error = pthread_mutex_lock(&(bars->sl_phase_5).mutex);
#line 712
	if (Error != 0) {
#line 712
		printf("Error while trying to get lock in barrier.\n");
#line 712
		exit(-1);
#line 712
	}
#line 712

#line 712
	Cycle = (bars->sl_phase_5).cycle;
#line 712
	if (++(bars->sl_phase_5).counter != (nprocs)) {
#line 712
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 712
		while (Cycle == (bars->sl_phase_5).cycle) {
#line 712
			Error = pthread_cond_wait(&(bars->sl_phase_5).cv, &(bars->sl_phase_5).mutex);
#line 712
			if (Error != 0) {
#line 712
				break;
#line 712
			}
#line 712
		}
#line 712
		pthread_setcancelstate(Cancel, &Temp);
#line 712
	} else {
#line 712
		(bars->sl_phase_5).cycle = !(bars->sl_phase_5).cycle;
#line 712
		(bars->sl_phase_5).counter = 0;
#line 712
		Error = pthread_cond_broadcast(&(bars->sl_phase_5).cv);
#line 712
	}
#line 712
	pthread_mutex_unlock(&(bars->sl_phase_5).mutex);
#line 712
}
#else
   {
#line 714
	unsigned long	Error, Cycle;
#line 714
	int		Cancel, Temp;
#line 714

#line 714
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 714
	if (Error != 0) {
#line 714
		printf("Error while trying to get lock in barrier.\n");
#line 714
		exit(-1);
#line 714
	}
#line 714

#line 714
	Cycle = (bars->barrier).cycle;
#line 714
	if (++(bars->barrier).counter != (nprocs)) {
#line 714
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 714
		while (Cycle == (bars->barrier).cycle) {
#line 714
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 714
			if (Error != 0) {
#line 714
				break;
#line 714
			}
#line 714
		}
#line 714
		pthread_setcancelstate(Cancel, &Temp);
#line 714
	} else {
#line 714
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 714
		(bars->barrier).counter = 0;
#line 714
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 714
	}
#line 714
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 714
}
#endif

   zsim_stamp();
/*     *******************************************************

               s i x t h   p h a s e

       *******************************************************  */

   istart = 1;
   iend = istart + gp[procid].rel_num_y[numlev-1] - 1;
   jstart = 1;
   jend = jstart + gp[procid].rel_num_x[numlev-1] - 1;
   ist = istart;
   ien = iend;
   jst = jstart;
   jen = jend;

   if (gp[procid].neighbors[UP] == -1) {
     istart = 0;
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     jstart = 0;
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     iend = im-1;
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     jend = jm-1;
   }
   t2a = (double **) rhs_multi[procid][numlev-1];
   t2b = (double **) ga[procid];
   t2c = (double **) oldga[procid];
   t2d = (double **) q_multi[procid][numlev-1];
   for(i=istart;i<=iend;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     for(j=jstart;j<=jend;j++) {
       t1a[j] = t1b[j] * ressqr;
     }
   }

   if (gp[procid].neighbors[UP] == -1) {
     t1d = (double *) t2d[0];
     t1b = (double *) t2b[0];
     for(j=jstart;j<=jend;j++) {
       t1d[j] = t1b[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1d = (double *) t2d[im-1];
     t1b = (double *) t2b[im-1];
     for(j=jstart;j<=jend;j++) {
       t1d[j] = t1b[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(i=istart;i<=iend;i++) {
       t2d[i][0] = t2b[i][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(i=istart;i<=iend;i++) {
       t2d[i][jm-1] = t2b[i][jm-1];
     }
   }

   fac = 1.0 / (4.0 - ressqr*eig2);
   for(i=ist;i<=ien;i++) {
     t1d = (double *) t2d[i];
     t1c = (double *) t2c[i];
     for(j=jst;j<=jen;j++) {
       t1d[j] = t1c[j];
     }
   }

   if ((procid == MASTER) || (do_stats)) {
     {
#line 792
	struct timeval	FullTime;
#line 792

#line 792
	gettimeofday(&FullTime, NULL);
#line 792
	(multi_start) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 792
};
   }

   multig(procid);

   if ((procid == MASTER) || (do_stats)) {
     {
#line 798
	struct timeval	FullTime;
#line 798

#line 798
	gettimeofday(&FullTime, NULL);
#line 798
	(multi_end) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 798
};
     gp[procid].multi_time += (multi_end - multi_start);
   }

/* the shared sum variable psiai is initialized to 0 at
   every time-step  */

   if (procid == MASTER) {
     global->psiai=0.0;
   }

/*  copy the solution for use as initial guess in next time-step  */

   for(i=istart;i<=iend;i++) {
     t1b = (double *) t2b[i];
     t1c = (double *) t2c[i];
     t1d = (double *) t2d[i];
     for(j=jstart;j<=jend;j++) {
       t1b[j] = t1d[j];
       t1c[j] = t1d[j];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 821
	unsigned long	Error, Cycle;
#line 821
	int		Cancel, Temp;
#line 821

#line 821
	Error = pthread_mutex_lock(&(bars->sl_phase_6).mutex);
#line 821
	if (Error != 0) {
#line 821
		printf("Error while trying to get lock in barrier.\n");
#line 821
		exit(-1);
#line 821
	}
#line 821

#line 821
	Cycle = (bars->sl_phase_6).cycle;
#line 821
	if (++(bars->sl_phase_6).counter != (nprocs)) {
#line 821
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 821
		while (Cycle == (bars->sl_phase_6).cycle) {
#line 821
			Error = pthread_cond_wait(&(bars->sl_phase_6).cv, &(bars->sl_phase_6).mutex);
#line 821
			if (Error != 0) {
#line 821
				break;
#line 821
			}
#line 821
		}
#line 821
		pthread_setcancelstate(Cancel, &Temp);
#line 821
	} else {
#line 821
		(bars->sl_phase_6).cycle = !(bars->sl_phase_6).cycle;
#line 821
		(bars->sl_phase_6).counter = 0;
#line 821
		Error = pthread_cond_broadcast(&(bars->sl_phase_6).cv);
#line 821
	}
#line 821
	pthread_mutex_unlock(&(bars->sl_phase_6).mutex);
#line 821
}
#else
   {
#line 823
	unsigned long	Error, Cycle;
#line 823
	int		Cancel, Temp;
#line 823

#line 823
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 823
	if (Error != 0) {
#line 823
		printf("Error while trying to get lock in barrier.\n");
#line 823
		exit(-1);
#line 823
	}
#line 823

#line 823
	Cycle = (bars->barrier).cycle;
#line 823
	if (++(bars->barrier).counter != (nprocs)) {
#line 823
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 823
		while (Cycle == (bars->barrier).cycle) {
#line 823
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 823
			if (Error != 0) {
#line 823
				break;
#line 823
			}
#line 823
		}
#line 823
		pthread_setcancelstate(Cancel, &Temp);
#line 823
	} else {
#line 823
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 823
		(bars->barrier).counter = 0;
#line 823
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 823
	}
#line 823
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 823
}
#endif

   zsim_stamp();
/*     *******************************************************

                s e v e n t h   p h a s e

       *******************************************************

   every process computes the running sum for its assigned portion
   in a private variable psiaipriv   */

   psiaipriv=0.0;
   t2a = (double **) ga[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     psiaipriv = psiaipriv + 0.25*(t2a[0][0]);
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     psiaipriv = psiaipriv + 0.25*(t2a[0][jm-1]);
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     psiaipriv=psiaipriv+0.25*(t2a[im-1][0]);
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     psiaipriv=psiaipriv+0.25*(t2a[im-1][jm-1]);
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     for(j=firstcol;j<=lastcol;j++) {
       psiaipriv = psiaipriv + 0.5*t1a[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       psiaipriv = psiaipriv + 0.5*t1a[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       psiaipriv = psiaipriv + 0.5*t2a[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       psiaipriv = psiaipriv + 0.5*t2a[j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       psiaipriv = psiaipriv + t1a[iindex];
     }
   }

/* after computing its private sum, every process adds that to the
   shared running sum psiai  */

   {pthread_mutex_lock(&(locks->psiailock));}
   global->psiai = global->psiai + psiaipriv;
   {pthread_mutex_unlock(&(locks->psiailock));}
#if defined(MULTIPLE_BARRIERS)
   {
#line 886
	unsigned long	Error, Cycle;
#line 886
	int		Cancel, Temp;
#line 886

#line 886
	Error = pthread_mutex_lock(&(bars->sl_phase_7).mutex);
#line 886
	if (Error != 0) {
#line 886
		printf("Error while trying to get lock in barrier.\n");
#line 886
		exit(-1);
#line 886
	}
#line 886

#line 886
	Cycle = (bars->sl_phase_7).cycle;
#line 886
	if (++(bars->sl_phase_7).counter != (nprocs)) {
#line 886
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 886
		while (Cycle == (bars->sl_phase_7).cycle) {
#line 886
			Error = pthread_cond_wait(&(bars->sl_phase_7).cv, &(bars->sl_phase_7).mutex);
#line 886
			if (Error != 0) {
#line 886
				break;
#line 886
			}
#line 886
		}
#line 886
		pthread_setcancelstate(Cancel, &Temp);
#line 886
	} else {
#line 886
		(bars->sl_phase_7).cycle = !(bars->sl_phase_7).cycle;
#line 886
		(bars->sl_phase_7).counter = 0;
#line 886
		Error = pthread_cond_broadcast(&(bars->sl_phase_7).cv);
#line 886
	}
#line 886
	pthread_mutex_unlock(&(bars->sl_phase_7).mutex);
#line 886
}
#else
   {
#line 888
	unsigned long	Error, Cycle;
#line 888
	int		Cancel, Temp;
#line 888

#line 888
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 888
	if (Error != 0) {
#line 888
		printf("Error while trying to get lock in barrier.\n");
#line 888
		exit(-1);
#line 888
	}
#line 888

#line 888
	Cycle = (bars->barrier).cycle;
#line 888
	if (++(bars->barrier).counter != (nprocs)) {
#line 888
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 888
		while (Cycle == (bars->barrier).cycle) {
#line 888
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 888
			if (Error != 0) {
#line 888
				break;
#line 888
			}
#line 888
		}
#line 888
		pthread_setcancelstate(Cancel, &Temp);
#line 888
	} else {
#line 888
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 888
		(bars->barrier).counter = 0;
#line 888
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 888
	}
#line 888
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 888
}
#endif
   zsim_stamp();
/*      *******************************************************

                e i g h t h   p h a s e

        *******************************************************

   augment ga(i,j) with [-psiai/psibi]*psib(i,j) */

   f4 = (-global->psiai)/(global->psibi);

   t2a = (double **) ga[procid];
   t2b = (double **) psib[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0] = t2a[0][0]+f4*t2b[0][0];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0] = t2a[im-1][0]+f4*t2b[im-1][0];
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1] = t2a[0][jm-1]+f4*t2b[0][jm-1];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1] = t2a[im-1][jm-1] +
			      f4*t2b[im-1][jm-1];
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) t2b[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1a[j]+f4*t1b[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) t2b[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1a[j]+f4*t1b[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = t2a[j][0]+f4*t2b[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = t2a[j][jm-1]+f4*t2b[j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       t1a[iindex] = t1a[iindex]+f4*t1b[iindex];
     }
   }

   t2a = (double **) rhs_multi[procid][numlev-1];
   t2b = (double **) gb[procid];
   t2c = (double **) oldgb[procid];
   t2d = (double **) q_multi[procid][numlev-1];
   for(i=istart;i<=iend;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     for(j=jstart;j<=jend;j++) {
       t1a[j] = t1b[j] * ressqr;
     }
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1d = (double *) t2d[0];
     t1b = (double *) t2b[0];
     for(j=jstart;j<=jend;j++) {
       t1d[j] = t1b[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1d = (double *) t2d[im-1];
     t1b = (double *) t2b[im-1];
     for(j=jstart;j<=jend;j++) {
       t1d[j] = t1b[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(i=istart;i<=iend;i++) {
       t2d[i][0] = t2b[i][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(i=istart;i<=iend;i++) {
       t2d[i][jm-1] = t2b[i][jm-1];
     }
   }

   fac = 1.0 / (4.0 - ressqr*eig2);
   for(i=ist;i<=ien;i++) {
     t1d = (double *) t2d[i];
     t1c = (double *) t2c[i];
     for(j=jst;j<=jen;j++) {
       t1d[j] = t1c[j];
     }
   }

   if ((procid == MASTER) || (do_stats)) {
     {
#line 994
	struct timeval	FullTime;
#line 994

#line 994
	gettimeofday(&FullTime, NULL);
#line 994
	(multi_start) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 994
};
   }

   multig(procid);

   if ((procid == MASTER) || (do_stats)) {
     {
#line 1000
	struct timeval	FullTime;
#line 1000

#line 1000
	gettimeofday(&FullTime, NULL);
#line 1000
	(multi_end) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 1000
};
     gp[procid].multi_time += (multi_end - multi_start);
   }

   for(i=istart;i<=iend;i++) {
     t1b = (double *) t2b[i];
     t1c = (double *) t2c[i];
     t1d = (double *) t2d[i];
     for(j=jstart;j<=jend;j++) {
       t1b[j] = t1d[j];
       t1c[j] = t1d[j];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 1014
	unsigned long	Error, Cycle;
#line 1014
	int		Cancel, Temp;
#line 1014

#line 1014
	Error = pthread_mutex_lock(&(bars->sl_phase_8).mutex);
#line 1014
	if (Error != 0) {
#line 1014
		printf("Error while trying to get lock in barrier.\n");
#line 1014
		exit(-1);
#line 1014
	}
#line 1014

#line 1014
	Cycle = (bars->sl_phase_8).cycle;
#line 1014
	if (++(bars->sl_phase_8).counter != (nprocs)) {
#line 1014
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 1014
		while (Cycle == (bars->sl_phase_8).cycle) {
#line 1014
			Error = pthread_cond_wait(&(bars->sl_phase_8).cv, &(bars->sl_phase_8).mutex);
#line 1014
			if (Error != 0) {
#line 1014
				break;
#line 1014
			}
#line 1014
		}
#line 1014
		pthread_setcancelstate(Cancel, &Temp);
#line 1014
	} else {
#line 1014
		(bars->sl_phase_8).cycle = !(bars->sl_phase_8).cycle;
#line 1014
		(bars->sl_phase_8).counter = 0;
#line 1014
		Error = pthread_cond_broadcast(&(bars->sl_phase_8).cv);
#line 1014
	}
#line 1014
	pthread_mutex_unlock(&(bars->sl_phase_8).mutex);
#line 1014
}
#else
   {
#line 1016
	unsigned long	Error, Cycle;
#line 1016
	int		Cancel, Temp;
#line 1016

#line 1016
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 1016
	if (Error != 0) {
#line 1016
		printf("Error while trying to get lock in barrier.\n");
#line 1016
		exit(-1);
#line 1016
	}
#line 1016

#line 1016
	Cycle = (bars->barrier).cycle;
#line 1016
	if (++(bars->barrier).counter != (nprocs)) {
#line 1016
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 1016
		while (Cycle == (bars->barrier).cycle) {
#line 1016
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 1016
			if (Error != 0) {
#line 1016
				break;
#line 1016
			}
#line 1016
		}
#line 1016
		pthread_setcancelstate(Cancel, &Temp);
#line 1016
	} else {
#line 1016
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 1016
		(bars->barrier).counter = 0;
#line 1016
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 1016
	}
#line 1016
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 1016
}
#endif
   zsim_stamp();
/*      *******************************************************

                n i n t h   p h a s e

        *******************************************************

   put appropriate linear combinations of ga and gb in work2 and work3;
   note that here (as in most cases) the constant multipliers are made
   private variables; the specific order in which things are done is
   chosen in order to hopefully reuse things brought into the cache

   note that here again we choose to have all processes share the work
   on both matrices despite the fact that the work done per element
   is the same, because the operand matrices are the same in both cases */

   t2a = (double **) ga[procid];
   t2b = (double **) gb[procid];
   t2c = (double **) work2[procid];
   t2d = (double **) work3[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2c[0][0] = t2b[0][0]-hh1*t2a[0][0];
     t2d[0][0] = t2b[0][0]+hh3*t2a[0][0];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2c[im-1][0] = t2b[im-1][0]-hh1*t2a[im-1][0];
     t2d[im-1][0] = t2b[im-1][0]+hh3*t2a[im-1][0];
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2c[0][jm-1] = t2b[0][jm-1]-hh1*t2a[0][jm-1];
     t2d[0][jm-1] = t2b[0][jm-1]+hh3*t2a[0][jm-1];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2c[im-1][jm-1] = t2b[im-1][jm-1] -
				 hh1*t2a[im-1][jm-1];
     t2d[im-1][jm-1] = t2b[im-1][jm-1] +
				 hh3*t2a[im-1][jm-1];
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) t2b[0];
     t1c = (double *) t2c[0];
     t1d = (double *) t2d[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1d[j] = t1b[j]+hh3*t1a[j];
       t1c[j] = t1b[j]-hh1*t1a[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) t2b[im-1];
     t1c = (double *) t2c[im-1];
     t1d = (double *) t2d[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1d[j] = t1b[j]+hh3*t1a[j];
       t1c[j] = t1b[j]-hh1*t1a[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2d[j][0] = t2b[j][0]+hh3*t2a[j][0];
       t2c[j][0] = t2b[j][0]-hh1*t2a[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2d[j][jm-1] = t2b[j][jm-1]+hh3*t2a[j][jm-1];
       t2c[j][jm-1] = t2b[j][jm-1]-hh1*t2a[j][jm-1];
     }
   }

   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     t1c = (double *) t2c[i];
     t1d = (double *) t2d[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       t1d[iindex] = t1b[iindex] + hh3*t1a[iindex];
       t1c[iindex] = t1b[iindex] - hh1*t1a[iindex];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 1100
	unsigned long	Error, Cycle;
#line 1100
	int		Cancel, Temp;
#line 1100

#line 1100
	Error = pthread_mutex_lock(&(bars->sl_phase_9).mutex);
#line 1100
	if (Error != 0) {
#line 1100
		printf("Error while trying to get lock in barrier.\n");
#line 1100
		exit(-1);
#line 1100
	}
#line 1100

#line 1100
	Cycle = (bars->sl_phase_9).cycle;
#line 1100
	if (++(bars->sl_phase_9).counter != (nprocs)) {
#line 1100
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 1100
		while (Cycle == (bars->sl_phase_9).cycle) {
#line 1100
			Error = pthread_cond_wait(&(bars->sl_phase_9).cv, &(bars->sl_phase_9).mutex);
#line 1100
			if (Error != 0) {
#line 1100
				break;
#line 1100
			}
#line 1100
		}
#line 1100
		pthread_setcancelstate(Cancel, &Temp);
#line 1100
	} else {
#line 1100
		(bars->sl_phase_9).cycle = !(bars->sl_phase_9).cycle;
#line 1100
		(bars->sl_phase_9).counter = 0;
#line 1100
		Error = pthread_cond_broadcast(&(bars->sl_phase_9).cv);
#line 1100
	}
#line 1100
	pthread_mutex_unlock(&(bars->sl_phase_9).mutex);
#line 1100
}
#else
   {
#line 1102
	unsigned long	Error, Cycle;
#line 1102
	int		Cancel, Temp;
#line 1102

#line 1102
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 1102
	if (Error != 0) {
#line 1102
		printf("Error while trying to get lock in barrier.\n");
#line 1102
		exit(-1);
#line 1102
	}
#line 1102

#line 1102
	Cycle = (bars->barrier).cycle;
#line 1102
	if (++(bars->barrier).counter != (nprocs)) {
#line 1102
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 1102
		while (Cycle == (bars->barrier).cycle) {
#line 1102
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 1102
			if (Error != 0) {
#line 1102
				break;
#line 1102
			}
#line 1102
		}
#line 1102
		pthread_setcancelstate(Cancel, &Temp);
#line 1102
	} else {
#line 1102
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 1102
		(bars->barrier).counter = 0;
#line 1102
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 1102
	}
#line 1102
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 1102
}
#endif
   zsim_stamp();
/*      *******************************************************

                t e n t h    p h a s e

        *******************************************************/


   timst = 2*dtau;

/* update the psi{1,3} matrices by adding 2*dtau*work3 to each */

   t2a = (double **) psi[procid][0];
   t2b = (double **) work3[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0] = t2a[0][0] + timst*t2b[0][0];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0] = t2a[im-1][0] +
			       timst*t2b[im-1][0];
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1] = t2a[0][jm-1] +
			       timst*t2b[0][jm-1];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1] = t2a[im-1][jm-1] +
				  timst*t2b[im-1][jm-1];
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) t2b[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1a[j] + timst*t1b[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) t2b[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1a[j] + timst*t1b[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = t2a[j][0] + timst*t2b[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = t2a[j][jm-1] +
				 timst*t2b[j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = t1a[iindex] + timst*t1b[iindex];
     }
   }

   t2a = (double **) psi[procid][1];
   t2b = (double **) work2[procid];
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[0][0] = t2a[0][0] + timst*t2b[0][0];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[LEFT] == -1)) {
     t2a[im-1][0] = t2a[im-1][0] +
			       timst*t2b[im-1][0];
   }
   if ((gp[procid].neighbors[UP] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[0][jm-1] = t2a[0][jm-1] +
			       timst*t2b[0][jm-1];
   }
   if ((gp[procid].neighbors[DOWN] == -1) && (gp[procid].neighbors[RIGHT] == -1)) {
     t2a[im-1][jm-1] = t2a[im-1][jm-1] +
				  timst*t2b[im-1][jm-1];
   }
   if (gp[procid].neighbors[UP] == -1) {
     t1a = (double *) t2a[0];
     t1b = (double *) t2b[0];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1a[j] + timst*t1b[j];
     }
   }
   if (gp[procid].neighbors[DOWN] == -1) {
     t1a = (double *) t2a[im-1];
     t1b = (double *) t2b[im-1];
     for(j=firstcol;j<=lastcol;j++) {
       t1a[j] = t1a[j] + timst*t1b[j];
     }
   }
   if (gp[procid].neighbors[LEFT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][0] = t2a[j][0] + timst*t2b[j][0];
     }
   }
   if (gp[procid].neighbors[RIGHT] == -1) {
     for(j=firstrow;j<=lastrow;j++) {
       t2a[j][jm-1] = t2a[j][jm-1] +
				 timst*t2b[j][jm-1];
     }
   }

   for(i=firstrow;i<=lastrow;i++) {
     t1a = (double *) t2a[i];
     t1b = (double *) t2b[i];
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
         t1a[iindex] = t1a[iindex] + timst*t1b[iindex];
     }
   }

#if defined(MULTIPLE_BARRIERS)
   {
#line 1218
	unsigned long	Error, Cycle;
#line 1218
	int		Cancel, Temp;
#line 1218

#line 1218
	Error = pthread_mutex_lock(&(bars->sl_phase_10).mutex);
#line 1218
	if (Error != 0) {
#line 1218
		printf("Error while trying to get lock in barrier.\n");
#line 1218
		exit(-1);
#line 1218
	}
#line 1218

#line 1218
	Cycle = (bars->sl_phase_10).cycle;
#line 1218
	if (++(bars->sl_phase_10).counter != (nprocs)) {
#line 1218
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 1218
		while (Cycle == (bars->sl_phase_10).cycle) {
#line 1218
			Error = pthread_cond_wait(&(bars->sl_phase_10).cv, &(bars->sl_phase_10).mutex);
#line 1218
			if (Error != 0) {
#line 1218
				break;
#line 1218
			}
#line 1218
		}
#line 1218
		pthread_setcancelstate(Cancel, &Temp);
#line 1218
	} else {
#line 1218
		(bars->sl_phase_10).cycle = !(bars->sl_phase_10).cycle;
#line 1218
		(bars->sl_phase_10).counter = 0;
#line 1218
		Error = pthread_cond_broadcast(&(bars->sl_phase_10).cv);
#line 1218
	}
#line 1218
	pthread_mutex_unlock(&(bars->sl_phase_10).mutex);
#line 1218
}
#else
   {
#line 1220
	unsigned long	Error, Cycle;
#line 1220
	int		Cancel, Temp;
#line 1220

#line 1220
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 1220
	if (Error != 0) {
#line 1220
		printf("Error while trying to get lock in barrier.\n");
#line 1220
		exit(-1);
#line 1220
	}
#line 1220

#line 1220
	Cycle = (bars->barrier).cycle;
#line 1220
	if (++(bars->barrier).counter != (nprocs)) {
#line 1220
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 1220
		while (Cycle == (bars->barrier).cycle) {
#line 1220
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 1220
			if (Error != 0) {
#line 1220
				break;
#line 1220
			}
#line 1220
		}
#line 1220
		pthread_setcancelstate(Cancel, &Temp);
#line 1220
	} else {
#line 1220
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 1220
		(bars->barrier).counter = 0;
#line 1220
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 1220
	}
#line 1220
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 1220
}
#endif

zsim_stamp();
//zsim_PIM_function_end();
}
