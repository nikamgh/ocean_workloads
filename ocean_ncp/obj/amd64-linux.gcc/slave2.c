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

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "decs.h"
#include "/mnt/panzer/mnika/Workload_Pool/zsim_hooks.h"
#if defined(ZSIM_TRACE_1) || defined(ZSIM_TRACE_2) || defined(ZSIM_TRACE_3) || defined(ZSIM_TRACE_4)
// #include "zsim_hooks.h"
#endif
int roi_included_2 = 0;
void slave2(long procid, long firstrow, long lastrow, long numrows, long firstcol, long lastcol, long numcols)
{
    zsim_stamp();
//#ifdef ZSIM_TRACE_4
  //         zsim_roi_begin();
//	zsim_PIM_function_begin();
//#endif
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
   double timst;
   double f4;
   long psiindex;
   double psiaipriv;
   long multi_start;
   long multi_end;

   ressqr = lev_res[numlev-1] * lev_res[numlev-1];
   zsim_stamp();
/*   ***************************************************************

          f i r s t     p h a s e   (of timestep calculation)

     ***************************************************************/

   if (procid == MASTER) {
     wrk1->ga[0][0]=0.0;
   }
   if (procid == nprocs-xprocs) {
     wrk1->ga[im-1][0]=0.0;
   }
   if (procid == xprocs-1) {
     wrk1->ga[0][jm-1]=0.0;
   }
   if (procid == nprocs-1) {
     wrk1->ga[im-1][jm-1]=0.0;
   }
   if (firstrow == 1) {
     for(j=firstcol;j<=lastcol;j++) {
       wrk1->ga[0][j] = 0.0;
     }
   }
   if ((firstrow+numrows) == im-1) {
     for(j=firstcol;j<=lastcol;j++) {
       wrk1->ga[im-1][j] = 0.0;
     }
   }
   if (firstcol == 1) {
     for(j=firstrow;j<=lastrow;j++) {
       wrk1->ga[j][0] = 0.0;
     }
   }
   if ((firstcol+numcols) == jm-1) {
     for(j=firstrow;j<=lastrow;j++) {
       wrk1->ga[j][jm-1] = 0.0;
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
         wrk1->ga[i][iindex] = 0.0;
     }
   }

   if (procid == MASTER) {
     wrk1->gb[0][0]=0.0;
   }
   if (procid == nprocs-xprocs) {
     wrk1->gb[im-1][0]=0.0;
   }
   if (procid == xprocs-1) {
     wrk1->gb[0][jm-1]=0.0;
   }
   if (procid == nprocs-1) {
     wrk1->gb[im-1][jm-1]=0.0;
   }
   if (firstrow == 1) {
     for(j=firstcol;j<=lastcol;j++) {
       wrk1->gb[0][j] = 0.0;
     }
   }
   if ((firstrow+numrows) == im-1) {
     for(j=firstcol;j<=lastcol;j++) {
       wrk1->gb[im-1][j] = 0.0;
     }
   }
   if (firstcol == 1) {
     for(j=firstrow;j<=lastrow;j++) {
       wrk1->gb[j][0] = 0.0;
     }
   }
   if ((firstcol+numcols) == jm-1) {
     for(j=firstrow;j<=lastrow;j++) {
       wrk1->gb[j][jm-1] = 0.0;
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       wrk1->gb[i][iindex] = 0.0;
     }
   }

/* put the laplacian of psi{1,3} in work1{1,2}
   note that psi(i,j,2) represents the psi3 array in
   the original equations  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     if (procid == MASTER) {
       wrk3->work1[psiindex][0][0] = 0;
     }
     if (procid == nprocs-xprocs) {
       wrk3->work1[psiindex][im-1][0] = 0;
     }
     if (procid == xprocs-1) {
       wrk3->work1[psiindex][0][jm-1] = 0;
     }
     if (procid == nprocs-1) {
       wrk3->work1[psiindex][im-1][jm-1] = 0;
     }
     laplacalc(fields->psi[psiindex],
	       wrk3->work1[psiindex],
	       firstrow,lastrow,firstcol,lastcol,numrows,numcols);

   }


   if (procid == MASTER) {
     wrk3->work2[0][0] = fields->psi[0][0][0]-fields->psi[1][0][0];
   }
   if (procid == nprocs-xprocs) {
     wrk3->work2[im-1][0] = fields->psi[0][im-1][0]-fields->psi[1][im-1][0];
   }
   if (procid == xprocs-1) {
     wrk3->work2[0][jm-1] = fields->psi[0][0][jm-1]-fields->psi[1][0][jm-1];
   }
   if (procid == nprocs-1) {
     wrk3->work2[im-1][jm-1] = fields->psi[0][im-1][jm-1]-fields->psi[1][im-1][jm-1];
   }
   if (firstrow == 1) {
     for(j=firstcol;j<=lastcol;j++) {
       wrk3->work2[0][j] = fields->psi[0][0][j]-fields->psi[1][0][j];
     }
   }
   if ((firstrow+numrows) == im-1) {
     for(j=firstcol;j<=lastcol;j++) {
       wrk3->work2[im-1][j] = fields->psi[0][im-1][j]-fields->psi[1][im-1][j];
     }
   }
   if (firstcol == 1) {
     for(j=firstrow;j<=lastrow;j++) {
       wrk3->work2[j][0] = fields->psi[0][j][0]-fields->psi[1][j][0];
     }
   }
   if ((firstcol+numcols) == jm-1) {
     for(j=firstrow;j<=lastrow;j++) {
       wrk3->work2[j][jm-1] = fields->psi[0][j][jm-1]-fields->psi[1][j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
         wrk3->work2[i][iindex] = fields->psi[0][i][iindex]-fields->psi[1][i][iindex];
     }
   }

/* set values of work3 array to h3/h * psi1 + h1/h * psi3  */

   hh3 = h3/h;
   hh1 = h1/h;
   if (procid == MASTER) {
     wrk2->work3[0][0] = hh3*fields->psi[0][0][0]+hh1*fields->psi[1][0][0];
   }
   if (procid == nprocs-xprocs) {
     wrk2->work3[im-1][0] = hh3*fields->psi[0][im-1][0]+hh1*fields->psi[1][im-1][0];
   }
   if (procid == xprocs-1) {
     wrk2->work3[0][jm-1] = hh3*fields->psi[0][0][jm-1]+hh1*fields->psi[1][0][jm-1];
   }
   if (procid == nprocs-1) {
     wrk2->work3[im-1][jm-1] = hh3*fields->psi[0][im-1][jm-1]+hh1*fields->psi[1][im-1][jm-1];
   }
   if (firstrow == 1) {
     for(j=firstcol;j<=lastcol;j++) {
       wrk2->work3[0][j] = hh3*fields->psi[0][0][j]+hh1*fields->psi[1][0][j];
     }
   }
   if ((firstrow+numrows) == im-1) {
     for(j=firstcol;j<=lastcol;j++) {
       wrk2->work3[im-1][j] = hh3*fields->psi[0][im-1][j]+hh1*fields->psi[1][im-1][j];
     }
   }
   if (firstcol == 1) {
     for(j=firstrow;j<=lastrow;j++) {
       wrk2->work3[j][0] = hh3*fields->psi[0][j][0]+hh1*fields->psi[1][j][0];
     }
   }
   if ((firstcol+numcols) == jm-1) {
     for(j=firstrow;j<=lastrow;j++) {
       wrk2->work3[j][jm-1] = hh3*fields->psi[0][j][jm-1]+hh1*fields->psi[1][j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
        wrk2->work3[i][iindex] = hh3*fields->psi[0][i][iindex]+hh1*fields->psi[1][i][iindex];
     }
   }

/* set values of temparray{1,3} to psim{1,3}  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     if (procid == MASTER) {
       wrk5->temparray[psiindex][0][0] = fields->psi[psiindex][0][0];
     }
     if (procid == nprocs-xprocs) {
       wrk5->temparray[psiindex][im-1][0] = fields->psi[psiindex][im-1][0];
     }
     if (procid == xprocs-1) {
       wrk5->temparray[psiindex][0][jm-1] = fields->psi[psiindex][0][jm-1];
     }
     if (procid == nprocs-1) {
       wrk5->temparray[psiindex][im-1][jm-1] = fields->psi[psiindex][im-1][jm-1];
     }
     if (firstrow == 1) {
       for(j=firstcol;j<=lastcol;j++) {
         wrk5->temparray[psiindex][0][j] = fields->psi[psiindex][0][j];
       }
     }
     if ((firstrow+numrows) == im-1) {
       for(j=firstcol;j<=lastcol;j++) {
         wrk5->temparray[psiindex][im-1][j] = fields->psi[psiindex][im-1][j];
       }
     }
     if (firstcol == 1) {
       for(j=firstrow;j<=lastrow;j++) {
         wrk5->temparray[psiindex][j][0] = fields->psi[psiindex][j][0];
       }
     }
     if ((firstcol+numcols) == jm-1) {
       for(j=firstrow;j<=lastrow;j++) {
         wrk5->temparray[psiindex][j][jm-1] = fields->psi[psiindex][j][jm-1];
       }
     }

     for(i=firstrow;i<=lastrow;i++) {
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         wrk5->temparray[psiindex][i][iindex] = fields->psi[psiindex][i][iindex];
       }
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 293
	unsigned long	Error, Cycle;
#line 293
	int		Cancel, Temp;
#line 293

#line 293
	Error = pthread_mutex_lock(&(bars->sl_phase_1).mutex);
#line 293
	if (Error != 0) {
#line 293
		printf("Error while trying to get lock in barrier.\n");
#line 293
		exit(-1);
#line 293
	}
#line 293

#line 293
	Cycle = (bars->sl_phase_1).cycle;
#line 293
	if (++(bars->sl_phase_1).counter != (nprocs)) {
#line 293
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 293
		while (Cycle == (bars->sl_phase_1).cycle) {
#line 293
			Error = pthread_cond_wait(&(bars->sl_phase_1).cv, &(bars->sl_phase_1).mutex);
#line 293
			if (Error != 0) {
#line 293
				break;
#line 293
			}
#line 293
		}
#line 293
		pthread_setcancelstate(Cancel, &Temp);
#line 293
	} else {
#line 293
		(bars->sl_phase_1).cycle = !(bars->sl_phase_1).cycle;
#line 293
		(bars->sl_phase_1).counter = 0;
#line 293
		Error = pthread_cond_broadcast(&(bars->sl_phase_1).cv);
#line 293
	}
#line 293
	pthread_mutex_unlock(&(bars->sl_phase_1).mutex);
#line 293
}
#else
   {
#line 295
	unsigned long	Error, Cycle;
#line 295
	int		Cancel, Temp;
#line 295

#line 295
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 295
	if (Error != 0) {
#line 295
		printf("Error while trying to get lock in barrier.\n");
#line 295
		exit(-1);
#line 295
	}
#line 295

#line 295
	Cycle = (bars->barrier).cycle;
#line 295
	if (++(bars->barrier).counter != (nprocs)) {
#line 295
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 295
		while (Cycle == (bars->barrier).cycle) {
#line 295
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 295
			if (Error != 0) {
#line 295
				break;
#line 295
			}
#line 295
		}
#line 295
		pthread_setcancelstate(Cancel, &Temp);
#line 295
	} else {
#line 295
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 295
		(bars->barrier).counter = 0;
#line 295
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 295
	}
#line 295
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 295
}
#endif
   zsim_stamp();
/*     *******************************************************

              s e c o n d   p h a s e

       *******************************************************

   set values of psi{1,3} to psim{1,3}   */

   for(psiindex=0;psiindex<=1;psiindex++) {
     if (procid == MASTER) {
       fields->psi[psiindex][0][0] = fields->psim[psiindex][0][0];
     }
     if (procid == xprocs-1) {
       fields->psi[psiindex][0][jm-1] = fields->psim[psiindex][0][jm-1];
     }
     if (procid == nprocs-xprocs) {
       fields->psi[psiindex][im-1][0] = fields->psim[psiindex][im-1][0];
     }
     if (procid == nprocs-1) {
       fields->psi[psiindex][im-1][jm-1] = fields->psim[psiindex][im-1][jm-1];
     }
     if (firstrow == 1) {
       for(j=firstcol;j<=lastcol;j++) {
         fields->psi[psiindex][0][j] = fields->psim[psiindex][0][j];
       }
     }
     if ((firstrow+numrows) == im-1) {
       for(j=firstcol;j<=lastcol;j++) {
         fields->psi[psiindex][im-1][j] = fields->psim[psiindex][im-1][j];
       }
     }
     if (firstcol == 1) {
       for(j=firstrow;j<=lastrow;j++) {
         fields->psi[psiindex][j][0] = fields->psim[psiindex][j][0];
       }
     }
     if ((firstcol+numcols) == jm-1) {
       for(j=firstrow;j<=lastrow;j++) {
         fields->psi[psiindex][j][jm-1] = fields->psim[psiindex][j][jm-1];
       }
     }

     for(i=firstrow;i<=lastrow;i++) {
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         fields->psi[psiindex][i][iindex] = fields->psim[psiindex][i][iindex];
       }
     }
   }

/* put the laplacian of the psim array
   into the work7 array; first part of a three-laplacian
   calculation to compute the friction terms  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     if (procid == MASTER) {
       wrk5->work7[psiindex][0][0] = 0;
     }
     if (procid == nprocs-xprocs) {
       wrk5->work7[psiindex][im-1][0] = 0;
     }
     if (procid == xprocs-1) {
       wrk5->work7[psiindex][0][jm-1] = 0;
     }
     if (procid == nprocs-1) {
       wrk5->work7[psiindex][im-1][jm-1] = 0;
     }
     laplacalc(fields->psim[psiindex],wrk5->work7[psiindex],firstrow,lastrow,firstcol,lastcol,numrows,numcols);
   }

/* to the values of the work1{1,2} arrays obtained from the
   laplacians of psi{1,2} in the previous phase, add to the
   elements of every column the corresponding value in the
   one-dimenional f array  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     if (procid == MASTER) {
       wrk3->work1[psiindex][0][0] = wrk3->work1[psiindex][0][0] + wrk2->f[0];
     }
     if (procid == nprocs-xprocs) {
       wrk3->work1[psiindex][im-1][0] = wrk3->work1[psiindex][im-1][0] + wrk2->f[0];
     }
     if (procid == xprocs-1) {
       wrk3->work1[psiindex][0][jm-1] = wrk3->work1[psiindex][0][jm-1] + wrk2->f[jm-1];
     }
     if (procid == nprocs-1) {
       wrk3->work1[psiindex][im-1][jm-1] = wrk3->work1[psiindex][im-1][jm-1] + wrk2->f[jm-1];
     }
     if (firstrow == 1) {
       for(j=firstcol;j<=lastcol;j++) {
         wrk3->work1[psiindex][0][j] = wrk3->work1[psiindex][0][j] + wrk2->f[j];
       }
     }
     if ((firstrow+numrows) == im-1) {
       for(j=firstcol;j<=lastcol;j++) {
         wrk3->work1[psiindex][im-1][j] = wrk3->work1[psiindex][im-1][j] + wrk2->f[j];
       }
     }
     if (firstcol == 1) {
       for(j=firstrow;j<=lastrow;j++) {
         wrk3->work1[psiindex][j][0] = wrk3->work1[psiindex][j][0] + wrk2->f[j];
       }
     }
     if ((firstcol+numcols) == jm-1) {
       for(j=firstrow;j<=lastrow;j++) {
         wrk3->work1[psiindex][j][jm-1] = wrk3->work1[psiindex][j][jm-1] + wrk2->f[j];
       }
     }
     for(i=firstrow;i<=lastrow;i++) {
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         wrk3->work1[psiindex][i][iindex] = wrk3->work1[psiindex][i][iindex] +
					   wrk2->f[iindex];
       }
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 413
	unsigned long	Error, Cycle;
#line 413
	int		Cancel, Temp;
#line 413

#line 413
	Error = pthread_mutex_lock(&(bars->sl_phase_2).mutex);
#line 413
	if (Error != 0) {
#line 413
		printf("Error while trying to get lock in barrier.\n");
#line 413
		exit(-1);
#line 413
	}
#line 413

#line 413
	Cycle = (bars->sl_phase_2).cycle;
#line 413
	if (++(bars->sl_phase_2).counter != (nprocs)) {
#line 413
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 413
		while (Cycle == (bars->sl_phase_2).cycle) {
#line 413
			Error = pthread_cond_wait(&(bars->sl_phase_2).cv, &(bars->sl_phase_2).mutex);
#line 413
			if (Error != 0) {
#line 413
				break;
#line 413
			}
#line 413
		}
#line 413
		pthread_setcancelstate(Cancel, &Temp);
#line 413
	} else {
#line 413
		(bars->sl_phase_2).cycle = !(bars->sl_phase_2).cycle;
#line 413
		(bars->sl_phase_2).counter = 0;
#line 413
		Error = pthread_cond_broadcast(&(bars->sl_phase_2).cv);
#line 413
	}
#line 413
	pthread_mutex_unlock(&(bars->sl_phase_2).mutex);
#line 413
}
#else
   {
#line 415
	unsigned long	Error, Cycle;
#line 415
	int		Cancel, Temp;
#line 415

#line 415
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 415
	if (Error != 0) {
#line 415
		printf("Error while trying to get lock in barrier.\n");
#line 415
		exit(-1);
#line 415
	}
#line 415

#line 415
	Cycle = (bars->barrier).cycle;
#line 415
	if (++(bars->barrier).counter != (nprocs)) {
#line 415
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 415
		while (Cycle == (bars->barrier).cycle) {
#line 415
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 415
			if (Error != 0) {
#line 415
				break;
#line 415
			}
#line 415
		}
#line 415
		pthread_setcancelstate(Cancel, &Temp);
#line 415
	} else {
#line 415
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 415
		(bars->barrier).counter = 0;
#line 415
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 415
	}
#line 415
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 415
}
#endif
   zsim_stamp();
/* 	*******************************************************

                 t h i r d   p h a s e

 	*******************************************************

   put the jacobian of the work1{1,2} and psi{1,3} arrays
   (the latter currently in temparray) in the work5{1,2} arrays  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     jacobcalc(wrk3->work1[psiindex],wrk5->temparray[psiindex],
               wrk4->work5[psiindex],procid,firstrow,lastrow,firstcol,lastcol,numrows,numcols);
   }


/* set values of psim{1,3} to temparray{1,3}  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     if (procid == MASTER) {
       fields->psim[psiindex][0][0] = wrk5->temparray[psiindex][0][0];
     }
     if (procid == nprocs-xprocs) {
       fields->psim[psiindex][im-1][0] = wrk5->temparray[psiindex][im-1][0];
     }
     if (procid == xprocs-1) {
       fields->psim[psiindex][0][jm-1] = wrk5->temparray[psiindex][0][jm-1];
     }
     if (procid == nprocs-1) {
       fields->psim[psiindex][im-1][jm-1] = wrk5->temparray[psiindex][im-1][jm-1];
     }
     if (firstrow == 1) {
       for(j=firstcol;j<=lastcol;j++) {
         fields->psim[psiindex][0][j] = wrk5->temparray[psiindex][0][j];
       }
     }
     if ((firstrow+numrows) == im-1) {
       for(j=firstcol;j<=lastcol;j++) {
         fields->psim[psiindex][im-1][j] = wrk5->temparray[psiindex][im-1][j];
       }
     }
     if (firstcol == 1) {
       for(j=firstrow;j<=lastrow;j++) {
         fields->psim[psiindex][j][0] = wrk5->temparray[psiindex][j][0];
       }
     }
     if ((firstcol+numcols) == jm-1) {
       for(j=firstrow;j<=lastrow;j++) {
         fields->psim[psiindex][j][jm-1] = wrk5->temparray[psiindex][j][jm-1];
       }
     }
     for(i=firstrow;i<=lastrow;i++) {
       for(iindex=firstcol;iindex<=lastcol;iindex++) {
         fields->psim[psiindex][i][iindex] = wrk5->temparray[psiindex][i][iindex];
       }
     }
   }

/* put the laplacian of the work7{1,2} arrays in the work4{1,2}
   arrays; second step in the three-laplacian friction calculation  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     laplacalc(wrk5->work7[psiindex],
	       wrk4->work4[psiindex],
               firstrow,lastrow,firstcol,lastcol,numrows,numcols);
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 484
	unsigned long	Error, Cycle;
#line 484
	int		Cancel, Temp;
#line 484

#line 484
	Error = pthread_mutex_lock(&(bars->sl_phase_3).mutex);
#line 484
	if (Error != 0) {
#line 484
		printf("Error while trying to get lock in barrier.\n");
#line 484
		exit(-1);
#line 484
	}
#line 484

#line 484
	Cycle = (bars->sl_phase_3).cycle;
#line 484
	if (++(bars->sl_phase_3).counter != (nprocs)) {
#line 484
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 484
		while (Cycle == (bars->sl_phase_3).cycle) {
#line 484
			Error = pthread_cond_wait(&(bars->sl_phase_3).cv, &(bars->sl_phase_3).mutex);
#line 484
			if (Error != 0) {
#line 484
				break;
#line 484
			}
#line 484
		}
#line 484
		pthread_setcancelstate(Cancel, &Temp);
#line 484
	} else {
#line 484
		(bars->sl_phase_3).cycle = !(bars->sl_phase_3).cycle;
#line 484
		(bars->sl_phase_3).counter = 0;
#line 484
		Error = pthread_cond_broadcast(&(bars->sl_phase_3).cv);
#line 484
	}
#line 484
	pthread_mutex_unlock(&(bars->sl_phase_3).mutex);
#line 484
}
#else
   {
#line 486
	unsigned long	Error, Cycle;
#line 486
	int		Cancel, Temp;
#line 486

#line 486
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 486
	if (Error != 0) {
#line 486
		printf("Error while trying to get lock in barrier.\n");
#line 486
		exit(-1);
#line 486
	}
#line 486

#line 486
	Cycle = (bars->barrier).cycle;
#line 486
	if (++(bars->barrier).counter != (nprocs)) {
#line 486
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 486
		while (Cycle == (bars->barrier).cycle) {
#line 486
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 486
			if (Error != 0) {
#line 486
				break;
#line 486
			}
#line 486
		}
#line 486
		pthread_setcancelstate(Cancel, &Temp);
#line 486
	} else {
#line 486
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 486
		(bars->barrier).counter = 0;
#line 486
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 486
	}
#line 486
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 486
}
#endif
   zsim_stamp();
/*     *******************************************************

                f o u r t h   p h a s e

       *******************************************************

   put the jacobian of the work2 and work3 arrays in the work6
   array  */

   jacobcalc(wrk3->work2,wrk2->work3,wrk6->work6,procid,firstrow,
             lastrow,firstcol,lastcol,numrows,numcols);

/* put the laplacian of the work4{1,2} arrays in the work7{1,2}
   arrays; third step in the three-laplacian friction calculation  */

   for(psiindex=0;psiindex<=1;psiindex++) {
     laplacalc(wrk4->work4[psiindex],
               wrk5->work7[psiindex],
               firstrow,lastrow,firstcol,lastcol,numrows,numcols);
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 510
	unsigned long	Error, Cycle;
#line 510
	int		Cancel, Temp;
#line 510

#line 510
	Error = pthread_mutex_lock(&(bars->sl_phase_4).mutex);
#line 510
	if (Error != 0) {
#line 510
		printf("Error while trying to get lock in barrier.\n");
#line 510
		exit(-1);
#line 510
	}
#line 510

#line 510
	Cycle = (bars->sl_phase_4).cycle;
#line 510
	if (++(bars->sl_phase_4).counter != (nprocs)) {
#line 510
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 510
		while (Cycle == (bars->sl_phase_4).cycle) {
#line 510
			Error = pthread_cond_wait(&(bars->sl_phase_4).cv, &(bars->sl_phase_4).mutex);
#line 510
			if (Error != 0) {
#line 510
				break;
#line 510
			}
#line 510
		}
#line 510
		pthread_setcancelstate(Cancel, &Temp);
#line 510
	} else {
#line 510
		(bars->sl_phase_4).cycle = !(bars->sl_phase_4).cycle;
#line 510
		(bars->sl_phase_4).counter = 0;
#line 510
		Error = pthread_cond_broadcast(&(bars->sl_phase_4).cv);
#line 510
	}
#line 510
	pthread_mutex_unlock(&(bars->sl_phase_4).mutex);
#line 510
}
#else
   {
#line 512
	unsigned long	Error, Cycle;
#line 512
	int		Cancel, Temp;
#line 512

#line 512
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 512
	if (Error != 0) {
#line 512
		printf("Error while trying to get lock in barrier.\n");
#line 512
		exit(-1);
#line 512
	}
#line 512

#line 512
	Cycle = (bars->barrier).cycle;
#line 512
	if (++(bars->barrier).counter != (nprocs)) {
#line 512
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 512
		while (Cycle == (bars->barrier).cycle) {
#line 512
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 512
			if (Error != 0) {
#line 512
				break;
#line 512
			}
#line 512
		}
#line 512
		pthread_setcancelstate(Cancel, &Temp);
#line 512
	} else {
#line 512
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 512
		(bars->barrier).counter = 0;
#line 512
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 512
	}
#line 512
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 512
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

   if (procid == MASTER) {
     wrk1->ga[0][0] = wrk4->work5[0][0][0]-wrk4->work5[1][0][0]+eig2*wrk6->work6[0][0]+h1inv*
                frcng->tauz[0][0]+lf*wrk5->work7[0][0][0]-lf*wrk5->work7[1][0][0];
     wrk1->gb[0][0] = hh1*wrk4->work5[0][0][0]+hh3*wrk4->work5[1][0][0]+hinv*frcng->tauz[0][0]+
                lf*hh1*wrk5->work7[0][0][0]+lf*hh3*wrk5->work7[1][0][0];
   }
   if (procid == nprocs-xprocs) {
     wrk1->ga[im-1][0] = wrk4->work5[0][im-1][0]-wrk4->work5[1][im-1][0]+eig2*wrk6->work6[im-1][0]+h1inv*
                   frcng->tauz[im-1][0]+lf*wrk5->work7[0][im-1][0]-lf*wrk5->work7[1][im-1][0];
     wrk1->gb[im-1][0] = hh1*wrk4->work5[0][im-1][0]+hh3*wrk4->work5[1][im-1][0]+hinv*frcng->tauz[im-1][0]+
                   lf*hh1*wrk5->work7[0][im-1][0]+lf*hh3*wrk5->work7[1][im-1][0];
   }
   if (procid == xprocs-1) {
     wrk1->ga[0][jm-1] = wrk4->work5[0][0][jm-1]-wrk4->work5[1][0][jm-1]+eig2*wrk6->work6[0][jm-1]+h1inv*
                   frcng->tauz[0][jm-1]+lf*wrk5->work7[0][0][jm-1]-lf*wrk5->work7[1][0][jm-1];
     wrk1->gb[0][jm-1] = hh1*wrk4->work5[0][0][jm-1]+hh3*wrk4->work5[1][0][jm-1]+hinv*frcng->tauz[0][jm-1]+
                   lf*hh1*wrk5->work7[0][0][jm-1]+lf*hh3*wrk5->work7[1][0][jm-1];
   }
   if (procid == nprocs-1) {
     wrk1->ga[im-1][jm-1] = wrk4->work5[0][im-1][jm-1]-wrk4->work5[1][im-1][jm-1]+eig2*wrk6->work6[im-1][jm-1]+
                      h1inv*frcng->tauz[im-1][jm-1]+lf*wrk5->work7[0][im-1][jm-1]-lf*wrk5->work7[1][im-1][jm-1];
     wrk1->gb[im-1][jm-1] = hh1*wrk4->work5[0][im-1][jm-1]+hh3*wrk4->work5[1][im-1][jm-1]+hinv*
		    frcng->tauz[im-1][jm-1]+lf*hh1*wrk5->work7[0][im-1][jm-1]+lf*hh3*wrk5->work7[1][im-1][jm-1];
   }
   if (firstrow == 1) {
     for(j=firstcol;j<=lastcol;j++) {
       wrk1->ga[0][j] = wrk4->work5[0][0][j]-wrk4->work5[1][0][j]+eig2*
                    wrk6->work6[0][j]+h1inv*frcng->tauz[0][j]+lf*wrk5->work7[0][0][j]-lf*wrk5->work7[0][0][j];
       wrk1->gb[0][j] = hh1*wrk4->work5[0][0][j]+hh3*wrk4->work5[1][0][j]+hinv*
                    frcng->tauz[0][j]+lf*hh1*wrk5->work7[0][0][j]+lf*hh3*wrk5->work7[1][0][j];
     }
   }
   if ((firstrow+numrows) == im-1) {
     for(j=firstcol;j<=lastcol;j++) {
       wrk1->ga[im-1][j] = wrk4->work5[0][im-1][j]-wrk4->work5[1][im-1][j]+eig2*
                       wrk6->work6[im-1][j]+h1inv*frcng->tauz[im-1][j]+
                       lf*wrk5->work7[0][im-1][j]-lf*wrk5->work7[1][im-1][j];
       wrk1->gb[im-1][j] = hh1*wrk4->work5[0][im-1][j]+hh3*wrk4->work5[1][im-1][j]+hinv*
                       frcng->tauz[im-1][j]+lf*hh1*wrk5->work7[0][im-1][j]+
                       lf*hh3*wrk5->work7[1][im-1][j];
     }
   }
   if (firstcol == 1) {
     for(j=firstrow;j<=lastrow;j++) {
       wrk1->ga[j][0] = wrk4->work5[0][j][0]-wrk4->work5[1][j][0]+eig2*
                  wrk6->work6[j][0]+h1inv*frcng->tauz[j][0]+lf*wrk5->work7[0][j][0]-lf*wrk5->work7[1][j][0];
       wrk1->gb[j][0] = hh1*wrk4->work5[0][j][0]+hh3*wrk4->work5[1][j][0]+hinv*
                  frcng->tauz[j][0]+lf*hh1*wrk5->work7[0][j][0]+lf*hh3*wrk5->work7[1][j][0];
     }
   }
   if ((firstcol+numcols) == jm-1) {
     for(j=firstrow;j<=lastrow;j++) {
       wrk1->ga[j][jm-1] = wrk4->work5[0][j][jm-1]-wrk4->work5[1][j][jm-1]+eig2*
                       wrk6->work6[j][jm-1]+h1inv*frcng->tauz[j][jm-1]+
                       lf*wrk5->work7[0][j][jm-1]-lf*wrk5->work7[1][j][jm-1];
       wrk1->gb[j][jm-1] = hh1*wrk4->work5[0][j][jm-1]+hh3*wrk4->work5[1][j][jm-1]+hinv*
                       frcng->tauz[j][jm-1]+lf*hh1*wrk5->work7[0][j][jm-1]+
                       lf*hh3*wrk5->work7[1][j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       wrk1->ga[i][iindex] = wrk4->work5[0][i][iindex]-wrk4->work5[1][i][iindex]+eig2*
                          wrk6->work6[i][iindex]+h1inv*frcng->tauz[i][iindex]+
                          lf*wrk5->work7[0][i][iindex]-lf*wrk5->work7[1][i][iindex];
       wrk1->gb[i][iindex] = hh1*wrk4->work5[0][i][iindex]+hh3*wrk4->work5[1][i][iindex]+hinv*
                          frcng->tauz[i][iindex]+lf*hh1*wrk5->work7[0][i][iindex]+
                          lf*hh3*wrk5->work7[1][i][iindex];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 599
	unsigned long	Error, Cycle;
#line 599
	int		Cancel, Temp;
#line 599

#line 599
	Error = pthread_mutex_lock(&(bars->sl_phase_5).mutex);
#line 599
	if (Error != 0) {
#line 599
		printf("Error while trying to get lock in barrier.\n");
#line 599
		exit(-1);
#line 599
	}
#line 599

#line 599
	Cycle = (bars->sl_phase_5).cycle;
#line 599
	if (++(bars->sl_phase_5).counter != (nprocs)) {
#line 599
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 599
		while (Cycle == (bars->sl_phase_5).cycle) {
#line 599
			Error = pthread_cond_wait(&(bars->sl_phase_5).cv, &(bars->sl_phase_5).mutex);
#line 599
			if (Error != 0) {
#line 599
				break;
#line 599
			}
#line 599
		}
#line 599
		pthread_setcancelstate(Cancel, &Temp);
#line 599
	} else {
#line 599
		(bars->sl_phase_5).cycle = !(bars->sl_phase_5).cycle;
#line 599
		(bars->sl_phase_5).counter = 0;
#line 599
		Error = pthread_cond_broadcast(&(bars->sl_phase_5).cv);
#line 599
	}
#line 599
	pthread_mutex_unlock(&(bars->sl_phase_5).mutex);
#line 599
}
#else
   {
#line 601
	unsigned long	Error, Cycle;
#line 601
	int		Cancel, Temp;
#line 601

#line 601
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 601
	if (Error != 0) {
#line 601
		printf("Error while trying to get lock in barrier.\n");
#line 601
		exit(-1);
#line 601
	}
#line 601

#line 601
	Cycle = (bars->barrier).cycle;
#line 601
	if (++(bars->barrier).counter != (nprocs)) {
#line 601
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 601
		while (Cycle == (bars->barrier).cycle) {
#line 601
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 601
			if (Error != 0) {
#line 601
				break;
#line 601
			}
#line 601
		}
#line 601
		pthread_setcancelstate(Cancel, &Temp);
#line 601
	} else {
#line 601
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 601
		(bars->barrier).counter = 0;
#line 601
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 601
	}
#line 601
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 601
}
#endif
   zsim_stamp();
/*     *******************************************************

               s i x t h   p h a s e

       *******************************************************  */

   istart = gp[procid].rel_start_y[numlev-1];
   iend = istart + gp[procid].rel_num_y[numlev-1] - 1;
   jstart = gp[procid].rel_start_x[numlev-1];
   jend = jstart + gp[procid].rel_num_x[numlev-1] - 1;
   ist = istart;
   ien = iend;
   jst = jstart;
   jen = jend;
   if (istart == 1) {
     istart = 0;
   }
   if (jstart == 1) {
     jstart = 0;
   }
   if (iend == im-2) {
     iend = im-1;
   }
   if (jend == jm-2) {
     jend = jm-1;
   }
   for(i=istart;i<=iend;i++) {
     for(j=jstart;j<=jend;j++) {
       multi->rhs_multi[numlev-1][i][j] = wrk1->ga[i][j] * ressqr;
     }
   }
   if (istart == 0) {
     for(j=jstart;j<=jend;j++) {
       multi->q_multi[numlev-1][0][j] = wrk1->ga[0][j];
     }
   }
   if (iend == im-1) {
     for(j=jstart;j<=jend;j++) {
       multi->q_multi[numlev-1][im-1][j] = wrk1->ga[im-1][j];
     }
   }
   if (jstart == 0) {
     for(i=istart;i<=iend;i++) {
       multi->q_multi[numlev-1][i][0] = wrk1->ga[i][0];
     }
   }
   if (jend == jm-1) {
     for(i=istart;i<=iend;i++) {
       multi->q_multi[numlev-1][i][jm-1] = wrk1->ga[i][jm-1];
     }
   }

   fac = 1.0 / (4.0 - ressqr*eig2);
   for(i=ist;i<=ien;i++) {
     for(j=jst;j<=jen;j++) {
       multi->q_multi[numlev-1][i][j] = guess->oldga[i][j];
     }
   }

   if ((procid == MASTER) || (do_stats)) {
     {
#line 664
	struct timeval	FullTime;
#line 664

#line 664
	gettimeofday(&FullTime, NULL);
#line 664
	(multi_start) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 664
};
   }

   multig(procid);

   if ((procid == MASTER) || (do_stats)) {
     {
#line 670
	struct timeval	FullTime;
#line 670

#line 670
	gettimeofday(&FullTime, NULL);
#line 670
	(multi_end) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 670
};
     gp[procid].multi_time += (multi_end - multi_start);
   }

   if (procid == MASTER) {
     global->psiai=0.0;
   }

/*  copy the solution for use as initial guess in next time-step  */

   for(i=istart;i<=iend;i++) {
     for(j=jstart;j<=jend;j++) {
       wrk1->ga[i][j] = multi->q_multi[numlev-1][i][j];
       guess->oldga[i][j] = multi->q_multi[numlev-1][i][j];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 687
	unsigned long	Error, Cycle;
#line 687
	int		Cancel, Temp;
#line 687

#line 687
	Error = pthread_mutex_lock(&(bars->sl_phase_6).mutex);
#line 687
	if (Error != 0) {
#line 687
		printf("Error while trying to get lock in barrier.\n");
#line 687
		exit(-1);
#line 687
	}
#line 687

#line 687
	Cycle = (bars->sl_phase_6).cycle;
#line 687
	if (++(bars->sl_phase_6).counter != (nprocs)) {
#line 687
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 687
		while (Cycle == (bars->sl_phase_6).cycle) {
#line 687
			Error = pthread_cond_wait(&(bars->sl_phase_6).cv, &(bars->sl_phase_6).mutex);
#line 687
			if (Error != 0) {
#line 687
				break;
#line 687
			}
#line 687
		}
#line 687
		pthread_setcancelstate(Cancel, &Temp);
#line 687
	} else {
#line 687
		(bars->sl_phase_6).cycle = !(bars->sl_phase_6).cycle;
#line 687
		(bars->sl_phase_6).counter = 0;
#line 687
		Error = pthread_cond_broadcast(&(bars->sl_phase_6).cv);
#line 687
	}
#line 687
	pthread_mutex_unlock(&(bars->sl_phase_6).mutex);
#line 687
}
#else
   {
#line 689
	unsigned long	Error, Cycle;
#line 689
	int		Cancel, Temp;
#line 689

#line 689
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 689
	if (Error != 0) {
#line 689
		printf("Error while trying to get lock in barrier.\n");
#line 689
		exit(-1);
#line 689
	}
#line 689

#line 689
	Cycle = (bars->barrier).cycle;
#line 689
	if (++(bars->barrier).counter != (nprocs)) {
#line 689
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 689
		while (Cycle == (bars->barrier).cycle) {
#line 689
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 689
			if (Error != 0) {
#line 689
				break;
#line 689
			}
#line 689
		}
#line 689
		pthread_setcancelstate(Cancel, &Temp);
#line 689
	} else {
#line 689
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 689
		(bars->barrier).counter = 0;
#line 689
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 689
	}
#line 689
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 689
}
#endif
   zsim_stamp();
/*     *******************************************************

                s e v e n t h   p h a s e

       *******************************************************

   every process computes the running sum for its assigned portion
   in a private variable psiaipriv   */



   psiaipriv=0.0;
   if (procid == MASTER) {
     psiaipriv = psiaipriv + 0.25*(wrk1->ga[0][0]);
   }
   if (procid == xprocs - 1) {
     psiaipriv = psiaipriv + 0.25*(wrk1->ga[0][jm-1]);
   }
   if (procid == nprocs-xprocs) {
     psiaipriv=psiaipriv+0.25*(wrk1->ga[im-1][0]);
   }
   if (procid == nprocs-1) {
     psiaipriv=psiaipriv+0.25*(wrk1->ga[im-1][jm-1]);
   }
   if (firstrow == 1) {
     for(j=firstcol;j<=lastcol;j++) {
       psiaipriv = psiaipriv + 0.5*wrk1->ga[0][j];
     }
   }
   if ((firstrow+numrows) == im-1) {
     for(j=firstcol;j<=lastcol;j++) {
       psiaipriv = psiaipriv + 0.5*wrk1->ga[im-1][j];
     }
   }
   if (firstcol == 1) {
     for(j=firstrow;j<=lastrow;j++) {
       psiaipriv = psiaipriv + 0.5*wrk1->ga[j][0];
     }
   }
   if ((firstcol+numcols) == jm-1) {
     for(j=firstrow;j<=lastrow;j++) {
       psiaipriv = psiaipriv + 0.5*wrk1->ga[j][jm-1];
     }
   }
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
   for(i=firstrow;i<=lastrow;i++) {
       psiaipriv = psiaipriv + wrk1->ga[i][iindex];
     }
   }

/* after computing its private sum, every process adds that to the
   shared running sum psiai  */

   {pthread_mutex_lock(&(locks->psibilock));}
   global->psiai = global->psiai + psiaipriv;
   {pthread_mutex_unlock(&(locks->psibilock));}
#if defined(MULTIPLE_BARRIERS)
   {
#line 749
	unsigned long	Error, Cycle;
#line 749
	int		Cancel, Temp;
#line 749

#line 749
	Error = pthread_mutex_lock(&(bars->sl_phase_7).mutex);
#line 749
	if (Error != 0) {
#line 749
		printf("Error while trying to get lock in barrier.\n");
#line 749
		exit(-1);
#line 749
	}
#line 749

#line 749
	Cycle = (bars->sl_phase_7).cycle;
#line 749
	if (++(bars->sl_phase_7).counter != (nprocs)) {
#line 749
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 749
		while (Cycle == (bars->sl_phase_7).cycle) {
#line 749
			Error = pthread_cond_wait(&(bars->sl_phase_7).cv, &(bars->sl_phase_7).mutex);
#line 749
			if (Error != 0) {
#line 749
				break;
#line 749
			}
#line 749
		}
#line 749
		pthread_setcancelstate(Cancel, &Temp);
#line 749
	} else {
#line 749
		(bars->sl_phase_7).cycle = !(bars->sl_phase_7).cycle;
#line 749
		(bars->sl_phase_7).counter = 0;
#line 749
		Error = pthread_cond_broadcast(&(bars->sl_phase_7).cv);
#line 749
	}
#line 749
	pthread_mutex_unlock(&(bars->sl_phase_7).mutex);
#line 749
}
#else
   {
#line 751
	unsigned long	Error, Cycle;
#line 751
	int		Cancel, Temp;
#line 751

#line 751
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 751
	if (Error != 0) {
#line 751
		printf("Error while trying to get lock in barrier.\n");
#line 751
		exit(-1);
#line 751
	}
#line 751

#line 751
	Cycle = (bars->barrier).cycle;
#line 751
	if (++(bars->barrier).counter != (nprocs)) {
#line 751
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 751
		while (Cycle == (bars->barrier).cycle) {
#line 751
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 751
			if (Error != 0) {
#line 751
				break;
#line 751
			}
#line 751
		}
#line 751
		pthread_setcancelstate(Cancel, &Temp);
#line 751
	} else {
#line 751
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 751
		(bars->barrier).counter = 0;
#line 751
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 751
	}
#line 751
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 751
}
#endif
   zsim_stamp();
/*      *******************************************************

                e i g h t h   p h a s e

        *******************************************************

   augment ga(i,j) with [-psiai/psibi]*psib(i,j)

   %%%%%%%%%%%%%%% f4 should be private  */

   f4 = (-global->psiai)/(global->psibi);

   if (procid == MASTER) {
     wrk1->ga[0][0] = wrk1->ga[0][0]+f4*wrk1->psib[0][0];
   }
   if (procid == nprocs-xprocs) {
     wrk1->ga[im-1][0] = wrk1->ga[im-1][0]+f4*wrk1->psib[im-1][0];
   }
   if (procid == xprocs-1) {
     wrk1->ga[0][jm-1] = wrk1->ga[0][jm-1]+f4*wrk1->psib[0][jm-1];
   }
   if (procid == nprocs-1) {
     wrk1->ga[im-1][jm-1] = wrk1->ga[im-1][jm-1]+f4*wrk1->psib[im-1][jm-1];
   }
   if (firstrow == 1) {
     for(j=firstcol;j<=lastcol;j++) {
       wrk1->ga[0][j] = wrk1->ga[0][j]+f4*wrk1->psib[0][j];
     }
   }
   if ((firstrow+numrows) == im-1) {
     for(j=firstcol;j<=lastcol;j++) {
       wrk1->ga[im-1][j] = wrk1->ga[im-1][j]+f4*wrk1->psib[im-1][j];
     }
   }
   if (firstcol == 1) {
     for(j=firstrow;j<=lastrow;j++) {
       wrk1->ga[j][0] = wrk1->ga[j][0]+f4*wrk1->psib[j][0];
     }
   }
   if ((firstcol+numcols) == jm-1) {
     for(j=firstrow;j<=lastrow;j++) {
       wrk1->ga[j][jm-1] = wrk1->ga[j][jm-1]+f4*wrk1->psib[j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       wrk1->ga[i][iindex] = wrk1->ga[i][iindex]+f4*wrk1->psib[i][iindex];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 804
	unsigned long	Error, Cycle;
#line 804
	int		Cancel, Temp;
#line 804

#line 804
	Error = pthread_mutex_lock(&(bars->sl_phase_8).mutex);
#line 804
	if (Error != 0) {
#line 804
		printf("Error while trying to get lock in barrier.\n");
#line 804
		exit(-1);
#line 804
	}
#line 804

#line 804
	Cycle = (bars->sl_phase_8).cycle;
#line 804
	if (++(bars->sl_phase_8).counter != (nprocs)) {
#line 804
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 804
		while (Cycle == (bars->sl_phase_8).cycle) {
#line 804
			Error = pthread_cond_wait(&(bars->sl_phase_8).cv, &(bars->sl_phase_8).mutex);
#line 804
			if (Error != 0) {
#line 804
				break;
#line 804
			}
#line 804
		}
#line 804
		pthread_setcancelstate(Cancel, &Temp);
#line 804
	} else {
#line 804
		(bars->sl_phase_8).cycle = !(bars->sl_phase_8).cycle;
#line 804
		(bars->sl_phase_8).counter = 0;
#line 804
		Error = pthread_cond_broadcast(&(bars->sl_phase_8).cv);
#line 804
	}
#line 804
	pthread_mutex_unlock(&(bars->sl_phase_8).mutex);
#line 804
}
#else
   {
#line 806
	unsigned long	Error, Cycle;
#line 806
	int		Cancel, Temp;
#line 806

#line 806
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 806
	if (Error != 0) {
#line 806
		printf("Error while trying to get lock in barrier.\n");
#line 806
		exit(-1);
#line 806
	}
#line 806

#line 806
	Cycle = (bars->barrier).cycle;
#line 806
	if (++(bars->barrier).counter != (nprocs)) {
#line 806
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 806
		while (Cycle == (bars->barrier).cycle) {
#line 806
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 806
			if (Error != 0) {
#line 806
				break;
#line 806
			}
#line 806
		}
#line 806
		pthread_setcancelstate(Cancel, &Temp);
#line 806
	} else {
#line 806
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 806
		(bars->barrier).counter = 0;
#line 806
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 806
	}
#line 806
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 806
}
#endif
   for(i=istart;i<=iend;i++) {
     for(j=jstart;j<=jend;j++) {
       multi->rhs_multi[numlev-1][i][j] = wrk1->gb[i][j] * ressqr;
     }
   }
   if (istart == 0) {
     for(j=jstart;j<=jend;j++) {
       multi->q_multi[numlev-1][0][j] = wrk1->gb[0][j];
     }
   }
   if (iend == im-1) {
     for(j=jstart;j<=jend;j++) {
       multi->q_multi[numlev-1][im-1][j] = wrk1->gb[im-1][j];
     }
   }
   if (jstart == 0) {
     for(i=istart;i<=iend;i++) {
       multi->q_multi[numlev-1][i][0] = wrk1->gb[i][0];
     }
   }
   if (jend == jm-1) {
     for(i=istart;i<=iend;i++) {
       multi->q_multi[numlev-1][i][jm-1] = wrk1->gb[i][jm-1];
     }
   }

   fac = 1.0 / (4.0 - ressqr*eig2);
   for(i=ist;i<=ien;i++) {
     for(j=jst;j<=jen;j++) {
       multi->q_multi[numlev-1][i][j] = guess->oldgb[i][j];
     }
   }

   if ((procid == MASTER) || (do_stats)) {
     {
#line 842
	struct timeval	FullTime;
#line 842

#line 842
	gettimeofday(&FullTime, NULL);
#line 842
	(multi_start) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 842
};
   }

   multig(procid);

   if ((procid == MASTER) || (do_stats)) {
     {
#line 848
	struct timeval	FullTime;
#line 848

#line 848
	gettimeofday(&FullTime, NULL);
#line 848
	(multi_end) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 848
};
     gp[procid].multi_time += (multi_end - multi_start);
   }

   for(i=istart;i<=iend;i++) {
     for(j=jstart;j<=jend;j++) {
       wrk1->gb[i][j] = multi->q_multi[numlev-1][i][j];
       guess->oldgb[i][j] = multi->q_multi[numlev-1][i][j];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 859
	unsigned long	Error, Cycle;
#line 859
	int		Cancel, Temp;
#line 859

#line 859
	Error = pthread_mutex_lock(&(bars->sl_phase_8).mutex);
#line 859
	if (Error != 0) {
#line 859
		printf("Error while trying to get lock in barrier.\n");
#line 859
		exit(-1);
#line 859
	}
#line 859

#line 859
	Cycle = (bars->sl_phase_8).cycle;
#line 859
	if (++(bars->sl_phase_8).counter != (nprocs)) {
#line 859
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 859
		while (Cycle == (bars->sl_phase_8).cycle) {
#line 859
			Error = pthread_cond_wait(&(bars->sl_phase_8).cv, &(bars->sl_phase_8).mutex);
#line 859
			if (Error != 0) {
#line 859
				break;
#line 859
			}
#line 859
		}
#line 859
		pthread_setcancelstate(Cancel, &Temp);
#line 859
	} else {
#line 859
		(bars->sl_phase_8).cycle = !(bars->sl_phase_8).cycle;
#line 859
		(bars->sl_phase_8).counter = 0;
#line 859
		Error = pthread_cond_broadcast(&(bars->sl_phase_8).cv);
#line 859
	}
#line 859
	pthread_mutex_unlock(&(bars->sl_phase_8).mutex);
#line 859
}
#else
   {
#line 861
	unsigned long	Error, Cycle;
#line 861
	int		Cancel, Temp;
#line 861

#line 861
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 861
	if (Error != 0) {
#line 861
		printf("Error while trying to get lock in barrier.\n");
#line 861
		exit(-1);
#line 861
	}
#line 861

#line 861
	Cycle = (bars->barrier).cycle;
#line 861
	if (++(bars->barrier).counter != (nprocs)) {
#line 861
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 861
		while (Cycle == (bars->barrier).cycle) {
#line 861
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 861
			if (Error != 0) {
#line 861
				break;
#line 861
			}
#line 861
		}
#line 861
		pthread_setcancelstate(Cancel, &Temp);
#line 861
	} else {
#line 861
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 861
		(bars->barrier).counter = 0;
#line 861
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 861
	}
#line 861
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 861
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

   if (procid == MASTER) {
     wrk3->work2[0][0] = wrk1->gb[0][0]-hh1*wrk1->ga[0][0];
     wrk2->work3[0][0] = wrk1->gb[0][0]+hh3*wrk1->ga[0][0];
   }
   if (procid == nprocs-xprocs) {
     wrk3->work2[im-1][0] = wrk1->gb[im-1][0]-hh1*wrk1->ga[im-1][0];
     wrk2->work3[im-1][0] = wrk1->gb[im-1][0]+hh3*wrk1->ga[im-1][0];
   }
   if (procid == xprocs-1) {
     wrk3->work2[0][jm-1] = wrk1->gb[0][jm-1]-hh1*wrk1->ga[0][jm-1];
     wrk2->work3[0][jm-1] = wrk1->gb[0][jm-1]+hh3*wrk1->ga[0][jm-1];
   }
   if (procid == nprocs-1) {
     wrk3->work2[im-1][jm-1] = wrk1->gb[im-1][jm-1]-hh1*wrk1->ga[im-1][jm-1];
     wrk2->work3[im-1][jm-1] = wrk1->gb[im-1][jm-1]+hh3*wrk1->ga[im-1][jm-1];
   }
   if (firstrow == 1) {
     for(j=firstcol;j<=lastcol;j++) {
       wrk2->work3[0][j] = wrk1->gb[0][j]+hh3*wrk1->ga[0][j];
       wrk3->work2[0][j] = wrk1->gb[0][j]-hh1*wrk1->ga[0][j];
     }
   }
   if ((firstrow+numrows) == im-1) {
     for(j=firstcol;j<=lastcol;j++) {
       wrk2->work3[im-1][j] = wrk1->gb[im-1][j]+hh3*wrk1->ga[im-1][j];
       wrk3->work2[im-1][j] = wrk1->gb[im-1][j]-hh1*wrk1->ga[im-1][j];
     }
   }
   if (firstcol == 1) {
     for(j=firstrow;j<=lastrow;j++) {
       wrk2->work3[j][0] = wrk1->gb[j][0]+hh3*wrk1->ga[j][0];
       wrk3->work2[j][0] = wrk1->gb[j][0]-hh1*wrk1->ga[j][0];
     }
   }
   if ((firstcol+numcols) == jm-1) {
     for(j=firstrow;j<=lastrow;j++) {
       wrk2->work3[j][jm-1] = wrk1->gb[j][jm-1]+hh3*wrk1->ga[j][jm-1];
       wrk3->work2[j][jm-1] = wrk1->gb[j][jm-1]-hh1*wrk1->ga[j][jm-1];
     }
   }

   for(i=firstrow;i<=lastrow;i++) {
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
       wrk2->work3[i][iindex] = wrk1->gb[i][iindex]+hh3*wrk1->ga[i][iindex];
       wrk3->work2[i][iindex] = wrk1->gb[i][iindex]-hh1*wrk1->ga[i][iindex];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 927
	unsigned long	Error, Cycle;
#line 927
	int		Cancel, Temp;
#line 927

#line 927
	Error = pthread_mutex_lock(&(bars->sl_phase_9).mutex);
#line 927
	if (Error != 0) {
#line 927
		printf("Error while trying to get lock in barrier.\n");
#line 927
		exit(-1);
#line 927
	}
#line 927

#line 927
	Cycle = (bars->sl_phase_9).cycle;
#line 927
	if (++(bars->sl_phase_9).counter != (nprocs)) {
#line 927
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 927
		while (Cycle == (bars->sl_phase_9).cycle) {
#line 927
			Error = pthread_cond_wait(&(bars->sl_phase_9).cv, &(bars->sl_phase_9).mutex);
#line 927
			if (Error != 0) {
#line 927
				break;
#line 927
			}
#line 927
		}
#line 927
		pthread_setcancelstate(Cancel, &Temp);
#line 927
	} else {
#line 927
		(bars->sl_phase_9).cycle = !(bars->sl_phase_9).cycle;
#line 927
		(bars->sl_phase_9).counter = 0;
#line 927
		Error = pthread_cond_broadcast(&(bars->sl_phase_9).cv);
#line 927
	}
#line 927
	pthread_mutex_unlock(&(bars->sl_phase_9).mutex);
#line 927
}
#else
   {
#line 929
	unsigned long	Error, Cycle;
#line 929
	int		Cancel, Temp;
#line 929

#line 929
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 929
	if (Error != 0) {
#line 929
		printf("Error while trying to get lock in barrier.\n");
#line 929
		exit(-1);
#line 929
	}
#line 929

#line 929
	Cycle = (bars->barrier).cycle;
#line 929
	if (++(bars->barrier).counter != (nprocs)) {
#line 929
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 929
		while (Cycle == (bars->barrier).cycle) {
#line 929
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 929
			if (Error != 0) {
#line 929
				break;
#line 929
			}
#line 929
		}
#line 929
		pthread_setcancelstate(Cancel, &Temp);
#line 929
	} else {
#line 929
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 929
		(bars->barrier).counter = 0;
#line 929
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 929
	}
#line 929
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 929
}
#endif
   zsim_stamp();
/*      *******************************************************

                t e n t h    p h a s e

        *******************************************************/


   timst = 2*dtau;

/* update the psi{1,3} matrices by adding 2*dtau*work3 to each */

   if (procid == MASTER) {
     fields->psi[0][0][0] = fields->psi[0][0][0] + timst*wrk2->work3[0][0];
   }
   if (procid == nprocs-xprocs) {
     fields->psi[0][im-1][0] = fields->psi[0][im-1][0] + timst*wrk2->work3[im-1][0];
   }
   if (procid == xprocs-1) {
     fields->psi[0][0][jm-1] = fields->psi[0][0][jm-1] + timst*wrk2->work3[0][jm-1];
   }
   if (procid == nprocs-1) {
     fields->psi[0][im-1][jm-1] = fields->psi[0][im-1][jm-1] + timst*wrk2->work3[im-1][jm-1];
   }
   if (firstrow == 1) {
     for(j=firstcol;j<=lastcol;j++) {
       fields->psi[0][0][j] = fields->psi[0][0][j] + timst*wrk2->work3[0][j];
     }
   }
   if ((firstrow+numrows) == im-1) {
     for(j=firstcol;j<=lastcol;j++) {
       fields->psi[0][im-1][j] = fields->psi[0][im-1][j] + timst*wrk2->work3[im-1][j];
     }
   }
   if (firstcol == 1) {
     for(j=firstrow;j<=lastrow;j++) {
       fields->psi[0][j][0] = fields->psi[0][j][0] + timst*wrk2->work3[j][0];
     }
   }
   if ((firstcol+numcols) == jm-1) {
     for(j=firstrow;j<=lastrow;j++) {
       fields->psi[0][j][jm-1] = fields->psi[0][j][jm-1] + timst*wrk2->work3[j][jm-1];
     }
   }
   for(i=firstrow;i<=lastrow;i++) {
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
         fields->psi[0][i][iindex] = fields->psi[0][i][iindex] + timst*wrk2->work3[i][iindex];
     }
   }

   if (procid == MASTER) {
     fields->psi[1][0][0] = fields->psi[1][0][0] + timst*wrk3->work2[0][0];
   }
   if (procid == nprocs-xprocs) {
     fields->psi[1][im-1][0] = fields->psi[1][im-1][0] + timst*wrk3->work2[im-1][0];
   }
   if (procid == xprocs-1) {
     fields->psi[1][0][jm-1] = fields->psi[1][0][jm-1] + timst*wrk3->work2[0][jm-1];
   }
   if (procid == nprocs-1) {
     fields->psi[1][im-1][jm-1] = fields->psi[1][im-1][jm-1] + timst*wrk3->work2[im-1][jm-1];
   }
   if (firstrow == 1) {
     for(j=firstcol;j<=lastcol;j++) {
       fields->psi[1][0][j] = fields->psi[1][0][j] + timst*wrk3->work2[0][j];
     }
   }
   if ((firstrow+numrows) == im-1) {
     for(j=firstcol;j<=lastcol;j++) {
       fields->psi[1][im-1][j] = fields->psi[1][im-1][j] + timst*wrk3->work2[im-1][j];
     }
   }
   if (firstcol == 1) {
     for(j=firstrow;j<=lastrow;j++) {
       fields->psi[1][j][0] = fields->psi[1][j][0] + timst*wrk3->work2[j][0];
     }
   }
   if ((firstcol+numcols) == jm-1) {
     for(j=firstrow;j<=lastrow;j++) {
       fields->psi[1][j][jm-1] = fields->psi[1][j][jm-1] + timst*wrk3->work2[j][jm-1];
     }
   }

   for(i=firstrow;i<=lastrow;i++) {
     for(iindex=firstcol;iindex<=lastcol;iindex++) {
         fields->psi[1][i][iindex] = fields->psi[1][i][iindex] + timst*wrk3->work2[i][iindex];
     }
   }
#if defined(MULTIPLE_BARRIERS)
   {
#line 1020
	unsigned long	Error, Cycle;
#line 1020
	int		Cancel, Temp;
#line 1020

#line 1020
	Error = pthread_mutex_lock(&(bars->sl_phase_10).mutex);
#line 1020
	if (Error != 0) {
#line 1020
		printf("Error while trying to get lock in barrier.\n");
#line 1020
		exit(-1);
#line 1020
	}
#line 1020

#line 1020
	Cycle = (bars->sl_phase_10).cycle;
#line 1020
	if (++(bars->sl_phase_10).counter != (nprocs)) {
#line 1020
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 1020
		while (Cycle == (bars->sl_phase_10).cycle) {
#line 1020
			Error = pthread_cond_wait(&(bars->sl_phase_10).cv, &(bars->sl_phase_10).mutex);
#line 1020
			if (Error != 0) {
#line 1020
				break;
#line 1020
			}
#line 1020
		}
#line 1020
		pthread_setcancelstate(Cancel, &Temp);
#line 1020
	} else {
#line 1020
		(bars->sl_phase_10).cycle = !(bars->sl_phase_10).cycle;
#line 1020
		(bars->sl_phase_10).counter = 0;
#line 1020
		Error = pthread_cond_broadcast(&(bars->sl_phase_10).cv);
#line 1020
	}
#line 1020
	pthread_mutex_unlock(&(bars->sl_phase_10).mutex);
#line 1020
}
#else
   {
#line 1022
	unsigned long	Error, Cycle;
#line 1022
	int		Cancel, Temp;
#line 1022

#line 1022
	Error = pthread_mutex_lock(&(bars->barrier).mutex);
#line 1022
	if (Error != 0) {
#line 1022
		printf("Error while trying to get lock in barrier.\n");
#line 1022
		exit(-1);
#line 1022
	}
#line 1022

#line 1022
	Cycle = (bars->barrier).cycle;
#line 1022
	if (++(bars->barrier).counter != (nprocs)) {
#line 1022
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 1022
		while (Cycle == (bars->barrier).cycle) {
#line 1022
			Error = pthread_cond_wait(&(bars->barrier).cv, &(bars->barrier).mutex);
#line 1022
			if (Error != 0) {
#line 1022
				break;
#line 1022
			}
#line 1022
		}
#line 1022
		pthread_setcancelstate(Cancel, &Temp);
#line 1022
	} else {
#line 1022
		(bars->barrier).cycle = !(bars->barrier).cycle;
#line 1022
		(bars->barrier).counter = 0;
#line 1022
		Error = pthread_cond_broadcast(&(bars->barrier).cv);
#line 1022
	}
#line 1022
	pthread_mutex_unlock(&(bars->barrier).mutex);
#line 1022
}
#endif
//#ifdef ZSIM_TRACE_4
  //      zsim_PIM_function_end();
  //  zsim_roi_end();
//#endif

   zsim_stamp();

}
