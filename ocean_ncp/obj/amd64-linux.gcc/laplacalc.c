#line 233 "/mnt/panzer/mnika/Workload_Pool/parsec-3.0/bin/../pkgs/libs/parmacs/inst/amd64-linux.gcc/m4/parmacs.pthreads.c.m4"

#line 1 "laplacalc.C"
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

/*      **************************************************************
                       end of subroutine jacobcalc
        **************************************************************

   performs the laplacian calculation for a subblock. */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "decs.h"
#if defined(ZSIM_TRACE_1) || defined(ZSIM_TRACE_2) || defined(ZSIM_TRACE_3) || defined(ZSIM_TRACE_4)
 #include "zsim_hooks.h"
#endif
int roi_included_1 = 0;
void laplacalc(double x[IMAX][JMAX], double z[IMAX][JMAX], long firstrow, long lastrow, long firstcol, long lastcol, long numrows, long numcols)
{
//#ifdef ZSIM_TRACE_3
  //         zsim_roi_begin();
//	zsim_PIM_function_begin();
//#endif
   long iindex;
   long indexp1;
   long indexm1;
   long ip1;
   long im1;
   long i;
   long j;

   for (iindex=firstcol;iindex<=lastcol;iindex++) {
     indexp1 = iindex+1;
     indexm1 = iindex-1;
     for (i=firstrow;i<=lastrow;i++) {
       ip1 = i+1;
       im1 = i-1;
       z[i][iindex] = factlap*(x[ip1][iindex]+x[im1][iindex]+x[i][indexp1]+
                              x[i][indexm1]-4.*x[i][iindex]);
     }
   }
   if (firstrow == 1) {
     for (j=firstcol;j<=lastcol;j++) {
       z[0][j] = 0.0;
     }
   }
   if ((firstrow+numrows) == im-1) {
     for (j=firstcol;j<=lastcol;j++) {
       z[im-1][j] = 0.0;
     }
   }
   if (firstcol == 1) {
     for (j=firstrow;j<=lastrow;j++) {
       z[j][0] = 0.0;
     }
   }
   if ((firstcol+numcols) == jm-1) {
     for (j=firstrow;j<=lastrow;j++) {
       z[j][jm-1] = 0.0;
     }
   }
//#ifdef ZSIM_TRACE_3
  //      zsim_PIM_function_end();
  //  zsim_roi_end();
//#endif

}
