/************************************************************
 * weight.c: initializes weights for weighted simulation    * 
 *                                                          *
 * last modification: 23.05.2003                            *
 ************************************************************/ 

# include "element.h" 
# include <stdlib.h>
# include <stdio.h>
# include <math.h>


/************************************************************
 * weight.c: initializes array w[i][j]= log (e.g. of P(i,j) *  
 *           i = number of polymersA, j=number of polymersB)*
 *                                                          *
 *           'w' exists in working directory ?              *
 *           no:  -> simulation will not be weighted,       *
 *                   all w[i][j]=0,                         *
 *           yes: -> simulation will be weighted,           *
 *                   w[i][j] according to 'w'-file.         *
 *                                                          *
 * returns:  pointer to 'w'-array                           *
 *                                                          * 
 * idea:     w[i][j] is a guess of (the log of the) prob.   *
 *           distribution P(#polyA,#polyB). In cbgc.c       *
 *           (grand-canonical MC step) the Hamiltonian is   *
 *           modified such that in the ideal case (guess =  *
 *           real distribution) P=constant for all          *
 *           (NPOLYMAXA,NPOLYMAXB).                         *
 *           This enables simulations below the critical    *
 *           point.                                         *
 *                                                          *
 * remark:   One can obtain a 'w'-file by extrapolation of  *
 *           an existing data-set.                          *
 ************************************************************/ 


void init_weight(void)
{
extern double w[NPOLYMAXA+1][NPOLYMAXB+1];
extern double H[NPOLYMAXA+1][NPOLYMAXB+1];
extern double g;
extern int    adjust_reweight;
extern int    iteration_step;
extern long   _MPP_MY_PE;

FILE *dp=NULL; /* file-pointer to w-file */
char w_file_name[30];
int ctr1, ctr2;
int x,y;


  /*---------------------------------------*
   | initialize 'w'-array: set all w[][]=0 |
   *---------------------------------------*/

  for (ctr1=0;ctr1<=NPOLYMAXA;ctr1++)
  { 
    for (ctr2=0;ctr2<=NPOLYMAXB;ctr2++)
      { w[ctr1][ctr2]=0;
      }   
  }  
   
 
  /*---------------------------------*
   | simulation will not be weighted |
   | -> all w[][] remain equal to 0  |
   *---------------------------------*/
  
#ifdef T3D  
  sprintf(w_file_name,"w_%d",_MPP_MY_PE);
  if ((dp =fopen(w_file_name,"r"))==NULL)  
    /* 'w'-file does not exist -> simulation will not be */
    /* weighted                                          */

#else
  if ((dp =fopen("w","r"))==NULL)  
#endif
 
    {
      fprintf(stderr,"No weights available !\n"
                     "Continue with unweighted simulation\n");
    }


  /*---------------------------------------------------*
   | simulation will be weighted                       |
   | -> overwrite 'w'-array with weights from 'w'-file |
   *---------------------------------------------------*/

  else
    {
      for(ctr1=0;ctr1<=NPOLYMAXA;ctr1++)
        for (ctr2=0;ctr2<=NPOLYMAXB;ctr2++)
          { 
            fscanf(dp,"%d %d",&x,&y);
            fscanf(dp,"%lg\n",&w[x][y]);
          }
     }
}
