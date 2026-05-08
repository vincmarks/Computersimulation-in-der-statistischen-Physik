/***********************************************
 * pressure.c: measures pressure in the system *
 *                                             *
 * last modification: 23.05.2003               *
 ***********************************************/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "element.h"


/* global variables */

extern int        nrchainsA, nrmonA, polA;
extern int        nrchainsB, nrmonB, polB;
extern double     LSX,LSY,LSZ;
extern MYVEC      *posA,*poslatA;
extern MYVEC      *posB,*poslatB;
extern double     T,Ewall,Fwall;
extern DATA       data;


/*****************************************
 * calc_p: calculates pressure using the *
 *         virial expression             *
 *                                       *
 * P  = (nkT)/V + 1/(3V) *               *
 *      sum (Fx*rx + Fy*ry + Fz*rz)      *
 *                                       *
 * Fx = -dE/dx, Fy = -dE/dy, Fz = -dE/dZ *
 *****************************************/

int calc_p(void)
{
int     i1,ii1,i2,ii2,type1,type1a,type2,type2a;
double  V,dummy,rad2,invrad2,invrad4,invrad8,invrad14;
double  p0,pljx,pljy,pljz,pfex,pfey,pfez,pwall1,pwall2,z1,z2;
double  z1_10,z1_4,z2_10,z2_4,px,py,pz,virialx,virialy,virialz;
MYVEC   *mon1,*mon2,diff;

   
   /*--------------------*
    | set variables to 0 |
    *--------------------*/  

   V      = LSX*LSY*LSZ;
   pljx   = 0;
   pljy   = 0;
   pljz   = 0;
   pfex   = 0;
   pfey   = 0;
   pfez   = 0;
   pwall1 = 0;
   pwall2 = 0;
   virialx =0;
   virialy =0;
   virialz =0;

#ifdef HARDWALLS
   for( i1=0,mon1=pos ; i1<nrchains ; i1++){
     for ( ii1=0 ; ii1< pol ; ii1++){
       z1    = 1.0/mon1->z;
       z2    = 1.0/(LSZ-mon1->z);
       z1_4  = z1*z1*z1*z1;
       z2_4  = z2*z2*z2*z2;
       z1_10 = z1_4*z1_4*z1*z1;
       z2_10 = z2_4*z2_4*z2*z2;
       dummy = 3*Ewall/(LSX*LSY)*
               (3*z1_10-Fwall*z1_4-3*z2_10+Fwall*z2_4);
       if (mon1->z<LSZ/2.0)
         pwall1 += dummy;
       else
         pwall1 -= dummy;                     
         /* modified pwall1 is the pressure the wall */
         /* exerts on the system                     */

       pwall2 += (mon1->z-0.5*LSZ)*dummy/LSZ; 
       /* pwall2 is the viral of the external forces */
       /* (check again)*/ 

       mon1++;
     }
   }
#endif

   
   /*----------------------------------*
    | step 1: consider AA interactions |
    *----------------------------------*/

   for (i1=0,mon1=poslatA;i1<nrchainsA;i1++)  
   /* loop over all polymers */

   for (ii1=0;ii1<polA;ii1++,mon1++)          
   /* loop over all monomers */

   for (i2=0,mon2=poslatA;i2<nrchainsA;i2++)  
   /* loop over all interaction partner polymers */

   for (ii2=0;ii2<polA;ii2++,mon2++) {        
   /* loop over all interaction partner monomers */
   /* -> all monomers interact with each other   */

     if (i1*polA+ii1>i2*polA+ii2) {           
     /* each interaction is only counted once      */

       /* determine distance between interacting monomers */ 
       diff.x = mon1->x - mon2->x;
       diff.y = mon1->y - mon2->y;                  
       diff.z = mon1->z - mon2->z; 

       /* Minimum Image Convention */
       diff.x = diff.x - LSX*((int)(2.0*diff.x/LSX));         
       diff.y = diff.y - LSY*((int)(2.0*diff.y/LSY));

#ifndef HARDWALLS
       diff.z = diff.z - LSZ*((int)(2.0*diff.z/LSZ));
#endif

        
       /*---------------------------*
        | compute LJ part of virial |
        *---------------------------*/

       rad2 = SQUARE(diff.x)+SQUARE(diff.y)+SQUARE(diff.z);
      
       if (rad2 < RANGE2AA) {                                  
          invrad2  = 1.0/rad2;    
          invrad4  = invrad2*invrad2;
          invrad8  = invrad4*invrad4;
          invrad14 = invrad8*invrad4*invrad2;  
          p0       = 24*(2*invrad14-invrad8);    
          /* compare with LJ potential definition  */
          /* in element.h EPSILONA and SIGMAA := 1 */

 
          pljx    += diff.x*diff.x*p0; /* = -(dE(LJ)/dx)*x */
          pljy    += diff.y*diff.y*p0; /* = -(dE(LJ)/dy)*y */
          pljz    += diff.z*diff.z*p0; /* = -(dE(LJ)/dz)*z */
       }

       
       /*-----------------------------*
        | compute FENE part of virial |          
        *-----------------------------*/

       if ((i1==i2)&&(ii1==ii2+1)) {                           
          p0       = -67.5*RNULL2INV/(1-rad2*RNULL2INV);  
          /* compare with FENE potential definition */ 
          /* in element.h SIGMAA in RNULL2INV :=1   */

          pfex    += diff.x*diff.x*p0; /* = -(dE(FENE)/dx)*x */
          pfey    += diff.y*diff.y*p0; /* = -(dE(FENE)/dy)*y */
          pfez    += diff.z*diff.z*p0; /* = -(dE(FENE)/dz)*z */
       }
     } /* each interaction is only counted once */
   } /* all AA interactions are considered */


   /*----------------------------------*
    | step 2: consider BB interactions |
    *----------------------------------*/

   for (i1=0,mon1=poslatB;i1<nrchainsB;i1++)  
   /* loop over all polymers */

   for (ii1=0;ii1<polB;ii1++,mon1++)          
   /* loop over all monomers */

   for (i2=0,mon2=poslatB;i2<nrchainsB;i2++)  
   /* loop over all interaction partner polymers */

   for (ii2=0;ii2<polB;ii2++,mon2++)  {       
   /* loop over all interaction partner monomers */
   /* -> all monomers interact with each other   */

     if (i1*polB+ii1>i2*polB+ii2) {
     /* each interaction is only counted once      */

       /* determine distance between interacting monomers */
       diff.x = mon1->x - mon2->x;
       diff.y = mon1->y - mon2->y;                   
       diff.z = mon1->z - mon2->z; 

       /* Minimum Image Convention */
       diff.x = diff.x - LSX*((int)(2.0*diff.x/LSX));         
       diff.y = diff.y - LSY*((int)(2.0*diff.y/LSY));

#ifndef HARDWALLS
       diff.z = diff.z - LSZ*((int)(2.0*diff.z/LSZ));
#endif
        
       /*---------------------------*
        | compute LJ part of virial |
        *---------------------------*/

       rad2 = SQUARE(diff.x)+SQUARE(diff.y)+SQUARE(diff.z);
      
       if (rad2 < RANGE2BB) {                                  
          invrad2  = 1.0/rad2;    
          invrad4  = invrad2*invrad2;
          invrad8  = invrad4*invrad4;
          invrad14 = invrad8*invrad4*invrad2;  
          p0       = 24*EPSILONB*(2*SQUARE(SIGMAB6)*invrad14-
                     SIGMAB6*invrad8);
          /* compare with LJ potential definition */
          /* in element.h                         */
 
          pljx    += diff.x*diff.x*p0;  /* = -(dE(LJ)/dx)*x */
          pljy    += diff.y*diff.y*p0;  /* = -(dE(LJ)/dy)*y */
          pljz    += diff.z*diff.z*p0;  /* = -(dE(LJ)/dz)*z */ 
       }

       
       /*-----------------------------*
        | compute FENE part of virial |          
        *-----------------------------*/

       if ((i1==i2)&&(ii1==ii2+1)) { 
       /* has to be adjusted to 2nd polymer species ! */
       /* -> change RNULL2INV to RNULL2INVB, EPSILONB */       
       /* right now it only works for N=1 */
    
          p0       = -67.5*RNULL2INV/(1-rad2*RNULL2INV);  
          /* compare with FENE potential definition */
          /* in element.h                           */
 
          pfex    += diff.x*diff.x*p0; /* = -(dE(FENE)/dx)*x */
          pfey    += diff.y*diff.y*p0; /* = -(dE(FENE)/dy)*y */
          pfez    += diff.z*diff.z*p0; /* = -(dE(FENE)/dz)*z */
       }
     } /* each interaction is only counted once */
   } /* all BB interactions are considered */


   /*----------------------------------*
    | step 3: consider AB interactions |
    *----------------------------------*/

   for (i1=0,mon1=poslatA;i1<nrchainsA;i1++)  
   /* loop over all polymers */

   for (ii1=0;ii1<polA;ii1++,mon1++)          
   /* loop over all monomers */

   for (i2=0,mon2=poslatB;i2<nrchainsB;i2++)  
   /* loop over all interaction partner polymers */

   for (ii2=0;ii2<polB;ii2++,mon2++) {        
   /* loop over all interaction partner monomers */
   /* -> all monomers interact with each other   */

        /* determine distance between interacting monomers */
       diff.x = mon1->x - mon2->x;
       diff.y = mon1->y - mon2->y;                  
       diff.z = mon1->z - mon2->z; 

       /* Minimum Image Convention */
       diff.x = diff.x - LSX*((int)(2.0*diff.x/LSX));         
       diff.y = diff.y - LSY*((int)(2.0*diff.y/LSY));

#ifndef HARDWALLS
       diff.z = diff.z - LSZ*((int)(2.0*diff.z/LSZ));
#endif

        
       /*---------------------------*
        | compute LJ part of virial |
        *---------------------------*/

       rad2 = SQUARE(diff.x)+SQUARE(diff.y)+SQUARE(diff.z);
      
       if (rad2 < RANGE2AB) {      
          invrad2  = 1.0/rad2;    
          invrad4  = invrad2*invrad2;
          invrad8  = invrad4*invrad4;
          invrad14 = invrad8*invrad4*invrad2;  
          p0       = 24*EPSILONAB*(2*SQUARE(SIGMAAB6)*invrad14-
                     SIGMAAB6*invrad8);   
          /* compare with LJ potential definition */
          /* in element.h                         */
 
          pljx    += diff.x*diff.x*p0;  /* = -(dE(LJ)/dx)*x */
          pljy    += diff.y*diff.y*p0;  /* = -(dE(LJ)/dy)*y */
          pljz    += diff.z*diff.z*p0;  /* = -(dE(LJ)/dz)*z */ 
       
          virialx += diff.x*diff.x*p0;
          virialy += diff.y*diff.y*p0;
          virialz += diff.z*diff.z*p0;
     }      
   } /* all AB interactions are considered */


   /*------------------------------------------*
    | step 4: calculate pressure in the system |
    *------------------------------------------*/

   px = (T*(double)(nrmonA+nrmonB) + pljx+pfex)/V;  
   /* px, py and pz each account for */

   py = (T*(double)(nrmonA+nrmonB) + pljy+pfey)/V;  
   /* total pressure in the system   */

   pz = (T*(double)(nrmonA+nrmonB) + pljz+pfez)/V;  
   /* - in a homogeneous system they */
   /* should yield the same result   */

   data.px     = px;  /* save results in data-array     */
   data.py     = py;
   data.pz     = pz;
   data.virialx = virialx/V;
   data.virialy = virialy/V;
   data.virialz = virialz/V;
   data.pwall1 = pwall1;
   data.pwall2 = pwall2;

return(OK);
}
