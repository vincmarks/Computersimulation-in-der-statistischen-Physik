/*******************************************************
 * sysout.c: writes sim. specifications and monomer    *
 *           positions for next run to file "conf_out" *
 *                                                     *
 * last modification: 23/05/2003                       *
 *******************************************************/


#include "element.h"
#include <stdio.h>

extern double    T,CPA,CPB,Ewall,Fwall,J,LSX,LSY,LSZ;
extern int       nrchainsA, nrmonA, polA;
extern int       nrchainsB, nrmonB, polB;
extern MYVEC     *posA,*posB;

#ifdef STRIPE
extern double    SEwall,SWidth;
#endif

int sysout(char *name,int mcs)
{
   int i, ii;
   MYVEC *mon;
   FILE *out;


   /*------------------------------------------------------*
    | write sim. specifications for next run to "conf_out" |
    *------------------------------------------------------*/

   out = fopen(name,"w");

#ifdef STRIPE
   fprintf(out,"T=%g\nCPA=%g\nCPB=%g\nEW=%g\nFW=%g\n"
               "SEW=%g\nSW=%g\nJ=%g\nLS=%g %g %g\nt=%d\n\n"
               "%d %d\n%d %d\n\n",T,CPA,CPB,Ewall,Fwall,\
               SEwall,SWidth,J,LSX,LSY,LSZ,mcs,nrchainsA,\
               nrmonA,nrchainsB,nrmonB);
#else
   fprintf(out,"T=%g\nCPA=%g\nCPB=%g\nEW=%g\nFW=%g\n"
               "J=%g\nLS=%g %g %g\nt=%d\n\n%d %d\n%d %d\n\n",\
               T,CPA,CPB,Ewall,Fwall,J,LSX,LSY,LSZ,mcs,\
               nrchainsA,nrmonA,nrchainsB,nrmonB);
#endif


   /*------------------------------------------------*
    | write type A - monomer positions to "conf_out" |
    *------------------------------------------------*/

   for( i=0,mon=posA ; i<nrchainsA ; i++){
      fprintf( out,"%d\n", polA);
      for( ii=0 ; ii<polA ; ii++){
         fprintf( out,"%f\t%f\t%f\n", mon->x, mon->y, mon->z);
         fprintf( out,"0.000 0.000 0.000\n"); 
         mon++ ; 
      }
      fprintf( out,"\n");
   }

   /*------------------------------------------------*
    | write type B - monomer positions to "conf_out" |
    *------------------------------------------------*/

   for( i=0,mon=posB ; i<nrchainsB ; i++){
      fprintf( out,"%d\n", polB);
      for( ii=0 ; ii<polB ; ii++){
         fprintf( out,"%f\t%f\t%f\n", mon->x, mon->y, mon->z);
         fprintf( out,"0.000 0.000 0.000\n"); 
         mon++ ; 
      }
      fprintf( out,"\n");
   }

   fclose(out);
   return(OK);
}
