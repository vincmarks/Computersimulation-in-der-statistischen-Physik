/*******************************************************
 * sysin.c: reads simulation specifications            *
 *          and monomer positions from file "conf_out" *
 *                                                     *
 * last modification: 23/05/2003                       *
 *******************************************************/


#include "element.h"
#include <stdio.h>
#include <stdlib.h>

extern int       nrchainsA, nrmonA, polA; 
extern int       nrchainsB, nrmonB, polB, mcs_start;
extern MYVEC     *posA, *poslatA,*posB, *poslatB;
extern double    T,CPA,CPB,Ewall,Fwall,J;
extern double    LSX,LSY,LSZ,LSSUBX,LSSUBY,LSSUBZ;
extern MYVEC     *cm0A,*cm0B;

#ifdef STRIPE
extern double    SEwall,SWidth;
#endif

int sysin(char *name)
{
int    i,ii;
MYVEC  *mon,*mon2;
FILE   *in;
double dummy;

   in = fopen( name,"r");


   /*------------------------------------------------*
    | read simulation specifications from "conf_out" |
    *------------------------------------------------*/ 

   fscanf(in,"T=%lf\n",&T);
   fscanf(in,"CPA=%lf\n",&CPA);
   fscanf(in,"CPB=%lf\n",&CPB);
   fscanf(in,"EW=%lf\n",&Ewall);
   fscanf(in,"FW=%lf\n",&Fwall);

#ifdef STRIPE
   fscanf(in,"SEW=%lf\n",&SEwall);
   fscanf(in,"SW=%lf\n",&SWidth);
#endif

   fscanf(in,"J=%lf\n",&J);
   fscanf(in,"LS=%lf %lf %lf\n",&LSX,&LSY,&LSZ);
   fscanf(in,"t=%d\n\n",&mcs_start);
   fscanf(in,"%d %d", &nrchainsA, &nrmonA);
   fscanf(in,"%d %d\n\n", &nrchainsB, &nrmonB);

   LSSUBX = LSX/((double)(NRSUBX));
   LSSUBY = LSY/((double)(NRSUBY));
   LSSUBZ = LSZ/((double)(NRSUBZ));


   /*-----------------*
    | allocate memory |
    *-----------------*/

   posA           = (MYVEC *) calloc(NMONOMAXA, sizeof(MYVEC));
   poslatA        = (MYVEC *) calloc(NMONOMAXA, sizeof(MYVEC));
   cm0A           = (MYVEC *) calloc(NPOLYMAXA, sizeof(MYVEC));
   posB           = (MYVEC *) calloc(NMONOMAXB, sizeof(MYVEC));
   poslatB        = (MYVEC *) calloc(NMONOMAXB, sizeof(MYVEC));
   cm0B           = (MYVEC *) calloc(NPOLYMAXB, sizeof(MYVEC));

   if (posA==NULL || poslatA==NULL || cm0A==NULL || 
       posB==NULL || poslatB==NULL ||cm0B==NULL) {
      fprintf(stderr,"ERROR: memory allocation denied"
                     "in function sysin\n");
      exit(ERROR);
   }

   if (nrchainsA==0) {
      polA = NMONOMAXA/NPOLYMAXA;
   } 
   if (nrchainsB==0) {
      polB = NMONOMAXB/NPOLYMAXB;
   } 


   /*--------------------------------------------------*
    | read typeA monomer positions from "conf_out" and |
    | update posA array                                |
    *--------------------------------------------------*/

   for( i=0,mon=posA ; i<nrchainsA ; i++){
      fscanf( in,"%d", &(polA));
      for ( ii=0 ; ii< polA ; ii++){
         fscanf( in,"%lf%lf%lf",&(mon->x),&(mon->y),&(mon->z));
         fscanf( in,"%lf%lf%lf", &dummy, &dummy, &dummy);  
         /* space for velocities -                  */
         /* makes conf_out compatible to MD version */
         mon ++; 
      }
   }


   /*-------------------------------*
    | update poslatA and cm0A array |
    *-------------------------------*/

   for(i=0;i<nrchainsA;i++) {
     mon  = posA    + i*polA;
     mon2 = poslatA + i*polA;

     calc_cm(polA,mon,(cm0A+i));  /* see wach.c */
  
     for (ii=0;ii<polA;ii++,mon++,mon2++) {
      mon2->x = mon->x-LSX*((int)(mon->x/LSX));
      mon2->y = mon->y-LSY*((int)(mon->y/LSY));
      mon2->z = mon->z-LSZ*((int)(mon->z/LSZ));

      if ( mon2->x < 0.0 ) mon2->x += LSX ;
      if ( mon2->y < 0.0 ) mon2->y += LSY ;
      if ( mon2->z < 0.0 ) mon2->z += LSZ ;
     }
   }

   
   /*--------------------------------------------------* 
    | read typeB monomer positions from "conf_out" and | 
    | update posA array                                |
    *--------------------------------------------------*/

   for( i=0,mon=posB ; i<nrchainsB ; i++){
      fscanf( in,"%d", &(polB));
      for ( ii=0 ; ii< polB ; ii++){
         fscanf( in,"%lf%lf%lf",&(mon->x),&(mon->y),&(mon->z));
         fscanf( in,"%lf%lf%lf", &dummy, &dummy, &dummy);
         mon ++;
      }
   }
   fclose( in);
   

   /*-------------------------------*
    | update poslatB and cm0B array |
    *-------------------------------*/   

   for(i=0;i<nrchainsB;i++) {
     mon  = posB    + i*polB;
     mon2 = poslatB + i*polB;

     calc_cm(polB,mon,(cm0B+i));

     for (ii=0;ii<polB;ii++,mon++,mon2++) {
      mon2->x = mon->x-LSX*((int)(mon->x/LSX));
      mon2->y = mon->y-LSY*((int)(mon->y/LSY));
      mon2->z = mon->z-LSZ*((int)(mon->z/LSZ));

      if ( mon2->x < 0.0 ) mon2->x += LSX ;
      if ( mon2->y < 0.0 ) mon2->y += LSY ;
      if ( mon2->z < 0.0 ) mon2->z += LSZ ;
     }
   }

   return(0);

}/* end sysin */
