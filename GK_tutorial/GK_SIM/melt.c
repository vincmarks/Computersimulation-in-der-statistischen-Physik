/***********************************************
 * melt.c: contains main()                     *
 *                                             *
 * last modification : 23/05/2003              *
 ***********************************************/


/* include header files */
#include "element.h"
#include "r250.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>


/* global variables */
int         nrchainsA, nrmonA, polA;
int         nrchainsB, nrmonB, polB;
MYVEC       *posA, *poslatA;
MYVEC       *posB, *poslatB;

double      T, Tinv,J,CPA,CPB,Ewall,Fwall;
double      LSX,LSY,LSZ,LSSUBX,LSSUBY,LSSUBZ;

#ifdef STRIPE
double      SEwall,SWidth=0;
#endif

int         mcs_start;

MYVEC       *cm0A,*cm0B,pressure;        
/* used in wach.c */

double      cenergie=0,e_start;          
/* required for a check of the total energy */ 

double      mx=0,my=0,mz=0,manzahl=0;    
/* used in mcmove.c */

ELEMENT     *element_liste;              
/* used in element.c, local.c */

BBOX        boxliste[NRSUB3];
double      pcenergie;                   

DATA        data;                        
/* used in wach.c to transfer results to main */

double      Jwallwtab[JwallN_MAX]  ;     
/* used in switch_Jwall.c, keep global for i/o */

double      w[NPOLYMAXA+1][NPOLYMAXB+1]; 
/* array to store weighting factors */ 


int main(int argc, char *argv[]) {

   int     outstep, countstep, seed, j;
   int     mcs, mcsmax, nreptation, local, ncbgc=0, placc;
   int     cbgc_acc0=0,cbgc_acc1=0,cbgc_count0=1,cbgc_count1=1;
   int     lakzeptiert=0,lgesamt=1,rep_acc=0,rep_count=1,\
           out_ctr=0;   
   int     pol_print_ctr=360;
   double  e_end;
   double  zbondA,zbondB;
   double  config_ctr=0;
   char    name[80], filenameout[30]="conf_out",dname[30],\
           out_name[30];
   FILE    *eingabe,*log_file;
   time_t  now;

   extern long _MPP_MY_PE, _MPP_N_PES;  
   /* only used on T3E machines - number of processors */

   /* test variables */
   double energy1, energy2, lj_energy3=0, lj_energy4=0;
   int i=0,ctr1,ctr2;
   static double *rndtabA=NULL, *rnd_ptrA=NULL;
   int a [2000];
   int dummy1;
   FILE *prob;

   
   /*-----------------*
    | read parameters |
    *-----------------*/

   time(&now);
   PRINTFS("program starts at %s\n",ctime(&now));
      
   if (argc>1) {


     /*--------------------------------------------------*
      | a) for use on a regular single processor machine |
      *--------------------------------------------------*/

#ifndef T3D
     sprintf(name,"%s",argv[1]);              /* 'eingabe' */
     if ((eingabe=fopen(name,"r"))==NULL) 
       {
         fprintf(stderr,"couldn't find %s \n",name);
         exit(1);
       }

     /* read parameters from 'eingabe' file */

     fscanf(eingabe,"Steps=%d \n",&mcsmax);   
     fscanf(eingabe,"out=%d \n",&outstep);
     fscanf(eingabe,"seed=%d \n",&seed);
     fscanf(eingabe,"datei=%s \n",dname);     /* 'conf_out' */
     fscanf(eingabe,"reptation=%d\n",&nreptation);
     fscanf(eingabe,"local=%d\n",&local);
     fscanf(eingabe,"cbgc=%d\n",&ncbgc);
     fclose(eingabe);

     if(seed<0){
              time_t seconds = time(NULL);
              seed = (int) seconds;
     }  
     /* get seed for rnd number generator r250 from cpu-time */
     r250_srandom(seed);                 
     /* assign seed to rnd number generator */
         
     sysin(dname);           
     /* read sim. specifications and monomer vectors */
     /* from file "conf_out" */

     init_weight();          
     /* initializes weights for reweighted simulation */ 
   
     /*------------------------*
      | b) for use on T3E only |
      *------------------------*/

#else   /* if T3D is defined */

     sprintf(filenameout,"%s_%d",filenameout,_MPP_MY_PE);  
     sprintf(name,"%s_%d",argv[1],_MPP_MY_PE); /* 'eingabe' */
     if ((eingabe=fopen(name,"r"))==NULL) 
       {            
         fprintf(stderr,"couldn't find %s \n",name);
         exit(1);
       }
     
     /* read parameters from 'eingabe' file */

     fscanf(eingabe,"Steps=%d \n",&mcsmax);  
     fscanf(eingabe,"out=%d \n",&outstep);
     fscanf(eingabe,"seed=%d \n",&seed);
     fscanf(eingabe,"datei=%s \n",dname);     /* 'conf_out' */
     fscanf(eingabe,"reptation=%d\n",&nreptation);
     fscanf(eingabe,"local=%d\n",&local);
     fscanf(eingabe,"cbgc=%d\n",&ncbgc);
     fclose(eingabe);

     sprintf(dname,"%s_%d",dname,_MPP_MY_PE);

     if(seed<0){
              time_t seconds = time(NULL);
              seed = (int) seconds;
     } 
     /* get seed for rnd number generator r250 from cpu-time */
     seed += 1711*_MPP_MY_PE;            
     /* assign different seed for each processor */ 
     r250_srandom(seed);

     sysin(dname);           
     /* read sim. specifications and monomer vectors */ 
     /* from file "conf_out" */

     init_weight();          
     /* initializes weights for reweighted simulation */
#endif
 
  
     /*-----------all parameters read !-----------------*/


     Tinv  = -1.0/T;  

     calcIntegral(&zbondA,&zbondB);  
     /* update rtable[] and itable[] which are used */ 
     /* in cbgc.c (choose_length)                   */

     populate();                     
     /* populate boxes with start configuration */
    
     /* initialize check for total energy */
     e_start= eenergie();
     cenergie=0;    


     for (mcs=0;mcs<mcsmax;) {    

       /*----------*
        | MC steps |
        *----------*/
     
       /*-----------------------------------------------*
        | grand-canonical configurational-bias MC moves | 
        *-----------------------------------------------*/     

       /* adjust CP for configurational-bias move */
       CPA += T*(log(LSX*LSY*LSZ)+(polA-1)*log(zbondA));   
       CPB += T*(log(LSX*LSY*LSZ/(SIGMAB*SIGMAB*SIGMAB))+\
                 (polB-1)*log(zbondB));  /* SIGMAA := 1 !!! */

       for(i=0;i<ncbgc;i++) {
         cbgc_acc0   += cbgc(0); 
         /* cbgc_acc1   += cbgc(1); */
         cbgc_count0 += 1;
         /* cbgc_count1 += 1; */
       }

       /* reajust CP */ 
       CPA -= T*(log(LSX*LSY*LSZ)+(polA-1)*log(zbondA));    
       CPB -= T*(log(LSX*LSY*LSZ/(SIGMAB*SIGMAB*SIGMAB))+\
                 (polB-1)*log(zbondB)); 


       /*----------------*  
        | local MC moves |
        *----------------*/

       for(i=0;i<local;i++) {
         lakzeptiert += lmcmove();
         lgesamt     += nrmonA+nrmonB;
       }


       /*--------------------*
        | reptation MC moves |
        *--------------------*/

       for(i=0;i<nreptation;i++) {
         rep_acc   += reptation();       
         
         if (polA==1)
           rep_count += nrchainsB;
         else if (polB==1)        
           rep_count += nrchainsA; 
         else 
           rep_count += nrchainsA+nrchainsB; 
       }
           
       mcs++;   /* MC step counter++ */
     
       
       /*------------end of MC circle---------------*/       


       /*--------------------------------*
        | analysis after 'outstep' steps |
        *--------------------------------*/      

       if (mcs%outstep==0) {
 
         /* calculate total energy in two ways:  */
         /* 1) eenergie in wach.c                */
         /* 2) sum of energy fcts in element.c   */ 
        
         energy1 = eenergie();
         analyse (mcs_start+mcs,1); 
         energy2 = data.energy_ljA + data.energy_ljB +
                   data.energy_ljAB + data.energy_feneA + 
                   data.energy_feneB;

         if ( (energy1 - energy2) > 1e-6) {
           fprintf(stderr,"error in energy calculation:"
                          "energy(element.c)=%g energy"
                          "(wach.c)=%g\n", energy1, energy2);} 

         /* check if acc. energy differences (cenergie) */
         /* = change in total energy (energie1-e_start) */

         e_end = cenergie-(energy1-e_start);
         e_end = sqrt(e_end*e_end);

         if (e_end > 1e-6) {
           fprintf(stderr,"accumulated energy difference=%g" 
                          "total change (e_start-e_end)=%g\n",\
                           cenergie,e_start-energy1);}
       } /* mcstep */  
     }  


     /*------------------------------------------*
      | write acceptance rate data to 'log' file |
      *------------------------------------------*/
      
    if ((log_file=fopen("log_file","a"))==NULL) 
       {            
         fprintf(stderr,"couldn't find %s \n",name);
         exit(1);
       }

     fprintf(log_file,"gccb: p:%g\tm:%g\t local: %g\trep: %g\n"
             ,cbgc_acc0/(double)cbgc_count0\
             ,cbgc_acc1/(double)cbgc_count1\
             ,lakzeptiert/(double)lgesamt\
             ,rep_acc/(double)rep_count);

     fclose(log_file);

     sysout(dname,mcs_start+mcs);  
     /* write sim. specifications and monomer vectors */
     /* for next run to file "conf_out"               */
   }
   

   free(posA);
   free(posB);
   free(poslatA);
   free(poslatB);
   free(cm0A);
   free(cm0B);
   free(element_liste);

   time(&now);
   PRINTFS("program ends at %s\n",ctime(&now));
   
   return (OK);
}






