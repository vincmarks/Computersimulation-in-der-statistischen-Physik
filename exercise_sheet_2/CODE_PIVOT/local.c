/*************************************
 * local.c: local MC move            *
 *          (tries to move monomers) *
 *                                   *
 * last modification: 23/05/2003     *
 *************************************/


#include "element.h"
#include "r250.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NRND_MAX     5*(NMONOMAXA+NMONOMAXB)   
/* (choose monomer, 3 coordinates and metropolis)
   times (total number of monomers)               */

extern int       nrmonA, polA;
extern int       nrmonB, polB;
extern MYVEC     *posA, *poslatA;
extern MYVEC     *posB, *poslatB;
extern double    Tinv,T;
extern double    LSX,LSY,LSZ,LSSUBX,LSSUBY,LSSUBZ;
extern long double    cenergie;

#ifdef STRIPE
extern double    SEwall;
#endif


/*******************************************
 * lmtrial: creates a trial vector for     *
 *          local MC moves                 *
 *                                         *
 * alt = monomer vector to be substituted  *
 * new = new trial vector coordinates      *
 *                                         *
 * rnd_ptr_adr = adress of random ptr      *
 * ident       = monomer position of 'alt' *
 *               in pos or poslat array    *
 *             = A (=0) or B (=1)          *
 *******************************************/

void lmtrial(MYVEC *new, MYVEC *alt, double ** rnd_ptr_adr)
{
double *rnd_ptr; 
double rad;
double ctheta,stheta,phi;

  rnd_ptr = *rnd_ptr_adr; 
  /* get pointer to rnd number table */

  rad     = MAX_STEP*(*rnd_ptr++); 
  /* determine distance from 'alt'   */

  /* determine arbitrary position on */  
  /* sphere around 'alt'             */

  ctheta = 2*(*rnd_ptr++) - 1;                 
  stheta = 2*(*rnd_ptr++) - 1;

  phi    = ctheta*ctheta + stheta*stheta;
  while (phi >= 1 || phi == 0) {
    ctheta = 2*double_r250() - 1;
    stheta = 2*double_r250() - 1;
    phi    = ctheta*ctheta + stheta*stheta;
  }

  new->z = alt->z + rad*(1.0 - 2*phi);
  phi    = 2*sqrt(1-phi);
  rad   *= phi;
  new->x = alt->x + rad*ctheta;
  new->y = alt->y + rad*stheta;

  *rnd_ptr_adr = rnd_ptr;                      
  /* set pointer back */
}


/*******************************************
 * lmcmove: executes nrmonA+nrmonB         *
 *          local moves, monomers are      *
 *          chosen at random, fct. returns *
 *          number of accepted moves       *
 *******************************************/

int lmcmove(void){

static double *rndtab = NULL;
static double *rnd_ptr;
int    ident,box_new,box_old,i_mono,i_poly,acc,type,pol;
MYVEC  *alt,*neu,*mon,*pos,*poslat,neud;
double dE;


  /*------------------------------------------*
   | initialize or update random number array |
   *------------------------------------------*/   

  if (rndtab==NULL) {                                        
  /* first call -> initialization    */

        rndtab = (double*) calloc(NRND_MAX,sizeof(double));  
        /* allocate memory for NRND_MAX #s */

        double_r250_vector(rndtab,NRND_MAX);                 
        /* generate all rndnumbers         */

        rnd_ptr = rndtab;  
        /* set rnd_ptr to first position   */
     }
 
  acc = (int)(rnd_ptr-rndtab);
  if ( (acc<0)||(acc>NRND_MAX) ) {
    fprintf(stderr,"ERROR: must generate too many RNDs" 
                   "in local.c (%d)\n",acc);
  }
 
  else {               /* update random number array      */   

    double_r250_vector(rndtab,acc); 
    /* refresh rnd #s which have already been used        */ 

    rnd_ptr = rndtab;  /* set rnd_ptr to first position   */
  }

  neu = &neud;

  
  /*----------------------------------*
   | loop over nrmonA+nrmonB monomers |
   | which are chosen at random       | 
   *----------------------------------*/

  for (i_mono=0,acc=0;i_mono<nrmonA+nrmonB;i_mono++) {  
     ident  = (int) ( ((double)(nrmonA+nrmonB))*(*rnd_ptr++) );
     /* choose monomer */

     if(ident==(nrmonA+nrmonB)) ident-=1;  /* if rnd_number is exactly 1 -> possible segmentation fault ! */

     /*-----------------------------------------*
      | set variables according to monomer type |
      *-----------------------------------------*/

     if (ident<nrmonA) {
       type   = 0;
       pol    = polA;
       pos    = posA;
       poslat = poslatA;
     }
     else {
       type   = 1;
       pol    = polB; 
       pos    = posB;
       poslat = poslatB;
       ident  = ident-nrmonA;
     }
     i_poly = ident/pol;  
     /* determine position of monomer in polymer */

     if ((i_poly >= NPOLYBRUSH) || (ident%pol > 0) || 
         (type==1) ) { /* not first monomer in a brush */


        
        /*----------------------------------------*
         | step 1: create trial vector 'neu' with |
         |         lmtrial-fct and calculate      |
         |         energy difference between old  |
         |         and new configuration          |
         *----------------------------------------*/       

        alt = pos + ident;

        lmtrial(neu,alt,&rnd_ptr);        
        /* generate trial vector 'neu' */

        dE  = fene(ident,type,neu)-fene(ident,type,alt);  
        dE += elj_energie(ident,type,neu)- 
              elj_energie(ident,type,alt);  

#ifdef HARDWALLS
        if (i_poly>=NPOLYBRUSH) { 
        /* monomer's not in a brush */

          dE += wenergy(neu)-wenergy(alt);
        }
        else {
          dE += wbenergy(neu)-wbenergy(alt);
        }
#endif


        /*-----------------------*
         | step 2: local MC step |
         *-----------------------*/   

        if ( (*rnd_ptr++) < exp(Tinv*dE)) {       


          /*-----------------------------*
           | accepted ? -> update arrays |
           *-----------------------------*/          

          box_old = VBOXNR(poslat+ident);
   
          alt->x = neu->x;            /* update pos-array */
          alt->y = neu->y;
          alt->z = neu->z;
    
          mon = poslat+ident;

          mon->x = alt->x-LSX*((int)(alt->x/LSX));
          mon->y = alt->y-LSY*((int)(alt->y/LSY));
          mon->z = alt->z-LSZ*((int)(alt->z/LSZ));

          /* update poslat-array */
          if ( mon->x < 0 ) mon->x += LSX ;        
          if ( mon->y < 0 ) mon->y += LSY ;
          if ( mon->z < 0 ) mon->z += LSZ ;

          box_new = VBOXNR(poslat+ident);    

          /* update double-linked subbox list if necessary */ 
          if (box_new !=box_old) { 
            if (type==0) {   
              del_element(box_old,ident);
              put_element(box_new,ident);
            }
            else {
              del_element(box_old,NMONOMAXA+ident);
              put_element(box_new,NMONOMAXA+ident);
            }
          }

          cenergie += dE; /* update energy control variable  */
          acc++;          /* number of accepted moves        */
        }
     } /* not first monomer in brush */
  }/* i_mono */

  return(acc);
}


/* undefine local definitions */
#undef NRND_MAX
