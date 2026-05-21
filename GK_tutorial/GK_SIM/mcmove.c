/***********************************************
 * mcmove.c: defines help-functions for cbgc.c *
 *                                             *
 * last modification : 23/05/2003              *
 ***********************************************/


#include "element.h"
#include "r250.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* substitution macros for choose_length() */

#define XSTART       (0.85)  /* minimum bondlength */
#define XEND         (1.15)  /* maximum bondlength */
#define MAX_INTEGRAL 1000    /* # steps for integral calc. */

#define DX           ((XEND-XSTART)/MAX_INTEGRAL)  
/* step-size for integral calc. */

extern int       nrmonA;
extern int       nrmonB;
extern MYVEC     *poslatA;
extern MYVEC     *poslatB;
extern double    Tinv,J;
extern double    LSX,LSY,LSZ,LSSUBX,LSSUBY,LSSUBZ;
extern BBOX      boxliste[NRSUB3];

/* rtable[i] = (XSTART+0.5*DX) + DX*i */
static double rtableA[MAX_INTEGRAL+5], rtableB[MAX_INTEGRAL+5];

/* itable[] stores area below P(l) from l=0 to rtable[i] */ 
static double itableA[MAX_INTEGRAL+5], itableB[MAX_INTEGRAL+5];

static int       wwlisteA[MAX_ANZAHL];     
/* interaction list for type A monomers */
static int       wwlistA;

static int       wwlisteB[MAX_ANZAHL];     
/* interaction list for type B monomers */
static int       wwlistB;


/************************************************
 * bondenergy: calculates  LJ+FENE energy of    *
 *             two monomers (same type) at      *
 *             given distance                   *
 *                                              *
 * type        =  A (=0) or B (=1)              *
 * length      =  distance between monomers     *
 ************************************************/

double bondenergy(int type, double length) {
  double benergie = 0;
  double rad2,rad6,rad6inv;

  rad2 = SQUARE(length);

  /* calc. LJ-energy for type A monomers */ 
  if (type==0 && rad2 < RANGE2AA) {        
    rad6      = CUBE(rad2);
    rad6inv   = 1.0/rad6;
    benergie += 4.0*LENJAA(rad6inv);
  }

  /* calc. LJ-energy for type B monomers */
  else if (type==1 && rad2 < RANGE2BB) {   
    rad6      = CUBE(rad2);
    rad6inv   = 1.0/rad6;
    benergie += 4.0*LENJBB(rad6inv);
  }

  else if (type!=0 && type!=1) {
  fprintf(stderr,"Wrong 'type' declaration in bondenergy() "
                 "... EXITING program\n");
  exit(1);
  }  

  benergie += BOND(rad2);  /* calc. FENE-energy */

  return benergie;

}


/***********************************************
 * calcIntegral: calculates area below         *
 *               probability(length) curve,    *
 *               stores data in arrays         *
 *               rtable[i] (=center of ith x-  *
 *               interval) and itable[i] (area *
 *               from 0 to interval i)         * 
 ***********************************************/
 
double calcIntegral(double *flaecheA, double *flaecheB) {
  int i;
  double f0,f1,x;
  double trapez,flaeche;
 
 
  /*--------------------------------*
   | update rtableA[] and itableA[] |
   *--------------------------------*/

  x       = XSTART;
  flaeche = 0.0;
  f0      = x*x*exp(Tinv*bondenergy(0,x));

  for (i=0;i<MAX_INTEGRAL;i++) {

    x += DX;    /* increment length by DX    */
    f1 = x*x*exp(Tinv*bondenergy(0,x)); 
    /* is proportional to P(l) ! */

    trapez    = (0.5*DX)*(f0+f1); 
    /* remark: P(l) is prop. to l^2 * Boltzmann-factor */
    
    flaeche  += trapez;
    rtableA[i] = x + (0.5*DX); /* save length */

    itableA[i] = flaeche;      
    /* save area up to this length */

    f0 = f1;
  }

  for(i=0;i<MAX_INTEGRAL;i++)  /* normalize area to one */
    itableA[i] /= flaeche;
 
  *flaecheA = (flaeche/(10*DX));        
  /* division by DX is not necessary ! - shifts CP */

  /*--------------------------------*
   | update rtableB[] and itableB[] |
   *--------------------------------*/

  x       = XSTART;
  flaeche = 0.0;
  f0      = x*x*exp(Tinv*bondenergy(1,x));

  for (i=0;i<MAX_INTEGRAL;i++) {

    x += DX;  /* increment length by DX    */
    f1 = x*x*exp(Tinv*bondenergy(1,x)); 
    /* is proportional to P(l) ! */

    trapez    = (0.5*DX)*(f0+f1);       
    /* remark: P(l) is prop. to l^2 * Boltzmann-factor */
    
    flaeche  += trapez;
    rtableB[i] = x + (0.5*DX); /* save length */

    itableB[i] = flaeche;     
    /* save area up to this length */

    f0 = f1;
  }

  for(i=0;i<MAX_INTEGRAL;i++)  /* normalize area to one */
    itableB[i] /= flaeche;  
 
  *flaecheB = (flaeche/(10*DX));

}


/***********************************************
 * choose_length: chooses bondlength such that *
 *                length is P(l) distributed   *
 *                                             *
 * type           = A (=0) or B (=1)           *
 * rnd            = random number (0..1)       *        
 ***********************************************/

/*---------------------------------------------------------*
 | The normalized area from l=l1 to l=l2 below the         |
 | probability distribution P(length) is equal to the      |
 | probability of choosing l in interval [l1,l2] when l is |
 | P(l)- distributed.                                      |
 *---------------------------------------------------------*/

double choose_length(int type, double rnd) {
  int i;

#ifdef FIXLENGTH
  return(FIXLENGTH);
#else

  if (type==0) {   /* choose bondlength for type A polymer */
    i=0;

    /* -> an interval is chosen according */
    /* to its area below P(l)             */

    while (rnd>itableA[i]) {    
      i++;                        
    }                                   
      
    return rtableA[i];  /* return chosen bondlength */
  }


  else if (type==1) { 
  /* choose bondlength for type B polymer */

    i=0;
    while (rnd>itableB[i]) {   
      i++;                     
    }                                   
      
    return rtableB[i];          
  }
                     
  else {
  fprintf(stderr,"Wrong 'type' declaration in "
                 "choose_bondlength() ... EXITING program\n");
  exit(1);
  }  
#endif

}


/************************************************
 * make_list: generates interaction list by     *
 *            saving the interaction monomers'  *
 *            type and position (in pos-arrays) *
 *            in wwliste-arrays                 *
 *                                              *
 * ident1   = position of monomer (in posA,B)   *
 *            to which wwliste refers           *
 * type     = A (=0) or B (=1)                  *
 * ort      = vector of monomer to which        *
 *            wwliste refers                    * 
 * radius   = interaction radius                *
 ************************************************/

/**********************************************************
 * The function exists in two versions. The first version *
 * is used for systems with a simulation box size of      *
 * 5*5*5 or greater ( -> LARGE is defined). The second    *
 * version is intended for smaller systems ( -> LARGE is  *
 * not defined). It creates an interaction list for all   *
 * W_ALL trial vectors in cbgc.c . These vectors are      *
 * located on a sphere around monomer "ort" and have an   *
 * interaction range determined by RANGE. To consider all *
 * vectors on the sphere the total interaction radius     *
 * adds up to (radius of the sphere + RANGE) <= 2*size of *
 * a subbox. Therefore all possible interaction partners  *
 * are located in a 5*5*5 box around the box which        *
 * contains "ort".                                        *    
 *                                                        *
 * 1. The first version loops over all 5*5*5 boxes around *
 *    the central box which contains "ort". It is only    *
 *    possible for simulation boxes >= 5*5*5. For large   *
 *    systems this method will speed up the program       *
 *    considerably.                                       *
 * 2. The second version puts all monomers (!) in the     *
 *    wwliste array. This is the only correct method for  *
 *    system sizes < 5*5*5.                               *
 **********************************************************/ 


/*---------------*
 | first version |
 *---------------*/
 
# ifdef LARGE
void make_list(int ident1, int type, MYVEC *ort, double radius)
{
  int i,j,i_box1,i_box2;
  MYVEC ortp,dist,*mon,*poslat;
  double rad2;
  double wwrad2; 
  ELEMENT *akt;
  BBOX *p_box1,*p_box2;
  
  wwlistA = 0;
  wwlistB = 0;
  wwrad2  = SQUARE(radius);

  if (radius>1.732*RANGE) 
    fprintf(stderr,"ERROR %g %g\n", radius,RANGE); 

  /*--------------------------------*
   | calculate periodic coordinates |
   *--------------------------------*/

  /* "ort" may be located outside the sim. box */
  /* -> put monomer back in box !              */

  ortp.x = ort->x - LSX*((int)(ort->x/LSX));  
  ortp.y = ort->y - LSY*((int)(ort->y/LSY));   
  ortp.z = ort->z - LSZ*((int)(ort->z/LSZ));
  
  if ( ortp.x < 0 ) ortp.x += LSX ;
  if ( ortp.y < 0 ) ortp.y += LSY ;
  if ( ortp.z < 0 ) ortp.z += LSZ ; 


  /*---------------------------*
   | generate interaction list |          
   *---------------------------*/

  i_box1 = FBOXNR(ortp.x,ortp.y,ortp.z);      
  /* determine box number of "ort" */

  p_box1 = boxliste+i_box1;

  for (i=0;i<125;i++) {  
  /* loop over inner and outer neighbor-boxes            */ 
  /* and the box itself (= 5*5*5 = 125 boxes altogether) */    

    i_box2 = p_box1->neighbours[i];  /* determine box    */

#ifdef HARDWALLS
    if (i_box2>=0) {
#endif

    p_box2 = boxliste + i_box2;
    akt    = p_box2->first;

    for (j=0;j<p_box2->population;j++,akt=akt->after) {  
      /* loop over double-linked list which */ 
      /* contains all monomers in box       */

      if ((akt->monomer != ident1) || (type!=akt->type)) { 
      /* "ort" is not part of the list */

        if (akt->type==0) {
          poslat = poslatA;
        }
        else {
          poslat = poslatB;
        }

        mon = poslat + akt->monomer; 

        dist.x = ortp.x - mon->x;                    
        dist.y = ortp.y - mon->y;
        dist.z = ortp.z - mon->z;

         /* minimum image convention */     
        dist.x = dist.x - LSX*((int)(2.0*dist.x/LSX));   
        dist.y = dist.y - LSY*((int)(2.0*dist.y/LSY));

#ifndef HARDWALLS
        dist.z = dist.z - LSZ*((int)(2.0*dist.z/LSZ));
#endif

        rad2 = SQUARE(dist.x)+SQUARE(dist.y) + SQUARE(dist.z);

        if (rad2<wwrad2){      /* check if distance <        */
          if (akt->type==0) {  /* total interaction radius   */
            wwlisteA[wwlistA] = akt->monomer;  /* and type=A */
            wwlistA++;  /* save monomer-position in wwlisteA */
          }
          else {
            wwlisteB[wwlistB] = akt->monomer;
            wwlistB++;
          }
        }
      }/* monomer "ort" excluded */
    }/* double-linked list o.k. */

#ifdef HARDWALLS
  }/* i_box2 o.k. */
#endif

  }/* loop over all neighboring boxes o.k. */ 
}


/*----------------*
 | second version |
 *----------------*/

#else
void make_list(int ident1, int type, MYVEC *ort,double radius) 
{
  int i;
  MYVEC ortp,dist,*mon;
  double rad2,rad6,rad6inv;
  double wwrad2;
  
  wwlistA = 0;
  wwlistB = 0;
  wwrad2  = SQUARE(radius);


  /*--------------------------------*
   | calculate periodic coordinates |
   *--------------------------------*/

  ortp.x = ort->x - LSX*((int)(ort->x/LSX));
  ortp.y = ort->y - LSY*((int)(ort->y/LSY));
  ortp.z = ort->z - LSZ*((int)(ort->z/LSZ));
  
  if ( ortp.x < 0 ) ortp.x += LSX ;
  if ( ortp.y < 0 ) ortp.y += LSY ;
  if ( ortp.z < 0 ) ortp.z += LSZ ; 
  

  /*--------------------------------------------------------*
   | generate interaction list: type A interaction monomers | 
   *--------------------------------------------------------*/
  
  mon = poslatA;
  for(i=0;i<nrmonA;i++,mon++){
    if( (type==1) || (type==0 && i!=ident1) ) {
  
      dist.x = ortp.x - mon->x;   
      dist.y = ortp.y - mon->y;
      dist.z = ortp.z - mon->z;
      
      /* minimum image convention */
      dist.x = dist.x - LSX*((int)(2.0*dist.x/LSX));  
      dist.y = dist.y - LSY*((int)(2.0*dist.y/LSY));

#ifndef HARDWALLS
      dist.z = dist.z - LSZ*((int)(2.0*dist.z/LSZ));
#endif
      
      rad2 = SQUARE(dist.x) + SQUARE(dist.y) + SQUARE(dist.z);
      if (rad2<wwrad2){
        wwlisteA[wwlistA] = i;
        wwlistA++;
      }
    }
  }


  /*--------------------------------------------------------*
   | generate interaction list: type B interaction monomers | 
   *--------------------------------------------------------*/

  mon = poslatB;
  for(i=0;i<nrmonB;i++,mon++){
    if( (type==0) || (type==1 && i!=ident1) ) {

      dist.x = ortp.x - mon->x;
      dist.y = ortp.y - mon->y;
      dist.z = ortp.z - mon->z;
      
      /* minimum image convention */
      dist.x = dist.x - LSX*((int)(2.0*dist.x/LSX));  
      dist.y = dist.y - LSY*((int)(2.0*dist.y/LSY));

#ifndef HARDWALLS
      dist.z = dist.z - LSZ*((int)(2.0*dist.z/LSZ));
#endif
      
      rad2 = SQUARE(dist.x) + SQUARE(dist.y) + SQUARE(dist.z);
      if (rad2<wwrad2){
        wwlisteB[wwlistB] = i;
        wwlistB++;
      }
    }
  }
}
#endif


/************************************************
 * lj_energy: calculates LJ interaction of      *
 *            monomer "ort" with all monomers   *
 *            in interaction lists "wwlisteA"+B *
 *                                              *
 * type       =  A (=0) or B (=1)               *
 * ort        =  monomer vector                 *
 ************************************************/
 
double lj_energy(int type, MYVEC *ort) {
  int i;
  MYVEC ortp,dist,*mon;
  double evalue,rad2,rad6,rad6inv;

  evalue = 0.0; 

  /*-----------------------------------------*
   | calculate periodic coordinates of "ort" |  
   *-----------------------------------------*/
  
  ortp.x = ort->x - LSX*((int)(ort->x/LSX));  
  ortp.y = ort->y - LSY*((int)(ort->y/LSY));
  ortp.z = ort->z - LSZ*((int)(ort->z/LSZ));

  if ( ortp.x < 0 ) ortp.x += LSX ;
  if ( ortp.y < 0 ) ortp.y += LSY ;
  if ( ortp.z < 0 ) ortp.z += LSZ ; 


  /*----------------------------------------*
   | calculate interaction of monomer "ort" |
   | with type A monomers in wwlisteA       |
   *----------------------------------------*/

  for(i=0;i<wwlistA;i++){           
  /* loop over all monomers in wwlisteA */

    /* calculate distance between "ort"   */
    /* and monomer from wwlisteA          */
    mon = poslatA + wwlisteA[i];

    dist.x = ortp.x - mon->x;       
    dist.y = ortp.y - mon->y;      
    dist.z = ortp.z - mon->z;

    /* minimum image convention */
    dist.x = dist.x - LSX*((int)(2.0*dist.x/LSX));  
    dist.y = dist.y - LSY*((int)(2.0*dist.y/LSY));

#ifndef HARDWALLS
    dist.z = dist.z - LSZ*((int)(2.0*dist.z/LSZ));
#endif

    rad2 = SQUARE(dist.x) + SQUARE(dist.y) + SQUARE(dist.z);

    if (type==0) {  /* monomer "ort" is of type A */

      if (rad2<RANGE2AA){ 
      /* check if distance < interaction range */
      /* ok -> calculate AA interaction        */

        rad6    = CUBE(rad2); 
        rad6inv = 1.0/rad6;         
        evalue += LENJAA(rad6inv);
      }
    }

    else {  /* monomer "ort" is of type B            */

      if (rad2<RANGE2AB){           
      /* check if distance < interaction range */
      /* ok -> calculate AB interaction        */

        rad6    = CUBE(rad2);       
        rad6inv = 1.0/rad6;      
        evalue += LENJAB(rad6inv);
      }
    }
  }


  /*----------------------------------------*
   | calculate interaction of monomer "ort" |
   | with type B monomers in wwlisteB       |
   *----------------------------------------*/

  for(i=0;i<wwlistB;i++){
  /* loop over all monomers in wwlisteB */

    mon = poslatB + wwlisteB[i];

    dist.x = ortp.x - mon->x;      
    dist.y = ortp.y - mon->y;
    dist.z = ortp.z - mon->z;

    /* minimum image convention */
    dist.x = dist.x - LSX*((int)(2.0*dist.x/LSX));  
    dist.y = dist.y - LSY*((int)(2.0*dist.y/LSY));

#ifndef HARDWALLS
    dist.z = dist.z - LSZ*((int)(2.0*dist.z/LSZ));
#endif

    rad2 = SQUARE(dist.x) + SQUARE(dist.y) + SQUARE(dist.z);

    if (type==0) {                  /* AB interaction */
      if (rad2<RANGE2AB){
        rad6    = CUBE(rad2);
        rad6inv = 1.0/rad6;
        evalue += LENJAB(rad6inv);
      }
    }
    else {                          /* BB interaction */
      if (rad2<RANGE2BB){
        rad6    = CUBE(rad2);
        rad6inv = 1.0/rad6;
        evalue += LENJBB(rad6inv);
      }
    }
  }
  evalue *= 4.0;

  return(evalue);

}
