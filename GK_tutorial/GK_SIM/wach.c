/*********************************************************
 * wach.c: analysis fcts-update data-array before output *
 *         file "histo.dat" is written                   *
 *                                                       *
 * WARNING : wall routines have not been adjusted !      * 
 * last modification: 23/05/2003                         *  
 *********************************************************/


#include "element.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern int       nrchainsA, nrmonA, polA;
extern MYVEC     *posA, *poslatA, *cm0A;
extern int       nrchainsB, nrmonB, polB;
extern MYVEC     *posB, *poslatB, *cm0B;
extern double    w[NPOLYMAXA+1][NPOLYMAXB+1];

extern double    LSX,LSY,LSZ;
extern double    J,Ewall;
extern long      _MPP_MY_PE,_MPP_N_PES;
extern DATA      data;
extern double    Jwallwtab[JwallN_MAX];   
extern int       switch_iwall;



/******************************************
 * analyse: updates data array by calling *
 *          most analysis functions,      *
 *          updates "histo.dat" output    *
 *          file                          *
 ******************************************/

void analyse(int zeit, int lauf){      
/* zeit = number of MC steps */

  FILE   *datei;
  char   fname[80];
  double weightA, weightB;


  /*-------------------------*
   | call analysis functions |
   *-------------------------*/

  bondwinkel();        /* calc bondangles                   */
  bondlength();        /* bondlength, b**2, bmax            */
  calc_polymer();      /* end-to-end and radius of gyration */

  calc_nrchainsleft();     
  /* calc number of chains in the left half of the sim. box */

  calc_msd();          /* calc mean square displacement     */
  calc_energy();       /* calc energy                       */
 
#ifdef PRESSURE
  calc_p();            /* calc pressure via virial          */
#endif

#ifdef SWITCH
  weight=Jwallwtab[switch_iwall];
#else
  weightA=0;
  weightB=0;
#endif

# ifndef T3D
  sprintf(fname,"histo.dat");
# else
  sprintf(fname,"histo.dat_%d",_MPP_MY_PE);
# endif

  if ((datei=fopen(fname,"a"))==NULL) {
    fprintf(stderr,"ERROR: Couldn't open %s\n",fname);
    exit(ERROR); }  


  /*------------------------------*
   | write results to "histo.dat" |
   *------------------------------*/

  /* fprintf(datei,"%d %d %d %g %g %g %g %g\n",
          zeit,nrchainsA,nrchainsB,data.energy_ljA,
          data.energy_ljB,data.energy_ljAB,data.energy_feneA,
          data.energy_feneB); */

  fprintf(datei,"%d %d %d %g %g %g %g %g %g %g %g %g %g %g %g"
                "\n",zeit,nrchainsA,nrchainsB,data.energy_ljA,
                data.energy_ljB,data.energy_ljAB,
                data.energy_feneA,data.energy_feneB,data.px,
                data.py,data.pz,data.virialx,data.virialy,
                data.virialz,w[nrchainsA][nrchainsB]); 

  fclose(datei);

  /* # ifndef T3D
    sprintf(fname,"analyse"); 
  # else
    sprintf(fname,"analyse_%d",_MPP_MY_PE);
  # endif

  if ((datei=fopen(fname,"a"))==NULL) {
     fprintf(stderr,"ERROR: Couldn't open %s\n",fname);
     exit(ERROR); }

    fprintf(datei,"%d %d %d %g %g %g %g %g %g\n",
                zeit,nrchainsA,nrchainsB,data.bondangleA,
                data.bondlength1A,data.endA,
                data.gyrA,data.msdA,
                w[nrchainsA][nrchainsB]); 

  fclose(datei); */

}


/*****************************************
 * bondwinkel: calculates cosine of      *
 *             bond-angle and averages   *
 *             over all bonds            *
 *****************************************/

void bondwinkel() {
  int i,j;
  double awinkel;
  MYVEC *mon;
  MYVEC bond[2];

  
  /*-----------------------------------------------------*
   | calculate average cos(bondangle) for polymer type A |  
   *-----------------------------------------------------*/

  awinkel=0;

  for(i=0;i<nrchainsA;i++){
    
    mon=posA+i*polA;

    for(j=0;j<(polA-2);j++,mon++) {
      DIFF(mon+1,mon,bond);
      DIFF(mon+1,mon+2,bond+1);
      awinkel+= CWINKEL(bond,bond+1);  
      /* cos(bondangle) ! - see def. element.h */
    }
  }

  if (nrchainsA>0)
    awinkel /= ((double)(nrmonA-2*nrchainsA));
  else
    awinkel=0;
  
  /* update data array */
  data.bondangleA = awinkel;


  /*-----------------------------------------------------*
   | calculate average cos(bondangle) for polymer type B |  
   *-----------------------------------------------------*/

  awinkel=0;  

  for(i=0;i<nrchainsB;i++){
    
    mon=posB+i*polB;

    for(j=0;j<(polB-2);j++,mon++) {
      DIFF(mon+1,mon,bond);
      DIFF(mon+1,mon+2,bond+1);
      awinkel+= CWINKEL(bond,bond+1);
    }
  }

  if (nrchainsB>0)
    awinkel /= ((double)(nrmonB-2*nrchainsB));
  else
    awinkel=0;

  data.bondangleB = awinkel;

}


/*****************************************
 * bondlength: calculates av.bondlenght, *
 *             b**2, bmax                *
 *****************************************/

void bondlength() {    
  int i,j;
  double length,abslength,maxlength,distance;
  MYVEC *mon;


  /*----------------------------------------*
   | calculate average bondlength, b**2 and |
   | bmax for type A polymer                |
   *----------------------------------------*/

  length=0;    /* average bondlength */
  abslength=0; /* (average bondlength)^2 */ 
  maxlength=0; /* maximum bondlength */

  for(i=0;i<nrchainsA;i++){
    
    mon=posA+i*polA;

    for(j=0;j<(polA-1);j++,mon++) {

      distance   = DIST(mon,mon+1);
      length    += (distance);
      abslength += SQUARE(distance);
      if (ABS(distance)>maxlength) maxlength=ABS(distance);

      if (ABS(distance)>2.0) {                             
      /* check if chain is broken */                          

        fprintf(stderr,"WARNING: typeA: chain %d mon "
                       "%i\nlaenge=%f\n",i,j,ABS(distance));

        fprintf(stderr,"WARNING: monomer 1   %f %f %f\n",
                       mon->x,mon->y,mon->z);
        fprintf(stderr,"WARNING: monomer 2   %f %f %f\n",
                       (mon+1)->x,(mon+1)->y,(mon+1)->z);
      }
    }
  }
  
  if ((nrchainsA>0)&&(polA>1))
    data.bondlength1A = length/((double)(nrmonA-nrchainsA));
  else
    data.bondlength1A = 0;
    
  if ((nrchainsA>0)&&(polA>1))
    data.bondlength2A = abslength/((double)(nrmonA-nrchainsA));
  else
    data.bondlength2A = 0;

  data.bondlengthmaxA = maxlength;


  /*------------------------------------*
   | calculate average bondlength, b**2 |
   | and bmax for type B polymer        |
   *------------------------------------*/

  length=0; 
  abslength=0;  
  maxlength=0;

  for(i=0;i<nrchainsB;i++){
    
    mon=posB+i*polB;

    for(j=0;j<(polB-1);j++,mon++) {

      distance   = DIST(mon,mon+1);
      length    += (distance);
      abslength += SQUARE(distance);
      if (ABS(distance)>maxlength) maxlength=ABS(distance);
      if (ABS(distance)>2.0) {

        fprintf(stderr,"WARNING: typeB: chain %d mon "
                       "%i\nlaenge=%f\n",i,j,ABS(distance));

        fprintf(stderr,"WARNING: monomer 1   %f %f %f\n",
                       mon->x,mon->y,mon->z);

        fprintf(stderr,"WARNING: monomer 2   %f %f %f\n",
                       (mon+1)->x,(mon+1)->y,(mon+1)->z);
      }
    }
  }
   
  if ((nrchainsB>0)&&(polB>1))
    data.bondlength1B = length/((double)(nrmonB-nrchainsB));
  else
    data.bondlength1B = 0;
    
  if ((nrchainsB>0)&&(polB>1))
    data.bondlength2B = abslength/((double)(nrmonB-nrchainsB));
  else
    data.bondlength2B = 0;

  data.bondlengthmaxB = maxlength;

}


/*****************************************
 * calc_polymer: calculates end-to-end   *
 *               distance and radius of  *
 *               gyration                *
 *****************************************/
  
void calc_polymer(){

  double end,gyration;
  MYVEC *mon1,*mon2;
  int i,j,k;


  /*---------------------------------------*
   | calculate end-to-end distance and     |
   | radius of gyration for type A polymer |
   *---------------------------------------*/

  end      = 0;
  gyration = 0;

  for(i=0;i<nrchainsA;i++){

    mon1 = posA + i*polA;
    end += DISTANCE2(mon1,mon1+polA-1);

    for(j=0;j<polA;j++,mon1++) {
      mon2 = posA +i*polA;
      for(k=0;k<polA;k++,mon2++)
        gyration += DISTANCE2(mon1,mon2);
    }
  }
  if (nrchainsA>0) {
    end      /= nrchainsA;
    gyration /= (2.0*SQUARE(polA)*nrchainsA);
  }
  else {
    end      = 0;
    gyration = 0;
  }
 
  data.endA = end;
  data.gyrA = gyration;


  /*---------------------------------------*
   | calculate end-to-end distance and     | 
   | radius of gyration for type B polymer |
   *---------------------------------------*/

  end      = 0;
  gyration = 0;

  for(i=0;i<nrchainsB;i++){

    mon1 = posB + i*polB;
    end += DISTANCE2(mon1,mon1+polB-1);

    for(j=0;j<polB;j++,mon1++) {
      mon2 = posB +i*polB;
      for(k=0;k<polB;k++,mon2++)
        gyration += DISTANCE2(mon1,mon2);
    }
  }
  if (nrchainsB>0) {
    end      /= nrchainsB;
    gyration /= (2.0*SQUARE(polB)*nrchainsB);
  }
  else {
    end      = 0;
    gyration = 0;
  }
 
  data.endB = end;
  data.gyrB = gyration;

}


/*****************************************
 * calc_nrchainsleft: calculates #chains *
 *                    in "left" half of  *
 *                    simulation box     *
 *****************************************/

void calc_nrchainsleft() {
  int i_mono,i_poly,nrchainsleft;
  double z;
  MYVEC  *mon;


  /*---------------------------------------------------*
   | calculate #chains in left half for type A polymer |
   *---------------------------------------------------*/

  nrchainsleft=0; 
  for (i_poly=0;i_poly<nrchainsA;i_poly++) {
  
    mon = posA + i_poly*polA;
    z   = 0;
  
    for (i_mono=0;i_mono<polA;i_mono++,mon++) {
      z += mon->z;
    }

    z  /= polA;
    if (2.0*z < LSZ) nrchainsleft++;
  }
  data.nrchainsleftA = nrchainsleft;
 

  /*---------------------------------------------------*
   | calculate #chains in left half for type B polymer |
   *---------------------------------------------------*/

  nrchainsleft=0; 
  for (i_poly=0;i_poly<nrchainsB;i_poly++) {
  
    mon = posB + i_poly*polB;
    z   = 0;
  
    for (i_mono=0;i_mono<polB;i_mono++,mon++) {
      z += mon->z;
    }

    z  /= polB;
    if (2.0*z < LSZ) nrchainsleft++;
  }
  data.nrchainsleftB = nrchainsleft;

}


/******************************************
 * calc_cm: calculates center-of-mass of  *
 *          polymer at position "start"   *
 *          in pos-array                  *
 ******************************************/

void calc_cm(int pol, MYVEC *start, MYVEC *result){

  int j;
  MYVEC *mon;

  result->x = 0;
  result->y = 0;
  result->z = 0;

  mon = start;

  for(j=0;j<pol;j++,mon++) {
    result->x += mon->x;
    result->y += mon->y;
    result->z += mon->z;
  }

  result->x /= pol;
  result->y /= pol;
  result->z /= pol;

}


/******************************************
 * calc_msd: calculates mean-square-      *
 *           displacement                 *
 ******************************************/

void calc_msd() {

  int i;
  MYVEC cm, *mon, *vhilf;
  double rmove;


  /*-------------------------------------------------------*
   | calculate mean-square-displacement for type A polymer |
   *-------------------------------------------------------*/ 

  vhilf= cm0A;
  rmove=0;

  for(i=0;i<nrchainsA;i++,vhilf++) {

    mon = posA + i*polA;
    calc_cm(polA,mon,&cm);

    rmove += (SQUARE((vhilf->x - cm.x))+\
              SQUARE((vhilf->y - cm.y))+\
              SQUARE((vhilf->z - cm.z)));
  }
 
  if (nrchainsA>0)
    data.msdA = rmove/nrchainsA;
  else
    data.msdA = 0;


  /*-------------------------------------------------------*
   | calculate mean-square-displacement for type B polymer |
   *-------------------------------------------------------*/  

  vhilf= cm0B;
  rmove=0;

  for(i=0;i<nrchainsB;i++,vhilf++) {

    mon = posB + i*polB;
    calc_cm(polB,mon,&cm);
    
    rmove += (SQUARE((vhilf->x - cm.x))+\
              SQUARE((vhilf->y - cm.y))+\
              SQUARE((vhilf->z - cm.z)));
  }

  if (nrchainsB>0)
    data.msdB = rmove/nrchainsB;
  else
    data.msdB = 0;

}


/*****************************************
 * calc_energy: calculates total FENE,   *
 *              LJ and wall energy       *
 *****************************************/

void calc_energy(){

  int i,j;
  double poten,poly,zcount,wall,rad2,rad6,rad6inv;
  MYVEC *mon1,*mon2,diff;


  /*------------------------------------------------------*
   | step 1a: calculate FENE potential for type A polymer |
   *------------------------------------------------------*/
 
  poly = 0.0;  

  for(i=0;i<nrchainsA;i++){

    mon1 = posA + i*polA;
    mon2 = mon1 + 1;
    
    for(j=0;j<polA-1;j++,mon1++,mon2++){
      rad2 = DISTANCE2(mon1,mon2);
      poly += BOND(rad2);             /* FENE potential */
    }
  }
  
  data.energy_feneA = poly;
  

  /*------------------------------------------------------*
   | step 1b: calculate FENE potential for type B polymer |
   *------------------------------------------------------*/
 
  poly = 0.0;  

  for(i=0;i<nrchainsB;i++){

    mon1 = posB + i*polB;
    mon2 = mon1 + 1;
    
    for(j=0;j<polB-1;j++,mon1++,mon2++){
      rad2 = DISTANCE2(mon1,mon2);
      poly += BOND(rad2);          
    }
  }
  
  data.energy_feneB = poly;


  /*------------------------------------------------------*
   | step 2a: calculate LJ energy between type A polymers |
   *------------------------------------------------------*/

  poten  = 0.0;
  zcount = 0.0;
  mon1 = poslatA;

  for(i=0;i<nrmonA;i++,mon1++){
   
    mon2 = mon1+1;
  
    for(j=i+1;j<nrmonA;j++,mon2++){ 
    /* +2 because of bond energy */

      diff.x = mon1->x - mon2->x;
      diff.y = mon1->y - mon2->y;
      diff.z = mon1->z - mon2->z;

      /* Minimum Image Convention */
      
      diff.x = diff.x - LSX*((int)(2.0*diff.x/LSX));
      diff.y = diff.y - LSY*((int)(2.0*diff.y/LSY));

#ifndef HARDWALLS
      diff.z = diff.z - LSZ*((int)(2.0*diff.z/LSZ));  
#endif
      
      rad2 = SQUARE(diff.x)+SQUARE(diff.y)+SQUARE(diff.z);

      if (rad2 < RANGE2AA) {  
      /* check if distance between monomers < */
      /* interaction range                    */

        rad6    = CUBE(rad2);
        rad6inv = 1.0/rad6;
        poten  += LENJAA(rad6inv);
        zcount += 1-0.5*rad6;
      }
    }
  }
  poten *= 4.0;
  
  data.energy_ljA   = poten;
  data.energy_sqA   = zcount;


  /*------------------------------------------------------*
   | step 2b: calculate LJ energy between type B polymers |
   *------------------------------------------------------*/

  poten  = 0.0;
  zcount = 0.0;
  mon1 = poslatB;

  for(i=0;i<nrmonB;i++,mon1++){
   
    mon2 = mon1+1;
  
    for(j=i+1;j<nrmonB;j++,mon2++){ 

      diff.x = mon1->x - mon2->x;
      diff.y = mon1->y - mon2->y;
      diff.z = mon1->z - mon2->z;

      /* Minimum Image Convention */
      
      diff.x = diff.x - LSX*((int)(2.0*diff.x/LSX));
      diff.y = diff.y - LSY*((int)(2.0*diff.y/LSY));

#ifndef HARDWALLS
      diff.z = diff.z - LSZ*((int)(2.0*diff.z/LSZ));  
#endif
      
      rad2 = SQUARE(diff.x)+SQUARE(diff.y)+SQUARE(diff.z);

      if (rad2 < RANGE2BB) {
        rad6    = CUBE(rad2);
        rad6inv = 1.0/rad6;
        poten  += LENJBB(rad6inv);
        zcount += 1-0.5*rad6;
      }
    }
  }
  poten *= 4.0;
  
  data.energy_ljB   = poten;
  data.energy_sqB   = zcount;


  /*-----------------------------------------------------*
   | step 2c: calculate (type A - type B) LJ interaction |
   *-----------------------------------------------------*/

  poten  = 0.0;
  zcount = 0.0;
  mon1 = poslatA;

  for(i=0;i<nrmonA;i++,mon1++){

    mon2 = poslatB;   

    for(j=0;j<nrmonB;j++,mon2++){ 

      diff.x = mon1->x - mon2->x;
      diff.y = mon1->y - mon2->y;
      diff.z = mon1->z - mon2->z;

      /* Minimum Image Convention */
      
      diff.x = diff.x - LSX*((int)(2.0*diff.x/LSX));
      diff.y = diff.y - LSY*((int)(2.0*diff.y/LSY));

#ifndef HARDWALLS
      diff.z = diff.z - LSZ*((int)(2.0*diff.z/LSZ));  
#endif
      
      rad2 = SQUARE(diff.x)+SQUARE(diff.y)+SQUARE(diff.z);

      if (rad2 < RANGE2AB && rad2 != 0) {
        rad6    = CUBE(rad2);
        rad6inv = 1.0/rad6;
        poten  += LENJAB(rad6inv);
        zcount += 1-0.5*rad6;
      }
    }
  }
  poten *= 4.0;
  
  data.energy_ljAB   = poten;
  data.energy_sqAB   = zcount;


  /*------------------------------------*
   | step 3: calculate wall energy      | 
   |                                    |
   | WARNING: this routine has not been |
   |          adjusted to two species ! |
   *------------------------------------*/
  
  wall   = 0.0;

/*  mon1 = poslat;

    for (i=0;i<nrmon;i++,mon1++){
      if (i<NPOLYBRUSH*pol) wall += wbenergy(mon1);
      else                  wall += wenergy(mon1);
    }
*/

  data.energy_wall = wall;

}
