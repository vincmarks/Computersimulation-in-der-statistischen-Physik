/***************************************
 * element.c: box and energy functions *
 *                                     *
 * last modifications: 23/05/2003      *
 ***************************************/

#include"element.h"
#include<stdio.h>
#include<math.h>
#include <stdlib.h>

extern MYVEC     *posA, *poslatA;
extern int       nrmonA,polA;
extern MYVEC     *posB, *poslatB;
extern int       nrmonB,polB;

extern double    LSX,LSY,LSZ,LSSUBX,LSSUBY,LSSUBZ;
extern double    J,T;
extern BBOX      boxliste[NRSUB3];
extern ELEMENT   *element_liste;

static int box_periodic(int ix, int iy, int iz); 
/* help-function for populate() */


/***************************************
 * put_element: put monomer in box by  *
 *              setting pointers in    *
 *              elementliste           *
 *              -> double linked list  *
 *                 for each subbox     *                     
 *                                     *
 * ibox       = box number             *
 * element    = position of monomer in *
 *              elementliste           *       
 ***************************************/

void put_element(int ibox, int element) {
  ELEMENT *before, *particle;
  BBOX *pbox;

  particle = element_liste + element;
  pbox     = boxliste + ibox;

  /* monomer is not first particle in box */
  if (pbox->population != 0){        
    before=(pbox->first)->before;
    before->after    = particle;
    (pbox->first)->before=particle;
    particle->before = before;
    particle->after  = pbox->first;
    pbox->population++;
   }
  /* monomer is first particle in box */
  else {                            
    pbox->population++;
    pbox->first      = particle;
    particle->before = particle;
    particle->after  = particle;  
   }

}
    
 
/***************************************
 * del_element: delete element from    *
 *              box by setting ptrs in *
 *              elementliste           *
 *                                     *
 * ibox       = box number             *
 * element    = position of monomer    *
 *              in elementliste        *
 ***************************************/

void del_element(int ibox, int element) {
  ELEMENT *before,*after,*particle;
  BBOX *pbox;

  particle = element_liste + element;
  pbox     = boxliste + ibox;

  before = particle->before;
  after  = particle->after;

  before->after = after;
  after->before = before;
  pbox->first   = before;
  pbox->population--;

}
  

/*****************************************
 * box_periodic: determines box number   *
 *               for a set of normalized,*
 *               non-periodic, integer   *
 *               box coordinates         *  
 *                                       *
 * ix, iy, iz :  box coordinates         *
 *               (0,1, .. (NRSUBX,Y,Z)-1)* 
 *                                       *
 * remark:       this function is only   *
 *               called in populate()    *
 *****************************************/

static int box_periodic(int ix, int iy, int iz){
  
  /* non-periodic -> periodic coordinates */
  if (ix<0) ix += NRSUBX;          
  if (iy<0) iy += NRSUBY;
  if (iz<0) iz += NRSUBZ;
  
  if (ix>=NRSUBX) ix -= NRSUBX;
  if (iy>=NRSUBY) iy -= NRSUBY;

#ifdef HARDWALLS                   
/* HARDWALLS   rturn -1  invalid number */
  if (iz>=NRSUBZ) return(-1);
#else    
  if (iz>=NRSUBZ) iz -= NRSUBZ;
  return (BOXNR(ix,iy,iz));    
#endif
     
}


/****************************************
 * populate: populates boxes with       *
 *           start configuration        *
 *                                      *
 * remark:   called at the beginning of *
 *           each run                   *
 ****************************************/

void populate(void){
  int     i,ix,iy,iz,ibox;
  int     nx,ny,nz,n;
  MYVEC   *mon;
  BBOX    *pbox;
  ELEMENT *act;

  if ((LSSUBX<RANGE)||(LSSUBY<RANGE)||(LSSUBZ<RANGE)) {
     fprintf(stderr,"ERROR: length of subbox %g %g "
                    "%gsmaller than interaction range %g\n",\
                    LSSUBX,LSSUBY,LSSUBZ,RANGE);
     exit(ERROR);
  }


  /*--------------------------------------------*
   | step 1: set up neighbor list in pbox-array |
   *--------------------------------------------*/

  /* loop over all boxes (determines "specified box") */  
  for(ix=0;ix<NRSUBX;ix++)            
    for(iy=0;iy<NRSUBY;iy++)
      for(iz=0;iz<NRSUBZ;iz++){

        ibox = BOXNR(ix,iy,iz);
        pbox = boxliste+ibox;
        pbox->population=0;
        
        n=0;

        for(nx= -1;nx<=1;nx++)         
        /* inner boxes: 3*3*3 around "specified box" */  
          for(ny= -1;ny<=1;ny++)       
          /* neighbor-array position 0-26 */                
            for (nz= -1;nz<=1;nz++) {
              (pbox->neighbours)[n]=\
              box_periodic(ix+nx,iy+ny,iz+nz);
              n++;
            }

        for(nx=-2;nx<=2;nx++)         
        /* outer boxes: 5*5*5 around "specified box" */
          for(ny=-2;ny<=2;ny++)       
          /* neighbor array position 27-125 */
            for (nz= -2;nz<=2;nz++) {
              if ( (abs(nx)==2)||(abs(ny)==2)||(abs(nz)==2) ) {
                (pbox->neighbours)[n]=\
                box_periodic(ix+nx,iy+ny,iz+nz);
                n++;
               }
            }
      }
            

  /*-------------------------------------------------------*
   | step 2a: - set monomer# and monomer type              |
   |            in elementliste-array for type A particles |
   |          - put monomers of type A in boxes by calling |
   |            put_element()                              |
   *-------------------------------------------------------*/
  
  element_liste  =\
  (ELEMENT *)calloc(NMONOMAXA+NMONOMAXB,sizeof(ELEMENT));

  act = element_liste;
  mon = poslatA;

  for(i=0;i<nrmonA;i++,mon++,act++){
    act->monomer = i;  /* set monomer#          */
    act->type    = 0;  /* set monomer type (A!) */

    /* determine box number */   
    ibox = FBOXNR(mon->x,mon->y,mon->z); 
 
    if ((ibox<0) || (ibox>=NRSUB3)){
      fprintf(stderr,"ERROR: ibox out of range (%f %f %f)\n", 
                     mon->x,mon->y,mon->z);
      exit(ERROR);
    }
    put_element(ibox,i);  /* put element in box */
  }


  /*-------------------------------------------------------*
   | step 2b: - set monomer# and monomer type              |
   |            in elementliste-array for type B particles |
   |          - put monomers of type B in boxes by calling |
   |            put_element()                              |
   *-------------------------------------------------------*/

  act = element_liste+NMONOMAXA; 
  /* pointer to first typeB particle */

  mon = poslatB;

  for(i=0;i<nrmonB;i++,mon++,act++){
    act->monomer = i;  /* set monomer#              */
   /* monomer# gives pos. in pos & poslat arrays    */
   /* (two particles (A,B) may have the same number */

    act->type    = 1;      /* set monomer type (B!) */

     /* determine box number */
    ibox = FBOXNR(mon->x,mon->y,mon->z); 
 
    if ((ibox<0) || (ibox>=NRSUB3)){
      fprintf(stderr,"ERROR: ibox out of range (%f %f %f)\n",
              mon->x,mon->y,mon->z);
      exit(ERROR);
    }
    put_element(ibox,i+NMONOMAXA);  /* put element in box */
  }


  /*----------------------------------------------------*
   | step 3: verify that all particles are put in boxes | 
   *----------------------------------------------------*/

  for (ibox=0,i=0;ibox<NRSUB3;ibox++) 
      i += (boxliste+ibox)->population;
 
  if ( i != nrmonA+nrmonB) {
    fprintf(stderr,"ERROR: number of particles in box %d"
                   "and total number %d differ\n"\
                   ,i,nrmonA,nrmonB);
    exit(ERROR);
  }

}


/****************************************
 * elj_energie: calculates LJ energy    *
 *              of a single monomer     *
 *              "box" method            *
 *                                      *
 * ident = position of monomer in       *
 *         pos and poslat arrays        *
 * type  = A (=0) or B (=1)             *      
 * mon   = monomer vector               *
 ****************************************/

/*---------------------------------------------------------*
 | The program calculates the LJ interaction of particle   |
 | "mon" with all particles whose distance to "mon" <      |
 | interaction range. They are all located in a 3*3*3      | 
 | box with the central subbox containing "mon" (size of   |
 | a subbox = interaction range).                          |
 *---------------------------------------------------------*/  

double elj_energie(int ident, int type, MYVEC *mon) {
  int     ibox=0,i=0,j=0,type2=0;
  BBOX    *pbox1=NULL,*pbox2=NULL;
  ELEMENT *akt=NULL;
  MYVEC   *mon1=NULL,*mon2=NULL,diff,vhilf;
  double  rad2=0,rad6inv=0,rad6=0,ljenergie=0;
  
  ljenergie=0.0;
  mon1 = &vhilf;

  mon1->x = mon->x;          /* -> periodic coordinates */
  mon1->y = mon->y;
  mon1->z = mon->z;

  mon1->x -= LSX*((int)(mon1->x/LSX));
  mon1->y -= LSY*((int)(mon1->y/LSY));
  mon1->z -= LSZ*((int)(mon1->z/LSZ));
  
  if (mon1->x<0.0) mon1->x+=LSX;
  if (mon1->y<0.0) mon1->y+=LSY;
  if (mon1->z<0.0) mon1->z+=LSZ;


  ibox  = FBOXNR(mon1->x,mon1->y,mon1->z);    
  /* determine box number of "mon" particle */

  pbox1 = boxliste+ibox;


  /*----------------------------------------------*
   | scan each of the 3*3*3 subboxes around "mon" |
   *----------------------------------------------*/
  
  for(i=0;i<27;i++) {                       

#ifdef HARDWALLS
    if ((pbox1->neighbours)[i]>=0) {
#endif
   
    pbox2 = (pbox1->neighbours)[i]+boxliste;  
    akt = pbox2->first;


    /*------------------------------------------------*
     | scan double-linked list which contains         |
     | possible interaction monomers in specified box |
     *------------------------------------------------*/
    
    for(j=0;j<pbox2->population;j++,akt=akt->after){    

      type2 = akt->type;                                
      /* determine interact. monomer type */

      if ((akt->monomer != ident) || (type2 != type)){  
      /* exclude interaction of monomer with itself */

        if (type2==0)
          mon2 = poslatA + akt->monomer;
        else
           mon2 = poslatB + akt->monomer;

        /* determine distance between "mon" */
        /* and interaction particle         */
        diff.x = mon2->x - mon1->x;  
        diff.y = mon2->y - mon1->y;                     
        diff.z = mon2->z - mon1->z;
       
        /* Minimum Image Convention */
        diff.x = diff.x - LSX*((int)(2.0*diff.x/LSX));  
        diff.y = diff.y - LSY*((int)(2.0*diff.y/LSY));

#ifndef HARDWALLS
        diff.z = diff.z - LSZ*((int)(2.0*diff.z/LSZ));  
#endif

        rad2 = SQUARE(diff.x)+SQUARE(diff.y)+SQUARE(diff.z);
       

        /*---------------------*
         | calculate LJ energy |
         *---------------------*/

        if ((type==0)&&(type2==0)){  /* AA interaction */
          /* check if distance < interaction range */
          if (rad2 < RANGE2AA) {    
            rad6    = CUBE(rad2); 
            rad6inv = 1.0/rad6;
            ljenergie += LENJAA(rad6inv);
          }
        }

        else if ((type==1)&&(type2==1)){  /* BB interaction */
          if (rad2 < RANGE2BB) {
            rad6    = CUBE(rad2); 
            rad6inv = 1.0/rad6;
            ljenergie += LENJBB(rad6inv);
          }
        }

        else {  /* AB interaction */
                /* exclude interaction between particles of */
                /* different type at the same position      */
          if ((rad2 < RANGE2AB)) {  
            rad6    = CUBE(rad2); 
            rad6inv = 1.0/rad6;
            ljenergie += LENJAB(rad6inv);

          }
        } /* in WW range */
      }/* akt o.k. */
    }/* j */    

#ifdef HARDWALLS
    }/* pbox2 o.k. */
#endif

  } /* i */

  return (4.0*ljenergie);
}


/****************************************
 * fene: calculates additional FENE-    *
 *       potential of an inserted       *
 *       monomer                        *
 *                                      *
 * ident = position of monomer in       *
 *         pos and poslat arrays        *
 * type  = A (=0) or B (=1)             *      
 * neu   = monomer vector               *
 ****************************************/

double fene(int ident, int type, MYVEC *neu)
{
int mnr;
double energie,rad2,pol;
MYVEC *mon,*pos;

  if (type==0) {
     pol = polA;
     pos = posA;
  }

  else{
     pol = polB;
     pos = posB;
  }

  if (pol==1) return(0.0); 
  /* works also for the monomer system */

  energie=0.0;

  /* calculate monomer position in polymer */
  mnr = ident - pol*((int)((double)(ident)/((double)(pol)))); 

  mon = pos + ident;


  /*-----------------------*
   | calculate FENE energy |
   *-----------------------*/

  if (mnr == 0){               /* first monomer in chain */
     rad2 = DISTANCE2(neu,mon+1); 
     if (rad2*RNULL2INV<1.0) energie = BOND(rad2);           
     else     {              energie = T*INFENERGY; }
  } 

  else if (mnr == (pol-1)) {   /* last monomer in chain */
     rad2 = DISTANCE2(neu,(mon-1));
     if (rad2*RNULL2INV<1.0) energie = BOND(rad2);
     else     {              energie = T*INFENERGY; }
  }

  else {                       /* monomer in middle of chain */
     rad2 = DISTANCE2(neu,(mon-1));
     if (rad2*RNULL2INV<1.0) energie = BOND(rad2);
     else     {              energie = T*INFENERGY; }

     rad2 = DISTANCE2(neu,(mon+1));
     if (rad2*RNULL2INV<1.0) energie += BOND(rad2);
     else     {              energie = T*INFENERGY; }
  }
 
  return(energie);

}


/***************************************
 * eenergie: calculates total energy   *
 *           using functions defined in*
 *           element.c                 *
 *                                     *
 * remark:   also checks if monomers   *
 *           are in correct box        *     
 ***************************************/

double eenergie(void) {
  int i,j,ix,iy,iz,number=0;
  double evalue;
  MYVEC *mon;
  BBOX *pbox;
  ELEMENT *philf;


  /*--------------------------------------*
   | check if monomers are in correct box |
   *--------------------------------------*/

  for(ix=0;ix<NRSUBX;ix++)                   
    for(iy=0;iy<NRSUBY;iy++)
      for(iz=0;iz<NRSUBZ;iz++){
        pbox    = boxliste+BOXNR(ix,iy,iz);
        number += pbox->population;
        philf   = pbox->first;
        for(j=0;j<pbox->population;j++,philf=philf->after) { 
          if (philf->type==0)  mon = poslatA+philf->monomer;
          else                 mon = poslatB+philf->monomer;
          if ((mon->x < (((double)(ix))*LSSUBX))        ||
              (mon->y < (((double)(iy))*LSSUBY))        ||
              (mon->z < (((double)(iz))*LSSUBZ))        ||
              (mon->x > (((double)(ix))*LSSUBX+LSSUBX)) ||
              (mon->y > (((double)(iy))*LSSUBY+LSSUBY)) ||
              (mon->z > (((double)(iz))*LSSUBZ+LSSUBZ))   )
          {fprintf(stderr,"ERROR: monomer number %d : "
                         "%d %d %d\n",philf->monomer,ix,iy,iz);
           VPRINT(mon);}
         }
      }


  /*--------------------------------------*
   | check if number of monomers in box = |
   | total number of monomers             |  
   *--------------------------------------*/

  if (number != nrmonA+nrmonB) 
     fprintf(stderr,"ERROR: found %d monomer in boxes, "
                    "should be %d\n",number,nrmonA+nrmonB);

  
  /*------------------------*
   | calculate total energy |
   *------------------------*/

  for(i=0,evalue=0;i<nrmonA;i++) 
#ifdef HARDWALLS
    if (i<NPOLYBRUSH*polA)
      evalue += elj_energie(i,0,posA+i) +\
                fene(i,0,posA+i) + 2*wbenergy((posA+i) );
    else
      evalue += elj_energie(i,0,posA+i) + fene(i,0,posA+i) +
                2*wenergy( (posA+i) );
#else
      evalue += elj_energie(i,0,posA+i) + fene(i,0,posA+i);
#endif

  for(i=0;i<nrmonB;i++) 
#ifdef HARDWALLS
    evalue += elj_energie(i,1,posB+i) + fene(i,1,posB+i) +
              2*wenergy( (posB+i) );
#else
    evalue += elj_energie(i,1,posB+i) + fene(i,1,posB+i); 
#endif

  return (evalue/2); /* each interaction is counted twice ! */
}
