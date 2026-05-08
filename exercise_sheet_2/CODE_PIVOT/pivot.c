/*********************************
 * pivot.c: pivot MC move        *
 *                               *
 * last modification: 08/11/2004 *
 *********************************/


#include "element.h"
#include "r250.h"
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>


#define NRND_MAX 6*(NPOLYMAXA+NPOLYMAXB)
/* determine chain, determine pivot center, determine angle (2), determine side to be rotated, MC step */

extern int       nrchainsA, polA;
extern int       nrchainsB, polB;
extern MYVEC     *posA, *poslatA;
extern MYVEC     *posB, *poslatB;
extern MYVEC     *pos_tempA, *pos_tempB; /* temporary arrays only used in pivot.c */
extern double    Tinv,J;
extern double    LSX,LSY,LSZ,LSSUBX,LSSUBY,LSSUBZ;
extern long double    cenergie;
extern BBOX      boxliste[BOX_LENGTH_3];
extern DATA      data;
extern double    time1,time2;

/*********************
 * pivot: calls piv  *
 *********************/

int pivot(void) {
int piv_ctr=0;

  if (polA !=1)
     piv_ctr += piv(0); 
  if (polB !=1)
     piv_ctr += piv(1);
  return (piv_ctr);      
  /* returns # of successful pivot attempts */
}


/****************************************
 * piv: tries nrchains times to rotate  *
 *      one arbitrarily chosen end of a *
 *      polymer around an arbitrarily   *
 *      chosen pivot center             *
 *                                      *
 * type: A(=0) or B(=1)                 *
 *                                      *
 * remark: no bondlengths are changed ! *
 ****************************************/

int piv(int type){

static double *rndtab = NULL;
static double *rnd_ptr;

int    rwert = 0;
int    kette, ident_kette, pivot_center, ident_pivot_center, start, stop, dx; 
int    an, ab , acc, box_new,box_old;
int    i_poly,i,  ident2;
MYVEC  trials[polA],dummy[polA],dummy1[polA];
MYVEC  *pos_temp;
double energy_old=0,energy_new=0,energy_old1=0,energy_new1=0,energy_old2=0,energy_new2=0;
double prob, length;

int    nrchains, pol, flag;
MYVEC  *pos, *poslat;
MYVEC  bond[2];
double phi,cphi,sphi,theta,ctheta,stheta,psi,cpsi,spsi;
double x_old,y_old,z_old;
double time0;
time_t now;

  /* time0=time(&now); */

  /*-----------------------------------------*
   | set variables according to monomer type |
   *-----------------------------------------*/

  if (type==0) {
    pol      = polA;
    pos      = posA;
    pos_temp = pos_tempA;
    poslat   = poslatA;
    nrchains = nrchainsA;

    if (nrchainsA > NPOLYMAXA) {
      fprintf(stderr,"ERROR in rep.c: nrchainsA>=NPOLYMAXA\n");
      return(0);
    }
  }

  else if (type==1) {
    pol      = polB; 
    pos      = posB;
    pos_temp = pos_tempB;
    poslat   = poslatB; 
    nrchains = nrchainsB;
  
    if (nrchainsB > NPOLYMAXB) {
      fprintf(stderr,"ERROR in rep.c: nrchainsB>=NPOLYMAXB\n");
      return(0);
    }
  }

  else {
    fprintf(stderr,"Incorrect 'type' declaration in piv"
                   " - EXITING program ...\n");  
    exit(1);
  }


  /*------------------------------------------*
   | initialize or update random number array |
   *------------------------------------------*/ 

  if (rndtab==NULL) {  /* first call -> initialization */

    rndtab = (double*) calloc(NRND_MAX,sizeof(double));   
    /* allocate memory for NRND_MAX #s */

    double_r250_vector(rndtab,NRND_MAX);                  
    /* assign random #s to array       */

    rnd_ptr = rndtab;                                     
    /* set rnd_ptr to first position   */
  }

  acc = (int)(rnd_ptr-rndtab);
  if ( (acc<0)||(acc>NRND_MAX) ) {
    fprintf(stderr,"ERROR: must generate too many RNDs" 
                   "in rep.c (%d)\n",acc);
  }

  else {   /* update random number array      */

    double_r250_vector(rndtab,acc);                       
    /* refresh rnd #s which have already been used */

    rnd_ptr = rndtab;                                     
    /* set rnd_ptr to first position */
  }

  /*-----------------------------------------*
   | try nrchains pivot moves -              |
   | each time a polymer is chosen at random |
   *-----------------------------------------*/

  for (i_poly=0;i_poly<nrchains;i_poly++) {  

    energy_new=0;
    energy_old=0;

    /*---------------------------------------------------------------------------------------*
     | step 1: choose chain, pivot center, side of polymer to be rotated and rotation angles |  
     *---------------------------------------------------------------------------------------*/

    /* energy_old2=eenergie(); */

    /*------------------*
     | 1.1 choose chain |
     *------------------*/

    kette = (int)((nrchains) * (*rnd_ptr++) );
            /* choose chain */

    if(kette==nrchains) kette-=1;            /* if rnd_number is exactly 1 -> possible segmentation fault ! */
    ident_kette = pol * kette;                                       
    /*-------------------------*
     | 1.2 choose pivot center |
     *-------------------------*/

    pivot_center = (int)((pol) * (*rnd_ptr++) );   
    if(pivot_center==pol) pivot_center-=1;   /* same as above */
    ident_pivot_center = ident_kette + pivot_center;

    /* printf("%d\n",ident_pivot_center); */ 
    
    /*--------------------------*
     | 1.3 choose rotation side |
     *--------------------------*/ 

    if ( (*rnd_ptr++) < 0.5){      
      start = ident_kette;               
      stop  = ident_pivot_center;
      flag=1;
      /* printf("1\n"); */
    } 
    else {                        
      start = ident_pivot_center + 1;
      stop  = ident_kette + pol;
      flag=2;
      /* printf("2\n"); */
    }

    /* printf("%d\t%d\n\n",start,stop); */
   
    /*--------------------------------------------------------------------------*
     | 1.4 choose Euler angles (Phi,Theta,Psi convention)                       |
     |     - rotate vector by Psi around the z-axis (0<Phi<=2Pi; matrix D)      |
     |     - rotate res. vector by Theta around x-axis (0<=Theta<=Pi; matrix C) |
     |     - rotate res. vector by Psi around z-axis (0<=Psi<=2Pi; matrix B)    |
     |     - resulting matrix: A=BCD                                            |
     *--------------------------------------------------------------------------*/

    phi    = 2*PI*(*rnd_ptr++);
    cphi   = cos(phi);
    sphi   = sin(phi);

    theta  = PI*(*rnd_ptr++);
    ctheta = cos(theta);
    stheta = sin(theta);

    psi    = 2*PI*(*rnd_ptr);
    cpsi   = cos(psi);
    spsi   = sin(psi);

   /*------------------------------------------------------------------------------------*
    | step 2: rotate monomers, check interactions and store positions in temporary array |
    |         early rejection criterium: overlaps are rejected immediately               |
    *------------------------------------------------------------------------------------*/

    for(i=start;i<stop;i++)
    {
        /*--------------------------------------------------------------------------------------*
         | 2.1 determine new position of monomer bead and copy it to temporary array (pos_temp) |
         *--------------------------------------------------------------------------------------*/

    	length = DISTANCE(pos+i,pos+ident_pivot_center); 
        (trials+i)->monomer_number = (pos+i)->monomer_number;

        /* "move" pivot center to origin  */
        x_old = (pos+i)->x - ((pos+ident_pivot_center)->x);
        y_old = (pos+i)->y - ((pos+ident_pivot_center)->y);
        z_old = (pos+i)->z - ((pos+ident_pivot_center)->z);

        /* rotate */
        /* compare f.i. with http://mathworld.wolfram.com/EulerAngles.html */
    	(trials+i)->x = (cpsi*cphi-ctheta*sphi*spsi)*x_old  + (cpsi*sphi+ctheta*cphi*spsi)*y_old  + (spsi*stheta)*z_old; 
    	(trials+i)->y = (-spsi*cphi-ctheta*sphi*cpsi)*x_old + (-spsi*sphi+ctheta*cphi*cpsi)*y_old + (cpsi*stheta)*z_old;
    	(trials+i)->z = (stheta*sphi)*x_old                 + (-stheta*cphi)*y_old                + (ctheta)*z_old;  

        /* "move" pivot center back */
        (trials+i)->x += ((pos+ident_pivot_center)->x);
        (trials+i)->y += ((pos+ident_pivot_center)->y);
        (trials+i)->z += ((pos+ident_pivot_center)->z);

        /*----------------------------------------------------------*
         | 2.2 calculate LJ interaction at old and new position,    |
         |     exclude interactions with rotated section of polymer |
         *----------------------------------------------------------*/

        if(flag==1){ /* compare with 1.3 */
    	energy_old  += elj_energie_pivot(i,type,pos+i,ident_kette,(ident_pivot_center)); 
    	energy_new  += elj_energie_pivot(i,type,trials+i,ident_kette,(ident_pivot_center)); 
   
        /* printf("%d\t%lf\t",i,elj_energie_pivot(i,type,pos+i,ident_kette,(ident_pivot_center)) );
        printf("%d\t%lf\t",i,elj_energie_pivot(i,type,trials+i,ident_kette,(ident_pivot_center)) );
        printf("%lf\n",energy_new-energy_old); */
        }

        if(flag==2){ /* compare with 1.3 */
        energy_old  += elj_energie_pivot(i,type,pos+i,(ident_pivot_center),(ident_kette + pol-1));
        energy_new  += elj_energie_pivot(i,type,trials+i,(ident_pivot_center),(ident_kette + pol-1));
        }

        /*-------------------------------------------------------------------------* 
         | 2.3 early rejection criterium : overlap -> reject whole move right away |
         *-------------------------------------------------------------------------*/

         if( (energy_new-energy_old) > INFENERGY) {
           return(rwert);
         }
    }

    /*------------------------------*
     | step 3: Metropolis criterium |
     *------------------------------*/
    
    prob = exp(Tinv*(energy_new-energy_old));   
    if ( (prob>1) || ((*rnd_ptr++) <prob)) {         

    /* printf("%.10lf\t",energy_new-energy_old);*/
    /* MC move accepted ! */
      rwert++;

      /*--------------------------------------------------------*
       | accepted: YES: Transfer positions from temporary array |
       *--------------------------------------------------------*/

      for(i=start;i<stop;i++)  {

      VCOPY(trials+i, pos+i);  

      (trials+i)->x = (trials+i)->x - LSX*((int)((trials+i)->x/LSX));
      (trials+i)->y = (trials+i)->y - LSY*((int)((trials+i)->y/LSY));
      (trials+i)->z = (trials+i)->z - LSZ*((int)((trials+i)->z/LSZ));
  
      if ( (trials+i)->x < 0 ) (trials+i)->x += LSX ;
      if ( (trials+i)->y < 0 ) (trials+i)->y += LSY ;
      if ( (trials+i)->z < 0 ) (trials+i)->z += LSZ ;  /* ? */
    
      box_old = VBOXNR(poslat+i); 
      VCOPY(trials+i, poslat+i);       /* update poslat array */
      box_new = VBOXNR(poslat+i);

      /* update double-linked subbox list if necessary */ 
      if (box_old!=box_new && type==0) {   
        del_element(box_old, i); 
        put_element(box_new, i);
      }
    
      else if (box_old!=box_new && type==1) {
        del_element(box_old, NMONOMAXA+i);  
        put_element(box_new, NMONOMAXA+i);
      } 
      }
 
      /* energy_new2=eenergie();
      printf("%.10lf\n",energy_new2-energy_old2); */

      cenergie += (energy_new-energy_old);   
      /* printf("!%lf\t%lf\n",cenergie,(energy_new - energy_old)); */
      /* update energy control variable */
    }
  }

  /* if(flag==1){
  time1+=time(&now)-time0;
  }
 
  if(flag==2){
  time2+=time(&now)-time0;
  } */


  return(rwert);   
  /* returns number of successful pivot moves */
}


/***************************************************
 * elj_energy_pivot: similar to elj_energie but    *
 *                  excludes interactions of       *
 *                  monomers from exclude_start to *
 *                  exclude_stop                   *
 *                                                 *
 * function is only used in piv() !                *
 *                                                 *
 * type           =  A (=0) or B (=1)              *
 * ort            =  monomer vector                *
 * exclude_start  =  starts list of monomers to be *
 *                   excluded in interactions      *
 * exclude_stop   =  end of list                   * 
 ***************************************************/

double elj_energie_pivot(int ident, int type, MYVEC *mon, int exclude_start, int exclude_stop)
{

/*---------------------------------------------------------*
 | The program calculates the LJ interaction of particle   |
 | "mon" with all particles whose distance to "mon" <      |
 | interaction range. They are all located in a 3*3*3      |
 | box with the central subbox containing "mon" (size of   |
 | a subbox = interaction range).                          |
 |                                                         |
 | Monomers from pos+exclude_start to pos+exclude_stop are |
 | excluded.                                               |
 *---------------------------------------------------------*/
 
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

      /* printf("!%d\n",akt->monomer); */
      type2 = akt->type;
      /* determine interact. monomer type */
 
      if ((akt->monomer != ident) || (type2 != type)){
      /* exclude interaction of monomer with itself */
 
      if ( ((akt->monomer) < exclude_start) || ((akt->monomer) > exclude_stop) || (type2!=type) ) {
      /* exclude interactions with monomers from excluded list */
      /* only difference to elj_energie !                      */

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
        }
	/* printf("%d\t%d\n",ident,akt->monomer);*/
      } /* not from excluded list */ 
      } /* akt o.k. */
    }/* j */
 
#ifdef HARDWALLS
    }/* pbox2 o.k. */
#endif
 
  } /* i */

  return (4.0*ljenergie);
}
