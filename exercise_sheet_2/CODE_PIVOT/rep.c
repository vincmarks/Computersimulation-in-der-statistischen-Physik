/***********************************************
 * rep.c: reptation MC move                    *
 *                                             *
 * last modification: 30/08/2000               *
 ***********************************************/


#include "element.h"
#include "r250.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NRND_MAX 5*(NPOLYMAXA+NPOLYMAXB)    
/* (choose chain and monomer which will be cut off,         */
/* 2* create new monomer, MC step) * total number of chains */ 

extern int       nrchainsA, polA;
extern int       nrchainsB, polB;
extern MYVEC     *posA, *poslatA;
extern MYVEC     *posB, *poslatB;
extern double    Tinv,J;
extern double    LSX,LSY,LSZ,LSSUBX,LSSUBY,LSSUBZ;
extern long double  cenergie;


/***************************************
 * reptation: calls rep                *
 ***************************************/

int reptation(void) {
int rep_ctr=0;

  if (polA !=1)
     rep_ctr += rep(0); 
  if (polB !=1)
     rep_ctr += rep(1);

  return (rep_ctr);      
  /* returns # of succesful reptation moves */
}


/****************************************
 * rep: tries nrchains times to cut off *
 *      one end of a polymer and to     *
 *      attach it to the opposite side  *
 *                                      *
 * type: A(=0) or B(=1)                 *
 *                                      *
 * remark: no bondlengths are changed ! *
 ****************************************/

int rep(int type){

static double *rndtab = NULL;
static double *rnd_ptr;

int    rwert = 0;
int    flag;
int    kette, ident, start, stop, dx; 
int    an, ab , acc, box_new,box_old;
int    i_poly,i,  ident2;
MYVEC  trials, *trial= &trials;
double energy_old,energy_new,energy_old1,energy_new1,ctheta,stheta,phi;
double prob,prob1, length;
double random_number;

int    nrchains, pol;
MYVEC  *pos, *poslat;

  /*-----------------------------------------*
   | set variables according to monomer type |
   *-----------------------------------------*/

  if (type==0) {
    pol      = polA;
    pos      = posA;
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
    poslat   = poslatB; 
    nrchains = nrchainsB;
  
    if (nrchainsB > NPOLYMAXB) {
      fprintf(stderr,"ERROR in rep.c: nrchainsB>=NPOLYMAXB\n");
      return(0);
    }
  }

  else {
    fprintf(stderr,"Incorrect 'type' declaration in rep"
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
   | try nrchains reptation moves -          |
   | each time a polymer is chosen at random |
   *-----------------------------------------*/

  for (i_poly=NPOLYBRUSH;i_poly<nrchains;i_poly++) {  

    /*------------------------------------------------------*
     | step 1: choose chain and reptation monomer at random |  
     *------------------------------------------------------*/

    kette = NPOLYBRUSH +
            (int)((nrchains-NPOLYBRUSH) * (*rnd_ptr++) );
            /* choose chain */

    if(kette==nrchains) kette-=1;            /* if rnd_number is exactly 1 -> possible segmentation fault ! */
    ident = pol * kette;                                       
  
    
    /*---------------------------------------*
     | decide whether first or last monomer  |
     | of the chosen polymer will be cut off |
     *---------------------------------------*/ 

    if ( (*rnd_ptr++) < 0.5){      
      an     = ident;               
      /* an = monomer on which a new monomer will be put  */

      ab     = ident + pol - 1;     
      /* ab = monomer which will be cut off - */
      /* here last monomer of chain           */

      ident2 = ab - 1;
      dx     = -1;            /* compare with step 3 */ 
    } 
    else {                        
      an     = ident + pol - 1;
      ab     = ident;
      /* first monomer of chain will be cut off */

      ident2 = ab + 1;
      dx     = 1;
    }
    start  = ab;
    stop   = an;

   
    /*--------------------------------------------------*
     | step 2: create a new monomer on opposite side of |
     |         chain- bondlength remains the same       |
     *--------------------------------------------------*/

    length = DISTANCE(pos+ab,pos+ident2);            
    /* bondlength between "ab" and its neighbor =  */
    /* bondlength between new monomer and "an"     */
    /* -> function doesn't change bondlengths !    */

    /* determine position of new monomer on sphere */
    /* around "an" with radius = old bondlength    */
    ctheta = 2*(*rnd_ptr++) - 1;                        

    stheta = 2*(*rnd_ptr++) - 1;                     

    phi    = ctheta*ctheta + stheta*stheta;        
    while (phi >= 1 || phi == 0) {
      ctheta = 2*double_r250() - 1;
      stheta = 2*double_r250() - 1;
      phi    = ctheta*ctheta + stheta*stheta;
    }
    trial->z = (pos+an)->z + length*(1.0 - 2*phi);
    phi      = 2*sqrt(1-phi)*length;
    trial->x = (pos+an)->x + phi*ctheta;
    trial->y = (pos+an)->y + phi*stheta;

    energy_old  = elj_energie(ab,type,pos+ab);  
    /* old energy at position ab without counting monomer ab */

    energy_new  = elj_energie(ab,type,trial);   
    /* new energy at position trial without counting ab */

    /*******************************************************************************/

    energy_old1 = eenergie();

    /*******************************************************************************/

#ifdef HARDWALLS
    energy_old += wenergy((pos+ab));
    energy_new += wenergy(trial);
#endif
  

    /*---------------------------*
     | step 3: reptation MC step |
     *---------------------------*/

    random_number = (*rnd_ptr++);
    
    prob = exp(Tinv*(energy_new-energy_old));    
    if ( (prob>1) || (random_number < prob) ) {         
    /* MC move accepted ! */
    flag=1;   
    /* printf("%.10lf\n",prob);*/ 

      rwert++;


      /*---------------------------------------*
       | move along chain and copy monomers to |
       | new positions in chain                |
       *---------------------------------------*/

      for(i=start;i!=stop;i+=dx){          
      /* copy old monomers to new positions */

        box_old = VBOXNR(poslat+i);       
        VCOPY(pos+i+dx, pos+i); 
        /* in pos array (cp pos+i+dx to pos+i)*/

        VCOPY(poslat+i+dx, poslat+i);      
        /* in poslat array                    */

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
    
      VCOPY(trial, pos+an);  
      /* copy newly created monomer to first or    */
      /* last position in chain - update pos array */

      trial->x = trial->x - LSX*((int)(trial->x/LSX));
      trial->y = trial->y - LSY*((int)(trial->y/LSY));
      trial->z = trial->z - LSZ*((int)(trial->z/LSZ));
  
      if ( trial->x < 0 ) trial->x += LSX ;
      if ( trial->y < 0 ) trial->y += LSY ;
      if ( trial->z < 0 ) trial->z += LSZ ; 
    
      box_old = VBOXNR(poslat+an); 
      VCOPY(trial, poslat+an);       /* update poslat array */
      box_new = VBOXNR(poslat+an);

      /* update double-linked subbox list if necessary */ 
      if (box_old!=box_new && type==0) {   
        del_element(box_old, an); 
        put_element(box_new, an);
      }
    
      else if (box_old!=box_new && type==1) {
        del_element(box_old, NMONOMAXA+i);  
        put_element(box_new, NMONOMAXA+i);
      }

      cenergie += (energy_new - energy_old);   
      /* update energy control variable */

      energy_new1=eenergie();
      prob1 = exp(Tinv*(energy_new1-energy_old1));
      /* printf("%.10lf\n",prob1);*/
      if ( (prob1>1) || (random_number < prob1) ) {
      /* MC move accepted ! */
      }
      //else{printf("1\n"); }

      //printf("%.10lf\t%.10lf\n",(energy_new-energy_old)-(energy_new1-energy_old1),prob-prob1); 
    }
  }


  return(rwert);   
  /* return number of succesful reptation moves */
}


