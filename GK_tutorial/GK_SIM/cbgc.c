/***********************************************************
 * cbgc.c: grand-canonical configurational-bias MC move    *
 *         (puts chains in box - deletes chains from box)  *
 *                                                         *
 * last modification: 23/05/2003                           * 
 ***********************************************************/


#include "element.h"
#include "r250.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NRND_MAX 10*((2*W_ALL+2)*pol+4)

extern int       nrchainsA, nrmonA, polA;
extern int       nrchainsB, nrmonB, polB;
extern MYVEC     *posA, *poslatA;
extern MYVEC     *posB, *poslatB;
extern double    Tinv,CPA,CPB;
extern double    LSX,LSY,LSZ,LSSUBX,LSSUBY,LSSUBZ;
extern double    cenergie;
extern double    w[NPOLYMAXA+1][NPOLYMAXB+1];
extern ELEMENT   *element_liste;


/**********************************************
 * cbgc: attemps to either insert (50% prob)  *
 *       or delete (50% prob) a chain from    *
 *       simulation box                       *
 *                                            *
 * type: A (=0) or B (=1)                     *            
 **********************************************/

double cbgc(int type) {
  int count=0; 
 
  if (r250_random()&1)    
    {
      if (nrchainsA< MAXA) 
         count=insert(type);           
         /* try to insert chain with 50% prob. */
    }

  else 
    {
      if (nrchainsA> MINA)
         count=delete(type);           
         /* try to delete chain with 50% prob. */
    }

  return(count);
}


/**********************************************
 * insert: attemps to insert a chain into     *
 *         simulation box                     *
 *                                            *
 * type:   A (=0) or B (=1)                   *
 **********************************************/

int insert(int type) {
  static double *rndtab = NULL;
  static double *rnd_ptr;
  static int    *box;
  static MYVEC  *new;
  double        ctheta,stheta,phi,length;

  double        energy_new,znew,x,y,z,prn;             
                /* znew = Rosenbluth-weight */

  double        energy[W_ALL],prob[W_ALL+1];           
                /* energy_new = energy of new chain */

  int           i_mono,i,acc,ident,i_chosen,i_break; 
  MYVEC         *mono,trial[W_ALL];

  int           nrchains,pol;
  MYVEC         *pos, *poslat;
  double        CP;
  double        weight_factor;

   /* initialize arrays for insert and delete at first call */
  if(rndtab==0) {       
    new     = (MYVEC*) calloc(polA+polB,sizeof(MYVEC));  
    box     = (int*)   calloc(polA+polB,sizeof(int));    
    /* stores box number of created monomers  */
  }

  /*---------------------------------------*
   | set variables according to chain-type |
   *---------------------------------------*/

  if (type==0) { 
    nrchains = nrchainsA;
    pos      = posA;
    poslat   = poslatA;
    pol      = polA;
    CP       = CPA;
 
    if (nrchainsA >= NPOLYMAXA) {
      fprintf(stderr,"ERROR in cbgc.c: nrchainsA>="
                     "NPOLYMAXA\n");
      return(0);
    }
  }    

  else if (type==1) { 
    nrchains = nrchainsB;
    pos      = posB;
    poslat   = poslatB;
    pol      = polB;
    CP       = CPB;

    if (nrchainsB >= NPOLYMAXB) {
      fprintf(stderr,"ERROR in cbgc.c: nrchainsB>="
                     "NPOLYMAXB\n");
      return(0);
    }
  }

  else {
    fprintf(stderr,"Incorrect 'type' declaration in insert()." 
                   "EXITING program ...\n");  
    exit(1);
  }


  /*------------------------------------------*
   | initialize or update random number array |
   *------------------------------------------*/ 

  if (rndtab==NULL) {    /* first call -> initialization    */

    rndtab = (double*) calloc(NRND_MAX,sizeof(double));   
    /* allocate memory for NRND_MAX #s */

    double_r250_vector(rndtab,NRND_MAX);                  
    /* assign random #s to array       */

    rnd_ptr = rndtab;                                     
    /* set rnd_ptr to first position   */
  }

  acc = (int)(rnd_ptr-rndtab);
  if ( (acc<0)||(acc>NRND_MAX) ) {
    fprintf(stderr,"ERROR: must generate too many RNDs "
                   "in cbgc.c (%d)\n",acc);
  }

  else {              /* update random number array    */
    double_r250_vector(rndtab,acc);                       
    /* refresh rnd #s which have already been used     */

    rnd_ptr = rndtab; /* set rnd_ptr to first position */
  }


  /*-------------------------------------------------*
   | step 1: create first monomer at random position |
   *-------------------------------------------------*/

  ident = nrchains*pol;                       
  /* assign first free position to first monomer */

  mono  = new;                                

  mono->x = LSX*(*rnd_ptr++); /* create at random position */
  mono->y = LSY*(*rnd_ptr++);
  mono->z = LSZ*(*rnd_ptr++);

  energy_new = elj_energie(ident,type,mono);  
  /* update energy of new chain */  

#ifdef HARDWALLS
  energy_new += wenergy(mono);
#endif

  znew    = exp(Tinv*energy_new);        
  /* store Boltzmann-weight of first monomer          */

  VCOPY(mono,pos+ident);                 
  /* copy monomer vector to last pos. of pos-array    */     

  VCOPY(mono,poslat+ident);              
  /* copy monomer vector to last pos. of poslat-array */

  box[0]  = VBOXNR(poslat+ident);

  if (type==0) {
  
    /* update double-linked subbox list */
    put_element(box[0],ident);            
    
    /* update elementliste-array        */
    (element_liste+ident)->monomer = ident;
    (element_liste+ident)->type    = 0; 

    nrmonA++;    /* increment # of monomers by one    */
  }

  else {
    put_element(box[0],ident+NMONOMAXA);
    (element_liste+NMONOMAXA+ident)->monomer = ident;
    (element_liste+NMONOMAXA+ident)->type = 1;
    nrmonB++;
  }
                                           
  i_break=pol; 


  /*----------------------------------*
   | step 2: create rest of the chain |
   *----------------------------------*/

  for (i_mono=1;i_mono<pol;i_mono++) {        
  /* loop to create rest of the chain */
    length = choose_length(type,*rnd_ptr++);  
    /* choose bondlength */


    
    /*---------------------------------------------*
     | 2.1.: create W_ALL trial-monomers on sphere |
     |       around first monomer                  |
     *---------------------------------------------*/

    for (i=0;i<W_ALL;i++) {                     

      /* determine trial-monomer pos. on sphere */
      /* trial-vectors are evenly distributed ! */
      ctheta = 2*(*rnd_ptr++) - 1;              
      stheta = 2*(*rnd_ptr++) - 1;              
      phi    = ctheta*ctheta + stheta*stheta;
      while (phi >= 1 || phi == 0) {
        ctheta = 2*double_r250() - 1;
        stheta = 2*double_r250() - 1;
        phi    = ctheta*ctheta + stheta*stheta;
      }
      (trial+i)->z = mono->z + length*(1.0 - 2*phi);
      phi    = 2*sqrt(1-phi)*length;
      (trial+i)->x = mono->x + phi*ctheta;      
      (trial+i)->y = mono->y + phi*stheta;
      /* trial vectors are stored in trial[]    */
    }

    /* set up interaction list to speed up      */ 
    /* energy calculations                      */
    /* remark: start-monomer "mono" is not part */ 
    /* of the interaction list                  */
    make_list(ident,type,mono,RANGE+length+0.05);  



    /*--------------------------------------------------*
     | 2.2.: choose one of the trial vectors according  |
     |       to its Boltzmann-weight                    |
     *--------------------------------------------------*/

    for (i=0,z=0;i<W_ALL;i++) {                 
      energy[i]  = lj_energy(type,trial+i);    
      /* calculate LJ energy of each trial vector */

#ifdef HARDWALLS  
      energy[i] += wenergy((trial+i));     
#endif

      z         += exp(Tinv*energy[i]);    
      /* calc. contribution to Rosenbluth-weight  */
      /* z is also used as a normalization factor */
      prob[i]    = z;
    }
 
    prob[W_ALL]=2;  
    /* always accept this W_ALL+1 possibility (z==0) ,*/ 
    /* but write warning */ 

    znew *= z;  /* update Rosenbluth-weight for MC step */

    i_chosen = 0;

    
    /*-----------------------------------------------------*
     | Each trial vector is assigned a probabilty of being |
     | chosen proportional to its Boltzmann-factor.        |
     | A random number is drawn and normalized. The while  |
     | loop determines to which trial vector the number    |
     | corresponds                                         |
     *-----------------------------------------------------*/ 

    prn      = z*(*rnd_ptr++);             
    /* draw random number and normalize it */

    while(prn>prob[i_chosen]) i_chosen++;  
    /* choose trial vector                 */  
                                           
    if (i_chosen>=W_ALL) {
      fprintf(stderr,"ERROR: run out of trials => rnd==1 "
                     "or z=%lg  resolved to i_chosen=0 and "
                     "break\n",z);
      i_chosen = 0;
      znew     = 0;
      i_break  = i_mono; 
      break;
    }

    
    /*---------------------*
     | 2.3.: update arrays |
     *---------------------*/

    energy_new += energy[i_chosen]+bondenergy(type,length);  
    /* update energy of new chain */

    mono++;
    ident++;

    x = (trial+i_chosen)->x;
    y = (trial+i_chosen)->y;
    z = (trial+i_chosen)->z;

    mono->x           = x;
    mono->y           = y;
    mono->z           = z;

    (pos+ident)->x    = x;        /* update pos-array    */
    (pos+ident)->y    = y;
    (pos+ident)->z    = z;

    x = x - LSX*((int)(x/LSX));   /* update poslat-array */
    if(x<0) x += LSX; 
    y = y - LSY*((int)(y/LSY)); 
    if(y<0) y += LSY; 
    z = z - LSZ*((int)(z/LSZ)); 
    if(z<0) z += LSZ; 

    (poslat+ident)->x = x;         
    (poslat+ident)->y = y;
    (poslat+ident)->z = z;

    box[i_mono]       = VBOXNR(poslat+ident);

    if (type==0) {

      /* update double-linked subbox list */
      put_element(box[i_mono],ident);              

      /* update element_liste             */
      (element_liste+ident)->monomer = ident;     
      (element_liste+ident)->type    = 0;
      nrmonA++;      /* increment # of monomers by one   */
    }

    else {
      put_element(box[i_mono],ident+NMONOMAXA);   
      (element_liste+NMONOMAXA+ident)->monomer = ident;
      (element_liste+NMONOMAXA+ident)->type = 1;
      nrmonB++;
    }
                                    
  } /* a new chain has been created ! */


  /*----------------------*
   | step 3: cbgc MC step |   
   *----------------------*/   

  if (type==0)                                                
    weight_factor=w[nrchainsA][nrchainsB]-\
                  w[nrchainsA+1][nrchainsB]; 
  else 
    weight_factor=w[nrchainsA][nrchainsB]-\
                  w[nrchainsA][nrchainsB+1]; 


  if (znew*exp( -Tinv*CP + weight_factor )\
              /(nrchains+1) > (*rnd_ptr++) ) {
  /* MC move is accepted ! */

    cenergie += energy_new;        
    if (type==0) nrchainsA++;
    else         nrchainsB++;           
   
    return(1);
  }

  else {  /* MC move is rejected */
    for (i_mono=0,ident=nrchains*pol;i_mono<i_break;\
         i_mono++,ident++) {

    /* delete monomers from boxes */
    if (type==0) del_element(box[i_mono],ident);      
    else         del_element(box[i_mono],ident+NMONOMAXA);
    }    
 
    if (type==0) nrmonA -= i_break;
    else         nrmonB -= i_break;        
    return(0);
  }

}


/**********************************************
 * delete: attemps to delete a chain from     *
 *         simulation box                     *
 *                                            *
 * type:   A (=0) or B (=1)                   *
 **********************************************/

int delete(int type) {
  static double *rndtab = NULL;
  static double *rnd_ptr;
  static int    *box;
  MYVEC         *mon0,trial[W_ALL];
  int           i_poly,i_mono,ident,ident2,box1,box2,i,acc;
  double        length,ctheta,stheta,phi,x,y,z,e;

  double        zold,energy_old;        
  /* zold = Rosenbluth-weight                   */
  /* energy_old = energy of chain to be deleted */

  int           nrchains, pol;
  MYVEC         *pos, *poslat;
  double        CP;
  double        weight_factor;

  if(rndtab==0) {                                       
    box     = (int*)calloc(polA+polB,sizeof(int));    
    /* stores box number of created monomers  */
  }


  /*---------------------------------------*
   | set variables according to chain-type |
   *---------------------------------------*/

  if (type==0) { 
    nrchains = nrchainsA;
    pos      = posA;
    poslat   = poslatA;
    pol      = polA;
    CP       = CPA;
 
    if (nrchains >= NPOLYMAXA) {
      fprintf(stderr,"ERROR in cbgc.c: nrchainsA=NPOLYMAXA\n");
      return(0);
    }
  }

  else if (type==1) { 
    nrchains = nrchainsB;
    pos      = posB;
    poslat   = poslatB;
    pol      = polB;
    CP       = CPB;

    if (nrchains >= NPOLYMAXB) {
      fprintf(stderr,"ERROR in cbgc.c: nrchainsB=NPOLYMAXB\n");
      return(0);
    }
  }

  else {
    fprintf(stderr,"Incorrect 'type' declaration in delete(). "
                   "EXITING program ...\n");  
    exit(1);
  }

  /*------------------------------------------*
   | initialize or update random number array |
   *------------------------------------------*/ 

  if (rndtab==NULL) {  /* first call -> initialization    */

    rndtab = (double*) calloc(NRND_MAX,sizeof(double));   
    /* allocate memory for NRND_MAX #s */

    double_r250_vector(rndtab,NRND_MAX);                  
    /* assign random #s to array       */

    rnd_ptr = rndtab;                                     
    /* set rnd_ptr to first position   */
  }

  acc = (int)(rnd_ptr-rndtab);
  if ( (acc<0)||(acc>NRND_MAX) ) {
    fprintf(stderr,"ERROR: must generate too many RNDs "
                   "in cbgc.c (%d)\n",acc);
  }

  else {    /* update random number array      */

    double_r250_vector(rndtab,acc);                       
    /* refresh rnd #s which have already been used  */

    rnd_ptr = rndtab;  /* set rnd_ptr to first position   */
  }


  /*----------------------------------------------*
   | step 1.1: choose polymer at random which fct |
   |           will attempt to delete             | 
   *----------------------------------------------*/

  if (nrchains==NPOLYBRUSH) return(0);
  else i_poly = NPOLYBRUSH+\
                r250_random()%(nrchains-NPOLYBRUSH);  
  ident = i_poly*pol;    
  /* chosen polymer's first monomer position */
  /* in pos and poslat arrays                */

 
  /*-------------------------------------------------------*
   | step 1.2: swap positions with last polymer:           |
   |           newly inserted polymers are always stored   |
   |           at the end of pos, poslat and elementliste- |
   |           arrays                                      |
   |            -> chosen polymer must be swaped with      |
   |               last polymer                            |
   *-------------------------------------------------------*/ 

  if (i_poly!=nrchains-1) {                    
  /* chosen polymer is not last polymer !       */

    ident2 = (nrchains-1)*pol;                 
    /* position of last polymer in pos and poslat */

    for (i_mono=0;i_mono<pol;i_mono++) {           
    /* swap positions of all monomers with        */
    /*  monomers of last polymer                  */

       box1 = VBOXNR(poslat+(ident +i_mono));  
       box2 = VBOXNR(poslat+(ident2+i_mono));


       if (type==0) {
       /* swap positions in double-linked subbox lists */
       del_element(box1,ident+ i_mono);       
       put_element(box2,ident+ i_mono);        
       del_element(box2,ident2+i_mono);
       put_element(box1,ident2+i_mono);
       }

       else {
       del_element(box1,ident+ i_mono + NMONOMAXA);       
       put_element(box2,ident+ i_mono + NMONOMAXA);        
       del_element(box2,ident2+i_mono + NMONOMAXA);
       put_element(box1,ident2+i_mono + NMONOMAXA);
       }       

       x = (pos+(ident+ i_mono))->x;   
       /* swap positions in pos-array */
       (pos+(ident+ i_mono))->x = (pos+(ident2+i_mono))->x;
       (pos+(ident2+i_mono))->x = x;
       y                        = (pos+(ident+ i_mono))->y;
       (pos+(ident+ i_mono))->y = (pos+(ident2+i_mono))->y;
       (pos+(ident2+i_mono))->y = y;
       z                        = (pos+(ident+ i_mono))->z;
       (pos+(ident+ i_mono))->z = (pos+(ident2+i_mono))->z;
       (pos+(ident2+i_mono))->z = z;

       /* swap positions in poslat-array */
       x = (poslat+(ident +i_mono))->x; 
       (poslat+(ident+ i_mono))->x =
               (poslat+(ident2+i_mono))->x;
       (poslat+(ident2+i_mono))->x = x;

       y =(poslat+(ident+ i_mono))->y;
       (poslat+(ident+ i_mono))->y =
               (poslat+(ident2+i_mono))->y;
       (poslat+(ident2+i_mono))->y = y;

       z = (poslat+(ident+ i_mono))->z;
       (poslat+(ident+ i_mono))->z = 
               (poslat+(ident2+i_mono))->z;
       (poslat+(ident2+i_mono))->z = z;

    } /* swap process succesful ! */
  } /* chosen polymer is not last polymer !*/

  
  /*-----------------------------*
   | step 2: delete chosen chain |
   *-----------------------------*/ 

  ident = (nrchains-1)*pol;
  zold=1;
  energy_old=0;

  /* delete all monomers of chain except first monomer  */
  for (i_mono=pol-2;i_mono>=0;i_mono--) {              
    box[i_mono+1] = VBOXNR(poslat+(ident+i_mono+1));

    if (type==0) {
    del_element(box[i_mono+1],(ident+i_mono+1));       
    /* update double-linked subbox list */
    nrmonA--;                                          
    /* decrement # of monomers by one   */
    }
   
    else {
    del_element(box[i_mono+1],(ident+i_mono+1+NMONOMAXA));
    nrmonB--;
    }

    mon0   = pos+(ident+i_mono);

    length =  DISTANCE(mon0,(mon0+1));
    /* determine bondlength */

    energy_old += bondenergy(type,length); 
    /* update total energy  */

    VCOPY(mon0+1,trial);
    /* copy monomer vector to first position of trial array */
  

    /* create W_ALL-1 trial vectors   */
    /* trial vectors are evenly       */
    /* distributed on a sphere around */
    /* mon0                           */   
    /* trial spheres are created for  */
    /* all monomers except last !     */

    for (i=1;i<W_ALL;i++) {                            
       ctheta = 2*(*rnd_ptr++) - 1;
       stheta = 2*(*rnd_ptr++) - 1;                    
       phi    = ctheta*ctheta + stheta*stheta;         
       while (phi >= 1 || phi == 0) {
         ctheta = 2*double_r250() - 1;                 
         stheta = 2*double_r250() - 1;                
         phi    = ctheta*ctheta + stheta*stheta;
       }
      (trial+i)->z = mon0->z + length*(1.0 - 2*phi);
      phi    = 2*sqrt(1-phi)*length;
      (trial+i)->x = mon0->x + phi*ctheta;
      (trial+i)->y = mon0->y + phi*stheta;
    }

    make_list(ident+i_mono,type,mon0,RANGE+length+0.05);    
    /* set up interaction list         */ 
    /* to speed up energy calculations */               
                                                         
    for (i=0,z=0;i<W_ALL;i++) {
      e  = lj_energy(type,trial+i);                    
     /* calc. LJ energy of each trial vector */

#ifdef HARDWALLS
      e += wenergy((trial+i)); 
#endif

      if (i==0) {                                 
      /* monomer belongs to chosen chain  */

      energy_old += e;  /* add energy to total energy       */
      }

      z += exp(Tinv*e);
    }

    zold *= z;  /* update Rosenbluth-weight */
  }/* all but first monomer have been deleted */

  if (type==0) nrmonA--;
  else         nrmonB--;

  e = elj_energie(ident,type,poslat+ident);      
  /* calc. energy of first monomer */

#ifdef HARDWALLS
  e += wenergy((poslat+ident));
#endif

  energy_old += e;     /* total energy of removed chain */
  zold *=exp(Tinv*e);  /* update Rosenbluth-weight      */
  box[0] = VBOXNR(poslat+(ident));     

  /* delete first monomer of chain */
  if (type==0)  del_element(box[0],ident);    
  else          del_element(box[0],ident+NMONOMAXA);


  /*----------------------*
   | step 3: cbgc MC step |
   *----------------------*/

  if (type==0)
    weight_factor=w[nrchainsA][nrchainsB]-\
                  w[nrchainsA-1][nrchainsB];
  else
    weight_factor=w[nrchainsA][nrchainsB]-\
                  w[nrchainsA][nrchainsB-1];    

  if ( nrchains/(zold*exp(-Tinv*CP - weight_factor) )\
        > (*rnd_ptr++) ) {   /* MC move is accepted */

    cenergie -= energy_old;                    
    if(type==0)  nrchainsA--;
    else         nrchainsB--;                                
    return(1);
  }

  else {  /* MC move is rejected */
    for (i_mono=0;i_mono<pol;i_mono++) {       
      if (type==0)  put_element(box[i_mono],ident+i_mono);
      else  put_element(box[i_mono],ident+i_mono+NMONOMAXA);
    }
    if (type==0)  nrmonA+=pol;
    else          nrmonB+=pol;
    return(0);
  }

}/* delete */
