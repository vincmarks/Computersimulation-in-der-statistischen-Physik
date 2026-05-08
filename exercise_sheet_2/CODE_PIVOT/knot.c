/*********************************************************
 * knot.c: functions to determine knots in configuration *
 *                                                       *
 * last modification: 03/25/2004                         *
 *********************************************************/


#include "element.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern int     mcs, polA, nrmonA;
extern MYVEC   *posA;
extern int     nrmon_knotA, nrmon_knotB;
extern MYVEC   *knotA, *knotB;
extern DATA    data;

/************************************************************
 * check_for_knots(): - reduces complexity of configuration *
 *                      by calling reduce_complexity        *
 *                    - analyzes remaining structure with   *
 *                      knot polynomials                    *
 *                                                          *    
 * returns:	      YES: found a knot                     *
 * 		      NO : didn't find a knot               *
 *							    *		
 * last modification: 11/30/2004                            *
 ************************************************************/


  int check_for_knots(){

  int i;
  MYVEC *mon, *mon_knotB;
  char filenameout[30];
  double t1,t2;
  double polynomial1, polynomial2, polynomial3, polynomial4, polynomial5, polynomial6, polynomial7, polynomial8;
  int number_of_bonds_left;
 
  /* allocate memory for arrays which are needed to analyze knots */
  knotA           = (MYVEC *) calloc(NMONOMAXA+10, sizeof(MYVEC));
  knotB           = (MYVEC *) calloc(NMONOMAXA+10, sizeof(MYVEC));
 
  copy_MYVEC_array(posA,knotB,0,nrmonA);
  nrmon_knotB=nrmonA;
  closure(0);
  number_of_bonds_left=reduce_complexity_2();
  t1=-1.1;
  t2=1/t1;
  /* polynomial1  = alexander_polynomial(-1); */
  polynomial2 = alexander_polynomial(t1)*alexander_polynomial(t2);
 
  /* copy_MYVEC_array(posA,knotB,0,nrmonA);
  nrmon_knotB=nrmonA;
  closure(2);
  number_of_bonds_left=reduce_complexity_2();
  t1=-1.1;
  t2=1/t1;
  polynomial3  = alexander_polynomial(-1);
  polynomial4 = alexander_polynomial(t1)*alexander_polynomial(t2); */
 
  copy_MYVEC_array(posA,knotB,0,nrmonA);
  nrmon_knotB=nrmonA;
  number_of_bonds_left=reduce_complexity_2();
  closure(0);
  t1=-1.1;
  t2=1/t1;
  /* polynomial5  = alexander_polynomial(-1); */
  polynomial6 = alexander_polynomial(t1)*alexander_polynomial(t2);
 
  /* copy_MYVEC_array(posA,knotB,0,nrmonA);
  nrmon_knotB=nrmonA;
  number_of_bonds_left=reduce_complexity_2();
  closure(2);
  t1=-1.1;
  t2=1/t1;
  polynomial7  = alexander_polynomial(-1);
  polynomial8 = alexander_polynomial(t1)*alexander_polynomial(t2); */

  printf("%.5lf \t%.5lf \n",polynomial2,polynomial6);

  free(knotA);
  free(knotB);

  if(sqrt((polynomial2-1)*(polynomial2-1)) > 0.000001 ||
     sqrt((polynomial6-1)*(polynomial6-1)) > 0.000001 ){
   
     return(YES);
  }

  else { 
     return(NO);
  }
}



/****************************************************************
 * check_for_knots_1(): as above, but calls reduce_complexity_1 *
 *                                                              *    
 * returns:	      YES: found a knot                         *
 * 		      NO : didn't find a knot                   *
 *							        *		
 * last modification: 05/31/2004                                *
 ****************************************************************/

int check_for_knots_1()
{
  /* allocate memory for arrays which are needed to analyze knots */
  knotA           = (MYVEC *) calloc(NMONOMAXA+1, sizeof(MYVEC));
  knotB           = (MYVEC *) calloc(NMONOMAXA+1, sizeof(MYVEC));

  if(reduce_complexity_1()==2){  /* if (# of remaining monomers in reduced configuration == 2) */
    free(knotA);               /* -> no knots ! */
    free(knotB);
    return(NO);              
  }

  else{                       /* else: check remaining configuration */                       
    alexander_polynomial(-1);    
    free(knotA);
    free(knotB);
    return(YES);   /* will eventually become more sophisticated ! */
  }
}



/*********************************************************************************
 * determine_knot_size: dito                                                     *
 *									         *
 * input:   p value for the product of the alexander polynomials (A(t1)*A(1/t1)) *
 * returns: knot size (in monomer units)                                         *
 *                                                                               *
 * last modification: 11/30/2004                                                 *
 *********************************************************************************/

int determine_knot_size(double p)
{
  double polynomial;
  int i=0;
  double t1,t2;
  int start=0,end=0;
  
  char filenameout[30]="conf_red";

  /*---------------------------------------------------------------*
   | The first and the last particle are always part of the chain. |
   | The programs starts with one end of the chain and removes the |
   | monomer next to the end particle. Then it applies Taylor's    |
   | reduction and calculates the product of the Alexander         |
   | polynomials at value t1 and 1/t1 (A(t1)*A(1/t1)).             |
   | If the structure stays the same, variable "start" is set to   |
   | the monomer number which has been removed                     |
   | The iteration stops when the structure becomes an unknot.     |
   | The proceedure is repeated from the other side of the chain.  |
   | This time variable "end" is determined.                       |
   | The size of the knot is defined as "end"-"start".             |
   *---------------------------------------------------------------*/

  t1=-1.1;      /* t1 and t2 should be the same as in check_for_knots */
  t2=1/t1; 
  polynomial=p;

  /*--------------------------------*
   | 1. reduce chain from beginning |
   *--------------------------------*/

  for(i=1;polynomial-1>0.000001;i++){
    copy_MYVEC_array_2(posA,knotB,i,nrmonA-1);  /* prepare reduced structure */
    nrmon_knotB=nrmonA-i+1;

    closure(0);
    reduce_complexity_2();                      /* apply Taylor reduction    */
    /* closure(10); */                               /* first and last particle are always part of the chain */
    polynomial=alexander_polynomial(t1)*alexander_polynomial(t2);  /* determine alexander polynomial */       

    if( sqrt((polynomial-p)*(polynomial-p))  < 0.000001) {   /* set variable "start" if "polynomial==p" */
      start=i;
    }

    /* printf("%d %d %.5lf\n",start,0,polynomial); */
  }

  polynomial=p;

  /*--------------------------* 
   | 2. reduce chain from end | 
   *--------------------------*/

  for(i=1;polynomial-1>0.000001;i++){
    copy_MYVEC_array_2(posA,knotB,start,nrmonA-i);
    nrmon_knotB=nrmonA-(start)-i+2;

    closure(0);
    reduce_complexity_2();
    /* closure(10); */
    polynomial=alexander_polynomial(t1)*alexander_polynomial(t2);

    if( sqrt((polynomial-p)*(polynomial-p)) < 0.000001) {
      end=nrmonA-i;
    }

    /* printf("%d %d %.5lf\n",start,end,polynomial); */
  } 

  /*------------------------------------------*
   | 3. print out final reduced configuration |
   *------------------------------------------*/ 

  /* copy_MYVEC_array_2(posA,knotB,start,end);
  nrmon_knotB=end-start+2;
  reduce_complexity_2();
  closure(10); 

  sprintf(filenameout,"c_0");
  sysout(filenameout,knotB,nrmon_knotB,mcs); */


  printf("%d\t%d\t",start,end);
  return(end-start);
}


void determine_bond_size_distribution() 
{
  int i;
  double distance;

for(i=0;i<nrmon_knotB-1;i++) {
    printf("1\t%lf\n",sqrt( ((knotB+i)->x - (knotB+i+1)->x) * ((knotB+i)->x - (knotB+i+1)->x)
                  + ((knotB+i)->y - (knotB+i+1)->y) * ((knotB+i)->y - (knotB+i+1)->y) ));
    printf("2\t%lf\n",DISTANCE(knotB+i+1,knotB+i)); 
    printf("3\t%d\n",(knotB+i+1)->monomer_number - (knotB+i)->monomer_number); 
   }
}


/***************************************************************
 * copy_MYVEC_array: copies array_a of struct MYVEC to array_b *
 *                   starting at position "start" and ending   *
 *                   at position "end"                         *
 *                                                             *
 * last modification: 09/29/2004                               *
 ***************************************************************/

void copy_MYVEC_array (MYVEC *array_A, MYVEC *array_B, int start, int end)
{
int i,j=0;

  for(i=start;i<end;i++,j++) {
    VCOPY(array_A+i,array_B+j);
    /* printf("%d\t",j);
    VPRINT(array_B+j); */
  }
}

void copy_MYVEC_array_2 (MYVEC *array_A, MYVEC *array_B, int start, int end)
{
int i,j=0;
                       
  VCOPY(posA,array_B);
  j++;  
                                                                                                    
  for(i=start;i<end;i++,j++) {
    VCOPY(array_A+i,array_B+j);
  }

  VCOPY(posA+nrmonA-1,array_B+j);
  (array_B+j)->monomer_number=(posA+nrmonA-1)->monomer_number; 
}


/********************************************************************
 * closure: closes open polymer such that ring is formed            *
 *                                                                  *
 * int method: 0 close loop by connecting first and last particle   *
 *             1 " " " " first and last particle to infinity        *
 *               (along x-axis)                                     *
 *             2 like 1 but connection along line through endpoints *
 *                                                                  *
 * last modification: 10/12/2004                                    *
 ********************************************************************/

void closure (int method)
{
MYVEC mon[1],vhelp[1];
double large_value1, large_value2;

  if(method==0) {
    VCOPY(posA,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+1;
    nrmon_knotB++;
  }

  if(method==10) {
    VCOPY(knotB,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+1;
    nrmon_knotB++;
  }
 
  if(method==1) {

    if( (knotB+nrmon_knotB-1)->x > knotB->x) {
    large_value1=100000000;
    large_value2=100000000;
    }
    else {
    large_value1=-100000000;
    large_value2=100000000;
    }

    mon->x = ((knotB+nrmon_knotB-1)->x)+large_value1;
    mon->y = ((knotB+nrmon_knotB-1)->y);
    mon->z = ((knotB+nrmon_knotB-1)->z); 
       
    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+1;
    nrmon_knotB++;

    mon->x = ((knotB+nrmon_knotB-2)->x)+large_value1;
    mon->y = ((knotB+nrmon_knotB-2)->y)-large_value2;
    mon->z = ((knotB+nrmon_knotB-2)->z);
                                                                                                                             
    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+2;
    nrmon_knotB++;

    mon->x = ((knotB)->x)-large_value1;
    mon->y = ((knotB)->y)-large_value2;
    mon->z = ((knotB)->z);
                                                                                                                             
    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+3;
    nrmon_knotB++;

    mon->x = ((knotB)->x)-large_value1;
    mon->y = ((knotB)->y);
    mon->z = ((knotB)->z);
                                                                                                                             
    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+4;
    nrmon_knotB++;

    mon->x = ((knotB)->x);
    mon->y = ((knotB)->y);
    mon->z = ((knotB)->z);
                                                                                                                            
    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+5;
    nrmon_knotB++;
  }

  if(method==2){
  
    large_value1 = 100000000;
    large_value2 = large_value1/DIST(knotB+nrmon_knotB-1,knotB);

    mon->x = knotB->x + large_value2 * ((knotB+nrmon_knotB-1)->x - knotB->x);
    mon->y = knotB->y + large_value2 * ((knotB+nrmon_knotB-1)->y - knotB->y);
    mon->z = knotB->z + large_value2 * ((knotB+nrmon_knotB-1)->z - knotB->z);

    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+1;
    nrmon_knotB++;

    DIFF(knotB+nrmon_knotB-2,knotB,mon); 
    VDOT(posA+1,mon,vhelp);
    mon->x = knotB->x + large_value1 * (vhelp->x)/BETRAG(vhelp);
    mon->y = knotB->y + large_value1 * (vhelp->y)/BETRAG(vhelp);
    mon->z = knotB->z + large_value1 * (vhelp->z)/BETRAG(vhelp);  
  
    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+2;
    nrmon_knotB++;

    mon->x = knotB->x - large_value2 * ((knotB+nrmon_knotB-3)->x - knotB->x);
    mon->y = knotB->y - large_value2 * ((knotB+nrmon_knotB-3)->y - knotB->y);
    mon->z = knotB->z - large_value2 * ((knotB+nrmon_knotB-3)->z - knotB->z);

    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+3;
    nrmon_knotB++;
 
    VCOPY(knotB,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+4;
    nrmon_knotB++;
  }
 
  if(method==3){

    calc_cm(nrmon_knotB, knotB, vhelp);  /* save center of mass in vhelp */
    large_value1 = 100000000;
  
    large_value2 = large_value1/DIST(vhelp,knotB+nrmon_knotB-1);
    mon->x = vhelp->x + large_value2 * ((knotB+nrmon_knotB-1)->x - vhelp->x);
    mon->y = vhelp->y + large_value2 * ((knotB+nrmon_knotB-1)->y - vhelp->y);
    mon->z = vhelp->z + large_value2 * ((knotB+nrmon_knotB-1)->z - vhelp->z); 

    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+1;
    nrmon_knotB++;

    large_value2 = large_value1/DIST(vhelp,knotB);
    mon->x = vhelp->x + large_value2 * (knotB->x - vhelp->x);
    mon->y = vhelp->y + large_value2 * (knotB->y - vhelp->y);
    mon->z = vhelp->z + large_value2 * (knotB->z - vhelp->z);
 
    VCOPY(mon,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+2;
    nrmon_knotB++;

    VCOPY(knotB,knotB+nrmon_knotB);
    (knotB+nrmon_knotB)->monomer_number=nrmonA+3;
    nrmon_knotB++;
  }
}
