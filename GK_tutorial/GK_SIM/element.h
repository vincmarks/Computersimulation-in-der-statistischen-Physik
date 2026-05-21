/*************************************************************
 * element.h: header file                                    *
 *            contains substitution macros, function macros, *
 *            definitions of data structures and function    *
 *            prototypes                                     *
 *                                                           *
 * last modification : 23/05/2003                            *
 *************************************************************/


/*#define PAR*/
/*#define SWITCH*/
/*#define NOSYNC*/
/*#define PRESSURE*/
/*#define T3D */
/*#define FIXLENGTH 0.974732*/
/*#define STRIPE     1 */
/*#define HARDWALLS*/
/*#define REPULSIVE*/

#define PRESSURE

/* #define TEST  */
#define NPOLYMAXA  500    /* A: POLYMER (N=5) */ 
#define NMONOMAXA  500
#define NPOLYMAXB  2    /* B: CO2 (N=1)     */
#define NMONOMAXB  2

/* number of trial vectors in ccb/rep */
#define W_ALL      25  

/* maximal number of monomers in wwlist*/
#define MAX_ANZAHL (NMONOMAXA+NMONOMAXB) 

#define NRSUBX     3      /* defines box size      */
#define NRSUBY     3      /* must be larger than 2 */      
#define NRSUBZ     3     

/* #define LARGE */   
/* speeds up calculations for box size>5 */
/* cannot be used for system sizes<4     */                   

#define MINA 0            /* defines simulation interval */
#define MAXA 500

#define INFENERGY  200
#define JwallN_MAX 1024
#define NPOLYBRUSH 0

#define NBOXX    4
#define NBOXY    4




/*--------------------------------------*
 | define energy interaction parameters |
 | (EPSILONA and SIGMAA :=1 !)          |
 *--------------------------------------*/

/* Epsilon: depth of LJ potential */
#define EPSILONB   0.726         /* B:  CO2 (N=1)   */
#define EPSILONAB  0.754921913   /* AB: POLYMER-CO2 */

/* Sigma: zero of LJ-potential   */
#define SIGMAB     0.816   
#define SIGMAAB    0.908         
#define SIGMAB6    0.295216721   /* Sigma^6 */
#define SIGMAAB6   0.56042189


#define OK       0
#define ERROR    1
#define TRUE     1
#define FALSE    0
#define OUT     -1
#define MEASURE  1

#ifdef PAR
#define NOSYNC
#endif

# ifndef T3D
# define PRINTF     printf
# define FPRINTF    fprintf
# define PRINTFS    printf
# define FPRINTFS   fprintf
# else
# define PRINTF     printf
# define FPRINTF    fprintf
# define PRINTFS    printf
# define FPRINTFS   fprintf
# endif

/*-----------------------------------------------------------*/

/*------------------------*
 | Makros fuer Ringlisten |
 *------------------------*/

#define BOXNR(a,b,c) ((a)+NRSUBX*(b)+NRSUB2*(c))
#define IFLOAT(a,L) ((int)((a)/L))

#define FBOXNR(ax,ay,az) BOXNR(IFLOAT(ax,LSSUBX),\
IFLOAT(ay,LSSUBY),IFLOAT(az,LSSUBZ))

#define VBOXNR(a) FBOXNR((a)->x,(a)->y,(a)->z)


/*------------------------------*
 | Makros zur Energieberechnung |
 *------------------------------*/

#define RNULL2INV (1.0/SQUARE(1.5))    

/* Offset fuer Bondenergie, diese hat Nullpunkt bei ca. 20 */
#define BOND_OFFSET (-20.0)         

#define BOND(x) (-33.75*log(1.0- (x)*RNULL2INV))   

/* sets LJ_potential = 0 at RANGE */
#define NULLPUNKT (0.007751464836) 
  
#define LENJAA(x)   ((SQUARE(x)-x)+NULLPUNKT)     
#define RANGEAA     (2.0*1.12246204830937298)           
#define RANGE2AA    SQUARE(RANGEAA)

#define LENJBB(x)   EPSILONB*(SQUARE(SIGMAB6*x)-\
                             (SIGMAB6*x)+NULLPUNKT)

#define RANGEBB     (2.0*1.12246204830937298*SIGMAB)           
#define RANGE2BB    SQUARE(RANGEBB)

#define LENJAB(x)   EPSILONAB*(SQUARE(SIGMAAB6*x)-\
                              (SIGMAAB6*x)+NULLPUNKT)
#define RANGEAB     (2.0*1.12246204830937298*SIGMAAB)         
#define RANGE2AB    SQUARE(RANGEAB)
#define RANGE       (MAX(MAX(RANGEAA,RANGEBB),RANGEAB))

/* maximale Verschiebung fuer local displacement*/
#define MAX_STEP    0.3 

#define NRSUB2  (NRSUBX*NRSUBY)
#define NRSUB3  (NRSUBX*NRSUBY*NRSUBZ)


/*--------------------------------------*
 | Makros fuer mathematische Funktionen |
 *--------------------------------------*/

#define PI 3.141592
#define ABS(x) (( (x)>0 )? (x) : -(x) )
#define MAX(x,y) (( (x)>(y) )? (x) : (y) )
#define MIN(x,y) (( (x)<(y) )? (x) : (y) )
#define SQUARE(x) ( (x)*(x) )
#define CUBE(x) ( (x)*(x)*(x) )
#define ASGN(x) (( (x)>0 )? (-1.0) : (1.0) )


/*---------------------------------------------------* 
 | Makros fuer mathematische Funktionen mit Vektoren |
 *---------------------------------------------------*/ 

#define BETRAG(a)    (sqrt(SQUARE((a)->x)+SQUARE((a)->y)+\
                                          SQUARE((a)->z)))

#define BETRAG2(a)   (SQUARE((a)->x)+SQUARE((a)->y)+\
                                     SQUARE((a)->z))

#define DOT(a,b)     (((a)->x)*((b)->x) + ((a)->y) *\
((b)->y)+ ((a)->z) * ((b)->z))

#define CWINKEL(a,b) ((DOT(a,b))/(BETRAG(a)*BETRAG(b)))

#define DIST(a,b)    (sqrt(SQUARE((a)->x - (b)->x)+\
                           SQUARE((a)->y - (b)->y)+\
                           SQUARE((a)->z - (b)->z)))

#define DISTANCE(a,b) (sqrt(SQUARE((a)->x - (b)->x)+\
                            SQUARE((a)->y - (b)->y)+\
                            SQUARE((a)->z - (b)->z)))
#define DISTANCE2(a,b) (SQUARE((a)->x - (b)->x)+\
                        SQUARE((a)->y - (b)->y)+\
                        SQUARE((a)->z - (b)->z))

#define VCOPY(a,b)     (b)->x = (a)->x;\
                       (b)->y = (a)->y;\
                       (b)->z = (a)->z;

#define DIFF(a,b,c)    (c)->x = (a)->x - (b)->x;\
                       (c)->y = (a)->y - (b)->y;\
                       (c)->z = (a)->z - (b)->z

#define ADD(a,b,c)     (c)->x = (a)->x + (b)->x;\
                       (c)->y = (a)->y + (b)->y;\
                       (c)->z = (a)->z + (b)->z

#define SDOT(n,a,b)    (b)->x = (a)->x*(n);\
                       (b)->y = (a)->y*(n);\
                       (b)->z = (a)->z*(n)

#define VDOT(a,b,c)    (c)->x = (a)->y*\
                                (b)->z-(a)->z * (b)->y;\
                       (c)->y = (a)->z*\
                                (b)->x-(a)->x * (b)->z;\
                       (c)->z = (a)->x*\
                                (b)->y-(a)->y * (b)->x



#define VPRINT(a) printf("%f %f %f\n",(a)->x,(a)->y,(a)->z)
     

/*----------------------*
 | Strukturdefinitionen |
 *----------------------*/

typedef struct o_element ELEMENT;
typedef struct o_box     BBOX;

typedef struct {
  double x;
  double y;
  double z;
} MYVEC;

struct o_element {
  int monomer;               /* monomer number */
  int type;                  /* monomer type   */
  struct o_element *before;  /* pointer to previous element  */
  struct o_element *after;   /* pointer to next element      */
};


/*                                             element 0     */
/*                                             after         */
/*                      element 1       <------before        */
/*                      after           ------>              */
/*  element 2     <-----before                               */
/*  after         ----->                                     */
/*  before                                                   */

struct o_box {
  struct o_element *first;  /* points somewhere in the list, */
                            /* there is no beginning         */
  int population;           /* number of list elements ,     */
                            /* number of particles in box    */
  int neighbours[125];      /* indices of neighboring boxes  */
                            /* - 27 inner and (125-27) outer */
};

typedef struct DATA {
  double bondangleA;
  double bondangleB;
  double bondlength1A;
  double bondlength1B;
  double bondlength2A;
  double bondlength2B;
  double bondlengthmaxA;
  double bondlengthmaxB;
  double gyrA;
  double gyrB;
  double endA;
  double endB;
  double energy_ljA;
  double energy_ljB;
  double energy_ljAB;
  double energy_feneA;
  double energy_feneB;
  double energy_sqA;
  double energy_sqB;
  double energy_sqAB;
  double energy_wall;
  double msdA;
  double msdB;
  double px;
  double py;
  double pz;
  double pwall1;
  double pwall2;
  double virialx;
  double virialy;
  double virialz;
  int    nrchainsleftA;
  int    nrchainsleftB;
} DATA;


#ifndef FUNKTIONEN
#define FUNKTIONEN

/*-----------------------------------------------------------*/

/*-----------------------*
 | Funktionendeklaration |
 *-----------------------*/


/*-------------------------------------*
 | element.c: box and energy functions |
 *-------------------------------------*/

void put_element(int ibox, int ident);
/* puts monomer in box */
             
void del_element(int ibox,int ident);                 
/* deletes element from box */  

void populate(void);                                  
/* fills boxes with monomers from starting configuration */

double elj_energie(int ident, int type, MYVEC *mon);  
/* calculates LJ interactions with BOX method */

double fene(int ident, int type, MYVEC *neu);  
/* calculates bond energy with all neighbors  */

double eenergie(void);                                
/* calculates total energy of the system using elj_energie,  */
/* fene, und wenergy, also checks box scheme                */

                    
/*----------------------------*
 | wach.c: analysis functions |
 *----------------------------*/ 

void analyse(int zeit,int lauf);                      
/* updates data-array by calling analysis functions 
  - writes results to "histo.dat" */
                                                      
void bondwinkel();   /* calc. average cos(bondangle)         */
void bondlength();   /* calc. average bondlength, b**2, bmax */
void calc_polymer(); /* calc. (R_e^2) and (R_g^2)            */

void calc_nrchainsleft(void);  
/* calc. #chains in the left half of sim. box */

void calc_cm(int pol, MYVEC *start,MYVEC *result);    
/* calc. center-of-mass of polymer at pos. start */

void calc_msd();     /* calc. mean-square displacement */
void calc_energy();  /* calc. total LJ, FENE and wall energy */

/* pressure.c */
int  calc_p(void);   /* calc. pressure via virial expression */


/*--------------------------------*
 | mcmove.c: functions for cbgc.c |
 *--------------------------------*/

double bondenergy(int type, double length);  
/* calc. LJ+FENE energy of two monomers (same type) */ 
/* at distance "length"                             */

double calcIntegral(double *, double *);              
/* updates rtable[] and itable[] in mcmove.c        */
/* - prerequesit for choose_length()                */     
   
double choose_length(int,double rnd);          
/* chooses Boltzmann-distributed bondlength         */

void   make_list(int, int, MYVEC *, double);   
/* creates interaction-arrays wwlisteA and wwlisteB */

double lj_energy(int, MYVEC *);                       
/* calculates LJ interaction of single monomer      */
/* with all monomers in wwlisteA+B                  */       


/*---------------------------------------------------*
 | cbgc.c: grand-canonical configurational-bias move |
 *---------------------------------------------------*/

double cbgc(int type);  
/* tries to insert (50%) or delete (50%) chains from sim box */

int insert (int type);
int delete (int type);


/*------------------------*
 | local.c: local MC move |
 *------------------------*/

void lmtrial(MYVEC *new, MYVEC *alt, double ** rnd_ptr_adr);  
/* creates trial vector for lmcmove */

int  lmcmove(void);                                           
/* tries nrmonA+nrmonB local moves, */ 
/* monomers are chosen at random    */


/*--------------------------*
 | rep.c: reptation MC move |
 *--------------------------*/

int reptation(void);              /* calls rep fct */
int rep      (int type);          
/* tries nrmonA (if type==0) or nrmonB (if type==1) */
/* reptation moves, monomers are chosen at random   */


/*--------------*
 | io-functions |
 *--------------*/

int sysout(char *,int );      
/* reads sim. specifications and monomer positions  */ 
/* from "conf_out"                                  */

int sysin(char *);            
/* writes sim. specifications and monomer positions */  
/* to "conf_out"                                    */

void init_weight(void);      
/* initializes weights for reweighted simulation    */


/*----------------*
 | wall functions |   -these functions have not been adjusted 
 *----------------*/

extern double (*wenergy)(MYVEC *pos);
extern double (*wbenergy)(MYVEC *pos);

double wall_energy(MYVEC *pos);   
/* calculates interaction energy with walls */

double bwall(MYVEC *pos);         
/* purely repulsive interaction SEwall */ 

double piston(MYVEC *pos);        
/* calculates interaction energy with walls and piston */

void   profile(int status);
#endif

/*
  
Erklaerungen (wichtigste Makros und Datenstrukturen)
====================================================

BOXNR        : Ermittelt aus (Integer) 
               Koordinaten die zugehoerige Boxnummer

IFLOAT       : Rechnet Fliesskommakoordinaten in 
               zugehoerige Boxkoordinaten um

BOXNR        : Wie BOXNR, akzeptiert aber 
               Fliesskommakoordinaten

BOX_PERIODIC : Ermittelt, mittels Funktion box_periodic, 
               fuer Teilchen die zugehoerige Boxennummer

BOND         : FENE-potential aus rad2 berechnen
LENJ         : Lenard Jones Potential aus rad6inv berechnen
DIST         : Abstand zwischen zwei Punkten
BETRAG       : Betrag eines (Pointer auf) Vektors
PRODUKT/DOT  : Vektorprodukt
CWINKEL      : Winkel zwischen zwei Vektoren
DISTANCE2    : Quadrat der Distanz zwischen zwei Punkten
BETRAG2      : Betragsquadrat eines Vektors


(o_)vector   : Datenstruktur fuer Vektoren
(o_)element  : Datenstruktur fuer Ringlisten, 
               deren Elemente auf Teilchen zeigen
(o_)box      : Datenstruktur fuer Boxen, 
               first zeigt auf ein Element der Ringliste,
               welche die Teilchen enthaelt, 
               die sich in der Box befinden

pos          : absolute positions of the monomers
poslat       : periodic boundary conditions applied
*/ 
