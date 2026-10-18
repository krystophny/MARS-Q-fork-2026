C*COMDECK COMCON
C
***********************************************************************
*                                                                     *
*     C.04         CONTROL VARIABLES                                  *
*                                                                     *
*     NANAL  : ANALYTIC SOLUTION FOR SOLOVEV CASE (NSURF = 1)         *
*              NANAL = 0 -------> NUMERICAL SOLOVEV EQUILIBRIUM       *
*              NANAL = 1 -------> ANALYTICAL SOLOVEV EQUILIBRIUM      *
*     NBAL   : BALOONING STABILITY FLAG                               *
*              NBAL = 0 -------> NO BALOONING TEST                    *
*              NBAL = 1 -------> PERFORM BALOONING TEST               *
*     NBLOPT : BALOONING OPTIMIZATION                                 *
*              NBLOPT = 0 -------> NO BALOONING OPTIMIZATION          *
*              NBLOPT = 1 -------> BALOONING OPTIMIZATION             *
*     NBSFUN : SELECT SHAPE OF BOOTSTRAP CURRENT                      *
*              NBSFUN = 1 -------> POLYNOMIAL                         *
*              NBSFUN = 2 -------> POLYNOMIAL IN 3 SECTIONS           *
*              NBSFUN = 3 -------> PRINCETON DEFINITION               *
*     NBSOPT : BOOTSTRAP CURRENT SPECIFIED                            *
*              NBSOPT = 0 -------> NO SPECIFICATION                   *
*              NBSOPT = 1 -------> BOOTSTRAP CURRENT SPECIFIED        *
*     NBSTRP : SELECT RADIAL DEPENDENCE OF BOOTSTRAP CURRENT          *
*              NBSTRP = 1 -------> INDEPENDENT OF S                   *
*              NBSTRP = 2 -------> DEPENDENT OF S, PROFILE SELECTED   *
*                                  WITH NBSFUN                        *
*     NCSCAL : 1 ----> Q IS SPECIFIED BY QSPEC AT S = CSSPEC          *
*              2 ----> TOTAL TOROIDAL CURRENT IS SPECIFIED BY CURRT   *
*              3 ----> Q IS SPECIFIED BY QSPEC AT RHO = CSSPEC        *
*     NDIFPS : 0 ----> USUAL S-MESH ; 1----> S-MESH WITH CONSTANT     *
*              FLUX AREA BETWEEN SURFACES                             *
*     NDIFT  : 0 ----> USUAL THETA-MESH ; 1----> THETA-MESH WITH      *
*              CONSTANT FLUX AREA IN THETA INTERVAL; 2---> THETA-MESH *
*              WITH CONSTANT ARC-LENGHT INTERVAL AT PLASMA SURFACE    *
*     NEGP   : EXPONENT OF GRAD-PSI OF JACOBIAN DEFINING CHI-MESH     *
*     NEQDSK : 0 ----> NORMAL INPUT AND NO EQDSK.OUT OUTPUT (DEFAULT) *
*              1 .AND. NSURF.EQ.6: READ EQDSK TYPE EQUILIBRIUM FROM   *
*                     EXPEQ FILE AND NORMALIZE IT USING RMAJOR AND B0 *
*              2 .AND. NSURF.EQ.6: READ EQDSK TYPE EQUILIBRIUM FROM   *
*                     EXPEQ FILE AND NORMALIZE IT ACCORDING TO NTMF0  *
*              .NE. 0 ---> WRITE EQDSK.OUT FILE                       *
*              -2: FORCES EQDSK.OUT TO BE NORMALIZED AS NTMF0=0       *
*     NER    : EXPONENT OF R OF JACOBIAN DEFINING CHI-MESH            *
*     NFFTOPT: 0 -> NUMERICAL INTEGRATION OF FOURIER FOR MARS         *
*              1 -> FAST FFT FOR FOURIER NEEDED FOR MARS (DEFAULT)    *
*     NFUNC  : SELECT FUNCTIONAL FORM OF TT', I* OR I_PARRALEL        *
*              NFUNC = 1 -------> POLYNOMIAL                          *
*              NFUNC = 2 -------> POLYNOMIAL IN 3 SECTIONS            *
*              NFUNC = 3 -------> PRINCETON DEFINITION                *
*              NFUNC = 4 -------> EXPERIMENTAL DATA                   *
*     NIDEAL : 0 ---> MARS ; 1 ----> ERATO ; 2 ----> LION             *
*              3 ---> NOVA-W AND PEST ; 4 ----> PENN ;                *
*              5 ---> XTOR ; 6 ---> EQDSK (OR USE NEQDSK.NE.0)        *
*     NIPR   : 1 ----> I-PRIME GIVEN BY 1 PROFILE                     *
*              2 ----> I-PRIME GIVEN BY ADDITION OF 2 PROFILES        *
*              3 ----> I-PRIME GIVEN BY ADDITION OF 3 PROFILES        *
*              4 ----> I-PRIME GIVEN BY ADDITION OF 4 PROFILES        *
*     NOPT   : 0 ----> EQUILIBRIUM IS COMPUTED                        *
*              1 ----> EQUILIBRIUM IS READ FROM FILE NIN              *
*             -1 ----> GUESS       IS READ FROM FILE NIN              *
*             -2 ----> GUESS FROM NIN, USE DIRECTLY LARGE DIMENSION   *
*     NPLOT  : 0 ----> NO PLOTS                                       *
*              1 ----> SAVE QUANTITIES NEEDED BY PLOT PROGRAM PCHEASE *
*                      ON FILE NUPLO                                  *
*     NPP    : 1 ----> P-PRIME GIVEN BY 1 PROFILE                     *
*              2 ----> P-PRIME GIVEN BY ADDITION OF 2 PROFILES        *
*     NPPFUN : SELECT FUNCTIONAL FORM OF P'                           *
*              NPPFUN = 1 -------> POLYNOMIAL                         *
*              NPPFUN = 2 -------> POLYNOMIAL IN 3 SECTIONS           *
*              NPPFUN = 3 -------> PRINCETON DEFINITION               *
*              NPPFUN = 4 -------> EXPERIMENTAL DATA                  *
*     NPRPSI : 0 ----> NO PRINT OF PSI                                *
*              1 ----> PRINT OUT PSI, DPSI/DS, DPSI/DT AND D2PSI/DSDT *
*                      ON (SIGMA, THETA) NODES                        *
*     NPROFZ : 0 ----> PROFILE DEFINITION WITH RESPECT TO S           *
*              1 ----> PROFILE DEFINITION IN CARTESIAN COORDINATE     *
*     NQMIN  : > 0 --> QSPEC AT QMINIMUM                              *
*     NRFP   : 0 ----> TOKAMAK EQUILIBRIUM                            *
*              1 ----> REVERSED FIELD PINCH EQUILIBRIUM               *
*     NRSCAL : 0 ----> NO NORMALIZATION                               *
*              1 ----> NORMALIZATION AGAINST RMAG                     *
*     NSMOOTH: 0 ----> NO BICUBIC SPLINE SMOOTHING                    *
*              1 ----> BICUBIC SPLINE SMOOTHING OF BICUBIC HERMITE    *
*                      SOLUTION                                       *
*     NSTTP  : 1 ----> T-TPRIME IS GIVEN                              *
*              2 ----> I* IS GIVEN                                    *
*              3 ----> I_PARALLEL IS GIVEN                            *
*     NSURF  : 1 ----> SOLOVEV PLASMA SHAPE                           *
*              2 ----> INTOR-LIKE PLASMA SHAPE                        *
*              3 ----> RACETRACK LIKE PLASMA SHAPE (FOR THE FUTUR)    *
*              4 ----> UP-DOWN ASYMETRIC PLASMA SHAPE                 *
*              5 ----> OCTOPOLE-LIKE PLASMA SHAPE                     *
*              6 ----> EXPERIMENTAL PLASMA BOUNDARY                   *
*              7 ----> EXPERIMENTAL PLASMA BOUNDARY GIVEN BY FOURIER  *
*                      COMPONENTS                                     *
*     NSYM   : 0 ----> USE OF ASYMETRIC ERATO VERSION FOR INSTABILITY *
*                      CALCULATION                                    *
*              1 ----> USE OF SYMETRIC ERATO VERSION FOR INSTABILITY  *
*                      CALCULATION                                    *
*     NTCASE : 0 ----> CASE IS SETTED UP WITH NAMELIST                *
*              1 ----> SOLOVEV TEST CASE                              *
*              2 ----> JET TEST CASE                                  *
*              3 ----> NET TEST CASE                                  *
*     NTEST  : 0 ----> NO RESIDU COMPUTATION FOR SOLOVEV CASE         *
*              1 ----> WITH RESIDU COMPUTATION FOR SOLOVEV CASE       *
*     NTMF0  : 0 ----> TSURF = 1                                      *
*              1 ----> T0    = 1                                      *
*     NVEXP  : 0 ----> EQUIDISTANT VACUUM MESH FOR MARS               *
*              1 ----> MARS VACUUM MESH WITH EXPONENTIAL DECAYING MESH*
*                      DENSITY SUCH THAT THE MESH SIZE IS CONTINUOUS  *
*                      ON THE PLASMA BOUNDARY                         *
*                                                                     *
***********************************************************************
C
         COMMON /COMCON/
     I   NANAL,   NBAL,    NBLOPT,  NBSFUN,  NBSOPT,  
     I   NBSTRP,  NCSCAL,  NDIFPS,  NDIFT,   NEGP,    
     I   NEQDSK,  NER,     NFFTOPT, NFUNC,   NIDEAL,  
     I   NIPR,    NOPT,    NPLOT,   NPP,     NPPFUN,  
     I   NPRPSI,  NPROFZ,  NPROPT,  NQMIN,   NRFP,    NRSCAL,  
     I   NSMOOTH, NSTTP,   NSURF,   NSYM,    NTCASE,  
     I   NTEST,   NTMF0,   NVEXP,   NINWALL

