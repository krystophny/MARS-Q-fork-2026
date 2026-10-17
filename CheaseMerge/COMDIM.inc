C*COMDECK COMDIM 
C  
***********************************************************************
*                                                                     *
*     DEFINE ADJUSTABLE DIMENSION STARING FROM:                       *
*                                                                     *
*        NPS    = NUMBER OF RADIAL INTERVALS                          *
*        NPT    = NUMBER OF POLOIDAL INTERVALS                        *
*        NPCHI  = NUMBER OF CHI ANGLES FOR MAPPING                    *
*        NPPSI  = NUMBER OF PSI INTERVALS FOR MAPPING                 *
*        NPISO  = NUMBER OF PSI INTERVALS IN GENERAL                  *
*        MPSMAX = MAXIMUM NUMBER OF POLOIDAL MODES FOR LINEAR CODE    *
*        NPSMAX = MAXIMUM NUMBER OF TOROIDAL MODES FOR LINEAR CODE    *
*        NPTURN = 2 * MAXIMUM NUMBER OF POLOIDAL PERIODS FOR BALOONING*
*                 INTEGRATION                                         *
*        NPPSBAL= NB OF PSI SURFACES FOR WHICH THE BALLOONING OPT.    *
*                 IS COMPUTED IN PARALLEL (TO REDUCE MEMORY IN ABAL)  *
*        NPPSCUB= NB OF ISO SURFACES USED AT ONE TIME IN CUBRT        *
*        NPV    = NUMBER OF RADIAL VACUUM INTERVALS                   *
*        NPMGS  = # OF GAUSS QUADRATURE POINTS USED FOR THE LINE      *
*                 INTEGRALS COMPUTED IN THE MAPPINGS                  *
*        NPSGS  = # OF GAUSS QUADRATURE POINTS USED FOR THE VARIATIO- *
*                 NALS PROBLEM IN THE SIGMA DIRECTION                 *
*        NPTGS  = # OF GAUSS QUADRATURE POINTS USED FOR THE VARIATIO- *
*                 NALS PROBLEM IN THE THETA DIRECTION                 *
*                                                                     *
*                                                                     *
*     INTRODUCE MEMORY FLAGS, TO REDUCE MEMORY FOR UNUSED MAPPING     *
*     IF A FLAG IS SET TO ZERO, SOME ARRAYS USED ONLY FOR THAT MAPPING*
*     ARE SET TO MINIMUM. IF SET TO ONE, STAYS SAME AS BEFORE         *
*                                                                     *
*        MFLGERL : FLAG FOR ERATO AND LION (NIDEAL=1,2)               *
*                  => NDEQ=NDEQ*MFLGERL+1*(1-MFLGERL)                 *
*                                                                     *
*        MFLGMAR : FLAG FOR MARS (NIDEAL=0)                           *
*                  => MPSMAX = MPSMAX*MFLGMAR + 1*(1-MFLGMAR)         *
*                  => NPSMAX = NPSMAX*MFLGMAR + 1*(1-MFLGMAR)         *
*                                                                     *
*        MFLGNVW : FLAG FOR NOVA-PEST (NIDEAL=3)                      *
*                  => ARRAYS IN OUTNVW ARE REDUCED TO (.,1)           *
*                     (NPPSNVW=NPPSI1*MFLGNVW+1*(1-MFLGNVW))          *
*                                                                     *
*        MFLGPEN : FLAG FOR PENN (NIDEAL=4)                           *
*                  => NPISO IS INCREASED                              *
*            ( NPISO=NPPSI1*(((NPMGS+1)/2+1)*MFLGPEN+1*(1-MFLGPEN)) ) *
*                                                                     *
***********************************************************************
C
         PARAMETER(
     I             NPS   = 320,
     I             NPT   = 640,
     I             NPCHI =  601,
     I             NPPSI =  1026,
     I             MFLGERL = 1,
     I             MFLGMAR = 1,
     I             MFLGNVW = 1,
     I             MFLGPEN = 0)
         PARAMETER(
     I             NPBPS  = 5000,
     I             NWBPS0 = 10,
     I             NPPACK = 300,
     I             NPBLC0 = 16,
     I             NPTURN = 20,
     I             NPPSBAL= 20,                   !      <=  NPPSI + 1
     I             NPPSCUB= 10,                   !      <=  NPPSI + 1
     I             MPSMAX = 80*MFLGMAR + 1*(1-MFLGMAR),
     I             NPSMAX =  1*MFLGMAR + 1*(1-MFLGMAR),
     I             NPV    = 360,
     I             NPMGS  = 4,
     I             NPSGS  = 4,
     I             NPTGS  = 4)
C
C
***********************************************************************
*                                                                     *
* AUXILLIARY PARAMETERS                                               *
*     ERATO ---> NDEQ = 25                                            *
*     LION  ---> NDEQ = 29                                            *
*                                                                     *
***********************************************************************
C
         PARAMETER(
     I             NDEQ   = 25*MFLGERL+1*(1-MFLGERL),
     I             NPBAND = 4 * NPT + 12,
     I             NPCHI1 = NPCHI + 1,
     I             NPGAUS = NPSGS * NPTGS,
     I             NPPSI1 = NPPSI + 1,
     I             NPPSNVW= NPPSI1*MFLGNVW+1*(1-MFLGNVW),
     I             NPISO  = NPPSI1*( ((NPMGS+1)/2+1)*MFLGPEN
     I                               + 1*(1-MFLGPEN) ),
     I             NPV1   = NPV + 1,
     I             NSNT   = (NPS + 1) * NPT,
     I             NSP1   = NPS + 1,
     I             NTP1   = NPT + 1,
     I             NTP2   = NPT + 2,
     I             NP4NST = 4 * NSNT,
     I             N4NT   = 4 * NPT,
     I             N8NT   = 8 * NPT,
     I             N16NT  = 16 * NPT)
C
***********************************************************************
*                                                                     *
* PARAMETERS NEEDED TO INITIATE COMMON VARIABLES TO ZERO INTO SUB-    *
* ROUTINE PRESET                                                      *
*                                                                     *
***********************************************************************
C
        PARAMETER(
     I             NBAL1 = NPISO * (NPBLC0 + 7) + 4,
     I             NBAL2 = 36 * NPISO + NPBLC0 + 
     I                     4 * NPPSBAL * (NPTURN * NPCHI + 1) + 5,
     I             NBLA  = (NPBAND + 5) * NP4NST,
     I             NBND1 = 2*14 + 6*NPBPS + 7*NPBPS*NWBPS0 + 1*5 +
     I                     5 + NWBPS0,
     I             NBND2 = 8 + NWBPS0,
     I             NERA  = (NDEQ + 9) * NPCHI * NPPSI1,
     I             NESH1 = 13,
     I             NESH2 = 12*NPPACK+NPCHI*(2*(NPMGS+1)+1)+2*NPCHI1 +
     I                     2*NPPSI1 + 2*NSP1 + 3*NPISO + NPV1 + NPV +
     I                     4*NTP1 + 19,
     I             NETA1 = (MPSMAX+1)*NPSMAX + (32+22+17)*NPMGS*NTP1 +1,
     I             NETA2 = NPISO * MPSMAX * (29 + 12 * NPSMAX),
     I             NETA3 = 64 * NPISO * MPSMAX,
     I             NETA4 = 34 * NPISO * MPSMAX,
     I             NICON = 33,
     I             NINT1 = NSNT + 20 * NPT,
     I             NINT2 = NPGAUS * (3+NPS+3*NPT+16*NPS*NPT),
     I             NIOD  = 21)
        PARAMETER(
     I             NIS2  = 2 * NPMGS * NTP1 * (16 * NPISO + 1),
     I             NISUR = (2*56 + 9) * NPISO + 1,
     I             NMAP2 = 12 * NTP2 * NPISO + 5 * NPCHI1 + 9 * NPCHI + 
     I                     2 * NPCHI * ((NPMGS + 1) + NPPSI1),
     I             NNUM  = 5 + 26,
     I             NPHY1 = 1,
     I             NPHY2 = 68 + 8*15,
     I             NPLO1 = 1,
     I             NPLO2 = NPPSI1 * (NPCHI + 8),
     I             NSOL1 = 3,
     I             NSOL2 = 10,
     I             NCVAC = 6 * MPSMAX * (NPV + NPV1),
     I             NVEV  = 11 * NSNT + 4,
     I             NEQD1 = NPISO*NPISO + 3,
     I             NEQD2 = 2)
