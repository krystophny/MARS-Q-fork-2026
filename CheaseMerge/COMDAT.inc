C*COMDECK COMDAT
C
***********************************************************************
*                                                                     *
*     C.20   EQUILBRIUM NAMELIST                                      *
*                                                                     *
***********************************************************************
C
         NAMELIST /EQDATA/
     R   AFBS,    AFBS2,   AP,      APLACE,  AP2,     ASPCT,   
     R   AT,      AT2,     AT3,     AT4,     AWIDTH,  BEANS,   
     R   BPLACE,  BSFRAC,  BWIDTH,  B0EXP,   CETA,    CFBAL,   CFNRESS,
     R   CFNRESSO,COMPTYP, CPLACE,  CPRESS,  CPRESSO, CQ0,     CSSPEC,
     R   CURRT,   CWIDTH,  DELTA,   DPLACE,  DWIDTH,  ELONG,   EPLACE, 
     R   ETAEI,   EWIDTH,  EPSLON,  GAMMA,   HEND,    PANGLE,  PREDGE,  
     R   PSISCL,  QSPEC,   QPLACE,  QWIDTH,  QWIDTH0, RBOXLEN, RBOXLFT, 
     R   RC,      RCOIL,   RELAX,   REXT,    RNU,     RPEOP,   RZION,   
     R   R0,      R0EXP,   R0W,     RZ0,     RZ0W,    SCALNE,  
     R   SCEXP,   SGMA,    SOLPDA,  SOLPDB,  SOLPDC,  SOLPDD,  
     R   SOLPDE,  THETA0,  TRIANG,  TRIPLT,  XI,      ZBOXLEN,
     R   RMJ,     WIDTHCOILN,       RCCOILN, ZCCOILN, WCCOILN,
     R   RICOILN, ZICOILN, WICOILN, QSHAVE,  ROTE,    VZ2GP,
     R   ZDEL,    ZFRC,    ZPTS,    PSISHAVE,
     I   MSMAX,   NANAL,   NBAL,    NBLC0,   NBLOPT,  NBPSOUT, 
     I   NBSFUN,  NBSOPT,  NBSTRP,  NCHI,    NCSCAL,  NDIFPS,  
     I   NDIFT,   NEGP,    NEQDSK,  NER,     NFFTOPT, NFUNC,   
     I   NIDEAL,  NINMAP,  NINSCA,  NINWALL,  NIPR,    NISO,
     I   NMESHA,  
     I   NMESHB,  NMESHC,  NMESHD,  NMESHE,  NMGAUS,  NOPT,    
     I   NPLOT,   NPOIDA,  NPOIDB,  NPOIDC,  NPOIDD,  NPOIDE,  
     I   NPOIDQ,  NPP,     NPPFUN,  NPPR,    NPROFZ,  NPROPT,  
     I   NPRPSI,  NPSI,    NQMIN,   NRBOX,   NRFP,    NRSCAL,  NS,      
     I   NSGAUS,  NSMOOTH, NSOUR,   NSTTP,   NSURF,   NSYM,    
     I   NT,      NTCASE,  NTGAUS,  NTEST,   NTMF0,   NTNOVA,  
     I   NTOR,    NTURN,   NV,      NVEXP,   NZBOX,   IDIIICOIL,
     I   NVACUUMRNW,       KMETHOD, NPTS
 

