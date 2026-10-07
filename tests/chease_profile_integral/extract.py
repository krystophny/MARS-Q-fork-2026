"""Extract executed ISOFUN branches and native SPLINE/TRIDAG, without rewriting them."""
from pathlib import Path
import hashlib,json,sys
src=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True);text=src.read_text();isofun=text.split('         SUBROUTINE ISOFUN(KN)',1)[1].split('C*DECK',1)[0]
body=isofun.split('         IF (NTMF0.EQ.0) THEN',1)[1].split('C           WRITE(*,1001)',1)[0]
header='''         SUBROUTINE NATIVE_INTEGRATE(KN,NTMF0,NSTTP,NRFP,
     &    AT,PSIISO,TTP,D2TTP,TMF,CIPR,CIDQ)
         IMPLICIT NONE
         INTEGER KN,NTMF0,NSTTP,NRFP,J7,J7P1,J8,J8P1
         DOUBLE PRECISION AT(10),PSIISO(KN),TTP(KN),D2TTP(KN),
     &    TMF(KN),CIPR(KN),CIDQ(KN),CPI
         CPI=ACOS(-1.D0)
         IF (NTMF0.EQ.0) THEN'''
code=header+body+'         RETURN\n         END\n'
for name in ('SPLINE(X,Y,N,YP2,A,B)','TRIDAG(A,B,R,N,EPS)'):
 unit=text.split('         SUBROUTINE '+name,1)[1].split('C*DECK',1)[0];code+='         SUBROUTINE '+name+unit
code+='''         SUBROUTINE NATIVE_SPLINE(N,X,Y,M)
         INCLUDE 'DECLAR.inc'
         DOUBLE PRECISION X(N),Y(N),M(N),A(N),B(N)
         RC1M14=1.D-14
         CALL SPLINE(X,Y,N,M,A,B)
         END
'''
(out/'native_extracted.f').write_text(code)
(out/'extraction.json').write_text(json.dumps({'source':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'extracted_sha256':hashlib.sha256(code.encode()).hexdigest(),'scope':'Literal ISOFUN NTMF0 branches; SPLINE/TRIDAG native routines; wrapper supplies independent polynomial jets. No PDE/production executable invoked.'},indent=2)+'\n')
