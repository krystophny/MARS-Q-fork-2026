"""Native prescribed-q current must obey jphi=-R*pprime-FFprime/R."""
from pathlib import Path
import subprocess
import math
from support import NATIVE, NativeTest

def kernel(source):
    text=(source/'chease.f').read_text()
    begin=text.index('                  ZX      =-1./ZH*CID2(I1(J3))')
    end=text.index('\nC                 IF (ZS.LT.',begin)
    return text[begin:end]


def build(source,folder):
    folder.mkdir(parents=True,exist_ok=False)
    code='''      PROGRAM CURRENT_ORACLE
      IMPLICIT NONE
      INTEGER J3,I1(1),K,IOS
      DOUBLE PRECISION CID2(2),D2CID2(2),ZCID2(1),ZFUNC(1)
      DOUBLE PRECISION ZFUNCD(1),ZPPRIM(1),PR(1),PJIPHI(1)
      DOUBLE PRECISION SPSIM,ZS,ZH,ZA,ZB,ZX,CPI
      CPI=ACOS(-1.D0)
      J3=1
      I1=1
      DO
         READ(*,*,IOSTAT=IOS) SPSIM,ZS,ZH,ZA,ZB,PR(1),
     &          ZPPRIM(1),ZFUNC(1),ZFUNCD(1),CID2,D2CID2,ZCID2(1)
         IF (IOS.NE.0) EXIT
KERNELPLACE
         WRITE(*,'(ES26.17)') PJIPHI(1)
      ENDDO
      END
'''.replace('KERNELPLACE',kernel(source))
    file=folder/'current.f';file.write_text(code)
    binary=folder/'current'
    r=subprocess.run(['gfortran','-O0','-fcheck=all','-ffixed-line-length-none',str(file),'-o',str(binary)],capture_output=True,text=True)
    if r.returncode:raise RuntimeError(r.stderr[-2000:])
    return binary



class PrescribedQCurrent(NativeTest):
    def test_analytic_current_for_both_q_signs_and_elongations(self):
        rows=[];inputs=[]
        for kappa in [1.,1.7]:
            for q0 in [1.5,-1.5]:
                for rho in [.15,.4,.8]:
                    for pp in [0.,.13]:
                        # Analytic circular/elliptical psi contours in native R,Z.
                        # CIDQ=int dl/(R|gradpsi|), F=2*pi*q/CIDQ.
                        psi_axis=-.1;minor=.1;h=.0002;x=[rho-h*.37,rho+h*.63]
                        q=q0;qprime=0.
                        c=-psi_axis/(math.pi*kappa*minor**2)
                        d=[c*math.sqrt(1-minor**2*v*v) for v in x]
                        dd=[-c*minor**2/(1-minor**2*v*v)**1.5 for v in x]
                        a=(x[1]-rho)/h;b=(rho-x[0])/h
                        dc=a*d[0]+b*d[1]+((a+1)*(a-1)*h*(x[1]-rho)*dd[0]+(b+1)*(b-1)*h*(rho-x[0])*dd[1])/6
                        R=1+minor*rho*.31
                        vals=[psi_axis,rho,h,a,b,R,pp,q,0.,*d,*dd,dc]
                        inputs.append(' '.join(f'{v:.17e}'for v in vals))
                        # Independent exact F^2=Fedge^2+2*FFprime*psi.
                        ff=2*psi_axis*q*q/(kappa**2*minor**2)
                        expected=-R*pp-ff/R
                        rows.append(dict(shape='circle'if kappa==1 else 'ellipse',kappa=kappa,q=q,rho=rho,
                            pprime=pp,R=R,FFprime_exact=ff,current_exact=expected,
                            flux_native=psi_axis*(1-rho*rho),
                            GS_operator_exact=-2*psi_axis/minor**2*(1/kappa**2+1/R)))
        binary=build(NATIVE,self.work/'current')
        result=subprocess.run([str(binary)], input='\n'.join(inputs)+'\n', capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode,0,result.stderr)
        values=[float(s) for s in result.stdout.split()]
        expected=[row['current_exact'] for row in rows]
        self.assertEqual(len(values),len(expected))
        for value,truth in zip(values,expected):
            self.assertLessEqual(abs(value-truth),1e-8+2e-8*abs(truth))
