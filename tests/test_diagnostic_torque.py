"""Compile actual native writers and production torque against a real-field oracle.

This regression targets the corrected-production staging revision: density work
uses 2*pi**2 and production -n Im(DWK)/(4*pi**2). It must not be transplanted to
an older production normalization without auditing that consuming revision.
"""
from pathlib import Path
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import unittest


def routine(text, name):
    start=text.index('      SUBROUTINE '+name+'(')
    end=text.index('      END SUBROUTINE '+name,start)
    return text[start:end]+'      END SUBROUTINE '+name+'\n'


def run(source):
    source=Path(source);text=(source/'kinetic.f').read_text()
    body=routine(text,'CALCDWKCOMP')
    start=body.find('      TNORM =')
    if start < 0:
        start=body.index('      TORQUENTVI = 0.')
    # Include the complete native aggregation, including per-species output
    # and the closing KNTV conditional; the next section is not torque code.
    end=re.search(r'^C\s+(?:SAME-DRIVE WORK|OUTPUT PROFILES)\b',body[start:],re.M)
    if end is None:raise RuntimeError('Native torque section boundary missing')
    production=body[start:start+end.start()]
    # Execute both exact native radial-fold assignments, independently of each
    # diagnostic writer. The rest of CALCDWKCOMP is unnecessary for constants.
    folds=[]
    for target in ('DWPPARY','DWPPERY'):
        start=body.index('         '+target+'(IS,:) =')
        lines=body[start:].splitlines()
        folds.append('\n'.join(lines[:2]))
    native=routine(text,'WRITEDWKBREAKDOWN')+routine(text,'WRITEDWKDRIVELEDGER')
    driver=r'''
      PROGRAM FIXTURE
      USE DIMENSIM
      USE GLOBALM
      USE RCOMDM
      USE KINETICM
      USE TORQUEM
      USE TOOLBOX
      IMPLICIT NONE
      INTEGER IS,I,INDX,U,KP,IH,NSP,FID,FIDS
      REAL*8 PI2,TH(2),TD(2),TW(2),TNORM,TSPEC
      CHARACTER*16 ARG
      COMPLEX*16 DWPPARX(3,1),DWPPERX(3,1),DWPPARY(3,1),DWPPERY(3,1)
      COMPLEX*16 PYD(3,1,5),EYD(3,1,5),DWHERM(3,1)
      NR=2
      NRP1=3
      NSPECIES=2
      ISPECIES_F0=0
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      READ(ARG,*) RNTOR
      ALLOCATE(CSM(2),CSH(2),INDXDWKC(2,10),INDXHWKC(2,5))
      ALLOCATE(TORQUENTVI(2),TORQUENTVE(2),TORQUENTV(2))
      CSM=(/0.25D0,0.75D0/)
      CSH=(/0.2D0,0.8D0/)
      INDXDWKC=-1
      INDXDWKC(1,5)=1
      INDXDWKC(1,10)=1
      INDXHWKC=-1
      INDXHWKC(1,5)=1
      IFOWP=0
      IFOWT=0
      DWPPARX=(1.0D0,2.0D0)
      DWPPERX=(3.0D0,4.0D0)
      DWPPARY=(5.0D0,6.0D0)
      DWPPERY=(7.0D0,8.0D0)
      CALL WRITEDWKBREAKDOWN(DWPPARX,DWPPERX,DWPPARY,DWPPERY,1)
      PI2=2.0D0*PI*PI
      DO IS=1,NR
FOLD
      ENDDO
C     MANUFACTURED IDEAL-INDUCTION CASE: SAME-DRIVE AND WEAK WORK AGREE.
C     THE REAL-FIELD FORCE QUADRATURE BELOW FIXES THE TORQUE INDEPENDENTLY.
      DWHERM=DWPPARY+DWPPERY
      IF (.TRUE.) THEN
PRODUCTION
      PYD=0
      EYD=0
      PYD(:,:,2)=DWPPARY
      EYD(:,:,2)=DWPPERY
      CALL WRITEDWKDRIVELEDGER(PYD,EYD,1,0.0D0,1.0D0)
      OPEN(NEWUNIT=U,FILE='PRODUCTION.OUT')
      DO IS=1,NR
         WRITE(U,'(5(1X,ES26.17E3))') DWPPARY(IS,1),
     &      DWPPERY(IS,1),TORQUENTVI(IS)
      ENDDO
      CLOSE(U)
      IF (ANY(DWPPARX.NE.(1.0D0,2.0D0))) STOP 2
      IF (ANY(DWPPERX.NE.(3.0D0,4.0D0))) STOP 3
      IF (ANY(PYD(:,:,2).NE.DWPPARY)) STOP 4
      IF (ANY(EYD(:,:,2).NE.DWPPERY)) STOP 5
      END
'''.replace('FOLD','\n'.join(folds)).replace('PRODUCTION\n',production+'\n')
    toolbox='''module toolbox
contains
integer function ASSIGNFREEFILEUNIT()
integer unit
logical opened
do unit=20,99
inquire(unit=unit,opened=opened)
if (.not.opened) then
ASSIGNFREEFILEUNIT=unit
return
endif
enddo
error stop 'No free fixture file unit'
end function
end module
'''
    def rows(path):
        return [[float(x) for x in line.split()] for line in path.read_text().splitlines()
                if line.strip() and not line.lstrip().startswith('%')]
    result=[]
    with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR','/var/tmp'),prefix='mars-torque-') as tmp:
        d=Path(tmp);(d/'native.f').write_text(native);(d/'driver.f').write_text(driver);(d/'toolbox.f90').write_text(toolbox)
        modules=['dimensim.f','globalm.f','rcomdm.f','kineticm.f','torquem.f']
        cmd=['gfortran','-fcheck=all','-fdefault-real-8','-fdefault-double-8','-ffixed-line-length-none',
             *[str(source/m) for m in modules],'toolbox.f90','native.f','driver.f','-o','check.x']
        p=subprocess.run(cmd,cwd=d,text=True,capture_output=True)
        if p.returncode:raise RuntimeError(p.stdout+p.stderr)
        for n in (-3,3):
            subprocess.run([str(d/'check.x'),str(n)],cwd=d,check=True,capture_output=True,text=True)
            # H=cos(n*phi), f=20*sin(n*phi); -<f*dH/dphi> by real quadrature.
            expected=math.fsum(20*n*math.sin(n*2*math.pi*j/4096)**2 for j in range(4096))/4096
            prod=rows(d/'PRODUCTION.OUT');breakdown=rows(d/'DWK_BREAKDOWN.OUT');ledger=rows(d/'DWK_DRIVE_LEDGER.OUT')
            err_prod=max(abs(row[-1]/expected-1) for row in prod)
            err_break=max(abs(row[-1]/expected-1) for row in breakdown)
            err_drive=max(abs(sum(row[-1] for row in ledger if int(row[0])==i)/expected-1) for i in (1,2))
            # Independently known complex density work: 2*pi²*(6+8i,10+12i).
            work=[2*math.pi**2*x for x in (6,8,10,12)]
            work_error=max(abs(row[4+k]/work[k]-1) for row in ledger if int(row[2])==2 for k in range(4))
            work_error=max(work_error,max(abs(row[12+k]/work[k]-1) for row in breakdown for k in range(4)))
            result.append(dict(n=n,real_field_oracle=expected,production_relative_error=err_prod,
                               breakdown_relative_error=err_break,drive_ledger_relative_error=err_drive,
                               unchanged_work_relative_error=work_error))
    return dict(passed=all(max(r[k] for k in ('production_relative_error','breakdown_relative_error',
                    'drive_ledger_relative_error','unchanged_work_relative_error'))<1e-7 for r in result),rows=result)


@unittest.skipUnless(shutil.which('gfortran'), 'native fixture requires gfortran')
class DiagnosticTorqueReports(unittest.TestCase):
    def test_native_reports_match_independent_real_field_force(self):
        here=Path(__file__).resolve()
        roots=[here.parents[1]/'MarsQ_2FK',here.parents[1]/'source/MarsQ_2FK',here.parent/'source/MarsQ_2FK']
        source=next((p for p in roots if (p/'kinetic.f').exists()),None)
        if source is None:raise RuntimeError('Owning native source required')
        result=run(source)
        self.assertTrue(result['passed'],json.dumps(result,indent=2))


if __name__=='__main__':
    import sys
    if len(sys.argv)==1:
        unittest.main()
        raise SystemExit(0)
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source',type=Path);p.add_argument('--output',type=Path)
    a=p.parse_args();result=run(a.source)
    text=json.dumps(result,indent=2)+'\n'
    if a.output:a.output.write_text(text)
    print(text,end='')
    raise SystemExit(0 if result['passed'] else 1)
