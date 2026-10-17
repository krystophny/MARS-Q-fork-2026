"""The native COMMON-size report must retain totals above signed int32."""
import re
import os
import subprocess

from support import NATIVE, NativeTest, routine


class CheaseCommonSizes(NativeTest):
    def test_total_matches_independent_sum_of_components(self):
        output = routine("OUTPUT")
        declaration = re.search(r"^ +INTEGER\(KIND=SELECTED_INT_KIND\(12\)\) ITOT$",
                                output, re.M).group()
        block = output.split("  120    CONTINUE", 1)[1].split("         RETURN", 1)[0]
        driver = ("      PROGRAM COMMON_SIZES\n" + declaration + "\n"
                  "      INCLUDE 'COMDIM.inc'\n" + block +
                  "      WRITE(6,'(1X,A,I12)') 'COMEQD =',NEQD1+NEQD2\n" +
                  " 2400 FORMAT('COMMON SIZES')\n      END\n")
        (self.work / "native.f").write_text(driver + routine("IVAR"))
        for include in ("COMDIM.inc", "DECLAR.inc"):
            (self.work / include).write_text((NATIVE / include).read_text())
        # This scalar report fixture needs no large static COMMON memory model.
        result = subprocess.run([os.environ.get("FC", "gfortran"),
                                 "-fcheck=all", "-ffixed-line-length-none",
                                 "native.f", "-o", "oracle"], cwd=self.work,
                                capture_output=True, text=True, timeout=300)
        self.assertEqual(result.returncode, 0, result.stderr)
        exe = self.work / "oracle"
        printed = dict(re.findall(r"^\s*(COM\w+|TOTAL)\s*=\s*(-?\d+)\s*$",
                                  self.execute(exe).stdout, re.M))
        total = int(printed.pop("TOTAL"))
        # The legacy COMVAC line prints complex elements; TOTAL counts real words.
        expected = sum(map(int, printed.values())) + int(printed["COMVAC"])
        self.assertEqual(total, expected)
        self.assertGreater(total, 2**31 - 1)
