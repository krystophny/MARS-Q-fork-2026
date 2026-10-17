MarsQ: source code for MARS-F/K/Q, run "make" under Linux to compile, update makefile for proper compiler
CheaseMerge: source code for CHEASE, solve fixed boundary Grad-Shafranov equation to provide input to MARS-*
EXAMPLE: various test examples
RZplot: post-processing Matlab scripts
Manual: MARS-F/K/Q user manual

Portable builds write executables under build/ and leave the historical
tracked executable untouched:

    make -C MarsQ_2FK gnu
    make -C MarsQ_2FK gnu-debug
    make -C MarsQ_2FK ifx
    make -C MarsQ_2FK nvhpc

Configure the selected compiler's MPI wrapper on PATH first. Each MARS target
cleans compiler intermediates, then writes a binary and its .provenance.json
record with compiler flags, Git source state and artifact SHA-256. Serialize
builds that share a source checkout. Verify an artifact without executing it:

    python3 tools/build_with_provenance.py --verify-manifest /path/to/marsq-gnu.x.provenance.json

The default keeps the dependency-free PAMS backend. Enable the certified
global steady-response backend with a GNU profile and the fortsparse prefix:

    make -C MarsQ_2FK gnu PAMS_CSC=1 FORTSPARSE_PREFIX=/path/to/fortsparse

FORTSPARSE_LIBS supplies any additional platform link flags. The provenance
record also binds the fortsparse modules, library and solver helper.
CheaseMerge supports gnu and ifx targets, retaining the legacy source's
zero-initialization and large-memory-model compiler flags.





