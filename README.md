# CMAQv5.4-CUDAv3.0

## Getting started
To prepare for compiling CMAQ-CUDA, first compile CMAQ as normal for your chosen chemical mechanism using the ros3 solver.
Ensure you have the NVIDIA HPC SDK installed for .cuf file compilation down the line.

## Compiling CMAQ-CUDA
### What you'll need to have
In addition to the NVHPC SDK, you need the rbsolver.F and rbkernel.cuf files found in this repository. Simply copy them into your build directory (note: rbsolver.F will overwrite the normal rbsolver.F--it may be best to back up that file).

### What you'll need to do
Although there is no need to modify either rbsolver.F or rbkernel.cuf from this repository, you will need to modify the Makefile your build uses to compile. There is an example Makefile (Makefile.intel) included in this repository.

For your compilation, there are three key modifications to make to the Makefile. Refer to Makefile.intel in this repository to complete these steps:
1. Lines 337-343: Editing the GAS list to include rbkernel and omit rbfeval, rbjacob, rbdecomp, and rbsolve as these subroutines are absorbed into rbkernel.
2. Lines 442-472: Editing the compilation steps to compile .cuf files using your nvfortran compiler.
3. Lines 733-744: Adjusting dependencies to reflect rbsolver's new dependence on rbkernel.
