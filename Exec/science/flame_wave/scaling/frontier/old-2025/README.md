Frontier scaling on 2025-04-09

This uses ROCm 6.3.1

A sample submission script is included here.

These runs disabled the EXTRACXXFLAGS that disabled inlining with HIP,
since those are no longer necessary.  But note that has little effect
on the runtime.

Note that for 2048
