
void dertexact
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void derterror
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);
