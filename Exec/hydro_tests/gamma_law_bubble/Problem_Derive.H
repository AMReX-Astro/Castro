
void ca_derpi
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derpioverp0
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derrhopert
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dertpert
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

