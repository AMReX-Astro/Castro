
void ca_dergradpoverp
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dergradpoverp1
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dergradpx
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dergradpy
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dergradrhooverrho
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);
