
void ca_derextheating
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_deradinvariant
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derenthalpyfluxy
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derkinengfluxy
  (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
   const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
   amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);
