
void ca_derinertialmomentumx
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derinertialmomentumy
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derinertialmomentumz
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derinertialangmomx
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derinertialangmomy
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derinertialangmomz
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derinertialradmomx
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derinertialradmomy
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derinertialradmomz
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derphieff
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derphieffpm_p
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derphieffpm_s
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derrhophiGrav
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derrhophiRot
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_derprimarymask
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);

void ca_dersecondarymask
(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
 const amrex::FArrayBox& datfab, const amrex::Geometry& geomdata,
 amrex::Real /*time*/, const int* /*bcrec*/, int /*level*/);
