function sph_length(lon1, lat1, lon2, lat2)
   ! Computes the geodesic length in radians between geographical positions
   ! (lon1,lat1) and (lon2,lat2)

   use types

   implicit none

   real (r8), intent(in) :: lon1, lat1, lon2, lat2
   real (r8) :: sph_length

   real (r8), parameter :: rad = 1.74532925199432958e-02_r8

   real (r8) :: lambda, phi, x1, y1, z1, x2, y2, z2

   phi = lon1*rad
   lambda = lat1*rad
   x1 = cos(lambda)*cos(phi) 
   y1 = cos(lambda)*sin(phi) 
   z1 = sin(lambda)

   phi = lon2*rad
   lambda = lat2*rad
   x2 = cos(lambda)*cos(phi) 
   y2 = cos(lambda)*sin(phi) 
   z2 = sin(lambda)

   sph_length = acos(min(1._r8, max(- 1._r8, x1*x2 + y1*y2 + z1*z2)))

end function sph_length
