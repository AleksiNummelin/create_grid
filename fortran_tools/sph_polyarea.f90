function sph_polyarea(lon, lat, n)
   ! This function computes the area in radians squared of a spherical polygon
   ! spanned out by n geographical vertices
   ! Created 20.01.2012 by mats.bentsen@nersc.no

   use types

   implicit none

   integer, intent(in) :: n
   real (r8), intent(in) :: lon(n), lat(n)
   real (r8) :: sph_polyarea

   real (r8) :: a, b, c
   integer :: i, j

   real (r8) :: sph_length

   sph_polyarea = 0._r8
   c = sph_length(lon(1), lat(1), lon(2), lat(2))
   do i = 2, n - 1
      j = i + 1
      a = c
      b = sph_length(lon(i), lat(i), lon(j), lat(j))
      c = sph_length(lon(1), lat(1), lon(j), lat(j))
      sph_polyarea = sph_polyarea &
         + 4._r8*atan(sqrt(max(0._r8,tan(.25_r8*(  a + b + c)) &
                                    *tan(.25_r8*(- a + b + c)) &
                                    *tan(.25_r8*(  a - b + c)) &
                                    *tan(.25_r8*(  a + b - c)))))
   enddo

end function sph_polyarea
