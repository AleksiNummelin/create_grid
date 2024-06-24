subroutine ncerr(status)
! ------------------------------------------------------------------------------
! If a netCDF error is encountered, display error message and abort the
! program, else return quietly.
! ------------------------------------------------------------------------------

   use netcdf
   
   implicit none

   integer, intent(in) :: status

   if (status /= nf90_noerr) then
      write (*,*) 'netCDF error: ', nf90_strerror(status)
      stop
   endif

end subroutine ncerr
