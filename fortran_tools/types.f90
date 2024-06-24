module types

   implicit none

   integer, parameter :: i1 = selected_int_kind(2)
   integer, parameter :: i2 = selected_int_kind(4)
   integer, parameter :: i4 = selected_int_kind(9)
   integer, parameter :: i8 = selected_int_kind(18)
   integer, parameter :: r4 = selected_real_kind(6,37)
   integer, parameter :: r8 = selected_real_kind(13,307)

end module types
