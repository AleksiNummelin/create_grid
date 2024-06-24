subroutine micom_fill(idm, jdm, inw, jnw, ine, jne, ins, jns, inn, jnn, depth)
   ! Fill ocean points that will not be used by MICOM

   use types

   implicit none

   integer, intent(in) :: idm, jdm
   integer, dimension(idm,jdm), intent(in) :: &
      inw, jnw, ine, jne, ins, jns, inn, jnn
   real (r8), dimension(idm,jdm), intent(inout) :: depth

   integer n, i, j, nzero
   logical changed

   n = 0
   changed = .true.

   do while(changed)
      changed = .false.
      n = n + 1
      do j = 1, jdm
        do i = 1, idm
          if (depth(i,j) > 0._r8) then
            nzero = 0
            if (depth(inw(i,j),jnw(i,j)) <= 0._r8) nzero = nzero + 1
            if (depth(ine(i,j),jne(i,j)) <= 0._r8) nzero = nzero + 1
            if (depth(ins(i,j),jns(i,j)) <= 0._r8) nzero = nzero + 1
            if (depth(inn(i,j),jnn(i,j)) <= 0._r8) nzero = nzero + 1
            if (nzero.ge.3) then
              write (*,'(a,i4,a,i4,a)') ' depth(',i,',',j,') set to zero'
              depth(i,j) = 0._r8
              changed = .true.
            endif
          endif
        enddo
      enddo
   enddo

end subroutine micom_fill
