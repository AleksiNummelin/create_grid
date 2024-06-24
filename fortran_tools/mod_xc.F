      module mod_xc
c
c --- dummy mod_xc
c
      implicit none
c
      public
c
#include "dimensions.h"
c
      integer mnproc
c
      contains
c
      subroutine xcspmd
c
      mnproc=1
      ii=idm
      jj=jdm
      lp=6
c
      end subroutine xcspmd
c
      subroutine xcaget(aa,a,mnflg)
c
      real,    intent(out)   :: aa(itdm,jtdm)
      real,    intent(in)    :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    :: mnflg
c
      aa=a
c
      end subroutine xcaget
c
      subroutine xcaput(aa,a,mnflg)
c
      real,    intent(inout) :: aa(itdm,jtdm)
      real,    intent(out)   :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    :: mnflg
c
      a=aa
c
      end subroutine xcaput
c
      subroutine xcmaxr(a)
c
      real, intent(in) :: a
c
      end subroutine xcmaxr
c
      subroutine xcminr(a)
c
      real, intent(in) :: a
c
      end subroutine xcminr
c
      subroutine xchalt(cerror)
c
      character*(*), intent(in) :: cerror
c
      if     (cerror.ne.' ') then
        write(lp,*) '**************************************************'
        write(lp,*) cerror
      endif
      write(lp,*) '**************************************************'
      write(lp,*) 'XCHALT CALLED ON PROC = ',mnproc
      write(lp,*) '**************************************************'
      call flush(lp)
c
      stop '(xchalt)'
c
      end subroutine xchalt
c
      end module mod_xc
