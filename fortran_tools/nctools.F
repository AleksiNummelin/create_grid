c ----------------------------------------------------------------------
c --- NetCDF tools -----------------------------------------------------
c ----------------------------------------------------------------------
c
c Description:
c    nctools.f is a boundle of fortran subroutines designed to write
c    output from the mpi version of MICOM to NetCDF data-format. 
c
c    Prerequisites are the linkage of the f90 netcdf library and the
c    access to the module netcdf.mod.
c
c Comments:
c    The order of actions is relevant: i.e. dimensions should be
c    declared first, then fields not containing the unlimited (time)
c    dimension, and finally the fields containing the time dimension.
c    Attributes may be defined any time.
c
c Contents:
c
c    ncfopn        - creates a file or opens existing file for reading
c    ncfcls        - closes an opened nc-file
c
c    ncdims        - defines a simple axis-dimension
c    ncdimc        - defines a dimension for compressed data storage 
c                    and optionally stores the index information that 
c                    is necessary to de-compress data
c    nctime        - defines a time dimension and stores a time value
c
c    ncattr        - adds a text attribute
c
c    ncputr        - writes a skalar or vector in real format 
c    ncputi        - writes a skalar or vector in integer format  
c    ncgetr        - reads a skalar or vector in real format
c    ncgeti        - reads a skalar or vector in integer format 
c
c    ncread        - reads a 2d or 3d field 
c
c    ncwrti        - writes a 2d or 3d field as int4 
c    ncwrtr        - writes a 2d or 3d field as real8 
c    ncpack        - writes a 2d or 3d field in packed format as int2 
c                    with scale factor and offset
c    nccomp        - writes a 2d or 3d field in compressed (skipping  
c                    land points) in real8 format 
c    nccopa        - writes a 2d or 3d field in compressed/packed format   
c                    as int2 with scale factor and offset  
c
c    ncerro        - displays error massage
c    ncsevl        - evaluates strings with deliminators ' ',':' and '-'
c    ncslen        - corrects intrinsic function len for white spaces
c    ncdnum        - computes a reference number for a given date
c
c Revision history:
c    may2008       - switched from F77-API to F90-API
c    may2008       - removed precision argument from ncwrti and ncwrtr 
c    apr2008       - added header padding
c    apr2008       - added ncputr, ncputi, ncgetr, ncgeti 
c    mar2008       - added nccomp, nccopa, ncread, ncdimc  
c    mar2008       - added MPI support 
c
c Contact:
c    Ingo Bethke (ingo.bethke@nersc.no)
c
c ----------------------------------------------------------------------



      subroutine ncfopn(fnm,faccess)
c
c ----------------------------------------------------------------------
c --- Description:
c       Creates NetCDF file
c
c --- Arguments:
c       char(*) fnm     (in) -  file name
c       char(*) faccess (in) -  'w' for writing and 'r' for reading 
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      character*(*) fnm,faccess
      integer old_mode
c
      if(mnproc.eq.1)then 
c 
c --- - Open file
         if(faccess(1:1).eq.'w')then 
c            status=nf90_create(fnm,nf90_clobber,ncid)
            status=nf90_create(fnm,
     .      cmode=or(nf90_clobber,nf90_64bit_offset),ncid=ncid)
            if(status.ne.nf90_noerr) call ncerro(status,'ncfopn')
         else 
            status=nf90_open(fnm,nf90_nowrite,ncid)
            if(status.ne.nf90_noerr) call ncerro(status,'ncfopn')
         end if
c
c --- - Set variable id to global (needed for global attributes)
         rhid=nf90_global
c
c --- - Initialise header padding
         flgpad=.false. 
c
      end if 
c
      end

      subroutine ncfcls
c
c ----------------------------------------------------------------------
c --- Description:
c       Closes NetCDF file
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      if(mnproc.eq.1)then 
         status=nf90_close(ncid)
         if(status.ne.nf90_noerr) call ncerro(status,'ncfcls')
      end if 
c
      end

      subroutine ncdims(dnm,dim)
c
c ----------------------------------------------------------------------
c --- Description:
c       Creates a simple dimension
c
c --- Arguments:
c       char(*) dnm  (in) -  name of the dimension
c       integer dim  (in) -  number of values on axis
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      integer dim,dimid
      character*(*) dnm 
c
c --- define dimension
      if(mnproc.eq.1)then 
         status=nf90_def_dim(ncid,dnm,dim,dimid)
         if(status.ne.nf90_noerr) call ncerro(status,'ncdims')
      end if 
c
      end


      subroutine ncdimc(dnm,msk,flg)
c
c ----------------------------------------------------------------------
c --- Description:
c       Creates a dimension for compressed data storage and optionally
c       stores the index information that is neccessary to de-compress
c       data.
c
c --- Arguments:
c       char(*)      dnm (in) -  name of dimension, e.g. 'pcomp'
c       integer(...) msk (in) -  2d mask with dimensions
c                                (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) 
c       integer      flg (in) -  1 if index variable should be written 
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      integer n,i,j,dimid,flg 
      integer, dimension(itdm*jtdm) :: ivec    
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: msk   
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rmsk   
      real, dimension(itdm,jtdm) :: rmskt
      character*(*) dnm 
c 
      do j=1,jj
         do i=1,ii
            rmsk(i,j)=msk(i,j)
         end do
      end do 
      call xcaget(rmskt,rmsk,1)
      if(mnproc.eq.1)then 
         n=0
         if(flg.eq.1)then 
            do j=1,jtdm
               do i=1,itdm
                  if(rmskt(i,j).gt.0.5)then
                     n=n+1
                     ivec(n)=i+(j-1)*itdm
                  end if
               end do 
            end do  
         else 
            do j=1,jtdm
               do i=1,itdm
                  if(rmskt(i,j).gt.0.5) n=n+1
               end do 
            end do  
         end if  
         status=nf90_def_dim(ncid,dnm,n,dimid)
         if(status.ne.nf90_noerr) call ncerro(status,'ncdimc')
         if(flg.eq.1)then
            status=nf90_def_var(ncid,dnm,nf90_int,dimid,rhid)
            if(status.ne.nf90_noerr) call ncerro(status,'ncdimc')
            status=nf90_enddef(ncid)
            if(status.ne.nf90_noerr) call ncerro(status,'ncdimc')
            status=nf90_put_var(ncid,rhid,ivec(1:n))
            if(status.ne.nf90_noerr) call ncerro(status,'ncdimc')
            status=nf90_redef(ncid)
            if(status.ne.nf90_noerr) call ncerro(status,'ncdimc')
         end if 
      end if
c
      end



      subroutine ncattr(attname,attrib)
c
c ----------------------------------------------------------------------
c --- Description:
c       Adds a single text attribute to the variable that has been 
c       created last. If no variable has been created then a global 
c       attribute is added.
c
c --- Arguments:
c       char(*) attname (in) -  attribute name
c       char(*) attrib  (in) -  attribute text
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      integer ncslen
      character*(*) attname,attrib
      external function ncslen
c
      if(mnproc.eq.1)then 
         status=nf90_put_att(ncid,rhid,attname,attrib)
         if(status.ne.nf90_noerr) call ncerro(status,'ncattr')
      end if 
c
      end



      subroutine ncputr(vnm,rval,dim)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes a skalar as a global attribute in real8 format. 
c
c --- Arguments:
c       char(*) vnm    (in) -  variable name
c       real rval(dim) (in) -  real skalar or vector
c       integer dim    (in) -  variable dimension 
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      character*(*) vnm
      integer dim
      real rval(dim) 
c
      if(mnproc.eq.1)then
         status=nf90_put_att(ncid,nf90_global,vnm,rval)
         if(status.ne.nf90_noerr) call ncerro(status,'ncputr')
      end if
c
      end



      subroutine ncputi(vnm,ival,dim)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes a skalar as a global attribute in real8 format.
c
c --- Arguments:
c       char(*) vnm       (in) -  variable name
c       integer ival(dim) (in) -  integer skalar or vector
c       integer dim       (in) -  variable dimension 
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      character*(*) vnm
      integer dim,ival(dim)
c
      if(mnproc.eq.1)then
         status=nf90_put_att(ncid,nf90_global,vnm,ival)
         if(status.ne.nf90_noerr) call ncerro(status,'ncputr')
      end if
c
      end



      subroutine ncgetr(vnm,rval)
c
c ----------------------------------------------------------------------
c --- Description:
c      Reads a skalar or vector in real8 format.
c
c --- Arguments:
c       char(*) vnm  (in)  -  variable name
c       real(*) rval (out) -  real skalar or vector
c ----------------------------------------------------------------------
c 
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      character*(*) vnm
      integer dim 
      real rval(*)
c
      if(mnproc.eq.1)then
         status=nf90_inquire_attribute(ncid,nf90_global,vnm,len=dim)
         if(status.ne.nf90_noerr) call ncerro(status,'ncgetr')
         status=nf90_get_att(ncid,nf90_global,vnm,rval(1:dim))
         if(status.ne.nf90_noerr) call ncerro(status,'ncgetr')
      end if
c
      end



      subroutine ncgeti(vnm,ival)
c
c ----------------------------------------------------------------------
c --- Description:
c      Reads a skalar or vector in integer format.
c 
c --- Arguments:
c       char(*)    vnm  (in)  -  variable name
c       integer(*) ival (out) -  integer skalar or vector
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      character*(*) vnm
      integer dim,ival(*) 
c
      if(mnproc.eq.1)then
         status=nf90_inquire_attribute(ncid,nf90_global,vnm,len=dim)
         if(status.ne.nf90_noerr) call ncerro(status,'ncgetr')
         status=nf90_get_att(ncid,nf90_global,vnm,ival(1:dim))
         if(status.ne.nf90_noerr) call ncerro(status,'ncgeti')
      end if
c
      end



      subroutine nctime(datenum,calendar,units,startdate)
c
c ----------------------------------------------------------------------
c --- Description:
c       Create time dimension and store time value.
c       Valid calendars are '360 days', 'noleap',
c       'all_leap', 'julian','gregorian', 'proleptic_gregorian'.
c
c --- Arguments:
c       real    datenum    (in) -  time value
c       char(*) calendar   (in) -  calendar choice
c       char(*) units      (in) -  time units, 
c                                  e.g. 'hours since 0001-01-01'
c       char(*) startdatde (in) -  start date, e.g. '1997-01-01'
c
c --- Comment:
c       The script assumes that the precision of startdate is not
c	higher than the unit specified in the units string. E.g. if
c       units='years since 0001-01-01' and startdate='1997-06-01' then
c       the resulting offset will be rounded to 1996 yrs.
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      integer dimid,rec,ncslen,
     .   bstrn,sstrn,strind(2,10),days,by,bm,bd,bh,sy,sm,sd,sh
      character*(*) calendar,units,startdate
      real datenum,caloffset
      real*8 datenum8
      external function ncslen
      data rec/1/,by/1/,bm/1/,bd/1/,bh/0/,sy/1/,sm/1/,sd/1/,sh/0/
c
      if(mnproc.eq.1)then 
c --- define time dimension
      status=nf90_def_dim(ncid,'time',nf90_unlimited,dimid)
      if(status.ne.nf90_noerr) call ncerro(status,'nctime')
c
c --- define time coordinate variable
      status=nf90_def_var(ncid,'time',nf90_double,dimid,rhid)
      if(status.ne.nf90_noerr) call ncerro(status,'nctime')
      status=nf90_put_att(ncid,rhid,'long_name','time')
      if(status.ne.nf90_noerr) call ncerro(status,'nctime')
      status=nf90_put_att(ncid,rhid,'units',units)
      if(status.ne.nf90_noerr) call ncerro(status,'nctime')
      status=nf90_put_att(ncid,rhid,'calendar',calendar)
      if(status.ne.nf90_noerr) call ncerro(status,'nctime')
c
c --- analyse units string
      call ncsevl(units,bstrn,strind)
      read(units(strind(1,3):strind(2,3)),*) by
      if(bstrn.ge.4) read(units(strind(1,4):strind(2,4)),*) bm
      if(bstrn.ge.5) read(units(strind(1,5):strind(2,5)),*) bd
      if(bstrn.ge.6) read(units(strind(1,6):strind(2,6)),*) bh
c
c --- analyse startdate string
      call ncsevl(startdate,sstrn,strind)
      read(startdate(strind(1,1):strind(2,1)),*) sy
      if(sstrn.ge.2) read(startdate(strind(1,2):strind(2,2)),*) sm
      if(sstrn.ge.3) read(startdate(strind(1,3):strind(2,3)),*) sd
      if(sstrn.ge.4) read(startdate(strind(1,4):strind(2,4)),*) sh
c
c --- calculate calendar offset
      if(units(1:1).eq.'y')then
         caloffset=real(sy-by)
      elseif(units(1:1).eq.'m')then
         caloffset=real(sy-by)*12.+real(sm-bm)
      elseif(units(1:1).eq.'d')then
         call nccaln(days,calendar,by,bm,bd,sy,sm,sd)
         caloffset=real(days)
      elseif(units(1:1).eq.'h')then
         call nccaln(days,calendar,by,bm,bd,sy,sm,sd)
         caloffset=real(days)*24.+real(sh-bh)
      end if
      datenum8=datenum+caloffset
c
c --- Change to data mode and write time
      if(flgpad)then 
         status=nf90_enddef(ncid)
      else
         status=nf90_enddef(ncid,20480,4,0,4)
         flgpad=.true.
      end if 
      if(status.ne.nf90_noerr) call ncerro(status,'nctime')
      status=nf90_put_var(ncid,rhid,datenum8,(/rec/))
      if(status.ne.nf90_noerr) call ncerro(status,'nctime')
c
c --- Switch back to define mode
      status=nf90_redef(ncid)
      if(status.ne.nf90_noerr) call ncerro(status,'nctime')
c
      end if 
c
      end


      
      subroutine ncread(vnm,fld,msk,mskflg,fill) 
c
c ----------------------------------------------------------------------
c --- Description:
c       Reads 2d or 3d field.
c
c --- Arguments:
c       char(*)      vnm    (in)  -  variable name
c       real(...)    fld    (out) -  field array with dimension
c                                    (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,:)
c       integer(...) msk    (in)  -  field mask, if mskflg=1 then the 
c                                    dimensions assumed to be 
c                                    (idm+2*nbdy)*(jdm+2*nbdy)
c       integer      mskflg (in)  -  1 if mask exists  
c       real         fill   (in)  -  fill value to be set on land   
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad
      save /io/
c
      character*(*) vnm 
      character*100 dimname 
      integer i,j,ij,k,kd,n,ndm
      integer, parameter :: maxdm=5, ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
      integer ndims,dimid,dimids(maxdm),dimlen,strn,ncslen,
     .        strind(2,maxdm),msk(*),mskflg
      integer, dimension(maxdm) :: start,count
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk
      real, dimension(itdm,jtdm) :: rfldt,rmskt
      real, dimension(itdm*jtdm) :: rfldtcmp
      real rkd,fill,fld(*),ofs,scf,rcmpflg
      external function ncslen
c
c --- Initialise fields
      rcmpflg=0.
      rkd=1.
      do n=1,maxdm
         start(n)=1
         count(n)=1
      end do  
c       
c --- Inquire dimensions, ect.
      if(mnproc.eq.1)then      
         status=nf90_inq_varid(ncid,vnm,rhid)
         if(status.ne.nf90_noerr) then 
            write(lp,*) 'Problems reading variable ',trim(vnm) 
            call ncerro(status,'ncread')
         end if
         status=nf90_inquire_variable(ncid,rhid,ndims=ndims)
         if(status.ne.nf90_noerr) call ncerro(status,'ncread')
         status=nf90_inquire_variable(ncid,rhid,dimids=dimids)
         if(status.ne.nf90_noerr) call ncerro(status,'ncread')
         ndm=1
         do n=1,ndims
            status=nf90_inquire_dimension(ncid,dimids(n),len=dimlen)
            if(status.ne.nf90_noerr) call ncerro(status,'ncread')
            status=nf90_inquire_dimension(ncid,dimids(n),name=dimname)
            if(status.ne.nf90_noerr) call ncerro(status,'ncread')
            if(dimname(2:5).eq.'comp') rcmpflg=1. 
            if(nint(n+rcmpflg).eq.3) rkd=dimlen
            ndm=ndm*dimlen
         end do
      end if 
      call xcmaxr(rcmpflg)
      call xcmaxr(rkd)
      kd=nint(rkd)
c
c --- Get attributes
      scf=-abs(nf90_fill_double)
      ofs=-abs(nf90_fill_double)
      if(mnproc.eq.1)then      
         status=nf90_get_att(ncid,rhid,'scale_factor',scf)
         if(status.ne.nf90_noerr) scf=1.
         status=nf90_get_att(ncid,rhid,'add_offset',ofs)
         if(status.ne.nf90_noerr) ofs=0.
      end if 
      call xcmaxr(scf)
      call xcmaxr(ofs)
c
c --- get mask for the full domain  
      if(mskflg.eq.1)then
         do j=1,jj
            do i=1,ii
               rmsk(i,j)=msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
            end do
         end do
         call xcaget(rmskt,rmsk,1)
      end if
c
c --- Get data
      do k=1,kd
         if(mnproc.eq.1)then 
            if(rcmpflg.gt.0.5)then 
               start(2)=k
               count(1)=nint(dble(ndm)/dble(kd)) 
            else 
               start(3)=k 
               count(1)=itdm
               count(2)=jtdm
            end if
            if(rcmpflg.gt.0.5)then 
               status=nf90_get_var(ncid,rhid,rfldtcmp,start,count)
               if(status.ne.nf90_noerr) call ncerro(status,'ncread')
               n=0 
               do j=1,jtdm
                  do i=1,itdm
                      if(rmskt(i,j).gt.0.5) then
                         n=n+1
                         rfldt(i,j)=rfldtcmp(n) 
                      end if 
                  end do
               end do  
            else
               status=nf90_get_var(ncid,rhid,rfldt,start,count)
               if(status.ne.nf90_noerr) call ncerro(status,'ncread')
            end if  
         end if
         call xcaput(rfldt,rfld,1)
         if(mskflg.eq.1)then 
            if(scf.eq.1..and.ofs.eq.0)then 
               do j=1,jj
                  do i=1,ii
                     ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                     if(msk(ij).eq.1)then
                        fld(ij+(k-1)*ijdm)=rfld(i,j)
                     else
                        fld(ij+(k-1)*ijdm)=fill
                     end if
                  end do
               end do
            else
               do j=1,jj
                  do i=1,ii
                     ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                     if(msk(ij).eq.1)then
                        fld(ij+(k-1)*ijdm)=rfld(i,j)*scf+ofs
                     else
                        fld(ij+(k-1)*ijdm)=fill
                     end if
                  end do
               end do
            end if  
         else
            if(scf.eq.1..and.ofs.eq.0)then
               do j=1,jj
                  do i=1,ii
                     ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                     fld(ij+(k-1)*ijdm)=rfld(i,j)*scf+ofs 
                  end do
               end do     
            else
               do j=1,jj
                  do i=1,ii
                     ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                     fld(ij+(k-1)*ijdm)=rfld(i,j)
                  end do
               end do
            end if
         end if        
      end do 
c
      end


 
      subroutine ncpack(vnm,dims,fld,msk,sfac,offs,mskflg)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes real field to nc-file in packed format using scale factor
c       and offset.
c
c --- Arguments:
c       char(*)  vnm    (in) -  variable name
c       char(*)  dims   (in) -  axes string, e.g. 'x y z time'
c       real(*)  fld    (in) -  field array with dimension
c       int(*)   msk    (in) -  field mask, if mskflg=1 the dimension
c                               is assumed (idm+2*nbdy)*(jdm+2*nbdy)
c       real     sfac   (in) -  additional scale factor
c       real     offs   (in) -  additional offset
c       integer  mskflg (in) -  flag indicating the presence of a mask
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      character*(*) vnm,dims
      integer i,j,ij,k,kd,n
      integer, parameter :: maxdm=5, ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
      integer, dimension(maxdm) :: start,count
      integer dimid,dimids(maxdm),strn,ncslen,
     .        strind(2,maxdm),msk(*),mskflg     
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk
      real, dimension(itdm,jtdm) :: rfldt,rmskt
      real scf,sfac,ofs,offs,arng(2),fldmin,fldmax,fld(*),rkd,ruvflg
      integer*2, dimension(itdm,jtdm) :: fldout
      integer*2 fill,i2max,i2min,vrng(2)
      parameter(fill=-32768,i2max=32767,i2min=-i2max)
      external function ncslen
      data vrng /i2min, i2max/, start /maxdm*1/, count /itdm,jtdm,1,1,1/
c
c --- Initialise fields
      ruvflg=0.
      rkd=1.
      kd=1
      start(3)=1   
c
c --- define variable 
      if(mnproc.eq.1)then
         call ncsevl(dims,strn,strind)
         do n=1,strn
            status=nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .        dimid)
            if(status.ne.nf90_noerr) call ncerro(status,'ncpack')
            dimids(n)=dimid
            if(n.eq.3)then
               status=nf90_inquire_dimension(ncid,dimid,len=kd)
               if(status.ne.nf90_noerr) call ncerro(status,'ncpack')
            end if 
            if(dims(1:1).eq.'u'.or.dims(1:1).eq.'v') ruvflg=1.
         end do
         rkd=kd
         status=nf90_def_var(ncid,vnm,nf90_int2,dimids(1:strn),rhid)
         if(status.ne.nf90_noerr) call ncerro(status,'ncpack')
      end if
      call xcmaxr(ruvflg)
      call xcmaxr(rkd)
      kd=nint(rkd) 
c
c --- compute scale factor and offset
      fldmin=abs(nf90_fill_double)
      fldmax=-abs(nf90_fill_double)
      if(mskflg.eq.1)then 
c$OMP PARALLEL DO REDUCTION(min:fldmin) REDUCTION(max:fldmax)
c$OMP+  SCHEDULE(STATIC) PRIVATE(ij) 
         do k=1,kd
            do j=1,jj
               do i=1,ii
                  ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  if(msk(ij).eq.1) then
                     fldmin=min(fldmin,(fld(ij+(k-1)*ijdm)*sfac)+offs)
                     fldmax=max(fldmax,(fld(ij+(k-1)*ijdm)*sfac)+offs)
                  end if
               end do
            end do
         end do
c$OMP END PARALLEL DO
      else
c$OMP PARALLEL DO REDUCTION(min:fldmin) REDUCTION(max:fldmax)
c$OMP+  SCHEDULE(STATIC) PRIVATE(ij) 
         do k=1,kd
            do j=1,jj
               do i=1,ii
                  ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
                  fldmin=min(fldmin,(fld(ij+(k-1)*ijdm)*sfac)+offs)
                  fldmax=max(fldmax,(fld(ij+(k-1)*ijdm)*sfac)+offs)
               end do
            end do
         end do
c$OMP END PARALLEL DO
      end if
      call xcminr(fldmin)
      call xcmaxr(fldmax)
      if(ruvflg.gt.0.5)then 
         if(fldmin.ge.fldmax)then
            scf=1.d0
         else
            scf=max(abs(fldmax),abs(fldmin))/dble(i2max)
         end if
         ofs=0.d0
      else
         if(fldmin.ge.fldmax)then
            scf=1.d0
         else
            scf=(fldmax-fldmin)/dble(2*i2max)
         end if
         ofs=0.5*(fldmin+fldmax)
      end if 
      arng(1)=fldmin
      arng(2)=fldmax
c
c --- Define attributes
      if(mnproc.eq.1)then 
         status=nf90_put_att(ncid,rhid,'_FillValue',fill)
         if(status.ne.nf90_noerr) call ncerro(status,'ncpack')
         status=nf90_put_att(ncid,rhid,'valid_range',vrng)
         if(status.ne.nf90_noerr) call ncerro(status,'ncpack')
         status=nf90_put_att(ncid,rhid,'actual_range',arng)
         if(status.ne.nf90_noerr) call ncerro(status,'ncpack')
         status=nf90_put_att(ncid,rhid,'scale_factor',scf)
         if(status.ne.nf90_noerr) call ncerro(status,'ncpack')
         status=nf90_put_att(ncid,rhid,'add_offset',ofs)
         if(status.ne.nf90_noerr) call ncerro(status,'ncpack')
         if(flgpad)then 
            status=nf90_enddef(ncid)
         else
            status=nf90_enddef(ncid,20480,4,0,4)
            flgpad=.true.
         end if 
         if(status.ne.nf90_noerr) call ncerro(status,'ncpack')
      end if
c
c --- get mask for the full domain
      if(mskflg.eq.1)then
         do j=1,jj
            do i=1,ii
               rmsk(i,j)=msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
            end do
         end do
      else
         do j=1,jj
            do i=1,ii
               rmsk(i,j)=1.
            end do
         end do
      end if 
      call xcaget(rmskt,rmsk,1)
c
c --- Prepare output field
      do k=1,kd
         start(3)=k
         do j=1,jj
            do i=1,ii 
               ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
               rfld(i,j)=((fld(ij+(k-1)*ijdm)*sfac)+offs-ofs)/scf
            end do 
         end do     
         call xcaget(rfldt,rfld,1)
         if(mnproc.eq.1)then
            do j=1,jtdm
               do i=1,itdm
                  if(rmskt(i,j).gt..5) then
                     fldout(i,j)=nint(rfldt(i,j))
                  else 
                     fldout(i,j)=fill
                  end if
               end do 
            end do 
            status=nf90_put_var(ncid,rhid,fldout,start,count)
            if(status.ne.nf90_noerr) call ncerro(status,'ncpack')
         end if 
      end do 
c
c --- Put file back to define mode
      if(mnproc.eq.1)then
        status=nf90_redef(ncid)
        if(status.ne.nf90_noerr) call ncerro(status,'ncpack')
        write(lp,*) vnm,arng
      end if
c
      end



      subroutine nccomp(vnm,dims,fld,msk)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes real field to nc-file in compressed format.
c       The horizontal dimensions are replaced with a single, compressed
c       dimension and only ocean points are written.  
c
c --- Arguments:
c       char(*)  vnm  (in) -  variable name
c       char(*)  dims (in) -  axes string, e.g. 'pcomp z time'
c       real(*)  fld  (in) -  field array with dimension
c                             (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kd)
c       int(*)   msk  (in) -  field mask with dimensions 
c                             (idm+2*nbdy)*(jdm+2*nbdy)  
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      character*(*) vnm,dims
      integer omode
      integer i,j,ij,k,kd,n
      integer, parameter :: maxdm=5, ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
      integer, dimension(maxdm) :: start,count
      integer dimid,dimids(maxdm),dimlen(maxdm),strn,ncslen,
     .        strind(2,maxdm),msk(*)     
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk   
      real, dimension(itdm,jtdm) :: rfldt,rmskt
      real, dimension(itdm*jtdm) :: fldout
      real fld(*),rkd
      external function ncslen
      data start /maxdm*1/, count /maxdm*1/ 
c 
c --- Initialise fields
      rkd=1.
      kd=1
      start(2)=1   
      count(1)=1 
c  
c --- define variable 
      if(mnproc.eq.1)then
         call ncsevl(dims,strn,strind)
         do n=1,strn
            status=nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .        dimid)
            if(status.ne.nf90_noerr) call ncerro(status,'nccomp')
            dimids(n)=dimid
            if(n.eq.2)then
               status=nf90_inquire_dimension(ncid,dimid,len=kd)
               if(status.ne.nf90_noerr) call ncerro(status,'nccomp')
            end if 
         end do
         rkd=kd
         status=nf90_def_var(ncid,vnm,nf90_double,dimids(1:strn),rhid)
         if(status.ne.nf90_noerr) call ncerro(status,'nccomp')
         status=nf90_put_att(ncid,rhid,'compress','x y')
         if(status.ne.nf90_noerr) call ncerro(status,'nccomp')
         status=nf90_set_fill(ncid,nf90_nofill,omode)
         if(status.ne.nf90_noerr) call ncerro(status,'nccomp')
         if(flgpad)then 
            status=nf90_enddef(ncid)
         else
            status=nf90_enddef(ncid,20480,4,0,4)
            flgpad=.true.
         end if 
         if(status.ne.nf90_noerr) call ncerro(status,'nccomp')
      end if
      call xcmaxr(rkd)
      kd=nint(rkd) 
c
c --- get mask for the full domain
      do j=1,jj
         do i=1,ii
            rmsk(i,j)=msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
         end do
      end do
      call xcaget(rmskt,rmsk,1)
c
c --- Prepare output field    
      do k=1,kd
         start(2)=k 
         do j=1,jj
            do i=1,ii 
               ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
               rfld(i,j)=fld(ij+(k-1)*ijdm)
            end do 
         end do       
         call xcaget(rfldt,rfld,1) 
         if(mnproc.eq.1)then
            n=0
            do j=1,jtdm
               do i=1,itdm
                  if(rmskt(i,j).gt..5) then
                     n=n+1
                     fldout(n)=rfldt(i,j)
                  end if 
               end do 
            end do 
            count(1)=n
            status=nf90_put_var(ncid,rhid,fldout,start,count)
            if(status.ne.nf90_noerr) call ncerro(status,'nccomp')
         end if 
      end do 
c
c --- Put file back to define mode
      if(mnproc.eq.1)then
        status=nf90_redef(ncid)
        if(status.ne.nf90_noerr) call ncerro(status,'nccomp')
      end if
c
      end



      subroutine nccopa(vnm,dims,fld,msk,sfac,offs)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes real field to nc-file in compressed+packed format using 
c       scale factor and offset.
c       The horizontal dimensions are replaced with a single, compressed
c       dimension and only ocean points are written.  
c
c --- Arguments:
c       char(*)  vnm  (in) -  variable name
c       char(*)  dims (in) -  axes string, e.g. 'pcomp z time'
c       real(*)  fld  (in) -  field array with dimension
c                             (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kd)
c       int(*)   msk  (in) -  field mask with dimensions 
c                             (idm+2*nbdy)*(jdm+2*nbdy)  
c       real     sfac (in) -  additional scale factor
c       real     offs (in) -  additional offset
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      character*(*) vnm,dims
      integer i,j,ij,ijdm,k,kd,maxdm,n
      parameter(maxdm=5)
      integer, dimension(maxdm) :: start,count
      integer dimid,dimids(maxdm),strn,ncslen,
     .        strind(2,maxdm),msk(*)     
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk   
      real, dimension(itdm,jtdm) :: rfldt,rmskt
      real scf,sfac,ofs,offs,arng(2),fldmin,fldmax,fld(*),rkd,ruvflg
      integer*2, dimension(itdm*jtdm) :: fldout
      integer*2 fill,i2max,i2min,vrng(2)
      parameter(fill=-32768,i2max=32767,i2min=-i2max)
      external function ncslen
      data vrng /i2min, i2max/, start /maxdm*1/, count /maxdm*1/
c
c --- Initialise fields
      ruvflg=0.
      rkd=1.
      kd=1
      start(2)=1   
      count(1)=1
c 
c --- define variable 
      if(mnproc.eq.1)then
         call ncsevl(dims,strn,strind)
         do n=1,strn
            status=nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .        dimid)
            if(status.ne.nf90_noerr) call ncerro(status,'nccopa')
            dimids(n)=dimid
            if(n.eq.2)then
               status=nf90_inquire_dimension(ncid,dimid,len=kd)
               if(status.ne.nf90_noerr) call ncerro(status,'nccopa')
            end if 
            if(dims(1:1).eq.'u'.or.dims(1:1).eq.'v') ruvflg=1.
         end do
         rkd=kd
         status=nf90_def_var(ncid,vnm,nf90_int2,dimids(1:strn),rhid)
         if(status.ne.nf90_noerr) call ncerro(status,'nccopa')
      end if
      call xcmaxr(ruvflg)
      call xcmaxr(rkd)
      kd=nint(rkd) 
c
c --- compute scale factor and offset
      fldmin=abs(nf90_fill_double)
      fldmax=-abs(nf90_fill_double)
c$OMP PARALLEL DO REDUCTION(min:fldmin) REDUCTION(max:fldmax)
c$OMP+  SCHEDULE(STATIC)
      ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
      do k=1,kd
         do j=1,jj
            do i=1,ii
               ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1) 
               if(msk(ij).eq.1) then
                  fldmin=min(fldmin,(fld(ij+(k-1)*ijdm)*sfac)+offs)
                  fldmax=max(fldmax,(fld(ij+(k-1)*ijdm)*sfac)+offs)
               end if
            end do
         end do
      end do
c$OMP END PARALLEL DO
      call xcminr(fldmin)
      call xcmaxr(fldmax)
      if(ruvflg.gt.0.5)then 
         if(fldmin.ge.fldmax)then
            scf=1.d0
         else
            scf=max(abs(fldmax),abs(fldmin))/dble(i2max)
         end if
         ofs=0.d0
      else
         if(fldmin.ge.fldmax)then
            scf=1.d0
         else
            scf=(fldmax-fldmin)/dble(2*i2max)
         end if
         ofs=0.5*(fldmin+fldmax)
      end if 
      arng(1)=fldmin
      arng(2)=fldmax
c
c --- Define attributes
      if(mnproc.eq.1)then
         status=nf90_put_att(ncid,rhid,'valid_range',vrng)
         if(status.ne.nf90_noerr) call ncerro(status,'nccopa')
         status=nf90_put_att(ncid,rhid,'actual_range',arng)
         if(status.ne.nf90_noerr) call ncerro(status,'nccopa')
         status=nf90_put_att(ncid,rhid,'scale_factor',scf)
         if(status.ne.nf90_noerr) call ncerro(status,'nccopa')
         status=nf90_put_att(ncid,rhid,'add_offset',ofs)
         if(status.ne.nf90_noerr) call ncerro(status,'nccopa')
         status=nf90_put_att(ncid,rhid,'compress','x y')
         if(status.ne.nf90_noerr) call ncerro(status,'nccopa')
         if(flgpad)then 
            status=nf90_enddef(ncid)
         else
            status=nf90_enddef(ncid,20480,4,0,4)
            flgpad=.true.
         end if 
         if(status.ne.nf90_noerr) call ncerro(status,'nccopa')
      end if
c
c --- get mask for the full domain
      do j=1,jj
         do i=1,ii
            rmsk(i,j)=msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
         end do
      end do
      call xcaget(rmskt,rmsk,1)
c
c --- Prepare output field    
      do k=1,kd
         do j=1,jj
            do i=1,ii 
               ij=i+nbdy+(idm+2*nbdy)*(j+nbdy-1)
               rfld(i,j)=((fld(ij+(k-1)*ijdm)*sfac)+offs-ofs)/scf
            end do 
         end do       
         call xcaget(rfldt,rfld,1) 
         if(mnproc.eq.1)then
            n=0
            do j=1,jtdm
               do i=1,itdm
                  if(rmskt(i,j).gt.0.5) then
                     n=n+1
                     fldout(n)=nint(rfldt(i,j))
                  end if 
               end do 
            end do 
            start(2)=k 
            count(1)=n 
            status=nf90_put_var(ncid,rhid,fldout,start,count)
            if(status.ne.nf90_noerr) call ncerro(status,'nccopa')
         end if 
      end do 
c
c --- Put file back to define mode
      if(mnproc.eq.1)then
        status=nf90_redef(ncid)
        if(status.ne.nf90_noerr) call ncerro(status,'nccopa')
        write(lp,*) vnm,arng
      end if
c
      end 



      subroutine ncwrtr(vnm,dims,fld,msk,mskflg)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes field to nc-file as real8. 
c
c --- Arguments:
c       char(*)  vnm       (in) -  variable name
c       char(*)  dims      (in) -  axes string, e.g. 'pcomp z time'
c       real(*)  fld       (in) -  field array with dimension
c                                  (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kd)
c       int(*)   msk       (in) -  field mask with dimensions 
c                                  (idm+2*nbdy)*(jdm+2*nbdy) 
c       integer  mskflg    (in) -  set to 1 if mask is used, 0 else   
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      character*(*) vnm,dims 
      integer i,j,ij,k,kd,n
      integer, parameter :: maxdm=5,ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
      integer dimid,dimids(maxdm),strn,ncslen,
     .        strind(2,maxdm),msk(*),mskflg     
      integer, dimension(maxdm) :: start,count
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk   
      real, dimension(itdm,jtdm) :: rfldt,rmskt
      real rkd,fill,fld(*)
      external function ncslen
      data start /maxdm*1/, count /itdm,jtdm,1,1,1/
c
c --- Initialise fields
      rkd=1.
      kd=1
      start(3)=1  
c  
c --- define variable 
      fill=nf90_fill_double
      if(mnproc.eq.1)then
         call ncsevl(dims,strn,strind)
         do n=1,strn
            status=nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .        dimid)
            if(status.ne.nf90_noerr) call ncerro(status,'ncwrtr')
            dimids(n)=dimid
            if(n.eq.3)then
               status=nf90_inquire_dimension(ncid,dimid,len=kd)
               if(status.ne.nf90_noerr) call ncerro(status,'ncwrtr')
            end if 
         end do
         rkd=kd
         status=nf90_def_var(ncid,vnm,nf90_double,dimids(1:strn),rhid)
         if(status.ne.nf90_noerr) call ncerro(status,'ncwrtr')
         if(mskflg.eq.1)then
            status=nf90_put_att(ncid,rhid,'_FillValue',fill)
            if(status.ne.nf90_noerr) call ncerro(status,'ncwrtr')
         end if 
         if(flgpad)then 
            status=nf90_enddef(ncid)
         else
            status=nf90_enddef(ncid,20480,4,0,4)
            flgpad=.true.
         end if 
         if(status.ne.nf90_noerr) call ncerro(status,'ncwrtr')
      end if
      call xcmaxr(rkd)
      kd=nint(rkd) 
c
c --- get mask for the full domain
      if(mskflg.eq.1)then 
         do j=1,jj
            do i=1,ii
               rmsk(i,j)=msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
            end do
         end do
      else 
         do j=1,jj
            do i=1,ii
               rmsk(i,j)=1.
            end do
         end do
      end if 
      call xcaget(rmskt,rmsk,1)
c
c --- Prepare output field   
      do k=1,kd
         start(3)=k 
         do j=1,jj
            do i=1,ii
               rfld(i,j)=
     .           fld(i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm)
            end do
         end do  
         call xcaget(rfldt,rfld,1)  
         if(mnproc.eq.1)then
            do j=1,jtdm
               do i=1,itdm
                  if(rmskt(i,j).lt.0.5) rfldt(i,j)=fill
               end do 
            end do 
            status=nf90_put_var(ncid,rhid,rfldt,start,count)
            if(status.ne.nf90_noerr) call ncerro(status,'ncwrtr')
         end if 
      end do 
c
c --- Put file back to define mode
      if(mnproc.eq.1)then
        status=nf90_redef(ncid)
        if(status.ne.nf90_noerr) call ncerro(status,'ncwrtr')
      end if
c
      end


      subroutine ncwrti(vnm,dims,fld,msk,mskflg)
c
c ----------------------------------------------------------------------
c --- Description:
c       Writes field to nc-file as int4. 
c
c --- Arguments:
c       char(*)    vnm       (in) -  variable name
c       char(*)    dims      (in) -  axes string, e.g. 'pcomp z time'
c       integer(*) fld       (in) -  field array with dimension
c                                   (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kd)
c       integer(*) msk       (in) -  field mask with dimensions 
c                                   (idm+2*nbdy)*(jdm+2*nbdy) 
c       integer    mskflg    (in) -  set to 1 if mask is used, 0 else   
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      logical flgpad 
      integer ncid,rhid,status
      common /io/ ncid,rhid,status,flgpad 
      save /io/
c
      character*(*) vnm,dims 
      integer i,j,ij,k,kd,n
      integer, parameter :: maxdm=5,ijdm=(idm+2*nbdy)*(jdm+2*nbdy)
      integer dimid,dimids(maxdm),strn,ncslen,
     .        strind(2,maxdm),msk(*),mskflg,fld(*)     
      integer, dimension(maxdm) :: start,count
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rfld,rmsk   
      real, dimension(itdm,jtdm) :: rfldt,rmskt
      real rkd,fill
      external function ncslen
      data start /maxdm*1/, count /itdm,jtdm,1,1,1/, rkd /1./, kd /1/
c
c --- Initialise fields
      rkd=1.
      kd=1
      start(3)=1   
c  
c --- define variable 
      fill=nf90_fill_int
      if(mnproc.eq.1)then
         call ncsevl(dims,strn,strind)
         do n=1,strn
            status=nf90_inq_dimid(ncid,dims(strind(1,n):strind(2,n)),
     .        dimid)
            if(status.ne.nf90_noerr) call ncerro(status,'ncwrti')
            dimids(n)=dimid
            if(n.eq.3)then
               status=nf90_inquire_dimension(ncid,dimid,len=kd)
               if(status.ne.nf90_noerr) call ncerro(status,'ncwrti')
            end if 
         end do
         rkd=kd
         status=nf90_def_var(ncid,vnm,nf90_int,dimids(1:strn),rhid)
         if(status.ne.nf90_noerr) call ncerro(status,'ncwrti')
         if(mskflg.eq.1)then
            status=nf90_put_att(ncid,rhid,'_FillValue',fill)
         end if
         if(status.ne.nf90_noerr) call ncerro(status,'ncwrti')
         if(flgpad)then 
            status=nf90_enddef(ncid)
         else
            status=nf90_enddef(ncid,20480,4,0,4)
            flgpad=.true.
         end if 
         if(status.ne.nf90_noerr) call ncerro(status,'ncwrti')
      end if
      call xcmaxr(rkd)
      kd=nint(rkd) 
c
c --- get mask for the full domain
      if(mskflg.eq.1)then
         do j=1,jj
            do i=1,ii
               rmsk(i,j)=msk(i+nbdy+(idm+2*nbdy)*(j+nbdy-1))
            end do
         end do
      else
         do j=1,jj
            do i=1,ii
               rmsk(i,j)=1.
            end do
         end do
      end if
      call xcaget(rmskt,rmsk,1)
c
c --- Prepare output field    
      do k=1,kd
         start(3)=k
         do j=1,jj
            do i=1,ii
               rfld(i,j)=
     .           fld(i+nbdy+(idm+2*nbdy)*(j+nbdy-1)+(k-1)*ijdm)
            end do
         end do
         call xcaget(rfldt,rfld,1)
         if(mnproc.eq.1)then
            do j=1,jtdm
               do i=1,itdm
                  if(rmskt(i,j).lt.0.5) rfldt(i,j)=fill
               end do
            end do
            status=nf90_put_var(ncid,rhid,rfldt,start,count)
            if(status.ne.nf90_noerr) call ncerro(status,'ncwrti')
         end if
      end do
c
c --- Put file back to define mode
      if(mnproc.eq.1)then
        status=nf90_redef(ncid)
        if(status.ne.nf90_noerr) call ncerro(status,'ncwrti')
      end if
c
      end


c ----------------------------------------------------------------------
c --- auxilary subroutines ---------------------------------------------
c ----------------------------------------------------------------------

 

      subroutine ncerro(status,strg)
c
c ----------------------------------------------------------------------
c --- Description:
c       Display error massages related to the netcdf library routines
c
c --- Arguments:
c       char(*) strg (in) -  name of subroutine where the error occured
c ----------------------------------------------------------------------
c
      use mod_xc
      use netcdf
c
      implicit none
c
      integer status
      character*(*) strg
c
      write(lp,*) strg,': ',nf90_strerror(status)
      call xchalt('(ncerro)')
             stop '(ncerro)'
c
      end



      subroutine ncsevl(strg,strgn,strgind)
c
c ----------------------------------------------------------------------
c --- Description:
c       Finds the number and the locations of sub-strings. 
c       Valid deliminators are ' ', '-' and ':'.  
c
c --- Arguments:
c       char(*) strg   (in)  -  input string
c       int strgn      (out) -  number of sub-strings
c       int strgind(*) (out) -  start/end locations of sub-strings (the
c                              dimension must at least equal to strgn*2)
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) strg
      character charold,charnew
      integer strgn,strgind(*),i
c
      charold=' '
      strgn=0
      do i=1,len(strg)
         charnew=strg(i:i)
         if((charold.eq.' '.or.charold.eq.'-'.or.charold.eq.':').and.
     .      (charnew.ne.' '.and.charnew.ne.'-'.and.charnew.ne.':'))then
            strgn=strgn+1
            strgind(strgn)=i
         elseif((charnew.eq.' '.or.charnew.eq.'-'.or.charnew.eq.':')
     .    .and.(charold.ne.' '.and.charold.ne.'-'.and.charold.ne.':'))
     .   then
            strgn=strgn+1
            strgind(strgn)=i-1
         endif
         charold=charnew
      enddo
      if(mod(strgn,2).eq.1)then
         strgn=strgn+1
         strgind(strgn)=len(strg)
      endif
      strgn=int(strgn/2.)
c
      end



      integer function ncslen(str)
c
c ----------------------------------------------------------------------
c --- Description:
c       ncslen(str) returns the length minus the number of white spaces
c       of the string str.
c
c --- Arguments:
c       character*(*) str (in) - string to evaluate
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) str
      integer n
c
      ncslen=0
      do n=len(str),1,-1
         if(str(n:n).ne.' '.and.ncslen.eq.0) ncslen=n
      end do
c
      end



      subroutine nccaln(days,calendar,y1,m1,d1,y2,m2,d2)
c
c ----------------------------------------------------------------------
c --- Description:
c       nccaln calculates the amount of days between two
c       dates. Valid calendars are '360 days', 'noleap', 'all_leap',
c       'julian', 'proleptic_gregorian' and 'standard' (udunits mixed
c       julian/gregorian).
c
c --- Arguments:
c       int days         (out) -  number of days
c       char(*) calendar (in)  -  choice of calendar
c       int y1           (in)  -  start year
c       int m1           (in)  -  start month
c       int d1           (in)  -  start day
c       int y2           (in)  -  last year
c       int m2           (in)  -  last month
c       int d2           (in)  -  last day
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) calendar
      integer n1,y1,m1,d1,n2,y2,m2,d2,days,nd1,nd2,nd3
c
      if(calendar(1:3).eq.'sta'.or.calendar(1:3).eq.'mix'.or.
     .   calendar(1:3).eq.'gre')then
         if(y1.lt.1582.or.(y1.eq.1582.and.m1.lt.10).or.
     .      (y1.eq.1582.and.m1.eq.10.and.d1.le.15))then
            call ncdnum('julian',y1,m1,d1,n1)
         else
            call ncdnum('julian',1582,10,15,nd1)
            call ncdnum('proleptic',1582,10,15,nd2)
            call ncdnum('proleptic',y1,m1,d1,nd3)
            n1=nd1-nd2+nd3-10
         endif
         if(y2.lt.1582.or.(y2.eq.1582.and.m2.lt.10).or.
     .      (y2.eq.1582.and.m2.eq.10.and.d2.le.15))then
            call ncdnum('julian',y2,m2,d2,n2)
         else
            call ncdnum('julian',1582,10,15,nd1)
            call ncdnum('proleptic',1582,10,15,nd2)
            call ncdnum('proleptic',y2,m2,d2,nd3)
            n2=nd1-nd2+nd3-10
         endif
      else
         call ncdnum(calendar,y1,m1,d1,n1)
         call ncdnum(calendar,y2,m2,d2,n2)
      endif
      days=n2-n1
c
      end



      subroutine ncdnum(calendar,year,month,day,ndays)
c
c ----------------------------------------------------------------------
c --- Description:
c       Auxilary subroutine for nccaln which calculates the number
c       of days between a given date and a reference date. Thereby the
c       reference date varies from calendar to calendar and is not
c       standerized. The purpose is more of differential nature rather
c       than to create an absolute time axis.
c ----------------------------------------------------------------------
c
      implicit none
c
      character*(*) calendar
      integer y,m,year,month,day,ndays,yr1,yr4,yr100,yr400
      integer acc_leap(12),acc_noleap(12),yoffset
      parameter(yoffset=1000000,yr1=365,yr4=1461,yr100=36524,
     .   yr400=146097)
      data acc_leap   /0,31,60,91,121,152,182,213,244,274,305,335/
      data acc_noleap /0,31,59,90,120,151,181,212,243,273,304,334/
c
      y=year+yoffset
      if(calendar(1:3).eq.'pro')then
         ndays=int((y-1)/400.)*yr400
     .      + int(mod((y-1),400)/100.)*yr100
     .      + int(mod(mod((y-1),400),100)/4.)*yr4
     .      + mod(mod(mod((y-1),400),100),4)*yr1
         if(month.gt.2.and.((mod(y,4).eq.0.and.mod(y,100).ne.0).or.
     .      mod(y,400).eq.0))then
            ndays=ndays+acc_leap(month)+day
         else
            ndays=ndays+acc_noleap(month)+day
         end if
      elseif(calendar(1:3).eq.'jul')then
         ndays=int((y-1)/4.)*yr4+mod(y,4)*yr1
         if(month.gt.2.and.mod(y,4).eq.0)then
            ndays=ndays+acc_leap(month)+day
         else
            ndays=ndays+acc_noleap(month)+day
         end if
      else if(calendar(1:3).eq.'360')then
         ndays=(y-1)*360+30*(month-1)+day
        else if(calendar(1:3).eq.'nol')then
         ndays=(y-1)*365+acc_noleap(month)+day
        else if(calendar(1:3).eq.'all')then
         ndays=(y-1)*366+acc_leap(month)+day
      end if
c
      end



c ----------------------------------------------------------------------
c ----------------------------------------------------------------------
c ----------------------------------------------------------------------
