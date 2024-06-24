program micom_grid
! ------------------------------------------------------------------------------
! Produce a MICOM grid file using extended grid information provided in the file
! lonlatang.bin
! ------------------------------------------------------------------------------

   use types, only: i1, i4, r8
   use mod_xc
   use netcdf

   implicit none

   ! Grid dimensions of extended grid
   integer, parameter :: iedm = 2*idm + 2, jedm = 2*jdm + 2

   real (r8), parameter :: &
     pi       = 3.14159265358979324_r8, &
     rad2deg  = 180._r8/pi, &
     deg2rad  = pi/180._r8, &
     rearth   = 6.37122e6_r8, & ! NCAR CCCM constant
     mindepth = 25._r8, &
     area_rme_limit = 1.e-9_r8
   character (len = 120), parameter :: &
     gridid  = 'tnx0.25'
   character (len = 120), parameter :: &
     gridinfo  = 'Tripolar grid with 0.25 degree resolution along the equator'
   character (len = 120), parameter :: &
     datestr  = 'DOTPaleo'
     

   integer (i4), dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
      insw, jnsw, ins, jns, inse, jnse, ine, jne, inne, jnne, inn, jnn, &
      innw, jnnw, jnw, inw, ip, iu, iv, iq, iuu, ivv, iqq, cplmsk
   integer (i4), dimension(:,:), allocatable :: &
      mask_pop, elementConn
   integer (i4), dimension(:), allocatable :: &
      grid_imask, elementMask
   integer (i4) :: ioerr, reclength, i, j, nreg, ie, je, idm_pop, jdm_pop, &
      grid_size, grid_size_dimid, grid_rank_dimid, grid_corners_dimid, &
      grid_dims_varid, grid_center_lat_varid, grid_center_lon_varid, &
      grid_area_varid, grid_imask_varid, grid_corner_lat_varid, &
      grid_corner_lon_varid, &
      nodeCount, elementCount, maxNodePElement, coordDim, nn, ne, &
      ncid, ncvarid, status, nodeCount_dimid, elementCount_dimid, &
      maxNodePElement_dimid, coordDim_dimid, nodeCoords_varid, &
      elementConn_varid, numElementConn_varid, centerCoords_varid, &
      elementArea_varid, elementMask_varid
   integer (i1), dimension(:), allocatable :: &
      numElementConn
   real (r8), dimension (iedm,jedm) :: &
      lon_e, lat_e, angle_e
   real (r8), dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
      plon, plat, ulat, ulon, vlat, vlon, qlat, qlon, &
      pdx, pdy, udx, udy, vdx, vdy, qdx, qdy, &
      parea, uarea, varea, qarea, angle, pdepth, udepth, vdepth, &
      sarea, maxfac
   real (r8), dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,4) :: &
      pclon, pclat, uclat, uclon, vclat, vclon, qclat, qclon
   real (r8), dimension(:,:), allocatable :: &
      ulat_pop, ulon_pop, htn_pop, hte_pop, hus_pop, huw_pop, angle_pop, &
      tlat_pop, tlon_pop, grid_corner_lat, grid_corner_lon, &
      nodeCoords, centerCoords
   real (r8), dimension(:), allocatable :: &
      grid_center_lat, grid_center_lon, grid_area, elementArea
   real (r8) :: pi2, lonofs, r, r2, area_rme, maxfacp, maxfacm
   logical exists

   ! External functions
   real (r8) :: sph_length, sph_polyarea

   call xcspmd

   ! ---------------------------------------------------------------------------
   ! Read extended grid information
   ! ---------------------------------------------------------------------------

   !inquire (iolength = reclength) lon_e
   !open (unit = 9, file = 'lonlatang.bin', status = 'old', &
   !      form = 'unformatted', access = 'direct', recl = reclength, &
   !      iostat = ioerr)
   !
   !if (ioerr /= 0) then
   !   write (*,*) 'Error opening ''lonlatang.bin'''
   !   stop
   !endif
   !
   !read (9, rec = 1 , iostat = ioerr) lon_e
   !read (9, rec = 2 , iostat = ioerr) lat_e
   !read (9, rec = 3 , iostat = ioerr) angle_e
   !
   !close (unit = 9)
   ! open netcdf4
   status = nf90_open("../grid_out/lonlatang_python.nc", nf90_nowrite,ncid)
   !status = nf90_open("lonlatang_matlab.nc", nf90_nowrite,ncid)
   write(*,*) status
   ! read variables
   status = nf90_inq_varid(ncid, "lat", ncvarid)
   write(*,*) status
   status = nf90_get_var(ncid,ncvarid,lat_e)
   write(*,*) status
   status = nf90_inq_varid(ncid, "lon", ncvarid)
   write(*,*) status
   status = nf90_get_var(ncid,ncvarid,lon_e)
   write(*,*) status
   status = nf90_inq_varid(ncid, "ang", ncvarid)
   write(*,*) status
   status = nf90_get_var(ncid,ncvarid,angle_e)
   write(*,*) status
   ! close netcdf4
   status = nf90_close(ncid)
   !
   write (*,*) lon_e(10,10), lat_e(10,10), angle_e(10,10)
   ! ---------------------------------------------------------------------------
   ! Read depth array if it exists
   ! ---------------------------------------------------------------------------

   !inquire(file = 'depth.uf', exist = exists)
   inquire(file = '../grid_out/depth.nc', exist = exists)
   if (exists) then
     !open (10, file = 'depth.uf', form = 'unformatted')
     status = nf90_open("../grid_out/depth.nc", nf90_nowrite,ncid)
     status = nf90_inq_varid(ncid, "z", ncvarid)
     status = nf90_get_var(ncid,ncvarid,pdepth)
     status = nf90_close(ncid)
     !
     !read (10) i, j
     !if (i /= idm .or. j /= jdm) then
     !   write (*,*) 'Inconsistent dimensions in depth.uf!'
     !   stop
     !endif
     !read (10) pdepth
     !close (10)
     write (*,*) 'Depth array is read'
   else
     pdepth = mindepth
     write (*,*) 'Depth array is not read'
   endif

   ! ---------------------------------------------------------------------------
   ! Define grid positions and angle
   ! ---------------------------------------------------------------------------
   !
   do j = 1, jdm
      je = 2*(j - 1) + 2
      do i = 1, idm
         ie = 2*(i - 1) + 2

         qlon(i,j) = lon_e(ie  ,je  )
         qlat(i,j) = lat_e(ie  ,je  )

         plon(i,j) = lon_e(ie+1,je+1)
         plat(i,j) = lat_e(ie+1,je+1)

         ulon(i,j) = lon_e(ie  ,je+1)
         ulat(i,j) = lat_e(ie  ,je+1)

         vlon(i,j) = lon_e(ie+1,je  )
         vlat(i,j) = lat_e(ie+1,je  )

         qclon(i,j,1) = lon_e(ie-1,je-1)
         qclon(i,j,2) = lon_e(ie+1,je-1)
         qclon(i,j,3) = lon_e(ie+1,je+1)
         qclon(i,j,4) = lon_e(ie-1,je+1)
         qclat(i,j,1) = lat_e(ie-1,je-1)
         qclat(i,j,2) = lat_e(ie+1,je-1)
         qclat(i,j,3) = lat_e(ie+1,je+1)
         qclat(i,j,4) = lat_e(ie-1,je+1)

         pclon(i,j,1) = lon_e(ie  ,je  )
         pclon(i,j,2) = lon_e(ie+2,je  )
         pclon(i,j,3) = lon_e(ie+2,je+2)
         pclon(i,j,4) = lon_e(ie  ,je+2)
         pclat(i,j,1) = lat_e(ie  ,je  )
         pclat(i,j,2) = lat_e(ie+2,je  )
         pclat(i,j,3) = lat_e(ie+2,je+2)
         pclat(i,j,4) = lat_e(ie  ,je+2)

         uclon(i,j,1) = lon_e(ie-1,je  )
         uclon(i,j,2) = lon_e(ie+1,je  )
         uclon(i,j,3) = lon_e(ie+1,je+2)
         uclon(i,j,4) = lon_e(ie-1,je+2)
         uclat(i,j,1) = lat_e(ie-1,je  )
         uclat(i,j,2) = lat_e(ie+1,je  )
         uclat(i,j,3) = lat_e(ie+1,je+2)
         uclat(i,j,4) = lat_e(ie-1,je+2)

         vclon(i,j,1) = lon_e(ie  ,je-1)
         vclon(i,j,2) = lon_e(ie+2,je-1)
         vclon(i,j,3) = lon_e(ie+2,je+1)
         vclon(i,j,4) = lon_e(ie  ,je+1)
         vclat(i,j,1) = lat_e(ie  ,je-1)
         vclat(i,j,2) = lat_e(ie+2,je-1)
         vclat(i,j,3) = lat_e(ie+2,je+1)
         vclat(i,j,4) = lat_e(ie  ,je+1)

         angle(i,j) = angle_e(ie+1,je+1)
      enddo
   enddo

   ! Make sure longitudes of grid points are in the range [-180,180) and that
   ! longitudes of cell corners are numerically close to the corresponding cell
   ! point.
   do j = 1, jdm
      do i = 1, idm

         qlon(i,j) = modulo(qlon(i,j) + 180._r8, 360._r8) - 180._r8
         lonofs = qlon(i,j) - 180._r8
         qclon(i,j,1) = modulo(qclon(i,j,1) - lonofs, 360._r8) + lonofs
         qclon(i,j,2) = modulo(qclon(i,j,2) - lonofs, 360._r8) + lonofs
         qclon(i,j,3) = modulo(qclon(i,j,3) - lonofs, 360._r8) + lonofs
         qclon(i,j,4) = modulo(qclon(i,j,4) - lonofs, 360._r8) + lonofs

         plon(i,j) = modulo(plon(i,j) + 180._r8, 360._r8) - 180._r8
         lonofs = plon(i,j) - 180._r8
         pclon(i,j,1) = modulo(pclon(i,j,1) - lonofs, 360._r8) + lonofs
         pclon(i,j,2) = modulo(pclon(i,j,2) - lonofs, 360._r8) + lonofs
         pclon(i,j,3) = modulo(pclon(i,j,3) - lonofs, 360._r8) + lonofs
         pclon(i,j,4) = modulo(pclon(i,j,4) - lonofs, 360._r8) + lonofs

         ulon(i,j) = modulo(ulon(i,j) + 180._r8, 360._r8) - 180._r8
         lonofs = ulon(i,j) - 180._r8
         uclon(i,j,1) = modulo(uclon(i,j,1) - lonofs, 360._r8) + lonofs
         uclon(i,j,2) = modulo(uclon(i,j,2) - lonofs, 360._r8) + lonofs
         uclon(i,j,3) = modulo(uclon(i,j,3) - lonofs, 360._r8) + lonofs
         uclon(i,j,4) = modulo(uclon(i,j,4) - lonofs, 360._r8) + lonofs

         vlon(i,j) = modulo(vlon(i,j) + 180._r8, 360._r8) - 180._r8
         lonofs = vlon(i,j) - 180._r8
         vclon(i,j,1) = modulo(vclon(i,j,1) - lonofs, 360._r8) + lonofs
         vclon(i,j,2) = modulo(vclon(i,j,2) - lonofs, 360._r8) + lonofs
         vclon(i,j,3) = modulo(vclon(i,j,3) - lonofs, 360._r8) + lonofs
         vclon(i,j,4) = modulo(vclon(i,j,4) - lonofs, 360._r8) + lonofs

     enddo
   enddo

   ! ---------------------------------------------------------------------------
   ! Check target grid type and create index arrays of grid neighbors
   ! ---------------------------------------------------------------------------

   if     (sph_length(pclon(idm/2,  1,1),pclat(idm/2,  1,1), &
                      pclon(idm/2,jdm,4),pclat(idm/2,jdm,4)) < 1.e-6) then
      write (*,*) 'Target domain is cyclic in j-direction!'
      nreg = 4
      do j = 1, jdm
         do i = 1, idm
            insw(i,j) = max(i - 1,   1)
            jnsw(i,j) = mod(j - 2 + jdm, jdm) + 1
            ins (i,j) =     i
            jns (i,j) = mod(j - 2 + jdm, jdm) + 1
            inse(i,j) = min(i + 1, idm)
            jnse(i,j) = mod(j - 2 + jdm, jdm) + 1
            ine (i,j) = min(i + 1, idm)
            jne (i,j) =     j
            inne(i,j) = min(i + 1, idm)
            jnne(i,j) = mod(j, jdm) + 1
            inn (i,j) =     i
            jnn (i,j) = mod(j, jdm) + 1
            innw(i,j) = max(i - 1,   1)
            jnnw(i,j) = mod(j, jdm) + 1
            inw (i,j) = max(i - 1,   1)
            jnw (i,j) =     j
         enddo
      enddo
   elseif (sph_length(pclon(  1,jdm/2,1),pclat(  1,jdm/2,1), &
                      pclon(idm,jdm/2,2),pclat(idm,jdm/2,2)) < 1.e-6) then
      i=idm/4
      if  (sph_length(plon(i      ,jdm  ),plat(i      ,jdm  ), &
                      plon(idm-i+1,jdm-1),plat(idm-i+1,jdm-1)) < 1.e-6) then
         write (*,*) 'Target domain is periodic in i-index with arctic patch!'
         nreg = 2
         do j = 1, jdm - 1
            do i = 1, idm
               insw(i,j) = mod(i - 2 + idm, idm) + 1
               jnsw(i,j) = max(j - 1,   1)
               ins (i,j) =     i
               jns (i,j) = max(j - 1,   1)
               inse(i,j) = mod(i, idm) + 1
               jnse(i,j) = max(j - 1,   1)
               ine (i,j) = mod(i, idm) + 1
               jne (i,j) =     j
               inne(i,j) = mod(i, idm) + 1
               jnne(i,j) =     j + 1
               inn (i,j) =     i
               jnn (i,j) =     j + 1
               innw(i,j) = mod(i - 2 + idm, idm) + 1
               jnnw(i,j) =     j + 1
               inw (i,j) = mod(i - 2 + idm, idm) + 1
               jnw (i,j) =     j
            enddo
         enddo
         j = jdm
         do i = 1, idm
            insw(i,j) = mod(i - 2 + idm, idm) + 1
            jnsw(i,j) =     j - 1
            ins (i,j) =     i
            jns (i,j) =     j - 1
            inse(i,j) = mod(i, idm) + 1
            jnse(i,j) =     j - 1
            ine (i,j) = mod(i, idm) + 1
            jne (i,j) =     j
            inne(i,j) = max(idm - i, 1)
            jnne(i,j) =     j - 2
            inn (i,j) = idm - i + 1
            jnn (i,j) =     j - 2
            innw(i,j) = min(idm - i + 2, idm)
            jnnw(i,j) =     j - 2
            inw (i,j) = mod(i - 2 + idm, idm) + 1
            jnw (i,j) =     j
         enddo
      else
         write (*,*) 'Target domain is cyclic in i-direction!'
         nreg = 1
         do j = 1, jdm
            do i = 1, idm
               insw(i,j) = mod(i - 2 + idm, idm) + 1
               jnsw(i,j) = max(j - 1,   1)
               ins (i,j) =     i
               jns (i,j) = max(j - 1,   1)
               inse(i,j) = mod(i, idm) + 1
               jnse(i,j) = max(j - 1,   1)
               ine (i,j) = mod(i, idm) + 1
               jne (i,j) =     j
               inne(i,j) = mod(i, idm) + 1
               jnne(i,j) = min(j + 1, jdm)
               inn (i,j) =     i
               jnn (i,j) = min(j + 1, jdm)
               innw(i,j) = mod(i - 2 + idm, idm) + 1
               jnnw(i,j) = min(j + 1, jdm)
               inw (i,j) = mod(i - 2 + idm, idm) + 1
               jnw (i,j) =     j
            enddo
         enddo
      endif
   else
      write (*,*) 'Target domain is closed!'
      nreg = 0
      do j = 1, jdm
         do i = 1, idm
            insw(i,j) = max(i - 1,   1)
            jnsw(i,j) = max(j - 1,   1)
            ins (i,j) =     i
            jns (i,j) = max(j - 1,   1)
            inse(i,j) = min(i + 1, idm)
            jnse(i,j) = max(j - 1,   1)
            ine (i,j) = min(i + 1, idm)
            jne (i,j) =     j
            inne(i,j) = min(i + 1, idm)
            jnne(i,j) = min(j + 1, jdm)
            inn (i,j) =     i
            jnn (i,j) = min(j + 1, jdm)
            innw(i,j) = max(i - 1,   1)
            jnnw(i,j) = min(j + 1, jdm)
            inw (i,j) = max(i - 1,   1)
            jnw (i,j) =     j
         enddo
      enddo
   endif

   ! ---------------------------------------------------------------------------
   ! Compute the various grid cell lenghts and areas.
   ! ---------------------------------------------------------------------------

   r = rearth
   r2 = r*r

   do j = 1, jdm
      do i = 1, idm
         qdx(i,j) = sph_length(uclon(i,j,1), uclat(i,j,1), &
                               uclon(i,j,2), uclat(i,j,2))*r
         qdy(i,j) = sph_length(vclon(i,j,1), vclat(i,j,1), &
                               vclon(i,j,4), vclat(i,j,4))*r
         pdx(i,j) = sph_length(vclon(i,j,4), vclat(i,j,4), &
                               vclon(i,j,3), vclat(i,j,3))*r
         pdy(i,j) = sph_length(uclon(i,j,2), uclat(i,j,2), &
                               uclon(i,j,3), uclat(i,j,3))*r
         udx(i,j) = sph_length(qclon(i,j,4), qclat(i,j,4), &
                               qclon(i,j,3), qclat(i,j,3))*r
         udy(i,j) = sph_length(pclon(i,j,1), pclat(i,j,1), &
                               pclon(i,j,4), pclat(i,j,4))*r
         vdx(i,j) = sph_length(pclon(i,j,1), pclat(i,j,1), &
                               pclon(i,j,2), pclat(i,j,2))*r
         vdy(i,j) = sph_length(qclon(i,j,2), qclat(i,j,2), &
                               qclon(i,j,3), qclat(i,j,3))*r
      enddo
   enddo
   !
   write(*,*) pclon(1,1,:), pclat(1,1,:)
   ! Compute grid cell area consistent with the POP/CICE definition

   idm_pop=idm
   if (nreg == 2) then
      jdm_pop=jdm - 1
   else
      jdm_pop=jdm
   endif

   do j = 1, jdm
      do i = 1, idm
         sarea(i,j) = sph_polyarea(pclon(i,j,:), pclat(i,j,:), 4)*r2
      enddo
   enddo
   do
      area_rme = 0._r8
      do j = 1, jdm_pop
         do i = 1, idm
            parea(i,j) = .25_r8*(vdx(inn(i,j),jnn(i,j)) + vdx(i,j)) &
                               *(udy(ine(i,j),jne(i,j)) + udy(i,j))
            !if (i==1) then
            !if (isnan(sarea(i,j))) .or. (isnan(parea(j,i))) then
            !    write (*,*) i, j, sarea(i,j), parea(i,j)
            !endif
            maxfac(i,j) = sqrt(sarea(i,j)/parea(i,j)) !check if nans come from there
            area_rme = area_rme + (1._r8 - maxfac(i,j))**2
         enddo
      enddo
      area_rme = sqrt(area_rme/(idm*jdm)) !check if nans come from here
      write (*,*) 'area_rme',area_rme
      if (area_rme < area_rme_limit) exit
      do j = 1, jdm_pop
         do i = 1, idm
            maxfacp = maxfac(i,j)
            maxfacm = maxfac(inw(i,j),jnw(i,j))
            udy(i,j) = .5_r8*(maxfacm + maxfacp)*udy(i,j)
            maxfacm = maxfac(ins(i,j),jns(i,j))
            vdx(i,j) = .5_r8*(maxfacm + maxfacp)*vdx(i,j)
         enddo
      enddo
      if (nreg == 2) then
         do i = idm/2 + 1, idm
            maxfacp = maxfac(i,jdm - 1)
            maxfacm = maxfac(idm - i + 1,jdm - 1)
            vdx(i,jdm) = .5_r8*(maxfacm + maxfacp)*vdx(i,j)
         enddo
         do i = 1, idm/2
            vdx(i,jdm) = vdx(idm-i+1,jdm)
         enddo
      endif
   enddo

   do j = 1, jdm
      do i = 1, idm
         qarea(i,j) = .25_r8*(udx(i,j) + udx(ins(i,j),jns(i,j))) &
                            *(vdy(i,j) + vdy(inw(i,j),jnw(i,j)))
         parea(i,j) = .25_r8*(vdx(inn(i,j),jnn(i,j)) + vdx(i,j)) &
                            *(udy(ine(i,j),jne(i,j)) + udy(i,j))
         uarea(i,j) = .25_r8*(qdx(inn(i,j),jnn(i,j)) + qdx(i,j)) &
                            *(pdy(i,j) + pdy(inw(i,j),jnw(i,j)))
         varea(i,j) = .25_r8*(pdx(i,j) + pdx(ins(i,j),jns(i,j))) &
                            *(qdy(ine(i,j),jne(i,j)) + qdy(i,j))
      enddo
   enddo

   ! If the grid type has a artic patch, make sure overlapping grid cells are
   ! consistent
   if (nreg == 2) then
      do i = 2, idm/2
         qdx  (i,jdm) = qdx  (idm-i+2,jdm)
         qdy  (i,jdm) = qdy  (idm-i+2,jdm)
         qarea(i,jdm) = qarea(idm-i+2,jdm)
      enddo
      do i = 1, idm
         pdx  (i,jdm) = pdx  (idm-i+1,jdm-1)
         pdy  (i,jdm) = pdy  (idm-i+1,jdm-1)
         parea(i,jdm) = parea(idm-i+1,jdm-1)
      enddo
      do i = 2, idm
         udx  (i,jdm) = udx  (idm-i+2,jdm-1)
         udy  (i,jdm) = udy  (idm-i+2,jdm-1)
         uarea(i,jdm) = uarea(idm-i+2,jdm-1)
      enddo
      do i = 1, idm/2
         vdx  (i,jdm) = vdx  (idm-i+1,jdm)
         vdy  (i,jdm) = vdy  (idm-i+1,jdm)
         varea(i,jdm) = varea(idm-i+1,jdm)
      enddo
   endif

   ! ---------------------------------------------------------------------------
   ! Modify depth array and create ocean mask
   ! ---------------------------------------------------------------------------
   ! fill in points not used by micom
   if (exists) then
       call micom_fill(idm, jdm, inw, jnw, ine, jne, ins, jns, inn, jnn, pdepth)
   endif
   !
   do j = 1, jdm
      do i = 1, idm
         if (pdepth(i,j) > 1._r8) then
            pdepth(i,j) = max(mindepth, pdepth(i,j))
            ip(i,j) = 1
         else
            pdepth(i,j) = 0._r8
            ip(i,j) = 0
         endif
      enddo
   enddo

   do j = 1, jdm
      do i = 1, idm
         iu(i,j) = 0
         iv(i,j) = 0
         iq(i,j) = 0
         if (ip(inw(i,j),jnw(i,j)) > 0 .and. ip(i,j) > 0) then
            iu(i,j) = 1
         endif
         if (ip(ins(i,j),jns(i,j)) > 0 .and. ip(i,j) > 0) then
            iv(i,j) = 1
         endif
         if (min(ip(i,j), ip(inw(i,j),jnw(i,j)), &
                 ip(ins(i,j),jns(i,j)), ip(insw(i,j),jnsw(i,j))) > 0) then
            iq(i,j) = 1
         elseif ((ip(i,j) > 0 .and. ip(insw(i,j),jnsw(i,j)) > 0) .or. &
            (ip(inw(i,j),jnw(i,j)) > 0 .and. ip(ins(i,j),jns(i,j)) > 0)) then
            iq(i,j) = 1
         endif
         udepth(i,j) = min(pdepth(i,j), pdepth(inw(i,j),jnw(i,j)))
         vdepth(i,j) = min(pdepth(i,j), pdepth(ins(i,j),jns(i,j)))
      enddo
   enddo

   do j = 1, jdm
      do i = 1, idm
         if ((ip(i,j) + ip(inw(i,j),jnw(i,j))) >= 1) then
            iuu(i,j) = 1
         else 
            iuu(i,j) = 0
         endif  
         if ((ip(i,j) + ip(ins(i,j),jns(i,j))) >= 1) then
            ivv(i,j) = 1
         else 
            ivv(i,j) = 0
         endif
         if ((iu(i,j) + iv(i,j) + &
              iu(ins(i,j),jns(i,j)) + iv(inw(i,j),jnw(i,j))) >= 1) then
            iqq(i,j) = 1
         else
            iqq(i,j) = 0
         endif
      enddo 
   enddo

   cplmsk = ip

   ! ---------------------------------------------------------------------------
   ! Write MICOM grid file
   ! ---------------------------------------------------------------------------

   call ncfopn('../grid_out/grid.nc', 'w')

   call ncdims('x', idm)
   call ncdims('y', jdm)
   call ncdims('nv', 4)
   call ncdimc('pcomp', ip, 1)
   call ncdimc('qcomp', iqq, 1)
   call ncdimc('ucomp', iuu, 1)
   call ncdimc('vcomp', ivv, 1)

   call ncputi('nreg', nreg, 1)

   call ncwrti('insw', 'x y', insw, ip, 0)
   call ncattr('long_name', 'i-index neighbor in logical south-west direction')

   call ncwrti('jnsw', 'x y', jnsw, ip, 0)
   call ncattr('long_name', 'j-index neighbor in logical south-west direction')

   call ncwrti('ins', 'x y', ins, ip, 0)
   call ncattr('long_name', 'i-index neighbor in logical south direction')

   call ncwrti('jns', 'x y', jns, ip, 0)
   call ncattr('long_name', 'j-index neighbor in logical south direction')

   call ncwrti('inse', 'x y', inse, ip, 0)
   call ncattr('long_name', 'i-index neighbor in logical south-east direction')

   call ncwrti('jnse', 'x y', jnse, ip, 0)
   call ncattr('long_name', 'j-index neighbor in logical south-east direction')

   call ncwrti('ine', 'x y', ine, ip, 0)
   call ncattr('long_name', 'i-index neighbor in logical east direction')

   call ncwrti('jne', 'x y', jne, ip, 0)
   call ncattr('long_name', 'j-index neighbor in logical east direction')

   call ncwrti('inne', 'x y', inne, ip, 0)
   call ncattr('long_name', 'i-index neighbor in logical north-east direction')

   call ncwrti('jnne', 'x y', jnne, ip, 0)
   call ncattr('long_name', 'j-index neighbor in logical north-east direction')

   call ncwrti('inn', 'x y', inn, ip, 0)
   call ncattr('long_name', 'i-index neighbor in logical north direction')

   call ncwrti('jnn', 'x y', jnn, ip, 0)
   call ncattr('long_name', 'j-index neighbor in logical north direction')

   call ncwrti('innw', 'x y', innw, ip, 0)
   call ncattr('long_name', 'i-index neighbor in logical north-west direction')

   call ncwrti('jnnw', 'x y', jnnw, ip, 0)
   call ncattr('long_name', 'j-index neighbor in logical north-west direction')

   call ncwrti('inw', 'x y', inw, ip, 0)
   call ncattr('long_name', 'i-index neighbor in logical west direction')

   call ncwrti('jnw', 'x y', jnw, ip, 0)
   call ncattr('long_name', 'j-index neighbor in logical west direction')

   call ncwrtr('plon', 'x y', plon, ip, 0)
   call ncattr('units', 'degrees_east')
   call ncattr('long_name', 'Longitude at p-points')
   call ncattr('corners', 'pclon')

   call ncwrtr('plat', 'x y', plat, ip, 0)
   call ncattr('units', 'degrees_north')
   call ncattr('long_name', 'Latitude at p-points')
   call ncattr('corners', 'pclat')

   call ncwrtr('ulon', 'x y', ulon, iu, 0)
   call ncattr('units', 'degrees_east')
   call ncattr('long_name', 'Longitude at u-points')
   call ncattr('corners', 'uclon')

   call ncwrtr('ulat', 'x y', ulat, iu, 0)
   call ncattr('units', 'degrees_north')
   call ncattr('long_name', 'Latitude at u-points')
   call ncattr('corners', 'uclat')

   call ncwrtr('vlon', 'x y', vlon, iv, 0)
   call ncattr('units', 'degrees_east')
   call ncattr('long_name', 'Longitude at v-points')
   call ncattr('corners', 'vclon')

   call ncwrtr('vlat', 'x y', vlat, iv, 0)
   call ncattr('units', 'degrees_north')
   call ncattr('long_name', 'Latitude at v-points')
   call ncattr('corners', 'vclat')

   call ncwrtr('qlon', 'x y', qlon, iq, 0)
   call ncattr('units', 'degrees_east')
   call ncattr('long_name', 'Longitude at q-points')
   call ncattr('corners', 'qclon')

   call ncwrtr('qlat', 'x y', qlat, iq, 0)
   call ncattr('units', 'degrees_north')
   call ncattr('long_name', 'Latitude at q-points')
   call ncattr('corners', 'qclat')

   call ncwrtr('pdx', 'x y', pdx, ip, 0)
   call ncattr('units', 'm')
   call ncattr('long_name', 'Grid scale at p-points in x-direction')

   call ncwrtr('pdy', 'x y', pdy, ip, 0)
   call ncattr('units', 'm')
   call ncattr('long_name', 'Grid scale at p-points in y-direction')

   call ncwrtr('udx', 'x y', udx, iu, 0)
   call ncattr('units', 'm')
   call ncattr('long_name', 'Grid scale at u-points in x-direction')

   call ncwrtr('udy', 'x y', udy, iu, 0)
   call ncattr('units', 'm')
   call ncattr('long_name', 'Grid scale at u-points in y-direction')

   call ncwrtr('vdx', 'x y', vdx, iv, 0)
   call ncattr('units', 'm')
   call ncattr('long_name', 'Grid scale at v-points in x-direction')

   call ncwrtr('vdy', 'x y', vdy, iv, 0)
   call ncattr('units', 'm')
   call ncattr('long_name', 'Grid scale at v-points in y-direction')

   call ncwrtr('qdx', 'x y', qdx, iq, 0)
   call ncattr('units', 'm')
   call ncattr('long_name', 'Grid scale at q-points in x-direction')

   call ncwrtr('qdy', 'x y', qdy, iq, 0)
   call ncattr('units', 'm')
   call ncattr('long_name', 'Grid scale at q-points in y-direction')

   call ncwrtr('parea', 'x y', parea, ip, 0)
   call ncattr('units', 'm2')
   call ncattr('long_name', 'Area of p-cells')

   call ncwrtr('uarea', 'x y', uarea, iu, 0)
   call ncattr('units', 'm2')
   call ncattr('long_name', 'Area of u-cells')

   call ncwrtr('varea', 'x y', varea, iv, 0)
   call ncattr('units', 'm2')
   call ncattr('long_name', 'Area of v-cells')

   call ncwrtr('qarea', 'x y', qarea, iq, 0)
   call ncattr('units', 'm2')
   call ncattr('long_name', 'Area of q-cells')

   call ncwrtr('angle', 'x y', angle, ip, 0)
   call ncattr('units', 'radians')
   call ncattr('long_name', &
      'Local angle between x-direction and west-east direction at p-points')

   call ncwrtr('pclon', 'x y nv', pclon, ip, 0)
   call ncattr('units', 'degrees_east')
   call ncattr('long_name', 'Longitude of p-cell corners')

   call ncwrtr('pclat', 'x y nv', pclat, ip, 0)
   call ncattr('units', 'degrees_north')
   call ncattr('long_name', 'Latitude of p-cell corners')

   call ncwrtr('uclon', 'x y nv', uclon, iu, 0)
   call ncattr('units', 'degrees_east')
   call ncattr('long_name', 'Longitude of u-cell corners')

   call ncwrtr('uclat', 'x y nv', uclat, iu, 0)
   call ncattr('units', 'degrees_north')
   call ncattr('long_name', 'Latitude of u-cell corners')

   call ncwrtr('vclon', 'x y nv', vclon, iv, 0)
   call ncattr('units', 'degrees_east')
   call ncattr('long_name', 'Longitude of v-cell corners')

   call ncwrtr('vclat', 'x y nv', vclat, iv, 0)
   call ncattr('units', 'degrees_north')
   call ncattr('long_name', 'Latitude of v-cell corners')

   call ncwrtr('qclon', 'x y nv', qclon, iq, 0)
   call ncattr('units', 'degrees_east')
   call ncattr('long_name', 'Longitude of q-cell corners')

   call ncwrtr('qclat', 'x y nv', qclat, iq, 0)
   call ncattr('units', 'degrees_north')
   call ncattr('long_name', 'Latitude of q-cell corners')

   call ncwrtr('pdepth', 'x y', pdepth, ip, 0)
   call ncattr('units', 'm')
   call ncattr('long_name', 'Depth of p-cells')

   call ncwrtr('udepth', 'x y', udepth, iu, 0)
   call ncattr('units', 'm')
   call ncattr('long_name', 'Depth of u-cells')

   call ncwrtr('vdepth', 'x y', vdepth, iv, 0)
   call ncattr('units', 'm')
   call ncattr('long_name', 'Depth of v-cells')

   call ncwrti('pmask', 'x y', ip, ip, 0)
   call ncattr('long_name', 'Ocean mask of p-cells')

   call ncwrti('umask', 'x y', iu, iu, 0)
   call ncattr('long_name', 'Ocean mask of u-cells')

   call ncwrti('vmask', 'x y', iv, iv, 0)
   call ncattr('long_name', 'Ocean mask of v-cells')

   call ncwrti('qmask', 'x y', iq, iq, 0)
   call ncattr('long_name', 'Ocean mask of q-cells')

   call ncwrti('cplmask', 'x y', cplmsk, cplmsk, 0)
   call ncattr('long_name', 'Ocean mask of CESM')

   call ncfcls

   ! ---------------------------------------------------------------------------
   ! Make POP grid information
   ! ---------------------------------------------------------------------------

   allocate(mask_pop(idm_pop,jdm_pop))
   allocate(ulat_pop(idm_pop,jdm_pop))
   allocate(ulon_pop(idm_pop,jdm_pop))
   allocate(htn_pop(idm_pop,jdm_pop))
   allocate(hte_pop(idm_pop,jdm_pop))
   allocate(hus_pop(idm_pop,jdm_pop))
   allocate(huw_pop(idm_pop,jdm_pop))
   allocate(angle_pop(idm_pop,jdm_pop))
   allocate(tlat_pop(idm_pop,jdm_pop))
   allocate(tlon_pop(idm_pop,jdm_pop))

   do j = 1, jdm_pop
      je = 2*(j - 1) + 2
      do i = 1, idm_pop
         ie = 2*(i - 1) + 2
         mask_pop(i,j) = ip(i,j)
         ulat_pop(i,j) = lat_e(ie+2,je+2)*deg2rad
         ulon_pop(i,j) = lon_e(ie+2,je+2)*deg2rad
         tlat_pop(i,j) = plat(i,j)
         tlon_pop(i,j) = plon(i,j)
         htn_pop(i,j) = vdx(inn(i,j),jnn(i,j))*1.e2_r8
         hte_pop(i,j) = udy(ine(i,j),jne(i,j))*1.e2_r8
         hus_pop(i,j) = udx(ine(i,j),jne(i,j))*1.e2_r8
         huw_pop(i,j) = vdy(inn(i,j),jnn(i,j))*1.e2_r8
         angle_pop(i,j)=angle_e(ie+2,je+2)
      enddo
   enddo
   !
   ! NETCDF VERSION
   !
   write (*,*) 'Write cice grid file'
   write (*,*) 'dims', idm_pop, jdm_pop
   !
   !call ncfopn_cice('grid_cice.nc', 'w')

   !call ncdims_cice('x', idm_pop)
   !call ncdims_cice('y', jdm_pop)
   
   !call ncwrti_cice('mask', 'x y', mask_pop, 0.0*mask_pop, 0)
   !
   !call ncattr('long_name', 'Ocean mask of CICE')
   
   !call ncwrtr_cice('ANGLE', 'x y', angle_pop, mask_pop, 0)
   !call ncattr('long_name', 'grid angle of CICE')
   
   !call ncwrtr_cice('ULAT', 'x y', ulat_pop, mask_pop, 0)
   !call ncattr('long_name', 'Latitude at U/V point of CICE')
   
   !call ncwrtr_cice('ULON', 'x y', ulon_pop, mask_pop, 0)
   !call ncattr('long_name', 'Longitude at U/V point of CICE')
   
   !call ncwrtr_cice('TLAT', 'x y', tlat_pop, mask_pop, 0)
   !call ncattr('long_name', 'Latitude at T point of CICE')
   
   !call ncwrtr_cice('TLON', 'x y', tlon_pop, mask_pop, 0)
   !call ncattr('long_name', 'Longitude at T point of CICE')
   
   !call ncfcls_cice
   !write(*,*) 'done'
   !
   ! BINARY VERSION
   inquire (iolength = reclength) mask_pop
   open (unit = 9, file = '../grid_out/topography_' // trim(datestr) // '.ieeei4', &
         form = 'unformatted', access = 'direct', recl = reclength, &
         iostat = ioerr)
   if (ioerr /= 0) then
      write (*,*) 'Error opening ''topography_' // trim(datestr) // '.ieeei4'''
      stop
   endif
   write (9, rec = 1 , iostat = ioerr) mask_pop
   close (unit = 9)

   inquire (iolength = reclength) ulat_pop
   open (unit = 9, file = '../grid_out/horiz_grid_' // trim(datestr) // '.ieeer8', &
         form = 'unformatted', access = 'direct', recl = reclength, &
         iostat = ioerr)

   if (ioerr /= 0) then
      write (*,*) 'Error opening ''horiz_grid_' // trim(datestr) // '.ieeer8'''
      stop
   endif

   write (9, rec = 1 , iostat = ioerr) ulat_pop
   write (9, rec = 2 , iostat = ioerr) ulon_pop
   write (9, rec = 3 , iostat = ioerr) htn_pop
   write (9, rec = 4 , iostat = ioerr) hte_pop
   write (9, rec = 5 , iostat = ioerr) hus_pop
   write (9, rec = 6 , iostat = ioerr) huw_pop
   write (9, rec = 7 , iostat = ioerr) angle_pop

   close (unit = 9)

   deallocate(mask_pop)
   deallocate(ulat_pop)
   deallocate(ulon_pop)
   deallocate(htn_pop)
   deallocate(hte_pop)
   deallocate(hus_pop)
   deallocate(huw_pop)
   deallocate(angle_pop)

   ! ---------------------------------------------------------------------------
   ! Make SCRIP remapping grid file
   ! ---------------------------------------------------------------------------

   grid_size = idm_pop*jdm_pop

   allocate(grid_center_lat(grid_size), grid_center_lon(grid_size), &
            grid_area(grid_size), grid_imask(grid_size), &
            grid_corner_lat(4, grid_size), grid_corner_lon(4, grid_size))

   nn = 0
   do j = 1, jdm_pop
      do i = 1, idm_pop
         nn = nn + 1
         grid_center_lat(nn) = plat(i,j)
         grid_center_lon(nn) = plon(i,j)
         grid_area(nn) = parea(i,j)/r2
         grid_imask(nn) = ip(i,j)
         grid_corner_lat(:,nn) = pclat(i,j,:)
         grid_corner_lon(:,nn) = pclon(i,j,:)
      enddo
   enddo

   ! Make sure longitudes of grid cell centers are in the range [0,360) and that
   ! longitudes of cell corners are numerically close to the corresponding cell
   ! center.
   do nn = 1, grid_size
      grid_center_lon(nn) = modulo(grid_center_lon(nn), 360._r8)
      lonofs = grid_center_lon(nn) - 180._r8
      grid_corner_lon(1,nn) = &
         modulo(grid_corner_lon(1,nn) - lonofs, 360._r8) + lonofs
      grid_corner_lon(2,nn) = &
         modulo(grid_corner_lon(2,nn) - lonofs, 360._r8) + lonofs
      grid_corner_lon(3,nn) = &
         modulo(grid_corner_lon(3,nn) - lonofs, 360._r8) + lonofs
      grid_corner_lon(4,nn) = &
         modulo(grid_corner_lon(4,nn) - lonofs, 360._r8) + lonofs
   enddo

   ! Create netCDF dataset and enter define mode.
   call ncerr(nf90_create('../grid_out/remap_grid_' // trim(gridid) // '_' // &
                          trim(datestr) // '.nc', nf90_clobber, ncid))

   ! Define dimensions.
   call ncerr(nf90_def_dim(ncid, 'grid_size', grid_size, grid_size_dimid))
   call ncerr(nf90_def_dim(ncid, 'grid_rank', 2, grid_rank_dimid))
   call ncerr(nf90_def_dim(ncid, 'grid_corners', 4, grid_corners_dimid))

   ! Define variables and assign attributes.

   call ncerr(nf90_def_var(ncid, 'grid_dims', nf90_int, grid_rank_dimid, &
                           grid_dims_varid))

   call ncerr(nf90_def_var(ncid, 'grid_center_lat', nf90_double, &
                           grid_size_dimid, grid_center_lat_varid))
   call ncerr(nf90_put_att(ncid, grid_center_lat_varid, "units", "degrees"))

   call ncerr(nf90_def_var(ncid, 'grid_center_lon', nf90_double, &
                           grid_size_dimid, grid_center_lon_varid))
   call ncerr(nf90_put_att(ncid, grid_center_lon_varid, "units", "degrees"))

   call ncerr(nf90_def_var(ncid, 'grid_area', nf90_double, grid_size_dimid, &
                           grid_area_varid))
   call ncerr(nf90_put_att(ncid, grid_area_varid, "units", "radians^2"))

   call ncerr(nf90_def_var(ncid, 'grid_imask', nf90_int, grid_size_dimid, &
                           grid_imask_varid))
   call ncerr(nf90_put_att(ncid, grid_imask_varid, "units", "unitless"))

   call ncerr(nf90_def_var(ncid, 'grid_corner_lat', nf90_double, &
                           (/grid_corners_dimid, grid_size_dimid/), &
                           grid_corner_lat_varid))
   call ncerr(nf90_put_att(ncid, grid_corner_lat_varid, "units", "degrees"))

   call ncerr(nf90_def_var(ncid, 'grid_corner_lon', nf90_double, &
                           (/grid_corners_dimid, grid_size_dimid/), &
                           grid_corner_lon_varid))
   call ncerr(nf90_put_att(ncid, grid_corner_lon_varid, "units", "degrees"))

   ! Global attributes
   call ncerr(nf90_put_att(ncid, nf90_global, 'title', trim(gridinfo)))

   ! End definitions and leave define mode.
   call ncerr(nf90_enddef(ncid))

   ! Provide values for variables.
   call ncerr(nf90_put_var(ncid, grid_dims_varid, (/idm_pop, jdm_pop/)))
   call ncerr(nf90_put_var(ncid, grid_center_lat_varid, grid_center_lat))
   call ncerr(nf90_put_var(ncid, grid_center_lon_varid, grid_center_lon))
   call ncerr(nf90_put_var(ncid, grid_area_varid, grid_area))
   call ncerr(nf90_put_var(ncid, grid_imask_varid, grid_imask))
   call ncerr(nf90_put_var(ncid, grid_corner_lat_varid, grid_corner_lat))
   call ncerr(nf90_put_var(ncid, grid_corner_lon_varid, grid_corner_lon))

   ! Close and save netCDF dataset.
   call ncerr(nf90_close(ncid))

   deallocate(grid_center_lat, grid_center_lon, grid_area, grid_imask, &
              grid_corner_lat, grid_corner_lon)

   ! ---------------------------------------------------------------------------
   ! Make ESMF unstructured grid file
   ! ---------------------------------------------------------------------------

   if (nreg /= 1 .and. nreg /= 2) then
     write (*,*) 'ESMF grid file only made for cyclic i-index'
     stop
   endif

   if (nreg == 2) then
      nodeCount = idm_pop*jdm_pop + idm_pop/2 + 1
   else
      nodeCount = idm_pop*(jdm_pop + 1)
   endif
   elementCount = idm_pop*jdm_pop
   maxNodePElement = 4
   coordDim = 2

   allocate(nodeCoords(coordDim,nodeCount))
   allocate(elementConn(maxNodePElement,elementCount))
   allocate(numElementConn(elementCount))
   allocate(centerCoords(coordDim,elementCount))
   allocate(elementArea(elementCount))
   allocate(elementMask(elementCount))

   nn = 0
   do j = 1, jdm_pop
      je = 2*(j - 1) + 2
      do i = 1, idm_pop
         ie = 2*(i - 1) + 2
         nn = nn + 1
         nodeCoords(1,nn) = lon_e(ie,je)
         nodeCoords(2,nn) = lat_e(ie,je)
      enddo
   enddo
   if (nreg == 2) then
      je = 2*jdm_pop + 2
      do i = 1, idm_pop/2 + 1
         ie = 2*(i - 1) + 2
         nn = nn + 1
         nodeCoords(1,nn) = lon_e(ie,je)
         nodeCoords(2,nn) = lat_e(ie,je)
      enddo
   else
      je = 2*jdm_pop + 2
      do i = 1, idm_pop
         ie = 2*(i - 1) + 2
         nn = nn + 1
         nodeCoords(1,nn) = lon_e(ie,je)
         nodeCoords(2,nn) = lat_e(ie,je)
      enddo
   endif

   ne = 0
   if (nreg == 2) then
      do j = 1, jdm_pop - 1
         do i = 1, idm_pop - 1
            ne = ne + 1
            elementConn(1,ne) = ne
            elementConn(2,ne) = ne + 1
            elementConn(3,ne) = ne + idm_pop + 1
            elementConn(4,ne) = ne + idm_pop
         enddo
         ne = ne + 1
         elementConn(1,ne) = ne
         elementConn(2,ne) = ne - idm_pop + 1
         elementConn(3,ne) = ne + 1
         elementConn(4,ne) = ne + idm_pop
      enddo
      do i = 1, idm_pop/2
         ne = ne + 1
         elementConn(1,ne) = ne
         elementConn(2,ne) = ne + 1
         elementConn(3,ne) = ne + idm_pop + 1
         elementConn(4,ne) = ne + idm_pop
      enddo
      do i = idm_pop/2 + 1, idm_pop - 1
         ne = ne + 1
         elementConn(1,ne) = ne
         elementConn(2,ne) = ne + 1
         elementConn(3,ne) = idm_pop*(jdm_pop + 1) - i + 1
         elementConn(4,ne) = idm_pop*(jdm_pop + 1) - i + 2
      enddo
      ne = ne + 1
      elementConn(1,ne) = ne
      elementConn(2,ne) = ne - idm_pop + 1
      elementConn(3,ne) = idm_pop*jdm_pop + 1
      elementConn(4,ne) = idm_pop*jdm_pop + 2
   else
      do j = 1, jdm_pop
         do i = 1, idm_pop - 1
            ne = ne + 1
            elementConn(1,ne) = ne
            elementConn(2,ne) = ne + 1
            elementConn(3,ne) = ne + idm_pop + 1
            elementConn(4,ne) = ne + idm_pop
         enddo
         ne = ne + 1
         elementConn(1,ne) = ne
         elementConn(2,ne) = ne - idm_pop + 1
         elementConn(3,ne) = ne + 1
         elementConn(4,ne) = ne + idm_pop
      enddo
   endif

   numElementConn(:) = 4

   ne = 0
   do j = 1, jdm_pop
      je = 2*(j - 1) + 2
      do i = 1, idm_pop
         ie = 2*(i - 1) + 2
         ne = ne + 1
         centerCoords(1,ne) = lon_e(ie+1,je+1)
         centerCoords(2,ne) = lat_e(ie+1,je+1)
         elementArea(ne) = parea(i,j)/r2
         elementMask(ne) = ip(i,j)
      enddo
   enddo

   write (*,*) 'nodes   ', nodeCount, nn
   write (*,*) 'elements', elementCount, ne
   write (*,*) 'max node in elements', maxval(elementConn(:,:))

   ! Create netCDF dataset and enter define mode.
   call ncerr(nf90_create('../grid_out/esmf_grid_' // trim(gridid) // '_' // &
                          trim(datestr) // '.nc', nf90_clobber, ncid))

   ! Define dimensions.
   call ncerr(nf90_def_dim(ncid, 'nodeCount', nodeCount, nodeCount_dimid))
   call ncerr(nf90_def_dim(ncid, 'elementCount', elementCount, &
                           elementCount_dimid))
   call ncerr(nf90_def_dim(ncid, 'maxNodePElement', maxNodePElement, &
                           maxNodePElement_dimid))
   call ncerr(nf90_def_dim(ncid, 'coordDim', coordDim, coordDim_dimid))

   ! Define variables and assign attributes.

   call ncerr(nf90_def_var(ncid, 'nodeCoords', nf90_double, &
                           (/coordDim_dimid, nodeCount_dimid/), &
                           nodeCoords_varid))
   call ncerr(nf90_put_att(ncid, nodeCoords_varid, "units", "degrees"))

   call ncerr(nf90_def_var(ncid, 'elementConn', nf90_int, &
                           (/maxNodePElement_dimid, elementCount_dimid/), &
                           elementConn_varid))
   call ncerr(nf90_put_att(ncid, elementConn_varid, "long_name", &
                           "Node indices that define the element connectivity"))
   call ncerr(nf90_put_att(ncid, elementConn_varid, "_FillValue", -1))

   call ncerr(nf90_def_var(ncid, 'numElementConn', nf90_byte, &
                           elementCount_dimid, numElementConn_varid))
   call ncerr(nf90_put_att(ncid, numElementConn_varid, "long_name", &
                           "Number of nodes per element"))

   call ncerr(nf90_def_var(ncid, 'centerCoords', nf90_double, &
                           (/coordDim_dimid, elementCount_dimid/), &
                           centerCoords_varid))
   call ncerr(nf90_put_att(ncid, centerCoords_varid, "units", "degrees"))

   call ncerr(nf90_def_var(ncid, 'elementArea', nf90_double, &
                           elementCount_dimid, elementArea_varid))
   call ncerr(nf90_put_att(ncid, elementArea_varid, "units", "radians^2"))
   call ncerr(nf90_put_att(ncid, elementArea_varid, "long_name", &
                           "Area weights"))

   call ncerr(nf90_def_var(ncid, 'elementMask', nf90_int, elementCount_dimid, &
                           elementMask_varid))
   call ncerr(nf90_put_att(ncid, elementMask_varid, "_FillValue", -9999))

   ! Global attributes
   call ncerr(nf90_put_att(ncid, nf90_global, 'title', trim(gridinfo)))

   ! End definitions and leave define mode.
   call ncerr(nf90_enddef(ncid))

   ! Provide values for variables.
   call ncerr(nf90_put_var(ncid, nodeCoords_varid, nodeCoords))
   call ncerr(nf90_put_var(ncid, elementConn_varid, elementConn))
   call ncerr(nf90_put_var(ncid, numElementConn_varid, numElementConn))
   call ncerr(nf90_put_var(ncid, centerCoords_varid, centerCoords))
   call ncerr(nf90_put_var(ncid, elementArea_varid, elementArea))
   call ncerr(nf90_put_var(ncid, elementMask_varid, elementMask))

   ! Close and save netCDF dataset.
   call ncerr(nf90_close(ncid))

   deallocate(nodeCoords)
   deallocate(elementConn)
   deallocate(numElementConn)
   deallocate(centerCoords)
   deallocate(elementArea)
   deallocate(elementMask)

end program micom_grid
