!***********************************************************************
      module m_sfcdrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/07/03
!     Modification: 2002/07/26, 2002/09/02, 2002/09/09, 2003/04/30,
!                   2003/05/19, 2003/07/15, 2003/12/12, 2004/04/10,
!                   2004/04/15, 2004/05/07, 2004/05/31, 2004/06/10,
!                   2004/08/01, 2004/09/01, 2005/01/14, 2005/02/10,
!                   2005/04/04, 2006/01/10, 2006/09/21, 2006/11/06,
!                   2007/01/20, 2007/08/24, 2008/05/02, 2008/07/01,
!                   2008/08/19, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2011/11/10, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures to interpolate the surface data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chkmxn
      use m_chkstd
      use m_comindx
      use m_comkind
      use m_commpi
      use m_cpondsfc
      use m_currpe
      use m_destroy
      use m_estimsfc
      use m_getcname
      use m_getdate
      use m_getiname
      use m_getrij
      use m_getrname
      use m_hint2d
      use m_hintlnd
      use m_inichar
      use m_outsfc
      use m_outsst
      use m_outstd05
      use m_outstd12
      use m_rdice
      use m_rdland
      use m_rdsst
      use m_setproj
      use m_sststep
      use m_undefice
      use m_undefsst

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: sfcdrv, s_sfcdrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface sfcdrv

        module procedure s_sfcdrv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_sfcdrv(fpsfcdat,fpidate,                             &
     &                    fpintopt_lnd,fplnduse_lnd,fpsstitv,           &
     &                    fpalbe_lnd,fpbeta_lnd,fpz0m_lnd,fpz0h_lnd,    &
     &                    fpcap_lnd,fpnuu_lnd,ni,nj,ri,rj,land,albe,    &
     &                    beta,z0m,z0h,cap,nuu,sst,kai,di,dj,tmp1,      &
     &                    nid_lnd,njd_lnd,landat,ltmp1,nid_sst,njd_sst, &
     &                    sstdat,stmp1,nid_ice,njd_ice,icedat,itmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsfcdat
                       ! Formal parameter of unique index of sfcdat

      integer, intent(in) :: fpidate
                       ! Formal parameter of unique index of idate

      integer, intent(in) :: fpintopt_lnd
                       ! Formal parameter of unique index of intopt_lnd

      integer, intent(in) :: fplnduse_lnd
                       ! Formal parameter of unique index of lnduse_lnd

      integer, intent(in) :: fpsstitv
                       ! Formal parameter of unique index of sstitv

      integer, intent(in) :: fpalbe_lnd
                       ! Formal parameter of unique index of albe_lnd

      integer, intent(in) :: fpbeta_lnd
                       ! Formal parameter of unique index of beta_lnd

      integer, intent(in) :: fpz0m_lnd
                       ! Formal parameter of unique index of z0m_lnd

      integer, intent(in) :: fpz0h_lnd
                       ! Formal parameter of unique index of z0h_lnd

      integer, intent(in) :: fpcap_lnd
                       ! Formal parameter of unique index of cap_lnd

      integer, intent(in) :: fpnuu_lnd
                       ! Formal parameter of unique index of nuu_lnd

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nid_lnd
                       ! Land use data dimension in x direction

      integer, intent(in) :: njd_lnd
                       ! Land use data dimension in y direction

      integer, intent(in) :: nid_sst
                       ! Sea surface temperature data dimension
                       ! in x direction

      integer, intent(in) :: njd_sst
                       ! Sea surface temperature data dimension
                       ! in y direction

      integer, intent(in) :: nid_ice
                       ! Sea ice distribution data dimension
                       ! in x direction

      integer, intent(in) :: njd_ice
                       ! Sea ice distribution data dimension
                       ! in y direction

! Internal shared variables

      character(len=108) sfcdat
                       ! Control flag of input surface data type

      character(len=108) idate
                       ! Forecast start date
                       ! with Gregorian calendar, yyyymmddhhmm

      character(len=12) cdate
                       ! Current forecast date
                       ! with Gregorian calendar, yyyymmddhhmm

      integer intopt_lnd
                       ! Option for
                       ! linear interpolating for land use data

      integer iid      ! Index of land use categories data

      integer(kind=i8) it
                       ! Index of main do loop

      integer(kind=i8) nstp0
                       ! Start index of main do loop

      integer(kind=i8) nstp1
                       ! End index of main do loop

      integer(kind=i8) ctime
                       ! Model current forecast time

      integer(kind=i8) sst103
                       ! 1000 x int(sstitv + 0.1)

      integer stat     ! Runtime status

      integer lnduse_lnd(1:100)
                       ! User specified land use data

      integer, intent(inout) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      integer, intent(inout) :: landat(1:nid_lnd,1:njd_lnd)
                       ! Land use in data

      real sstitv      ! Time interval of
                       ! sea surface temperature data file

      real x0          ! x origin of model grid
      real y0          ! y origin of model grid

      real cpj(1:7)    ! Map projection parameters of model grid

      real x0lnd       ! x origin of land use data grid
      real y0lnd       ! y origin of land use data grid

      real cpjlnd(1:7) ! Map projection parameters of land use data grid

      real x0sst       ! x origin of sea surface temperature data grid
      real y0sst       ! y origin of sea surface temperature data grid

      real cpjsst(1:7) ! Map projection parameters
                       ! of sea surface temperature data grid

      real x0ice       ! x origin of sea ice distribution data grid
      real y0ice       ! y origin of sea ice distribution data grid

      real cpjice(1:7) ! Map projection parameters
                       ! of sea ice distribution data grid

      real albe_lnd(1:100)
                       ! User specified albedo data

      real beta_lnd(1:100)
                       ! User specified
                       ! evapotranspiration efficiency data

      real z0m_lnd(1:100)
                       ! User specified
                       ! roughness length data for velocity

      real z0h_lnd(1:100)
                       ! User specified
                       ! roughness length data for scalar

      real cap_lnd(1:100)
                       ! User specified thermal capacity data

      real nuu_lnd(1:100)
                       ! User specified thermal diffusivity data

      real, intent(inout) :: ri(0:ni+1,0:nj+1)
                       ! Real indices in data region in x direction

      real, intent(inout) :: rj(0:ni+1,0:nj+1)
                       ! Real indices in data region in y direction

      real, intent(inout) :: albe(0:ni+1,0:nj+1)
                       ! Albedo

      real, intent(inout) :: beta(0:ni+1,0:nj+1)
                       ! Evapotranspiration efficiency

      real, intent(inout) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(inout) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(inout) :: cap(0:ni+1,0:nj+1)
                       ! Thermal capacity

      real, intent(inout) :: nuu(0:ni+1,0:nj+1)
                       ! Thermal diffusivity

      real, intent(inout) :: sst(0:ni+1,0:nj+1)
                       ! Sea surface temperature

      real, intent(inout) :: kai(0:ni+1,0:nj+1)
                       ! Sea ice distribution

      real, intent(inout) :: di(0:ni+1,0:nj+1)
                       ! Distance between data and model grid points
                       ! in x direction

      real, intent(inout) :: dj(0:ni+1,0:nj+1)
                       ! Distance between data and model grid points
                       ! in y direction

      real, intent(inout) :: sstdat(1:nid_sst,1:njd_sst)
                       ! Sea surface temperature in data

      real, intent(inout) :: icedat(1:nid_ice,1:njd_ice)
                       ! Sea ice distribution in data

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: ltmp1(1:nid_lnd,1:njd_lnd)
                       ! Temporary array

      real, intent(inout) :: stmp1(1:nid_sst,1:njd_sst)
                       ! Temporary array

      real, intent(inout) :: itmp1(1:nid_ice,1:njd_ice)
                       ! Temporary array

! Remark

!     kai: This variable is also temporary.

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(sfcdat)
      call inichar(idate)

! -----

! Get the required namelist variables.

      call getcname(fpsfcdat,sfcdat)
      call getcname(fpidate,idate)
      call getiname(fpintopt_lnd,intopt_lnd)
      call getrname(fpsstitv,sstitv)

      do iid=0,99

        call getiname(fplnduse_lnd+iid,lnduse_lnd(iid+1))
        call getrname(fpalbe_lnd+iid,albe_lnd(iid+1))
        call getrname(fpbeta_lnd+iid,beta_lnd(iid+1))
        call getrname(fpz0m_lnd+iid,z0m_lnd(iid+1))
        call getrname(fpz0h_lnd+iid,z0h_lnd(iid+1))
        call getrname(fpcap_lnd+iid,cap_lnd(iid+1))
        call getrname(fpnuu_lnd+iid,nuu_lnd(iid+1))

      end do

! -----

! Set the map projection parameters of the model grid.

      call setproj(idmpopt,idnspol,iddx,iddy,idulat,idulon,idriu,idrju, &
     &             idtlat1,idtlat2,idtlon,'solver  ',6,x0,y0,cpj)

! -----

!!!!! Create the model input interpolated soil and ice surface file.

      if(sfcdat(1:1).eq.'o'.or.sfcdat(3:3).eq.'o') then

! Set the map projection parameters of the land use data grid and read
! out the data from the land use data file.

        if(sfcdat(1:1).eq.'o') then

          stat=0

          call setproj(idmpopt_lnd,idnspol_lnd,iddx_lnd,iddy_lnd,       &
     &                 idulat_lnd,idulon_lnd,idriu_lnd,idrju_lnd,       &
     &                 idtlat1_lnd,idtlat2_lnd,idtlon_lnd,              &
     &                 'surface ',7,x0lnd,y0lnd,cpjlnd)

          call rdland(iddatdir,idncdat,idwlngth,stat,nid_lnd,njd_lnd,   &
     &                landat)

          call chkerr(stat)

          if(stat.lt.0) then

            call destroy('rdland  ',6,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          call chkstd(root)

        end if

! -----

! Set the map projection parameters of the sea ice distribution data
! grid and read out the data from the sea ice distribution data file.

        if(sfcdat(3:3).eq.'o') then

          stat=0

          call setproj(idmpopt_ice,idnspol_ice,iddx_ice,iddy_ice,       &
     &                 idulat_ice,idulon_ice,idriu_ice,idrju_ice,       &
     &                 idtlat1_ice,idtlat2_ice,idtlon_ice,              &
     &                 'surface ',7,x0ice,y0ice,cpjice)

          call rdice(iddatdir,idncdat,idwlngth,stat,nid_ice,njd_ice,    &
     &               icedat)

          call chkerr(stat)

          if(stat.lt.0) then

            call destroy('rdice   ',5,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          call chkstd(root)

          call outstd12(0,'    ',4,'       ',0_i8,0,0,1,0.e0,0,0,1,0.e0)

          call s_chkmxn('kai ',3,'[%]    ','undef',0_i8,stat,           &
     &                 0.e0,100.e0,0.e0,100.e0,nid_ice,njd_ice,1,icedat)

          call chkerr(stat)

          if(stat.lt.0) then

            call destroy('chkmxn  ',6,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          call chkstd(root)

          call undefice(stat,nid_ice,njd_ice,icedat,itmp1)

          call chkerr(stat)

          if(stat.lt.0) then

            call destroy('undefice',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          call chkstd(root)

        end if

! -----

!!!! Create the model input file for each processor element.

        do mype=0,npe-1

! Get the current processor element number.

          call currpe('all     ',3,'unset')

! -----

!!! For the soil surface data.

          if(sfcdat(1:1).eq.'o') then

! Calculate the real indices at the model grid points in the data
! region.

            call getrij(idmpopt_lnd,idnspol_lnd,idtlon_lnd,             &
     &                  iddxiv_lnd,iddyiv_lnd,'surface ',7,'xx',        &
     &                  x0,y0,cpj,x0lnd,y0lnd,cpjlnd,ni,nj,ri,rj,di,dj, &
     &                  tmp1,kai)

! -----

! Interpolate the land use data to the model grid.

            call hintlnd(idmpopt_lnd,ni,nj,ri,rj,land,nid_lnd,njd_lnd,  &
     &                   landat)

! -----

!! Correspond the soil surface data to the land use categories.

            if(intopt_lnd.eq.0) then

! Correspond the albedo.

              call cpondsfc(idnumctg_lnd,'  ',lnduse_lnd,albe_lnd,      &
     &                      0,ni+1,0,nj+1,land,albe)

! -----

! Correspond the evapotranspiration efficiency.

              call cpondsfc(idnumctg_lnd,'  ',lnduse_lnd,beta_lnd,      &
     &                      0,ni+1,0,nj+1,land,beta)

! -----

! Correspond the roughness length for velocity.

              call cpondsfc(idnumctg_lnd,'  ',lnduse_lnd,z0m_lnd,       &
     &                      0,ni+1,0,nj+1,land,z0m)

! -----

! Correspond the roughness length for scalar.

              call cpondsfc(idnumctg_lnd,'  ',lnduse_lnd,z0h_lnd,       &
     &                      0,ni+1,0,nj+1,land,z0h)

! -----

! Correspond the thermal capacity.

              call cpondsfc(idnumctg_lnd,'  ',lnduse_lnd,cap_lnd,       &
     &                      0,ni+1,0,nj+1,land,cap)

! -----

! Correspond the thermal diffusivity.

              call cpondsfc(idnumctg_lnd,'  ',lnduse_lnd,nuu_lnd,       &
     &                      0,ni+1,0,nj+1,land,nuu)

! -----

!! -----

!! Interpolate the soil surface data to the model grid horizontally.

            else

! Interpolate the albedo.

              call cpondsfc(idnumctg_lnd,'oo',lnduse_lnd,albe_lnd,      &
     &                      1,nid_lnd,1,njd_lnd,landat,ltmp1)

              call hint2d(idmpopt_lnd,idoneopt,'  ','cal ',             &
     &                    ni,nj,ri,rj,di,dj,albe,nid_lnd,njd_lnd,ltmp1)

              call estimsfc(idnumctg_lnd,lnduse_lnd,albe_lnd,ni,nj,land,&
     &                      albe)

! -----

! Interpolate the evapotranspiration efficiency.

              call cpondsfc(idnumctg_lnd,'oo',lnduse_lnd,beta_lnd,      &
     &                      1,nid_lnd,1,njd_lnd,landat,ltmp1)

              call hint2d(idmpopt_lnd,idoneopt,'  ','skip',             &
     &                    ni,nj,ri,rj,di,dj,beta,nid_lnd,njd_lnd,ltmp1)

              call estimsfc(idnumctg_lnd,lnduse_lnd,beta_lnd,ni,nj,land,&
     &                      beta)

! -----

! Interpolate the roughness length for velocity.

              call cpondsfc(idnumctg_lnd,'oo',lnduse_lnd,z0m_lnd,       &
     &                      1,nid_lnd,1,njd_lnd,landat,ltmp1)

              call hint2d(idmpopt_lnd,idoneopt,'  ','skip',             &
     &                    ni,nj,ri,rj,di,dj,z0m,nid_lnd,njd_lnd,ltmp1)

              call estimsfc(idnumctg_lnd,lnduse_lnd,z0m_lnd,ni,nj,land, &
     &                      z0m)

! -----

! Interpolate the roughness length for scalar.

              call cpondsfc(idnumctg_lnd,'oo',lnduse_lnd,z0h_lnd,       &
     &                      1,nid_lnd,1,njd_lnd,landat,ltmp1)

              call hint2d(idmpopt_lnd,idoneopt,'  ','skip',             &
     &                    ni,nj,ri,rj,di,dj,z0h,nid_lnd,njd_lnd,ltmp1)

              call estimsfc(idnumctg_lnd,lnduse_lnd,z0h_lnd,ni,nj,land, &
     &                      z0h)

! -----

! Interpolate the thermal capacity.

              call cpondsfc(idnumctg_lnd,'oo',lnduse_lnd,cap_lnd,       &
     &                      1,nid_lnd,1,njd_lnd,landat,ltmp1)

              call hint2d(idmpopt_lnd,idoneopt,'  ','skip',             &
     &                    ni,nj,ri,rj,di,dj,cap,nid_lnd,njd_lnd,ltmp1)

              call estimsfc(idnumctg_lnd,lnduse_lnd,cap_lnd,ni,nj,land, &
     &                      cap)

! -----

! Interpolate the thermal diffusivity.

              call cpondsfc(idnumctg_lnd,'oo',lnduse_lnd,nuu_lnd,       &
     &                      1,nid_lnd,1,njd_lnd,landat,ltmp1)

              call hint2d(idmpopt_lnd,idoneopt,'  ','skip',             &
     &                    ni,nj,ri,rj,di,dj,nuu,nid_lnd,njd_lnd,ltmp1)

              call estimsfc(idnumctg_lnd,lnduse_lnd,nuu_lnd,ni,nj,land, &
     &                      nuu)

! -----

            end if

!! -----

          end if

!!! -----

!! For the sea ice distribution data.

          if(sfcdat(3:3).eq.'o') then

! Calculate the real indices at the model grid points in the data
! region.

            call getrij(idmpopt_ice,idnspol_ice,idtlon_ice,             &
     &                  iddxiv_ice,iddyiv_ice,'surface ',7,'  ',        &
     &                  x0,y0,cpj,x0ice,y0ice,cpjice,ni,nj,ri,rj,di,dj, &
     &                  tmp1,kai)

! -----

! Interpolate the sea ice distribution data to the model grid
! horizontally.

            call hint2d(idmpopt_ice,idoneopt,'  ','cal ',               &
     &                  ni,nj,ri,rj,di,dj,kai,nid_ice,njd_ice,icedat)

! -----

          end if

!! -----

! Read in the data to the interpolated surface file.

          call outsfc(idexprim,idcrsdir,idsfcdat,idncexp,idnccrs,       &
     &                idwlngth,ni,nj,land,albe,beta,z0m,z0h,cap,nuu,kai)

! -----

        end do

!!!! -----

! Read in the message to standard i/o.

        call outstd05(0)

! -----

      end if

!!!!! -----

!!! Create the model input interpolated sea surface temperature file.

      if(sfcdat(2:2).eq.'o') then

! Set the common used variable.

        sst103=1000_i8*int(sstitv+.1e0,i8)

! -----

! Calculate the number of steps of the main do loop.

        call sststep(idsfcopt,idsstitv,idstime,idetime,nstp0,nstp1)

! -----

        do it=nstp0,nstp1

! Calculate the current forecast date.

          ctime=sst103*(it-1_i8)

          call getdate(idate,ctime,cdate)

! -----

! Set the map projection parameters of the sea surface temperature data
! grid and read out the data from the sea surface temperature data file.

          stat=0

          call setproj(idmpopt_sst,idnspol_sst,iddx_sst,iddy_sst,       &
     &                 idulat_sst,idulon_sst,idriu_sst,idrju_sst,       &
     &                 idtlat1_sst,idtlat2_sst,idtlon_sst,              &
     &                 'surface ',7,x0sst,y0sst,cpjsst)

          call rdsst(iddatdir,idncdat,idwlngth,cdate,ctime,stat,        &
     &               nid_sst,njd_sst,sstdat)

          call chkerr(stat)

          if(stat.lt.0) then

            call destroy('rdsst   ',5,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          call chkstd(root)

          call outstd12(0,'    ',4,'       ',0_i8,0,0,1,0.e0,0,0,1,0.e0)

          call s_chkmxn('sst ',3,'[K]    ','undef',0_i8,stat,           &
     &                  268.16e0,323.16e0,268.16e0,323.16e0,            &
     &                  nid_sst,njd_sst,1,sstdat)

          call chkerr(stat)

          if(stat.lt.0) then

            call destroy('chkmxn  ',6,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          call chkstd(root)

          call undefsst(stat,nid_sst,njd_sst,sstdat,stmp1)

          call chkerr(stat)

          if(stat.lt.0) then

            call destroy('undefsst',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          call chkstd(root)

! -----

!! Create the model input file for each processor element.

          do mype=0,npe-1

! Get the current processor element number.

            call currpe('all     ',3,'unset')

! -----

! Calculate the real indices at the model grid points in the data
! region.

            call getrij(idmpopt_sst,idnspol_sst,idtlon_sst,             &
     &                  iddxiv_sst,iddyiv_sst,'surface ',7,'  ',        &
     &                  x0,y0,cpj,x0sst,y0sst,cpjsst,ni,nj,ri,rj,di,dj, &
     &                  tmp1,kai)

! -----

! Interpolate the sea surface temperature data to the model grid
! horizontally.

            call hint2d(idmpopt_sst,idoneopt,'  ','cal ',               &
     &                  ni,nj,ri,rj,di,dj,sst,nid_sst,njd_sst,sstdat)

! -----

! Read in the data to the interpolated sea surface temperature file.

            call outsst(idexprim,idcrsdir,idncexp,idnccrs,idwlngth,     &
     &                  it,nstp0,ctime,ni,nj,sst)

! -----

          end do

!! -----

! Read in the message to standard i/o.

          call outstd05(0)

! -----

        end do

      end if

!!! -----

      end subroutine s_sfcdrv

!-----7--------------------------------------------------------------7--

      end module m_sfcdrv
