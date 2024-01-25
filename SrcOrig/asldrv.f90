!***********************************************************************
      module m_asldrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2011/11/10, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures to interpolate the aerosol data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_aslstep
      use m_chkasl
      use m_chkerr
      use m_chkstd
      use m_comindx
      use m_comkind
      use m_commpi
      use m_copy3d
      use m_currpe
      use m_destroy
      use m_getcname
      use m_getdate
      use m_getrname
      use m_getxy
      use m_getz
      use m_getzph
      use m_inichar
      use m_intrpasl
      use m_outasl
      use m_outstd05
      use m_rdaero
      use m_setproj
      use m_vint31a
      use m_xy2ll

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: asldrv, s_asldrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface asldrv

        module procedure s_asldrv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic int

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_asldrv(fpidate,fpaslitv,ni,nj,nk,nqa,                &
     &                  z,zph,qasl,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,       &
     &                  nid_asl,njd_asl,nkd_asl,km_asl,zdat,qadat,dtmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpidate
                       ! Formal parameter of unique index of idate

      integer, intent(in) :: fpaslitv
                       ! Formal parameter of unique index of aslitv

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

      integer, intent(in) :: nid_asl
                       ! Aerosol data dimension in x direction

      integer, intent(in) :: njd_asl
                       ! Aerosol data dimension in y direction

      integer, intent(in) :: nkd_asl
                       ! Aerosol data dimension in z direction

      integer, intent(in) :: km_asl
                       ! Dimension of max(nk, nkd_asl)

! Internal shared variables

      character(len=108) idate
                       ! Forecast start date
                       ! with Gregorian calendar, yyyymmddhhmm

      character(len=12) cdate
                       ! Current forecast date
                       ! with Gregorian calendar, yyyymmddhhmm

      integer n        ! Array index in 4th direction

      integer(kind=i8) it
                       ! Index of main do loop

      integer(kind=i8) nstp0
                       ! Start index of main do loop

      integer(kind=i8) nstp1
                       ! End index of main do loop

      integer(kind=i8) ctime
                       ! Model current forecast time

      integer(kind=i8) asl103
                       ! 1000 x int(aslitv + 0.1)

      integer stat     ! Runtime status

      real aslitv      ! Time interval of aerosol data file

      real x0          ! x origin of model grid
      real y0          ! y origin of model grid

      real cpj(1:7)    ! Map projection parameters of model grid

      real x0asl       ! x origin of aerosol data grid
      real y0asl       ! y origin of aerosol data grid

      real cpjasl(1:7) ! Map projection parameters of aerosol data grid

      real, intent(inout) :: z(1:nk)
                       ! zeta coordinates

      real, intent(inout) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(inout) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp6(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: zdat(1:nid_asl,1:njd_asl,1:nkd_asl)
                       ! z physical coordinates in data

      real, intent(inout) ::                                            &
     &                   qadat(1:nid_asl,1:njd_asl,1:km_asl,1:nqa(0))
                       ! Aerosol mixing ratio in data

      real, intent(inout) :: dtmp1(1:nid_asl,1:njd_asl,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Initialize character variable.

      call inichar(idate)

! -----

! Get the required namelist variables.

      call getcname(fpidate,idate)
      call getrname(fpaslitv,aslitv)

! -----

! Read in the message to standard i/o.

      call outstd05(0)

! -----

! Set the map projection parameters of the model grid.

      call setproj(idmpopt,idnspol,iddx,iddy,idulat,idulon,idriu,idrju, &
     &             idtlat1,idtlat2,idtlon,'solver  ',6,x0,y0,cpj)

! -----

! Set the map projection parameters of the data grid.

      call setproj(idmpopt_asl,idnspol_asl,iddx_asl,iddy_asl,idulat_asl,&
     &           idulon_asl,idriu_asl,idrju_asl,idtlat1_asl,idtlat2_asl,&
     &           idtlon_asl,'asldata ',7,x0asl,y0asl,cpjasl)

! -----

! Get the x and the y coordinates at the data grid points.

      call s_getxy(iddx_asl,iddy_asl,'oo',1,nid_asl,1,njd_asl,          &
     &             dtmp1(1,1,1),dtmp1(1,1,2))

! -----

! Get the zeta coordinates at the model grid points.

      call getz(iddz,idzsfc,nk,z)

! -----

! Calculate the number of steps of the main do loop.

      call aslstep(idnggopt,idexbopt,idlspopt,idvspopt,idaslitv,        &
     &             idstime,idetime,nstp0,nstp1)

! -----

!!! Create the model input files.

! Set the common used variable.

      asl103=1000_i8*int(aslitv+.1e0,i8)

! -----

      do it=nstp0,nstp1

! Calculate the current forecast date.

        ctime=asl103*(it-1_i8)

        call getdate(idate,ctime,cdate)

! -----

! Read out the data from the aerosol data file.

        stat=0

        call rdaero(iddatdir,idncdat,idwlngth,cdate,ctime,stat,         &
     &              nid_asl,njd_asl,nkd_asl,nqa,zdat,qadat)

        call chkerr(stat)

        if(stat.lt.0) then

          call destroy('rdaero  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        call chkstd(root)

! -----

! Check the aerosol data variables.

        stat=0

        call chkasl(ctime,stat,nid_asl,njd_asl,nkd_asl,nqa,zdat,qadat)

        call chkerr(stat)

        if(stat.lt.0) then

          call destroy('chkasl  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        call chkstd(root)

! -----

! Interpolate the aerosol data to the flat plane vertically.

        do n=1,nqa(0)

          call s_vint31a(nid_asl,njd_asl,nkd_asl,zdat,qadat(1,1,1,n),   &
     &                   nk,z,dtmp1)

          call s_copy3d(1,nid_asl,1,njd_asl,1,nk,dtmp1,qadat(1,1,1,n))

        end do

! -----

!! Create the model input file for each processor element.

        do mype=0,npe-1

! Calculate the current processor element number.

         call currpe('all     ',3,'unset')

! -----

! Calculate the z physical coordinates.

         call s_getzph(idtrnopt,idexbopt,ni,nj,nk,z,zph,tmp1,tmp5,tmp6)

! -----

! Interpolate the aerosol data to the model grid points horizontally.

         call intrpasl(x0,y0,cpj,x0asl,y0asl,cpjasl,ni,nj,nk,nqa,       &
     &                 z,zph,qasl,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,        &
     &                 nid_asl,njd_asl,qadat)

! -----

! Read in the data to the interpolated aerosol file.

         call outasl(idexprim,idcrsdir,idncexp,idnccrs,idwlngth,        &
     &               it,nstp0,ctime,ni,nj,nk,nqa,qasl)

! -----

        end do

!! -----

! Read in the message to standard i/o.

        call outstd05(0)

! -----

      end do

!!! -----

      end subroutine s_asldrv

!-----7--------------------------------------------------------------7--

      end module m_asldrv
