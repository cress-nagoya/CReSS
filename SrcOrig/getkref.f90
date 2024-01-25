!***********************************************************************
      module m_getkref
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/08/01
!     Modification: 2004/08/31, 2004/09/10, 2005/02/10, 2007/01/31,
!                   2007/06/27, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2008/12/11, 2009/02/27, 2009/11/13, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the index of the base state reference pressure.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getkref, s_getkref

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getkref

        module procedure s_getkref

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getkref(fprefsfc_gpv,idstr,idend,jdstr,jdend,kref,   &
     &                     nid,njd,nkd,zlow,zdat,altmin,nk,z)
!***********************************************************************

! Input variables

      integer, intent(in) :: fprefsfc_gpv
                       ! Formal parameter of unique index of refsfc_gpv

      integer, intent(in) :: idstr
                       ! Minimum index of model grid in data region
                       ! in x direction

      integer, intent(in) :: idend
                       ! Maximum index of model grid in data region
                       ! in x direction

      integer, intent(in) :: jdstr
                       ! Minimum index of model grid in data region
                       ! in y direction

      integer, intent(in) :: jdend
                       ! Maximum index of model grid in data region
                       ! in y direction

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zlow(1:nid,1:njd)
                       ! Lowest z physical coordinates in data

      real, intent(in) :: zdat(1:nid,1:njd,1:nkd)
                       ! z physical coordinates in data

      real, intent(in) :: z(1:nk)
                       ! zeta coordinates

! Output variable

      integer, intent(out) :: kref
                       ! Reference index

! Internal shared variables

      integer refsfc_gpv
                       ! Option for
                       ! surface data reference in interpolating

      integer kdbot    ! Data index of lowest interpolated plane

      integer kbot     ! Index of lowest interpolated plane
      integer ktop     ! Index of highest interpolated plane

      real zdmax       ! Highest height at lowest data plane
      real zdmin       ! Lowest height at highest data plane

      real, intent(inout) :: altmin(1:nkd)
                       ! Lowest altitude above ground level
                       ! for each data plane

! Internal private variables

      integer id       ! Data array index in x direction
      integer jd       ! Data array index in y direction
      integer kd       ! Data array index in z direction

      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fprefsfc_gpv,refsfc_gpv)

! -----

! Initialize the data index of lowest interpolated plane.

      if(refsfc_gpv.eq.0) then
        kdbot=1
      else
        kdbot=2
      end if

! -----

! Initialize the processed variables.

      kbot=nk
      ktop=1

      zdmin=lim36
      zdmax=lim36n

! -----

!! Get the index of lowest and highest interpolated plane.

!$omp parallel default(shared) private(kd)

! Reset the data index of lowest interpolated plane.

      if(refsfc_gpv.eq.1) then

!$omp do schedule(runtime)

        do kd=2,nkd-1
          altmin(kd)=lim36
        end do

!$omp end do

        if(idstr.le.idend) then

          do kd=2,nkd-1

!$omp do schedule(runtime) private(id,jd)

            do jd=jdstr,jdend
            do id=idstr,idend
              altmin(kd)=min(altmin(kd),zdat(id,jd,kd)-zlow(id,jd))
            end do
            end do

!$omp end do

          end do

        else

          do kd=2,nkd-1

!$omp do schedule(runtime) private(id,jd)

            do jd=jdstr,jdend
            do id=idstr,nid
              altmin(kd)=min(altmin(kd),zdat(id,jd,kd)-zlow(id,jd))
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(id,jd)

            do jd=jdstr,jdend
            do id=1,idend
              altmin(kd)=min(altmin(kd),zdat(id,jd,kd)-zlow(id,jd))
            end do
            end do

!$omp end do

          end do

        end if

!$omp single

        do_kd: do kd=2,nkd-1

          if(altmin(kd).gt.0.e0) then

            kdbot=kd

            exit do_kd

          end if

        end do do_kd

!$omp end single

      end if

! -----

! Get the highest height at lowest data plane and the lowest height at
! highest data plane.

      if(idstr.le.idend) then

!$omp do schedule(runtime) private(id,jd)                               &
!$omp&   reduction(min: zdmin) reduction(max: zdmax)

        do jd=jdstr,jdend
        do id=idstr,idend
          zdmin=min(zdat(id,jd,nkd),zdmin)
          zdmax=max(zdat(id,jd,kdbot),zdmax)
        end do
        end do

!$omp end do

      else

!$omp do schedule(runtime) private(id,jd)                               &
!$omp&   reduction(min: zdmin) reduction(max: zdmax)

        do jd=jdstr,jdend
        do id=idstr,nid
          zdmin=min(zdat(id,jd,nkd),zdmin)
          zdmax=max(zdat(id,jd,kdbot),zdmax)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(id,jd)                               &
!$omp&   reduction(min: zdmin) reduction(max: zdmax)

        do jd=jdstr,jdend
        do id=1,idend
          zdmin=min(zdat(id,jd,nkd),zdmin)
          zdmax=max(zdat(id,jd,kdbot),zdmax)
        end do
        end do

!$omp end do

      end if

! -----

! Calculate the index of lowest and highest interpolated plane.

!$omp single private(k)

      do_k_1: do k=1,nk-1

        if(z(k).gt.zdmax) then

          kbot=k

          exit do_k_1

        end if

      end do do_k_1

      do_k_2: do k=2,nk

        if(z(k).gt.zdmin) then

          ktop=k-1

          exit do_k_2

        end if

      end do do_k_2

!$omp end single

! -----

!$omp end parallel

!! -----

! Finally get the index of the base state reference pressure.

      if(z(nk).lt.zdmin) then
        ktop=nk
      end if

      if(kbot.gt.ktop) then
        kref=1
      else
        kref=kbot
      end if

! -----

      end subroutine s_getkref

!-----7--------------------------------------------------------------7--

      end module m_getkref
