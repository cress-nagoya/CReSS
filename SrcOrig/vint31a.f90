!***********************************************************************
      module m_vint31a
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     interpolate the 3 dimensional input variable to the 1 dimensional
!     flat plane.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vint31a, s_vint31a

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vint31a

        module procedure s_vint31a

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
      subroutine s_vint31a(nid,njd,nkd,zdat,vardat,nk,z,varef)
!***********************************************************************

! Input variables

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zdat(1:nid,1:njd,1:nkd)
                       ! z physical coordinates in data

      real, intent(in) :: vardat(1:nid,1:njd,1:nkd)
                       ! Optional variable in data

      real, intent(in) :: z(1:nk)
                       ! zeta coordinates

! Output variable

      real, intent(out) :: varef(1:nid,1:njd,1:nk)
                       ! Optional vertically interpolated variable

! Internal shared variable

      integer nkdm1    ! nkd - 1

! Internal private variables

      integer id       ! Data array index in x direction
      integer jd       ! Data array index in y direction
      integer kd       ! Data array index in z direction

      integer k        ! Array index in z direction

      real dk          ! Distance in z direction
                       ! between flat plane and data points

!-----7--------------------------------------------------------------7--

! Set the common used variable.

      nkdm1=nkd-1

! -----

!! Interpolate the variable to the flat plane vertically.

!$omp parallel default(shared) private(k,kd)

! Extrapolate the variable.

      do k=1,nk

!$omp do schedule(runtime) private(id,jd,dk)

        do jd=1,njd
        do id=1,nid

          if(zdat(id,jd,1).gt.z(k)) then

            varef(id,jd,k)=vardat(id,jd,1)

          end if

          if(zdat(id,jd,nkd).le.z(k)) then

            dk=(z(k)-zdat(id,jd,nkdm1))                                 &
     &        /(zdat(id,jd,nkd)-zdat(id,jd,nkdm1))

            varef(id,jd,k)                                              &
     &        =(1.e0-dk)*vardat(id,jd,nkdm1)+dk*vardat(id,jd,nkd)

          end if

        end do
        end do

!$omp end do

      end do

! -----

! Interpolate the variable.

      do kd=1,nkd-1

        do k=1,nk

!$omp do schedule(runtime) private(id,jd,dk)

          do jd=1,njd
          do id=1,nid

            if(zdat(id,jd,kd).le.z(k).and.zdat(id,jd,kd+1).gt.z(k)) then

              dk=(z(k)-zdat(id,jd,kd))/(zdat(id,jd,kd+1)-zdat(id,jd,kd))

              varef(id,jd,k)                                            &
     &          =(1.e0-dk)*vardat(id,jd,kd)+dk*vardat(id,jd,kd+1)

            end if

          end do
          end do

!$omp end do

        end do

      end do

! -----

!$omp end parallel

!! -----

      end subroutine s_vint31a

!-----7--------------------------------------------------------------7--

      end module m_vint31a
