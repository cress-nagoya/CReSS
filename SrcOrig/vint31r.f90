!***********************************************************************
      module m_vint31r
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2003/04/30, 2003/05/19, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     interpolate the 3 dimensional input variable to the 1 dimensional
!     flat plane.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vint31r, s_vint31r

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vint31r

        module procedure s_vint31r

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
      subroutine s_vint31r(nid,njd,nkd,zdat,vardat,nk,z,varef)
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

! Internal private variables

      integer id       ! Data array index in x direction
      integer jd       ! Data array index in y direction
      integer kd       ! Data array index in z direction

      integer k        ! Array index in z direction

      real dk          ! Distance in z direction
                       ! between flat plane and data points

!-----7--------------------------------------------------------------7--

!! Interpolate the 3 dimensional input variable to the 1 dimensional
!! flat plane.

!$omp parallel default(shared) private(k,kd)

! Fill in the undifined value outside of the flat plane.

      do k=1,nk

!$omp do schedule(runtime) private(id,jd)

        do jd=1,njd
        do id=1,nid

          if(zdat(id,jd,1).gt.z(k).or.zdat(id,jd,nkd).le.z(k)) then

            varef(id,jd,k)=lim35n

          end if

        end do
        end do

!$omp end do

      end do

! -----

! Interpolate the variable to the flat plane vertically.

      do kd=1,nkd-1

        do k=1,nk

!$omp do schedule(runtime) private(id,jd,dk)

          do jd=1,njd
          do id=1,nid

            if(zdat(id,jd,kd).le.z(k).and.zdat(id,jd,kd+1).gt.z(k)) then

              if(vardat(id,jd,kd).gt.lim34n                             &
     &          .and.vardat(id,jd,kd+1).gt.lim34n) then

                dk=(z(k)-zdat(id,jd,kd))                                &
     &            /(zdat(id,jd,kd+1)-zdat(id,jd,kd))

                varef(id,jd,k)                                          &
     &            =(1.e0-dk)*vardat(id,jd,kd)+dk*vardat(id,jd,kd+1)

              else

                varef(id,jd,k)=lim35n

              end if

            end if

          end do
          end do

!$omp end do

        end do

      end do

! -----

!$omp end parallel

!! -----

      end subroutine s_vint31r

!-----7--------------------------------------------------------------7--

      end module m_vint31r
