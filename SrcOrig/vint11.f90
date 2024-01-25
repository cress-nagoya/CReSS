!***********************************************************************
      module m_vint11
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/01/15
!     Modification: 2001/04/15, 2001/05/29, 2001/06/29, 2002/04/02,
!                   2002/12/02, 2003/04/30, 2003/05/19, 2004/09/10,
!                   2007/01/31, 2007/06/27, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     interpolate the 1 dimensional data to the fine horizontally
!     averaged levels vertically.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vint11, s_vint11

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vint11

        module procedure s_vint11

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
      subroutine s_vint11(nlev,z1d,var1d,nkd,zref,varef)
!***********************************************************************

! Input variables

      integer, intent(in) :: nlev
                       ! Horizontally averaged levels

      integer, intent(in) :: nkd
                       ! Data dimension

      real, intent(in) :: z1d(1:nlev)
                       ! Horizontally averaged z physical coordinates

      real, intent(in) :: zref(1:nkd)
                       ! z physical coordinates in data

      real, intent(in) :: varef(1:nkd)
                       ! Optional variable in data

! Output variable

      real, intent(out) :: var1d(1:nlev)
                       ! Horizontally averaged interpolated variable

! Internal shared variables

      integer nkdm1    ! nkd - 1

      real dkivl       ! Inverse of distance in z direction
                       ! between data points in lowest layer

      real dkivh       ! Inverse of distance in z direction
                       ! between data points in highest layer

! Internal private variables

      integer kl       ! Averaged array index in z direction

      integer kd       ! Data array index in z direction

      real dk          ! Distance in z direction
                       ! between averaged levels and data points

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      nkdm1=nkd-1

      dkivl=1.e0/(zref(nkd)-zref(nkdm1))
      dkivh=1.e0/(zref(2)-zref(1))

! -----

!! Interpolate the 1 dimensional data to the fine horizontally averaged
!! levels vertically.

!$omp parallel default(shared) private(kd)

! Extrapolate the data to the averaged levels vertically.

!$omp do schedule(runtime) private(kl,dk)

      do kl=1,nlev

        if(zref(1).gt.z1d(kl)) then

          dk=dkivl*(z1d(kl)-zref(1))

          var1d(kl)=(1.e0-dk)*varef(1)+dk*varef(2)

        else if(zref(nkd).le.z1d(kl)) then

          dk=dkivh*(z1d(kl)-zref(nkdm1))

          var1d(kl)=(1.e0-dk)*varef(nkdm1)+dk*varef(nkd)

        end if

      end do

!$omp end do

! -----

! Interpolate the data to the averaged levels vertically.

      do kd=1,nkd-1

!$omp do schedule(runtime) private(kl,dk)

        do kl=1,nlev

          if(zref(kd+1).gt.z1d(kl).and.zref(kd).le.z1d(kl)) then

            dk=(z1d(kl)-zref(kd))/(zref(kd+1)-zref(kd))

            var1d(kl)=(1.e0-dk)*varef(kd)+dk*varef(kd+1)

          end if

        end do

!$omp end do

      end do

! -----

!$omp end parallel

!! -----

      end subroutine s_vint11

!-----7--------------------------------------------------------------7--

      end module m_vint11
