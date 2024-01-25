!***********************************************************************
      module m_vint13
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/03/24
!     Modification: 1999/04/06, 1999/05/10, 1999/05/20, 1999/06/28,
!                   1999/07/05, 1999/10/12, 1999/11/01, 2000/01/17,
!                   2001/01/15, 2001/04/15, 2001/05/29, 2001/06/29,
!                   2001/11/20, 2002/04/02, 2002/06/06, 2002/06/18,
!                   2002/09/09, 2003/04/30, 2003/05/19, 2004/01/09,
!                   2004/09/10, 2005/04/04, 2006/09/21, 2007/06/27,
!                   2007/10/19, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     interpolate the variable to the model or data grid vertically.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getindx

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vint13, s_vint13

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vint13

        module procedure s_vint13

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
      subroutine s_vint13(xo,imin,imax,jmin,jmax,kmin,kmax,zph,outvar,  &
     &                    nlev,z1d,var1d)
!***********************************************************************

! Input variables

      character(len=2), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: imin
                       ! Minimum array index in x direction

      integer, intent(in) :: imax
                       ! Maximum array index in x direction

      integer, intent(in) :: jmin
                       ! Minimum array index in y direction

      integer, intent(in) :: jmax
                       ! Maximum array index in y direction

      integer, intent(in) :: kmin
                       ! Minimum array index in z direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      integer, intent(in) :: nlev
                       ! Horizontally averaged vertical dimension

      real, intent(in) :: zph(imin:imax,jmin:jmax,kmin:kmax)
                       ! z physical coordinates in model or data

      real, intent(in) :: z1d(0:nlev)
                       ! Horizontally averaged z physical coordinates

      real, intent(in) :: var1d(0:nlev)
                       ! Horizontally averaged optional variable

! Output variable

      real, intent(out) :: outvar(imin:imax,jmin:jmax,kmin:kmax)
                       ! Optional interpolated variable

! Internal shared variables

      integer nlevm1   ! nlev - 1

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction
      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction

      real dkivl       ! Inverse of distance in z direction
                       ! between averaged levels in lowest layer

      real dkivh       ! Inverse of distance in z direction
                       ! between averaged levels in highest layer

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer kl       ! Averaged array index in z direction

      real dk          ! Distance in z direction
                       ! between model or data and averaged points

!-----7--------------------------------------------------------------7--

! Get the maximum and minimim indices of do loops.

      call getindx(xo,imin,imax,jmin,jmax,istr,iend,jstr,jend)

! -----

! Set the common used variables.

      nlevm1=nlev-1

      dkivl=1.e0/(z1d(1)-z1d(0))
      dkivh=1.e0/(z1d(nlev)-z1d(nlevm1))

! -----

!! Interpolate the variable to the model or data grid vertically.

!$omp parallel default(shared) private(k,kl)

! Extrapolate the variable.

      do k=kmin,kmax

!$omp do schedule(runtime) private(i,j,dk)

        do j=jstr,jend
        do i=istr,iend

          if(z1d(0).gt.zph(i,j,k)) then

            dk=dkivl*(zph(i,j,k)-z1d(0))

            outvar(i,j,k)=(1.e0-dk)*var1d(0)+dk*var1d(1)

          else if(z1d(nlev).le.zph(i,j,k)) then

            dk=dkivh*(zph(i,j,k)-z1d(nlevm1))

            outvar(i,j,k)=(1.e0-dk)*var1d(nlevm1)+dk*var1d(nlev)

          end if

        end do
        end do

!$omp end do

      end do

! -----

! Interpolate the variable.

      do kl=0,nlev-1

        do k=kmin,kmax

!$omp do schedule(runtime) private(i,j,dk)

          do j=jstr,jend
          do i=istr,iend

            if(z1d(kl).le.zph(i,j,k).and.z1d(kl+1).gt.zph(i,j,k)) then

              dk=(zph(i,j,k)-z1d(kl))/(z1d(kl+1)-z1d(kl))

              outvar(i,j,k)=(1.e0-dk)*var1d(kl)+dk*var1d(kl+1)

            end if

          end do
          end do

!$omp end do

        end do

      end do

! -----

!$omp end parallel

!! -----

      end subroutine s_vint13

!-----7--------------------------------------------------------------7--

      end module m_vint13
