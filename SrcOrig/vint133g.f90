!***********************************************************************
      module m_vint133g
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/07/15
!     Modification: 2003/04/30, 2003/05/19, 2003/12/12, 2004/09/10,
!                   2006/09/21, 2007/06/27, 2007/10/19, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     interpolate the variable to the model grid vertically.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getindx

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vint133g, s_vint133g

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vint133g

        module procedure s_vint133g

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
      subroutine s_vint133g(xo,ni,nj,nk,zph,outvar,nlev,z1d,invar)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nlev
                       ! Horizontally averaged vertical dimension

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: z1d(1:nlev)
                       ! Horizontally averaged z physical coordinates

      real, intent(in) :: invar(0:ni+1,0:nj+1,1:nlev)
                       ! Optional variable
                       ! at horizontally averraged plane

! Output variable

      real, intent(out) :: outvar(0:ni+1,0:nj+1,1:nk)
                       ! Optional interpolated variable

! Internal shared variables

      integer nlevm1   ! nlev - 1

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction
      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction
      integer kend     ! Maximum do loops index in z direction

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
                       ! between model and averaged points

!-----7--------------------------------------------------------------7--

! Get the maximum and minimim indices of do loops.

      call getindx(xo,0,ni+1,0,nj+1,istr,iend,jstr,jend)

      if(xo(3:3).eq.'o') then
        kend=nk-1
      else if(xo(3:3).eq.'x') then
        kend=nk-2
      end if

! -----

! Set the common used variables.

      nlevm1=nlev-1

      dkivl=1.e0/(z1d(2)-z1d(1))
      dkivh=1.e0/(z1d(nlev)-z1d(nlevm1))

! -----

!! Interpolate the variable to the model grid vertically.

!$omp parallel default(shared) private(k,kl)

! Extrapolate the variable.

      do k=2,kend

!$omp do schedule(runtime) private(i,j,dk)

        do j=jstr,jend
        do i=istr,iend

          if(z1d(1).gt.zph(i,j,k)) then

            dk=dkivl*(zph(i,j,k)-z1d(1))

            outvar(i,j,k)=(1.e0-dk)*invar(i,j,1)+dk*invar(i,j,2)

          else if(z1d(nlev).le.zph(i,j,k)) then

            dk=dkivh*(zph(i,j,k)-z1d(nlevm1))

            outvar(i,j,k)=(1.e0-dk)*invar(i,j,nlevm1)+dk*invar(i,j,nlev)

          end if

        end do
        end do

!$omp end do

      end do

! -----

! Interpolate the variable.

      do kl=1,nlev-1

        do k=2,kend

!$omp do schedule(runtime) private(i,j,dk)

          do j=jstr,jend
          do i=istr,iend

            if(z1d(kl).le.zph(i,j,k).and.z1d(kl+1).gt.zph(i,j,k)) then

              dk=(zph(i,j,k)-z1d(kl))/(z1d(kl+1)-z1d(kl))

              outvar(i,j,k)=(1.e0-dk)*invar(i,j,kl)+dk*invar(i,j,kl+1)

            end if

          end do
          end do

!$omp end do

        end do

      end do

! -----

!$omp end parallel

!! -----

      end subroutine s_vint133g

!-----7--------------------------------------------------------------7--

      end module m_vint133g
