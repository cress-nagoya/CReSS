!***********************************************************************
      module m_buoywsi
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/01
!     Modification: 1999/11/19, 2000/01/17, 2000/07/05, 2001/06/29,
!                   2001/11/20, 2001/12/11, 2002/04/02, 2002/08/15,
!                   2003/01/04, 2003/04/30, 2003/05/19, 2003/11/05,
!                   2003/12/26, 2004/04/15, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the buoyancy in the small time steps integration for the
!     horizontally explicit and vertically implicit method.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: buoywsi, s_buoywsi

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface buoywsi

        module procedure s_buoywsi

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
      subroutine s_buoywsi(fpgwmopt,fpweicoe,dts,ni,nj,nk,rbr,ptbr,rst, &
     &                     rcsq,pp,ptp,fp,fw,wb8s)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgwmopt
                       ! Formal parameter of unique index of gwmopt

      integer, intent(in) :: fpweicoe
                       ! Formal parameter of unique index of weicoe

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: rcsq(0:ni+1,0:nj+1,1:nk)
                       ! rbr x sound wave speed squared

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(in) :: fp(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in pressure equation

! Input and output variable

      real, intent(inout) :: fw(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

! Internal shared variables

      integer gwmopt   ! Option for gravity wave mode integration

      real weicoe      ! Weighting coefficient for implicit method

      real dtw         ! dts x weicoe

      real g05n        ! - 0.5 x g

      real, intent(inout) :: wb8s(0:ni+1,0:nj+1,1:nk)
                       ! 0.5 x buoyancy value
                       ! in small time steps at scalar points

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpgwmopt,gwmopt)
      call getrname(fpweicoe,weicoe)

! -----

! Set the common used variables.

      dtw=dts*weicoe

      g05n=-.5e0*g

! -----

! Calculate the buoyancy in the small time steps.

!$omp parallel default(shared) private(k)

      if(gwmopt.eq.0) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            wb8s(i,j,k)=(pp(i,j,k)*rst(i,j,k)+fp(i,j,k)*rbr(i,j,k)*dtw) &
     &        /rcsq(i,j,k)*g05n
          end do
          end do

!$omp end do

        end do

      else

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            wb8s(i,j,k)=((pp(i,j,k)*rst(i,j,k)+fp(i,j,k)*rbr(i,j,k)*dtw)&
     &        /rcsq(i,j,k)-ptp(i,j,k)*rst(i,j,k)/ptbr(i,j,k))*g05n
          end do
          end do

!$omp end do

        end do

      end if

      do k=3,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          fw(i,j,k)=fw(i,j,k)+(wb8s(i,j,k-1)+wb8s(i,j,k))
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_buoywsi

!-----7--------------------------------------------------------------7--

      end module m_buoywsi
