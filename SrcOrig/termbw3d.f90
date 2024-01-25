!***********************************************************************
      module m_termbw3d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/08/08
!     Modification: 2006/09/30, 2007/10/19, 2008/01/11, 2008/05/02,
!                   2008/08/25, 2008/10/10, 2009/02/27, 2011/03/18,
!                   2011/08/18, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the terminal velocity for water bin.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: termbw3d, s_termbw3d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface termbw3d

        module procedure s_termbw3d

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_termbw3d(ncp,ni,nj,nk,nqw,rbw,rrbw,rbr,rbv,t,ubwmax, &
     &                      c1,c2,c3,c4,c5,ubw)
!***********************************************************************

! Input variables

      integer, intent(in) :: ncp
                       ! Index of current processed bin

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      real, intent(in) :: rbw(1:nqw)
                       ! Standard radius
                       ! between adjacent water bins [cm]

      real, intent(in) :: rrbw(1:nqw,1:5)
                       ! Related parameters of rbw

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Bese state density [g/cm^3]

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of bese state density [cm^3/g]

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

! Input and output variables

      real, intent(inout) :: ubwmax(0:ni+1,0:nj+1,1:nk)
                       ! Maximum terminal velocity of water bin [cm/s]

      real, intent(inout) :: c1(0:ni+1,0:nj+1,1:nk)
                       ! Common used temporary array

      real, intent(inout) :: c2(0:ni+1,0:nj+1,1:nk)
                       ! Common used temporary array

      real, intent(inout) :: c3(0:ni+1,0:nj+1,1:nk)
                       ! Common used temporary array

      real, intent(inout) :: c4(0:ni+1,0:nj+1,1:nk)
                       ! Common used temporary array

      real, intent(inout) :: c5(0:ni+1,0:nj+1,1:nk)
                       ! Common used temporary array

! Output variable

      real, intent(out) :: ubw(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of water bin [cm/s]

! Internal shared variables

      real g163        ! 1600.0 x g / 3.0

      real rhow03      ! 0.001 x rhow

      real rbwm2       ! rbwmax x rbwmax
      real rbwmiv      ! 1.0 / rbwmax

      real eta05       ! 0.5 x eta0
      real etav48      ! 48.0 / eta0

      real geta29      ! 200.0 x g / (9.0 x eta0)
      real geta4v      ! 0.01 / (g x eta0 x eta0 x eta0 x eta0)

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real a           ! Temporary variable
      real b           ! Temporary variable
      real c           ! Temporary variable
      real d           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Set the common used variable.

      g163=1600.e0*oned3*g

      rhow03=.001e0*rhow

      rbwm2=rbwmax*rbwmax
      rbwmiv=1.e0/rbwmax

      eta05=.5e0*eta0
      etav48=48.e0/eta0

      geta29=200.e0*g/(9.e0*eta0)
      geta4v=.01e0/(g*eta0*eta0*eta0*eta0)

! -----

!! Calculate the terminal velocity for water bin.

!$omp parallel default(shared) private(k)

! Set the common used variable.

      if(ncp.eq.1) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,a,b,c,d)

          do j=1,nj-1
          do i=1,ni-1

            a=rhow03-rbr(i,j,k)
            b=76.1e0-.155e0*(t(i,j,k)-t0)

            c=exp(oned6*log(b*b*b*geta4v*rbr(i,j,k)*rbr(i,j,k)/a))

            c1(i,j,k)=geta29*a
            c2(i,j,k)=etav48*c1(i,j,k)*rbr(i,j,k)
            c3(i,j,k)=eta05*rbv(i,j,k)
            c4(i,j,k)=g163*a*c/b
            c5(i,j,k)=c*c3(i,j,k)

            a=log(rbwm2*c4(i,j,k))

            b=a*a
            c=a*b

            d=5.23778e0*a-5.00015e0+.475294e0*c-2.04914e0*b             &
     &        +.238449e-2*b*c-.542819e-1*b*b

            ubwmax(i,j,k)=rbwmiv*exp(d)*c5(i,j,k)

          end do
          end do

!$omp end do

        end do

      end if

! -----

! Get the terminal velocity.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j,a,b,c,d)

        do j=1,nj-1
        do i=1,ni-1

          if(rbw(ncp).lt.1.e-3) then

            ubw(i,j,k)=c1(i,j,k)*rrbw(ncp,4)

          else if(rbw(ncp).ge.1.e-3.and.rbw(ncp).lt.5.35e-2) then

            a=log(c2(i,j,k)*rrbw(ncp,2))

            b=a*a
            c=a*b

            d=.992696e0*a-3.18657e0-.987059e-3*c-.153193e-2*b           &
     &        +.855176e-4*b*c-.578878e-3*b*b-.327815e-5*c*c

            ubw(i,j,k)=exp(d)*c3(i,j,k)*rrbw(ncp,5)

          else if(rbw(ncp).ge.5.35e-2.and.rbw(ncp).lt..35e0) then

            a=log(c4(i,j,k)*rrbw(ncp,1))

            b=a*a
            c=a*b

            d=5.23778e0*a-5.00015e0+.475294e0*c-2.04914e0*b             &
     &        +.238449e-2*b*c-.542819e-1*b*b

            ubw(i,j,k)=exp(d)*c5(i,j,k)*rrbw(ncp,3)

          else

            ubw(i,j,k)=ubwmax(i,j,k)

          end if

        end do
        end do

!$omp end do

      end do

! -----

!$omp end parallel

!! -----

      end subroutine s_termbw3d

!-----7--------------------------------------------------------------7--

      end module m_termbw3d
