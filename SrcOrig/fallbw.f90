!***********************************************************************
      module m_fallbw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/09/30
!     Modification: 2007/05/14, 2007/10/19, 2008/01/11, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2009/11/05, 2011/08/18,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the fall out.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkfall
      use m_comindx
      use m_commath
      use m_getrname
      use m_termbw3d
      use m_upwnbin
      use m_upwmbin

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: fallbw, s_fallbw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface fallbw

        module procedure s_fallbw

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic int
      intrinsic min
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_fallbw(fpdz,dtb,ni,nj,nk,nqw,nnw,rbw,rrbw,           &
     &                    jcb,rbr,rst,rbv,t,mwbin,nwbin,prr,            &
     &                    ubw,ubwmax,c1,c2,c3,c4,c5,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentrations

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: rbw(1:nqw)
                       ! Standard radius
                       ! between adjacent water bins [cm]

      real, intent(in) :: rrbw(1:nqw,1:5)
                       ! Related parameters of rbw

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density [g/cm^3]

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian [g/cm^3]

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density [cm^3/g]

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

! Input and output variables

      real, intent(inout) :: mwbin(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Total water mass [g/cm^3]

      real, intent(inout) :: nwbin(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations [1/cm^3]

      real, intent(inout) :: prr(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation
                       ! for rain [cm/s], [cm]

! Internal shared variables

      integer n        ! Array index in water bin categories

      integer npstp    ! Number of steps for fall out time integration

      integer ipstp    ! Index of steps for fall out time integration

      real dz          ! Grid distance in z direction

      real dz2         ! 100.0 x dz

      real dtp         ! Time steps interval of fall out integration

      real, intent(inout) :: ubw(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of water bin

      real, intent(inout) :: ubwmax(0:ni+1,0:nj+1,1:nk)
                       ! Maximum terminal velocity of water bin

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

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpdz,dz)

! -----

! Set the common used variable.

      dz2=100.e0*dz

! -----

!! Perform the fall out.

      do n=1,nqw

! Calculate the terminal velocity for water bin.

        call termbw3d(n,ni,nj,nk,nqw,                                   &
     &                rbw,rrbw,rbr,rbv,t,ubwmax,c1,c2,c3,c4,c5,ubw)

! -----

! Initialize the processed variable, dtp.

        dtp=dtb

! -----

! Get the minimum time interval.

!$omp parallel default(shared) private(k)

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j) reduction(min: dtp)

          do j=1,nj-1
          do i=1,ni-1
            dtp=min(dz2*jcb(i,j,k)/(ubw(i,j,k)+eps),dtp)
          end do
          end do

!$omp end do

        end do

!$omp end parallel

! -----

! Get the time interval and number of steps.

        call chkfall(dtp)

        if(dtp.lt.dtb) then
          npstp=int((dtb+.001e0)/dtp)+1
          dtp=dtb/real(npstp)

        else
          npstp=1
          dtp=dtb

        end if

! -----

! Calculate the sedimentation and precipitation.

        do ipstp=1,npstp

          call upwmbin(idadvopt,iddziv,n,dtp,ni,nj,nk,nqw,rbr,rst,rbv,  &
     &                 ubw,mwbin,prr,tmp1)

          call upwnbin(iddziv,n,dtp,ni,nj,nk,nnw,rbr,rst,ubw,nwbin,tmp1)

        end do

! -----

      end do

!! -----

      end subroutine s_fallbw

!-----7--------------------------------------------------------------7--

      end module m_fallbw
