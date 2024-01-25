!***********************************************************************
      module m_fallqr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/01
!     Modification: 1999/11/19, 1999/11/24, 1999/12/15, 2000/01/17,
!                   2000/03/08, 2000/04/18, 2000/07/05, 2000/11/17,
!                   2001/06/29, 2001/12/11, 2002/01/15, 2002/04/02,
!                   2003/04/30, 2003/05/19, 2003/12/12, 2004/03/05,
!                   2004/03/22, 2004/04/15, 2004/06/10, 2004/09/01,
!                   2004/09/10, 2004/09/25, 2004/12/17, 2006/02/13,
!                   2006/04/03, 2006/05/12, 2006/09/30, 2007/05/14,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/11/05, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the fall out of the rain water.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkfall
      use m_comindx
      use m_commath
      use m_getrname
      use m_upwqp

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: fallqr, s_fallqr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface fallqr

        module procedure s_fallqr

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
      subroutine s_fallqr(fpdz,dtb,ni,nj,nk,jcb,rbr,rst,urq,qrf,prr,    &
     &                    tmp1)
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

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: urq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of rain water

! Input and output variables

      real, intent(inout) :: qrf(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio at future

      real, intent(inout) :: prr(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for rain

! Internal shared variables

      integer npstp    ! Number of steps for fall out time integration

      integer ipstp    ! Index of steps for fall out time integration

      real dz          ! Grid distance in z direction

      real dtp         ! Time steps interval of fall out integration

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

!! Perform the fall out.

! Initialize the processed variable, dtp.

      dtp=dtb

! -----

! Get the minimum time interval.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j) reduction(min: dtp)

        do j=1,nj-1
        do i=1,ni-1
          dtp=min(dz*jcb(i,j,k)/(urq(i,j,k)+eps),dtp)
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

! Calculate the sedimentation and precipitation of the rain water.

      do ipstp=1,npstp

        call upwqp(idadvopt,iddziv,dtp,ni,nj,nk,rbr,rst,urq,qrf,prr,    &
     &             tmp1)

      end do

! -----

!! -----

      end subroutine s_fallqr

!-----7--------------------------------------------------------------7--

      end module m_fallqr
