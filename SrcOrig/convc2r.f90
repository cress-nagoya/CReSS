!***********************************************************************
      module m_convc2r
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/01
!     Modification: 2000/01/17, 2000/03/08, 2000/06/01, 2000/07/05,
!                   2001/10/18, 2001/11/20, 2002/04/02, 2003/04/30,
!                   2003/05/19, 2003/12/12, 2004/08/01, 2004/09/25,
!                   2004/12/17, 2006/04/03, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the auto conversion rate from the cloud water to the
!     rain water.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: convc2r, s_convc2r

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface convc2r

        module procedure s_convc2r

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
      subroutine s_convc2r(dtb,ni,nj,nk,qcp,qcf,qrf)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: qcp(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio at past

! Input and output variables

      real, intent(inout) :: qcf(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio at future

      real, intent(inout) :: qrf(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio at future

! Internal shared variable

      real ccncr       ! 0.001 x dtb

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real cncr        ! Auto conversion rate
                       ! from cloud water to rain water

!-----7--------------------------------------------------------------7--

! Set the common used variable.

      ccncr=.001e0*dtb

! -----

! Calculate the auto conversion rate from the cloud water to the rain
! water.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j,cncr)

        do j=1,nj-1
        do i=1,ni-1

          if(qcp(i,j,k).gt..001e0) then

            cncr=ccncr*(qcp(i,j,k)-.001e0)

            if(qcf(i,j,k).gt.cncr) then

              qcf(i,j,k)=qcf(i,j,k)-cncr
              qrf(i,j,k)=qrf(i,j,k)+cncr

            else

              qrf(i,j,k)=qrf(i,j,k)+qcf(i,j,k)
              qcf(i,j,k)=0.e0

            end if

          end if

        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_convc2r

!-----7--------------------------------------------------------------7--

      end module m_convc2r
