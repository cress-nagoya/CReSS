!***********************************************************************
      module m_termqr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/01
!     Modification: 1999/11/19, 1999/11/24, 1999/12/15, 2000/01/17,
!                   2000/03/08, 2000/04/18, 2000/07/05, 2000/11/17,
!                   2001/06/29, 2001/12/11, 2002/01/15, 2002/04/02,
!                   2003/04/30, 2003/05/19, 2003/12/12, 2004/03/05,
!                   2004/03/22, 2004/04/15, 2004/06/10, 2004/09/01,
!                   2004/09/10, 2004/09/25, 2004/12/17, 2006/02/13,
!                   2006/04/03, 2006/09/30, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the terminal velocity of the rain water.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: termqr, s_termqr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface termqr

        module procedure s_termqr

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_termqr(fpthresq,ni,nj,nk,rbr,rbv,qr,urq)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpthresq
                       ! Formal parameter of unique index of thresq

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(in) :: qr(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio

! Output variable

      real, intent(out) :: urq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of rain water

! Internal shared variable

      real thresq      ! Minimum threshold value of mixing ratio

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getrname(fpthresq,thresq)

! -----

! Calculate the terminal velocity of the rain water.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

          if(qr(i,j,k).gt.thresq) then

            urq(i,j,k)=14.1640e0*sqrt(r0*rbv(i,j,k))                    &
     &        *exp(.1364e0*log(rbr(i,j,k)*qr(i,j,k)))

          else

            urq(i,j,k)=0.e0

          end if

        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_termqr

!-----7--------------------------------------------------------------7--

      end module m_termqr
