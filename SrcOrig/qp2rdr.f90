!***********************************************************************
      module m_qp2rdr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2002/10/31, 2003/04/30, 2003/05/19, 2003/11/28,
!                   2003/12/12, 2007/05/07, 2007/07/30, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2009/03/23,
!                   2009/11/13, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the analysis nudging to radar data of optional
!     precipitation mixing ratio.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: qp2rdr, s_qp2rdr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface qp2rdr

        module procedure s_qp2rdr

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_qp2rdr(fpngropt,ngrdmp,rtinc,ni,nj,nk,rst,qpp,       &
     &                    qprdr,qprtd,qpfrc)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpngropt
                       ! Formal parameter of unique index of ngropt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: ngrdmp(1:2)
                       ! Analysis nudging damping coefficient

      real, intent(in) :: rtinc(1:2)
                       ! Lapse of forecast time from radar data reading

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jabobian

      real, intent(in) :: qpp(0:ni+1,0:nj+1,1:nk)
                       ! Optional precipitation mixing ratio at past

      real, intent(in) :: qprdr(0:ni+1,0:nj+1,1:nk)
                       ! Optional precipitation mixing ratio
                       ! at marked time

      real, intent(in) :: qprtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! optional precipitation mixing ratio

! Input and output variable

      real, intent(inout) :: qpfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! of optional precipitation mixing ratio

! Internal shared variable

      integer ngropt   ! Option for analysis nudging to radar

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real cqprdr      ! Optional precipitation mixing ratio
                       ! at current forecast time

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpngropt,ngropt)

! -----

! Calculate the analysis nudging terms for optional precipitation mixing
! ratio.

!$omp parallel default(shared) private(k)

      if(ngropt.eq.1.and.ngrdmp(1).gt.0.e0) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j,cqprdr)

          do j=2,nj-2
          do i=2,ni-2

            if(qprdr(i,j,k).gt.lim34n.and.qprtd(i,j,k).gt.lim34n) then

              cqprdr=max(qprdr(i,j,k)+qprtd(i,j,k)*rtinc(1),0.e0)

              qpfrc(i,j,k)=qpfrc(i,j,k)                                 &
     &          +ngrdmp(1)*rst(i,j,k)*(cqprdr-qpp(i,j,k))

            end if

          end do
          end do

!$omp end do

        end do

      else if(ngropt.ge.2.and.ngrdmp(2).gt.0.e0) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2

            if(qprdr(i,j,k).gt.lim34n) then

              qpfrc(i,j,k)=qpfrc(i,j,k)                                 &
     &          +ngrdmp(2)*rst(i,j,k)*(qprdr(i,j,k)-qpp(i,j,k))

            end if

          end do
          end do

!$omp end do

        end do

      end if

!$omp end parallel

! -----

      end subroutine s_qp2rdr

!-----7--------------------------------------------------------------7--

      end module m_qp2rdr
