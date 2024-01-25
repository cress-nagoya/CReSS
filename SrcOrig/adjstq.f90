!***********************************************************************
      module m_adjstq
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/01
!     Modification: 2000/01/17, 2000/04/18, 2000/07/05, 2002/01/15,
!                   2002/04/02, 2003/04/30, 2003/05/19, 2004/06/10,
!                   2004/09/01, 2004/09/25, 2004/10/12, 2006/01/10,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/11/05, 2011/08/18, 2011/09/22, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     force the hydrometeor mixing ratio more than user specified value.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: adjstq, s_adjstq

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface adjstq

        module procedure s_adjstq

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_adjstq(fpcphopt,fphaiopt,ni,nj,nk,nqw,nqi,           &
     &                    qv,qwtr,qice)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

! Input and output variables

      real, intent(inout) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(inout) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor

      real, intent(inout) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)

! -----

! Force the hydrometeor mixing ratio more than user specified value.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          qv(i,j,k)=max(qv(i,j,k),0.e0)
        end do
        end do

!$omp end do

      end do

      if(abs(cphopt).lt.10) then

        if(abs(cphopt).ge.1) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              qwtr(i,j,k,1)=max(qwtr(i,j,k,1),0.e0)
              qwtr(i,j,k,2)=max(qwtr(i,j,k,2),0.e0)
            end do
            end do

!$omp end do

          end do

        end if

        if(abs(cphopt).ge.2) then

          if(haiopt.eq.0) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                qice(i,j,k,1)=max(qice(i,j,k,1),0.e0)
                qice(i,j,k,2)=max(qice(i,j,k,2),0.e0)
                qice(i,j,k,3)=max(qice(i,j,k,3),0.e0)
              end do
              end do

!$omp end do

            end do

          else

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                qice(i,j,k,1)=max(qice(i,j,k,1),0.e0)
                qice(i,j,k,2)=max(qice(i,j,k,2),0.e0)
                qice(i,j,k,3)=max(qice(i,j,k,3),0.e0)
                qice(i,j,k,4)=max(qice(i,j,k,4),0.e0)
              end do
              end do

!$omp end do

            end do

          end if

        end if

      end if

!$omp end parallel

! -----

      end subroutine s_adjstq

!-----7--------------------------------------------------------------7--

      end module m_adjstq
