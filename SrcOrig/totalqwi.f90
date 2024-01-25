!***********************************************************************
      module m_totalqwi
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/01/10
!     Modification: 2006/09/30, 2007/05/14, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/01/30, 2009/02/27, 2011/08/18,
!                   2011/09/22, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the total water and ice mixing ratio.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: totalqwi, s_totalqwi

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface totalqwi

        module procedure s_totalqwi

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_totalqwi(fpcphopt,fphaiopt,                          &
     &                      ni,nj,nk,nqw,nqi,qwtr,qice,qall)
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

      real, intent(in) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor

      real, intent(in) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor

! Output variable

      real, intent(out) :: qall(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in bin categories

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)

! -----

!! Get the total water and ice mixing ratio.

!$omp parallel default(shared) private(k,n)

      if(abs(cphopt).ge.1) then

! For the bulk categories.

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).eq.1) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                qall(i,j,k)=qwtr(i,j,k,1)+qwtr(i,j,k,2)
              end do
              end do

!$omp end do

            end do

          else if(abs(cphopt).ge.2) then

            if(haiopt.eq.0) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  qall(i,j,k)=qwtr(i,j,k,1)+qwtr(i,j,k,2)               &
     &              +qice(i,j,k,1)+qice(i,j,k,2)+qice(i,j,k,3)
                end do
                end do

!$omp end do

              end do

            else

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  qall(i,j,k)=qwtr(i,j,k,1)+qwtr(i,j,k,2)+qice(i,j,k,1) &
     &              +qice(i,j,k,2)+qice(i,j,k,3)+qice(i,j,k,4)
                end do
                end do

!$omp end do

              end do

            end if

          end if

! -----

! For the bin categories.

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).eq.11) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                qall(i,j,k)=qwtr(i,j,k,1)
              end do
              end do

!$omp end do

            end do

            do n=2,nqw

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  qall(i,j,k)=qall(i,j,k)+qwtr(i,j,k,n)
                end do
                end do

!$omp end do

              end do

            end do

          else if(abs(cphopt).eq.12) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                qall(i,j,k)=qwtr(i,j,k,1)
              end do
              end do

!$omp end do

            end do

            do n=2,nqw

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  qall(i,j,k)=qall(i,j,k)+qwtr(i,j,k,n)
                end do
                end do

!$omp end do

              end do

            end do

            do n=2,nqi

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1
                  qall(i,j,k)=qall(i,j,k)+qice(i,j,k,n)
                end do
                end do

!$omp end do

              end do

            end do

          end if

        end if

! -----

      end if

!$omp end parallel

!! -----

      end subroutine s_totalqwi

!-----7--------------------------------------------------------------7--

      end module m_totalqwi
