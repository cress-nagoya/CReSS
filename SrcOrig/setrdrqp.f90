!***********************************************************************
      module m_setrdrqp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/07/30
!     Modification: 2007/09/25, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2009/03/23, 2011/09/22, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the interpolated radar variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setrdrqp, s_setrdrqp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setrdrqp

        module procedure s_setrdrqp

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
      subroutine s_setrdrqp(fpcphopt,fphaiopt,fprdritv,ird,             &
     &                      ni,nj,nk,nqw,nqi,qwrdr,qwrtd,qirdr,qirtd)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fprdritv
                       ! Formal parameter of unique index of rdritv

      integer, intent(in) :: ird
                       ! Index of count to read out in rdrdrnxt

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

      real, intent(inout) :: qwrdr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of radar data
                       ! at marked time

      real, intent(inout) :: qwrtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of
                       ! water hydrometeor of radar data

      real, intent(inout) :: qirdr(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of radar data
                       ! at marked time

      real, intent(inout) :: qirtd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of
                       ! ice hydrometeor of radar data

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes

      real rdritv      ! Time interval of radar data

      real rdriv       ! Inverse of rdritv

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getrname(fprdritv,rdritv)

! -----

! Set the common used variable.

      rdriv=1.e0/rdritv

! -----

!!! Set the interpolated radar variables.

!$omp parallel default(shared) private(k)

!! Set the time tendency of variables at current marked time.

      if(ird.eq.1) then

        if(abs(cphopt).lt.10) then

! Set the time tendency of rain water at current marked time.

          if(abs(cphopt).ge.1) then

            do k=1,nk

!$omp do schedule(runtime) private(i,j)

              do j=1,nj
              do i=1,ni

                if(qwrdr(i,j,k,2).gt.lim34n                             &
     &            .and.qwrtd(i,j,k,2).gt.lim34n) then

                  qwrtd(i,j,k,2)=(qwrtd(i,j,k,2)-qwrdr(i,j,k,2))*rdriv

                else

                  qwrtd(i,j,k,2)=lim35n

                end if

              end do
              end do

!$omp end do

            end do

          end if

! -----

! Set the time tendency of snow, graupel and hail at current marked
! time.

          if(abs(cphopt).ge.2) then

            if(haiopt.eq.0) then

              do k=1,nk

!$omp do schedule(runtime) private(i,j)

                do j=1,nj
                do i=1,ni

                  if(qirdr(i,j,k,2).gt.lim34n                           &
     &              .and.qirtd(i,j,k,2).gt.lim34n) then

                    qirtd(i,j,k,2)=(qirtd(i,j,k,2)-qirdr(i,j,k,2))*rdriv

                  else

                    qirtd(i,j,k,2)=lim35n

                  end if

                  if(qirdr(i,j,k,3).gt.lim34n                           &
     &              .and.qirtd(i,j,k,3).gt.lim34n) then

                    qirtd(i,j,k,3)=(qirtd(i,j,k,3)-qirdr(i,j,k,3))*rdriv

                  else

                    qirtd(i,j,k,3)=lim35n

                  end if

                end do
                end do

!$omp end do

              end do

            else

              do k=1,nk

!$omp do schedule(runtime) private(i,j)

                do j=1,nj
                do i=1,ni

                  if(qirdr(i,j,k,2).gt.lim34n                           &
     &              .and.qirtd(i,j,k,2).gt.lim34n) then

                    qirtd(i,j,k,2)=(qirtd(i,j,k,2)-qirdr(i,j,k,2))*rdriv

                  else

                    qirtd(i,j,k,2)=lim35n

                  end if

                  if(qirdr(i,j,k,3).gt.lim34n                           &
     &              .and.qirtd(i,j,k,3).gt.lim34n) then

                    qirtd(i,j,k,3)=(qirtd(i,j,k,3)-qirdr(i,j,k,3))*rdriv

                  else

                    qirtd(i,j,k,3)=lim35n

                  end if

                  if(qirdr(i,j,k,4).gt.lim34n                           &
     &              .and.qirtd(i,j,k,4).gt.lim34n) then

                    qirtd(i,j,k,4)=(qirtd(i,j,k,4)-qirdr(i,j,k,4))*rdriv

                  else

                    qirtd(i,j,k,4)=lim35n

                  end if

                end do
                end do

!$omp end do

              end do

            end if

          end if

! -----

        end if

      end if

!! -----

!! Set the variables at current marked time.

      if(ird.eq.2) then

        if(abs(cphopt).lt.10) then

! Set the rain water at current marked time.

          if(abs(cphopt).ge.1) then

            do k=1,nk

!$omp do schedule(runtime) private(i,j)

              do j=1,nj
              do i=1,ni

                if(qwrtd(i,j,k,2).gt.lim34n) then

                  qwrdr(i,j,k,2)=qwrtd(i,j,k,2)

                else

                  qwrdr(i,j,k,2)=lim35n

                end if

              end do
              end do

!$omp end do

            end do

          end if

! -----

! Set the snow, graupel and hail at current marked time.

          if(abs(cphopt).ge.2) then

            if(haiopt.eq.0) then

              do k=1,nk

!$omp do schedule(runtime) private(i,j)

                do j=1,nj
                do i=1,ni

                  if(qirtd(i,j,k,2).gt.lim34n) then

                    qirdr(i,j,k,2)=qirtd(i,j,k,2)

                  else

                    qirdr(i,j,k,2)=lim35n

                  end if

                  if(qirtd(i,j,k,3).gt.lim34n) then

                    qirdr(i,j,k,3)=qirtd(i,j,k,3)

                  else

                    qirdr(i,j,k,3)=lim35n

                  end if

                end do
                end do

!$omp end do

              end do

            else

              do k=1,nk

!$omp do schedule(runtime) private(i,j)

                do j=1,nj
                do i=1,ni

                  if(qirtd(i,j,k,2).gt.lim34n) then

                    qirdr(i,j,k,2)=qirtd(i,j,k,2)

                  else

                    qirdr(i,j,k,2)=lim35n

                  end if

                  if(qirtd(i,j,k,3).gt.lim34n) then

                    qirdr(i,j,k,3)=qirtd(i,j,k,3)

                  else

                    qirdr(i,j,k,3)=lim35n

                  end if

                  if(qirtd(i,j,k,4).gt.lim34n) then

                    qirdr(i,j,k,4)=qirtd(i,j,k,4)

                  else

                    qirdr(i,j,k,4)=lim35n

                  end if

                end do
                end do

!$omp end do

              end do

            end if

          end if

! -----

        end if

      end if

!! -----

!$omp end parallel

!!! -----

      end subroutine s_setrdrqp

!-----7--------------------------------------------------------------7--

      end module m_setrdrqp
