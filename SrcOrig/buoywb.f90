!***********************************************************************
      module m_buoywb
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/05,
!                   1999/08/23, 1999/09/30, 1999/10/12, 1999/11/01,
!                   1999/11/19, 2000/01/17, 2000/04/12, 2000/04/18,
!                   2000/07/05, 2001/01/15, 2001/08/07, 2001/11/20,
!                   2001/12/11, 2002/04/02, 2002/08/15, 2002/12/02,
!                   2003/04/30, 2003/05/19, 2003/11/05, 2003/11/28,
!                   2003/12/12, 2004/09/10, 2006/01/10, 2006/11/06,
!                   2007/01/20, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the buoyancy in the large time steps integration.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: buoywb, s_buoywb

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface buoywb

        module procedure s_buoywb

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
      subroutine s_buoywb(fpgwmopt,fpcphopt,fmois,ni,nj,nk,ptbr,qvbr,   &
     &                    rst,ptp,qv,qall,wfrc,qvd,wb8s)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpgwmopt
                       ! Formal parameter of unique index of gwmopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: qall(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio

! Input and output variable

      real, intent(inout) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

! Internal shared variables

      integer gwmopt   ! Option for gravity wave mode integration
      integer cphopt   ! Option for cloud micro physics

      real g05         ! 0.5 x g

      real, intent(inout) :: qvd(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor perturbation

      real, intent(inout) :: wb8s(0:ni+1,0:nj+1,1:nk)
                       ! 0.5 x buoyacy value
                       ! in large time steps at scalar points

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpgwmopt,gwmopt)
      call getiname(fpcphopt,cphopt)

! -----

! Set the common used variable.

      g05=.5e0*g

! -----

!! Calculate the buoyancy in the large time steps integration.

!$omp parallel default(shared) private(k)

! For dry air case.

      if(fmois(1:3).eq.'dry') then

        if(gwmopt.eq.0) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              wb8s(i,j,k)=g05*rst(i,j,k)*ptp(i,j,k)/ptbr(i,j,k)
            end do
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              wb8s(i,j,k)=0.e0
            end do
            end do

!$omp end do

          end do

        end if

! -----

! For moist air case.

      else if(fmois(1:5).eq.'moist') then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            qvd(i,j,k)=qv(i,j,k)-qvbr(i,j,k)
          end do
          end do

!$omp end do

        end do

        if(abs(cphopt).eq.0) then

          if(gwmopt.eq.0) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                wb8s(i,j,k)=g05*rst(i,j,k)                              &
     &            *(qvd(i,j,k)/(epsva+qvbr(i,j,k))                      &
     &            -qvd(i,j,k)/(1.e0+qvbr(i,j,k))+ptp(i,j,k)/ptbr(i,j,k))
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                wb8s(i,j,k)=g05*rst(i,j,k)                              &
     &            *(qvd(i,j,k)/(epsva+qvbr(i,j,k))                      &
     &            -qvd(i,j,k)/(1.e0+qvbr(i,j,k)))
              end do
              end do

!$omp end do

            end do

          end if

        else

          if(gwmopt.eq.0) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                wb8s(i,j,k)=g05*rst(i,j,k)                              &
     &            *(qvd(i,j,k)/(epsva+qvbr(i,j,k))                      &
     &            -(qvd(i,j,k)+qall(i,j,k))/(1.e0+qvbr(i,j,k))          &
     &            +ptp(i,j,k)/ptbr(i,j,k))
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                wb8s(i,j,k)=g05*rst(i,j,k)                              &
     &            *(qvd(i,j,k)/(epsva+qvbr(i,j,k))                      &
     &            -(qvd(i,j,k)+qall(i,j,k))/(1.e0+qvbr(i,j,k)))
              end do
              end do

!$omp end do

            end do

          end if

        end if

      end if

! -----

! Finally be averaged vartically.

      do k=3,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          wfrc(i,j,k)=wfrc(i,j,k)+(wb8s(i,j,k-1)+wb8s(i,j,k))
        end do
        end do

!$omp end do

      end do

! -----

!$omp end parallel

!! -----

      end subroutine s_buoywb

!-----7--------------------------------------------------------------7--

      end module m_buoywb
