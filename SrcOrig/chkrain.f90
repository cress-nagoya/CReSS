!***********************************************************************
      module m_chkrain
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/08/01
!     Modification: 2004/09/01, 2006/01/10, 2006/04/03, 2007/01/20,
!                   2007/05/07, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/01/30, 2009/02/27, 2010/02/01, 2011/09/22,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the precipitation on the surface.

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

      public :: chkrain, s_chkrain

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkrain

        module procedure s_chkrain

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
      subroutine s_chkrain(fpcphopt,fphaiopt,fmois,ni,nj,nqw,nqi,       &
     &                     prwtr,price,fall)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      real, intent(in) :: prwtr(0:ni+1,0:nj+1,1:2,1:nqw)
                       ! Precipitation and accumulation for water

      real, intent(in) :: price(0:ni+1,0:nj+1,1:2,1:nqi)
                       ! Precipitation and accumulation for ice

! Output variable

      real, intent(out) :: fall(0:ni+1,0:nj+1)
                       ! Precipitation flag

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)

! -----

!!!! Check the precipitation on the surface.

!$omp parallel default(shared)

! Fill in the undefined value in the case of dry run.

      if(fmois(1:3).eq.'dry') then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          fall(i,j)=-1.e0
        end do
        end do

!$omp end do

! -----

!!! Check the precipitation in the case of moist run.

      else if(fmois(1:5).eq.'moist') then

! Fill in the undefined value in the case of no cloud physics.

        if(abs(cphopt).eq.0) then

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            fall(i,j)=-1.e0
          end do
          end do

!$omp end do

        else

! -----

!! Check the precipitation in the case of processing cloud physics.

! For the bulk method.

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).eq.1) then

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1

                if(prwtr(i,j,1,1).gt.prmin.or.                          &
     &             prwtr(i,j,1,2).gt.prmin) then

                  fall(i,j)=1.e0

                else

                  fall(i,j)=-1.e0

                end if

              end do
              end do

!$omp end do

            else if(abs(cphopt).ge.2) then

              if(haiopt.eq.0) then

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1

                  if(prwtr(i,j,1,1).gt.prmin.or.                        &
     &               prwtr(i,j,1,2).gt.prmin.or.                        &
     &               price(i,j,1,1).gt.prmin.or.                        &
     &               price(i,j,1,2).gt.prmin.or.                        &
     &               price(i,j,1,3).gt.prmin) then

                    fall(i,j)=1.e0

                  else

                    fall(i,j)=-1.e0

                  end if

                end do
                end do

!$omp end do

              else

!$omp do schedule(runtime) private(i,j)

                do j=1,nj-1
                do i=1,ni-1

                  if(prwtr(i,j,1,1).gt.prmin.or.                        &
     &               prwtr(i,j,1,2).gt.prmin.or.                        &
     &               price(i,j,1,1).gt.prmin.or.                        &
     &               price(i,j,1,2).gt.prmin.or.                        &
     &               price(i,j,1,3).gt.prmin.or.                        &
     &               price(i,j,1,4).gt.prmin) then

                    fall(i,j)=1.e0

                  else

                    fall(i,j)=-1.e0

                  end if

                end do
                end do

!$omp end do

              end if

            end if

! -----

! For the bin method.

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).eq.11) then

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1

                if(prwtr(i,j,1,1).gt.prmin) then
                  fall(i,j)=1.e0
                else
                  fall(i,j)=-1.e0
                end if

              end do
              end do

!$omp end do

            else if(abs(cphopt).eq.12) then

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1

                if(prwtr(i,j,1,1).gt.prmin.or.                          &
     &             price(i,j,1,1).gt.prmin) then

                  fall(i,j)=1.e0

                else

                  fall(i,j)=-1.e0

                end if

              end do
              end do

!$omp end do

            end if

          end if

! -----

!! -----

        end if

      end if

!!! -----

!$omp end parallel

!!!! -----

      end subroutine s_chkrain

!-----7--------------------------------------------------------------7--

      end module m_chkrain
