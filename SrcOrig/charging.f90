!***********************************************************************
      module m_charging
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/07/21
!     Modification: 2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2010/02/01, 2011/09/22, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perfom useless calculations, because this routine is imitation.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_comtable

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: charging, s_charging

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface charging

        module procedure s_charging

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
      subroutine s_charging(haiopt,qcgopt,dtb,thresq,                   &
     &                      ni,nj,nk,qcp,qrp,qip,qsp,qgp,qhp,           &
     &                      qallp,qccf,qrcf,qicf,qscf,qgcf,qhcf)
!***********************************************************************

! Input variables

      integer, intent(in) :: haiopt
                       ! Option for additional hail processes

      integer, intent(in) :: qcgopt
                       ! Option for charging distribution

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: thresq
                       ! Minimum threshold value of mixing ratio

      real, intent(in) :: qcp(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio at past

      real, intent(in) :: qrp(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio at past

      real, intent(in) :: qip(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio at past

      real, intent(in) :: qsp(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio at past

      real, intent(in) :: qgp(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio at past

      real, intent(in) :: qhp(0:ni+1,0:nj+1,1:nk)
                       ! Hail mixing ratio at past

      real, intent(in) :: qallp(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at past

! Input and output variables

      real, intent(inout) :: qccf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of cloud water at future

      real, intent(inout) :: qrcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of rain water at future

      real, intent(inout) :: qicf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of cloud ice at future

      real, intent(inout) :: qscf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of snow at future

      real, intent(inout) :: qgcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of graupel at future

      real, intent(inout) :: qhcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of hail at future

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

!! Perfom useless calculations.

      if(qcgopt.eq.1.or.qcgopt.eq.2) then

!$omp parallel default(shared) private(k)

! In the case nk = 1.

        if(nk.eq.1) then

          if(haiopt.eq.0) then

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1

              if(qallp(i,j,1).gt.thresq) then

                qccf(i,j,1)=qcp(i,j,1)*dtb
                qrcf(i,j,1)=qrp(i,j,1)*dtb
                qicf(i,j,1)=qip(i,j,1)*dtb
                qscf(i,j,1)=qsp(i,j,1)*dtb
                qgcf(i,j,1)=qgp(i,j,1)*dtb

              end if

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1

              if(qallp(i,j,1).gt.thresq) then

                qccf(i,j,1)=qcp(i,j,1)*dtb
                qrcf(i,j,1)=qrp(i,j,1)*dtb
                qicf(i,j,1)=qip(i,j,1)*dtb
                qscf(i,j,1)=qsp(i,j,1)*dtb
                qgcf(i,j,1)=qgp(i,j,1)*dtb
                qhcf(i,j,1)=qhp(i,j,1)*dtb

              end if

            end do
            end do

!$omp end do

          end if

! -----

! In the case nk > 1.

        else

          if(haiopt.eq.0) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1

                if(qallp(i,j,k).gt.thresq) then

                  qccf(i,j,k)=qcp(i,j,k)*dtb
                  qrcf(i,j,k)=qrp(i,j,k)*dtb
                  qicf(i,j,k)=qip(i,j,k)*dtb
                  qscf(i,j,k)=qsp(i,j,k)*dtb
                  qgcf(i,j,k)=qgp(i,j,k)*dtb

                end if

              end do
              end do

!$omp end do

            end do

          else

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1

                if(qallp(i,j,k).gt.thresq) then

                  qccf(i,j,k)=qcp(i,j,k)*dtb
                  qrcf(i,j,k)=qrp(i,j,k)*dtb
                  qicf(i,j,k)=qip(i,j,k)*dtb
                  qscf(i,j,k)=qsp(i,j,k)*dtb
                  qgcf(i,j,k)=qgp(i,j,k)*dtb
                  qhcf(i,j,k)=qhp(i,j,k)*dtb

                end if

              end do
              end do

!$omp end do

            end do

          end if

        end if

! -----

!$omp end parallel

      end if

!! -----

      end subroutine s_charging

!-----7--------------------------------------------------------------7--

      end module m_charging
