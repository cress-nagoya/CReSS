!***********************************************************************
      module m_steptund
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/16
!     Modification: 2001/12/11, 2002/01/15, 2002/04/02, 2002/07/03,
!                   2003/01/04, 2003/01/20, 2003/04/30, 2003/05/19,
!                   2003/07/15, 2003/10/31, 2003/12/12, 2004/04/15,
!                   2004/08/01, 2004/08/20, 2005/01/14, 2005/04/04,
!                   2006/04/03, 2007/07/30, 2007/10/19, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2009/02/27, 2009/08/20,
!                   2009/11/13, 2011/11/10, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     solve the soil and sea temperature to the next time step.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_comphy
      use m_gaussel
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: steptund, s_steptund

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface steptund

        module procedure s_steptund

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_steptund(fpsfcopt,fpdzgrd,fpdzsea,dtsoil,stinc,      &
     &                      ni,nj,nk,nund,t,land,cap,nuu,sst,sstd,      &
     &                      hs,le,rsd,rld,rlu,tundp,tundf,rr,ss,tt,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fpdzgrd
                       ! Formal parameter of unique index of dzgrd

      integer, intent(in) :: fpdzsea
                       ! Formal parameter of unique index of dzsea

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: dtsoil
                       ! Time interval of soil temperature calculation

      real, intent(in) :: stinc
                       ! Lapse of forecast time
                       ! from sea surface temperature data reading

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: cap(0:ni+1,0:nj+1)
                       ! Thermal capacity

      real, intent(in) :: nuu(0:ni+1,0:nj+1)
                       ! Thermal diffusivity

      real, intent(in) :: sst(0:ni+1,0:nj+1)
                       ! Sea surface temperature of external data
                       ! at marked time

      real, intent(in) :: sstd(0:ni+1,0:nj+1)
                       ! Time tendency of
                       ! sea surface temperature of external data

      real, intent(in) :: hs(0:ni+1,0:nj+1)
                       ! Sensible heat

      real, intent(in) :: le(0:ni+1,0:nj+1)
                       ! Latent heat

      real, intent(in) :: rsd(0:ni+1,0:nj+1)
                       ! Net downward short wave radiation

      real, intent(in) :: rld(0:ni+1,0:nj+1)
                       ! Downward long wave radiation

      real, intent(in) :: rlu(0:ni+1,0:nj+1)
                       ! Upward long wave radiation

! Input and output variable

      real, intent(inout) :: tundp(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at past

! Output variable

      real, intent(out) :: tundf(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at future

! Internal shared variables

      integer sfcopt   ! Option for surface physics

      integer nundm1   ! nund - 1

      real dzgrd       ! Grid distance in soil layers in z direction
      real dzsea       ! Grid distance in sea layers in z direction

      real ctg1        ! 1.0 / dzgrd x dtsoil
      real ctgm1       ! 1.0 / (dzgrd x dzgrd) x dtsoil

      real cts1        ! 1.0 / dzsea x dtsoil
      real ctsm1       ! 1.0 / (dzsea x dzsea) x dtsoil

      real s1g         ! 1.0 / (dzgrd x dzgrd) x dtsoil
      real skg         ! 2.0 / (dzgrd x dzgrd) x dtsoil
      real rkg         ! - 1.0 / (dzgrd x dzgrd) x dtsoil
      real tkg         ! - 1.0 / (dzgrd x dzgrd) x dtsoil

      real s1s         ! 1.0 / (dzsea x dzsea) x dtsoil
      real sks         ! 2.0 / (dzsea x dzsea) x dtsoil
      real rks         ! - 1.0 / (dzsea x dzsea) x dtsoil
      real tks         ! - 1.0 / (dzsea x dzsea) x dtsoil

      real, intent(inout) :: rr(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient matrix

      real, intent(inout) :: ss(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient matrix

      real, intent(inout) :: tt(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient matrix

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsfcopt,sfcopt)
      call getrname(fpdzgrd,dzgrd)
      call getrname(fpdzsea,dzsea)

! -----

! Set the common used variables.

      nundm1=nund-1

      ctg1=1.e0/dzgrd*dtsoil
      ctgm1=1.e0/(dzgrd*dzgrd)*dtsoil

      cts1=1.e0/dzsea*dtsoil
      ctsm1=1.e0/(dzsea*dzsea)*dtsoil

      s1g=1.e0/(dzgrd*dzgrd)*dtsoil
      skg=2.e0/(dzgrd*dzgrd)*dtsoil
      rkg=-1.e0/(dzgrd*dzgrd)*dtsoil
      tkg=-1.e0/(dzgrd*dzgrd)*dtsoil

      s1s=1.e0/(dzsea*dzsea)*dtsoil
      sks=2.e0/(dzsea*dzsea)*dtsoil
      rks=-1.e0/(dzsea*dzsea)*dtsoil
      tks=-1.e0/(dzsea*dzsea)*dtsoil

! -----

!!! Solve the soil and sea temperature to the next time step.

!! Set the top and bottom boundary conditions and coefficient matrix.

!$omp parallel default(shared) private(k)

! Set the ice and snow surface temperature.

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1

        if(land(i,j).ge.3.and.land(i,j).lt.10) then

          tundp(i,j,1)=min(.5e0*(t(i,j,1)+t(i,j,2)),t0)

        end if

      end do
      end do

!$omp end do

      do k=2,nund

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

          if(land(i,j).ge.3.and.land(i,j).lt.10) then

            tundp(i,j,k)=tundp(i,j,1)

          end if

        end do
        end do

!$omp end do

      end do

! -----

! Copy the past value to future and convert the unit from Kelvin to
! Celsius degrees.

      do k=1,nund

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          tundf(i,j,k)=tundp(i,j,k)-t0
        end do
        end do

!$omp end do

      end do

! -----

! Set the top and bottom boundary conditions.

      if(sfcopt.eq.1.or.sfcopt.eq.11) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

         if(land(i,j).lt.3) then

           tundf(i,j,1)=tundf(i,j,1)                                    &
     &       +cts1*(rsd(i,j)+rld(i,j)-rlu(i,j)-hs(i,j)-le(i,j))/cap(i,j)

           tundf(i,j,nundm1)=(1.e0+ctsm1*nuu(i,j))*tundf(i,j,nundm1)

         end if

         if(land(i,j).ge.10) then

           tundf(i,j,1)=tundf(i,j,1)                                    &
     &       +ctg1*(rsd(i,j)+rld(i,j)-rlu(i,j)-hs(i,j)-le(i,j))/cap(i,j)

           tundf(i,j,nundm1)                                            &
     &       =tundf(i,j,nundm1)+ctgm1*nuu(i,j)*tundf(i,j,nund)

         end if

        end do
        end do

!$omp end do

      else if(sfcopt.eq.2.or.sfcopt.eq.12) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

         if(land(i,j).ge.10) then

           tundf(i,j,1)=tundf(i,j,1)                                    &
     &       +ctg1*(rsd(i,j)+rld(i,j)-rlu(i,j)-hs(i,j)-le(i,j))/cap(i,j)

           tundf(i,j,nundm1)                                            &
     &       =tundf(i,j,nundm1)+ctgm1*nuu(i,j)*tundf(i,j,nund)

         end if

        end do
        end do

!$omp end do

      else if(sfcopt.eq.3.or.sfcopt.eq.13) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

         if(land(i,j).lt.3) then

           tundf(i,j,1)=(sst(i,j)+sstd(i,j)*stinc)-t0

         end if

         if(land(i,j).ge.10) then

           tundf(i,j,1)=tundf(i,j,1)                                    &
     &       +ctg1*(rsd(i,j)+rld(i,j)-rlu(i,j)-hs(i,j)-le(i,j))/cap(i,j)

           tundf(i,j,nundm1)                                            &
     &       =tundf(i,j,nundm1)+ctgm1*nuu(i,j)*tundf(i,j,nund)

         end if

        end do
        end do

!$omp end do

      end if

! -----

! Set the constant sea temperature.

      if(sfcopt.eq.3.or.sfcopt.eq.13) then

        do k=2,nund

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

            if(land(i,j).lt.3) then

              tundf(i,j,k)=tundf(i,j,1)

            end if

          end do
          end do

!$omp end do

        end do

      end if

! -----

! Set the coefficient matrix.

      if(sfcopt.eq.1.or.sfcopt.eq.11) then

        do k=1,nund-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

            if(land(i,j).lt.3) then

              rr(i,j,k)=rks*nuu(i,j)
              ss(i,j,k)=sks*nuu(i,j)+1.e0
              tt(i,j,k)=tks*nuu(i,j)

            else if(land(i,j).ge.10) then

              rr(i,j,k)=rkg*nuu(i,j)
              ss(i,j,k)=skg*nuu(i,j)+1.e0
              tt(i,j,k)=tkg*nuu(i,j)

            else

              rr(i,j,k)=0.e0
              ss(i,j,k)=1.e0
              tt(i,j,k)=0.e0

            end if

          end do
          end do

!$omp end do

        end do

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

          if(land(i,j).lt.3) then

            ss(i,j,1)=s1s*nuu(i,j)+1.e0

          else if(land(i,j).ge.10) then

            ss(i,j,1)=s1g*nuu(i,j)+1.e0

          else

            ss(i,j,1)=1.e0

          end if

        end do
        end do

!$omp end do

      else if(sfcopt.eq.2.or.sfcopt.eq.3.or.sfcopt.ge.12) then

        do k=1,nund-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

            if(land(i,j).ge.10) then

              rr(i,j,k)=rkg*nuu(i,j)
              ss(i,j,k)=skg*nuu(i,j)+1.e0
              tt(i,j,k)=tkg*nuu(i,j)

            else

              rr(i,j,k)=0.e0
              ss(i,j,k)=1.e0
              tt(i,j,k)=0.e0

            end if

          end do
          end do

!$omp end do

        end do

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

          if(land(i,j).ge.10) then

            ss(i,j,1)=s1g*nuu(i,j)+1.e0

          else

            ss(i,j,1)=1.e0

          end if

        end do
        end do

!$omp end do

      end if

! -----

!$omp end parallel

!! ------

! Solve the trifiagonal equation with the Gauss elimination.

      call gaussel(idoneopt,1,ni-1,1,nj-1,1,nund-1,ni,nj,nund,rr,ss,tt, &
     &             tundf,tmp1)

! -----

!! Set the bottom boundary conditions and convert the unit from Celsius
!! to Kelvin degrees.

!$omp parallel default(shared) private(k)

! Set the bottom boundary condition.

      if(sfcopt.eq.1.or.sfcopt.eq.11) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

          if(land(i,j).lt.3) then

            tundf(i,j,nund)=tundf(i,j,nundm1)

          end if

        end do
        end do

!$omp end do

      end if

! -----

! Convert the unit from Celsius to Kelvin degrees.

      do k=1,nund

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          tundf(i,j,k)=tundf(i,j,k)+t0
        end do
        end do

!$omp end do

      end do

! -----

!$omp end parallel

!! -----

!!! -----

      end subroutine s_steptund

!-----7--------------------------------------------------------------7--

      end module m_steptund
