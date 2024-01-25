!***********************************************************************
      module m_exbcpt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/06/07
!     Modification: 1999/07/05, 1999/08/03, 1999/08/09, 1999/09/30,
!                   1999/11/01, 2000/01/17, 2000/02/02, 2000/04/18,
!                   2001/01/15, 2001/03/13, 2001/04/15, 2001/05/29,
!                   2001/06/06, 2001/07/13, 2001/08/07, 2001/12/11,
!                   2002/04/02, 2002/07/23, 2002/08/15, 2002/10/31,
!                   2003/03/28, 2003/04/30, 2003/05/19, 2003/06/27,
!                   2003/11/05, 2003/11/28, 2003/12/12, 2004/04/15,
!                   2004/05/07, 2004/08/20, 2005/01/07, 2005/02/10,
!                   2006/04/03, 2006/09/21, 2006/12/04, 2007/01/05,
!                   2007/01/31, 2007/05/07, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2008/12/11, 2009/02/27, 2009/03/23,
!                   2011/09/22, 2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     force the lateral boundary value to the external boundary value
!     for potential temperature.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: exbcpt, s_exbcpt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface exbcpt

        module procedure s_exbcpt

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_exbcpt(fpexbvar,fpwbc,fpebc,fpadvopt,fpexnews,       &
     &                    ivstp,dt,gtinc,ni,nj,nk,ptp,ptpp,             &
     &                    ptcpx,ptcpy,ptpgpv,ptptd,ptpf)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexbvar
                       ! Formal parameter of unique index of exbvar

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpexnews
                       ! Formal parameter of unique index of exnews

      integer, intent(in) :: ivstp
                       ! Index of time steps
                       ! of vertical Cubic Lagrange advection

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dt
                       ! Time steps interval

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: ptcpx(1:nj,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on west and east boundary

      real, intent(in) :: ptcpy(1:ni,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on south and north boundary

      real, intent(in) :: ptpgpv(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation of GPV data
                       ! at marked time

      real, intent(in) :: ptptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! potential temperature perturbation of GPV data

! Input and output variable

      real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

! Internal shared variables

      character(len=108) exbvar
                       ! Control flag of
                       ! extrenal boundary forced variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer advopt   ! Option for advection scheme

      integer nim1     ! ni - 1
      integer nim2     ! ni - 2
      integer njm1     ! nj - 1
      integer njm2     ! nj - 2

      real exnews      ! Boundary damping coefficient

      real gtinc1      ! gtinc
      real gtinc2      ! gtinc + dt

      real dt2         ! 2.0 x dt

      real dmpdt       ! exnews x dt or 2.0 x exnews x dt

      real tpdt        ! gtinc + real(ivstp - 1) x dt

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real ptb1        ! Temporary variable
      real ptb2        ! Temporary variable
      real ptb2i       ! Temporary variable
      real ptb2j       ! Temporary variable

      real gamma       ! Temporary variable

      real radwe       ! Temporary variable
      real radsn       ! Temporary variable

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(exbvar)

! -----

! Get the required namelist variables.

      call getcname(fpexbvar,exbvar)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpadvopt,advopt)
      call getrname(fpexnews,exnews)

! -----

! Set the common used variables.

      nim1=ni-1
      nim2=ni-2
      njm1=nj-1
      njm2=nj-2

      gtinc1=gtinc
      gtinc2=gtinc+dt

      dt2=2.e0*dt

      if(advopt.le.3) then
        dmpdt=exnews*dt2
      else
        dmpdt=exnews*dt
      end if

      tpdt=gtinc+real(ivstp-1)*dt

! -----

!! Force the lateral boundary value to the external boundary value.

!$omp parallel default(shared)

! Force the boundary value to the external boundary value at the four
! corners.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(exbvar(5:5).eq.'-') then

          if(abs(wbc).ne.1.and.abs(ebc).ne.1) then

            if(advopt.le.3) then

              if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime) private(k,ptb1,ptb2i,ptb2j,radwe,radsn)

                do k=2,nk-2
                  ptb1=ptpgpv(1,1,k)+ptptd(1,1,k)*gtinc1
                  ptb2i=ptpgpv(2,1,k)+ptptd(2,1,k)*gtinc2
                  ptb2j=ptpgpv(1,2,k)+ptptd(1,2,k)*gtinc2

                  radwe=((ptp(2,1,k)-ptpp(1,1,k))-(ptb2i-ptb1))         &
     &              *ptcpx(1,k,1)/(1.e0-ptcpx(1,k,1))

                  radsn=((ptp(1,2,k)-ptpp(1,1,k))-(ptb2j-ptb1))         &
     &              *ptcpy(1,k,1)/(1.e0-ptcpy(1,k,1))

                  ptpf(1,1,k)=ptpp(1,1,k)+ptptd(1,1,k)*dt2              &
     &              -2.e0*(radwe+radsn)-dmpdt*(ptpp(1,1,k)-ptb1)

                end do

!$omp end do

              end if

              if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime) private(k,ptb1,ptb2i,ptb2j,radwe,radsn)

                do k=2,nk-2
                  ptb1=ptpgpv(nim1,1,k)+ptptd(nim1,1,k)*gtinc1
                  ptb2i=ptpgpv(nim2,1,k)+ptptd(nim2,1,k)*gtinc2
                  ptb2j=ptpgpv(nim1,2,k)+ptptd(nim1,2,k)*gtinc2

                  radwe=((ptp(nim2,1,k)-ptpp(nim1,1,k))-(ptb2i-ptb1))   &
     &              *ptcpx(1,k,2)/(1.e0+ptcpx(1,k,2))

                  radsn=((ptp(nim1,2,k)-ptpp(nim1,1,k))-(ptb2j-ptb1))   &
     &              *ptcpy(nim1,k,1)/(1.e0-ptcpy(nim1,k,1))

                  ptpf(nim1,1,k)=ptpp(nim1,1,k)+ptptd(nim1,1,k)*dt2     &
     &              +2.e0*(radwe-radsn)-dmpdt*(ptpp(nim1,1,k)-ptb1)

                end do

!$omp end do

              end if

            else

              if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime) private(k,ptb1,ptb2i,ptb2j,radwe,radsn)

                do k=2,nk-2
                  ptb1=ptpgpv(1,1,k)+ptptd(1,1,k)*tpdt
                  ptb2i=ptpgpv(2,1,k)+ptptd(2,1,k)*tpdt
                  ptb2j=ptpgpv(1,2,k)+ptptd(1,2,k)*tpdt

                  radwe=((ptpp(2,1,k)-ptpp(1,1,k))-(ptb2i-ptb1))        &
     &              *ptcpx(1,k,1)

                  radsn=((ptpp(1,2,k)-ptpp(1,1,k))-(ptb2j-ptb1))        &
     &              *ptcpy(1,k,1)

                  ptpf(1,1,k)=ptpp(1,1,k)+ptptd(1,1,k)*dt               &
     &              -(radwe+radsn)-dmpdt*(ptpp(1,1,k)-ptb1)

                end do

!$omp end do

              end if

              if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime) private(k,ptb1,ptb2i,ptb2j,radwe,radsn)

                do k=2,nk-2
                  ptb1=ptpgpv(nim1,1,k)+ptptd(nim1,1,k)*tpdt
                  ptb2i=ptpgpv(nim2,1,k)+ptptd(nim2,1,k)*tpdt
                  ptb2j=ptpgpv(nim1,2,k)+ptptd(nim1,2,k)*tpdt

                  radwe=((ptpp(nim2,1,k)-ptpp(nim1,1,k))-(ptb2i-ptb1))  &
     &              *ptcpx(1,k,2)

                  radsn=((ptpp(nim1,2,k)-ptpp(nim1,1,k))-(ptb2j-ptb1))  &
     &              *ptcpy(nim1,k,1)

                  ptpf(nim1,1,k)=ptpp(nim1,1,k)+ptptd(nim1,1,k)*dt      &
     &              +(radwe-radsn)-dmpdt*(ptpp(nim1,1,k)-ptb1)

                end do

!$omp end do

              end if

            end if

          end if

        end if

      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(exbvar(5:5).eq.'-') then

          if(abs(wbc).ne.1.and.abs(ebc).ne.1) then

            if(advopt.le.3) then

              if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime) private(k,ptb1,ptb2i,ptb2j,radwe,radsn)

                do k=2,nk-2
                  ptb1=ptpgpv(1,njm1,k)+ptptd(1,njm1,k)*gtinc1
                  ptb2i=ptpgpv(2,njm1,k)+ptptd(2,njm1,k)*gtinc2
                  ptb2j=ptpgpv(1,njm2,k)+ptptd(1,njm2,k)*gtinc2

                  radwe=((ptp(2,njm1,k)-ptpp(1,njm1,k))-(ptb2i-ptb1))   &
     &              *ptcpx(njm1,k,1)/(1.e0-ptcpx(njm1,k,1))

                  radsn=((ptp(1,njm2,k)-ptpp(1,njm1,k))-(ptb2j-ptb1))   &
     &              *ptcpy(1,k,2)/(1.e0+ptcpy(1,k,2))

                  ptpf(1,njm1,k)=ptpp(1,njm1,k)+ptptd(1,njm1,k)*dt2     &
     &              -2.e0*(radwe-radsn)-dmpdt*(ptpp(1,njm1,k)-ptb1)

                end do

!$omp end do

              end if

              if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime) private(k,ptb1,ptb2i,ptb2j,radwe,radsn)

                do k=2,nk-2
                 ptb1=ptpgpv(nim1,njm1,k)+ptptd(nim1,njm1,k)*gtinc1
                 ptb2i=ptpgpv(nim2,njm1,k)+ptptd(nim2,njm1,k)*gtinc2
                 ptb2j=ptpgpv(nim1,njm2,k)+ptptd(nim1,njm2,k)*gtinc2

                 radwe=((ptp(nim2,njm1,k)-ptpp(nim1,njm1,k))            &
     &             -(ptb2i-ptb1))*ptcpx(njm1,k,2)/(1.e0+ptcpx(njm1,k,2))

                 radsn=((ptp(nim1,njm2,k)-ptpp(nim1,njm1,k))            &
     &             -(ptb2j-ptb1))*ptcpy(nim1,k,2)/(1.e0+ptcpy(nim1,k,2))

                 ptpf(nim1,njm1,k)=ptpp(nim1,njm1,k)                    &
     &             +ptptd(nim1,njm1,k)*dt2+2.e0*(radwe+radsn)           &
     &             -dmpdt*(ptpp(nim1,njm1,k)-ptb1)

                end do

!$omp end do

              end if

            else

              if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime) private(k,ptb1,ptb2i,ptb2j,radwe,radsn)

                do k=2,nk-2
                  ptb1=ptpgpv(1,njm1,k)+ptptd(1,njm1,k)*tpdt
                  ptb2i=ptpgpv(2,njm1,k)+ptptd(2,njm1,k)*tpdt
                  ptb2j=ptpgpv(1,njm2,k)+ptptd(1,njm2,k)*tpdt

                  radwe=((ptpp(2,njm1,k)-ptpp(1,njm1,k))-(ptb2i-ptb1))  &
     &              *ptcpx(njm1,k,1)

                  radsn=((ptpp(1,njm2,k)-ptpp(1,njm1,k))-(ptb2j-ptb1))  &
     &              *ptcpy(1,k,2)

                  ptpf(1,njm1,k)=ptpp(1,njm1,k)+ptptd(1,njm1,k)*dt      &
     &              -(radwe-radsn)-dmpdt*(ptpp(1,njm1,k)-ptb1)

                end do

!$omp end do

              end if

              if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime) private(k,ptb1,ptb2i,ptb2j,radwe,radsn)

                do k=2,nk-2
                  ptb1=ptpgpv(nim1,njm1,k)+ptptd(nim1,njm1,k)*tpdt
                  ptb2i=ptpgpv(nim2,njm1,k)+ptptd(nim2,njm1,k)*tpdt
                  ptb2j=ptpgpv(nim1,njm2,k)+ptptd(nim1,njm2,k)*tpdt

                  radwe=((ptpp(nim2,njm1,k)-ptpp(nim1,njm1,k))          &
     &              -(ptb2i-ptb1))*ptcpx(njm1,k,2)

                  radsn=((ptpp(nim1,njm2,k)-ptpp(nim1,njm1,k))          &
     &              -(ptb2j-ptb1))*ptcpy(nim1,k,2)

                  ptpf(nim1,njm1,k)=ptpp(nim1,njm1,k)                   &
     &              +ptptd(nim1,njm1,k)*dt+(radwe+radsn)                &
     &              -dmpdt*(ptpp(nim1,njm1,k)-ptb1)

                end do

!$omp end do

              end if

            end if

          end if

        end if

      end if

! -----

! Force the west boundary value to the external boundary value.

      if(ebw.eq.1.and.isub.eq.0) then

        if(abs(wbc).ne.1) then

          if(advopt.le.3) then

            if(exbvar(5:5).eq.'-') then

!$omp do schedule(runtime) private(j,k,ptb1,ptb2,gamma)

              do k=2,nk-2
              do j=2,nj-2
                gamma=2.e0*ptcpx(j,k,1)/(1.e0-ptcpx(j,k,1))

                ptb1=ptpgpv(1,j,k)+ptptd(1,j,k)*gtinc1
                ptb2=ptpgpv(2,j,k)+ptptd(2,j,k)*gtinc2

                ptpf(1,j,k)=ptpp(1,j,k)+ptptd(1,j,k)*dt2                &
     &            -gamma*((ptp(2,j,k)-ptpp(1,j,k))-(ptb2-ptb1))         &
     &            -dmpdt*(ptpp(1,j,k)-ptb1)

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                ptpf(1,j,k)=ptpp(1,j,k)+ptptd(1,j,k)*dt2
              end do
              end do

!$omp end do

            end if

          else

            if(exbvar(5:5).eq.'-') then

!$omp do schedule(runtime) private(j,k,ptb1,ptb2)

              do k=2,nk-2
              do j=2,nj-2
                ptb1=ptpgpv(1,j,k)+ptptd(1,j,k)*tpdt
                ptb2=ptpgpv(2,j,k)+ptptd(2,j,k)*tpdt

                ptpf(1,j,k)=ptpp(1,j,k)+ptptd(1,j,k)*dt                 &
     &            -ptcpx(j,k,1)*((ptpp(2,j,k)-ptpp(1,j,k))-(ptb2-ptb1)) &
     &            -dmpdt*(ptpp(1,j,k)-ptb1)

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                ptpf(1,j,k)=ptpp(1,j,k)+ptptd(1,j,k)*dt
              end do
              end do

!$omp end do

            end if

          end if

        end if

      end if

! -----

! Force the east boundary value to the external boundary value.

      if(ebe.eq.1.and.isub.eq.nisub-1) then

        if(abs(ebc).ne.1) then

          if(advopt.le.3) then

            if(exbvar(5:5).eq.'-') then

!$omp do schedule(runtime) private(j,k,ptb1,ptb2,gamma)

              do k=2,nk-2
              do j=2,nj-2
                gamma=2.e0*ptcpx(j,k,2)/(1.e0+ptcpx(j,k,2))

                ptb1=ptpgpv(nim1,j,k)+ptptd(nim1,j,k)*gtinc1
                ptb2=ptpgpv(nim2,j,k)+ptptd(nim2,j,k)*gtinc2

                ptpf(nim1,j,k)=ptpp(nim1,j,k)+ptptd(nim1,j,k)*dt2       &
     &            +gamma*((ptp(nim2,j,k)-ptpp(nim1,j,k))-(ptb2-ptb1))   &
     &            -dmpdt*(ptpp(nim1,j,k)-ptb1)

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                ptpf(nim1,j,k)=ptpp(nim1,j,k)+ptptd(nim1,j,k)*dt2
              end do
              end do

!$omp end do

            end if

          else

            if(exbvar(5:5).eq.'-') then

!$omp do schedule(runtime) private(j,k,ptb1,ptb2)

              do k=2,nk-2
              do j=2,nj-2
                ptb1=ptpgpv(nim1,j,k)+ptptd(nim1,j,k)*tpdt
                ptb2=ptpgpv(nim2,j,k)+ptptd(nim2,j,k)*tpdt

                ptpf(nim1,j,k)=ptpp(nim1,j,k)+ptptd(nim1,j,k)*dt        &
     &            +ptcpx(j,k,2)*((ptpp(nim2,j,k)-ptpp(nim1,j,k))        &
     &            -(ptb2-ptb1))-dmpdt*(ptpp(nim1,j,k)-ptb1)

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                ptpf(nim1,j,k)=ptpp(nim1,j,k)+ptptd(nim1,j,k)*dt
              end do
              end do

!$omp end do

            end if

          end if

        end if

      end if

! -----

! Force the south boundary value to the external boundary value.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(exbvar(5:5).eq.'-') then

          if(advopt.le.3) then

!$omp do schedule(runtime) private(i,k,ptb1,ptb2,gamma)

            do k=2,nk-2
            do i=2,ni-2
              gamma=2.e0*ptcpy(i,k,1)/(1.e0-ptcpy(i,k,1))

              ptb1=ptpgpv(i,1,k)+ptptd(i,1,k)*gtinc1
              ptb2=ptpgpv(i,2,k)+ptptd(i,2,k)*gtinc2

              ptpf(i,1,k)=ptpp(i,1,k)+ptptd(i,1,k)*dt2                  &
     &          -gamma*((ptp(i,2,k)-ptpp(i,1,k))-(ptb2-ptb1))           &
     &          -dmpdt*(ptpp(i,1,k)-ptb1)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k,ptb1,ptb2)

            do k=2,nk-2
            do i=2,ni-2
              ptb1=ptpgpv(i,1,k)+ptptd(i,1,k)*tpdt
              ptb2=ptpgpv(i,2,k)+ptptd(i,2,k)*tpdt

              ptpf(i,1,k)=ptpp(i,1,k)+ptptd(i,1,k)*dt                   &
     &          -ptcpy(i,k,1)*((ptpp(i,2,k)-ptpp(i,1,k))-(ptb2-ptb1))   &
     &          -dmpdt*(ptpp(i,1,k)-ptb1)

            end do
            end do

!$omp end do

          end if

        else

          if(advopt.le.3) then

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=1,ni-1
              ptpf(i,1,k)=ptpp(i,1,k)+ptptd(i,1,k)*dt2
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=1,ni-1
              ptpf(i,1,k)=ptpp(i,1,k)+ptptd(i,1,k)*dt
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Force the north boundary value to the external boundary value.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(exbvar(5:5).eq.'-') then

          if(advopt.le.3) then

!$omp do schedule(runtime) private(i,k,ptb1,ptb2,gamma)

            do k=2,nk-2
            do i=2,ni-2
              gamma=2.e0*ptcpy(i,k,2)/(1.e0+ptcpy(i,k,2))

              ptb1=ptpgpv(i,njm1,k)+ptptd(i,njm1,k)*gtinc1
              ptb2=ptpgpv(i,njm2,k)+ptptd(i,njm2,k)*gtinc2

              ptpf(i,njm1,k)=ptpp(i,njm1,k)+ptptd(i,njm1,k)*dt2         &
     &          +gamma*((ptp(i,njm2,k)-ptpp(i,njm1,k))-(ptb2-ptb1))     &
     &          -dmpdt*(ptpp(i,njm1,k)-ptb1)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k,ptb1,ptb2)

            do k=2,nk-2
            do i=2,ni-2
              ptb1=ptpgpv(i,njm1,k)+ptptd(i,njm1,k)*tpdt
              ptb2=ptpgpv(i,njm2,k)+ptptd(i,njm2,k)*tpdt

              ptpf(i,njm1,k)=ptpp(i,njm1,k)+ptptd(i,njm1,k)*dt          &
     &          +ptcpy(i,k,2)*((ptpp(i,njm2,k)-ptpp(i,njm1,k))          &
     &          -(ptb2-ptb1))-dmpdt*(ptpp(i,njm1,k)-ptb1)

            end do
            end do

!$omp end do

          end if

        else

          if(advopt.le.3) then

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=1,ni-1
              ptpf(i,njm1,k)=ptpp(i,njm1,k)+ptptd(i,njm1,k)*dt2
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=1,ni-1
              ptpf(i,njm1,k)=ptpp(i,njm1,k)+ptptd(i,njm1,k)*dt
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_exbcpt

!-----7--------------------------------------------------------------7--

      end module m_exbcpt
