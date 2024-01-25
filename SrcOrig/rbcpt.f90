!***********************************************************************
      module m_rbcpt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/28,
!                   1999/08/03, 1999/08/18, 1999/08/23, 1999/09/01,
!                   1999/09/06, 1999/09/30, 1999/10/07, 1999/11/01,
!                   1999/12/06, 2000/01/17, 2000/03/17, 2001/09/13,
!                   2001/12/11, 2002/04/02, 2002/07/23, 2002/08/15,
!                   2002/10/31, 2003/03/21, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2003/06/27, 2003/11/05, 2003/11/28,
!                   2003/12/12, 2004/01/09, 2004/04/15, 2004/05/07,
!                   2004/08/20, 2005/01/07, 2006/04/03, 2006/09/21,
!                   2006/12/04, 2007/01/05, 2007/01/31, 2007/03/10,
!                   2007/05/07, 2007/05/21, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2008/12/11, 2009/02/27, 2009/03/23,
!                   2011/09/22, 2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the radiative lateral boundary conditions for potential
!     temperature.

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

      public :: rbcpt, s_rbcpt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rbcpt

        module procedure s_rbcpt

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rbcpt(fplbcvar,fpwbc,fpebc,fpsbc,fpnbc,fpnggopt,     &
     &                   fplspopt,fpvspopt,fpadvopt,fplbnews,           &
     &                   ivstp,dt,gtinc,ni,nj,nk,ptp,ptpp,              &
     &                   ptcpx,ptcpy,ptpgpv,ptptd,ptpf)
!***********************************************************************

! Input variables

      integer, intent(in) :: fplbcvar
                       ! Formal parameter of unique index of lbcvar

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: fpnggopt
                       ! Formal parameter of unique index of nggopt

      integer, intent(in) :: fplspopt
                       ! Formal parameter of unique index of lspopt

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fplbnews
                       ! Formal parameter of unique index of lbnews

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

      character(len=108) lbcvar
                       ! Control flag of
                       ! lateral boundary forced variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer nggopt   ! Option for analysis nudging to GPV
      integer lspopt   ! Option for lateral sponge damping
      integer vspopt   ! Option for vertical sponge damping
      integer advopt   ! Option for advection scheme

      integer nim1     ! ni - 1
      integer nim2     ! ni - 2
      integer njm1     ! nj - 1
      integer njm2     ! nj - 2

      real lbnews      ! Boundary damping coefficient

      real dmpdt       ! lbnews x dt or 2.0 x lbnews x dt

      real tpdt        ! gtinc + real(ivstp - 1) x dt

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real dgpv        ! Temporary variable

      real gamma       ! Temporary variable

      real radwe       ! Temporary variable
      real radsn       ! Temporary variable

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(lbcvar)

! -----

! Get the required namelist variables.

      call getcname(fplbcvar,lbcvar)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)
      call getiname(fpnggopt,nggopt)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getiname(fpadvopt,advopt)
      call getrname(fplbnews,lbnews)

! -----

! Set the common used variables.

      nim1=ni-1
      nim2=ni-2
      njm1=nj-1
      njm2=nj-2

      if(advopt.le.3) then

        if(lbcvar(5:5).eq.'o') then
          dmpdt=2.e0*lbnews*dt
        else
          dmpdt=0.e0
        end if

      else

        if(lbcvar(5:5).eq.'o') then
          dmpdt=lbnews*dt
        else
          dmpdt=0.e0
        end if

      end if

      tpdt=gtinc+real(ivstp-1)*dt

! -----

!! Set the radiative lateral boundary conditions.

!$omp parallel default(shared) private(k)

! Set the boundary conditions at the four corners.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(advopt.le.3) then

          if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.sbc.ge.4) then

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(dgpv,radwe,radsn)

              do k=2,nk-2
                dgpv=ptpp(1,1,k)-(ptpgpv(1,1,k)+ptptd(1,1,k)*gtinc)

                radwe=(ptp(2,1,k)-ptpp(1,1,k))                          &
     &            *ptcpx(1,k,1)/(1.e0-ptcpx(1,k,1))

                radsn=(ptp(1,2,k)-ptpp(1,1,k))                          &
     &            *ptcpy(1,k,1)/(1.e0-ptcpy(1,k,1))

                ptpf(1,1,k)=ptpp(1,1,k)-dmpdt*dgpv-2.e0*(radwe+radsn)

              end do

!$omp end do

            else

!$omp do schedule(runtime) private(radwe,radsn)

              do k=2,nk-2
                radwe=(ptp(2,1,k)-ptpp(1,1,k))                          &
     &            *ptcpx(1,k,1)/(1.e0-ptcpx(1,k,1))

                radsn=(ptp(1,2,k)-ptpp(1,1,k))                          &
     &            *ptcpy(1,k,1)/(1.e0-ptcpy(1,k,1))

                ptpf(1,1,k)=ptpp(1,1,k)-dmpdt*ptpp(1,1,k)               &
     &            -2.e0*(radwe+radsn)

              end do

!$omp end do

            end if

          end if

          if(ebe.eq.1.and.isub.eq.nisub-1                               &
     &      .and.ebc.ge.4.and.sbc.ge.4) then

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(dgpv,radwe,radsn)

              do k=2,nk-2
                dgpv=ptpp(nim1,1,k)                                     &
     &            -(ptpgpv(nim1,1,k)+ptptd(nim1,1,k)*gtinc)

                radwe=(ptp(nim2,1,k)-ptpp(nim1,1,k))                    &
     &            *ptcpx(1,k,2)/(1.e0+ptcpx(1,k,2))

                radsn=(ptp(nim1,2,k)-ptpp(nim1,1,k))                    &
     &            *ptcpy(nim1,k,1)/(1.e0-ptcpy(nim1,k,1))

                ptpf(nim1,1,k)=ptpp(nim1,1,k)-dmpdt*dgpv                &
     &            +2.e0*(radwe-radsn)

              end do

!$omp end do

            else

!$omp do schedule(runtime) private(radwe,radsn)

              do k=2,nk-2
                radwe=(ptp(nim2,1,k)-ptpp(nim1,1,k))                    &
     &            *ptcpx(1,k,2)/(1.e0+ptcpx(1,k,2))

                radsn=(ptp(nim1,2,k)-ptpp(nim1,1,k))                    &
     &            *ptcpy(nim1,k,1)/(1.e0-ptcpy(nim1,k,1))

                ptpf(nim1,1,k)=ptpp(nim1,1,k)-dmpdt*ptpp(nim1,1,k)      &
     &            +2.e0*(radwe-radsn)

              end do

!$omp end do

            end if

          end if

        else

          if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.sbc.ge.4) then

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(dgpv,radwe,radsn)

              do k=2,nk-2
                dgpv=ptpp(1,1,k)-(ptpgpv(1,1,k)+ptptd(1,1,k)*tpdt)

                radwe=(ptpp(2,1,k)-ptpp(1,1,k))*ptcpx(1,k,1)
                radsn=(ptpp(1,2,k)-ptpp(1,1,k))*ptcpy(1,k,1)

                ptpf(1,1,k)=ptpp(1,1,k)-dmpdt*dgpv-(radwe+radsn)

              end do

!$omp end do

            else

!$omp do schedule(runtime) private(radwe,radsn)

              do k=2,nk-2
                radwe=(ptpp(2,1,k)-ptpp(1,1,k))*ptcpx(1,k,1)
                radsn=(ptpp(1,2,k)-ptpp(1,1,k))*ptcpy(1,k,1)

                ptpf(1,1,k)=ptpp(1,1,k)-dmpdt*ptpp(1,1,k)-(radwe+radsn)

              end do

!$omp end do

            end if

          end if

          if(ebe.eq.1.and.isub.eq.nisub-1                               &
     &      .and.ebc.ge.4.and.sbc.ge.4) then

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(dgpv,radwe,radsn)

              do k=2,nk-2
                dgpv=ptpp(nim1,1,k)                                     &
     &            -(ptpgpv(nim1,1,k)+ptptd(nim1,1,k)*tpdt)

                radwe=(ptpp(nim2,1,k)-ptpp(nim1,1,k))*ptcpx(1,k,2)
                radsn=(ptpp(nim1,2,k)-ptpp(nim1,1,k))*ptcpy(nim1,k,1)

                ptpf(nim1,1,k)=ptpp(nim1,1,k)-dmpdt*dgpv+(radwe-radsn)

              end do

!$omp end do

            else

!$omp do schedule(runtime) private(radwe,radsn)

              do k=2,nk-2
                radwe=(ptpp(nim2,1,k)-ptpp(nim1,1,k))*ptcpx(1,k,2)
                radsn=(ptpp(nim1,2,k)-ptpp(nim1,1,k))*ptcpy(nim1,k,1)

                ptpf(nim1,1,k)=ptpp(nim1,1,k)-dmpdt*ptpp(nim1,1,k)      &
     &            +(radwe-radsn)

              end do

!$omp end do

            end if

          end if

        end if

      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(advopt.le.3) then

          if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.nbc.ge.4) then

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(dgpv,radwe,radsn)

              do k=2,nk-2
                dgpv=ptpp(1,njm1,k)                                     &
     &            -(ptpgpv(1,njm1,k)+ptptd(1,njm1,k)*gtinc)

                radwe=(ptp(2,njm1,k)-ptpp(1,njm1,k))                    &
     &            *ptcpx(njm1,k,1)/(1.e0-ptcpx(njm1,k,1))

                radsn=(ptp(1,njm2,k)-ptpp(1,njm1,k))                    &
     &            *ptcpy(1,k,2)/(1.e0+ptcpy(1,k,2))

                ptpf(1,njm1,k)=ptpp(1,njm1,k)-dmpdt*dgpv                &
     &            -2.e0*(radwe-radsn)

              end do

!$omp end do

            else

!$omp do schedule(runtime) private(radwe,radsn)

              do k=2,nk-2
                radwe=(ptp(2,njm1,k)-ptpp(1,njm1,k))                    &
     &            *ptcpx(njm1,k,1)/(1.e0-ptcpx(njm1,k,1))

                radsn=(ptp(1,njm2,k)-ptpp(1,njm1,k))                    &
     &            *ptcpy(1,k,2)/(1.e0+ptcpy(1,k,2))

                ptpf(1,njm1,k)=ptpp(1,njm1,k)-dmpdt*ptpp(1,njm1,k)      &
     &            -2.e0*(radwe-radsn)

              end do

!$omp end do

            end if

          end if

          if(ebe.eq.1.and.isub.eq.nisub-1                               &
     &      .and.ebc.ge.4.and.nbc.ge.4) then

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(dgpv,radwe,radsn)

              do k=2,nk-2
                dgpv=ptpp(nim1,njm1,k)                                  &
     &            -(ptpgpv(nim1,njm1,k)+ptptd(nim1,njm1,k)*gtinc)

                radwe=(ptp(nim2,njm1,k)-ptpp(nim1,njm1,k))              &
     &            *ptcpx(njm1,k,2)/(1.e0+ptcpx(njm1,k,2))

                radsn=(ptp(nim1,njm2,k)-ptpp(nim1,njm1,k))              &
     &            *ptcpy(nim1,k,2)/(1.e0+ptcpy(nim1,k,2))

                ptpf(nim1,njm1,k)=ptpp(nim1,njm1,k)-dmpdt*dgpv          &
     &           +2.e0*(radwe+radsn)

              end do

!$omp end do

            else

!$omp do schedule(runtime) private(radwe,radsn)

              do k=2,nk-2
                radwe=(ptp(nim2,njm1,k)-ptpp(nim1,njm1,k))              &
     &            *ptcpx(njm1,k,2)/(1.e0+ptcpx(njm1,k,2))

                radsn=(ptp(nim1,njm2,k)-ptpp(nim1,njm1,k))              &
     &            *ptcpy(nim1,k,2)/(1.e0+ptcpy(nim1,k,2))

                ptpf(nim1,njm1,k)=ptpp(nim1,njm1,k)                     &
     &            -dmpdt*ptpp(nim1,njm1,k)+2.e0*(radwe+radsn)

              end do

!$omp end do

            end if

          end if

        else

          if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.nbc.ge.4) then

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(dgpv,radwe,radsn)

              do k=2,nk-2
                dgpv=ptpp(1,njm1,k)                                     &
     &            -(ptpgpv(1,njm1,k)+ptptd(1,njm1,k)*tpdt)

                radwe=(ptpp(2,njm1,k)-ptpp(1,njm1,k))*ptcpx(njm1,k,1)
                radsn=(ptpp(1,njm2,k)-ptpp(1,njm1,k))*ptcpy(1,k,2)

                ptpf(1,njm1,k)=ptpp(1,njm1,k)-dmpdt*dgpv-(radwe-radsn)

              end do

!$omp end do

           else

!$omp do schedule(runtime) private(radwe,radsn)

              do k=2,nk-2
                radwe=(ptpp(2,njm1,k)-ptpp(1,njm1,k))*ptcpx(njm1,k,1)
                radsn=(ptpp(1,njm2,k)-ptpp(1,njm1,k))*ptcpy(1,k,2)

                ptpf(1,njm1,k)=ptpp(1,njm1,k)-dmpdt*ptpp(1,njm1,k)      &
     &            -(radwe-radsn)

              end do

!$omp end do

            end if

          end if

          if(ebe.eq.1.and.isub.eq.nisub-1                               &
     &      .and.ebc.ge.4.and.nbc.ge.4) then

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

!$omp do schedule(runtime) private(dgpv,radwe,radsn)

              do k=2,nk-2
                dgpv=ptpp(nim1,njm1,k)                                  &
     &            -(ptpgpv(nim1,njm1,k)+ptptd(nim1,njm1,k)*tpdt)

                radwe=(ptpp(nim2,njm1,k)-ptpp(nim1,njm1,k))             &
     &            *ptcpx(njm1,k,2)

                radsn=(ptpp(nim1,njm2,k)-ptpp(nim1,njm1,k))             &
     &            *ptcpy(nim1,k,2)

                ptpf(nim1,njm1,k)=ptpp(nim1,njm1,k)-dmpdt*dgpv          &
     &            +(radwe+radsn)

              end do

!$omp end do

            else

!$omp do schedule(runtime) private(radwe,radsn)

              do k=2,nk-2
                radwe=(ptpp(nim2,njm1,k)-ptpp(nim1,njm1,k))             &
     &            *ptcpx(njm1,k,2)

                radsn=(ptpp(nim1,njm2,k)-ptpp(nim1,njm1,k))             &
     &            *ptcpy(nim1,k,2)

                ptpf(nim1,njm1,k)=ptpp(nim1,njm1,k)                     &
     &            -dmpdt*ptpp(nim1,njm1,k)+(radwe+radsn)

              end do

!$omp end do

            end if

          end if

        end if

      end if

! -----

! Set the west boundary conditions.

      if(ebw.eq.1.and.isub.eq.0) then

        if(wbc.ge.4) then

          if(advopt.le.3) then

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

              do k=2,nk-2

!$omp do schedule(runtime) private(j,dgpv,gamma)

                do j=2,nj-2
                  dgpv=ptpp(1,j,k)-(ptpgpv(1,j,k)+ptptd(1,j,k)*gtinc)

                  gamma=2.e0*ptcpx(j,k,1)/(1.e0-ptcpx(j,k,1))

                  ptpf(1,j,k)=ptpp(1,j,k)-dmpdt*dgpv                    &
     &              -gamma*(ptp(2,j,k)-ptpp(1,j,k))

                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(j,gamma)

                do j=2,nj-2
                  gamma=2.e0*ptcpx(j,k,1)/(1.e0-ptcpx(j,k,1))

                  ptpf(1,j,k)=ptpp(1,j,k)-dmpdt*ptpp(1,j,k)             &
     &              -gamma*(ptp(2,j,k)-ptpp(1,j,k))

                end do

!$omp end do

              end do

            end if

          else

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

              do k=2,nk-2

!$omp do schedule(runtime) private(j,dgpv)

                do j=2,nj-2
                  dgpv=ptpp(1,j,k)-(ptpgpv(1,j,k)+ptptd(1,j,k)*tpdt)

                  ptpf(1,j,k)=ptpp(1,j,k)-dmpdt*dgpv                    &
     &              -ptcpx(j,k,1)*(ptpp(2,j,k)-ptpp(1,j,k))

                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(j)

                do j=2,nj-2
                  ptpf(1,j,k)=ptpp(1,j,k)-dmpdt*ptpp(1,j,k)             &
     &              -ptcpx(j,k,1)*(ptpp(2,j,k)-ptpp(1,j,k))

                end do

!$omp end do

              end do

            end if

          end if

        end if

      end if

! -----

! Set the east boundary conditions.

      if(ebe.eq.1.and.isub.eq.nisub-1) then

        if(ebc.ge.4) then

          if(advopt.le.3) then

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

              do k=2,nk-2

!$omp do schedule(runtime) private(j,dgpv,gamma)

                do j=2,nj-2
                  dgpv=ptpp(nim1,j,k)                                   &
     &              -(ptpgpv(nim1,j,k)+ptptd(nim1,j,k)*gtinc)

                  gamma=2.e0*ptcpx(j,k,2)/(1.e0+ptcpx(j,k,2))

                  ptpf(nim1,j,k)=ptpp(nim1,j,k)-dmpdt*dgpv              &
     &              +gamma*(ptp(nim2,j,k)-ptpp(nim1,j,k))

                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(j,gamma)

                do j=2,nj-2
                  gamma=2.e0*ptcpx(j,k,2)/(1.e0+ptcpx(j,k,2))

                  ptpf(nim1,j,k)=ptpp(nim1,j,k)-dmpdt*ptpp(nim1,j,k)    &
     &              +gamma*(ptp(nim2,j,k)-ptpp(nim1,j,k))

                end do

!$omp end do

              end do

            end if

          else

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

              do k=2,nk-2

!$omp do schedule(runtime) private(j,dgpv)

                do j=2,nj-2
                  dgpv=ptpp(nim1,j,k)                                   &
     &              -(ptpgpv(nim1,j,k)+ptptd(nim1,j,k)*tpdt)

                  ptpf(nim1,j,k)=ptpp(nim1,j,k)-dmpdt*dgpv              &
     &              +ptcpx(j,k,2)*(ptpp(nim2,j,k)-ptpp(nim1,j,k))

                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(j)

                do j=2,nj-2
                  ptpf(nim1,j,k)=ptpp(nim1,j,k)-dmpdt*ptpp(nim1,j,k)    &
     &              +ptcpx(j,k,2)*(ptpp(nim2,j,k)-ptpp(nim1,j,k))

                end do

!$omp end do

              end do

            end if

          end if

        end if

      end if

! -----

! Set the south boundary conditions.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(sbc.ge.4) then

          if(advopt.le.3) then

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,dgpv,gamma)

                do i=2,ni-2
                  dgpv=ptpp(i,1,k)-(ptpgpv(i,1,k)+ptptd(i,1,k)*gtinc)

                  gamma=2.e0*ptcpy(i,k,1)/(1.e0-ptcpy(i,k,1))

                  ptpf(i,1,k)=ptpp(i,1,k)-dmpdt*dgpv                    &
     &              -gamma*(ptp(i,2,k)-ptpp(i,1,k))

                end do

!$omp end do

              end do

            else

!$omp do schedule(runtime) private(i,gamma)

              do k=2,nk-2
              do i=2,ni-2
                gamma=2.e0*ptcpy(i,k,1)/(1.e0-ptcpy(i,k,1))

                ptpf(i,1,k)=ptpp(i,1,k)-dmpdt*ptpp(i,1,k)               &
     &            -gamma*(ptp(i,2,k)-ptpp(i,1,k))

              end do
              end do

!$omp end do

            end if

          else

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,dgpv)

                do i=2,ni-2
                  dgpv=ptpp(i,1,k)-(ptpgpv(i,1,k)+ptptd(i,1,k)*tpdt)

                  ptpf(i,1,k)=ptpp(i,1,k)-dmpdt*dgpv                    &
     &              -ptcpy(i,k,1)*(ptpp(i,2,k)-ptpp(i,1,k))

                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i)

                do i=2,ni-2
                  ptpf(i,1,k)=ptpp(i,1,k)-dmpdt*ptpp(i,1,k)             &
     &              -ptcpy(i,k,1)*(ptpp(i,2,k)-ptpp(i,1,k))

                end do

!$omp end do

              end do

            end if

          end if

        end if

      end if

! -----

! Set the north boundary conditions.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(nbc.ge.4) then

          if(advopt.le.3) then

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,dgpv,gamma)

                do i=2,ni-2
                  dgpv=ptpp(i,njm1,k)                                   &
     &              -(ptpgpv(i,njm1,k)+ptptd(i,njm1,k)*gtinc)

                  gamma=2.e0*ptcpy(i,k,2)/(1.e0+ptcpy(i,k,2))

                  ptpf(i,njm1,k)=ptpp(i,njm1,k)-dmpdt*dgpv              &
     &              +gamma*(ptp(i,njm2,k)-ptpp(i,njm1,k))

                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,gamma)

                do i=2,ni-2
                  gamma=2.e0*ptcpy(i,k,2)/(1.e0+ptcpy(i,k,2))

                  ptpf(i,njm1,k)=ptpp(i,njm1,k)-dmpdt*ptpp(i,njm1,k)    &
     &              +gamma*(ptp(i,njm2,k)-ptpp(i,njm1,k))

                end do

!$omp end do

              end do

            end if

          else

            if(nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,dgpv)

                do i=2,ni-2
                  dgpv=ptpp(i,njm1,k)                                   &
     &              -(ptpgpv(i,njm1,k)+ptptd(i,njm1,k)*tpdt)

                  ptpf(i,njm1,k)=ptpp(i,njm1,k)-dmpdt*dgpv              &
     &              +ptcpy(i,k,2)*(ptpp(i,njm2,k)-ptpp(i,njm1,k))

                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i)

                do i=2,ni-2
                  ptpf(i,njm1,k)=ptpp(i,njm1,k)-dmpdt*ptpp(i,njm1,k)    &
     &              +ptcpy(i,k,2)*(ptpp(i,njm2,k)-ptpp(i,njm1,k))

                end do

!$omp end do

              end do

            end if

          end if

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_rbcpt

!-----7--------------------------------------------------------------7--

      end module m_rbcpt
