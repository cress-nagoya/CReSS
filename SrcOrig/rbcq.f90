!***********************************************************************
      module m_rbcq
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
!                   2009/11/13, 2011/09/22, 2013/01/28, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the radiative lateral boundary conditions for optional mixing
!     ratio.

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

      public :: rbcq, s_rbcq

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rbcq

        module procedure s_rbcq

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max
      intrinsic mod
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rbcq(fpgpvvar,fplbcvar,fpwbc,fpebc,fpsbc,fpnbc,      &
     &                  fpnggopt,fplspopt,fpvspopt,fpadvopt,fplbnews,   &
     &                  apg,apl,ivstp,dt,gtinc,ni,nj,nk,q,qp,           &
     &                  qcpx,qcpy,qgpv,qtd,qf)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

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

      integer, intent(in) :: apg
                       ! Pointer of gpvvar

      integer, intent(in) :: apl
                       ! Pointer of lbcvar

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

      real, intent(in) :: q(0:ni+1,0:nj+1,1:nk)
                       ! Optional mixing ratio at present

      real, intent(in) :: qp(0:ni+1,0:nj+1,1:nk)
                       ! Optional mixing ratio at past

      real, intent(in) :: qcpx(1:nj,1:nk,1:2)
                       ! Phase speed of optional mixing ratio
                       ! on west and east boundary

      real, intent(in) :: qcpy(1:ni,1:nk,1:2)
                       ! Phase speed of optional mixing ratio
                       ! on south and north boundary

      real, intent(in) :: qgpv(0:ni+1,0:nj+1,1:nk)
                       ! Optional mixing ratio of GPV data
                       ! at marked time

      real, intent(in) :: qtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! optional mixing ratio of GPV data

! Input and output variable

      real, intent(inout) :: qf(0:ni+1,0:nj+1,1:nk)
                       ! Optional mixing ratio at future

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

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

      real gamma       ! Temporary variable

      real radwe       ! Temporary variable
      real radsn       ! Temporary variable

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(gpvvar)
      call inichar(lbcvar)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
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

        if(lbcvar(apl:apl).eq.'o') then
          dmpdt=2.e0*lbnews*dt
        else
          dmpdt=0.e0
        end if

      else

        if(lbcvar(apl:apl).eq.'o') then
          dmpdt=lbnews*dt
        else
          dmpdt=0.e0
        end if

      end if

      tpdt=gtinc+real(ivstp-1)*dt

! -----

!! Set the radiative lateral boundary conditions.

!$omp parallel default(shared)

! Set the boundary conditions at the four corners.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(advopt.le.3) then

         if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.sbc.ge.4) then

           if(gpvvar(apg:apg).eq.'o'.and.                               &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(q(2,1,k)-qp(1,1,k))                               &
     &           *qcpx(1,k,1)/(1.e0-qcpx(1,k,1))

               radsn=(q(1,2,k)-qp(1,1,k))                               &
     &           *qcpy(1,k,1)/(1.e0-qcpy(1,k,1))

               qf(1,1,k)=max(qp(1,1,k)-2.e0*(radwe+radsn)               &
     &           -dmpdt*(qp(1,1,k)-(qgpv(1,1,k)+qtd(1,1,k)*gtinc)),0.e0)

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(q(2,1,k)-qp(1,1,k))                               &
     &           *qcpx(1,k,1)/(1.e0-qcpx(1,k,1))

               radsn=(q(1,2,k)-qp(1,1,k))                               &
     &           *qcpy(1,k,1)/(1.e0-qcpy(1,k,1))

               qf(1,1,k)                                                &
     &           =max(qp(1,1,k)-2.e0*(radwe+radsn)-dmpdt*qp(1,1,k),0.e0)

             end do

!$omp end do

           end if

         end if

         if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.sbc.ge.4) then

           if(gpvvar(apg:apg).eq.'o'.and.                               &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(q(nim2,1,k)-qp(nim1,1,k))                         &
     &           *qcpx(1,k,2)/(1.e0+qcpx(1,k,2))

               radsn=(q(nim1,2,k)-qp(nim1,1,k))                         &
     &           *qcpy(nim1,k,1)/(1.e0-qcpy(nim1,k,1))

               qf(nim1,1,k)                                             &
     &           =max(0.e0,qp(nim1,1,k)+2.e0*(radwe-radsn)-dmpdt        &
     &           *(qp(nim1,1,k)-(qgpv(nim1,1,k)+qtd(nim1,1,k)*gtinc)))

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(q(nim2,1,k)-qp(nim1,1,k))                         &
     &           *qcpx(1,k,2)/(1.e0+qcpx(1,k,2))

               radsn=(q(nim1,2,k)-qp(nim1,1,k))                         &
     &           *qcpy(nim1,k,1)/(1.e0-qcpy(nim1,k,1))

               qf(nim1,1,k)=max(qp(nim1,1,k)+2.e0*(radwe-radsn)         &
     &           -dmpdt*qp(nim1,1,k),0.e0)

             end do

!$omp end do

           end if

         end if

        else

         if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.sbc.ge.4) then

           if(gpvvar(apg:apg).eq.'o'.and.                               &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qp(2,1,k)-qp(1,1,k))*qcpx(1,k,1)
               radsn=(qp(1,2,k)-qp(1,1,k))*qcpy(1,k,1)

               qf(1,1,k)=max(qp(1,1,k)-(radwe+radsn)                    &
     &           -dmpdt*(qp(1,1,k)-(qgpv(1,1,k)+qtd(1,1,k)*tpdt)),0.e0)

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qp(2,1,k)-qp(1,1,k))*qcpx(1,k,1)
               radsn=(qp(1,2,k)-qp(1,1,k))*qcpy(1,k,1)

               qf(1,1,k)                                                &
     &           =max(qp(1,1,k)-(radwe+radsn)-dmpdt*qp(1,1,k),0.e0)

             end do

!$omp end do

           end if

         end if

         if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.sbc.ge.4) then

           if(gpvvar(apg:apg).eq.'o'.and.                               &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qp(nim2,1,k)-qp(nim1,1,k))*qcpx(1,k,2)
               radsn=(qp(nim1,2,k)-qp(nim1,1,k))*qcpy(nim1,k,1)

               qf(nim1,1,k)=max(0.e0,qp(nim1,1,k)+(radwe-radsn)-dmpdt   &
     &           *(qp(nim1,1,k)-(qgpv(nim1,1,k)+qtd(nim1,1,k)*tpdt)))

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qp(nim2,1,k)-qp(nim1,1,k))*qcpx(1,k,2)
               radsn=(qp(nim1,2,k)-qp(nim1,1,k))*qcpy(nim1,k,1)

               qf(nim1,1,k)=max(qp(nim1,1,k)+(radwe-radsn)              &
     &           -dmpdt*qp(nim1,1,k),0.e0)

             end do

!$omp end do

           end if

         end if

        end if

      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(advopt.le.3) then

         if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.nbc.ge.4) then

           if(gpvvar(apg:apg).eq.'o'.and.                               &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(q(2,njm1,k)-qp(1,njm1,k))                         &
     &           *qcpx(njm1,k,1)/(1.e0-qcpx(njm1,k,1))

               radsn=(q(1,njm2,k)-qp(1,njm1,k))                         &
     &           *qcpy(1,k,2)/(1.e0+qcpy(1,k,2))

               qf(1,njm1,k)                                             &
     &           =max(0.e0,qp(1,njm1,k)-2.e0*(radwe-radsn)-dmpdt        &
     &           *(qp(1,njm1,k)-(qgpv(1,njm1,k)+qtd(1,njm1,k)*gtinc)))

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(q(2,njm1,k)-qp(1,njm1,k))                         &
     &           *qcpx(njm1,k,1)/(1.e0-qcpx(njm1,k,1))

               radsn=(q(1,njm2,k)-qp(1,njm1,k))                         &
     &           *qcpy(1,k,2)/(1.e0+qcpy(1,k,2))

               qf(1,njm1,k)=max(qp(1,njm1,k)-2.e0*(radwe-radsn)         &
     &           -dmpdt*qp(1,njm1,k),0.e0)

             end do

!$omp end do

           end if

         end if

         if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.nbc.ge.4) then

           if(gpvvar(apg:apg).eq.'o'.and.                               &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(q(nim2,njm1,k)-qp(nim1,njm1,k))                   &
     &           *qcpx(njm1,k,2)/(1.e0+qcpx(njm1,k,2))

               radsn=(q(nim1,njm2,k)-qp(nim1,njm1,k))                   &
     &           *qcpy(nim1,k,2)/(1.e0+qcpy(nim1,k,2))

               qf(nim1,njm1,k)=max(qp(nim1,njm1,k)                      &
     &           +2.e0*(radwe+radsn)-dmpdt*(qp(nim1,njm1,k)             &
     &           -(qgpv(nim1,njm1,k)+qtd(nim1,njm1,k)*gtinc)),0.e0)

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(q(nim2,njm1,k)-qp(nim1,njm1,k))                   &
     &           *qcpx(njm1,k,2)/(1.e0+qcpx(njm1,k,2))

               radsn=(q(nim1,njm2,k)-qp(nim1,njm1,k))                   &
     &           *qcpy(nim1,k,2)/(1.e0+qcpy(nim1,k,2))

               qf(nim1,njm1,k)=max(qp(nim1,njm1,k)+2.e0*(radwe+radsn)   &
     &           -dmpdt*qp(nim1,njm1,k),0.e0)

             end do

!$omp end do

           end if

         end if

        else

         if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.nbc.ge.4) then

           if(gpvvar(apg:apg).eq.'o'.and.                               &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qp(2,njm1,k)-qp(1,njm1,k))*qcpx(njm1,k,1)
               radsn=(qp(1,njm2,k)-qp(1,njm1,k))*qcpy(1,k,2)

               qf(1,njm1,k)=max(0.e0,qp(1,njm1,k)-(radwe-radsn)-dmpdt   &
     &           *(qp(1,njm1,k)-(qgpv(1,njm1,k)+qtd(1,njm1,k)*tpdt)))

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qp(2,njm1,k)-qp(1,njm1,k))*qcpx(njm1,k,1)
                radsn=(qp(1,njm2,k)-qp(1,njm1,k))*qcpy(1,k,2)

               qf(1,njm1,k)=max(qp(1,njm1,k)-(radwe-radsn)              &
     &           -dmpdt*qp(1,njm1,k),0.e0)

             end do

!$omp end do

           end if

         end if

         if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.nbc.ge.4) then

           if(gpvvar(apg:apg).eq.'o'.and.                               &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qp(nim2,njm1,k)-qp(nim1,njm1,k))*qcpx(njm1,k,2)
               radsn=(qp(nim1,njm2,k)-qp(nim1,njm1,k))*qcpy(nim1,k,2)

               qf(nim1,njm1,k)=max(qp(nim1,njm1,k)                      &
     &           +(radwe+radsn)-dmpdt*(qp(nim1,njm1,k)                  &
     &           -(qgpv(nim1,njm1,k)+qtd(nim1,njm1,k)*tpdt)),0.e0)

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qp(nim2,njm1,k)-qp(nim1,njm1,k))*qcpx(njm1,k,2)
               radsn=(qp(nim1,njm2,k)-qp(nim1,njm1,k))*qcpy(nim1,k,2)

               qf(nim1,njm1,k)=max(qp(nim1,njm1,k)+(radwe+radsn)        &
     &           -dmpdt*qp(nim1,njm1,k),0.e0)

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

            if(gpvvar(apg:apg).eq.'o'.and.                              &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(j,k,gamma)

              do k=2,nk-2
              do j=2,nj-2
                gamma=2.e0*qcpx(j,k,1)/(1.e0-qcpx(j,k,1))

                qf(1,j,k)=max(0.e0,qp(1,j,k)-gamma*(q(2,j,k)-qp(1,j,k)) &
     &            -dmpdt*(qp(1,j,k)-(qgpv(1,j,k)+qtd(1,j,k)*gtinc)))

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k,gamma)

              do k=2,nk-2
              do j=2,nj-2
                gamma=2.e0*qcpx(j,k,1)/(1.e0-qcpx(j,k,1))

                qf(1,j,k)=max(qp(1,j,k)-gamma*(q(2,j,k)-qp(1,j,k))      &
     &            -dmpdt*qp(1,j,k),0.e0)

              end do
              end do

!$omp end do

            end if

          else

            if(gpvvar(apg:apg).eq.'o'.and.                              &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                qf(1,j,k)                                               &
     &            =max(qp(1,j,k)-qcpx(j,k,1)*(qp(2,j,k)-qp(1,j,k))      &
     &            -dmpdt*(qp(1,j,k)-(qgpv(1,j,k)+qtd(1,j,k)*tpdt)),0.e0)
              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                qf(1,j,k)                                               &
     &            =max(qp(1,j,k)-qcpx(j,k,1)*(qp(2,j,k)-qp(1,j,k))      &
     &            -dmpdt*qp(1,j,k),0.e0)
              end do
              end do

!$omp end do

            end if

          end if

        end if

      end if

! -----

! Set the east boundary conditions.

      if(ebe.eq.1.and.isub.eq.nisub-1) then

        if(ebc.ge.4) then

          if(advopt.le.3) then

            if(gpvvar(apg:apg).eq.'o'.and.                              &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(j,k,gamma)

              do k=2,nk-2
              do j=2,nj-2
                gamma=2.e0*qcpx(j,k,2)/(1.e0+qcpx(j,k,2))

                qf(nim1,j,k)=max(0.e0,qp(nim1,j,k)                      &
     &            +gamma*(q(nim2,j,k)-qp(nim1,j,k))-dmpdt               &
     &            *(qp(nim1,j,k)-(qgpv(nim1,j,k)+qtd(nim1,j,k)*gtinc)))

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k,gamma)

              do k=2,nk-2
              do j=2,nj-2
                gamma=2.e0*qcpx(j,k,2)/(1.e0+qcpx(j,k,2))

                qf(nim1,j,k)=max(0.e0,qp(nim1,j,k)                      &
     &            +gamma*(q(nim2,j,k)-qp(nim1,j,k))-dmpdt*qp(nim1,j,k))

              end do
              end do

!$omp end do

            end if

          else

            if(gpvvar(apg:apg).eq.'o'.and.                              &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                qf(nim1,j,k)=max(0.e0,qp(nim1,j,k)                      &
     &            +qcpx(j,k,2)*(qp(nim2,j,k)-qp(nim1,j,k))-dmpdt        &
     &            *(qp(nim1,j,k)-(qgpv(nim1,j,k)+qtd(nim1,j,k)*tpdt)))
              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                qf(nim1,j,k)=max(qp(nim1,j,k)+qcpx(j,k,2)               &
     &            *(qp(nim2,j,k)-qp(nim1,j,k))-dmpdt*qp(nim1,j,k),0.e0)
              end do
              end do

!$omp end do

            end if

          end if

        end if

      end if

! -----

! Set the south boundary conditions.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(sbc.ge.4) then

          if(advopt.le.3) then

            if(gpvvar(apg:apg).eq.'o'.and.                              &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(i,k,gamma)

              do k=2,nk-2
              do i=2,ni-2
                gamma=2.e0*qcpy(i,k,1)/(1.e0-qcpy(i,k,1))

                qf(i,1,k)=max(0.e0,qp(i,1,k)-gamma*(q(i,2,k)-qp(i,1,k)) &
     &            -dmpdt*(qp(i,1,k)-(qgpv(i,1,k)+qtd(i,1,k)*gtinc)))

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(i,k,gamma)

              do k=2,nk-2
              do i=2,ni-2
                gamma=2.e0*qcpy(i,k,1)/(1.e0-qcpy(i,k,1))

                qf(i,1,k)=max(qp(i,1,k)-gamma*(q(i,2,k)-qp(i,1,k))      &
     &            -dmpdt*qp(i,1,k),0.e0)

              end do
              end do

!$omp end do

            end if

          else

            if(gpvvar(apg:apg).eq.'o'.and.                              &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(i,k)

              do k=2,nk-2
              do i=2,ni-2
                qf(i,1,k)                                               &
     &            =max(qp(i,1,k)-qcpy(i,k,1)*(qp(i,2,k)-qp(i,1,k))      &
     &            -dmpdt*(qp(i,1,k)-(qgpv(i,1,k)+qtd(i,1,k)*tpdt)),0.e0)
              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(i,k)

              do k=2,nk-2
              do i=2,ni-2
                qf(i,1,k)                                               &
     &            =max(qp(i,1,k)-qcpy(i,k,1)*(qp(i,2,k)-qp(i,1,k))      &
     &            -dmpdt*qp(i,1,k),0.e0)
              end do
              end do

!$omp end do

            end if

          end if

        end if

      end if

! -----

! Set the north boundary conditions.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(nbc.ge.4) then

          if(advopt.le.3) then

            if(gpvvar(apg:apg).eq.'o'.and.                              &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(i,k,gamma)

              do k=2,nk-2
              do i=2,ni-2
                gamma=2.e0*qcpy(i,k,2)/(1.e0+qcpy(i,k,2))

                qf(i,njm1,k)=max(0.e0,qp(i,njm1,k)                      &
     &            +gamma*(q(i,njm2,k)-qp(i,njm1,k))-dmpdt               &
     &            *(qp(i,njm1,k)-(qgpv(i,njm1,k)+qtd(i,njm1,k)*gtinc)))

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(i,k,gamma)

              do k=2,nk-2
              do i=2,ni-2
                gamma=2.e0*qcpy(i,k,2)/(1.e0+qcpy(i,k,2))

                qf(i,njm1,k)=max(0.e0,qp(i,njm1,k)                      &
     &            +gamma*(q(i,njm2,k)-qp(i,njm1,k))-dmpdt*qp(i,njm1,k))

              end do
              end do

!$omp end do

            end if

          else

            if(gpvvar(apg:apg).eq.'o'.and.                              &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(i,k)

              do k=2,nk-2
              do i=2,ni-2
                qf(i,njm1,k)=max(0.e0,qp(i,njm1,k)                      &
     &            +qcpy(i,k,2)*(qp(i,njm2,k)-qp(i,njm1,k))-dmpdt        &
     &            *(qp(i,njm1,k)-(qgpv(i,njm1,k)+qtd(i,njm1,k)*tpdt)))
              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(i,k)

              do k=2,nk-2
              do i=2,ni-2
                qf(i,njm1,k)=max(qp(i,njm1,k)+qcpy(i,k,2)               &
     &            *(qp(i,njm2,k)-qp(i,njm1,k))-dmpdt*qp(i,njm1,k),0.e0)
              end do
              end do

!$omp end do

            end if

          end if

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_rbcq

!-----7--------------------------------------------------------------7--

      end module m_rbcq
