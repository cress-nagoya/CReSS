!***********************************************************************
      module m_rbcqv
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
!                   2003/12/12, 2004/01/09, 2004/05/07, 2004/08/20,
!                   2005/01/07, 2006/04/03, 2006/09/21, 2006/12/04,
!                   2007/01/05, 2007/01/31, 2007/03/10, 2007/05/07,
!                   2007/05/21, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2008/12/11, 2009/02/27, 2009/03/23, 2009/11/13,
!                   2011/09/22, 2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the radiative lateral boundary conditions for the water vapor
!     mixing ratio.

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

      public :: rbcqv, s_rbcqv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rbcqv

        module procedure s_rbcqv

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
      subroutine s_rbcqv(fpgpvvar,fplbcvar,fpwbc,fpebc,fpsbc,fpnbc,     &
     &                   fpnggopt,fplspopt,fpvspopt,fpadvopt,fplbnews,  &
     &                   ivstp,dt,gtinc,ni,nj,nk,qvbr,qv,qvp,           &
     &                   qvcpx,qvcpy,qvgpv,qvtd,qvf)
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

      real, intent(in) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at present

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qvcpx(1:nj,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on west and east boundary

      real, intent(in) :: qvcpy(1:ni,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on south and north boundary

      real, intent(in) :: qvgpv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, intent(in) :: qvtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! water vapor mixing ratio of GPV data

! Input and output variable

      real, intent(inout) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

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

        if(lbcvar(6:6).eq.'o') then
          dmpdt=2.e0*lbnews*dt
        else
          dmpdt=0.e0
        end if

      else

        if(lbcvar(6:6).eq.'o') then
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

           if(gpvvar(2:2).eq.'o'.and.                                   &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qv(2,1,k)-qvp(1,1,k))                             &
     &           *qvcpx(1,k,1)/(1.e0-qvcpx(1,k,1))

               radsn=(qv(1,2,k)-qvp(1,1,k))                             &
     &           *qvcpy(1,k,1)/(1.e0-qvcpy(1,k,1))

               qvf(1,1,k)=max(0.e0,qvp(1,1,k)-2.e0*(radwe+radsn)        &
     &           -dmpdt*(qvp(1,1,k)-(qvgpv(1,1,k)+qvtd(1,1,k)*gtinc)))

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qv(2,1,k)-qvp(1,1,k))                             &
     &           *qvcpx(1,k,1)/(1.e0-qvcpx(1,k,1))

               radsn=(qv(1,2,k)-qvp(1,1,k))                             &
     &           *qvcpy(1,k,1)/(1.e0-qvcpy(1,k,1))

               qvf(1,1,k)=max(qvp(1,1,k)-2.e0*(radwe+radsn)             &
     &           -dmpdt*(qvp(1,1,k)-qvbr(1,1,k)),0.e0)

             end do

!$omp end do

           end if

         end if

         if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.sbc.ge.4) then

           if(gpvvar(2:2).eq.'o'.and.                                   &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qv(nim2,1,k)-qvp(nim1,1,k))                       &
     &           *qvcpx(1,k,2)/(1.e0+qvcpx(1,k,2))

               radsn=(qv(nim1,2,k)-qvp(nim1,1,k))                       &
     &           *qvcpy(nim1,k,1)/(1.e0-qvcpy(nim1,k,1))

               qvf(nim1,1,k)                                            &
     &          =max(0.e0,qvp(nim1,1,k)+2.e0*(radwe-radsn)-dmpdt        &
     &          *(qvp(nim1,1,k)-(qvgpv(nim1,1,k)+qvtd(nim1,1,k)*gtinc)))

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qv(nim2,1,k)-qvp(nim1,1,k))                       &
     &           *qvcpx(1,k,2)/(1.e0+qvcpx(1,k,2))

               radsn=(qv(nim1,2,k)-qvp(nim1,1,k))                       &
     &           *qvcpy(nim1,k,1)/(1.e0-qvcpy(nim1,k,1))

               qvf(nim1,1,k)=max(qvp(nim1,1,k)+2.e0*(radwe-radsn)       &
     &           -dmpdt*(qvp(nim1,1,k)-qvbr(nim1,1,k)),0.e0)

             end do

!$omp end do

           end if

         end if

        else

         if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.sbc.ge.4) then

           if(gpvvar(2:2).eq.'o'.and.                                   &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qvp(2,1,k)-qvp(1,1,k))*qvcpx(1,k,1)
               radsn=(qvp(1,2,k)-qvp(1,1,k))*qvcpy(1,k,1)

               qvf(1,1,k)=max(0.e0,qvp(1,1,k)-(radwe+radsn)             &
     &           -dmpdt*(qvp(1,1,k)-(qvgpv(1,1,k)+qvtd(1,1,k)*tpdt)))

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qvp(2,1,k)-qvp(1,1,k))*qvcpx(1,k,1)
               radsn=(qvp(1,2,k)-qvp(1,1,k))*qvcpy(1,k,1)

               qvf(1,1,k)=max(qvp(1,1,k)-(radwe+radsn)                  &
     &           -dmpdt*(qvp(1,1,k)-qvbr(1,1,k)),0.e0)

             end do

!$omp end do

           end if

         end if

         if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.sbc.ge.4) then

           if(gpvvar(2:2).eq.'o'.and.                                   &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qvp(nim2,1,k)-qvp(nim1,1,k))*qvcpx(1,k,2)
               radsn=(qvp(nim1,2,k)-qvp(nim1,1,k))*qvcpy(nim1,k,1)

               qvf(nim1,1,k)=max(0.e0,qvp(nim1,1,k)+(radwe-radsn)-dmpdt &
     &           *(qvp(nim1,1,k)-(qvgpv(nim1,1,k)+qvtd(nim1,1,k)*tpdt)))

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qvp(nim2,1,k)-qvp(nim1,1,k))*qvcpx(1,k,2)
               radsn=(qvp(nim1,2,k)-qvp(nim1,1,k))*qvcpy(nim1,k,1)

               qvf(nim1,1,k)=max(qvp(nim1,1,k)+(radwe-radsn)            &
     &           -dmpdt*(qvp(nim1,1,k)-qvbr(nim1,1,k)),0.e0)

             end do

!$omp end do

           end if

         end if

        end if

      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(advopt.le.3) then

         if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.nbc.ge.4) then

           if(gpvvar(2:2).eq.'o'.and.                                   &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qv(2,njm1,k)-qvp(1,njm1,k))                       &
     &           *qvcpx(njm1,k,1)/(1.e0-qvcpx(njm1,k,1))

               radsn=(qv(1,njm2,k)-qvp(1,njm1,k))                       &
     &           *qvcpy(1,k,2)/(1.e0+qvcpy(1,k,2))

               qvf(1,njm1,k)                                            &
     &          =max(0.e0,qvp(1,njm1,k)-2.e0*(radwe-radsn)-dmpdt        &
     &          *(qvp(1,njm1,k)-(qvgpv(1,njm1,k)+qvtd(1,njm1,k)*gtinc)))

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qv(2,njm1,k)-qvp(1,njm1,k))                       &
     &           *qvcpx(njm1,k,1)/(1.e0-qvcpx(njm1,k,1))

               radsn=(qv(1,njm2,k)-qvp(1,njm1,k))                       &
     &           *qvcpy(1,k,2)/(1.e0+qvcpy(1,k,2))

               qvf(1,njm1,k)=max(qvp(1,njm1,k)-2.e0*(radwe-radsn)       &
     &           -dmpdt*(qvp(1,njm1,k)-qvbr(1,njm1,k)),0.e0)

             end do

!$omp end do

           end if

         end if

         if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.nbc.ge.4) then

           if(gpvvar(2:2).eq.'o'.and.                                   &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qv(nim2,njm1,k)-qvp(nim1,njm1,k))                 &
     &           *qvcpx(njm1,k,2)/(1.e0+qvcpx(njm1,k,2))

               radsn=(qv(nim1,njm2,k)-qvp(nim1,njm1,k))                 &
     &           *qvcpy(nim1,k,2)/(1.e0+qvcpy(nim1,k,2))

               qvf(nim1,njm1,k)=max(qvp(nim1,njm1,k)                    &
     &           +2.e0*(radwe+radsn)-dmpdt*(qvp(nim1,njm1,k)            &
     &           -(qvgpv(nim1,njm1,k)+qvtd(nim1,njm1,k)*gtinc)),0.e0)

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qv(nim2,njm1,k)-qvp(nim1,njm1,k))                 &
     &           *qvcpx(njm1,k,2)/(1.e0+qvcpx(njm1,k,2))

               radsn=(qv(nim1,njm2,k)-qvp(nim1,njm1,k))                 &
     &           *qvcpy(nim1,k,2)/(1.e0+qvcpy(nim1,k,2))

               qvf(nim1,njm1,k)=max(qvp(nim1,njm1,k)+2.e0*(radwe+radsn) &
     &           -dmpdt*(qvp(nim1,njm1,k)-qvbr(nim1,njm1,k)),0.e0)

             end do

!$omp end do

           end if

         end if

        else

         if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.nbc.ge.4) then

           if(gpvvar(2:2).eq.'o'.and.                                   &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qvp(2,njm1,k)-qvp(1,njm1,k))*qvcpx(njm1,k,1)
               radsn=(qvp(1,njm2,k)-qvp(1,njm1,k))*qvcpy(1,k,2)

               qvf(1,njm1,k)=max(0.e0,qvp(1,njm1,k)-(radwe-radsn)-dmpdt &
     &           *(qvp(1,njm1,k)-(qvgpv(1,njm1,k)+qvtd(1,njm1,k)*tpdt)))

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qvp(2,njm1,k)-qvp(1,njm1,k))*qvcpx(njm1,k,1)
               radsn=(qvp(1,njm2,k)-qvp(1,njm1,k))*qvcpy(1,k,2)

               qvf(1,njm1,k)=max(qvp(1,njm1,k)-(radwe-radsn)            &
     &           -dmpdt*(qvp(1,njm1,k)-qvbr(1,njm1,k)),0.e0)

             end do

!$omp end do

           end if

         end if

         if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.nbc.ge.4) then

           if(gpvvar(2:2).eq.'o'.and.                                   &
     &       (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qvp(nim2,njm1,k)-qvp(nim1,njm1,k))                &
     &           *qvcpx(njm1,k,2)

               radsn=(qvp(nim1,njm2,k)-qvp(nim1,njm1,k))                &
     &           *qvcpy(nim1,k,2)

               qvf(nim1,njm1,k)=max(qvp(nim1,njm1,k)                    &
     &           +(radwe+radsn)-dmpdt*(qvp(nim1,njm1,k)                 &
     &           -(qvgpv(nim1,njm1,k)+qvtd(nim1,njm1,k)*tpdt)),0.e0)

             end do

!$omp end do

           else

!$omp do schedule(runtime) private(k,radwe,radsn)

             do k=2,nk-2
               radwe=(qvp(nim2,njm1,k)-qvp(nim1,njm1,k))                &
     &           *qvcpx(njm1,k,2)

               radsn=(qvp(nim1,njm2,k)-qvp(nim1,njm1,k))                &
     &           *qvcpy(nim1,k,2)

               qvf(nim1,njm1,k)=max(qvp(nim1,njm1,k)+(radwe+radsn)      &
     &           -dmpdt*(qvp(nim1,njm1,k)-qvbr(nim1,njm1,k)),0.e0)

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

            if(gpvvar(2:2).eq.'o'.and.                                  &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(j,k,gamma)

              do k=2,nk-2
              do j=2,nj-2
                gamma=2.e0*qvcpx(j,k,1)/(1.e0-qvcpx(j,k,1))

                qvf(1,j,k)                                              &
     &            =max(0.e0,qvp(1,j,k)-gamma*(qv(2,j,k)-qvp(1,j,k))     &
     &            -dmpdt*(qvp(1,j,k)-(qvgpv(1,j,k)+qvtd(1,j,k)*gtinc)))

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k,gamma)

              do k=2,nk-2
              do j=2,nj-2
                gamma=2.e0*qvcpx(j,k,1)/(1.e0-qvcpx(j,k,1))

                qvf(1,j,k)=max(qvp(1,j,k)-gamma*(qv(2,j,k)-qvp(1,j,k))  &
     &            -dmpdt*(qvp(1,j,k)-qvbr(1,j,k)),0.e0)

              end do
              end do

!$omp end do

            end if

          else

            if(gpvvar(2:2).eq.'o'.and.                                  &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                qvf(1,j,k)=max(0.e0,qvp(1,j,k)                          &
     &            -qvcpx(j,k,1)*(qvp(2,j,k)-qvp(1,j,k))                 &
     &            -dmpdt*(qvp(1,j,k)-(qvgpv(1,j,k)+qvtd(1,j,k)*tpdt)))
              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                qvf(1,j,k)                                              &
     &            =max(qvp(1,j,k)-qvcpx(j,k,1)*(qvp(2,j,k)-qvp(1,j,k))  &
     &            -dmpdt*(qvp(1,j,k)-qvbr(1,j,k)),0.e0)
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

            if(gpvvar(2:2).eq.'o'.and.                                  &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(j,k,gamma)

              do k=2,nk-2
              do j=2,nj-2
               gamma=2.e0*qvcpx(j,k,2)/(1.e0+qvcpx(j,k,2))

               qvf(nim1,j,k)=max(0.e0,qvp(nim1,j,k)                     &
     &          +gamma*(qv(nim2,j,k)-qvp(nim1,j,k))-dmpdt               &
     &          *(qvp(nim1,j,k)-(qvgpv(nim1,j,k)+qvtd(nim1,j,k)*gtinc)))

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k,gamma)

              do k=2,nk-2
              do j=2,nj-2
                gamma=2.e0*qvcpx(j,k,2)/(1.e0+qvcpx(j,k,2))

                qvf(nim1,j,k)=max(qvp(nim1,j,k)                         &
     &            +gamma*(qv(nim2,j,k)-qvp(nim1,j,k))                   &
     &            -dmpdt*(qvp(nim1,j,k)-qvbr(nim1,j,k)),0.e0)

              end do
              end do

!$omp end do

            end if

          else

            if(gpvvar(2:2).eq.'o'.and.                                  &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                qvf(nim1,j,k)=max(0.e0,qvp(nim1,j,k)                    &
     &           +qvcpx(j,k,2)*(qvp(nim2,j,k)-qvp(nim1,j,k))-dmpdt      &
     &           *(qvp(nim1,j,k)-(qvgpv(nim1,j,k)+qvtd(nim1,j,k)*tpdt)))
              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k,gamma)

              do k=2,nk-2
              do j=2,nj-2
                gamma=2.e0*qvcpx(j,k,2)/(1.e0+qvcpx(j,k,2))

                qvf(nim1,j,k)=max(qvp(nim1,j,k)                         &
     &            +qvcpx(j,k,2)*(qvp(nim2,j,k)-qvp(nim1,j,k))           &
     &            -dmpdt*(qvp(nim1,j,k)-qvbr(nim1,j,k)),0.e0)

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

            if(gpvvar(2:2).eq.'o'.and.                                  &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(i,k,gamma)

              do k=2,nk-2
              do i=2,ni-2
                gamma=2.e0*qvcpy(i,k,1)/(1.e0-qvcpy(i,k,1))

                qvf(i,1,k)                                              &
     &            =max(0.e0,qvp(i,1,k)-gamma*(qv(i,2,k)-qvp(i,1,k))     &
     &            -dmpdt*(qvp(i,1,k)-(qvgpv(i,1,k)+qvtd(i,1,k)*gtinc)))

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(i,k,gamma)

              do k=2,nk-2
              do i=2,ni-2
                gamma=2.e0*qvcpy(i,k,1)/(1.e0-qvcpy(i,k,1))

                qvf(i,1,k)=max(qvp(i,1,k)-gamma*(qv(i,2,k)-qvp(i,1,k))  &
     &            -dmpdt*(qvp(i,1,k)-qvbr(i,1,k)),0.e0)

              end do
              end do

!$omp end do

            end if

          else

            if(gpvvar(2:2).eq.'o'.and.                                  &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(i,k)

              do k=2,nk-2
              do i=2,ni-2
                qvf(i,1,k)=max(0.e0,qvp(i,1,k)                          &
     &            -qvcpy(i,k,1)*(qvp(i,2,k)-qvp(i,1,k))                 &
     &            -dmpdt*(qvp(i,1,k)-(qvgpv(i,1,k)+qvtd(i,1,k)*tpdt)))
              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(i,k)

              do k=2,nk-2
              do i=2,ni-2
                qvf(i,1,k)                                              &
     &            =max(qvp(i,1,k)-qvcpy(i,k,1)*(qvp(i,2,k)-qvp(i,1,k))  &
     &            -dmpdt*(qvp(i,1,k)-qvbr(i,1,k)),0.e0)
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

            if(gpvvar(2:2).eq.'o'.and.                                  &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(i,k,gamma)

              do k=2,nk-2
              do i=2,ni-2
               gamma=2.e0*qvcpy(i,k,2)/(1.e0+qvcpy(i,k,2))

               qvf(i,njm1,k)=max(0.e0,qvp(i,njm1,k)                     &
     &          +gamma*(qv(i,njm2,k)-qvp(i,njm1,k))-dmpdt               &
     &          *(qvp(i,njm1,k)-(qvgpv(i,njm1,k)+qvtd(i,njm1,k)*gtinc)))

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(i,k,gamma)

              do k=2,nk-2
              do i=2,ni-2
                gamma=2.e0*qvcpy(i,k,2)/(1.e0+qvcpy(i,k,2))

                qvf(i,njm1,k)=max(qvp(i,njm1,k)                         &
     &            +gamma*(qv(i,njm2,k)-qvp(i,njm1,k))                   &
     &            -dmpdt*(qvp(i,njm1,k)-qvbr(i,njm1,k)),0.e0)

              end do
              end do

!$omp end do

            end if

          else

            if(gpvvar(2:2).eq.'o'.and.                                  &
     &        (nggopt.eq.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

!$omp do schedule(runtime) private(i,k)

              do k=2,nk-2
              do i=2,ni-2
                qvf(i,njm1,k)=max(0.e0,qvp(i,njm1,k)                    &
     &           +qvcpy(i,k,2)*(qvp(i,njm2,k)-qvp(i,njm1,k))-dmpdt      &
     &           *(qvp(i,njm1,k)-(qvgpv(i,njm1,k)+qvtd(i,njm1,k)*tpdt)))
              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(i,k)

              do k=2,nk-2
              do i=2,ni-2
                qvf(i,njm1,k)=max(qvp(i,njm1,k)                         &
     &            +qvcpy(i,k,2)*(qvp(i,njm2,k)-qvp(i,njm1,k))           &
     &            -dmpdt*(qvp(i,njm1,k)-qvbr(i,njm1,k)),0.e0)
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

      end subroutine s_rbcqv

!-----7--------------------------------------------------------------7--

      end module m_rbcqv
