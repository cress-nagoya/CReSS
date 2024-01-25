!***********************************************************************
      module m_exbcq
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
!                   2004/05/07, 2004/08/20, 2005/01/07, 2006/04/03,
!                   2006/09/21, 2006/12/04, 2007/01/05, 2007/01/31,
!                   2007/05/07, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2008/12/11, 2009/02/27, 2009/03/23, 2009/11/13,
!                   2011/09/22, 2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     force the lateral boundary value to the external boundary value
!     for optional mixing ratio.

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

      public :: exbcq, s_exbcq

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface exbcq

        module procedure s_exbcq

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic max
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_exbcq(fpexbvar,fpwbc,fpebc,fpadvopt,fpexnews,        &
     &                   ape,ivstp,dt,gtinc,ni,nj,nk,q,qp,              &
     &                   qcpx,qcpy,qgpv,qtd,qf)
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

      integer, intent(in) :: ape
                       ! Pointer of exbvar

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

      real qb1         ! Temporary variable
      real qb2         ! Temporary variable
      real qb2i        ! Temporary variable
      real qb2j        ! Temporary variable

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

        if(exbvar(ape:ape).eq.'-') then

          if(abs(wbc).ne.1.and.abs(ebc).ne.1) then

            if(advopt.le.3) then

              if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)

                do k=2,nk-2
                  qb1=qgpv(1,1,k)+qtd(1,1,k)*gtinc1
                  qb2i=qgpv(2,1,k)+qtd(2,1,k)*gtinc2
                  qb2j=qgpv(1,2,k)+qtd(1,2,k)*gtinc2

                  radwe=((q(2,1,k)-qp(1,1,k))-(qb2i-qb1))               &
     &              *qcpx(1,k,1)/(1.e0-qcpx(1,k,1))

                  radsn=((q(1,2,k)-qp(1,1,k))-(qb2j-qb1))               &
     &              *qcpy(1,k,1)/(1.e0-qcpy(1,k,1))

                  qf(1,1,k)=max(qp(1,1,k)+qtd(1,1,k)*dt2                &
     &              -2.e0*(radwe+radsn)-dmpdt*(qp(1,1,k)-qb1),0.e0)

                end do

!$omp end do

              end if

              if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)

                do k=2,nk-2
                  qb1=qgpv(nim1,1,k)+qtd(nim1,1,k)*gtinc1
                  qb2i=qgpv(nim2,1,k)+qtd(nim2,1,k)*gtinc2
                  qb2j=qgpv(nim1,2,k)+qtd(nim1,2,k)*gtinc2

                  radwe=((q(nim2,1,k)-qp(nim1,1,k))-(qb2i-qb1))         &
     &              *qcpx(1,k,2)/(1.e0+qcpx(1,k,2))

                  radsn=((q(nim1,2,k)-qp(nim1,1,k))-(qb2j-qb1))         &
     &              *qcpy(nim1,k,1)/(1.e0-qcpy(nim1,k,1))

                  qf(nim1,1,k)=max(qp(nim1,1,k)+qtd(nim1,1,k)*dt2       &
     &              +2.e0*(radwe-radsn)-dmpdt*(qp(nim1,1,k)-qb1),0.e0)

                end do

!$omp end do

              end if

            else

              if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)

                do k=2,nk-2
                  qb1=qgpv(1,1,k)+qtd(1,1,k)*tpdt
                  qb2i=qgpv(2,1,k)+qtd(2,1,k)*tpdt
                  qb2j=qgpv(1,2,k)+qtd(1,2,k)*tpdt

                  radwe=((qp(2,1,k)-qp(1,1,k))-(qb2i-qb1))*qcpx(1,k,1)
                  radsn=((qp(1,2,k)-qp(1,1,k))-(qb2j-qb1))*qcpy(1,k,1)

                  qf(1,1,k)=max(qp(1,1,k)+qtd(1,1,k)*dt                 &
     &              -(radwe+radsn)-dmpdt*(qp(1,1,k)-qb1),0.e0)

                end do

!$omp end do

              end if

              if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)

                do k=2,nk-2
                  qb1=qgpv(nim1,1,k)+qtd(nim1,1,k)*tpdt
                  qb2i=qgpv(nim2,1,k)+qtd(nim2,1,k)*tpdt
                  qb2j=qgpv(nim1,2,k)+qtd(nim1,2,k)*tpdt

                  radwe=((qp(nim2,1,k)-qp(nim1,1,k))-(qb2i-qb1))        &
     &              *qcpx(1,k,2)

                  radsn=((qp(nim1,2,k)-qp(nim1,1,k))-(qb2j-qb1))        &
     &              *qcpy(nim1,k,1)

                  qf(nim1,1,k)=max(qp(nim1,1,k)+qtd(nim1,1,k)*dt        &
     &              +(radwe-radsn)-dmpdt*(qp(nim1,1,k)-qb1),0.e0)

                end do

!$omp end do

              end if

            end if

          end if

        end if

      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(exbvar(ape:ape).eq.'-') then

          if(abs(wbc).ne.1.and.abs(ebc).ne.1) then

            if(advopt.le.3) then

              if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)

                do k=2,nk-2
                  qb1=qgpv(1,njm1,k)+qtd(1,njm1,k)*gtinc1
                  qb2i=qgpv(2,njm1,k)+qtd(2,njm1,k)*gtinc2
                  qb2j=qgpv(1,njm2,k)+qtd(1,njm2,k)*gtinc2

                  radwe=((q(2,njm1,k)-qp(1,njm1,k))-(qb2i-qb1))         &
     &              *qcpx(njm1,k,1)/(1.e0-qcpx(njm1,k,1))

                  radsn=((q(1,njm2,k)-qp(1,njm1,k))-(qb2j-qb1))         &
     &              *qcpy(1,k,2)/(1.e0+qcpy(1,k,2))

                  qf(1,njm1,k)=max(qp(1,njm1,k)+qtd(1,njm1,k)*dt2       &
     &              -2.e0*(radwe-radsn)-dmpdt*(qp(1,njm1,k)-qb1),0.e0)

                end do

!$omp end do

              end if

              if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)

                do k=2,nk-2
                  qb1=qgpv(nim1,njm1,k)+qtd(nim1,njm1,k)*gtinc1
                  qb2i=qgpv(nim2,njm1,k)+qtd(nim2,njm1,k)*gtinc2
                  qb2j=qgpv(nim1,njm2,k)+qtd(nim1,njm2,k)*gtinc2

                  radwe=((q(nim2,njm1,k)-qp(nim1,njm1,k))-(qb2i-qb1))   &
     &              *qcpx(njm1,k,2)/(1.e0+qcpx(njm1,k,2))

                  radsn=((q(nim1,njm2,k)-qp(nim1,njm1,k))-(qb2j-qb1))   &
     &              *qcpy(nim1,k,2)/(1.e0+qcpy(nim1,k,2))

                  qf(nim1,njm1,k)                                       &
     &             =max(qp(nim1,njm1,k)+qtd(nim1,njm1,k)*dt2            &
     &             +2.e0*(radwe+radsn)-dmpdt*(qp(nim1,njm1,k)-qb1),0.e0)

                end do

!$omp end do

              end if

            else

              if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)

                do k=2,nk-2
                  qb1=qgpv(1,njm1,k)+qtd(1,njm1,k)*tpdt
                  qb2i=qgpv(2,njm1,k)+qtd(2,njm1,k)*tpdt
                  qb2j=qgpv(1,njm2,k)+qtd(1,njm2,k)*tpdt

                  radwe=((qp(2,njm1,k)-qp(1,njm1,k))-(qb2i-qb1))        &
     &              *qcpx(njm1,k,1)

                  radsn=((qp(1,njm2,k)-qp(1,njm1,k))-(qb2j-qb1))        &
     &              *qcpy(1,k,2)

                  qf(1,njm1,k)=max(qp(1,njm1,k)+qtd(1,njm1,k)*dt        &
     &              -(radwe-radsn)-dmpdt*(qp(1,njm1,k)-qb1),0.e0)

                end do

!$omp end do

              end if

              if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime) private(k,qb1,qb2i,qb2j,radwe,radsn)

                do k=2,nk-2
                  qb1=qgpv(nim1,njm1,k)+qtd(nim1,njm1,k)*tpdt
                  qb2i=qgpv(nim2,njm1,k)+qtd(nim2,njm1,k)*tpdt
                  qb2j=qgpv(nim1,njm2,k)+qtd(nim1,njm2,k)*tpdt

                  radwe=((qp(nim2,njm1,k)-qp(nim1,njm1,k))-(qb2i-qb1))  &
     &              *qcpx(njm1,k,2)

                  radsn=((qp(nim1,njm2,k)-qp(nim1,njm1,k))-(qb2j-qb1))  &
     &              *qcpy(nim1,k,2)

                  qf(nim1,njm1,k)                                       &
     &              =max(qp(nim1,njm1,k)+qtd(nim1,njm1,k)*dt            &
     &              +(radwe+radsn)-dmpdt*(qp(nim1,njm1,k)-qb1),0.e0)

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

            if(exbvar(ape:ape).eq.'-') then

!$omp do schedule(runtime) private(j,k,qb1,qb2,gamma)

              do k=2,nk-2
              do j=2,nj-2
                gamma=2.e0*qcpx(j,k,1)/(1.e0-qcpx(j,k,1))

                qb1=qgpv(1,j,k)+qtd(1,j,k)*gtinc1
                qb2=qgpv(2,j,k)+qtd(2,j,k)*gtinc2

                qf(1,j,k)=max(qp(1,j,k)+qtd(1,j,k)*dt2                  &
     &            -gamma*((q(2,j,k)-qp(1,j,k))-(qb2-qb1))               &
     &            -dmpdt*(qp(1,j,k)-qb1),0.e0)

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                qf(1,j,k)=max(qp(1,j,k)+qtd(1,j,k)*dt2,0.e0)
              end do
              end do

!$omp end do

            end if

          else

            if(exbvar(ape:ape).eq.'-') then

!$omp do schedule(runtime) private(j,k,qb1,qb2)

              do k=2,nk-2
              do j=2,nj-2
                qb1=qgpv(1,j,k)+qtd(1,j,k)*tpdt
                qb2=qgpv(2,j,k)+qtd(2,j,k)*tpdt

                qf(1,j,k)=max(qp(1,j,k)+qtd(1,j,k)*dt                   &
     &            -qcpx(j,k,1)*((qp(2,j,k)-qp(1,j,k))-(qb2-qb1))        &
     &            -dmpdt*(qp(1,j,k)-qb1),0.e0)

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                qf(1,j,k)=max(qp(1,j,k)+qtd(1,j,k)*dt,0.e0)
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

            if(exbvar(ape:ape).eq.'-') then

!$omp do schedule(runtime) private(j,k,qb1,qb2,gamma)

              do k=2,nk-2
              do j=2,nj-2
                gamma=2.e0*qcpx(j,k,2)/(1.e0+qcpx(j,k,2))

                qb1=qgpv(nim1,j,k)+qtd(nim1,j,k)*gtinc1
                qb2=qgpv(nim2,j,k)+qtd(nim2,j,k)*gtinc2

                qf(nim1,j,k)=max(qp(nim1,j,k)+qtd(nim1,j,k)*dt2         &
     &            +gamma*((q(nim2,j,k)-qp(nim1,j,k))-(qb2-qb1))         &
     &            -dmpdt*(qp(nim1,j,k)-qb1),0.e0)

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                qf(nim1,j,k)=max(qp(nim1,j,k)+qtd(nim1,j,k)*dt2,0.e0)
              end do
              end do

!$omp end do

            end if

          else

            if(exbvar(ape:ape).eq.'-') then

!$omp do schedule(runtime) private(j,k,qb1,qb2)

              do k=2,nk-2
              do j=2,nj-2
                qb1=qgpv(nim1,j,k)+qtd(nim1,j,k)*tpdt
                qb2=qgpv(nim2,j,k)+qtd(nim2,j,k)*tpdt

                qf(nim1,j,k)=max(qp(nim1,j,k)+qtd(nim1,j,k)*dt          &
     &            +qcpx(j,k,2)*((qp(nim2,j,k)-qp(nim1,j,k))-(qb2-qb1))  &
     &            -dmpdt*(qp(nim1,j,k)-qb1),0.e0)

              end do
              end do

!$omp end do

            else

!$omp do schedule(runtime) private(j,k)

              do k=2,nk-2
              do j=2,nj-2
                qf(nim1,j,k)=max(qp(nim1,j,k)+qtd(nim1,j,k)*dt,0.e0)
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

        if(exbvar(ape:ape).eq.'-') then

          if(advopt.le.3) then

!$omp do schedule(runtime) private(i,k,qb1,qb2,gamma)

            do k=2,nk-2
            do i=2,ni-2
              gamma=2.e0*qcpy(i,k,1)/(1.e0-qcpy(i,k,1))

              qb1=qgpv(i,1,k)+qtd(i,1,k)*gtinc1
              qb2=qgpv(i,2,k)+qtd(i,2,k)*gtinc2

              qf(i,1,k)=max(qp(i,1,k)+qtd(i,1,k)*dt2                    &
     &          -gamma*((q(i,2,k)-qp(i,1,k))-(qb2-qb1))                 &
     &          -dmpdt*(qp(i,1,k)-qb1),0.e0)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k,qb1,qb2)

            do k=2,nk-2
            do i=2,ni-2
              qb1=qgpv(i,1,k)+qtd(i,1,k)*tpdt
              qb2=qgpv(i,2,k)+qtd(i,2,k)*tpdt

              qf(i,1,k)=max(qp(i,1,k)+qtd(i,1,k)*dt                     &
     &          -qcpy(i,k,1)*((qp(i,2,k)-qp(i,1,k))-(qb2-qb1))          &
     &          -dmpdt*(qp(i,1,k)-qb1),0.e0)

            end do
            end do

!$omp end do

          end if

        else

          if(advopt.le.3) then

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=1,ni-1
              qf(i,1,k)=max(qp(i,1,k)+qtd(i,1,k)*dt2,0.e0)
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=1,ni-1
              qf(i,1,k)=max(qp(i,1,k)+qtd(i,1,k)*dt,0.e0)
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Force the north boundary value to the external boundary value.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(exbvar(ape:ape).eq.'-') then

          if(advopt.le.3) then

!$omp do schedule(runtime) private(i,k,qb1,qb2,gamma)

            do k=2,nk-2
            do i=2,ni-2
              gamma=2.e0*qcpy(i,k,2)/(1.e0+qcpy(i,k,2))

              qb1=qgpv(i,njm1,k)+qtd(i,njm1,k)*gtinc1
              qb2=qgpv(i,njm2,k)+qtd(i,njm2,k)*gtinc2

              qf(i,njm1,k)=max(qp(i,njm1,k)+qtd(i,njm1,k)*dt2           &
     &          +gamma*((q(i,njm2,k)-qp(i,njm1,k))-(qb2-qb1))           &
     &          -dmpdt*(qp(i,njm1,k)-qb1),0.e0)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k,qb1,qb2)

            do k=2,nk-2
            do i=2,ni-2
              qb1=qgpv(i,njm1,k)+qtd(i,njm1,k)*tpdt
              qb2=qgpv(i,njm2,k)+qtd(i,njm2,k)*tpdt

              qf(i,njm1,k)=max(qp(i,njm1,k)+qtd(i,njm1,k)*dt            &
     &          +qcpy(i,k,2)*((qp(i,njm2,k)-qp(i,njm1,k))-(qb2-qb1))    &
     &          -dmpdt*(qp(i,njm1,k)-qb1),0.e0)

            end do
            end do

!$omp end do

          end if

        else

          if(advopt.le.3) then

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=1,ni-1
              qf(i,njm1,k)=max(qp(i,njm1,k)+qtd(i,njm1,k)*dt2,0.e0)
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=1,ni-1
              qf(i,njm1,k)=max(qp(i,njm1,k)+qtd(i,njm1,k)*dt,0.e0)
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_exbcq

!-----7--------------------------------------------------------------7--

      end module m_exbcq
