!***********************************************************************
      module m_allocbuf
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/07/15, 2003/08/08, 2003/11/05, 2003/12/12,
!                   2004/03/05, 2004/05/31, 2004/08/20, 2005/02/10,
!                   2005/10/05, 2006/01/10, 2006/02/13, 2006/04/03,
!                   2006/07/21, 2006/11/06, 2006/12/04, 2007/01/20,
!                   2007/01/31, 2007/10/19, 2007/11/26, 2008/05/02,
!                   2008/08/25, 2009/01/05, 2009/01/30, 2009/02/27,
!                   2009/11/13, 2011/08/18, 2011/09/22, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the communication buffer.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_combuf
      use m_comgrp
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: allocbuf, s_allocbuf

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allocbuf

        module procedure s_allocbuf

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic max
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_allocbuf(fpwbc,fpebc,fpsbc,fpnbc,                    &
     &                      fpgwmopt,fpadvopt,fpsmtopt,fpcphopt,        &
     &                      fpqcgopt,fpaslopt,fptrkopt,fptubopt,        &
     &                      ni,nj,nk,nqw,nnw,nqi,nni,nqa)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: fpgwmopt
                       ! Formal parameter of unique index of gwmopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fpqcgopt
                       ! Formal parameter of unique index of qcgopt

      integer, intent(in) :: fpaslopt
                       ! Formal parameter of unique index of aslopt

      integer, intent(in) :: fptrkopt
                       ! Formal parameter of unique index of trkopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentrations

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      integer, intent(in) :: nni
                       ! Number of categories of ice concentrations

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer gwmopt   ! Option for gravity wave mode integration
      integer advopt   ! Option for advection scheme
      integer smtopt   ! Option for numerical smoothing
      integer cphopt   ! Option for cloud micro physics
      integer qcgopt   ! Option for charging distribution
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing

      integer nigrp1   ! nigrp + 1
      integer njgrp1   ! njgrp + 1

      integer nb       ! Number of exchanging variables

      integer igc      ! Array index in group domain in x direction
      integer jgc      ! Array index in group domain in y direction

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

      integer siz      ! Communication buffer size

! Internal private variables

      integer ijpe     ! Array index in sub domain

      integer ijsc     ! Index of serial number table

      integer igc_sub  ! Substitute for igc
      integer jgc_sub  ! Substitute for jgc

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)
      call getiname(fpgwmopt,gwmopt)
      call getiname(fpadvopt,advopt)
      call getiname(fpsmtopt,smtopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fpqcgopt,qcgopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)

! -----

! Set the common used variables.

      nigrp1=nigrp+1
      njgrp1=njgrp+1

! -----

! Count the maximum number of sending and receiving variables.

      if(mod(wbc,10).eq.6.or.mod(ebc,10).eq.6.or.                       &
     &   mod(sbc,10).eq.6.or.mod(nbc,10).eq.6.or.advopt.ge.2.or.        &
     &   tubopt.ge.1.or.mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

        nb=6

      else

        if(gwmopt.eq.0) then
          nb=2
        else
          nb=1
        end if

      end if

      if(abs(cphopt).lt.10) then

        if(abs(cphopt).ge.1) then
          nb=nb+nqw
        end if

        if(abs(cphopt).ge.2) then

          if(abs(cphopt).eq.2) then
            nb=nb+nqi+1
          else
            nb=nb+nqi+nni
          end if

        end if

        if(abs(cphopt).eq.4) then
          nb=nb+nnw
        end if

        if(cphopt.lt.0) then

          if(qcgopt.eq.2) then
            nb=nb+nqw
          end if

          nb=nb+nqi

        end if

      else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

        if(abs(cphopt).ge.11) then
          nb=nb+nqw+nnw
        end if

        if(abs(cphopt).eq.12) then
          nb=nb+nqi+nni
        end if

      end if

      if(aslopt.ge.1) then
        nb=nb+nqa(0)
      end if

      if(trkopt.ge.1) then
        nb=nb+1
      end if

      if(tubopt.ge.2) then
        nb=nb+1
      end if

! -----

! Set the communication buffer size.

      if(npe.eq.1) then

        siz=1

      else

        siz=max(2*(ni+2)*nk*nb,2*(nj+2)*nk*nb)

      end if

! -----

!! Allocate the communication buffer.

! Perform allocate.

      stat=0

      allocate(idxbuf(1:3,0:npe-1),mxnbuf(0:npe-1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(sbuf(1:siz),rbuf(1:siz),stat=cstat)

      stat=stat+abs(cstat)

      allocate(chrgrp(1:nigrp,1:njgrp),stat=cstat)

      stat=stat+abs(cstat)

      allocate(grpxy(0:nigrp+1,0:njgrp+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(xgrp(0:nsrl-1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ygrp(0:nsrl-1),stat=cstat)

      stat=stat+abs(cstat)

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('allocbuf',8,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('allocbuf',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

!! Initialize the table and set boundary conditions.

! Initiaze the character table.

      do jgc=1,njgrp
      do igc=1,nigrp

        write(chrgrp(igc,jgc)(1:1),'(a1)') '-'

      end do
      end do

! -----

! Initialize the other table and set boundary conditions.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(ijpe)

      do ijpe=0,npe-1
        idxbuf(1,ijpe)=0
        idxbuf(2,ijpe)=0
        idxbuf(3,ijpe)=0

        mxnbuf(ijpe)=0.e0

      end do

!$omp end do

!$omp do schedule(runtime) private(ijpe)

      do ijpe=1,siz
        sbuf(ijpe)=0.e0
        rbuf(ijpe)=0.e0
      end do

!$omp end do

!$omp do schedule(runtime) private(igc_sub,jgc_sub)

      do jgc_sub=0,njgrp+1
      do igc_sub=0,nigrp+1
        grpxy(igc_sub,jgc_sub)=-1
      end do
      end do

!$omp end do

      if(abs(wbc).eq.1.and.abs(ebc).eq.1) then

!$omp do schedule(runtime) private(jgc_sub)

        do jgc_sub=0,njgrp+1
          grpxy(0,jgc_sub)=nsrl
          grpxy(nigrp1,jgc_sub)=nsrl
        end do

!$omp end do

      end if

      if(abs(sbc).eq.1.and.abs(nbc).eq.1) then

!$omp do schedule(runtime) private(igc_sub)

        do igc_sub=0,nigrp+1
          grpxy(igc_sub,0)=nsrl
          grpxy(igc_sub,njgrp1)=nsrl
        end do

!$omp end do

      end if

!$omp do schedule(runtime) private(ijsc)

      do ijsc=0,nsrl-1
        xgrp(ijsc)=-1
        ygrp(ijsc)=-1
      end do

!$omp end do

!$omp end parallel

! -----

!! -----

      end subroutine s_allocbuf

!-----7--------------------------------------------------------------7--

      end module m_allocbuf
