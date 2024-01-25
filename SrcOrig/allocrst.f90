!***********************************************************************
      module m_allocrst
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/07/01
!     Modification: 2005/08/05, 2005/09/30, 2005/10/05, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/07/21, 2007/01/05,
!                   2007/01/20, 2007/05/07, 2007/07/30, 2007/11/26,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2008/12/11,
!                   2009/01/05, 2009/01/30, 2009/02/27, 2011/08/18,
!                   2011/09/22, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the array for rstruct.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commpi
      use m_comrst
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_getiname
      use m_inichar
      use m_setcst2d
      use m_setcst3d
      use m_setcst4d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: allocrst, s_allocrst

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allocrst

        module procedure s_allocrst

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_allocrst(fpdmpvar,fpsavmem,fpwbc,fpebc,fpsbc,fpnbc,  &
     &                      fpsfcopt,fpadvopt,fpcphopt,fpqcgopt,        &
     &                      fpaslopt,fptrkopt,fptubopt,fpdiaopt,        &
     &                      ni,nj,nk,nqw,nnw,nqi,nni,nqa,nund,          &
     &                      ni_rst,nj_rst)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdmpvar
                       ! Formal parameter of unique index of dmpvar

      integer, intent(in) :: fpsavmem
                       ! Formal parameter of unique index of savmem

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

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

      integer, intent(in) :: fpdiaopt
                       ! Formal parameter of unique index of diaopt

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

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

! Output variables

      integer, intent(out) :: ni_rst
                       ! Restructed files dimension in x direction

      integer, intent(out) :: nj_rst
                       ! Restructed files dimension in y direction

! Internal shared variables

      character(len=108) dmpvar
                       ! Control flag of dumped variables

      integer savmem   ! Option for memory saving

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer sfcopt   ! Option for surface physics
      integer advopt   ! Option for advection scheme
      integer cphopt   ! Option for cloud micro physics
      integer qcgopt   ! Option for charging distribution
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing
      integer diaopt   ! Option for diabatic calculation

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(dmpvar)

! -----

! Get the required namelist variables.

      call getcname(fpdmpvar,dmpvar)
      call getiname(fpsavmem,savmem)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)
      call getiname(fpsfcopt,sfcopt)
      call getiname(fpadvopt,advopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fpqcgopt,qcgopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)
      call getiname(fpdiaopt,diaopt)

! -----

! Get the model dimension of rstruct.

      ni_rst=((ni-3)*nisub)/nisub_rst+3
      nj_rst=((nj-3)*njsub)/njsub_rst+3

! -----

!! Allocate the array for rstruct.

! Perform allocate.

      stat=0

      allocate(ubr(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(vbr(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(pbr(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ptbr(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qvbr(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      if(savmem.eq.0.or.advopt.le.3) then

        allocate(u(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(v(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(w(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pp(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptp(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(u(0:0,0:0,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(v(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(w(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pp(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptp(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qv(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      allocate(up(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(vp(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(wp(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ppp(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ptpp(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qvp(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      if(savmem.eq.0.or.                                                &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        allocate(ucpx(1:nj,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ucpy(1:ni,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vcpx(1:nj,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vcpy(1:ni,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(wcpx(1:nj,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(wcpy(1:ni,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pcpx(1:nj,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pcpy(1:ni,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptcpx(1:nj,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptcpy(1:ni,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvcpx(1:nj,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvcpy(1:ni,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(ucpx(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ucpy(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vcpx(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vcpy(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(wcpx(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(wcpy(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pcpx(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pcpy(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptcpx(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptcpy(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvcpx(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvcpy(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.abs(cphopt).ge.1)) then

        allocate(qwtr(0:ni+1,0:nj+1,1:nk,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qwtr(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.                                                &
     &  (advopt.le.3.and.(abs(cphopt).ge.4.and.abs(cphopt).lt.20))) then

        allocate(nwtr(0:ni+1,0:nj+1,1:nk,1:nnw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(nwtr(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.abs(cphopt).ge.1) then

        allocate(qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qwtrp(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(abs(cphopt).ge.4.and.abs(cphopt).lt.20)) then

        allocate(nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(nwtrp(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.((abs(cphopt).ge.1.and.abs(cphopt).lt.20).and.  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(qwcpx(1:nj,1:nk,1:2,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qwcpy(1:ni,1:nk,1:2,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qwcpx(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qwcpy(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.((abs(cphopt).ge.4.and.abs(cphopt).lt.20).and.  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(nwcpx(1:nj,1:nk,1:2,1:nnw),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nwcpy(1:ni,1:nk,1:2,1:nnw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(nwcpx(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nwcpy(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.                               &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1))) then

        allocate(qice(0:ni+1,0:nj+1,1:nk,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qice(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.(abs(cphopt).ge.2.and.         &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11))) then

        if(abs(cphopt).eq.2) then

          allocate(nice(0:ni+1,0:nj+1,1:nk,1:1),stat=cstat)

          stat=stat+abs(cstat)

        else

          allocate(nice(0:ni+1,0:nj+1,1:nk,1:nni),stat=cstat)

          stat=stat+abs(cstat)

        end if

      else

        allocate(nice(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.                                                &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        allocate(qicep(0:ni+1,0:nj+1,1:nk,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qicep(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(abs(cphopt).ge.2.and.                          &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11)) then

        if(abs(cphopt).eq.2) then

          allocate(nicep(0:ni+1,0:nj+1,1:nk,1:1),stat=cstat)

          stat=stat+abs(cstat)

        else

          allocate(nicep(0:ni+1,0:nj+1,1:nk,1:nni),stat=cstat)

          stat=stat+abs(cstat)

        end if

      else

        allocate(nicep(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.((abs(cphopt).ge.2.and.                         &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(qicpx(1:nj,1:nk,1:2,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qicpy(1:ni,1:nk,1:2,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

        if(abs(cphopt).eq.2) then

          allocate(nicpx(1:nj,1:nk,1:2,1:1),stat=cstat)

          stat=stat+abs(cstat)

          allocate(nicpy(1:ni,1:nk,1:2,1:1),stat=cstat)

          stat=stat+abs(cstat)

        else

          allocate(nicpx(1:nj,1:nk,1:2,1:nni),stat=cstat)

          stat=stat+abs(cstat)

          allocate(nicpy(1:ni,1:nk,1:2,1:nni),stat=cstat)

          stat=stat+abs(cstat)

        end if

      else

        allocate(qicpx(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qicpy(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nicpx(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nicpy(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.abs(cphopt).ge.1) then

        if(abs(cphopt).lt.10) then

          allocate(prwtr(0:ni+1,0:nj+1,1:2,1:nqw),stat=cstat)

          stat=stat+abs(cstat)

        else

          allocate(prwtr(0:ni+1,0:nj+1,1:2,1:1),stat=cstat)

          stat=stat+abs(cstat)

        end if

      else

        allocate(prwtr(0:0,0:0,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.                                                &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        if(abs(cphopt).lt.10) then

          allocate(price(0:ni+1,0:nj+1,1:2,1:nqi),stat=cstat)

          stat=stat+abs(cstat)

        else

          allocate(price(0:ni+1,0:nj+1,1:2,1:1),stat=cstat)

          stat=stat+abs(cstat)

        end if

      else

        allocate(price(0:0,0:0,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.                                                &
     &  (advopt.le.3.and.cphopt.lt.0.and.qcgopt.eq.2)) then

        allocate(qcwtr(0:ni+1,0:nj+1,1:nk,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qcwtr(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.cphopt.lt.0)) then

        allocate(qcice(0:ni+1,0:nj+1,1:nk,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qcice(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(cphopt.lt.0.and.qcgopt.eq.2)) then

        allocate(qcwtrp(0:ni+1,0:nj+1,1:nk,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qcwtrp(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.cphopt.lt.0) then

        allocate(qcicep(0:ni+1,0:nj+1,1:nk,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qcicep(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(cphopt.lt.0.and.qcgopt.eq.2.and.               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(qcwcpx(1:nj,1:nk,1:2,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qcwcpy(1:ni,1:nk,1:2,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qcwcpx(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qcwcpy(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(cphopt.lt.0.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(qcicpx(1:nj,1:nk,1:2,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qcicpy(1:ni,1:nk,1:2,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qcicpx(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qcicpy(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.aslopt.ge.1)) then

        allocate(qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0)),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qasl(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.aslopt.ge.1) then

        allocate(qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0)),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qaslp(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(aslopt.ge.1.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(qacpx(1:nj,1:nk,1:2,1:nqa(0)),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qacpy(1:ni,1:nk,1:2,1:nqa(0)),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qacpx(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qacpy(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.trkopt.ge.1)) then

        allocate(qt(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qt(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.trkopt.ge.1) then

        allocate(qtp(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qtp(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(trkopt.ge.1.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(qtcpx(1:nj,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qtcpy(1:ni,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qtcpx(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qtcpy(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.tubopt.ge.2)) then

        allocate(tke(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(tke(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.tubopt.ge.2) then

        allocate(tkep(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(tkep(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(tubopt.ge.2.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(tkecpx(1:nj,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(tkecpy(1:ni,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(tkecpx(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(tkecpy(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(tubopt.ge.2.and.                               &
     &  (dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-'))) then

        allocate(maxvl(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(maxvl(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.diaopt.eq.1) then

        allocate(pdia(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(pdia(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.sfcopt.ge.1)) then

        allocate(tund(0:ni+1,0:nj+1,1:nund),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(tund(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.sfcopt.ge.1) then

        allocate(tundp(0:ni+1,0:nj+1,1:nund),stat=cstat)

        stat=stat+abs(cstat)

        allocate(z0m(0:ni+1,0:nj+1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(z0h(0:ni+1,0:nj+1),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(tundp(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(z0m(0:0,0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(z0h(0:0,0:0),stat=cstat)

        stat=stat+abs(cstat)

      end if

      allocate(ubr_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(vbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(pbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ptbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qvbr_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      if(savmem.eq.0.or.advopt.le.3) then

        allocate(u_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(v_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(w_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pp_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptp_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qv_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(u_rst(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(v_rst(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(w_rst(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pp_rst(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptp_rst(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qv_rst(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      allocate(up_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(vp_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(wp_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ppp_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ptpp_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qvp_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      if(savmem.eq.0.or.                                                &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        allocate(ucpx_rst(1:nj_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ucpy_rst(1:ni_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vcpx_rst(1:nj_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vcpy_rst(1:ni_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(wcpx_rst(1:nj_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(wcpy_rst(1:ni_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pcpx_rst(1:nj_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pcpy_rst(1:ni_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptcpx_rst(1:nj_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptcpy_rst(1:ni_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvcpx_rst(1:nj_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvcpy_rst(1:ni_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(ucpx_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ucpy_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vcpx_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vcpy_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(wcpx_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(wcpy_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pcpx_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pcpy_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptcpx_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptcpy_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvcpx_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvcpy_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.abs(cphopt).ge.1)) then

        allocate(qwtr_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qwtr_rst(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.                                                &
     &  (advopt.le.3.and.(abs(cphopt).ge.4.and.abs(cphopt).lt.20))) then

        allocate(nwtr_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nnw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(nwtr_rst(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.abs(cphopt).ge.1) then

        allocate(qwtrp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qwtrp_rst(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(abs(cphopt).ge.4.and.abs(cphopt).lt.20)) then

        allocate(nwtrp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nnw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(nwtrp_rst(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.((abs(cphopt).ge.1.and.abs(cphopt).lt.20).and.  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(qwcpx_rst(1:nj_rst,1:nk,1:2,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qwcpy_rst(1:ni_rst,1:nk,1:2,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qwcpx_rst(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qwcpy_rst(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.((abs(cphopt).ge.4.and.abs(cphopt).lt.20).and.  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(nwcpx_rst(1:nj_rst,1:nk,1:2,1:nnw),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nwcpy_rst(1:ni_rst,1:nk,1:2,1:nnw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(nwcpx_rst(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nwcpy_rst(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.                               &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1))) then

        allocate(qice_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qice_rst(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.(abs(cphopt).ge.2.and.         &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11))) then

        if(abs(cphopt).eq.2) then

          allocate(nice_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:1),            &
     &             stat=cstat)

          stat=stat+abs(cstat)

        else

          allocate(nice_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nni),          &
     &             stat=cstat)

          stat=stat+abs(cstat)

        end if

      else

        allocate(nice_rst(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.                                                &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        allocate(qicep_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qicep_rst(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(abs(cphopt).ge.2.and.                          &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11)) then

        if(abs(cphopt).eq.2) then

          allocate(nicep_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:1),           &
     &             stat=cstat)

          stat=stat+abs(cstat)

        else

          allocate(nicep_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nni),         &
     &             stat=cstat)

          stat=stat+abs(cstat)

        end if

      else

        allocate(nicep_rst(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.((abs(cphopt).ge.2.and.                         &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(qicpx_rst(1:nj_rst,1:nk,1:2,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qicpy_rst(1:ni_rst,1:nk,1:2,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

        if(abs(cphopt).eq.2) then

          allocate(nicpx_rst(1:nj_rst,1:nk,1:2,1:1),stat=cstat)

          stat=stat+abs(cstat)

          allocate(nicpy_rst(1:ni_rst,1:nk,1:2,1:1),stat=cstat)

          stat=stat+abs(cstat)

        else

          allocate(nicpx_rst(1:nj_rst,1:nk,1:2,1:nni),stat=cstat)

          stat=stat+abs(cstat)

          allocate(nicpy_rst(1:ni_rst,1:nk,1:2,1:nni),stat=cstat)

          stat=stat+abs(cstat)

        end if

      else

        allocate(qicpx_rst(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qicpy_rst(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nicpx_rst(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nicpy_rst(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.abs(cphopt).ge.1) then

        if(abs(cphopt).lt.10) then

          allocate(prwtr_rst(0:ni_rst+1,0:nj_rst+1,1:2,1:nqw),          &
     &             stat=cstat)

          stat=stat+abs(cstat)

        else

          allocate(prwtr_rst(0:ni_rst+1,0:nj_rst+1,1:2,1:1),            &
     &             stat=cstat)

          stat=stat+abs(cstat)

        end if

      else

        allocate(prwtr_rst(0:0,0:0,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.                                                &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        if(abs(cphopt).lt.10) then

          allocate(price_rst(0:ni_rst+1,0:nj_rst+1,1:2,1:nqi),          &
     &             stat=cstat)

          stat=stat+abs(cstat)

        else

          allocate(price_rst(0:ni_rst+1,0:nj_rst+1,1:2,1:1),            &
     &             stat=cstat)

          stat=stat+abs(cstat)

        end if

      else

        allocate(price_rst(0:0,0:0,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.                                                &
     &  (advopt.le.3.and.cphopt.lt.0.and.qcgopt.eq.2)) then

        allocate(qcwtr_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw),           &
     &           stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qcwtr_rst(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.cphopt.lt.0)) then

        allocate(qcice_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi),           &
     &           stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qcice_rst(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(cphopt.lt.0.and.qcgopt.eq.2)) then

        allocate(qcwtrp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqw),          &
     &           stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qcwtrp_rst(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.cphopt.lt.0) then

        allocate(qcicep_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqi),          &
     &           stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qcicep_rst(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(cphopt.lt.0.and.qcgopt.eq.2.and.               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(qcwcpx_rst(1:nj_rst,1:nk,1:2,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qcwcpy_rst(1:ni_rst,1:nk,1:2,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qcwcpx_rst(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qcwcpy_rst(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(cphopt.lt.0.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(qcicpx_rst(1:nj_rst,1:nk,1:2,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qcicpy_rst(1:ni_rst,1:nk,1:2,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qcicpx_rst(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qcicpy_rst(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.aslopt.ge.1)) then

        allocate(qasl_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqa(0)),         &
     &           stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qasl_rst(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.aslopt.ge.1) then

        allocate(qaslp_rst(0:ni_rst+1,0:nj_rst+1,1:nk,1:nqa(0)),        &
     &           stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qaslp_rst(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(aslopt.ge.1.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(qacpx_rst(1:nj_rst,1:nk,1:2,1:nqa(0)),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qacpy_rst(1:ni_rst,1:nk,1:2,1:nqa(0)),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qacpx_rst(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qacpy_rst(1:1,1:1,1:2,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.trkopt.ge.1)) then

        allocate(qt_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qt_rst(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.trkopt.ge.1) then

        allocate(qtp_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qtp_rst(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(trkopt.ge.1.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(qtcpx_rst(1:nj_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qtcpy_rst(1:ni_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qtcpx_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qtcpy_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.tubopt.ge.2)) then

        allocate(tke_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(tke_rst(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.tubopt.ge.2) then

        allocate(tkep_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(tkep_rst(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(tubopt.ge.2.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        allocate(tkecpx_rst(1:nj_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(tkecpy_rst(1:ni_rst,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(tkecpx_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(tkecpy_rst(1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(tubopt.ge.2.and.                               &
     &  (dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-'))) then

        allocate(maxvl_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(maxvl_rst(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.diaopt.eq.1) then

        allocate(pdia_rst(0:ni_rst+1,0:nj_rst+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(pdia_rst(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.sfcopt.ge.1)) then

        allocate(tund_rst(0:ni_rst+1,0:nj_rst+1,1:nund),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(tund_rst(0:0,0:0,1:nund),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.sfcopt.ge.1) then

        allocate(tundp_rst(0:ni_rst+1,0:nj_rst+1,1:nund),stat=cstat)

        stat=stat+abs(cstat)

        allocate(z0m_rst(0:ni_rst+1,0:nj_rst+1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(z0h_rst(0:ni_rst+1,0:nj_rst+1),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(tundp_rst(0:0,0:0,1:nund),stat=cstat)

        stat=stat+abs(cstat)

        allocate(z0m_rst(0:0,0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(z0h_rst(0:0,0:0),stat=cstat)

        stat=stat+abs(cstat)

      end if

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('allocrst',8,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('allocrst',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

! Fill in all array with 0.

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ubr)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vbr)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,pbr)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ptbr)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qvbr)

      if(savmem.eq.0.or.advopt.le.3) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,u)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,v)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,w)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,pp)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ptp)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qv)

      end if

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,up)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vp)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,wp)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ppp)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ptpp)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qvp)

      if(savmem.eq.0.or.                                                &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        call setcst3d(1,nj,1,nk,1,2,0.e0,ucpx)
        call setcst3d(1,ni,1,nk,1,2,0.e0,ucpy)

        call setcst3d(1,nj,1,nk,1,2,0.e0,vcpx)
        call setcst3d(1,ni,1,nk,1,2,0.e0,vcpy)

        call setcst3d(1,nj,1,nk,1,2,0.e0,wcpx)
        call setcst3d(1,ni,1,nk,1,2,0.e0,wcpy)

        call setcst3d(1,nj,1,nk,1,2,0.e0,pcpx)
        call setcst3d(1,ni,1,nk,1,2,0.e0,pcpy)

        call setcst3d(1,nj,1,nk,1,2,0.e0,ptcpx)
        call setcst3d(1,ni,1,nk,1,2,0.e0,ptcpy)

        call setcst3d(1,nj,1,nk,1,2,0.e0,qvcpx)
        call setcst3d(1,ni,1,nk,1,2,0.e0,qvcpy)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.abs(cphopt).ge.1)) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qwtr)

      end if

      if(savmem.eq.0.or.                                                &
     &  (advopt.le.3.and.(abs(cphopt).ge.4.and.abs(cphopt).lt.20))) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nnw,0.e0,nwtr)

      end if

      if(savmem.eq.0.or.abs(cphopt).ge.1) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qwtrp)

      end if

      if(savmem.eq.0.or.(abs(cphopt).ge.4.and.abs(cphopt).lt.20)) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nnw,0.e0,nwtrp)

      end if

      if(savmem.eq.0.or.((abs(cphopt).ge.1.and.abs(cphopt).lt.20).and.  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call setcst4d(1,nj,1,nk,1,2,1,nqw,0.e0,qwcpx)
        call setcst4d(1,ni,1,nk,1,2,1,nqw,0.e0,qwcpy)

      end if

      if(savmem.eq.0.or.((abs(cphopt).ge.4.and.abs(cphopt).lt.20).and.  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call setcst4d(1,nj,1,nk,1,2,1,nnw,0.e0,nwcpx)
        call setcst4d(1,ni,1,nk,1,2,1,nnw,0.e0,nwcpy)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.                               &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1))) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qice)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.(abs(cphopt).ge.2.and.         &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11))) then

        if(abs(cphopt).eq.2) then

          call setcst4d(0,ni+1,0,nj+1,1,nk,1,1,0.e0,nice)

        else

          call setcst4d(0,ni+1,0,nj+1,1,nk,1,nni,0.e0,nice)

        end if

      end if

      if(savmem.eq.0.or.                                                &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qicep)

      end if

      if(savmem.eq.0.or.(abs(cphopt).ge.2.and.                          &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11)) then

        if(abs(cphopt).eq.2) then

          call setcst4d(0,ni+1,0,nj+1,1,nk,1,1,0.e0,nicep)

        else

          call setcst4d(0,ni+1,0,nj+1,1,nk,1,nni,0.e0,nicep)

        end if

      end if

      if(savmem.eq.0.or.((abs(cphopt).ge.2.and.                         &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call setcst4d(1,nj,1,nk,1,2,1,nqi,0.e0,qicpx)
        call setcst4d(1,ni,1,nk,1,2,1,nqi,0.e0,qicpy)

        if(abs(cphopt).eq.2) then

          call setcst4d(1,nj,1,nk,1,2,1,1,0.e0,nicpx)
          call setcst4d(1,ni,1,nk,1,2,1,1,0.e0,nicpy)

        else

          call setcst4d(1,nj,1,nk,1,2,1,nni,0.e0,nicpx)
          call setcst4d(1,ni,1,nk,1,2,1,nni,0.e0,nicpy)

        end if

      end if

      if(savmem.eq.0.or.abs(cphopt).ge.1) then

        if(abs(cphopt).lt.10) then

          call setcst4d(0,ni+1,0,nj+1,1,2,1,nqw,0.e0,prwtr)

        else

          call setcst4d(0,ni+1,0,nj+1,1,2,1,1,0.e0,prwtr)

        end if

      end if

      if(savmem.eq.0.or.                                                &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        if(abs(cphopt).lt.10) then

          call setcst4d(0,ni+1,0,nj+1,1,2,1,nqi,0.e0,price)

        else

          call setcst4d(0,ni+1,0,nj+1,1,2,1,1,0.e0,price)

        end if

      end if

      if(savmem.eq.0.or.                                                &
     &  (advopt.le.3.and.cphopt.lt.0.and.qcgopt.eq.2)) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qcwtr)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.cphopt.lt.0)) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qcice)

      end if

      if(savmem.eq.0.or.(cphopt.lt.0.and.qcgopt.eq.2)) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qcwtrp)

      end if

      if(savmem.eq.0.or.cphopt.lt.0) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qcicep)

      end if

      if(savmem.eq.0.or.(cphopt.lt.0.and.qcgopt.eq.2.and.               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call setcst4d(1,nj,1,nk,1,2,1,nqw,0.e0,qcwcpx)
        call setcst4d(1,ni,1,nk,1,2,1,nqw,0.e0,qcwcpy)

      end if

      if(savmem.eq.0.or.(cphopt.lt.0.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call setcst4d(1,nj,1,nk,1,2,1,nqi,0.e0,qcicpx)
        call setcst4d(1,ni,1,nk,1,2,1,nqi,0.e0,qcicpy)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.aslopt.ge.1)) then

        call s_setcst4d(0,ni+1,0,nj+1,1,nk,1,nqa(0),0.e0,qasl)

      end if

      if(savmem.eq.0.or.aslopt.ge.1) then

        call s_setcst4d(0,ni+1,0,nj+1,1,nk,1,nqa(0),0.e0,qaslp)

      end if

      if(savmem.eq.0.or.(aslopt.ge.1.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call s_setcst4d(1,nj,1,nk,1,2,1,nqa(0),0.e0,qacpx)
        call s_setcst4d(1,ni,1,nk,1,2,1,nqa(0),0.e0,qacpy)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.trkopt.ge.1)) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qt)

      end if

      if(savmem.eq.0.or.trkopt.ge.1) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qtp)

      end if

      if(savmem.eq.0.or.(trkopt.ge.1.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call setcst3d(1,nj,1,nk,1,2,0.e0,qtcpx)
        call setcst3d(1,ni,1,nk,1,2,0.e0,qtcpy)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.tubopt.ge.2)) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tke)

      end if

      if(savmem.eq.0.or.tubopt.ge.2) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tkep)

      end if

      if(savmem.eq.0.or.(tubopt.ge.2.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call setcst3d(1,nj,1,nk,1,2,0.e0,tkecpx)
        call setcst3d(1,ni,1,nk,1,2,0.e0,tkecpy)

      end if

      if(savmem.eq.0.or.(tubopt.ge.2.and.                               &
     &  (dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-'))) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,maxvl)

      end if

      if(savmem.eq.0.or.diaopt.eq.1) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,pdia)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.sfcopt.ge.1)) then

        call setcst3d(0,ni+1,0,nj+1,1,nund,0.e0,tund)

      end if

      if(savmem.eq.0.or.sfcopt.ge.1) then

        call setcst3d(0,ni+1,0,nj+1,1,nund,0.e0,tundp)

        call setcst2d(0,ni+1,0,nj+1,0.e0,z0m)
        call setcst2d(0,ni+1,0,nj+1,0.e0,z0h)

      end if

      call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,ubr_rst)
      call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,vbr_rst)
      call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,pbr_rst)
      call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,ptbr_rst)
      call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,qvbr_rst)

      if(savmem.eq.0.or.advopt.le.3) then

        call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,u_rst)
        call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,v_rst)
        call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,w_rst)
        call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,pp_rst)
        call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,ptp_rst)
        call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,qv_rst)

      end if

      call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,up_rst)
      call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,vp_rst)
      call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,wp_rst)
      call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,ppp_rst)
      call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,ptpp_rst)
      call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,qvp_rst)

      if(savmem.eq.0.or.                                                &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4)) then

        call setcst3d(1,nj_rst,1,nk,1,2,0.e0,ucpx_rst)
        call setcst3d(1,ni_rst,1,nk,1,2,0.e0,ucpy_rst)

        call setcst3d(1,nj_rst,1,nk,1,2,0.e0,vcpx_rst)
        call setcst3d(1,ni_rst,1,nk,1,2,0.e0,vcpy_rst)

        call setcst3d(1,nj_rst,1,nk,1,2,0.e0,wcpx_rst)
        call setcst3d(1,ni_rst,1,nk,1,2,0.e0,wcpy_rst)

        call setcst3d(1,nj_rst,1,nk,1,2,0.e0,pcpx_rst)
        call setcst3d(1,ni_rst,1,nk,1,2,0.e0,pcpy_rst)

        call setcst3d(1,nj_rst,1,nk,1,2,0.e0,ptcpx_rst)
        call setcst3d(1,ni_rst,1,nk,1,2,0.e0,ptcpy_rst)

        call setcst3d(1,nj_rst,1,nk,1,2,0.e0,qvcpx_rst)
        call setcst3d(1,ni_rst,1,nk,1,2,0.e0,qvcpy_rst)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.abs(cphopt).ge.1)) then

        call setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,nqw,0.e0,qwtr_rst)

      end if

      if(savmem.eq.0.or.                                                &
     &  (advopt.le.3.and.(abs(cphopt).ge.4.and.abs(cphopt).lt.20))) then

        call setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,nnw,0.e0,nwtr_rst)

      end if

      if(savmem.eq.0.or.abs(cphopt).ge.1) then

        call setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,nqw,0.e0,qwtrp_rst)

      end if

      if(savmem.eq.0.or.(abs(cphopt).ge.4.and.abs(cphopt).lt.20)) then

        call setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,nnw,0.e0,nwtrp_rst)

      end if

      if(savmem.eq.0.or.((abs(cphopt).ge.1.and.abs(cphopt).lt.20).and.  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call setcst4d(1,nj_rst,1,nk,1,2,1,nqw,0.e0,qwcpx_rst)
        call setcst4d(1,ni_rst,1,nk,1,2,1,nqw,0.e0,qwcpy_rst)

      end if

      if(savmem.eq.0.or.((abs(cphopt).ge.4.and.abs(cphopt).lt.20).and.  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call setcst4d(1,nj_rst,1,nk,1,2,1,nnw,0.e0,nwcpx_rst)
        call setcst4d(1,ni_rst,1,nk,1,2,1,nnw,0.e0,nwcpy_rst)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.                               &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1))) then

        call setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,nqi,0.e0,qice_rst)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.(abs(cphopt).ge.2.and.         &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11))) then

        if(abs(cphopt).eq.2) then

          call setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,1,0.e0,nice_rst)

        else

          call setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,nni,0.e0,nice_rst)

        end if

      end if

      if(savmem.eq.0.or.                                                &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        call setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,nqi,0.e0,qicep_rst)

      end if

      if(savmem.eq.0.or.(abs(cphopt).ge.2.and.                          &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11)) then

        if(abs(cphopt).eq.2) then

          call setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,1,0.e0,nicep_rst)

        else

          call setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,nni,0.e0,nicep_rst)

        end if

      end if

      if(savmem.eq.0.or.((abs(cphopt).ge.2.and.                         &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11).and.                  &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call setcst4d(1,nj_rst,1,nk,1,2,1,nqi,0.e0,qicpx_rst)
        call setcst4d(1,ni_rst,1,nk,1,2,1,nqi,0.e0,qicpy_rst)

        if(abs(cphopt).eq.2) then

          call setcst4d(1,nj_rst,1,nk,1,2,1,1,0.e0,nicpx_rst)
          call setcst4d(1,ni_rst,1,nk,1,2,1,1,0.e0,nicpy_rst)

        else

          call setcst4d(1,nj_rst,1,nk,1,2,1,nni,0.e0,nicpx_rst)
          call setcst4d(1,ni_rst,1,nk,1,2,1,nni,0.e0,nicpy_rst)

        end if

      end if

      if(savmem.eq.0.or.abs(cphopt).ge.1) then

        if(abs(cphopt).lt.10) then

          call setcst4d(0,ni_rst+1,0,nj_rst+1,1,2,1,nqw,0.e0,prwtr_rst)

        else

          call setcst4d(0,ni_rst+1,0,nj_rst+1,1,2,1,1,0.e0,prwtr_rst)

        end if

      end if

      if(savmem.eq.0.or.                                                &
     &  (abs(cphopt).ge.2.and.mod(abs(cphopt),10).ne.1)) then

        if(abs(cphopt).lt.10) then

          call setcst4d(0,ni_rst+1,0,nj_rst+1,1,2,1,nqi,0.e0,price_rst)

        else

          call setcst4d(0,ni_rst+1,0,nj_rst+1,1,2,1,1,0.e0,price_rst)

        end if

      end if

      if(savmem.eq.0.or.                                                &
     &  (advopt.le.3.and.cphopt.lt.0.and.qcgopt.eq.2)) then

        call setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,nqw,0.e0,qcwtr_rst)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.cphopt.lt.0)) then

        call setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,nqi,0.e0,qcice_rst)

      end if

      if(savmem.eq.0.or.(cphopt.lt.0.and.qcgopt.eq.2)) then

        call setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,nqw,0.e0,qcwtrp_rst)

      end if

      if(savmem.eq.0.or.cphopt.lt.0) then

        call setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,nqi,0.e0,qcicep_rst)

      end if

      if(savmem.eq.0.or.(cphopt.lt.0.and.qcgopt.eq.2.and.               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call setcst4d(1,nj_rst,1,nk,1,2,1,nqw,0.e0,qcwcpx_rst)
        call setcst4d(1,ni_rst,1,nk,1,2,1,nqw,0.e0,qcwcpy_rst)

      end if

      if(savmem.eq.0.or.(cphopt.lt.0.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call setcst4d(1,nj_rst,1,nk,1,2,1,nqi,0.e0,qcicpx_rst)
        call setcst4d(1,ni_rst,1,nk,1,2,1,nqi,0.e0,qcicpy_rst)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.aslopt.ge.1)) then

        call s_setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,nqa(0),0.e0,       &
     &                  qasl_rst)

      end if

      if(savmem.eq.0.or.aslopt.ge.1) then

        call s_setcst4d(0,ni_rst+1,0,nj_rst+1,1,nk,1,nqa(0),0.e0,       &
     &                  qaslp_rst)

      end if

      if(savmem.eq.0.or.(aslopt.ge.1.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call s_setcst4d(1,nj_rst,1,nk,1,2,1,nqa(0),0.e0,qacpx_rst)
        call s_setcst4d(1,ni_rst,1,nk,1,2,1,nqa(0),0.e0,qacpy_rst)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.trkopt.ge.1)) then

        call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,qt_rst)

      end if

      if(savmem.eq.0.or.trkopt.ge.1) then

        call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,qtp_rst)

      end if

      if(savmem.eq.0.or.(trkopt.ge.1.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call setcst3d(1,nj_rst,1,nk,1,2,0.e0,qtcpx_rst)
        call setcst3d(1,ni_rst,1,nk,1,2,0.e0,qtcpy_rst)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.tubopt.ge.2)) then

        call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,tke_rst)

      end if

      if(savmem.eq.0.or.tubopt.ge.2) then

        call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,tkep_rst)

      end if

      if(savmem.eq.0.or.(tubopt.ge.2.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call setcst3d(1,nj_rst,1,nk,1,2,0.e0,tkecpx_rst)
        call setcst3d(1,ni_rst,1,nk,1,2,0.e0,tkecpy_rst)

      end if

      if(savmem.eq.0.or.(tubopt.ge.2.and.                               &
     &  (dmpvar(12:12).eq.'+'.or.dmpvar(12:12).eq.'-'))) then

        call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,maxvl_rst)

      end if

      if(savmem.eq.0.or.diaopt.eq.1) then

        call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nk,0.e0,pdia_rst)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.sfcopt.ge.1)) then

        call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nund,0.e0,tund_rst)

      end if

      if(savmem.eq.0.or.sfcopt.ge.1) then

        call setcst3d(0,ni_rst+1,0,nj_rst+1,1,nund,0.e0,tundp_rst)

        call setcst2d(0,ni_rst+1,0,nj_rst+1,0.e0,z0m_rst)
        call setcst2d(0,ni_rst+1,0,nj_rst+1,0.e0,z0h_rst)

      end if

! -----

      end subroutine s_allocrst

!-----7--------------------------------------------------------------7--

      end module m_allocrst
