!***********************************************************************
      module m_allocslv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/07/15, 2003/08/08, 2003/11/05, 2003/11/28,
!                   2003/12/12, 2004/01/09, 2004/02/01, 2004/04/15,
!                   2004/05/07, 2004/05/31, 2004/06/10, 2004/08/01,
!                   2004/08/20, 2005/01/14, 2005/02/10, 2005/08/05,
!                   2005/10/05, 2006/01/10, 2006/02/13, 2006/04/03,
!                   2006/06/21, 2006/07/21, 2006/09/30, 2006/11/06,
!                   2006/11/27, 2007/01/05, 2007/01/20, 2007/01/31,
!                   2007/03/10, 2007/05/07, 2007/05/21, 2007/07/30,
!                   2007/10/19, 2007/11/26, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2008/12/11, 2009/01/05, 2009/01/30,
!                   2009/02/27, 2011/08/18, 2011/09/22, 2011/11/10,
!                   2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the array for solver.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commpi
      use m_comslv
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_getiname
      use m_inichar
      use m_setcst1d
      use m_setcst2d
      use m_setcst3d
      use m_setcst4d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: allocslv, s_allocslv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allocslv

        module procedure s_allocslv

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
      subroutine s_allocslv(fpdmpvar,fpsavmem,                          &
     &                      fpwbc,fpebc,fpsbc,fpnbc,fpnggopt,           &
     &                      fpexbopt,fplspopt,fpvspopt,fpngropt,        &
     &                      fpiniopt,fpsfcopt,fpadvopt,fpcphopt,        &
     &                      fpqcgopt,fpaslopt,fptrkopt,fptubopt,        &
     &                      ni,nj,nk,nqw,nnw,nqi,nni,km,nqa,nund,nlev,  &
     &                      nid_rdr,njd_rdr,nkd_rdr,km_rdr)
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

      integer, intent(in) :: fpnggopt
                       ! Formal parameter of unique index of nggopt

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: fplspopt
                       ! Formal parameter of unique index of lspopt

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: fpngropt
                       ! Formal parameter of unique index of ngropt

      integer, intent(in) :: fpiniopt
                       ! Formal parameter of unique index of iniopt

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

      integer, intent(in) :: km
                       ! Dimension of max(nk, nqw, nqi)

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

      integer, intent(in) :: nlev
                       ! Horizontally averaged vertical dimension

      integer, intent(in) :: nid_rdr
                       ! Radar data dimension in x direction

      integer, intent(in) :: njd_rdr
                       ! Radar data dimension in y direction

      integer, intent(in) :: nkd_rdr
                       ! Radar data dimension in z direction

      integer, intent(in) :: km_rdr
                       ! Dimension of max(nk, nkd_rdr)

! Internal shared variables

      character(len=108) dmpvar
                       ! Control flag of dumped variables

      integer savmem   ! Option for memory saving

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer nggopt   ! Option for analysis nudging to GPV
      integer exbopt   ! Option for external boundary forcing
      integer lspopt   ! Option for lateral sponge damping
      integer vspopt   ! Option for vertical sponge damping
      integer ngropt   ! Option for analysis nudging to radar data
      integer iniopt   ! Option for model initialization
      integer sfcopt   ! Option for surface physics
      integer advopt   ! Option for advection scheme
      integer cphopt   ! Option for cloud micro physics
      integer qcgopt   ! Option for charging distribution
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

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
      call getiname(fpnggopt,nggopt)
      call getiname(fpexbopt,exbopt)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getiname(fpngropt,ngropt)
      call getiname(fpiniopt,iniopt)
      call getiname(fpsfcopt,sfcopt)
      call getiname(fpadvopt,advopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fpqcgopt,qcgopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)

! -----

!! Allocate the array for solver.

! Perform allocate.

      stat=0

      allocate(zph(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(zsth(1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(lat(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(lon(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(j31(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(j32(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(jcb(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(jcb8u(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(jcb8v(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(jcb8w(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(mf(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(mf8u(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(mf8v(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(rmf(0:ni+1,0:nj+1,1:4),stat=cstat)

      stat=stat+abs(cstat)

      allocate(rmf8u(0:ni+1,0:nj+1,1:3),stat=cstat)

      stat=stat+abs(cstat)

      allocate(rmf8v(0:ni+1,0:nj+1,1:3),stat=cstat)

      stat=stat+abs(cstat)

      allocate(fc(0:ni+1,0:nj+1,1:2),stat=cstat)

      stat=stat+abs(cstat)

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

      allocate(rbr(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(rst(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(rst8u(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(rst8v(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(rst8w(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(rcsq(0:ni+1,0:nj+1,1:nk),stat=cstat)

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

        allocate(u(0:0,0:0,1:1),stat=cstat)

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

      allocate(pdia(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp1(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp2(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp3(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp4(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp5(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp6(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp7(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp8(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp9(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp10(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp11(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp12(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp13(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp14(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp15(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp16(0:ni+1,0:nj+1,1:km),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp17(0:ni+1,0:nj+1,1:km),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp18(0:ni+1,0:nj+1,1:km),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp19(0:ni+1,0:nj+1,1:km),stat=cstat)

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

      if(savmem.eq.0.or.lspopt.ge.1) then

        allocate(rbcx(1:ni),stat=cstat)

        stat=stat+abs(cstat)

        allocate(rbcy(1:nj),stat=cstat)

        stat=stat+abs(cstat)

        allocate(rbcxy(1:ni,1:nj),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(rbcx(1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(rbcy(1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(rbcxy(1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.vspopt.ge.1) then

        allocate(rbct(1:ni,1:nj,1:nk,1:2),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(rbct(1:1,1:1,1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(nggopt.eq.1.or.                                &
     &   exbopt.ge.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

        allocate(ugpv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(utd(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vgpv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vtd(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(wgpv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(wtd(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ppgpv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pptd(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptpgpv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptptd(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvgpv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvtd(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(ugpv(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(utd(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vgpv(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vtd(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(wgpv(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(wtd(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ppgpv(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pptd(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptpgpv(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ptptd(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvgpv(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvtd(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.ngropt.ge.1) then

        allocate(urdr(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vrdr(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(wrdr(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(urdr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vrdr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(wrdr(0:0,0:0,1:1),stat=cstat)

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

        allocate(qwtmp(0:ni+1,0:nj+1,1:nk,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qwtrp(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qwtmp(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(abs(cphopt).ge.2.and.abs(cphopt).lt.20)) then

        allocate(nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(nwtrp(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(abs(cphopt).ge.4.and.abs(cphopt).lt.20)) then

        allocate(nwtmp(0:ni+1,0:nj+1,1:nk,1:nnw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(nwtmp(0:0,0:0,1:1,1:1),stat=cstat)

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

      if(savmem.eq.0.or.((nggopt.eq.1.or.                               &
     &   exbopt.ge.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1).and.        &
     &  (abs(cphopt).ge.1.and.abs(cphopt).lt.10))) then

        allocate(qwgpv(0:ni+1,0:nj+1,1:nk,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qwtd(0:ni+1,0:nj+1,1:nk,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qwgpv(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qwtd(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.                                                &
     &  (ngropt.ge.1.and.(abs(cphopt).ge.1.and.abs(cphopt).lt.10))) then

        allocate(qwrdr(0:ni+1,0:nj+1,1:nk,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qwrdr(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.                                                &
     &  (ngropt.eq.1.and.(abs(cphopt).ge.1.and.abs(cphopt).lt.10))) then

        allocate(qwrtd(0:ni+1,0:nj+1,1:nk,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qwrtd(0:0,0:0,1:1,1:1),stat=cstat)

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

        allocate(qitmp(0:ni+1,0:nj+1,1:nk,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qicep(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qitmp(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(abs(cphopt).ge.2.and.                          &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11)) then

        allocate(nicep(0:ni+1,0:nj+1,1:nk,1:nni),stat=cstat)

        stat=stat+abs(cstat)

        if(abs(cphopt).eq.2) then

          allocate(nitmp(0:ni+1,0:nj+1,1:nk,1:1),stat=cstat)

          stat=stat+abs(cstat)

        else

          allocate(nitmp(0:ni+1,0:nj+1,1:nk,1:nni),stat=cstat)

          stat=stat+abs(cstat)

        end if

      else

        allocate(nicep(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nitmp(0:0,0:0,1:1,1:1),stat=cstat)

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

      if(savmem.eq.0.or.((nggopt.eq.1.or.                               &
     &   exbopt.ge.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1).and.        &
     &  (abs(cphopt).ge.2.and.abs(cphopt).lt.10))) then

        allocate(qigpv(0:ni+1,0:nj+1,1:nk,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qitd(0:ni+1,0:nj+1,1:nk,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qigpv(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qitd(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.                                                &
     &  (ngropt.ge.1.and.(abs(cphopt).ge.2.and.abs(cphopt).lt.10))) then

        allocate(qirdr(0:ni+1,0:nj+1,1:nk,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qirdr(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.                                                &
     &  (ngropt.eq.1.and.(abs(cphopt).ge.2.and.abs(cphopt).lt.10))) then

        allocate(qirtd(0:ni+1,0:nj+1,1:nk,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qirtd(0:0,0:0,1:1,1:1),stat=cstat)

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

        allocate(qcwtmp(0:ni+1,0:nj+1,1:nk,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qcwtrp(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qcwtmp(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.cphopt.lt.0) then

        allocate(qcicep(0:ni+1,0:nj+1,1:nk,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qcitmp(0:ni+1,0:nj+1,1:nk,1:nqi),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qcicep(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qcitmp(0:0,0:0,1:1,1:1),stat=cstat)

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

        allocate(qatmp(0:ni+1,0:nj+1,1:nk,1:nqa(0)),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qaslp(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qatmp(0:0,0:0,1:1,1:1),stat=cstat)

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

      if(savmem.eq.0.or.((nggopt.eq.1.or.exbopt.ge.1.or.                &
     &   mod(lspopt,10).eq.1.or.vspopt.eq.1).and.aslopt.ge.1)) then

        allocate(qagpv(0:ni+1,0:nj+1,1:nk,1:nqa(0)),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qatd(0:ni+1,0:nj+1,1:nk,1:nqa(0)),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qagpv(0:0,0:0,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qatd(0:0,0:0,1:1,1:1),stat=cstat)

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

        allocate(qttmp(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qtp(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qttmp(0:0,0:0,1:1),stat=cstat)

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

        allocate(tketmp(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(tkep(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(tketmp(0:0,0:0,1:1),stat=cstat)

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

      if(savmem.eq.0.or.(advopt.le.3.and.abs(cphopt).ge.1)) then

        allocate(qall(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qall(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.abs(cphopt).ge.1) then

        allocate(qallp(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(qallp(0:0,0:0,1:1),stat=cstat)

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

        allocate(tutmp(0:ni+1,0:nj+1,1:nund),stat=cstat)

        stat=stat+abs(cstat)

        allocate(land(0:ni+1,0:nj+1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(albe(0:ni+1,0:nj+1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(beta(0:ni+1,0:nj+1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(z0m(0:ni+1,0:nj+1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(z0h(0:ni+1,0:nj+1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cap(0:ni+1,0:nj+1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nuu(0:ni+1,0:nj+1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(kai(0:ni+1,0:nj+1),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(tundp(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(tutmp(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(land(0:0,0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(albe(0:0,0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(beta(0:0,0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(z0m(0:0,0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(z0h(0:0,0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cap(0:0,0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nuu(0:0,0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(kai(0:0,0:0),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.(sfcopt.eq.3.or.sfcopt.eq.13)) then

        allocate(sst(0:ni+1,0:nj+1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(sstd(0:ni+1,0:nj+1),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(sst(0:0,0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(sstd(0:0,0:0),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.iniopt.eq.1) then

        allocate(z1d(0:nlev),stat=cstat)

        stat=stat+abs(cstat)

        allocate(u1d(0:nlev),stat=cstat)

        stat=stat+abs(cstat)

        allocate(v1d(0:nlev),stat=cstat)

        stat=stat+abs(cstat)

        allocate(p1d(0:nlev),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pt1d(0:nlev),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qv1d(0:nlev),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ltmp1(1:nlev),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ltmp2(1:nlev),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ltmp3(1:nlev),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(z1d(0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(u1d(0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(v1d(0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(p1d(0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pt1d(0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qv1d(0:0),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ltmp1(1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ltmp2(1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ltmp3(1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

      if(savmem.eq.0.or.ngropt.eq.12) then

        allocate(lon_rdr(1:nid_rdr,1:njd_rdr),stat=cstat)

        stat=stat+abs(cstat)

        allocate(z_rdr(1:nid_rdr,1:njd_rdr,1:nkd_rdr),stat=cstat)

        stat=stat+abs(cstat)

        allocate(u_rdr(1:nid_rdr,1:njd_rdr,1:km_rdr),stat=cstat)

        stat=stat+abs(cstat)

        allocate(v_rdr(1:nid_rdr,1:njd_rdr,1:km_rdr),stat=cstat)

        stat=stat+abs(cstat)

        allocate(w_rdr(1:nid_rdr,1:njd_rdr,1:km_rdr),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qp_rdr(1:nid_rdr,1:njd_rdr,1:km_rdr),stat=cstat)

        stat=stat+abs(cstat)

        allocate(tmp1_rdr(1:nid_rdr,1:njd_rdr,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(lon_rdr(1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(z_rdr(1:1,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(u_rdr(1:1,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(v_rdr(1:1,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(w_rdr(1:1,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qp_rdr(1:1,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(tmp1_rdr(1:1,1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('allocslv',8,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('allocslv',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

!! Fill in all array with 0.

! For the integer variable.

      if(savmem.eq.0.or.sfcopt.ge.1) then

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j)

        do j=0,nj+1
        do i=0,ni+1
          land(i,j)=0
        end do
        end do

!$omp end do

!$omp end parallel

      end if

! -----

! For the real variables.

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,zph)

      call setcst1d(1,nk,0.e0,zsth)

      call setcst2d(0,ni+1,0,nj+1,0.e0,lat)
      call setcst2d(0,ni+1,0,nj+1,0.e0,lon)

      call setcst2d(0,ni+1,0,nj+1,0.e0,mf)
      call setcst2d(0,ni+1,0,nj+1,0.e0,mf8u)
      call setcst2d(0,ni+1,0,nj+1,0.e0,mf8v)

      call setcst3d(0,ni+1,0,nj+1,1,4,0.e0,rmf)
      call setcst3d(0,ni+1,0,nj+1,1,3,0.e0,rmf8u)
      call setcst3d(0,ni+1,0,nj+1,1,3,0.e0,rmf8v)

      call setcst3d(0,ni+1,0,nj+1,1,2,0.e0,fc)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,j31)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,j32)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,jcb)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,jcb8u)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,jcb8v)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,jcb8w)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ubr)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vbr)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,pbr)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ptbr)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qvbr)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,rbr)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,rst)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,rst8u)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,rst8v)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,rst8w)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,rcsq)

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

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,pdia)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp1)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp2)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp3)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp4)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp5)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp6)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp7)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp8)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp9)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp10)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp11)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp12)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp13)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp14)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp15)

      call setcst3d(0,ni+1,0,nj+1,1,km,0.e0,tmp16)
      call setcst3d(0,ni+1,0,nj+1,1,km,0.e0,tmp17)
      call setcst3d(0,ni+1,0,nj+1,1,km,0.e0,tmp18)
      call setcst3d(0,ni+1,0,nj+1,1,km,0.e0,tmp19)

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

      if(savmem.eq.0.or.lspopt.ge.1) then

        call setcst1d(1,ni,0.e0,rbcx)
        call setcst1d(1,nj,0.e0,rbcy)

        call setcst2d(1,ni,1,nj,0.e0,rbcxy)

      end if

      if(savmem.eq.0.or.vspopt.ge.1) then

        call setcst4d(1,ni,1,nj,1,nk,1,2,0.e0,rbct)

      end if

      if(savmem.eq.0.or.(nggopt.eq.1.or.                                &
     &   exbopt.ge.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1)) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ugpv)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,utd)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vgpv)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vtd)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,wgpv)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,wtd)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ppgpv)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,pptd)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ptpgpv)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ptptd)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qvgpv)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qvtd)

      end if

      if(savmem.eq.0.or.ngropt.ge.1) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,urdr)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vrdr)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,wrdr)

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
        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qwtmp)

      end if

      if(savmem.eq.0.or.(abs(cphopt).ge.2.and.abs(cphopt).lt.20)) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nnw,0.e0,nwtrp)

      end if

      if(savmem.eq.0.or.(abs(cphopt).ge.4.and.abs(cphopt).lt.20)) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nnw,0.e0,nwtmp)

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

      if(savmem.eq.0.or.((nggopt.eq.1.or.                               &
     &   exbopt.ge.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1).and.        &
     &  (abs(cphopt).ge.1.and.abs(cphopt).lt.10))) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qwgpv)
        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qwtd)

      end if

      if(savmem.eq.0.or.                                                &
     &  (ngropt.ge.1.and.(abs(cphopt).ge.1.and.abs(cphopt).lt.10))) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qwrdr)

      end if

      if(savmem.eq.0.or.                                                &
     &  (ngropt.eq.1.and.(abs(cphopt).ge.1.and.abs(cphopt).lt.10))) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qwrtd)

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
        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qitmp)

      end if

      if(savmem.eq.0.or.(abs(cphopt).ge.2.and.                          &
     &   abs(cphopt).lt.20.and.abs(cphopt).ne.11)) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nni,0.e0,nicep)

        if(abs(cphopt).eq.2) then

          call setcst4d(0,ni+1,0,nj+1,1,nk,1,1,0.e0,nitmp)

        else

          call setcst4d(0,ni+1,0,nj+1,1,nk,1,nni,0.e0,nitmp)

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

      if(savmem.eq.0.or.((nggopt.eq.1.or.                               &
     &   exbopt.ge.1.or.mod(lspopt,10).eq.1.or.vspopt.eq.1).and.        &
     &  (abs(cphopt).ge.2.and.abs(cphopt).lt.10))) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qigpv)
        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qitd)

      end if

      if(savmem.eq.0.or.                                                &
     &  (ngropt.ge.1.and.(abs(cphopt).ge.2.and.abs(cphopt).lt.10))) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qirdr)

      end if

      if(savmem.eq.0.or.                                                &
     &  (ngropt.eq.1.and.(abs(cphopt).ge.2.and.abs(cphopt).lt.10))) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qirtd)

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
        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqw,0.e0,qcwtmp)

      end if

      if(savmem.eq.0.or.cphopt.lt.0) then

        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qcicep)
        call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqi,0.e0,qcitmp)

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
        call s_setcst4d(0,ni+1,0,nj+1,1,nk,1,nqa(0),0.e0,qatmp)

      end if

      if(savmem.eq.0.or.(aslopt.ge.1.and.                               &
     &  (wbc.ge.4.or.ebc.ge.4.or.sbc.ge.4.or.nbc.ge.4))) then

        call s_setcst4d(1,nj,1,nk,1,2,1,nqa(0),0.e0,qacpx)
        call s_setcst4d(1,ni,1,nk,1,2,1,nqa(0),0.e0,qacpy)

      end if

      if(savmem.eq.0.or.((nggopt.eq.1.or.exbopt.ge.1.or.                &
     &   mod(lspopt,10).eq.1.or.vspopt.eq.1).and.aslopt.ge.1)) then

        call s_setcst4d(0,ni+1,0,nj+1,1,nk,1,nqa(0),0.e0,qagpv)
        call s_setcst4d(0,ni+1,0,nj+1,1,nk,1,nqa(0),0.e0,qatd)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.trkopt.ge.1)) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qt)

      end if

      if(savmem.eq.0.or.trkopt.ge.1) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qtp)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qttmp)

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
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tketmp)

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

      if(savmem.eq.0.or.(advopt.le.3.and.abs(cphopt).ge.1)) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qall)

      end if

      if(savmem.eq.0.or.abs(cphopt).ge.1) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qallp)

      end if

      if(savmem.eq.0.or.(advopt.le.3.and.sfcopt.ge.1)) then

        call setcst3d(0,ni+1,0,nj+1,1,nund,0.e0,tund)

      end if

      if(savmem.eq.0.or.sfcopt.ge.1) then

        call setcst3d(0,ni+1,0,nj+1,1,nund,0.e0,tundp)
        call setcst3d(0,ni+1,0,nj+1,1,nund,0.e0,tutmp)

        call setcst2d(0,ni+1,0,nj+1,0.e0,albe)

        call setcst2d(0,ni+1,0,nj+1,0.e0,beta)

        call setcst2d(0,ni+1,0,nj+1,0.e0,z0m)
        call setcst2d(0,ni+1,0,nj+1,0.e0,z0h)

        call setcst2d(0,ni+1,0,nj+1,0.e0,cap)
        call setcst2d(0,ni+1,0,nj+1,0.e0,nuu)

        call setcst2d(0,ni+1,0,nj+1,0.e0,kai)

      end if

      if(savmem.eq.0.or.(sfcopt.eq.3.or.sfcopt.eq.13)) then

        call setcst2d(0,ni+1,0,nj+1,0.e0,sst)
        call setcst2d(0,ni+1,0,nj+1,0.e0,sstd)

      end if

      if(savmem.eq.0.or.iniopt.eq.1) then

        call setcst1d(0,nlev,0.e0,z1d)
        call setcst1d(0,nlev,0.e0,u1d)
        call setcst1d(0,nlev,0.e0,v1d)
        call setcst1d(0,nlev,0.e0,p1d)
        call setcst1d(0,nlev,0.e0,pt1d)
        call setcst1d(0,nlev,0.e0,qv1d)

        call setcst1d(1,nlev,0.e0,ltmp1)
        call setcst1d(1,nlev,0.e0,ltmp2)
        call setcst1d(1,nlev,0.e0,ltmp3)

      end if

      if(savmem.eq.0.or.ngropt.eq.12) then

        call setcst2d(1,nid_rdr,1,njd_rdr,0.e0,lon_rdr)

        call setcst3d(1,nid_rdr,1,njd_rdr,1,nkd_rdr,0.e0,z_rdr)

        call setcst3d(1,nid_rdr,1,njd_rdr,1,km_rdr,0.e0,u_rdr)
        call setcst3d(1,nid_rdr,1,njd_rdr,1,km_rdr,0.e0,v_rdr)
        call setcst3d(1,nid_rdr,1,njd_rdr,1,km_rdr,0.e0,w_rdr)
        call setcst3d(1,nid_rdr,1,njd_rdr,1,km_rdr,0.e0,qp_rdr)

        call setcst3d(1,nid_rdr,1,njd_rdr,1,nk,0.e0,tmp1_rdr)

      end if

! -----

!! -----

      end subroutine s_allocslv

!-----7--------------------------------------------------------------7--

      end module m_allocslv
