!***********************************************************************
      module m_allocmph
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/08/08, 2003/11/05, 2003/12/12, 2004/09/25,
!                   2004/10/12, 2005/04/04, 2006/01/10, 2006/09/30,
!                   2007/01/20, 2007/01/31, 2007/11/26, 2008/01/11,
!                   2008/05/02, 2008/08/25, 2009/01/05, 2009/01/30,
!                   2009/02/27, 2011/01/14, 2011/08/18, 2011/09/22,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the array for cloud physics.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_combin
      use m_comblk
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getiname
      use m_setcst1d
      use m_setcst2d
      use m_setcst3d
      use m_temparam

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: allocmph, s_allocmph

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allocmph

        module procedure s_allocmph

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_allocmph(fpsavmem,fpcphopt,fphaiopt,fpaslopt,        &
     &                      ni,nj,nk,nqw)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsavmem
                       ! Formal parameter of unique index of savmem

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fpaslopt
                       ! Formal parameter of unique index of aslopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

! Internal shared variables

      integer savmem   ! Option for memory saving
      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer aslopt   ! Option for aerosol processes

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

      integer nqw_sub  ! Substitute for nqw

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsavmem,savmem)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpaslopt,aslopt)

! -----

! Set the substituted variable.

      nqw_sub=nqw

! ----

!!! Allocate the array for bulk cold rain method.

!! Allocate the array for terminal velocity.

! Perform allocate.

      stat=0

      if(savmem.eq.0.or.(abs(cphopt).ge.2.and.abs(cphopt).lt.10)) then

        allocate(ucq(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(urq(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(uiq(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(usq(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ugq(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(uhq(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ucn(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(urn(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(uin(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(usn(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ugn(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(uhn(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(ucq(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(urq(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(uiq(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(usq(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ugq(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(uhq(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ucn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(urn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(uin(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(usn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ugn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(uhn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('allocmph',8,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('allocmph',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

!! Allocate the array for temperature and suchlike.

! Perform allocate.

      stat=0

      if(aslopt.ge.1.or.elemnt_opt.eq.1.or.                             &
     &  (savmem.eq.0.and.(abs(cphopt).ge.2.and.abs(cphopt).lt.10))) then

        allocate(t(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(tcel(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvsst0(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvsw(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvsi(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(lv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ls(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(lf(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(kp(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mu(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(dv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mi(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqc(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqr(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqi(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqs(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqg(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vntr(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vnts(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vntg(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nuvi(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nuci(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clcr(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clcs(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clcg(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clri(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrs(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrg(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clir(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clis(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clig(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clsr(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clsg(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrsg(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrin(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrsn(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clsrn(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clsgn(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(agcn(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(agrn(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(agin(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(agsn(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vdvr(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vdvi(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vdvs(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vdvg(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cncr(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cnis(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cnsg(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cnsgn(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(spsi(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(spgi(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mlic(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mlsr(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mlgr(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(frrg(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(frrgn(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(shsr(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(shgr(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pgwet(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ecs(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(t(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(tcel(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvsst0(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvsw(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvsi(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(lv(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ls(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(lf(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(kp(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mu(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(dv(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mi(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqc(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqi(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqs(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vntr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vnts(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vntg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nuvi(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nuci(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clcr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clcs(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clcg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clri(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrs(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clir(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clis(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clig(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clsr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clsg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrsg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrin(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrsn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clsrn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clsgn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(agcn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(agrn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(agin(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(agsn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vdvr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vdvi(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vdvs(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vdvg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cncr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cnis(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cnsg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cnsgn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(spsi(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(spgi(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mlic(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mlsr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mlgr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(frrg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(frrgn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(shsr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(shgr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pgwet(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ecs(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

! -----

! Check errors in allocating.

      call chkerr(stat)

      if(aslopt.ge.1.or.elemnt_opt.eq.1) then

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('allocmph',8,'cont',5,'              ',14,101, &
     &                   stat)

          end if

          call cpondpe

          call destroy('allocmph',8,'stop',1001,'              ',14,101,&
     &                 stat)

        else

          write(fmem(1:6),'(a6)') 'enough'

        end if

      else if(savmem.eq.0.and.                                          &
     &       (abs(cphopt).ge.2.and.abs(cphopt).lt.10)) then

        if(stat.lt.0) then

          write(fmem(1:6),'(a6)') 'full  '

        else

          write(fmem(1:6),'(a6)') 'enough'

        end if

      else

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('allocmph',8,'cont',5,'              ',14,101, &
     &                   stat)

          end if

          call cpondpe

          call destroy('allocmph',8,'stop',1001,'              ',14,101,&
     &                 stat)

        else

          write(fmem(1:6),'(a6)') 'save  '

        end if

      end if

! -----

! Perform deallocate.

      if(fmem(1:4).eq.'full') then

        cstat=0

        if(cstat.eq.0) then

          deallocate(t,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(tcel,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(qvsst0,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(qvsw,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(qvsi,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(lv,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(ls,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(lf,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(kp,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(mu,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(dv,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(mi,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(diaqc,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(diaqr,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(diaqi,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(diaqs,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(diaqg,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(vntr,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(vnts,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(vntg,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(nuvi,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(nuci,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clcr,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clcs,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clcg,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clri,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clrs,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clrg,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clir,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clis,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clig,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clsr,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clsg,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clrsg,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clrin,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clrsn,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clsrn,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(clsgn,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(agcn,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(agrn,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(agin,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(agsn,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(vdvr,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(vdvi,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(vdvs,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(vdvg,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(cncr,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(cnis,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(cnsg,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(cnsgn,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(spsi,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(spgi,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(mlic,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(mlsr,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(mlgr,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(frrg,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(frrgn,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(shsr,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(shgr,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(pgwet,stat=cstat)

        end if

        if(cstat.eq.0) then

          deallocate(ecs,stat=cstat)

        end if

      end if

! -----

! Perform allocate again.

      if(fmem(1:4).eq.'full') then

        stat=0

        allocate(t(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(tcel(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvsst0(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvsw(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(qvsi(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(lv(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ls(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(lf(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(kp(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mu(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(dv(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mi(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqc(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqi(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqs(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(diaqg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vntr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vnts(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vntg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nuvi(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(nuci(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clcr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clcs(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clcg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clri(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrs(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clir(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clis(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clig(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clsr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clsg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrsg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrin(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clrsn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clsrn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(clsgn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(agcn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(agrn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(agin(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(agsn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vdvr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vdvi(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vdvs(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vdvg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cncr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cnis(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cnsg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(cnsgn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(spsi(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(spgi(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mlic(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mlsr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(mlgr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(frrg(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(frrgn(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(shsr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(shgr(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(pgwet(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ecs(0:0,0:0,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

! -----

! If error occured, call the procedure destroy.

      if(fmem(1:4).eq.'full') then

        call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('allocmph',8,'cont',5,'              ',14,101, &
     &                   stat)

          end if

          call cpondpe

          call destroy('allocmph',8,'stop',1001,'              ',14,101,&
     &                 stat)

        else

          write(fmem(1:6),'(a6)') 'save  '

        end if

      end if

! -----

!! -----

!!! -----

!! Allocate the array for warm bin method.

! Perform allocate.

      stat=0

      if(savmem.eq.0.or.(abs(cphopt).gt.10.and.abs(cphopt).lt.20)) then

        allocate(rbw(1:nqw),stat=cstat)

        stat=stat+abs(cstat)

        allocate(rrbw(1:nqw,1:5),stat=cstat)

        stat=stat+abs(cstat)

        allocate(brw(1:nqw+1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(rbrw(1:nqw+1,1:8),stat=cstat)

        stat=stat+abs(cstat)

        allocate(bmw(1:nqw+1,1:3),stat=cstat)

        stat=stat+abs(cstat)

        allocate(rbmw(1:nqw,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(dbmw(1:nqw),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ewbw(1:nqw,1:nqw),stat=cstat)

        stat=stat+abs(cstat)

      else

        allocate(rbw(1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(rrbw(1:1,1:5),stat=cstat)

        stat=stat+abs(cstat)

        allocate(brw(1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(rbrw(1:2,1:8),stat=cstat)

        stat=stat+abs(cstat)

        allocate(bmw(1:2,1:3),stat=cstat)

        stat=stat+abs(cstat)

        allocate(rbmw(1:1,1:2),stat=cstat)

        stat=stat+abs(cstat)

        allocate(dbmw(1:1),stat=cstat)

        stat=stat+abs(cstat)

        allocate(ewbw(1:1,1:1),stat=cstat)

        stat=stat+abs(cstat)

      end if

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('allocmph',8,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('allocmph',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

!! Fill in all array with 0.

! For bulk cold rain method.

      if(savmem.eq.0.or.(abs(cphopt).ge.2.and.abs(cphopt).lt.10)) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ucq)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,urq)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,uiq)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,usq)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ugq)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,uhq)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ucn)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,urn)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,uin)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,usn)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ugn)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,uhn)

      end if

      if(fmem(1:6).eq.'enough') then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,t)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tcel)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qvsst0)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qvsw)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qvsi)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,lv)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ls)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,lf)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,kp)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,mu)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,dv)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,mi)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,diaqc)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,diaqr)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,diaqi)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,diaqs)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,diaqg)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vntr)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vnts)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vntg)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,nuvi)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,nuci)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clcr)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clcs)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clcg)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clri)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clrs)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clrg)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clir)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clis)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clig)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clsr)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clsg)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clrsg)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clrin)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clrsn)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clsrn)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,clsgn)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,agcn)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,agrn)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,agin)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,agsn)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vdvr)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vdvi)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vdvs)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vdvg)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,cncr)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,cnis)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,cnsg)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,cnsgn)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,spsi)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,spgi)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,mlic)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,mlsr)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,mlgr)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,frrg)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,frrgn)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,shsr)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,shgr)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,pgwet)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ecs)

      end if

! -----

! For warm bin method.

      if(savmem.eq.0.or.(abs(cphopt).gt.10.and.abs(cphopt).lt.20)) then

        call setcst1d(1,nqw,0.e0,rbw)
        call setcst2d(1,nqw,1,5,0.e0,rrbw)

        call setcst1d(1,nqw+1,0.e0,brw)
        call setcst2d(1,nqw+1,1,8,0.e0,rbrw)

        call setcst2d(1,nqw+1,1,3,0.e0,bmw)
        call setcst2d(1,nqw,1,2,0.e0,rbmw)

        call setcst1d(1,nqw,0.e0,dbmw)

        call setcst2d(1,nqw,1,nqw_sub,0.e0,ewbw)

      end if

! -----

!! -----

      end subroutine s_allocmph

!-----7--------------------------------------------------------------7--

      end module m_allocmph
