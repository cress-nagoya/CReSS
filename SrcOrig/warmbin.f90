!***********************************************************************
      module m_warmbin
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/09/30
!     Modification: 2007/01/20, 2007/07/30, 2008/01/11, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2009/01/30, 2009/02/27,
!                   2009/03/12, 2011/08/18

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the warm bin cloud physics.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_adjstbw
      use m_cg2mkbw
      use m_chkerr
      use m_coalbw
      use m_combin
      use m_comindx
      use m_commpi
      use m_copy3d
      use m_cpondpe
      use m_depsitbw
      use m_destroy
      use m_fallbw
      use m_getta2d
      use m_getta3d
      use m_inibinbw
      use m_mk2cgbw
      use m_termbw2d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: warmbin, s_warmbin

!-----7--------------------------------------------------------------7--

! Module variables

      real, allocatable :: rbr_sub(:,:,:)
                       ! Substitute for rbr

      real, allocatable :: rst_sub(:,:,:)
                       ! Substitute for rst

! Module procedure

      interface warmbin

        module procedure s_warmbin

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

!     none

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_warmbin(nclstp,dtb,dtcl,ni,nj,nk,nqw,nnw,km,         &
     &                     jcb,jcb8w,ptbr,rbr,rst,pi,wf,rbv,p,ptpf,qvf, &
     &                     qwbinf,nwbinf,prr,t,ubw,tmp1,tmp2,tmp3,      &
     &                     tmp4,tmp5,tmp6,tmp7,tmp8,tmp9)
!***********************************************************************

! Input variables

      integer, intent(in) :: nclstp(0:3)
                       ! Number of steps of cloud micro physics

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

      integer, intent(in) :: km
                       ! Dimension of max(nk, nqw, nqi)

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: dtcl(1:3)
                       ! Time interval of cloud micro physics

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exnar function

      real, intent(in) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

! Input and output variables

      real, intent(inout) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(inout) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(inout) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(inout) :: qwbinf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water mixing ratio at future

      real, intent(inout) :: nwbinf(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at future

      real, intent(inout) :: prr(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for rain

! Internal shared variables

      integer k        ! Array index in z direction

      integer iclstp   ! Index of steps of cloud micro physics

      integer iclact(1:2)
                       ! Activator of inner do loop

      integer iclint(1:2)
                       ! Interval of inner do loop

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

      real, intent(inout) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(inout) :: ubw(0:ni+1,0:nj+1,1:km)
                       ! Terminal velocity of water bin

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp6(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp7(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp8(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp9(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Allocate the substitute variables.

      stat=0

      allocate(rbr_sub(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(rst_sub(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('warmbin ',7,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('warmbin ',7,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      call copy3d(0,ni+1,0,nj+1,1,nk,rbr,rbr_sub)
      call copy3d(0,ni+1,0,nj+1,1,nk,rst,rst_sub)

! -----

! Change measurement from [m] [k] to [c] [g].

      call mk2cgbw(ni,nj,nk,nqw,nnw,rbr_sub,rst_sub,rbv,p,qwbinf,nwbinf,&
     &             prr)

! -----

! Adjust the mean water mixing ratio to be between their boundaries.

      call adjstbw(ni,nj,nk,nqw,nnw,bmw,qwbinf,nwbinf)

! -----

! Get the air temperature.

      call getta3d(ni,nj,nk,ptbr,pi,ptpf,t)

! -----

!!! Perform the deposition and coalescence.

! Initialize the activator of inner do loop.

      iclact(1)=1
      iclact(2)=1

! -----

! Set the interval of inner do loop.

      iclint(1)=nclstp(0)/nclstp(1)
      iclint(2)=nclstp(0)/nclstp(2)

! -----

!! Perform do loops in k and number of maximum steps among nclstp.

      do k=1,nk-1

        do iclstp=1,nclstp(0)

! Perform the depositon processes.

          if(iclact(1).eq.iclstp) then

            call s_depsitbw(k,dtcl(1),ni,nj,nk,nqw,nnw,                 &
     &                      rbv(0,0,k),pi(0,0,k),p(0,0,k),              &
     &                      t(0,0,k),brw,rbrw,bmw,rbmw,dbmw,            &
     &                      ptpf(0,0,k),qvf(0,0,k),qwbinf,nwbinf,       &
     &                      tmp1,tmp2,tmp3,tmp4(0,0,4),tmp5(0,0,4),     &
     &                      tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),        &
     &                      tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),        &
     &                      tmp6(0,0,1),tmp6(0,0,2),tmp6(0,0,3))

            call s_getta2d(ni,nj,ptbr(0,0,k),pi(0,0,k),ptpf(0,0,k),     &
     &                     t(0,0,k))

            iclact(1)=iclact(1)+iclint(1)

          end if

! -----

! Calculate the terminal velocity.

          call s_termbw2d(ni,nj,nqw,rbw,rrbw,                           &
     &                    rbr_sub(0,0,k),rbv(0,0,k),t(0,0,k),ubw,       &
     &                    tmp1(0,0,1),tmp1(0,0,2),tmp1(0,0,3),          &
     &                    tmp2(0,0,1),tmp2(0,0,2),tmp2(0,0,3))

! -----

! Perform the coalescence processes.

          if(iclact(2).eq.iclstp) then

            call s_coalbw(k,dtcl(2),ni,nj,nk,nqw,nnw,                   &
     &                    bmw,rbmw,dbmw,ewbw,ubw,qwbinf,nwbinf,         &
     &                    tmp5(0,0,1),tmp6(0,0,1),tmp6(0,0,3),          &
     &                    tmp5(0,0,2),tmp5(0,0,3),tmp1,tmp2,tmp3,tmp4)

            iclact(2)=iclact(2)+iclint(2)

          end if

! -----

        end do

      end do

!! -----

!!! -----

! Adjust the mean water mixing ratio to be between their boundaries.

      call adjstbw(ni,nj,nk,nqw,nnw,bmw,qwbinf,nwbinf)

! -----

! Get the air temperature.

      call getta3d(ni,nj,nk,ptbr,pi,ptpf,t)

! -----

! Perform the fall out.

      call fallbw(iddz,dtb,ni,nj,nk,nqw,nnw,rbw,rrbw,                   &
     &            jcb,rbr_sub,rst_sub,rbv,t,qwbinf,nwbinf,              &
     &            prr,ubw,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7)

! -----

! Adjust the mean water mixing ratio to be between their boundaries.

      call adjstbw(ni,nj,nk,nqw,nnw,bmw,qwbinf,nwbinf)

! -----

! Set the initial bin distribution of total water.

      call s_inibinbw(dtb,ni,nj,nk,nqw,nnw,                             &
     &                jcb8w,ptbr,rbv,pi,p,wf,t,brw,rbrw,                &
     &                ptpf,qvf,qwbinf,nwbinf,tmp3,tmp4,tmp5,tmp6,tmp7,  &
     &                tmp8(0,0,1),tmp8(0,0,2),tmp8(0,0,3),tmp8(0,0,4),  &
     &                tmp9(0,0,1),tmp9(0,0,2),tmp9(0,0,3),tmp1,tmp2)

! -----

! Adjust the mean water mixing ratio to be between their boundaries.

      call adjstbw(ni,nj,nk,nqw,nnw,bmw,qwbinf,nwbinf)

! -----

! Change measurement from [c] [g] to [m] [k].

      call cg2mkbw(ni,nj,nk,nqw,nnw,rbv,qwbinf,nwbinf,prr)

! -----

! Deallocate the substitute variables.

      stat=0

      deallocate(rbr_sub,stat=cstat)

      stat=stat+abs(cstat)

      deallocate(rst_sub,stat=cstat)

      stat=stat+abs(cstat)

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('warmbin ',7,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('warmbin ',7,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

      end subroutine s_warmbin

!-----7--------------------------------------------------------------7--

      end module m_warmbin
