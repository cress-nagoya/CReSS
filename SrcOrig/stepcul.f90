!***********************************************************************
      module m_stepcul
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/04/03
!     Modification: 2006/05/12, 2006/06/21, 2006/07/21, 2006/09/21,
!                   2006/09/30, 2006/11/06, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/05/07, 2007/07/30, 2007/11/26,
!                   2008/05/02, 2008/06/09, 2008/07/01, 2008/08/25,
!                   2009/01/30, 2009/02/27, 2011/07/15, 2011/08/18,
!                   2011/09/22, 2011/11/10, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     solve the variables to the next time step by the Cubic Lagrange
!     scheme.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc4news
      use m_bcycle
      use m_combuf
      use m_comindx
      use m_copy3d
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getiname
      use m_hculs
      use m_hculuvw
      use m_lbcs
      use m_lbcu
      use m_lbcv
      use m_lbcw
      use m_phy2cnt
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy
      use m_var8w8s
      use m_vbcp
      use m_vbcqcg
      use m_vbcs
      use m_vbcu
      use m_vbcv
      use m_vbcw
      use m_vculs
      use m_vculs0
      use m_vculuvw

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: stepcul, s_stepcul

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface stepcul

        module procedure s_stepcul

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
      subroutine s_stepcul(fpcphopt,fphaiopt,fpqcgopt,fpaslopt,fptrkopt,&
     &                     fptubopt,fmois,nvstp,ivstp,dtb,dtsep,        &
     &                     ni,nj,nk,nqw,nnw,nqi,nni,nqa,j31,j32,jcb8w,  &
     &                     mf,mf8u,mf8v,uf,vf,wf,ppf,ptpf,qvf,          &
     &                     qwtrf,nwtrf,qicef,nicef,qcwtrf,qcicef,       &
     &                     qaslf,qtf,tkef,up,vp,wp,wc,wc8s,             &
     &                     tmp1,tmp2,tmp3,tmp4,tmp5)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fpqcgopt
                       ! Formal parameter of unique index of qcgopt

      integer, intent(in) :: fpaslopt
                       ! Formal parameter of unique index of aslopt

      integer, intent(in) :: fptrkopt
                       ! Formal parameter of unique index of trkopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: nvstp
                       ! Number of steps
                       ! of vertical Cubic Lagrange advection

      integer, intent(in) :: ivstp
                       ! Index of time steps
                       ! of vertical Cubic Lagrange advection

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

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: dtsep
                       ! Time steps interval
                       ! of vertical Cubic Lagrange advection

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(in) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

! Input and output variables

      real, intent(inout) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(inout) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(inout) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

      real, intent(inout) :: ppf(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at future

      real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(inout) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(inout) :: qwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at future

      real, intent(inout) :: nwtrf(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at future

      real, intent(inout) :: qicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at future

      real, intent(inout) :: nicef(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at future

      real, intent(inout) :: qcwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at future

      real, intent(inout) :: qcicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at future

      real, intent(inout) :: qaslf(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at future

      real, intent(inout) :: qtf(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at future

      real, intent(inout) :: tkef(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at future

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer qcgopt   ! Option for charging distribution
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing

      integer n        ! Array index in 4th direction

      integer ib       ! Exchanging variables number
      integer nb       ! Number of exchanging variables

      integer ni_sub   ! Substitute for ni
      integer nj_sub   ! Substitute for nj

      real, intent(inout) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(inout) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(inout) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(inout) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

      real, intent(inout) :: wc8s(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity
                       ! at scalar points

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpqcgopt,qcgopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)

! -----

! Set the substituted variables.

      ni_sub=ni
      nj_sub=nj

! -----

!!! Solve the variables to the next time step.

!! At the 1st step.

      if(ivstp.eq.1) then

! Calculate the velocity advection horizontally.

        call copy3d(0,ni+1,0,nj+1,1,nk,uf,up)
        call copy3d(0,ni+1,0,nj+1,1,nk,vf,vp)
        call copy3d(0,ni+1,0,nj+1,1,nk,wf,wp)

        call s_hculuvw(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,dtb,     &
     &                 ni,nj,nk,mf8u,mf8v,up,vp,wp,tmp1,tmp2,tmp3,wc,   &
     &                 tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4), &
     &                 tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

! -----

! Calculate the zeta components of contravariant velocity.

        call s_phy2cnt(idsthopt,idtrnopt,idmpopt,idmfcopt,idadvopt,     &
     &                 ni,nj,nk,j31,j32,jcb8w,mf,tmp1,tmp2,tmp3,wc,     &
     &                 wc8s,tmp4,tmp5)

! -----

! Set the bottom and the top boundary conditions for the z components of
! velocity.

        call s_vbcw(idbbc,idtbc,idmpopt,idmfcopt,ni,nj,nk,j31,j32,jcb8w,&
     &              mf,tmp1,tmp2,wc,tmp3,wc8s,tmp4,tmp5)

! -----

! zeta components of contravariant velocity be averaged to scalar
! points.

        call var8w8s(ni,nj,nk,wc,wc8s)

! -----

! Calculate the velocity advection vertically.

        call vculuvw(iddziv,ivstp,dtsep,                                &
     &               ni,nj,nk,wc,wc8s,tmp1,tmp2,tmp3,uf,vf,wf)

! -----

! Calculate the pressure perturbation advection.

        call copy3d(0,ni+1,0,nj+1,1,nk,ppf,tmp1)

        call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,           &
     &               dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,       &
     &               tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),   &
     &               tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

        call vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,ppf)

! -----

! Solve the potential temperature perturbation to the next time step.

        call copy3d(0,ni+1,0,nj+1,1,nk,ptpf,tmp1)

        call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,           &
     &               dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,       &
     &               tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),   &
     &               tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

        call vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,ptpf)

! -----

! Solve the hydrometeor to the next time step.

        if(fmois(1:5).eq.'moist') then

          call copy3d(0,ni+1,0,nj+1,1,nk,qvf,tmp1)

          call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,         &
     &                 dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,     &
     &                 tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4), &
     &                 tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

          call vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,qvf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              call s_copy3d(0,ni+1,0,nj+1,1,nk,qwtrf(0,0,1,1),tmp1)

              call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,     &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

              call s_vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,      &
     &                      qwtrf(0,0,1,1))

              call s_copy3d(0,ni+1,0,nj+1,1,nk,qwtrf(0,0,1,2),tmp1)

              call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,     &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

              call s_vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,      &
     &                      qwtrf(0,0,1,2))

            end if

            if(abs(cphopt).eq.4) then

              call s_copy3d(0,ni+1,0,nj+1,1,nk,nwtrf(0,0,1,1),tmp1)

              call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,     &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,       &
     &                     nwtrf(0,0,1,1))

              call s_copy3d(0,ni+1,0,nj+1,1,nk,nwtrf(0,0,1,2),tmp1)

              call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,     &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,       &
     &                     nwtrf(0,0,1,2))

            end if

            if(abs(cphopt).ge.2) then

              call s_copy3d(0,ni+1,0,nj+1,1,nk,qicef(0,0,1,1),tmp1)

              call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,     &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

              call s_vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,      &
     &                      qicef(0,0,1,1))

              call s_copy3d(0,ni+1,0,nj+1,1,nk,qicef(0,0,1,2),tmp1)

              call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,     &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

              call s_vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,      &
     &                      qicef(0,0,1,2))

              call s_copy3d(0,ni+1,0,nj+1,1,nk,qicef(0,0,1,3),tmp1)

              call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,     &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

              call s_vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,      &
     &                      qicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_copy3d(0,ni+1,0,nj+1,1,nk,qicef(0,0,1,4),tmp1)

                call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,   &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

                call s_vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,    &
     &                        qicef(0,0,1,4))

              end if

            end if

            if(abs(cphopt).eq.2) then

              call s_copy3d(0,ni+1,0,nj+1,1,nk,nicef(0,0,1,1),tmp1)

              call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,     &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,       &
     &                     nicef(0,0,1,1))

            else if(abs(cphopt).ge.3) then

              call s_copy3d(0,ni+1,0,nj+1,1,nk,nicef(0,0,1,1),tmp1)

              call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,     &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,       &
     &                     nicef(0,0,1,1))

              call s_copy3d(0,ni+1,0,nj+1,1,nk,nicef(0,0,1,2),tmp1)

              call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,     &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,       &
     &                     nicef(0,0,1,2))

              call s_copy3d(0,ni+1,0,nj+1,1,nk,nicef(0,0,1,3),tmp1)

              call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,     &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,       &
     &                     nicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_copy3d(0,ni+1,0,nj+1,1,nk,nicef(0,0,1,4),tmp1)

                call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,   &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,     &
     &                       nicef(0,0,1,4))

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                call s_copy3d(0,ni+1,0,nj+1,1,nk,qcwtrf(0,0,1,1),tmp1)

                call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,   &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,     &
     &                       qcwtrf(0,0,1,1))

                call s_copy3d(0,ni+1,0,nj+1,1,nk,qcwtrf(0,0,1,2),tmp1)

                call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,   &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,     &
     &                       qcwtrf(0,0,1,2))

              end if

              call s_copy3d(0,ni+1,0,nj+1,1,nk,qcicef(0,0,1,1),tmp1)

              call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,     &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,       &
     &                     qcicef(0,0,1,1))

              call s_copy3d(0,ni+1,0,nj+1,1,nk,qcicef(0,0,1,2),tmp1)

              call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,     &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,       &
     &                     qcicef(0,0,1,2))

              call s_copy3d(0,ni+1,0,nj+1,1,nk,qcicef(0,0,1,3),tmp1)

              call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,     &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,       &
     &                     qcicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_copy3d(0,ni+1,0,nj+1,1,nk,qcicef(0,0,1,4),tmp1)

                call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,   &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,     &
     &                       qcicef(0,0,1,4))

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                call s_copy3d(0,ni+1,0,nj+1,1,nk,qwtrf(0,0,1,n),tmp1)

                call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,   &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,     &
     &                       qwtrf(0,0,1,n))

              end do

              do n=1,nnw

                call s_copy3d(0,ni+1,0,nj+1,1,nk,nwtrf(0,0,1,n),tmp1)

                call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,   &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,     &
     &                       nwtrf(0,0,1,n))

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                call s_copy3d(0,ni+1,0,nj+1,1,nk,qicef(0,0,1,n),tmp1)

                call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,   &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,     &
     &                       qicef(0,0,1,n))

              end do

              do n=1,nni

                call s_copy3d(0,ni+1,0,nj+1,1,nk,nicef(0,0,1,n),tmp1)

                call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,   &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,     &
     &                       nicef(0,0,1,n))

              end do

            end if

          end if

        end if

! -----

! Solve the aerosol to the next time step.

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            call s_copy3d(0,ni+1,0,nj+1,1,nk,qaslf(0,0,1,n),tmp1)

            call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,       &
     &                  dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,    &
     &                  tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4),&
     &                  tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

            call s_vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,        &
     &                    qaslf(0,0,1,n))

          end do

        end if

! -----

! Solve the tracer to the next time step.

        if(trkopt.ge.1) then

          call copy3d(0,ni+1,0,nj+1,1,nk,qtf,tmp1)

          call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,         &
     &                 dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,     &
     &                 tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4), &
     &                 tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

          call vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,qtf)

        end if

! -----

! Solve the turbulent kinetic energy to the next time step.

        if(tubopt.ge.2) then

          call copy3d(0,ni+1,0,nj+1,1,nk,tkef,tmp1)

          call s_hculs(idmpopt,idmfcopt,idadvopt,iddxiv,iddyiv,         &
     &                 dtb,ni,nj,nk,mf8u,mf8v,up,vp,tmp1,tmp2,tmp3,     &
     &                 tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4), &
     &                 tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4))

          call vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp2,tkef)

        end if

! -----

!! -----

!! At the other steps.

      else

! zeta components of contravariant velocity be averaged to scalar
! points.

        call var8w8s(ni,nj,nk,wc,wc8s)

! -----

! Calculate the velocity advection vertically.

        call vculuvw(iddziv,ivstp,dtsep,                                &
     &               ni,nj,nk,wc,wc8s,tmp1,tmp2,tmp3,uf,vf,wf)

! -----

! Calculate the pressure perturbation advection vertically.

        call vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,ppf)

! -----

! Solve the potential temperature perturbation to the next time step.

        call vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,ptpf)

! -----

! Solve the hydrometeor to the next time step.

        if(fmois(1:5).eq.'moist') then

          call vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,qvf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              call s_vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,      &
     &                      qwtrf(0,0,1,1))

              call s_vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,      &
     &                      qwtrf(0,0,1,2))

            end if

            if(abs(cphopt).eq.4) then

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,       &
     &                     nwtrf(0,0,1,1))

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,       &
     &                     nwtrf(0,0,1,2))

            end if

            if(abs(cphopt).ge.2) then

              call s_vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,      &
     &                      qicef(0,0,1,1))

              call s_vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,      &
     &                      qicef(0,0,1,2))

              call s_vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,      &
     &                      qicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,    &
     &                        qicef(0,0,1,4))

              end if

            end if

            if(abs(cphopt).eq.2) then

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,       &
     &                     nicef(0,0,1,1))

            else if(abs(cphopt).ge.3) then

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,       &
     &                     nicef(0,0,1,1))

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,       &
     &                     nicef(0,0,1,2))

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,       &
     &                     nicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,     &
     &                       nicef(0,0,1,4))

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,     &
     &                       qcwtrf(0,0,1,1))

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,     &
     &                       qcwtrf(0,0,1,2))

              end if

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,       &
     &                     qcicef(0,0,1,1))

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,       &
     &                     qcicef(0,0,1,2))

              call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,       &
     &                     qcicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,     &
     &                       qcicef(0,0,1,4))

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,     &
     &                       qwtrf(0,0,1,n))

              end do

              do n=1,nnw

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,     &
     &                       nwtrf(0,0,1,n))

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,     &
     &                       qicef(0,0,1,n))

              end do

              do n=1,nni

                call s_vculs(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,     &
     &                       nicef(0,0,1,n))

              end do

            end if

          end if

        end if

! -----

! Solve the aerosol to the next time step.

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            call s_vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,        &
     &                    qaslf(0,0,1,n))

          end do

        end if

! -----

! Solve the tracer to the next time step.

        if(trkopt.ge.1) then

          call vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,qtf)

        end if

! -----

! Solve the turbulent kinetic energy to the next time step.

        if(tubopt.ge.2) then

          call vculs0(iddziv,ivstp,dtsep,ni,nj,nk,wc8s,tmp1,tkef)

        end if

! -----

      end if

!! -----

!!! -----

! Set the lateral boundary conditions.

      if(ivstp.eq.nvstp) then

        call lbcu(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,uf)
        call lbcv(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,vf)
        call lbcw(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,wf)

        call lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,ppf)

        call lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,ptpf)

        if(fmois(1:5).eq.'moist') then

          call lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,qvf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qwtrf(0,0,1,1))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qwtrf(0,0,1,2))

            end if

            if(abs(cphopt).eq.4) then

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    nwtrf(0,0,1,1))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    nwtrf(0,0,1,2))

            end if

            if(abs(cphopt).ge.2) then

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qicef(0,0,1,1))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qicef(0,0,1,2))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,  &
     &                      qicef(0,0,1,4))

              end if

            end if

            if(abs(cphopt).eq.2) then

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    nicef(0,0,1,1))

            else if(abs(cphopt).ge.3) then

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    nicef(0,0,1,1))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    nicef(0,0,1,2))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    nicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,  &
     &                      nicef(0,0,1,4))

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,  &
     &                      qcwtrf(0,0,1,1))

                call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,  &
     &                      qcwtrf(0,0,1,2))

              end if

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qcicef(0,0,1,1))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qcicef(0,0,1,2))

              call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,    &
     &                    qcicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,  &
     &                      qcicef(0,0,1,4))

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,  &
     &                      qwtrf(0,0,1,n))

              end do

              do n=1,nnw

                call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,  &
     &                      nwtrf(0,0,1,n))

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,  &
     &                      qicef(0,0,1,n))

              end do

              do n=1,nni

                call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,  &
     &                      nicef(0,0,1,n))

              end do

            end if

          end if

        end if

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            call s_lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,      &
     &                  qaslf(0,0,1,n))

          end do

        end if

        if(trkopt.ge.1) then

          call lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,qtf)

        end if

        if(tubopt.ge.2) then

          call lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,tkef)

        end if

      end if

! -----

!!!! Exchange the value horizontally.

      if(ivstp.eq.nvstp) then

! Count the number of variables.

        nb=5

        if(fmois(1:5).eq.'moist') then

          nb=nb+1

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

!!! Exchange the value horizontally between sub domain.

!! in x direction.

! Put the exchanging buffer.

        ib=0

        ib=ib+1

        call s_putbufsx(idwbc,idebc,'all',3,ni-2,ni,nj,nk,uf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,vf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,wf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,ppf,          &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,ptpf,         &
     &                  ib,nb,sbuf)

        if(fmois(1:5).eq.'moist') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,qvf,        &
     &                    ib,nb,sbuf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,sbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          nicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,1),ib,nb,sbuf)

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qcicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qwtrf(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          nwtrf(0,0,1,n),ib,nb,sbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qicef(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          nicef(0,0,1,n),ib,nb,sbuf)

              end do

            end if

          end if

        end if

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qaslf(0,0,1,n),ib,nb,sbuf)

          end do

        end if

        if(trkopt.ge.1) then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,qtf,        &
     &                    ib,nb,sbuf)

        end if

        if(tubopt.ge.2) then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,tkef,       &
     &                    ib,nb,sbuf)

        end if

! -----

! Call the exchanger.

        if(nb.ne.0) then

          call s_shiftsx(idwbc,idebc,'all',nj,nk,nb,sbuf,rbuf)

        end if

! -----

! Get the exchanging buffer.

        ib=0

        ib=ib+1

        call s_getbufsx(idwbc,idebc,'all',1,ni_sub,ni,nj,nk,uf,         &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,vf,           &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,wf,           &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,ppf,          &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,ptpf,         &
     &                  ib,nb,rbuf)

        if(fmois(1:5).eq.'moist') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,qvf,        &
     &                    ib,nb,rbuf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,rbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          nicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,1),ib,nb,rbuf)

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qcicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qwtrf(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          nwtrf(0,0,1,n),ib,nb,rbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qicef(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          nicef(0,0,1,n),ib,nb,rbuf)

              end do

            end if

          end if

        end if

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qaslf(0,0,1,n),ib,nb,rbuf)

          end do

        end if

        if(trkopt.ge.1) then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,qtf,        &
     &                    ib,nb,rbuf)

        end if

        if(tubopt.ge.2) then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,tkef,       &
     &                    ib,nb,rbuf)

        end if

! -----

!! -----

!! In y direction.

! Put the exchanging buffer.

        ib=0

        ib=ib+1

        call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,uf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufsy(idsbc,idnbc,'all',3,nj-2,ni,nj,nk,vf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,wf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,ppf,          &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,ptpf,         &
     &                  ib,nb,sbuf)

        if(fmois(1:5).eq.'moist') then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,qvf,        &
     &                    ib,nb,sbuf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qicef(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          qicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,sbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          nicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,1),ib,nb,sbuf)

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          qcicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          qwtrf(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          nwtrf(0,0,1,n),ib,nb,sbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          qicef(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          nicef(0,0,1,n),ib,nb,sbuf)

              end do

            end if

          end if

        end if

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qaslf(0,0,1,n),ib,nb,sbuf)

          end do

        end if

        if(trkopt.ge.1) then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,qtf,        &
     &                    ib,nb,sbuf)

        end if

        if(tubopt.ge.2) then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,tkef,       &
     &                    ib,nb,sbuf)

        end if

! -----

! Call the exchanger.

        if(nb.ne.0) then

          call s_shiftsy(idsbc,idnbc,'all',ni,nk,nb,sbuf,rbuf)

        end if

! -----

! Get the exchanging buffer.

        ib=0

        ib=ib+1

        call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,uf,           &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufsy(idsbc,idnbc,'all',1,nj_sub,ni,nj,nk,vf,         &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,wf,           &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,ppf,          &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,ptpf,         &
     &                  ib,nb,rbuf)

        if(fmois(1:5).eq.'moist') then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,qvf,        &
     &                    ib,nb,rbuf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          qicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,rbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          nicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,1),ib,nb,rbuf)

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          qcicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          qwtrf(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          nwtrf(0,0,1,n),ib,nb,rbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          qicef(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          nicef(0,0,1,n),ib,nb,rbuf)

              end do

            end if

          end if

        end if

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qaslf(0,0,1,n),ib,nb,rbuf)

          end do

        end if

        if(trkopt.ge.1) then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,qtf,        &
     &                    ib,nb,rbuf)

        end if

        if(tubopt.ge.2) then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,tkef,       &
     &                    ib,nb,rbuf)

        end if

! -----

!! -----

!!! -----

!!! Exchange the value horizontally between group domain.

!! in x direction.

! Put the exchanging buffer.

        ib=0

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',3,ni-2,ni,nj,nk,uf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,vf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,wf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,ppf,          &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,ptpf,         &
     &                  ib,nb,sbuf)

        if(fmois(1:5).eq.'moist') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,qvf,        &
     &                    ib,nb,sbuf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,sbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          nicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,1),ib,nb,sbuf)

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qcicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qwtrf(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          nwtrf(0,0,1,n),ib,nb,sbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qicef(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          nicef(0,0,1,n),ib,nb,sbuf)

              end do

            end if

          end if

        end if

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qaslf(0,0,1,n),ib,nb,sbuf)

          end do

        end if

        if(trkopt.ge.1) then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,qtf,        &
     &                    ib,nb,sbuf)

        end if

        if(tubopt.ge.2) then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,tkef,       &
     &                    ib,nb,sbuf)

        end if

! -----

! Call the exchanger.

        if(nb.ne.0) then

          call s_shiftgx(idwbc,idebc,'all',nj,nk,nb,sbuf,rbuf)

        end if

! -----

! Get the exchanging buffer.

        ib=0

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni_sub,ni,nj,nk,uf,         &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,vf,           &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,wf,           &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,ppf,          &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,ptpf,         &
     &                  ib,nb,rbuf)

        if(fmois(1:5).eq.'moist') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,qvf,        &
     &                    ib,nb,rbuf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,rbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          nicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,1),ib,nb,rbuf)

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qcicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qwtrf(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          nwtrf(0,0,1,n),ib,nb,rbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qicef(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          nicef(0,0,1,n),ib,nb,rbuf)

              end do

            end if

          end if

        end if

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qaslf(0,0,1,n),ib,nb,rbuf)

          end do

        end if

        if(trkopt.ge.1) then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,qtf,        &
     &                    ib,nb,rbuf)

        end if

        if(tubopt.ge.2) then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,tkef,       &
     &                    ib,nb,rbuf)

        end if

! -----

!! -----

!! In y direction.

! Put the exchanging buffer.

        ib=0

        ib=ib+1

        call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,uf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufgy(idsbc,idnbc,'all',3,nj-2,ni,nj,nk,vf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,wf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,ppf,          &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,ptpf,         &
     &                  ib,nb,sbuf)

        if(fmois(1:5).eq.'moist') then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,qvf,        &
     &                    ib,nb,sbuf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qicef(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          qicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,sbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        nicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          nicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,1),ib,nb,sbuf)

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          qcicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          qwtrf(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          nwtrf(0,0,1,n),ib,nb,sbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          qicef(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,      &
     &                          nicef(0,0,1,n),ib,nb,sbuf)

              end do

            end if

          end if

        end if

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,          &
     &                      qaslf(0,0,1,n),ib,nb,sbuf)

          end do

        end if

        if(trkopt.ge.1) then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,qtf,        &
     &                    ib,nb,sbuf)

        end if

        if(tubopt.ge.2) then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,tkef,       &
     &                    ib,nb,sbuf)

        end if

! -----

! Call the exchanger.

        if(nb.ne.0) then

          call s_shiftgy(idsbc,idnbc,'all',ni,nk,nb,sbuf,rbuf)

        end if

! -----

! Get the exchanging buffer.

        ib=0

        ib=ib+1

        call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,uf,           &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufgy(idsbc,idnbc,'all',1,nj_sub,ni,nj,nk,vf,         &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,wf,           &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,ppf,          &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,ptpf,         &
     &                  ib,nb,rbuf)

        if(fmois(1:5).eq.'moist') then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,qvf,        &
     &                    ib,nb,rbuf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          qicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,rbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        nicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          nicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,1),ib,nb,rbuf)

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          qcicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          qwtrf(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          nwtrf(0,0,1,n),ib,nb,rbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          qicef(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,      &
     &                          nicef(0,0,1,n),ib,nb,rbuf)

              end do

            end if

          end if

        end if

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,          &
     &                      qaslf(0,0,1,n),ib,nb,rbuf)

          end do

        end if

        if(trkopt.ge.1) then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,qtf,        &
     &                    ib,nb,rbuf)

        end if

        if(tubopt.ge.2) then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,tkef,       &
     &                    ib,nb,rbuf)

        end if

! -----

!! -----

!! in x direction again.

! Put the exchanging buffer.

        ib=0

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',3,ni-2,ni,nj,nk,uf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,vf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,wf,           &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,ppf,          &
     &                  ib,nb,sbuf)

        ib=ib+1

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,ptpf,         &
     &                  ib,nb,sbuf)

        if(fmois(1:5).eq.'moist') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,qvf,        &
     &                    ib,nb,sbuf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,1),ib,nb,sbuf)

             ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nwtrf(0,0,1,2),ib,nb,sbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,sbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,sbuf)

             ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        nicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          nicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,1),ib,nb,sbuf)

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,2),ib,nb,sbuf)

              end if

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,1),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,2),ib,nb,sbuf)

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,        &
     &                        qcicef(0,0,1,3),ib,nb,sbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qcicef(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qwtrf(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          nwtrf(0,0,1,n),ib,nb,sbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          qicef(0,0,1,n),ib,nb,sbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,      &
     &                          nicef(0,0,1,n),ib,nb,sbuf)

              end do

            end if

          end if

        end if

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,          &
     &                      qaslf(0,0,1,n),ib,nb,sbuf)

          end do

        end if

        if(trkopt.ge.1) then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,qtf,        &
     &                    ib,nb,sbuf)

        end if

        if(tubopt.ge.2) then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,tkef,       &
     &                    ib,nb,sbuf)

        end if

! -----

! Call the exchanger.

        if(nb.ne.0) then

          call s_shiftgx(idwbc,idebc,'all',nj,nk,nb,sbuf,rbuf)

        end if

! -----

! Get the exchanging buffer.

        ib=0

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni_sub,ni,nj,nk,uf,         &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,vf,           &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,wf,           &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,ppf,          &
     &                  ib,nb,rbuf)

        ib=ib+1

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,ptpf,         &
     &                  ib,nb,rbuf)

        if(fmois(1:5).eq.'moist') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,qvf,        &
     &                    ib,nb,rbuf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).eq.4) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nwtrf(0,0,1,2),ib,nb,rbuf)

            end if

            if(abs(cphopt).ge.2) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(abs(cphopt).eq.2) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,rbuf)

            else if(abs(cphopt).ge.3) then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        nicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          nicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,1),ib,nb,rbuf)

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qcwtrf(0,0,1,2),ib,nb,rbuf)

              end if

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,1),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,2),ib,nb,rbuf)

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,        &
     &                        qcicef(0,0,1,3),ib,nb,rbuf)

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qcicef(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qwtrf(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nnw

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          nwtrf(0,0,1,n),ib,nb,rbuf)

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          qicef(0,0,1,n),ib,nb,rbuf)

              end do

              do n=1,nni

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,      &
     &                          nicef(0,0,1,n),ib,nb,rbuf)

              end do

            end if

          end if

        end if

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            ib=ib+1

            call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,          &
     &                      qaslf(0,0,1,n),ib,nb,rbuf)

          end do

        end if

        if(trkopt.ge.1) then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,qtf,        &
     &                    ib,nb,rbuf)

        end if

        if(tubopt.ge.2) then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,tkef,       &
     &                    ib,nb,rbuf)

        end if

! -----

!! -----

!!! -----

      end if

!!!! -----

! Set the periodic boundary conditions.

      if(ivstp.eq.nvstp) then

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              3,1,ni-2,ni_sub,2,1,nj-2,nj-1,ni,nj,nk,uf)

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              2,1,ni-2,ni-1,3,1,nj-2,nj_sub,ni,nj,nk,vf)

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,wf)

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,ppf)

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,ptpf)

        if(fmois(1:5).eq.'moist') then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,qvf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qwtrf(0,0,1,1))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qwtrf(0,0,1,2))

            end if

            if(abs(cphopt).eq.4) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,nwtrf(0,0,1,1))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,nwtrf(0,0,1,2))

            end if

            if(abs(cphopt).ge.2) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qicef(0,0,1,1))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qicef(0,0,1,2))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        2,1,ni-2,ni-1,2,1,nj-2,nj-1,              &
     &                        ni,nj,nk,qicef(0,0,1,4))

              end if

            end if

            if(abs(cphopt).eq.2) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,nicef(0,0,1,1))

            else if(abs(cphopt).ge.3) then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,nicef(0,0,1,1))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,nicef(0,0,1,2))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,nicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        2,1,ni-2,ni-1,2,1,nj-2,nj-1,              &
     &                        ni,nj,nk,nicef(0,0,1,4))

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        2,1,ni-2,ni-1,2,1,nj-2,nj-1,              &
     &                        ni,nj,nk,qcwtrf(0,0,1,1))

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        2,1,ni-2,ni-1,2,1,nj-2,nj-1,              &
     &                        ni,nj,nk,qcwtrf(0,0,1,2))

              end if

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qcicef(0,0,1,1))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qcicef(0,0,1,2))

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qcicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        2,1,ni-2,ni-1,2,1,nj-2,nj-1,              &
     &                        ni,nj,nk,qcicef(0,0,1,4))

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        2,1,ni-2,ni-1,2,1,nj-2,nj-1,              &
     &                        ni,nj,nk,qwtrf(0,0,1,n))

              end do

              do n=1,nnw

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        2,1,ni-2,ni-1,2,1,nj-2,nj-1,              &
     &                        ni,nj,nk,nwtrf(0,0,1,n))

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        2,1,ni-2,ni-1,2,1,nj-2,nj-1,              &
     &                        ni,nj,nk,qicef(0,0,1,n))

              end do

              do n=1,nni

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        2,1,ni-2,ni-1,2,1,nj-2,nj-1,              &
     &                        ni,nj,nk,nicef(0,0,1,n))

              end do

            end if

          end if

        end if

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            call s_bcycle(idwbc,idebc,idsbc,idnbc,                      &
     &              2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,qaslf(0,0,1,n))

          end do

        end if

        if(trkopt.ge.1) then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,qtf)

        end if

        if(tubopt.ge.2) then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,tkef)

        end if

      end if

! -----

! Set the boundary conditions at the four corners.

      if(ivstp.eq.nvstp) then

        call bc4news(idwbc,idebc,idsbc,idnbc,1,ni_sub,1,nj-1,           &
     &               ni,nj,nk,uf)

        call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj_sub,           &
     &               ni,nj,nk,vf)

        call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,             &
     &               ni,nj,nk,wf)

        call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,             &
     &               ni,nj,nk,ppf)

        call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,             &
     &               ni,nj,nk,ptpf)

        if(fmois(1:5).eq.'moist') then

          call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,           &
     &                 ni,nj,nk,qvf)

          if(abs(cphopt).lt.10) then

            if(abs(cphopt).ge.1) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,qwtrf(0,0,1,1))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,qwtrf(0,0,1,2))

            end if

            if(abs(cphopt).eq.4) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,nwtrf(0,0,1,1))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,nwtrf(0,0,1,2))

            end if

            if(abs(cphopt).ge.2) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,qicef(0,0,1,1))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,qicef(0,0,1,2))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,qicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,   &
     &                         ni,nj,nk,qicef(0,0,1,4))

              end if

            end if

            if(abs(cphopt).eq.2) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,nicef(0,0,1,1))

            else if(abs(cphopt).ge.3) then

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,nicef(0,0,1,1))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,nicef(0,0,1,2))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,nicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,   &
     &                         ni,nj,nk,nicef(0,0,1,4))

              end if

            end if

            if(cphopt.lt.0) then

              if(qcgopt.eq.2) then

                call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,   &
     &                         ni,nj,nk,qcwtrf(0,0,1,1))

                call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,   &
     &                         ni,nj,nk,qcwtrf(0,0,1,2))

              end if

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,qcicef(0,0,1,1))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,qcicef(0,0,1,2))

              call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,     &
     &                       ni,nj,nk,qcicef(0,0,1,3))

              if(haiopt.eq.1) then

                call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,   &
     &                         ni,nj,nk,qcicef(0,0,1,4))

              end if

            end if

          else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

            if(abs(cphopt).ge.11) then

              do n=1,nqw

                call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,   &
     &                         ni,nj,nk,qwtrf(0,0,1,n))

              end do

              do n=1,nnw

                call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,   &
     &                         ni,nj,nk,nwtrf(0,0,1,n))

              end do

            end if

            if(abs(cphopt).eq.12) then

              do n=1,nqi

                call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,   &
     &                         ni,nj,nk,qicef(0,0,1,n))

              end do

              do n=1,nni

                call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,   &
     &                         ni,nj,nk,nicef(0,0,1,n))

              end do

            end if

          end if

        end if

        if(aslopt.ge.1) then

          do n=1,nqa(0)

            call s_bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,       &
     &                     ni,nj,nk,qaslf(0,0,1,n))

          end do

        end if

        if(trkopt.ge.1) then

          call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,           &
     &                 ni,nj,nk,qtf)

        end if

        if(tubopt.ge.2) then

          call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,           &
     &                 ni,nj,nk,tkef)

        end if

      end if

! -----

! Calculate the zeta components of contravariant velocity.

      call s_phy2cnt(idsthopt,idtrnopt,idmpopt,idmfcopt,idadvopt,       &
     &               ni,nj,nk,j31,j32,jcb8w,mf,uf,vf,wf,wc,             &
     &               tmp1,tmp2,tmp3)

! -----

! Set the bottom and the top boundary conditions.

      call vbcu(ni,nj,nk,uf)
      call vbcv(ni,nj,nk,vf)

      call s_vbcw(idbbc,idtbc,idmpopt,idmfcopt,ni,nj,nk,j31,j32,jcb8w,  &
     &            mf,uf,vf,wc,wf,tmp1,tmp2,tmp3)

      call vbcp(idbbc,ni,nj,nk,ppf)

      call vbcs(ni,nj,nk,ptpf)

      if(fmois(1:5).eq.'moist') then

        call vbcs(ni,nj,nk,qvf)

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            call s_vbcs(ni,nj,nk,qwtrf(0,0,1,1))
            call s_vbcs(ni,nj,nk,qwtrf(0,0,1,2))

          end if

          if(abs(cphopt).eq.4) then

            call s_vbcs(ni,nj,nk,nwtrf(0,0,1,1))
            call s_vbcs(ni,nj,nk,nwtrf(0,0,1,2))

          end if

          if(abs(cphopt).ge.2) then

            call s_vbcs(ni,nj,nk,qicef(0,0,1,1))
            call s_vbcs(ni,nj,nk,qicef(0,0,1,2))
            call s_vbcs(ni,nj,nk,qicef(0,0,1,3))

            if(haiopt.eq.1) then

              call s_vbcs(ni,nj,nk,qicef(0,0,1,4))

            end if

          end if

          if(abs(cphopt).eq.2) then

            call s_vbcs(ni,nj,nk,nicef(0,0,1,1))

          else if(abs(cphopt).ge.3) then

            call s_vbcs(ni,nj,nk,nicef(0,0,1,1))
            call s_vbcs(ni,nj,nk,nicef(0,0,1,2))
            call s_vbcs(ni,nj,nk,nicef(0,0,1,3))

            if(haiopt.eq.1) then

              call s_vbcs(ni,nj,nk,nicef(0,0,1,4))

            end if

          end if

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              call s_vbcqcg(ni,nj,nk,qcwtrf(0,0,1,1))
              call s_vbcqcg(ni,nj,nk,qcwtrf(0,0,1,2))

            end if

            call s_vbcqcg(ni,nj,nk,qcicef(0,0,1,1))
            call s_vbcqcg(ni,nj,nk,qcicef(0,0,1,2))
            call s_vbcqcg(ni,nj,nk,qcicef(0,0,1,3))

            if(haiopt.eq.1) then

              call s_vbcqcg(ni,nj,nk,qcicef(0,0,1,4))

            end if

          end if

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              call s_vbcs(ni,nj,nk,qwtrf(0,0,1,n))

            end do

            do n=1,nnw

              call s_vbcs(ni,nj,nk,nwtrf(0,0,1,n))

            end do

          end if

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              call s_vbcs(ni,nj,nk,qicef(0,0,1,n))

            end do

            do n=1,nni

              call s_vbcs(ni,nj,nk,nicef(0,0,1,n))

            end do

          end if

        end if

      end if

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          call s_vbcs(ni,nj,nk,qaslf(0,0,1,n))

        end do

      end if

      if(trkopt.ge.1) then

        call vbcs(ni,nj,nk,qtf)

      end if

      if(tubopt.ge.2) then

        call vbcs(ni,nj,nk,tkef)

      end if

! -----

      end subroutine s_stepcul

!-----7--------------------------------------------------------------7--

      end module m_stepcul
