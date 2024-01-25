!***********************************************************************
      module m_phasevbc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/07/13
!     Modification: 2001/08/07, 2001/09/13, 2002/04/02, 2002/06/06,
!                   2002/07/23, 2002/08/15, 2002/10/31, 2002/12/02,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2003/10/31,
!                   2003/11/05, 2003/11/28, 2003/12/12, 2004/04/15,
!                   2006/01/10, 2006/04/03, 2006/05/12, 2006/09/21,
!                   2006/11/06, 2007/01/20, 2007/05/07, 2007/05/21,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2008/10/10,
!                   2009/02/27, 2009/03/23, 2011/08/18, 2011/09/22,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the diffrential phase speed term between the external
!     boundary and model grid variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_getcname
      use m_getiname
      use m_inichar
      use m_phvbcuvw
      use m_phvbcs

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: phasevbc, s_phasevbc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface phasevbc

        module procedure s_phasevbc

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
      subroutine s_phasevbc(fpexbvar,fpgwmopt,fpcphopt,fphaiopt,        &
     &                      fpaslopt,fmois,dtb,dts,dtsep,gtinc,ni,nj,nk,&
     &                      nqw,nqi,nqa,rmf,rmf8u,rmf8v,u,up,uf,v,vp,vf,&
     &                      w,wp,wf,pp,ppp,ppf,ptp,ptpp,ptpf,qv,qvp,qvf,&
     &                      qwtr,qwtrp,qwtrf,qice,qicep,qicef,          &
     &                      qasl,qaslp,qaslf,ugpv,utd,vgpv,vtd,         &
     &                      wgpv,wtd,ppgpv,pptd,ptpgpv,ptptd,           &
     &                      qvgpv,qvtd,qwgpv,qwtd,qigpv,qitd,           &
     &                      qagpv,qatd,ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,   &
     &                      pcpx,pcpy,ptcpx,ptcpy,qvcpx,qvcpy,          &
     &                      qwcpx,qwcpy,qicpx,qicpy,qacpx,qacpy,        &
     &                      tmp1,tmp2,tmp3,tmp4)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpexbvar
                       ! Formal parameter of unique index of exbvar

      integer, intent(in) :: fpgwmopt
                       ! Formal parameter of unique index of gwmopt

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

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: dtsep
                       ! Time steps interval
                       ! of vertical Cubic Lagrange advection

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(in) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(in) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(in) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(in) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(in) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at present

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: ppf(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at future

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at present

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(in) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at present

      real, intent(in) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(in) :: qwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at future

      real, intent(in) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at present

      real, intent(in) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(in) :: qicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at future

      real, intent(in) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at present

      real, intent(in) :: qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at past

      real, intent(in) :: qaslf(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at future

      real, intent(in) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: utd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! x components of velocity of GPV data

      real, intent(in) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: vtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! y components of velocity of GPV data

      real, intent(in) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: wtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! z components of velocity of GPV data

      real, intent(in) :: ppgpv(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation of GPV data
                       ! at marked time

      real, intent(in) :: pptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! pressure perturbation of GPV data

      real, intent(in) :: ptpgpv(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation of GPV data
                       ! at marked time

      real, intent(in) :: ptptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! potential temperature perturbation of GPV data

      real, intent(in) :: qvgpv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, intent(in) :: qvtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! water vapor mixing ratio of GPV data

      real, intent(in) :: qwgpv(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of GPV data at marked time

      real, intent(in) :: qwtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of water hydrometeor of GPV data

      real, intent(in) :: qigpv(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of GPV data at marked time

      real, intent(in) :: qitd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of ice hydrometeor of GPV data

      real, intent(in) :: qagpv(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Mixing ratio of aerosol data at marked time

      real, intent(in) :: qatd(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Time tendency of mixing ratio of aerosol data

! Output variables

      real, intent(out) :: ucpx(1:nj,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on west and east boundary

      real, intent(out) :: ucpy(1:ni,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on south and north boundary

      real, intent(out) :: vcpx(1:nj,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on west and east boundary

      real, intent(out) :: vcpy(1:ni,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on south and north boundary

      real, intent(out) :: wcpx(1:nj,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on west and east boundary

      real, intent(out) :: wcpy(1:ni,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on south and north boundary

      real, intent(out) :: pcpx(1:nj,1:nk,1:2)
                       ! Phase speed of pressure
                       ! on west and east boundary

      real, intent(out) :: pcpy(1:ni,1:nk,1:2)
                       ! Phase speed of pressure
                       ! on south and north boundary

      real, intent(out) :: ptcpx(1:nj,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on west and east boundary

      real, intent(out) :: ptcpy(1:ni,1:nk,1:2)
                       ! Phase speed of potential temperature
                       ! on south and north boundary

      real, intent(out) :: qvcpx(1:nj,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on west and east boundary

      real, intent(out) :: qvcpy(1:ni,1:nk,1:2)
                       ! Phase speed of water vapor mixing ratio
                       ! on south and north boundary

      real, intent(out) :: qwcpx(1:nj,1:nk,1:2,1:nqw)
                       ! Phase speed of water hydrometeor
                       ! on west and east boundary

      real, intent(out) :: qwcpy(1:ni,1:nk,1:2,1:nqw)
                       ! Phase speed of water hydrometeor
                       ! on south and north boundary

      real, intent(out) :: qicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on west and east boundary

      real, intent(out) :: qicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on south and north boundary

      real, intent(out) :: qacpx(1:nj,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on west and east boundary

      real, intent(out) :: qacpy(1:ni,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on south and north boundary

! Internal shared variables

      character(len=108) exbvar
                       ! Control flag of
                       ! extrenal boundary forced variables

      integer gwmopt   ! Option for gravity wave mode integration
      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer aslopt   ! Option for aerosol processes

      integer n        ! Array index in 4th direction

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(exbvar)

! -----

! Get the required namelist variables.

      call getcname(fpexbvar,exbvar)
      call getiname(fpgwmopt,gwmopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpaslopt,aslopt)

! -----

! Calculate the diffrential phase speed term between the external
! boundary and model grid velocity variables.

      call s_phvbcuvw(idexbvar,idwbc,idebc,idsbc,idnbc,                 &
     &               idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,            &
     &               dtb,dts,gtinc,ni,nj,nk,rmf,rmf8u,rmf8v,            &
     &               u,up,uf,v,vp,vf,w,wp,wf,ugpv,utd,vgpv,vtd,wgpv,wtd,&
     &               ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the diffrential phase speed term between the external
! boundary and model grid pressure.

      if(exbvar(4:4).eq.'-') then

        call s_phvbcs(idwbc,idebc,idsbc,idnbc,idoneopt,                 &
     &                idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,           &
     &                'sml',dtb,dts,dtsep,gtinc,ni,nj,nk,rmf,u,v,       &
     &                pp,ppp,ppf,ppgpv,pptd,pcpx,pcpy,tmp1,tmp2)

      end if

! -----

! Calculate the diffrential phase speed term between the external
! boundary and model grid potential temperature.

      if(exbvar(5:5).eq.'-') then

        if(gwmopt.eq.0) then

          call s_phvbcs(idwbc,idebc,idsbc,idnbc,idoneopt,               &
     &                  idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,         &
     &                  'big',dtb,dts,dtsep,gtinc,ni,nj,nk,rmf,u,v,     &
     &                  ptp,ptpp,ptpf,ptpgpv,ptptd,ptcpx,ptcpy,         &
     &                  tmp1,tmp2)

        else

          call s_phvbcs(idwbc,idebc,idsbc,idnbc,idoneopt,               &
     &                  idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,         &
     &                  'sml',dtb,dts,dtsep,gtinc,ni,nj,nk,rmf,u,v,     &
     &                  ptp,ptpp,ptpf,ptpgpv,ptptd,ptcpx,ptcpy,         &
     &                  tmp1,tmp2)

        end if

      end if

! -----

!!! Calculate the diffrential phase speed term between the external
!!! boundary and model grid hydrometeor.

      if(fmois(1:5).eq.'moist') then

! Calculate the diffrential phase speed term between the external
! boundary and model grid water vapor mixing ratio.

        if(exbvar(6:6).eq.'-') then

          call s_phvbcs(idwbc,idebc,idsbc,idnbc,idoneopt,               &
     &                  idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,         &
     &                  'big',dtb,dts,dtsep,gtinc,ni,nj,nk,rmf,u,v,     &
     &                  qv,qvp,qvf,qvgpv,qvtd,qvcpx,qvcpy,tmp1,tmp2)

        end if

! -----

!! For the bulk categories.

        if(abs(cphopt).lt.10) then

! Calculate the diffrential phase speed term between the external
! boundary and model water hydrometeor.

          if(abs(cphopt).ge.1) then

            if(exbvar(7:7).eq.'-') then

              call s_phvbcs(idwbc,idebc,idsbc,idnbc,idoneopt,           &
     &                      idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,     &
     &                      'big',dtb,dts,dtsep,gtinc,ni,nj,nk,         &
     &                      rmf,u,v,qwtr(0,0,1,1),qwtrp(0,0,1,1),       &
     &                      qwtrf(0,0,1,1),qwgpv(0,0,1,1),qwtd(0,0,1,1),&
     &                      qwcpx(1,1,1,1),qwcpy(1,1,1,1),tmp1,tmp2)

              call s_phvbcs(idwbc,idebc,idsbc,idnbc,idoneopt,           &
     &                      idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,     &
     &                      'big',dtb,dts,dtsep,gtinc,ni,nj,nk,         &
     &                      rmf,u,v,qwtr(0,0,1,2),qwtrp(0,0,1,2),       &
     &                      qwtrf(0,0,1,2),qwgpv(0,0,1,2),qwtd(0,0,1,2),&
     &                      qwcpx(1,1,1,2),qwcpy(1,1,1,2),tmp1,tmp2)

            end if

          end if

! -----

! Calculate the diffrential phase speed term between the external
! boundary and model ice hydrometeor.

          if(abs(cphopt).ge.2) then

            if(exbvar(7:7).eq.'-') then

              call s_phvbcs(idwbc,idebc,idsbc,idnbc,idoneopt,           &
     &                      idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,     &
     &                      'big',dtb,dts,dtsep,gtinc,ni,nj,nk,         &
     &                      rmf,u,v,qice(0,0,1,1),qicep(0,0,1,1),       &
     &                      qicef(0,0,1,1),qigpv(0,0,1,1),qitd(0,0,1,1),&
     &                      qicpx(1,1,1,1),qicpy(1,1,1,1),tmp1,tmp2)

              call s_phvbcs(idwbc,idebc,idsbc,idnbc,idoneopt,           &
     &                      idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,     &
     &                      'big',dtb,dts,dtsep,gtinc,ni,nj,nk,         &
     &                      rmf,u,v,qice(0,0,1,2),qicep(0,0,1,2),       &
     &                      qicef(0,0,1,2),qigpv(0,0,1,2),qitd(0,0,1,2),&
     &                      qicpx(1,1,1,2),qicpy(1,1,1,2),tmp1,tmp2)

              call s_phvbcs(idwbc,idebc,idsbc,idnbc,idoneopt,           &
     &                      idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,     &
     &                      'big',dtb,dts,dtsep,gtinc,ni,nj,nk,         &
     &                      rmf,u,v,qice(0,0,1,3),qicep(0,0,1,3),       &
     &                      qicef(0,0,1,3),qigpv(0,0,1,3),qitd(0,0,1,3),&
     &                      qicpx(1,1,1,3),qicpy(1,1,1,3),tmp1,tmp2)

              if(haiopt.eq.1) then

                call s_phvbcs(idwbc,idebc,idsbc,idnbc,idoneopt,         &
     &                      idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,     &
     &                      'big',dtb,dts,dtsep,gtinc,ni,nj,nk,         &
     &                      rmf,u,v,qice(0,0,1,4),qicep(0,0,1,4),       &
     &                      qicef(0,0,1,4),qigpv(0,0,1,4),qitd(0,0,1,4),&
     &                      qicpx(1,1,1,4),qicpy(1,1,1,4),tmp1,tmp2)

              end if

            end if

          end if

! -----

        end if

!! -----

      end if

!!! -----

! Calculate the diffrential phase speed term between the external
! boundary and model aerosol mixing ratio.

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          call s_phvbcs(idwbc,idebc,idsbc,idnbc,idoneopt,               &
     &                  idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,         &
     &                  'big',dtb,dts,dtsep,gtinc,ni,nj,nk,             &
     &                  rmf,u,v,qasl(0,0,1,n),qaslp(0,0,1,n),           &
     &                  qaslf(0,0,1,n),qagpv(0,0,1,n),qatd(0,0,1,n),    &
     &                  qacpx(1,1,1,n),qacpy(1,1,1,n),tmp1,tmp2)

        end do

      end if

! -----

      end subroutine s_phasevbc

!-----7--------------------------------------------------------------7--

      end module m_phasevbc
