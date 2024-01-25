!***********************************************************************
      module m_phasev
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/09/06
!     Modification: 1999/09/14, 1999/10/12, 1999/10/27, 1999/11/01,
!                   1999/12/06, 2000/01/17, 2000/03/17, 2000/07/05,
!                   2001/01/15, 2001/02/13, 2001/04/15, 2001/05/29,
!                   2001/06/06, 2001/06/29, 2001/07/13, 2001/08/07,
!                   2001/09/13, 2002/04/02, 2002/06/06, 2002/08/15,
!                   2002/10/31, 2002/12/02, 2003/04/30, 2003/05/19,
!                   2003/09/01, 2003/10/31, 2003/11/05, 2003/11/28,
!                   2003/12/12, 2004/04/15, 2004/05/31, 2005/10/05,
!                   2006/01/10, 2006/02/13, 2006/04/03, 2006/05/12,
!                   2006/07/21, 2006/11/06, 2007/01/20, 2007/05/21,
!                   2007/11/26, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2008/10/10, 2009/01/30, 2009/02/27, 2011/08/18,
!                   2011/09/22, 2011/11/10, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the phase speed for the open boundary conditions.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_getiname
      use m_phvuvw
      use m_phvs

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: phasev, s_phasev

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface phasev

        module procedure s_phasev

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
      subroutine s_phasev(fpgwmopt,fpcphopt,fphaiopt,                   &
     &                    fpqcgopt,fpaslopt,fptrkopt,fptubopt,fmois,    &
     &                    dtb,dts,dtsep,ni,nj,nk,nqw,nnw,nqi,nni,nqa,   &
     &                    rmf,rmf8u,rmf8v,u,up,uf,v,vp,vf,w,wp,wf,      &
     &                    pp,ppp,ppf,ptp,ptpp,ptpf,qv,qvp,qvf,          &
     &                    qwtr,qwtrp,qwtrf,nwtr,nwtrp,nwtrf,            &
     &                    qice,qicep,qicef,nice,nicep,nicef,            &
     &                    qcwtr,qcwtrp,qcwtrf,qcice,qcicep,qcicef,      &
     &                    qasl,qaslp,qaslf,qt,qtp,qtf,tke,tkep,tkef,    &
     &                    ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,pcpx,pcpy,      &
     &                    ptcpx,ptcpy,qvcpx,qvcpy,qwcpx,qwcpy,          &
     &                    nwcpx,nwcpy,qicpx,qicpy,nicpx,nicpy,          &
     &                    qcwcpx,qcwcpy,qcicpx,qcicpy,qacpx,qacpy,      &
     &                    qtcpx,qtcpy,tkecpx,tkecpy,tmp1,tmp2,tmp3,tmp4)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpgwmopt
                       ! Formal parameter of unique index of gwmopt

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

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: dtsep
                       ! Time steps interval
                       ! of vertical Cubic Lagrange advection

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

      real, intent(in) :: nwtr(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at present

      real, intent(in) :: nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at past

      real, intent(in) :: nwtrf(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at future

      real, intent(in) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at present

      real, intent(in) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(in) :: qicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at future

      real, intent(in) :: nice(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at present

      real, intent(in) :: nicep(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at past

      real, intent(in) :: nicef(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at future

      real, intent(in) :: qcwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at present

      real, intent(in) :: qcwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at past

      real, intent(in) :: qcwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at future

      real, intent(in) :: qcice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at present

      real, intent(in) :: qcicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at past

      real, intent(in) :: qcicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at future

      real, intent(in) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at present

      real, intent(in) :: qaslp(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at past

      real, intent(in) :: qaslf(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio at future

      real, intent(in) :: qt(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at present

      real, intent(in) :: qtp(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at past

      real, intent(in) :: qtf(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio at future

      real, intent(in) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at present

      real, intent(in) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

      real, intent(in) :: tkef(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at future

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

      real, intent(out) :: nwcpx(1:nj,1:nk,1:2,1:nnw)
                       ! Phase speed of water concentrations
                       ! on west and east boundary

      real, intent(out) :: nwcpy(1:ni,1:nk,1:2,1:nnw)
                       ! Phase speed of water concentrations
                       ! on south and north boundary

      real, intent(out) :: qicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on west and east boundary

      real, intent(out) :: qicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of ice hydrometeor
                       ! on south and north boundary

      real, intent(out) :: nicpx(1:nj,1:nk,1:2,1:nni)
                       ! Phase speed of ice concentrations
                       ! on west and east boundary

      real, intent(out) :: nicpy(1:ni,1:nk,1:2,1:nni)
                       ! Phase speed of ice concentrations
                       ! on south and north boundary

      real, intent(out) :: qcwcpx(1:nj,1:nk,1:2,1:nqw)
                       ! Phase speed of charging distribution for water
                       ! on west and east boundary

      real, intent(out) :: qcwcpy(1:ni,1:nk,1:2,1:nqw)
                       ! Phase speed of charging distribution for water
                       ! on south and north boundary

      real, intent(out) :: qcicpx(1:nj,1:nk,1:2,1:nqi)
                       ! Phase speed of charging distribution for ice
                       ! on west and east boundary

      real, intent(out) :: qcicpy(1:ni,1:nk,1:2,1:nqi)
                       ! Phase speed of charging distribution for ice
                       ! on south and north boundary

      real, intent(out) :: qacpx(1:nj,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on west and east boundary

      real, intent(out) :: qacpy(1:ni,1:nk,1:2,1:nqa(0))
                       ! Phase speed of aerosol mixing ratio
                       ! on south and north boundary

      real, intent(out) :: qtcpx(1:nj,1:nk,1:2)
                       ! Phase speed of tracer mixing ratio
                       ! on west and east boundary

      real, intent(out) :: qtcpy(1:ni,1:nk,1:2)
                       ! Phase speed of tracer mixing ratio
                       ! on south and north boundary

      real, intent(out) :: tkecpx(1:nj,1:nk,1:2)
                       ! Phase speed of turbulent kinetic energy
                       ! on west and east boundary

      real, intent(out) :: tkecpy(1:ni,1:nk,1:2)
                       ! Phase speed of turbulent kinetic energy
                       ! on south and north boundary

! Internal shared variables

      integer gwmopt   ! Option for gravity wave mode integration
      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer qcgopt   ! Option for charging distribution
      integer aslopt   ! Option for aerosol processes
      integer trkopt   ! Option for mixing ratio tracking
      integer tubopt   ! Option for turbulent mixing

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

! Get the required namelist variables.

      call getiname(fpgwmopt,gwmopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpqcgopt,qcgopt)
      call getiname(fpaslopt,aslopt)
      call getiname(fptrkopt,trkopt)
      call getiname(fptubopt,tubopt)

! -----

! Calculate the phase speed of the velocity.

      call s_phvuvw(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,          &
     &              idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,dtb,dts,     &
     &              ni,nj,nk,rmf,rmf8u,rmf8v,u,up,uf,v,vp,vf,w,wp,wf,   &
     &              ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the phase speed of the pressure.

      call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,            &
     &            idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,      &
     &            'sml',4,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,pp,ppp,ppf,    &
     &            pcpx,pcpy,tmp1,tmp2)

! -----

! Calculate the phase speed of the potential temperature.

      if(gwmopt.eq.0) then

        call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,          &
     &              idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,    &
     &              'big',5,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,             &
     &              ptp,ptpp,ptpf,ptcpx,ptcpy,tmp1,tmp2)

      else

        call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,          &
     &              idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,    &
     &              'sml',5,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,             &
     &              ptp,ptpp,ptpf,ptcpx,ptcpy,tmp1,tmp2)

      end if

! -----

!!! Calculate the phase speed of the hydrometeor.

      if(fmois(1:5).eq.'moist') then

! Calculate the phase speed of the water vapor mixing ratio.

        call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,          &
     &              idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,    &
     &              'big',6,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,qv,qvp,qvf,  &
     &              qvcpx,qvcpy,tmp1,tmp2)

! -----

!! For the bulk categories.

        if(abs(cphopt).lt.10) then

! Calculate the phase speed of the water hydrometeor.

          if(abs(cphopt).ge.1) then

            call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,      &
     &                  idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,        &
     &                  idgwave,'big',7,dtb,dts,dtsep,ni,nj,nk,rmf,u,v, &
     &                  qwtr(0,0,1,1),qwtrp(0,0,1,1),qwtrf(0,0,1,1),    &
     &                  qwcpx(1,1,1,1),qwcpy(1,1,1,1),tmp1,tmp2)

            call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,      &
     &                  idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,        &
     &                  idgwave,'big',7,dtb,dts,dtsep,ni,nj,nk,rmf,u,v, &
     &                  qwtr(0,0,1,2),qwtrp(0,0,1,2),qwtrf(0,0,1,2),    &
     &                  qwcpx(1,1,1,2),qwcpy(1,1,1,2),tmp1,tmp2)

          end if

! -----

! Calculate the phase speed of the water concentrations.

          if(abs(cphopt).eq.4) then

            call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,      &
     &                  idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,        &
     &                  idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v, &
     &                  nwtr(0,0,1,1),nwtrp(0,0,1,1),nwtrf(0,0,1,1),    &
     &                  nwcpx(1,1,1,1),nwcpy(1,1,1,1),tmp1,tmp2)

            call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,      &
     &                  idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,        &
     &                  idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v, &
     &                  nwtr(0,0,1,2),nwtrp(0,0,1,2),nwtrf(0,0,1,2),    &
     &                  nwcpx(1,1,1,2),nwcpy(1,1,1,2),tmp1,tmp2)

          end if

! -----

! Calculate the phase speed of the ice hydrometeor.

          if(abs(cphopt).ge.2) then

            call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,      &
     &                  idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,        &
     &                  idgwave,'big',7,dtb,dts,dtsep,ni,nj,nk,rmf,u,v, &
     &                  qice(0,0,1,1),qicep(0,0,1,1),qicef(0,0,1,1),    &
     &                  qicpx(1,1,1,1),qicpy(1,1,1,1),tmp1,tmp2)

            call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,      &
     &                  idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,        &
     &                  idgwave,'big',7,dtb,dts,dtsep,ni,nj,nk,rmf,u,v, &
     &                  qice(0,0,1,2),qicep(0,0,1,2),qicef(0,0,1,2),    &
     &                  qicpx(1,1,1,2),qicpy(1,1,1,2),tmp1,tmp2)

            call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,      &
     &                  idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,        &
     &                  idgwave,'big',7,dtb,dts,dtsep,ni,nj,nk,rmf,u,v, &
     &                  qice(0,0,1,3),qicep(0,0,1,3),qicef(0,0,1,3),    &
     &                  qicpx(1,1,1,3),qicpy(1,1,1,3),tmp1,tmp2)

            if(haiopt.eq.1) then

              call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,    &
     &                   idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,       &
     &                   idgwave,'big',7,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,&
     &                   qice(0,0,1,4),qicep(0,0,1,4),qicef(0,0,1,4),   &
     &                   qicpx(1,1,1,4),qicpy(1,1,1,4),tmp1,tmp2)

            end if

          end if

! -----

! Calculate the phase speed of the ice concentrations.

          if(abs(cphopt).eq.2) then

            call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,      &
     &                  idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,        &
     &                  idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v, &
     &                  nice(0,0,1,1),nicep(0,0,1,1),nicef(0,0,1,1),    &
     &                  nicpx(1,1,1,1),nicpy(1,1,1,1),tmp1,tmp2)

          else if(abs(cphopt).ge.3) then

            call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,      &
     &                  idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,        &
     &                  idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v, &
     &                  nice(0,0,1,1),nicep(0,0,1,1),nicef(0,0,1,1),    &
     &                  nicpx(1,1,1,1),nicpy(1,1,1,1),tmp1,tmp2)

            call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,      &
     &                  idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,        &
     &                  idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v, &
     &                  nice(0,0,1,2),nicep(0,0,1,2),nicef(0,0,1,2),    &
     &                  nicpx(1,1,1,2),nicpy(1,1,1,2),tmp1,tmp2)

            call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,      &
     &                  idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,        &
     &                  idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v, &
     &                  nice(0,0,1,3),nicep(0,0,1,3),nicef(0,0,1,3),    &
     &                  nicpx(1,1,1,3),nicpy(1,1,1,3),tmp1,tmp2)

            if(haiopt.eq.1) then

              call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,    &
     &                   idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,       &
     &                   idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,&
     &                   nice(0,0,1,4),nicep(0,0,1,4),nicef(0,0,1,4),   &
     &                   nicpx(1,1,1,4),nicpy(1,1,1,4),tmp1,tmp2)

            end if

          end if

! -----

! Calculate the phase speed of the charging distributions.

          if(cphopt.lt.0) then

            if(qcgopt.eq.2) then

              call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,    &
     &                   idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,       &
     &                   idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,&
     &                   qcwtr(0,0,1,1),qcwtrp(0,0,1,1),qcwtrf(0,0,1,1),&
     &                   qcwcpx(1,1,1,1),qcwcpy(1,1,1,1),tmp1,tmp2)

              call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,    &
     &                   idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,       &
     &                   idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,&
     &                   qcwtr(0,0,1,2),qcwtrp(0,0,1,2),qcwtrf(0,0,1,2),&
     &                   qcwcpx(1,1,1,2),qcwcpy(1,1,1,2),tmp1,tmp2)

            end if

            call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,      &
     &                  idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,        &
     &                  idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v, &
     &                  qcice(0,0,1,1),qcicep(0,0,1,1),qcicef(0,0,1,1), &
     &                  qcicpx(1,1,1,1),qcicpy(1,1,1,1),tmp1,tmp2)

            call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,      &
     &                  idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,        &
     &                  idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v, &
     &                  qcice(0,0,1,2),qcicep(0,0,1,2),qcicef(0,0,1,2), &
     &                  qcicpx(1,1,1,2),qcicpy(1,1,1,2),tmp1,tmp2)

            call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,      &
     &                  idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,        &
     &                  idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v, &
     &                  qcice(0,0,1,3),qcicep(0,0,1,3),qcicef(0,0,1,3), &
     &                  qcicpx(1,1,1,3),qcicpy(1,1,1,3),tmp1,tmp2)

            if(haiopt.eq.1) then

              call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,    &
     &                   idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,       &
     &                   idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,&
     &                   qcice(0,0,1,4),qcicep(0,0,1,4),qcicef(0,0,1,4),&
     &                   qcicpx(1,1,1,4),qcicpy(1,1,1,4),tmp1,tmp2)

            end if

          end if

! -----

!! -----

!! For the bin categories.

        else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

! Calculate the phase speed of the water hydrometeor.

          if(abs(cphopt).ge.11) then

            do n=1,nqw

              call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,    &
     &                   idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,       &
     &                   idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,&
     &                   qwtr(0,0,1,n),qwtrp(0,0,1,n),qwtrf(0,0,1,n),   &
     &                   qwcpx(1,1,1,n),qwcpy(1,1,1,n),tmp1,tmp2)

            end do

            do n=1,nnw

              call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,    &
     &                   idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,       &
     &                   idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,&
     &                   nwtr(0,0,1,n),nwtrp(0,0,1,n),nwtrf(0,0,1,n),   &
     &                   nwcpx(1,1,1,n),nwcpy(1,1,1,n),tmp1,tmp2)

            end do

          end if

! -----

! Calculate the phase speed of the ice hydrometeor.

          if(abs(cphopt).eq.12) then

            do n=1,nqi

              call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,    &
     &                   idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,       &
     &                   idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,&
     &                   qice(0,0,1,n),qicep(0,0,1,n),qicef(0,0,1,n),   &
     &                   qicpx(1,1,1,n),qicpy(1,1,1,n),tmp1,tmp2)

            end do

            do n=1,nni

              call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,    &
     &                   idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,       &
     &                   idgwave,'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,&
     &                   nice(0,0,1,n),nicep(0,0,1,n),nicef(0,0,1,n),   &
     &                   nicpx(1,1,1,n),nicpy(1,1,1,n),tmp1,tmp2)

            end do

          end if

! -----

        end if

!! -----

      end if

!!! -----

! Calculate the phase speed of the aerosol.

      if(aslopt.ge.1) then

        do n=1,nqa(0)

          call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,        &
     &                idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,  &
     &                'big',8,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,           &
     &                qasl(0,0,1,n),qaslp(0,0,1,n),qaslf(0,0,1,n),      &
     &                qacpx(1,1,1,n),qacpy(1,1,1,n),tmp1,tmp2)

        end do

      end if

! -----

! Calculate the phase speed of the tracer.

      if(trkopt.ge.1) then

        call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,          &
     &              idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,    &
     &              'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,qt,qtp,qtf,  &
     &              qtcpx,qtcpy,tmp1,tmp2)

      end if

! -----

! Calculate the phase speed of the turbulent kinetic energy.

      if(tubopt.ge.2) then

        call s_phvs(idexbvar,idexbopt,idwbc,idebc,idsbc,idnbc,          &
     &              idoneopt,idmpopt,idmfcopt,iddxiv,iddyiv,idgwave,    &
     &              'big',9,dtb,dts,dtsep,ni,nj,nk,rmf,u,v,             &
     &              tke,tkep,tkef,tkecpx,tkecpy,tmp1,tmp2)

      end if

! -----

      end subroutine s_phasev

!-----7--------------------------------------------------------------7--

      end module m_phasev
