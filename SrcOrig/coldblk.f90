!***********************************************************************
      module m_coldblk
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2000/12/18, 2001/04/15,
!                   2001/06/29, 2001/10/18, 2001/11/20, 2002/01/07,
!                   2002/01/15, 2002/04/02, 2002/06/06, 2002/09/09,
!                   2002/12/02, 2003/02/13, 2003/03/13, 2003/03/21,
!                   2003/03/28, 2003/04/30, 2003/05/19, 2003/11/05,
!                   2003/12/12, 2004/03/22, 2004/04/01, 2004/05/31,
!                   2004/06/10, 2004/08/01, 2004/09/01, 2004/09/25,
!                   2004/10/12, 2004/12/17, 2005/04/04, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/05/12, 2006/07/21,
!                   2006/09/30, 2007/01/20, 2007/05/07, 2007/11/26,
!                   2008/01/11, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2009/02/27, 2011/01/14, 2011/03/18, 2011/06/01,
!                   2011/09/22, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the source amounts of the bulk cold cloud micro physics.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_adjstnci
      use m_adjstni
      use m_adjstnp
      use m_adjstnsg
      use m_adjstnw
      use m_comblk
      use m_comindx
      use m_fallblk
      use m_getiname
      use m_getrname
      use m_srcblk
      use m_temparam
      use m_termblk

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: coldblk, s_coldblk

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface coldblk

        module procedure s_coldblk

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
      subroutine s_coldblk(fpcphopt,fphaiopt,fpqcgopt,fpthresq,         &
     &                     dtb,ni,nj,nk,jcb,ptbr,rbr,rst,rbv,           &
     &                     pi,p,ptpp,qvp,qcp,qrp,qip,qsp,qgp,qhp,       &
     &                     nccp,ncrp,ncip,ncsp,ncgp,nchp,qallp,ptpf,    &
     &                     qvf,qcf,qrf,qif,qsf,qgf,qhf,nccf,ncrf,ncif,  &
     &                     ncsf,ncgf,nchf,qccf,qrcf,qicf,qscf,qgcf,qhcf,&
     &                     prc,prr,pri,prs,prg,prh,tmp1,tmp2,tmp3,tmp4, &
     &                     tmp5,tmp6,tmp7,tmp8,tmp9)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fpqcgopt
                       ! Formal parameter of unique index of qcgopt

      integer, intent(in) :: fpthresq
                       ! Formal parameter of unique index of thresq

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(in) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exnar function

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qcp(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio at past

      real, intent(in) :: qrp(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio at past

      real, intent(in) :: qip(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio at past

      real, intent(in) :: qsp(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio at past

      real, intent(in) :: qgp(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio at past

      real, intent(in) :: qhp(0:ni+1,0:nj+1,1:nk)
                       ! Hail mixing ratio at past

      real, intent(in) :: nccp(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud water at past

      real, intent(in) :: ncrp(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of rain water at past

      real, intent(in) :: ncip(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud ice at past

      real, intent(in) :: ncsp(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow at past

      real, intent(in) :: ncgp(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of graupel at past

      real, intent(in) :: nchp(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of hail at past

      real, intent(in) :: qallp(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at past

! Input and output variables

      real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(inout) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(inout) :: qcf(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio at future

      real, intent(inout) :: qrf(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio at future

      real, intent(inout) :: qif(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio at future

      real, intent(inout) :: qsf(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio at future

      real, intent(inout) :: qgf(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio at future

      real, intent(inout) :: qhf(0:ni+1,0:nj+1,1:nk)
                       ! Hail mixing ratio at future

      real, intent(inout) :: nccf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud water at future

      real, intent(inout) :: ncrf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of rain water at future

      real, intent(inout) :: ncif(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud ice at future

      real, intent(inout) :: ncsf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow at future

      real, intent(inout) :: ncgf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of graupel at future

      real, intent(inout) :: nchf(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of hail at future

      real, intent(inout) :: qccf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of cloud water at future

      real, intent(inout) :: qrcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of rain water at future

      real, intent(inout) :: qicf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of cloud ice at future

      real, intent(inout) :: qscf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of snow at future

      real, intent(inout) :: qgcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of graupel at future

      real, intent(inout) :: qhcf(0:ni+1,0:nj+1,1:nk)
                       ! Charging distribution of hail at future

      real, intent(inout) :: prc(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for cloud water

      real, intent(inout) :: prr(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for rain

      real, intent(inout) :: pri(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for cloud ice

      real, intent(inout) :: prs(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for snow

      real, intent(inout) :: prg(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for graupel

      real, intent(inout) :: prh(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for hail

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer qcgopt   ! Option for charging distribution

      integer k        ! Array index in z direction

      real thresq      ! Minimum threshold value of mixing ratio

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

      real, intent(inout) :: tmp6(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp7(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp8(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp9(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpqcgopt,qcgopt)
      call getrname(fpthresq,thresq)

! -----

! Adjust the concentrations.

      if(abs(cphopt).eq.2) then

        call adjstnci(ni,nj,nk,qif,ncif)

      else if(abs(cphopt).eq.3) then

        call adjstni(idhaiopt,ni,nj,nk,                                 &
     &               rbr,rbv,qif,qsf,qgf,qhf,ncif,ncsf,ncgf,nchf)

      else if(abs(cphopt).eq.4) then

        call adjstnw(ni,nj,nk,rbr,rbv,qcf,qrf,nccf,ncrf)

        call adjstni(idhaiopt,ni,nj,nk,                                 &
     &               rbr,rbv,qif,qsf,qgf,qhf,ncif,ncsf,ncgf,nchf)

      end if

! -----

! Calculate the terminal velocity of the clod water, rain water, cloud
! ice, snow, graupel and hail.

      call termblk(idcphopt,idhaiopt,idthresq,ni,nj,nk,rbv,             &
     &            qcp,qrp,qip,qsp,qgp,qhp,nccp,ncrp,ncip,ncsp,ncgp,nchp,&
     &            ucq,urq,uiq,usq,ugq,uhq,ucn,urn,uin,usn,ugn,uhn)

! -----

!! Calculate the source amounts of the bulk cold cloud micro physics.

! If there is enough memory space, do loops perform 3 dimensionally.

      if(fmem(1:6).eq.'enough') then

        call srcblk(cphopt,haiopt,qcgopt,dtb,thresq,                    &
     &              ni,nj,nk,ptbr,rbr,rbv,pi,p,ptpp,qvp,qcp,qrp,        &
     &              qip,qsp,qgp,qhp,nccp,ncrp,ncip,ncsp,ncgp,nchp,      &
     &              qallp,urq,usq,ugq,uhq,urn,usn,ugn,uhn,ptpf,         &
     &              qvf,qcf,qrf,qif,qsf,qgf,qhf,nccf,ncrf,ncif,         &
     &              ncsf,ncgf,nchf,qccf,qrcf,qicf,qscf,qgcf,qhcf,       &
     &              t,tcel,qvsst0,qvsw,qvsi,lv,ls,lf,kp,mu,dv,mi,       &
     &              diaqc,diaqr,diaqi,diaqs,diaqg,vntr,vnts,vntg,       &
     &              nuvi,nuci,clcr,clcs,clcg,clri,clrs,clrg,            &
     &              clir,clis,clig,clsr,clsg,clrsg,clrin,               &
     &              clrsn,clsrn,clsgn,agcn,agrn,agin,agsn,              &
     &              vdvr,vdvi,vdvs,vdvg,cncr,cnis,cnsg,cnsgn,           &
     &              spsi,spgi,mlic,mlsr,mlgr,frrg,frrgn,                &
     &              shsr,shgr,pgwet,ecs)

! -----

! If there is not enough memory space, do loops perform 2 dimensionally.

      else if(fmem(1:4).eq.'save') then

        do k=1,nk-1

         call s_srcblk(cphopt,haiopt,qcgopt,dtb,thresq,ni,nj,1,         &
     &                 ptbr(0,0,k),rbr(0,0,k),rbv(0,0,k),pi(0,0,k),     &
     &                 p(0,0,k),ptpp(0,0,k),qvp(0,0,k),qcp(0,0,k),      &
     &                 qrp(0,0,k),qip(0,0,k),qsp(0,0,k),qgp(0,0,k),     &
     &                 qhp(0,0,k),nccp(0,0,k),ncrp(0,0,k),ncip(0,0,k),  &
     &                 ncsp(0,0,k),ncgp(0,0,k),nchp(0,0,k),qallp(0,0,k),&
     &                 urq(0,0,k),usq(0,0,k),ugq(0,0,k),uhq(0,0,k),     &
     &                 urn(0,0,k),usn(0,0,k),ugn(0,0,k),uhn(0,0,k),     &
     &                 ptpf(0,0,k),qvf(0,0,k),qcf(0,0,k),qrf(0,0,k),    &
     &                 qif(0,0,k),qsf(0,0,k),qgf(0,0,k),qhf(0,0,k),     &
     &                 nccf(0,0,k),ncrf(0,0,k),ncif(0,0,k),ncsf(0,0,k), &
     &                 ncgf(0,0,k),nchf(0,0,k),qccf(0,0,k),qrcf(0,0,k), &
     &                 qicf(0,0,k),qscf(0,0,k),qgcf(0,0,k),qhcf(0,0,k), &
     &                 tmp1(0,0,1),tmp1(0,0,2),tmp1(0,0,3),tmp1(0,0,4), &
     &                 tmp2(0,0,1),tmp2(0,0,2),tmp2(0,0,3),tmp2(0,0,4), &
     &                 tmp3(0,0,1),tmp3(0,0,2),tmp3(0,0,3),tmp3(0,0,4), &
     &                 tmp3(0,0,5),tmp3(0,0,6),tmp3(0,0,7),tmp3(0,0,8), &
     &                 tmp4(0,0,1),tmp4(0,0,2),tmp4(0,0,3),tmp4(0,0,4), &
     &                 tmp4(0,0,5),tmp4(0,0,6),tmp4(0,0,7),tmp4(0,0,8), &
     &                 tmp5(0,0,1),tmp5(0,0,2),tmp5(0,0,3),tmp5(0,0,4), &
     &                 tmp5(0,0,5),tmp5(0,0,6),tmp5(0,0,7),tmp5(0,0,8), &
     &                 tmp6(0,0,1),tmp6(0,0,2),tmp6(0,0,3),tmp6(0,0,4), &
     &                 tmp6(0,0,5),tmp6(0,0,6),tmp6(0,0,7),tmp6(0,0,8), &
     &                 tmp7(0,0,1),tmp7(0,0,2),tmp7(0,0,3),tmp7(0,0,4), &
     &                 tmp7(0,0,5),tmp7(0,0,6),tmp7(0,0,7),tmp7(0,0,8), &
     &                 tmp8(0,0,1),tmp8(0,0,2),tmp8(0,0,3),tmp8(0,0,4), &
     &                 tmp8(0,0,5),tmp8(0,0,6),tmp8(0,0,7),tmp8(0,0,8), &
     &                 tmp9(0,0,1),tmp9(0,0,2),tmp9(0,0,3),tmp9(0,0,4), &
     &                 tmp9(0,0,5))

        end do

      end if

! -----

!! -----

! Adjust the concentrations.

      if(abs(cphopt).eq.2) then

        call adjstnci(ni,nj,nk,qif,ncif)

      else if(abs(cphopt).eq.3) then

        call adjstni(idhaiopt,ni,nj,nk,                                 &
     &               rbr,rbv,qif,qsf,qgf,qhf,ncif,ncsf,ncgf,nchf)

      else if(abs(cphopt).eq.4) then

        call adjstnw(ni,nj,nk,rbr,rbv,qcf,qrf,nccf,ncrf)

        call adjstni(idhaiopt,ni,nj,nk,                                 &
     &               rbr,rbv,qif,qsf,qgf,qhf,ncif,ncsf,ncgf,nchf)

      end if

! -----

! Perform the fall out.

      call fallblk(idcphopt,idhaiopt,idqcgopt,iddz,dtb,ni,nj,nk,jcb,rbr,&
     &            rst,ucq,urq,uiq,usq,ugq,uhq,ucn,urn,uin,usn,ugn,uhn,  &
     &            qcf,qrf,qif,qsf,qgf,qhf,nccf,ncrf,ncif,ncsf,ncgf,nchf,&
     &            qccf,qrcf,qicf,qscf,qgcf,qhcf,prc,prr,pri,prs,prg,prh,&
     &            tmp1)

      if(flqcqi_opt.eq.0) then

        if(abs(cphopt).eq.3) then

          call adjstnsg(idhaiopt,ni,nj,nk,                              &
     &                  rbr,rbv,qsf,qgf,qhf,ncsf,ncgf,nchf)

        else if(abs(cphopt).eq.4) then

          call adjstnp(idhaiopt,ni,nj,nk,                               &
     &                 rbr,rbv,qrf,qsf,qgf,qhf,ncrf,ncsf,ncgf,nchf)

        end if

      else

        if(abs(cphopt).eq.2) then

          call adjstnci(ni,nj,nk,qif,ncif)

        else if(abs(cphopt).eq.3) then

          call adjstni(idhaiopt,ni,nj,nk,                               &
     &                 rbr,rbv,qif,qsf,qgf,qhf,ncif,ncsf,ncgf,nchf)

        else if(abs(cphopt).eq.4) then

          call adjstnw(ni,nj,nk,rbr,rbv,qcf,qrf,nccf,ncrf)

          call adjstni(idhaiopt,ni,nj,nk,                               &
     &                 rbr,rbv,qif,qsf,qgf,qhf,ncif,ncsf,ncgf,nchf)

        end if

      end if

! -----

      end subroutine s_coldblk

!-----7--------------------------------------------------------------7--

      end module m_coldblk
