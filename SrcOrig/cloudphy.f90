!***********************************************************************
      module m_cloudphy
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/01
!     Modification: 1999/11/19, 2000/01/17, 2000/03/08, 2000/04/18,
!                   2000/06/01, 2000/07/05, 2000/08/21, 2000/12/18,
!                   2001/01/15, 2001/06/29, 2001/10/18, 2002/01/07,
!                   2002/01/15, 2002/04/02, 2002/12/02, 2003/05/19,
!                   2003/12/12, 2004/03/22, 2004/04/01, 2004/06/10,
!                   2004/08/01, 2004/08/31, 2004/09/01, 2004/09/10,
!                   2004/09/25, 2004/10/12, 2004/12/17, 2005/04/04,
!                   2005/10/05, 2005/11/22, 2006/01/10, 2006/02/13,
!                   2006/04/03, 2006/05/12, 2006/07/21, 2006/09/30,
!                   2007/11/26, 2008/01/11, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2009/01/30, 2009/02/27, 2011/08/18,
!                   2011/09/22, 2012/06/19, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for the cloud micro physics.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_coldblk
      use m_comindx
      use m_diagnsg
      use m_diagnw
      use m_getexner
      use m_getiname
      use m_getvdens
      use m_warmbin
      use m_warmblk

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: cloudphy, s_cloudphy

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface cloudphy

        module procedure s_cloudphy

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
      subroutine s_cloudphy(fpadvopt,fpcphopt,fphaiopt,nclstp,          &
     &                      dtb,dtcl,ni,nj,nk,nqw,nnw,nqi,nni,km,       &
     &                      jcb,jcb8w,pbr,ptbr,rbr,rst,wf,ppp,          &
     &                      ptpp,qvp,qwtrp,qicep,qallp,ptpf,qvf,        &
     &                      qwtrf,nwtrf,nwtrp,qicef,nicef,nicep,        &
     &                      qcwtrf,qcicef,prwtr,price,rbv,pi,p,         &
     &                      tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,         &
     &                      tmp8,tmp9,tmp10,tmp11)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

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

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      integer, intent(in) :: nni
                       ! Number of categories of ice concentrations

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

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at past

      real, intent(in) :: qicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at past

      real, intent(in) :: qallp(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at past

! Input and output variables

      real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(inout) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(inout) :: qwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor at future

      real, intent(inout) :: nwtrf(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at future

      real, intent(inout) :: nwtrp(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations at past

      real, intent(inout) :: qicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor at future

      real, intent(inout) :: nicef(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at future

      real, intent(inout) :: nicep(0:ni+1,0:nj+1,1:nk,1:nni)
                       ! Ice concentrations at past

      real, intent(inout) :: qcwtrf(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at future

      real, intent(inout) :: qcicef(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at future

      real, intent(inout) :: prwtr(0:ni+1,0:nj+1,1:2,1:nqw)
                       ! Precipitation and accumulation for water

      real, intent(inout) :: price(0:ni+1,0:nj+1,1:2,1:nqi)
                       ! Precipitation and accumulation for ice

! Internal shared variables

      integer advopt   ! Option for advection scheme
      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes

      real dtb_sub     ! Substitute for dtb

      real, intent(inout) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(inout) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exnar function

      real, intent(inout) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

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

      real, intent(inout) :: tmp8(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

      real, intent(inout) :: tmp9(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

      real, intent(inout) :: tmp10(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

      real, intent(inout) :: tmp11(0:ni+1,0:nj+1,1:km)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpadvopt,advopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)

! -----

!! Perform the cloud physics.

      if(abs(cphopt).ge.1) then

! Reset the large time steps interval.

        if(advopt.le.3) then
          dtb_sub=2.e0*dtb
        else
          dtb_sub=dtb
        end if

! -----

! Get the inverse of the base state density.

        call getvdens(ni,nj,nk,rbr,rbv)

! -----

! Calculate the total pressure variable and Exner function.

        call getexner(ni,nj,nk,pbr,ppp,pi,p)

! -----

! Perform the bulk warm rain method.

        if(abs(cphopt).eq.1) then

          call s_warmblk(dtb_sub,ni,nj,nk,jcb,ptbr,rbr,rst,rbv,pi,p,    &
     &                   ptpp,qvp,qwtrp(0,0,1,1),qwtrp(0,0,1,2),        &
     &                   ptpf,qvf,qwtrf(0,0,1,1),qwtrf(0,0,1,2),        &
     &                   prwtr(0,0,1,2),tmp1,tmp2)

! -----

! Perform the bulk cold rain method.

        else if(abs(cphopt).eq.2) then

          call diagnw(ni,nj,nk,nqw,nnw,rbr,qwtrp,nwtrp)

          call diagnsg(idhaiopt,ni,nj,nk,nqi,nni,rbr,rbv,qicep,nicep)

          if(haiopt.eq.0) then

            call s_coldblk(idcphopt,idhaiopt,idqcgopt,idthresq,dtb_sub, &
     &                  ni,nj,nk,jcb,ptbr,rbr,rst,rbv,pi,p,ptpp,qvp,    &
     &                  qwtrp(0,0,1,1),qwtrp(0,0,1,2),qicep(0,0,1,1),   &
     &                  qicep(0,0,1,2),qicep(0,0,1,3),qicep(0,0,1,3),   &
     &                  nwtrp(0,0,1,1),nwtrp(0,0,1,2),nicep(0,0,1,1),   &
     &                  nicep(0,0,1,2),nicep(0,0,1,3),nicep(0,0,1,3),   &
     &                  qallp,ptpf,qvf,qwtrf(0,0,1,1),qwtrf(0,0,1,2),   &
     &                  qicef(0,0,1,1),qicef(0,0,1,2),qicef(0,0,1,3),   &
     &                  qicef(0,0,1,3),nwtrf(0,0,1,1),nwtrf(0,0,1,1),   &
     &                  nicef(0,0,1,1),nicef(0,0,1,2),nicef(0,0,1,3),   &
     &                  nicef(0,0,1,3),qcwtrf(0,0,1,1),qcwtrf(0,0,1,2), &
     &                  qcicef(0,0,1,1),qcicef(0,0,1,2),qcicef(0,0,1,3),&
     &                  qcicef(0,0,1,3),prwtr(0,0,1,1),prwtr(0,0,1,2),  &
     &                  price(0,0,1,1),price(0,0,1,2),price(0,0,1,3),   &
     &                  price(0,0,1,3),tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,   &
     &                  tmp7,tmp8,tmp9)

          else

            call s_coldblk(idcphopt,idhaiopt,idqcgopt,idthresq,dtb_sub, &
     &                  ni,nj,nk,jcb,ptbr,rbr,rst,rbv,pi,p,ptpp,qvp,    &
     &                  qwtrp(0,0,1,1),qwtrp(0,0,1,2),qicep(0,0,1,1),   &
     &                  qicep(0,0,1,2),qicep(0,0,1,3),qicep(0,0,1,4),   &
     &                  nwtrp(0,0,1,1),nwtrp(0,0,1,2),nicep(0,0,1,1),   &
     &                  nicep(0,0,1,2),nicep(0,0,1,3),nicep(0,0,1,4),   &
     &                  qallp,ptpf,qvf,qwtrf(0,0,1,1),qwtrf(0,0,1,2),   &
     &                  qicef(0,0,1,1),qicef(0,0,1,2),qicef(0,0,1,3),   &
     &                  qicef(0,0,1,4),nwtrf(0,0,1,1),nwtrf(0,0,1,1),   &
     &                  nicef(0,0,1,1),nicef(0,0,1,2),nicef(0,0,1,3),   &
     &                  nicef(0,0,1,4),qcwtrf(0,0,1,1),qcwtrf(0,0,1,2), &
     &                  qcicef(0,0,1,1),qcicef(0,0,1,2),qcicef(0,0,1,3),&
     &                  qcicef(0,0,1,4),prwtr(0,0,1,1),prwtr(0,0,1,2),  &
     &                  price(0,0,1,1),price(0,0,1,2),price(0,0,1,3),   &
     &                  price(0,0,1,4),tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,   &
     &                  tmp7,tmp8,tmp9)

          end if

! -----

! Perform the bulk cold rain method with ice concentrations.

        else if(abs(cphopt).eq.3) then

          call diagnw(ni,nj,nk,nqw,nnw,rbr,qwtrp,nwtrp)

          if(haiopt.eq.0) then

            call s_coldblk(idcphopt,idhaiopt,idqcgopt,idthresq,dtb_sub, &
     &                  ni,nj,nk,jcb,ptbr,rbr,rst,rbv,pi,p,ptpp,qvp,    &
     &                  qwtrp(0,0,1,1),qwtrp(0,0,1,2),qicep(0,0,1,1),   &
     &                  qicep(0,0,1,2),qicep(0,0,1,3),qicep(0,0,1,3),   &
     &                  nwtrp(0,0,1,1),nwtrp(0,0,1,2),nicep(0,0,1,1),   &
     &                  nicep(0,0,1,2),nicep(0,0,1,3),nicep(0,0,1,3),   &
     &                  qallp,ptpf,qvf,qwtrf(0,0,1,1),qwtrf(0,0,1,2),   &
     &                  qicef(0,0,1,1),qicef(0,0,1,2),qicef(0,0,1,3),   &
     &                  qicef(0,0,1,3),nwtrf(0,0,1,1),nwtrf(0,0,1,1),   &
     &                  nicef(0,0,1,1),nicef(0,0,1,2),nicef(0,0,1,3),   &
     &                  nicef(0,0,1,3),qcwtrf(0,0,1,1),qcwtrf(0,0,1,2), &
     &                  qcicef(0,0,1,1),qcicef(0,0,1,2),qcicef(0,0,1,3),&
     &                  qcicef(0,0,1,3),prwtr(0,0,1,1),prwtr(0,0,1,2),  &
     &                  price(0,0,1,1),price(0,0,1,2),price(0,0,1,3),   &
     &                  price(0,0,1,3),tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,   &
     &                  tmp7,tmp8,tmp9)

          else

            call s_coldblk(idcphopt,idhaiopt,idqcgopt,idthresq,dtb_sub, &
     &                  ni,nj,nk,jcb,ptbr,rbr,rst,rbv,pi,p,ptpp,qvp,    &
     &                  qwtrp(0,0,1,1),qwtrp(0,0,1,2),qicep(0,0,1,1),   &
     &                  qicep(0,0,1,2),qicep(0,0,1,3),qicep(0,0,1,4),   &
     &                  nwtrp(0,0,1,1),nwtrp(0,0,1,2),nicep(0,0,1,1),   &
     &                  nicep(0,0,1,2),nicep(0,0,1,3),nicep(0,0,1,4),   &
     &                  qallp,ptpf,qvf,qwtrf(0,0,1,1),qwtrf(0,0,1,2),   &
     &                  qicef(0,0,1,1),qicef(0,0,1,2),qicef(0,0,1,3),   &
     &                  qicef(0,0,1,4),nwtrf(0,0,1,1),nwtrf(0,0,1,1),   &
     &                  nicef(0,0,1,1),nicef(0,0,1,2),nicef(0,0,1,3),   &
     &                  nicef(0,0,1,4),qcwtrf(0,0,1,1),qcwtrf(0,0,1,2), &
     &                  qcicef(0,0,1,1),qcicef(0,0,1,2),qcicef(0,0,1,3),&
     &                  qcicef(0,0,1,4),prwtr(0,0,1,1),prwtr(0,0,1,2),  &
     &                  price(0,0,1,1),price(0,0,1,2),price(0,0,1,3),   &
     &                  price(0,0,1,4),tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,   &
     &                  tmp7,tmp8,tmp9)

          end if

! -----

! Perform the bulk cold rain method with water and ice concentrations.

        else if(abs(cphopt).eq.4) then

          if(haiopt.eq.0) then

            call s_coldblk(idcphopt,idhaiopt,idqcgopt,idthresq,dtb_sub, &
     &                  ni,nj,nk,jcb,ptbr,rbr,rst,rbv,pi,p,ptpp,qvp,    &
     &                  qwtrp(0,0,1,1),qwtrp(0,0,1,2),qicep(0,0,1,1),   &
     &                  qicep(0,0,1,2),qicep(0,0,1,3),qicep(0,0,1,3),   &
     &                  nwtrp(0,0,1,1),nwtrp(0,0,1,2),nicep(0,0,1,1),   &
     &                  nicep(0,0,1,2),nicep(0,0,1,3),nicep(0,0,1,3),   &
     &                  qallp,ptpf,qvf,qwtrf(0,0,1,1),qwtrf(0,0,1,2),   &
     &                  qicef(0,0,1,1),qicef(0,0,1,2),qicef(0,0,1,3),   &
     &                  qicef(0,0,1,3),nwtrf(0,0,1,1),nwtrf(0,0,1,2),   &
     &                  nicef(0,0,1,1),nicef(0,0,1,2),nicef(0,0,1,3),   &
     &                  nicef(0,0,1,3),qcwtrf(0,0,1,1),qcwtrf(0,0,1,2), &
     &                  qcicef(0,0,1,1),qcicef(0,0,1,2),qcicef(0,0,1,3),&
     &                  qcicef(0,0,1,3),prwtr(0,0,1,1),prwtr(0,0,1,2),  &
     &                  price(0,0,1,1),price(0,0,1,2),price(0,0,1,3),   &
     &                  price(0,0,1,3),tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,   &
     &                  tmp7,tmp8,tmp9)

          else

            call s_coldblk(idcphopt,idhaiopt,idqcgopt,idthresq,dtb_sub, &
     &                  ni,nj,nk,jcb,ptbr,rbr,rst,rbv,pi,p,ptpp,qvp,    &
     &                  qwtrp(0,0,1,1),qwtrp(0,0,1,2),qicep(0,0,1,1),   &
     &                  qicep(0,0,1,2),qicep(0,0,1,3),qicep(0,0,1,4),   &
     &                  nwtrp(0,0,1,1),nwtrp(0,0,1,2),nicep(0,0,1,1),   &
     &                  nicep(0,0,1,2),nicep(0,0,1,3),nicep(0,0,1,4),   &
     &                  qallp,ptpf,qvf,qwtrf(0,0,1,1),qwtrf(0,0,1,2),   &
     &                  qicef(0,0,1,1),qicef(0,0,1,2),qicef(0,0,1,3),   &
     &                  qicef(0,0,1,4),nwtrf(0,0,1,1),nwtrf(0,0,1,2),   &
     &                  nicef(0,0,1,1),nicef(0,0,1,2),nicef(0,0,1,3),   &
     &                  nicef(0,0,1,4),qcwtrf(0,0,1,1),qcwtrf(0,0,1,2), &
     &                  qcicef(0,0,1,1),qcicef(0,0,1,2),qcicef(0,0,1,3),&
     &                  qcicef(0,0,1,4),prwtr(0,0,1,1),prwtr(0,0,1,2),  &
     &                  price(0,0,1,1),price(0,0,1,2),price(0,0,1,3),   &
     &                  price(0,0,1,4),tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,   &
     &                  tmp7,tmp8,tmp9)

          end if

! -----

! Perform the warm bin method.

        else if(abs(cphopt).eq.11) then

          call s_warmbin(nclstp,dtb_sub,dtcl,ni,nj,nk,                  &
     &                   nqw,nnw,km,jcb,jcb8w,ptbr,rbr,rst,pi,wf,       &
     &                   rbv,p,ptpf,qvf,qwtrf,nwtrf,prwtr(0,0,1,1),     &
     &                   tmp7,tmp8,tmp9,tmp10,tmp11,tmp1,tmp2,tmp3,     &
     &                   tmp4,tmp5,tmp6)

        end if

! -----

      end if

!! -----

      end subroutine s_cloudphy

!-----7--------------------------------------------------------------7--

      end module m_cloudphy
