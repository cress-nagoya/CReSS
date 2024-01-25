!***********************************************************************
      module m_forceqcg
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/02/13
!     Modification: 2006/04/03, 2006/05/12, 2006/06/21, 2006/07/21,
!                   2006/09/21, 2006/11/06, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/08/18, 2011/09/22, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the forcing terms in the charging distribution equation.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_advs
      use m_comindx
      use m_getcname
      use m_getiname
      use m_inichar
      use m_lsps0
      use m_smoo2s
      use m_smoo4s
      use m_turbflx
      use m_turbs
      use m_vsps0

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: forceqcg, s_forceqcg

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface forceqcg

        module procedure s_forceqcg

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_forceqcg(fplspvar,fpvspvar,fplspopt,fpvspopt,        &
     &                    fpsmtopt,fpcphopt,fphaiopt,fpqcgopt,fptubopt, &
     &                    ksp0,ni,nj,nk,nqw,nqi,j31,j32,jcb,jcb8u,jcb8v,&
     &                    mf,rmf,rmf8u,rmf8v,rbr,rst,rstxu,rstxv,rstxwc,&
     &                    qcwtr,qcwtrp,qcice,qcicep,rkh8u,rkh8v,rkv8w,  &
     &                    rbcxy,rbct,qcwfrc,qcifrc,h3,tmp1,tmp2,tmp3,   &
     &                    tmp4,tmp5)
!***********************************************************************

! Input variables

      integer, intent(in) :: fplspvar
                       ! Formal parameter of unique index of lspvar

      integer, intent(in) :: fpvspvar
                       ! Formal parameter of unique index of vspvar

      integer, intent(in) :: fplspopt
                       ! Formal parameter of unique index of lspopt

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fpqcgopt
                       ! Formal parameter of unique index of qcgopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: ksp0(1:2)
                       ! Index of lowest vertical sponge level

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

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: rstxu(0:ni+1,0:nj+1,1:nk)
                       ! u x base state density x Jacobian

      real, intent(in) :: rstxv(0:ni+1,0:nj+1,1:nk)
                       ! v x base state density x Jacobian

      real, intent(in) :: rstxwc(0:ni+1,0:nj+1,1:nk)
                       ! wc x base state density x Jacobian

      real, intent(in) :: qcwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at present

      real, intent(in) :: qcwtrp(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Charging distribution for water at past

      real, intent(in) :: qcice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at present

      real, intent(in) :: qcicep(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Charging distribution for ice at past

      real, intent(in) :: rkh8u(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x rbr x horizontal eddy diffusivity / jcb
                       ! at u points

      real, intent(in) :: rkh8v(0:ni+1,0:nj+1,1:nk)
                       ! 2.0 x rbr x horizontal eddy diffusivity / jcb
                       ! at v points

      real, intent(in) :: rkv8w(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy diffusivity / jcb
                       ! at w points

      real, intent(in) :: rbcxy(1:ni,1:nj)
                       ! Relaxed lateral sponge damping coefficients

      real, intent(in) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

! Output variables

      real, intent(out) :: qcwfrc(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Forcing terms
                       ! in charging distribution for water equations

      real, intent(out) :: qcifrc(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Forcing terms
                       ! in charging distribution for ice equations

! Internal shared variables

      character(len=108) lspvar
                       ! Control flag of
                       ! lateral sponge damped variables

      character(len=108) vspvar
                       ! Control flag of
                       ! vertical sponge damped variables

      integer lspopt   ! Option for lateral sponge damping
      integer vspopt   ! Option for vertical sponge damping

      integer smtopt   ! Option for numerical smoothing
      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer qcgopt   ! Option for charging distribution
      integer tubopt   ! Option for turbulent mixing

      real, intent(inout) :: h3(0:ni+1,0:nj+1,1:nk)
                       ! z components of turbulent fluxes

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

! Initialize the character variables.

      call inichar(lspvar)
      call inichar(vspvar)

! -----

! Get the required namelist variables.

      call getcname(fplspvar,lspvar)
      call getcname(fpvspvar,vspvar)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getiname(fpsmtopt,smtopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpqcgopt,qcgopt)
      call getiname(fptubopt,tubopt)

! -----

!!!! Calculate the forcing terms in the charging distribution equation.

      if(cphopt.lt.0) then

!!! For the water hydrometeor.

        if(qcgopt.eq.2) then

!! For the cloud water.

! Calculate the advection.

          call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,       &
     &               iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,  &
     &               qcwtr(0,0,1,1),qcwfrc(0,0,1,1),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

          if(mod(smtopt,10).eq.1) then

            call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,                   &
     &                    rbr,qcwtrp(0,0,1,1),qcwfrc(0,0,1,1),tmp1)

          end if

! -----

! Calculate the 4th order smoothing.

          if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

            call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,   &
     &                   idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qcwtrp(0,0,1,1),&
     &                   qcwfrc(0,0,1,1),tmp1,tmp2,tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the turbulent mixing.

          if(tubopt.ge.1) then

            call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,      &
     &                     ni,nj,nk,j31,j32,jcb,qcwtrp(0,0,1,1),        &
     &                     qcwfrc(0,0,1,1),rkh8u,rkh8v,rkv8w,           &
     &                     tmp1,tmp2,h3,tmp3,tmp4,tmp5)

            call s_turbs(idtrnopt,idmpopt,idmfcopt,                     &
     &                   iddxiv,iddyiv,iddziv,ni,nj,nk,                 &
     &                   j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,        &
     &                   tmp1,tmp2,h3,qcwfrc(0,0,1,1),tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the lateral sponge damping.

          if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

            call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,ni,nj,nk,  &
     &                   rst,qcwtrp(0,0,1,1),rbcxy,qcwfrc(0,0,1,1),tmp1)

          end if

! -----

! Calculate the vertical sponge damping.

          if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

            call s_vsps0(ksp0,ni,nj,nk,rst,qcwtrp(0,0,1,1),rbct,        &
     &                   qcwfrc(0,0,1,1))

          end if

! -----

!! -----

!! For the rain water.

! Calculate the advection.

          call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,       &
     &               iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,  &
     &               qcwtr(0,0,1,2),qcwfrc(0,0,1,2),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

          if(mod(smtopt,10).eq.1) then

            call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,                   &
     &                    rbr,qcwtrp(0,0,1,2),qcwfrc(0,0,1,2),tmp1)

          end if

! -----

! Calculate the 4th order smoothing.

          if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

            call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,   &
     &                   idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qcwtrp(0,0,1,2),&
     &                   qcwfrc(0,0,1,2),tmp1,tmp2,tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the turbulent mixing.

          if(tubopt.ge.1) then

            call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,      &
     &                     ni,nj,nk,j31,j32,jcb,qcwtrp(0,0,1,2),        &
     &                     qcwfrc(0,0,1,2),rkh8u,rkh8v,rkv8w,           &
     &                     tmp1,tmp2,h3,tmp3,tmp4,tmp5)

            call s_turbs(idtrnopt,idmpopt,idmfcopt,                     &
     &                   iddxiv,iddyiv,iddziv,ni,nj,nk,                 &
     &                   j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,        &
     &                   tmp1,tmp2,h3,qcwfrc(0,0,1,2),tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the lateral sponge damping.

          if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

            call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,ni,nj,nk,  &
     &                   rst,qcwtrp(0,0,1,2),rbcxy,qcwfrc(0,0,1,2),tmp1)

          end if

! -----

! Calculate the vertical sponge damping.

          if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

            call s_vsps0(ksp0,ni,nj,nk,rst,qcwtrp(0,0,1,2),rbct,        &
     &                   qcwfrc(0,0,1,2))

          end if

! -----

!! -----

        end if

!!! -----

!!! For the ice hydrometeor.

!! For the cloud ice.

! Calculate the advection.

        call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,         &
     &              iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,   &
     &              qcice(0,0,1,1),qcifrc(0,0,1,1),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

        if(mod(smtopt,10).eq.1) then

          call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qcicep(0,0,1,1), &
     &                  qcifrc(0,0,1,1),tmp1)

        end if

! -----

! Calculate the 4th order smoothing.

        if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

          call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,     &
     &                  idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qcicep(0,0,1,1), &
     &                  qcifrc(0,0,1,1),tmp1,tmp2,tmp3,tmp4,tmp5)

        end if

! -----

! Calculate the turbulent mixing.

        if(tubopt.ge.1) then

          call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,        &
     &                   ni,nj,nk,j31,j32,jcb,qcicep(0,0,1,1),          &
     &                   qcifrc(0,0,1,1),rkh8u,rkh8v,rkv8w,             &
     &                   tmp1,tmp2,h3,tmp3,tmp4,tmp5)

          call s_turbs(idtrnopt,idmpopt,idmfcopt,                       &
     &                 iddxiv,iddyiv,iddziv,ni,nj,nk,                   &
     &                 j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,          &
     &                 tmp1,tmp2,h3,qcifrc(0,0,1,1),tmp3,tmp4,tmp5)

        end if

! -----

! Calculate the lateral sponge damping.

        if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

          call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,ni,nj,nk,    &
     &                 rst,qcicep(0,0,1,1),rbcxy,qcifrc(0,0,1,1),tmp1)

        end if

! -----

! Calculate the vertical sponge damping.

        if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

          call s_vsps0(ksp0,ni,nj,nk,rst,qcicep(0,0,1,1),rbct,          &
     &                 qcifrc(0,0,1,1))

        end if

! -----

!! -----

!! For the snow.

! Calculate the advection.

        call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,         &
     &              iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,   &
     &              qcice(0,0,1,2),qcifrc(0,0,1,2),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

        if(mod(smtopt,10).eq.1) then

          call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qcicep(0,0,1,2), &
     &                  qcifrc(0,0,1,2),tmp1)

        end if

! -----

! Calculate the 4th order smoothing.

        if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

          call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,     &
     &                  idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qcicep(0,0,1,2), &
     &                  qcifrc(0,0,1,2),tmp1,tmp2,tmp3,tmp4,tmp5)

        end if

! -----

! Calculate the turbulent mixing.

        if(tubopt.ge.1) then

          call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,        &
     &                   ni,nj,nk,j31,j32,jcb,qcicep(0,0,1,2),          &
     &                   qcifrc(0,0,1,2),rkh8u,rkh8v,rkv8w,             &
     &                   tmp1,tmp2,h3,tmp3,tmp4,tmp5)

          call s_turbs(idtrnopt,idmpopt,idmfcopt,                       &
     &                 iddxiv,iddyiv,iddziv,ni,nj,nk,                   &
     &                 j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,          &
     &                 tmp1,tmp2,h3,qcifrc(0,0,1,2),tmp3,tmp4,tmp5)

        end if

! -----

! Calculate the lateral sponge damping.

        if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

          call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,ni,nj,nk,    &
     &                 rst,qcicep(0,0,1,2),rbcxy,qcifrc(0,0,1,2),tmp1)

        end if

! -----

! Calculate the vertical sponge damping.

        if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

          call s_vsps0(ksp0,ni,nj,nk,rst,qcicep(0,0,1,2),rbct,          &
     &                 qcifrc(0,0,1,2))

        end if

! -----

!! -----

!! For the graupel.

! Calculate the advection.

        call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,         &
     &              iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,   &
     &              qcice(0,0,1,3),qcifrc(0,0,1,3),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

        if(mod(smtopt,10).eq.1) then

          call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qcicep(0,0,1,3), &
     &                  qcifrc(0,0,1,3),tmp1)

        end if

! -----

! Calculate the 4th order smoothing.

        if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

          call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,     &
     &                  idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qcicep(0,0,1,3), &
     &                  qcifrc(0,0,1,3),tmp1,tmp2,tmp3,tmp4,tmp5)

        end if

! -----

! Calculate the turbulent mixing.

        if(tubopt.ge.1) then

          call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,        &
     &                   ni,nj,nk,j31,j32,jcb,qcicep(0,0,1,3),          &
     &                   qcifrc(0,0,1,3),rkh8u,rkh8v,rkv8w,             &
     &                   tmp1,tmp2,h3,tmp3,tmp4,tmp5)

          call s_turbs(idtrnopt,idmpopt,idmfcopt,                       &
     &                 iddxiv,iddyiv,iddziv,ni,nj,nk,                   &
     &                 j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,          &
     &                 tmp1,tmp2,h3,qcifrc(0,0,1,3),tmp3,tmp4,tmp5)

        end if

! -----

! Calculate the lateral sponge damping.

        if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

          call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,ni,nj,nk,    &
     &                 rst,qcicep(0,0,1,3),rbcxy,qcifrc(0,0,1,3),tmp1)

        end if

! -----

! Calculate the vertical sponge damping.

        if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

          call s_vsps0(ksp0,ni,nj,nk,rst,qcicep(0,0,1,3),rbct,          &
     &                 qcifrc(0,0,1,3))

        end if

! -----

!! -----

!! For the hail.

        if(haiopt.eq.1) then

! Calculate the advection.

          call s_advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,       &
     &               iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,  &
     &               qcice(0,0,1,4),qcifrc(0,0,1,4),tmp1,tmp2,tmp3,tmp4)

! -----

! Calculate the 2nd order smoothing.

          if(mod(smtopt,10).eq.1) then

            call s_smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,                   &
     &                    rbr,qcicep(0,0,1,4),qcifrc(0,0,1,4),tmp1)

          end if

! -----

! Calculate the 4th order smoothing.

          if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

            call s_smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,   &
     &                   idsmhcoe,idsmvcoe,ni,nj,nk,rbr,qcicep(0,0,1,4),&
     &                   qcifrc(0,0,1,4),tmp1,tmp2,tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the turbulent mixing.

          if(tubopt.ge.1) then

            call s_turbflx(idtrnopt,idzeropt,iddxiv,iddyiv,iddziv,      &
     &                     ni,nj,nk,j31,j32,jcb,qcicep(0,0,1,4),        &
     &                     qcifrc(0,0,1,4),rkh8u,rkh8v,rkv8w,           &
     &                     tmp1,tmp2,h3,tmp3,tmp4,tmp5)

            call s_turbs(idtrnopt,idmpopt,idmfcopt,                     &
     &                   iddxiv,iddyiv,iddziv,ni,nj,nk,                 &
     &                   j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,        &
     &                   tmp1,tmp2,h3,qcifrc(0,0,1,4),tmp3,tmp4,tmp5)

          end if

! -----

! Calculate the lateral sponge damping.

          if(lspopt.ge.1.and.lspvar(7:7).eq.'o') then

            call s_lsps0(idlspopt,idwdnews,idlsnews,idlspsmt,ni,nj,nk,  &
     &                   rst,qcicep(0,0,1,4),rbcxy,qcifrc(0,0,1,4),tmp1)

          end if

! -----

! Calculate the vertical sponge damping.

          if(vspopt.ge.1.and.vspvar(7:7).eq.'o') then

            call s_vsps0(ksp0,ni,nj,nk,rst,qcicep(0,0,1,4),rbct,        &
     &                   qcifrc(0,0,1,4))

          end if

! -----

        end if

!! -----

!!! -----

      end if

!!!! -----

      end subroutine s_forceqcg

!-----7--------------------------------------------------------------7--

      end module m_forceqcg
