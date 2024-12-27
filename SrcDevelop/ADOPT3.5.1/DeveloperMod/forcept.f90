!***********************************************************************
      module m_forcept
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/05/20,
!                   1999/06/07, 1999/07/05, 1999/07/21, 1999/08/03,
!                   1999/08/09, 1999/08/18, 1999/08/23, 1999/09/16,
!                   1999/09/30, 1999/10/07, 1999/10/12, 1999/10/22,
!                   1999/11/01, 1999/11/19, 1999/11/24, 2000/01/17,
!                   2000/02/07, 2000/04/18, 2000/12/19, 2001/03/13,
!                   2001/04/15, 2001/05/29, 2001/06/06, 2001/07/13,
!                   2001/08/07, 2001/11/20, 2002/04/02, 2002/06/18,
!                   2002/07/23, 2002/08/15, 2002/09/09, 2002/10/31,
!                   2002/12/11, 2003/01/04, 2003/01/20, 2003/03/13,
!                   2003/03/21, 2003/04/30, 2003/05/19, 2003/07/15,
!                   2003/09/01, 2003/10/10, 2003/11/28, 2003/12/12,
!                   2004/02/01, 2004/03/05, 2004/04/15, 2004/05/31,
!                   2004/06/10, 2004/08/01, 2004/08/20, 2006/01/10,
!                   2006/02/13, 2006/04/03, 2006/05/12, 2006/06/21,
!                   2006/09/21, 2006/11/06, 2007/05/07, 2007/07/30,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/12/11,
!                   2009/02/27, 2009/03/23, 2011/09/22, 2013/01/28,
!                   2013/02/13, 2013/03/27

!     Author      : Satoki Tsujino
!     Modification: 2024/12/26

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the forcing term in the potential temperature equation.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_advbspt
      use m_advs
      use m_comindx
      use m_comkind
      use m_comtub
      use m_dmptub
      use m_getcname
      use m_getiname
      use m_inichar
      use m_lsps
      use m_nlsms
      use m_s2gpv
      use m_smoo2s
      use m_smoo4s
      use m_turbflx
      use m_turbs
      use m_vsps

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: forcept, s_forcept

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface forcept

        module procedure s_forcept

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
      subroutine s_forcept(fpnggvar,fplspvar,fpvspvar,                  &
     &                     fplspopt,fpvspopt,fpgwmopt,fpsmtopt,         &
     &                     fpadvopt,fptubopt,ksp0,nggdmp,gtinc,         &
     &                     ni,nj,nk,j31,j32,jcb,jcb8u,jcb8v,mf,rmf,     &
     &                     rmf8u,rmf8v,ptbr,rbr,rst,rstxu,rstxv,rstxwc, &
     &                     w,wp,ptp,ptpp,rkh8u,rkh8v,rkv8w,rbcxy,rbct,  &
     &                     ptpgpv,ptptd,ptfrc,pt,h3,tmp1,tmp2,          &
     &                     tmp3,tmp4,tmp5)
!***********************************************************************

      use m_comtub

! Input variables

      integer, intent(in) :: fpnggvar
                       ! Formal parameter of unique index of nggvar

      integer, intent(in) :: fplspvar
                       ! Formal parameter of unique index of lspvar

      integer, intent(in) :: fpvspvar
                       ! Formal parameter of unique index of vspvar

      integer, intent(in) :: fplspopt
                       ! Formal parameter of unique index of lspopt

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: fpgwmopt
                       ! Formal parameter of unique index of gwmopt

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

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

      real, intent(in) :: nggdmp
                       ! Analysis nudging damping coefficient for GPV

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

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

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

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

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(in) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at present

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

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

      real, intent(in) :: ptpgpv(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation of GPV data
                       ! at marked time

      real, intent(in) :: ptptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! potential temperature perturbation of GPV data

! Input and output variable

      real, intent(inout) :: ptfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in potential temperature equation

! Internal shared variables

      character(len=108) nggvar
                       ! Control flag of
                       ! analysis nudged variables to GPV

      character(len=108) lspvar
                       ! Control flag of
                       ! lateral sponge damped variables

      character(len=108) vspvar
                       ! Control flag of
                       ! vertical sponge damped variables

      integer lspopt   ! Option for lateral sponge damping
      integer vspopt   ! Option for vertical sponge damping
      integer gwmopt   ! Option for gravity wave mode integration

      integer smtopt   ! Option for numerical smoothing
      integer advopt   ! Option for advection scheme
      integer tubopt   ! Option for turbulent mixing

      real, intent(inout) :: pt(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature

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

      character(len=108) dmpvar
                       ! Control flag of dump variables

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(nggvar)
      call inichar(lspvar)
      call inichar(vspvar)

! -----

! Get the required namelist variables.

      call getcname(fpnggvar,nggvar)
      call getcname(fplspvar,lspvar)
      call getcname(fpvspvar,vspvar)
      call getiname(fplspopt,lspopt)
      call getiname(fpvspopt,vspopt)
      call getiname(fpgwmopt,gwmopt)
      call getiname(fpsmtopt,smtopt)
      call getiname(fpadvopt,advopt)
      call getiname(fptubopt,tubopt)

      call getcname(iddmpvar,dmpvar)

! -----

! Get the potential temperature.

      if(tubopt.ge.1) then

!$omp parallel default(shared) private(k)

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            pt(i,j,k)=ptbr(i,j,k)+ptpp(i,j,k)
          end do
          end do

!$omp end do

        end do

!$omp end parallel

      end if

! -----

! Calculate the advection.

      call advs(idadvopt,idiwest,idieast,idjsouth,idjnorth,             &
     &          iddxiv,iddyiv,iddziv,ni,nj,nk,rstxu,rstxv,rstxwc,       &
     &          ptp,ptfrc,tmp1,tmp2,tmp3,tmp4)

! -----

      if(dmpvar(16:16).eq.'+')then

        call s_dmptub( ni, nj, nk, ptfrc, 'numdpt', 'o', 1.0e0 )

      end if

! Calculate the 2nd order smoothing.

      if(mod(smtopt,10).eq.1) then

        call smoo2s(idsmhcoe,idsmvcoe,ni,nj,nk,rbr,ptpp,ptfrc,tmp1)

      end if

! -----

! Calculate the 4th order smoothing.

      if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

        call smoo4s(idsmtopt,idiwest,idieast,idjsouth,idjnorth,         &
     &              idsmhcoe,idsmvcoe,ni,nj,nk,rbr,ptpp,ptfrc,          &
     &              tmp1,tmp2,tmp3,tmp4,tmp5)

      end if

! -----

! Calculate the non linear smoothing.

      if(smtopt.ge.11) then

        call nlsms(idnlhcoe,idnlvcoe,1.e0,ni,nj,nk,rbr,ptpp,ptfrc,      &
     &             tmp1,tmp2,tmp3,tmp4)

      end if

      if(dmpvar(16:16).eq.'+')then

        call s_dmptub( ni, nj, nk, ptfrc, 'numdpt', 'm', 1.0e0,         &
     &                 fact3d=rst(0:ni+1,0:nj+1,1:nk) )

      end if

! -----

! Calculate the turbulent mixing.

      if(tubopt.ge.1) then

        if(dmpvar(16:16).eq.'o'.or.dmpvar(16:16).eq.'+')then

          call s_dmptub( ni, nj, nk, ptfrc, 'turbpt', 'o', 1.0e0 )

        end if

        call turbflx(idtrnopt,idsfcopt,iddxiv,iddyiv,iddziv,ni,nj,nk,   &
     &               j31,j32,jcb,pt,ptfrc,rkh8u,rkh8v,rkv8w,            &
     &               tmp1,tmp2,h3,tmp3,tmp4,tmp5)

        call turbs(idtrnopt,idmpopt,idmfcopt,iddxiv,iddyiv,iddziv,      &
     &             ni,nj,nk,j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,     &
     &             tmp1,tmp2,h3,ptfrc,tmp3,tmp4,tmp5)

        if(dmpvar(16:16).eq.'o'.or.dmpvar(16:16).eq.'+')then

          call s_dmptub( ni, nj, nk, ptfrc, 'turbpt', 'm', 1.0e0,         &
     &                   fact3d=rst(0:ni+1,0:nj+1,1:nk) )

        end if

      end if

! -----

! Perform the analysis nudging to GPV.

      if(nggdmp.gt.0.e0.and.nggvar(5:5).eq.'o') then

        call s2gpv(idgpvvar,9,nggdmp,gtinc,                             &
     &             ni,nj,nk,rst,ptpp,ptpgpv,ptptd,ptfrc)

      end if

! -----

! Calculate the lateral sponge damping.

      if(lspopt.ge.1.and.lspvar(5:5).eq.'o') then

        call lsps(idgpvvar,idlspopt,idwdnews,idlsnews,idlspsmt,9,gtinc, &
     &            ni,nj,nk,rst,ptpp,rbcxy,ptpgpv,ptptd,ptfrc,tmp1)

      end if

! -----

! Calculate the vertical sponge damping.

      if(vspopt.ge.1.and.vspvar(5:5).eq.'o') then

        call vsps(idgpvvar,idvspopt,9,ksp0,gtinc,ni,nj,nk,rst,ptpp,rbct,&
     &            ptpgpv,ptptd,ptfrc)

      end if

! -----

! Calculate the base state advection.

      if(gwmopt.eq.0) then

        if(advopt.le.3) then

          call advbspt(idgwmopt,iddziv,ni,nj,nk,ptbr,rbr,w,ptfrc,tmp1)

        else

          call advbspt(idgwmopt,iddziv,ni,nj,nk,ptbr,rbr,wp,ptfrc,tmp1)

        end if

      end if

! -----

      end subroutine s_forcept

!-----7--------------------------------------------------------------7--

      end module m_forcept
