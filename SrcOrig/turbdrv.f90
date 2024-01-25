!***********************************************************************
      module m_turbdrv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/04/18
!     Modification: 2000/12/19, 2001/04/15, 2001/06/06, 2001/08/07,
!                   2001/11/20, 2002/01/15, 2002/04/02, 2002/08/15,
!                   2003/01/04, 2003/01/20, 2003/02/13, 2003/03/13,
!                   2003/03/21, 2003/05/19, 2003/07/15, 2003/12/12,
!                   2004/02/01, 2004/03/05, 2004/04/15, 2004/05/31,
!                   2004/06/10, 2004/09/01, 2006/01/10, 2006/02/13,
!                   2006/04/03, 2006/05/12, 2006/10/20, 2006/11/06,
!                   2007/01/20, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/01/30, 2009/02/27, 2011/11/10, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the interior procedures to solve the turbulent mixing.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bruntv
      use m_comindx
      use m_defomten
      use m_defomssq
      use m_eddyvis
      use m_strsten
      use m_turbuvw

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: turbdrv, s_turbdrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface turbdrv

        module procedure s_turbdrv

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
      subroutine s_turbdrv(fmois,dtb,ni,nj,nk,                          &
     &                     zph,j31,j32,jcb,jcb8u,jcb8v,jcb8w,           &
     &                     mf,mf8u,mf8v,rmf,rmf8u,rmf8v,pbr,ptbr,rbr,   &
     &                     up,vp,wp,ppp,ptpp,qvp,tkep,qallp,            &
     &                     ufrc,vfrc,wfrc,rkh,rkv,priv,ssq,nsq8w,       &
     &                     t13,t23,t33,t11,t22,t12,t31,t32,tmp1)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

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

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(in) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(in) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(in) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

      real, intent(in) :: qallp(0:ni+1,0:nj+1,1:nk)
                       ! Total water and ice mixing ratio at past

! Input and output variables

      real, intent(inout) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(inout) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(inout) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

! Output variables

      real, intent(out) :: rkh(0:ni+1,0:nj+1,1:nk)
                       ! rbr x horizontal eddy viscosity

      real, intent(out) :: rkv(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy viscosity

      real, intent(out) :: priv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of turbulent Prandtl number

      real, intent(out) :: ssq(0:ni+1,0:nj+1,1:nk)
                       ! Magnitude of deformation squared

      real, intent(out) :: nsq8w(0:ni+1,0:nj+1,1:nk)
                       ! Half value of Brunt Vaisala frequency squared
                       ! at w points

      real, intent(out) :: t13(0:ni+1,0:nj+1,1:nk)
                       ! x-z components of stress tensor

      real, intent(out) :: t23(0:ni+1,0:nj+1,1:nk)
                       ! y-z components of stress tensor

      real, intent(out) :: t33(0:ni+1,0:nj+1,1:nk)
                       ! z-z components of stress tensor

! Internal shared variables

      real, intent(inout) :: t11(0:ni+1,0:nj+1,1:nk)
                       ! x-x components of stress tensor

      real, intent(inout) :: t22(0:ni+1,0:nj+1,1:nk)
                       ! y-y components of stress tensor

      real, intent(inout) :: t12(0:ni+1,0:nj+1,1:nk)
                       ! x-y components of stress tensor

      real, intent(inout) :: t31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of stress tensor

      real, intent(inout) :: t32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of stress tensor

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Remark

!     rkh,rkv,priv,ssq: These variables are also temporary.

!-----7--------------------------------------------------------------7--

! Calculate the Brunt-Vaisala frequency squared.

      call s_bruntv(idcphopt,iddziv,idthresq,fmois,ni,nj,nk,            &
     &              jcb8w,pbr,ptbr,ppp,ptpp,qvp,qallp,nsq8w,            &
     &              rkh,rkv,priv,tmp1(0,0,1))

! -----

! Calculate the deformation tensor.

      call defomten(idtrnopt,idmpopt,idmfcopt,                          &
     &              idiwest,idieast,idjsouth,idjnorth,                  &
     &              iddxiv,iddyiv,iddziv,ni,nj,nk,j31,j32,              &
     &              jcb,jcb8u,jcb8v,jcb8w,mf,mf8u,mf8v,up,vp,wp,        &
     &              t11,t22,t33,t12,t13,t23,t31,t32,rkh,rkv,priv,ssq)

! -----

! Calculate the deformation squared.

      call defomssq(ni,nj,nk,t11,t22,t33,t12,t31,t32,ssq)

! -----

! Calculate the eddy viscosity and the turbulent length scale.

      call eddyvis(idmpopt,idmfcopt,idsfcopt,idtubopt,idisoopt,         &
     &             iddx,iddy,iddz,dtb,ni,nj,nk,zph,jcb,rmf,rbr,         &
     &             tkep,ssq,nsq8w,rkh,rkv,priv)

! -----

! Calculate the stress tensor.

      call strsten(idsfcopt,ni,nj,nk,ufrc,vfrc,rkh,rkv,t11,t22,t33,     &
     &             t12,t13,t23,t31,t32)

! -----

! Calculate the velocity turbulent mixing.

      call turbuvw(idtrnopt,idmpopt,idmfcopt,idoneopt,                  &
     &             iddxiv,iddyiv,iddziv,ni,nj,nk,j31,j32,               &
     &             jcb,jcb8u,jcb8v,mf,mf8u,mf8v,rmf,rmf8u,rmf8v,        &
     &             t11,t22,t33,t12,t13,t23,t31,t32,ufrc,vfrc,wfrc,tmp1)

! -----

      end subroutine s_turbdrv

!-----7--------------------------------------------------------------7--

      end module m_turbdrv
