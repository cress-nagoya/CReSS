!***********************************************************************
      module m_masscon
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/06/01
!     Modification: 2000/07/05, 2000/12/18, 2001/02/13, 2001/04/15,
!                   2001/05/29, 2001/06/06, 2001/08/07, 2002/04/02,
!                   2002/06/06, 2002/07/15, 2002/08/15, 2002/10/15,
!                   2003/01/04, 2003/02/13, 2003/03/21, 2003/04/30,
!                   2003/05/19, 2003/07/15, 2003/11/05, 2004/03/05,
!                   2004/05/31, 2004/06/10, 2004/08/20, 2005/02/10,
!                   2006/01/10, 2006/04/03, 2006/06/21, 2006/09/30,
!                   2006/11/06, 2006/12/04, 2007/01/05, 2007/01/20,
!                   2007/07/30, 2007/10/19, 2008/05/02, 2008/06/09,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     fit the x, y and z components of velocity to the mass consistent
!     equation.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc4news
      use m_bcycle
      use m_bcyclex
      use m_bcycley
      use m_chkstd
      use m_cnt2phy
      use m_combuf
      use m_comindx
      use m_comkind
      use m_commpi
      use m_copy3d
      use m_diffequa
      use m_fituvwc
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getiname
      use m_getrname
      use m_outstd11
      use m_phy2cnt
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_setcst3d
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: masscon, s_masscon

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface masscon

        module procedure s_masscon

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
      subroutine s_masscon(fpadvopt,fpsmtopt,fptubopt,                  &
     &                     fpdxiv,fpdyiv,fpdziv,fpmaseps,               &
     &                     fpalpha1,fpalpha2,ctime,ni,nj,nk,j31,j32,    &
     &                     jcb8u,jcb8v,jcb8w,mf,u,up,v,vp,w,wp,         &
     &                     wc,lamb,known,tmp1,tmp2,tmp3)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: fpdxiv
                       ! Formal parameter of unique index of dxiv

      integer, intent(in) :: fpdyiv
                       ! Formal parameter of unique index of dyiv

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: fpmaseps
                       ! Formal parameter of unique index of maseps

      integer, intent(in) :: fpalpha1
                       ! Formal parameter of unique index of alpha1

      integer, intent(in) :: fpalpha2
                       ! Formal parameter of unique index of alpha2

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: ni
                       ! Model dimensionn in x direction

      integer, intent(in) :: nj
                       ! Model dimensionn in y direction

      integer, intent(in) :: nk
                       ! Model dimensionn in z direction

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

! Input and output variables

      real, intent(inout) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(inout) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(inout) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(inout) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(inout) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(inout) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

! Internal shared variables

      integer advopt   ! Option for advection scheme
      integer smtopt   ! Option for numerical smoothing
      integer tubopt   ! Option for turbulent mixing

      integer itcnt    ! Iteration count

      integer ni_sub   ! Substitute for ni
      integer nj_sub   ! Substitute for nj

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy
      real dziv        ! Inverse of dz

      real maseps      ! Value of convergence of iteration

      real alpha1      ! Weighting coeffient
      real alpha2      ! Weighting coeffient

      real adxiv2      ! 2.0 x dxiv x alpha1 x alpha1
      real adyiv2      ! 2.0 x dyiv x alpha1 x alpha1
      real adziv2      ! 2.0 x dziv x alpha1 x alpha1

      real, intent(inout) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

      real, intent(inout) :: lamb(0:ni+1,0:nj+1,1:nk)
                       ! Lagrange multiplier

      real, intent(inout) :: known(0:ni+1,0:nj+1,1:nk)
                       ! Known quantity

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in y direction

!-----7--------------------------------------------------------------7--

! Read in the message to standard i/o.

      itcnt=0

      if(mype.eq.root) then

        call outstd11('masscon ',7,'velocity                  ',8,itcnt,&
     &                1,1,ctime)

      end if

      call chkstd(root)

! -----

! Get the required namelist variables.

      call getiname(fpadvopt,advopt)
      call getiname(fpsmtopt,smtopt)
      call getiname(fptubopt,tubopt)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)
      call getrname(fpdziv,dziv)
      call getrname(fpmaseps,maseps)
      call getrname(fpalpha1,alpha1)
      call getrname(fpalpha2,alpha2)

! -----

! Set the common used variables.

      adxiv2=2.e0*dxiv*alpha1*alpha1
      adyiv2=2.e0*dyiv*alpha1*alpha1
      adziv2=2.e0*dziv*alpha1*alpha1

! -----

! Set the substituted variables.

      ni_sub=ni
      nj_sub=nj

! -----

! Calculate the zeta components of contravariant velocity.

      call s_phy2cnt(idsthopt,idtrnopt,idmpopt,idmfcopt,idoneopt,       &
     &               ni,nj,nk,j31,j32,jcb8w,mf,up,vp,wp,wc,             &
     &               tmp1,tmp2,tmp3)

! -----

! Calculate the known quantity in the right hand.

!$omp parallel default(shared) private(k)

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-1
          tmp1(i,j,k)=jcb8u(i,j,k)*up(i,j,k)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-2
          tmp2(i,j,k)=jcb8v(i,j,k)*vp(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          tmp3(i,j,k)=jcb8w(i,j,k)*wc(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          known(i,j,k)=(tmp1(i+1,j,k)-tmp1(i,j,k))*adxiv2               &
     &      +(tmp2(i,j+1,k)-tmp2(i,j,k))*adyiv2                         &
     &      +(tmp3(i,j,k+1)-tmp3(i,j,k))*adziv2
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

! Solve the parabolic pertial differential equation.

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,lamb)

      call diffequa(iddx,iddy,iddz,maseps,alpha1,alpha2,itcnt,ni,nj,nk, &
     &              known,lamb,tmp1)

! -----

! Fit the x and the y components of velocity and the zeta components of
! contravariant velocity to the mass consistent equation by the lamb.

      call fituvwc(idiwest,idieast,idjsouth,idjnorth,iddxiv,iddyiv,     &
     &             iddziv,idalpha1,idalpha2,ni,nj,nk,jcb8u,jcb8v,jcb8w, &
     &             lamb,up,vp,wc)

! -----

! Calculate the z components of velocity.

      call s_cnt2phy(idsthopt,idtrnopt,idmpopt,idmfcopt,ni,nj,nk,       &
     &               j31,j32,jcb8w,mf,up,vp,wc,wp,tmp1,tmp2,tmp3)

! -----

!! Exchange the value in the case the 4th order calcuration is
!! performed.

      if(advopt.ge.2.or.tubopt.ge.1                                     &
     &  .or.mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

! Exchange the value between sub domain.

        call s_putbufsx(idwbc,idebc,'all',4,ni-3,ni,nj,nk,up,1,2,sbuf)
        call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,wp,2,2,sbuf)

        call s_shiftsx(idwbc,idebc,'all',nj,nk,2,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'all',0,ni+1,ni,nj,nk,up,1,2,rbuf)
        call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,wp,2,2,rbuf)

        call s_putbufsy(idsbc,idnbc,'all',4,nj-3,ni,nj,nk,vp,1,2,sbuf)
        call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,wp,2,2,sbuf)

        call s_shiftsy(idsbc,idnbc,'all',ni,nk,2,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'all',0,nj+1,ni,nj,nk,vp,1,2,rbuf)
        call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,wp,2,2,rbuf)

! -----

! Exchange the value between group domain.

        call s_putbufgx(idwbc,idebc,'all',4,ni-3,ni,nj,nk,up,1,2,sbuf)
        call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,wp,2,2,sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,nk,2,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',0,ni+1,ni,nj,nk,up,1,2,rbuf)
        call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,wp,2,2,rbuf)

        call s_putbufgy(idsbc,idnbc,'all',4,nj-3,ni,nj,nk,vp,1,2,sbuf)
        call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,wp,2,2,sbuf)

        call s_shiftgy(idsbc,idnbc,'all',ni,nk,2,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'all',0,nj+1,ni,nj,nk,vp,1,2,rbuf)
        call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,wp,2,2,rbuf)

        call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,wp,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,wp,1,1,rbuf)

! -----

! Set the boundary conditions.

        call bcyclex(idwbc,idebc,4,0,ni-3,ni+1,ni,nj,nk,up)
        call bcycley(idsbc,idnbc,4,0,nj-3,nj+1,ni,nj,nk,vp)

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,nk,wp)

! -----

! Set the boundary conditions at the four corners.

        call bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub,         &
     &               ni,nj,nk,wp)

! -----

      end if

!! -----

! Finally copy the variables.

      if(advopt.le.3) then

        call copy3d(0,ni+1,0,nj+1,1,nk,up,u)
        call copy3d(0,ni+1,0,nj+1,1,nk,vp,v)
        call copy3d(0,ni+1,0,nj+1,1,nk,wp,w)

      end if

! -----

! Read in the message to standard i/o.

      if(mype.eq.root) then

        call outstd11('masscon ',7,'velocity                  ',8,itcnt,&
     &                2,1,ctime)

      end if

      call chkstd(root)

! -----

      end subroutine s_masscon

!-----7--------------------------------------------------------------7--

      end module m_masscon
