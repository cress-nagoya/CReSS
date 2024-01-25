!***********************************************************************
      module m_stepwi
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/12/20
!     Modification: 2000/01/17, 2000/02/02, 2000/03/17, 2000/03/23,
!                   2000/04/18, 2000/06/01, 2000/07/05, 2000/12/18,
!                   2001/01/04, 2001/01/15, 2001/02/24, 2001/03/13,
!                   2001/04/15, 2001/05/29, 2001/06/06, 2001/06/29,
!                   2001/07/13, 2001/08/07, 2001/09/13, 2001/10/18,
!                   2002/04/02, 2002/06/06, 2002/07/23, 2002/08/15,
!                   2002/10/31, 2003/01/04, 2003/03/21, 2003/04/30,
!                   2003/05/19, 2003/06/27, 2003/10/31, 2003/11/05,
!                   2003/12/12, 2004/01/09, 2004/03/05, 2004/05/07,
!                   2004/06/10, 2004/07/01, 2004/08/20, 2005/01/31,
!                   2005/02/10, 2005/08/05, 2006/02/13, 2006/04/03,
!                   2006/06/21, 2006/09/21, 2006/11/06, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2007/05/07, 2007/10/19,
!                   2008/05/02, 2008/06/09, 2008/08/25, 2008/12/11,
!                   2009/02/27, 2009/03/23, 2011/09/22, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     solve the z components of velocity to the next time step with the
!     horizontally explicit and vertically implicit method.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bbcw
      use m_bc4news
      use m_bcycle
      use m_combuf
      use m_comindx
      use m_comphy
      use m_exbcw
      use m_gaussel
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getcname
      use m_getiname
      use m_getrname
      use m_gseidel
      use m_inichar
      use m_lbcw
      use m_phy2cnt
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_rbcw
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy
      use m_vbcw

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: stepwi, s_stepwi

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface stepwi

        module procedure s_stepwi

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
      subroutine s_stepwi(fpexbvar,fpexbopt,fpimpopt,fpbuyopt,          &
     &                    fpdziv,fpweicoe,isstp,dts,gtinc,ni,nj,nk,     &
     &                    j31,j32,jcb,jcb8w,mf,rbr,rst,rst8w,rcsq,uf,vf,&
     &                    wcpx,wcpy,wgpv,wtd,fw,wf,wc,tmp1,tmp2,tmp3)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexbvar
                       ! Formal parameter of unique index of exbvar

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: fpimpopt
                       ! Formal parameter of unique index of impopt

      integer, intent(in) :: fpbuyopt
                       ! Formal parameter of unique index of buyopt

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: fpweicoe
                       ! Formal parameter of unique index of weicoe

      integer, intent(in) :: isstp
                       ! Index of small time steps integration

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: rst8w(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at w points

      real, intent(in) :: rcsq(0:ni+1,0:nj+1,1:nk)
                       ! rbr x sound wave speed squared

      real, intent(in) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(in) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(in) :: wcpx(1:nj,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on west and east boundary

      real, intent(in) :: wcpy(1:ni,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on south and north boundary

      real, intent(in) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: wtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! z components of velocity of GPV data

! Input and output variables

      real, intent(inout) :: fw(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

      real, intent(inout) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

! Output variable

      real, intent(out) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

! Internal shared variables

      character(len=108) exbvar
                       ! Control flag of
                       ! extrenal boundary forced variables

      integer exbopt   ! Option for external boundary forcing
      integer impopt   ! Option for vertical implicit method
      integer buyopt   ! Option for buoyancy calculation

      real dziv        ! Inverse of dz

      real weicoe      ! Weighting coefficient for impilicit method

      real g05         ! 0.5 x g

      real rstiv       ! 1.0 / rst

      real sbsqzi      ! dts^2 x weicoe^2 x dziv
      real sbsqg5      ! 0.5 x g x dts^2 x weicoe^2

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real a           ! Temporary variable
      real b           ! Temporary variable

      real mm          ! Temporary variable
      real nn          ! Temporary variable

! Remark

!     fw: This variable is also temporary, because it is not used again.

!     wc: This variable is also temporary.

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(exbvar)

! -----

! Get the required namelist variables.

      call getcname(fpexbvar,exbvar)
      call getiname(fpexbopt,exbopt)
      call getiname(fpimpopt,impopt)
      call getiname(fpbuyopt,buyopt)
      call getrname(fpdziv,dziv)
      call getrname(fpweicoe,weicoe)

! -----

! Force the lateral boundary value to the external boundary value in the
! case the lateral sponge damping is performed.

      if(exbopt.ge.1.and.exbvar(3:3).ne.'x') then

        call exbcw(idexbvar,idwbc,idebc,idexnews,isstp,dts,gtinc,       &
     &             ni,nj,nk,wcpx,wcpy,wgpv,wtd,wf)

! -----

! Set the radiative lateral boundary conditions.

      else

        call rbcw(idgpvvar,idlbcvar,idwbc,idebc,idsbc,idnbc,            &
     &            idnggopt,idlspopt,idvspopt,idlbnews,isstp,            &
     &            dts,gtinc,ni,nj,nk,wcpx,wcpy,wgpv,wtd,wf)

      end if

! -----

!! Solve the z components of velocity to the next time step.

! Set the common used variables.

      g05=.5e0*g

      sbsqzi=dts*dts*weicoe*weicoe*dziv
      sbsqg5=.5e0*g*dts*dts*weicoe*weicoe

! -----

! Calculate the solved vector.

!$omp parallel default(shared) private(k)

      do k=3,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          wf(i,j,k)=wf(i,j,k)+dts*fw(i,j,k)/rst8w(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          tmp1(i,j,k)=g05*rbr(i,j,k)/rcsq(i,j,k)
          tmp2(i,j,k)=dziv/jcb(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          fw(i,j,k)=tmp1(i,j,k)+tmp2(i,j,k)
          wc(i,j,k)=tmp1(i,j,k)-tmp2(i,j,k)
        end do
        end do

!$omp end do

      end do

      if(buyopt.eq.0) then

        do k=3,nk-2

!$omp do schedule(runtime) private(i,j,mm,nn,a)

          do j=2,nj-2
          do i=2,ni-2
            a=sbsqzi/rst8w(i,j,k)

            mm=a*rcsq(i,j,k-1)
            nn=a*rcsq(i,j,k)

            tmp1(i,j,k)=-mm*fw(i,j,k-1)

            tmp2(i,j,k)=1.e0+(nn*fw(i,j,k)-mm*wc(i,j,k-1))

            tmp3(i,j,k)=nn*wc(i,j,k)

          end do
          end do

!$omp end do

        end do

      else if(buyopt.eq.1) then

        do k=3,nk-2

!$omp do schedule(runtime) private(i,j,rstiv,mm,nn,a,b)

          do j=2,nj-2
          do i=2,ni-2
            rstiv=1.e0/rst8w(i,j,k)

            a=rstiv*sbsqzi
            b=rstiv*sbsqg5

            mm=b*rst(i,j,k-1)-a*rcsq(i,j,k-1)
            nn=b*rst(i,j,k)+a*rcsq(i,j,k)

            tmp1(i,j,k)=mm*fw(i,j,k-1)

            tmp2(i,j,k)=1.e0+(mm*wc(i,j,k-1)+nn*fw(i,j,k))

            tmp3(i,j,k)=nn*wc(i,j,k)

          end do
          end do

!$omp end do

        end do

      end if

!$omp end parallel

! -----

! Set the bottom boundary conditions.

      call s_bbcw(idmpopt,idmfcopt,ni,nj,nk,j31,j32,mf,tmp1,            &
     &            uf,vf,wf,fw,wc)

! -----

! Solve the trifiagonal equation with the Gauss elimination.

      if(impopt.le.2) then

        call gaussel(idimpopt,2,ni-2,2,nj-2,3,nk-2,ni,nj,nk,            &
     &               tmp1,tmp2,tmp3,wf,fw)

! -----

! Solve the trifiagonal equation with the Gauss-Seidel method.

      else if(impopt.eq.3) then

        call s_gseidel(idgsdeps,ni,nj,nk,tmp1,tmp2,tmp3,wf,fw,          &
     &                 wc(0,0,1),wc(0,0,2),wc(0,0,3))

      end if

! -----

!! -----

!!! Exchange the value horizontally.

!! Exchange the value horizontally between sub domain.

! In x direction.

      call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,wf,1,1,sbuf)

      call s_shiftsx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

      call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,wf,1,1,rbuf)

! -----

! In y direction.

      call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,wf,1,1,sbuf)

      call s_shiftsy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

      call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,wf,1,1,rbuf)

! -----

!! -----

!! Exchange the value horizontally between group domain.

! In x direction.

      call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,wf,1,1,sbuf)

      call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

      call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,wf,1,1,rbuf)

! -----

! In y direction.

      call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,wf,1,1,sbuf)

      call s_shiftgy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

      call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,wf,1,1,rbuf)

! -----

! In x direction again.

      call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,wf,1,1,sbuf)

      call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

      call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,wf,1,1,rbuf)

! -----

!! -----

!!! -----

! Set the periodic boundary conditions.

      call bcycle(idwbc,idebc,idsbc,idnbc,                              &
     &            2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,wf)

! -----

! Set the boundary conditions at the four corners.

      call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,ni,nj,nk,wf)

! -----

! Set the lateral boundary conditions.

      if(exbopt.eq.0.or.exbvar(3:3).eq.'x') then

        call lbcw(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,wf)

      end if

! -----

! Calculate the zeta components of contravariant velocity.

      call s_phy2cnt(idsthopt,idtrnopt,idmpopt,idmfcopt,idoneopt,       &
     &               ni,nj,nk,j31,j32,jcb8w,mf,uf,vf,wf,wc,             &
     &               tmp1,tmp2,tmp3)

! -----

! Set the bottom and the top boundary conditions.

      call s_vbcw(idbbc,idtbc,idmpopt,idmfcopt,ni,nj,nk,j31,j32,jcb8w,  &
     &            mf,uf,vf,wc,wf,tmp1,tmp2,tmp3)

! -----

      end subroutine s_stepwi

!-----7--------------------------------------------------------------7--

      end module m_stepwi
