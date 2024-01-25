!***********************************************************************
      module m_stepuv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/01/25, 1999/03/25, 1999/06/07,
!                   1999/07/05, 1999/07/21, 1999/08/03, 1999/08/09,
!                   1999/08/18, 1999/08/23, 1999/09/06, 1999/09/30,
!                   1999/10/07, 1999/10/12, 1999/10/27, 1999/11/01,
!                   1999/11/19, 1999/12/06, 2000/01/17, 2000/02/02,
!                   2000/03/17, 2000/04/18, 2000/12/18, 2001/01/15,
!                   2001/03/13, 2001/06/06, 2001/06/29, 2001/07/13,
!                   2001/08/07, 2001/09/13, 2002/04/02, 2002/06/06,
!                   2002/07/23, 2002/08/15, 2002/10/31, 2003/04/30,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2003/11/05,
!                   2003/12/12, 2004/01/09, 2004/05/07, 2004/08/20,
!                   2005/01/31, 2005/02/10, 2006/06/21, 2006/09/21,
!                   2006/09/30, 2006/12/04, 2007/01/05, 2007/01/20,
!                   2007/05/07, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2008/12/11, 2009/02/27, 2009/03/23, 2011/09/22,
!                   2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     solve the x components of velocity to the next time step.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc4news
      use m_bcycle
      use m_combuf
      use m_comindx
      use m_exbcu
      use m_exbcv
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getcname
      use m_getiname
      use m_inichar
      use m_lbcu
      use m_lbcv
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_rbcu
      use m_rbcv
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy
      use m_vbcu
      use m_vbcv

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: stepuv, s_stepuv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface stepuv

        module procedure s_stepuv

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
      subroutine s_stepuv(fpexbvar,fpexbopt,isstp,dts,                  &
     &                    gtinc,ni,nj,nk,ubr,vbr,rst8u,rst8v,           &
     &                    ufrc,vfrc,usml,vsml,ucpx,ucpy,vcpx,vcpy,      &
     &                    ugpv,utd,vgpv,vtd,uf,vf)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexbvar
                       ! Formal parameter of unique index of exbvar

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

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

      real, intent(in) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(in) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(in) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at u points

      real, intent(in) :: rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at v points

      real, intent(in) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(in) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(in) :: usml(0:ni+1,0:nj+1,1:nk)
                       ! Acoustic mood in u equation

      real, intent(in) :: vsml(0:ni+1,0:nj+1,1:nk)
                       ! Acoustic mood in v equation

      real, intent(in) :: ucpx(1:nj,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on west and east boundary

      real, intent(in) :: ucpy(1:ni,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on south and north boundary

      real, intent(in) :: vcpx(1:nj,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on west and east boundary

      real, intent(in) :: vcpy(1:ni,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on south and north boundary

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

! Input and output variables

      real, intent(inout) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(inout) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

! Internal shared variables

      character(len=108) exbvar
                       ! Control flag of
                       ! extrenal boundary forced variables

      integer exbopt   ! Option for external boundary forcing

      integer ni_sub   ! Substitute for ni
      integer nj_sub   ! Substitute for nj

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(exbvar)

! -----

! Get the required namelist variables.

      call getcname(fpexbvar,exbvar)
      call getiname(fpexbopt,exbopt)

! -----

! Set the substituted variables.

      ni_sub=ni
      nj_sub=nj

! -----

! Force the lateral boundary value to the external boundary value in the
! case the lateral sponge damping is performed.

      if(exbopt.ge.1) then

        if(exbvar(1:1).ne.'x') then

          call exbcu(idexbvar,idwbc,idebc,idexnews,idexnorm,isstp,dts,  &
     &               gtinc,ni,nj,nk,ucpx,ucpy,ugpv,utd,uf)

        end if

        if(exbvar(2:2).ne.'x') then

          call exbcv(idexbvar,idwbc,idebc,idexnews,idexnorm,isstp,dts,  &
     &               gtinc,ni,nj,nk,vcpx,vcpy,vgpv,vtd,vf)

        end if

      end if

! -----

! Set the radiative lateral boundary conditions.

      if(exbopt.eq.0.or.exbvar(1:1).eq.'x') then

        call rbcu(idlbcvar,idwbc,idebc,idsbc,idnbc,idnggopt,            &
     &            idlspopt,idvspopt,idlbnews,idlbnorm,isstp,            &
     &            dts,gtinc,ni,nj,nk,ubr,ucpx,ucpy,ugpv,utd,uf)

      end if

      if(exbopt.eq.0.or.exbvar(2:2).eq.'x') then

        call rbcv(idlbcvar,idwbc,idebc,idsbc,idnbc,idnggopt,            &
     &            idlspopt,idvspopt,idlbnews,idlbnorm,isstp,            &
     &            dts,gtinc,ni,nj,nk,vbr,vcpx,vcpy,vgpv,vtd,vf)

      end if

! -----

! Solve the x and y components of velocity to the next time step.

!$omp parallel default(shared) private(k)

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-1
          uf(i,j,k)=uf(i,j,k)+dts*(ufrc(i,j,k)+usml(i,j,k))/rst8u(i,j,k)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-2
          vf(i,j,k)=vf(i,j,k)+dts*(vfrc(i,j,k)+vsml(i,j,k))/rst8v(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

!!! Exchange the value horizontally.

!! Exchange the value horizontally between sub domain.

! In x direction.

      call s_putbufsx(idwbc,idebc,'all',3,ni-2,ni,nj,nk,uf,1,2,sbuf)
      call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,vf,2,2,sbuf)

      call s_shiftsx(idwbc,idebc,'all',nj,nk,2,sbuf,rbuf)

      call s_getbufsx(idwbc,idebc,'all',1,ni_sub,ni,nj,nk,uf,1,2,rbuf)
      call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,vf,2,2,rbuf)

! -----

! In y direction.

      call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,uf,1,2,sbuf)
      call s_putbufsy(idsbc,idnbc,'all',3,nj-2,ni,nj,nk,vf,2,2,sbuf)

      call s_shiftsy(idsbc,idnbc,'all',ni,nk,2,sbuf,rbuf)

      call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,uf,1,2,rbuf)
      call s_getbufsy(idsbc,idnbc,'all',1,nj_sub,ni,nj,nk,vf,2,2,rbuf)

! -----

!! -----

!! Exchange the value horizontally between group domain.

! In x direction.

      call s_putbufgx(idwbc,idebc,'all',3,ni-2,ni,nj,nk,uf,1,2,sbuf)
      call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,vf,2,2,sbuf)

      call s_shiftgx(idwbc,idebc,'all',nj,nk,2,sbuf,rbuf)

      call s_getbufgx(idwbc,idebc,'all',1,ni_sub,ni,nj,nk,uf,1,2,rbuf)
      call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,vf,2,2,rbuf)

! -----

! In y direction.

      call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,uf,1,2,sbuf)
      call s_putbufgy(idsbc,idnbc,'all',3,nj-2,ni,nj,nk,vf,2,2,sbuf)

      call s_shiftgy(idsbc,idnbc,'all',ni,nk,2,sbuf,rbuf)

      call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,uf,1,2,rbuf)
      call s_getbufgy(idsbc,idnbc,'all',1,nj_sub,ni,nj,nk,vf,2,2,rbuf)

! -----

! In x direction again.

      call s_putbufgx(idwbc,idebc,'all',3,ni-2,ni,nj,nk,uf,1,2,sbuf)
      call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,vf,2,2,sbuf)

      call s_shiftgx(idwbc,idebc,'all',nj,nk,2,sbuf,rbuf)

      call s_getbufgx(idwbc,idebc,'all',1,ni_sub,ni,nj,nk,uf,1,2,rbuf)
      call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,vf,2,2,rbuf)

! -----

!! -----

!!! -----

! Set the periodic boundary conditions.

      call bcycle(idwbc,idebc,idsbc,idnbc,                              &
     &            3,1,ni-2,ni_sub,2,1,nj-2,nj-1,ni,nj,nk,uf)

      call bcycle(idwbc,idebc,idsbc,idnbc,                              &
     &            2,1,ni-2,ni-1,3,1,nj-2,nj_sub,ni,nj,nk,vf)

! -----

! Set the boundary conditions at the four corners.

      call bc4news(idwbc,idebc,idsbc,idnbc,1,ni_sub,1,nj-1,ni,nj,nk,uf)
      call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj_sub,ni,nj,nk,vf)

! -----

! Set the lateral boundary conditions.

      if(exbopt.eq.0.or.exbvar(1:1).eq.'x') then

        call lbcu(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,uf)

      end if

      if(exbopt.eq.0.or.exbvar(2:2).eq.'x') then

        call lbcv(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,vf)

      end if

! -----

! Set the bottom and the top boundary conditions.

      call vbcu(ni,nj,nk,uf)
      call vbcv(ni,nj,nk,vf)

! -----

      end subroutine s_stepuv

!-----7--------------------------------------------------------------7--

      end module m_stepuv
