!***********************************************************************
      module m_steppi
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/12/17
!     Modification: 1999/12/20, 2000/01/17, 2000/02/02, 2000/03/17,
!                   2000/04/18, 2000/12/18, 2001/01/15, 2001/03/13,
!                   2001/05/29, 2001/06/06, 2001/06/29, 2001/07/13,
!                   2001/08/07, 2001/11/20, 2002/04/02, 2002/06/06,
!                   2002/07/23, 2002/08/15, 2002/10/31, 2003/04/30,
!                   2003/05/19, 2003/06/27, 2003/11/05, 2003/12/12,
!                   2004/01/09, 2004/04/15, 2004/05/07, 2004/08/20,
!                   2005/01/31, 2005/02/10, 2006/06/21, 2006/09/21,
!                   2006/12/04, 2007/01/05, 2007/01/20, 2007/05/07,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/12/11,
!                   2009/02/27, 2009/03/23, 2011/09/22, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     solve the pressure to the next time step with the horizontally
!     explicit and vertically implicit method.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc4news
      use m_bcycle
      use m_combuf
      use m_comindx
      use m_exbcss
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getcname
      use m_getiname
      use m_inichar
      use m_lbcs
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_rbcss
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy
      use m_vbcp

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: steppi, s_steppi

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface steppi

        module procedure s_steppi

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
      subroutine s_steppi(fpexbvar,fpexbopt,isstp,dts,gtinc,ni,nj,nk,   &
     &                    jcb,fp,pcpx,pcpy,ppgpv,pptd,ppf)
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

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: fp(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in pressure equation

      real, intent(in) :: pcpx(1:nj,1:nk,1:2)
                       ! Phase speed of pressure
                       ! on west and east boundary

      real, intent(in) :: pcpy(1:ni,1:nk,1:2)
                       ! Phase speed of pressure
                       ! on south and north boundary

      real, intent(in) :: ppgpv(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation of GPV data
                       ! at marked time

      real, intent(in) :: pptd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! pressure perturbation of GPV data

! Input and output variable

      real, intent(inout) :: ppf(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at future

! Internal shared variables

      character(len=108) exbvar
                       ! Control flag of
                       ! extrenal boundary forced variables

      integer exbopt   ! Option for external boundary forcing

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

! Force the lateral boundary value to the external boundary value in the
! case the lateral sponge damping is performed.

      if(exbopt.ge.1.and.exbvar(4:4).ne.'x') then

        call exbcss(idexbvar,idwbc,idebc,idexnews,4,isstp,dts,gtinc,    &
     &              ni,nj,nk,pcpx,pcpy,ppgpv,pptd,ppf)

! -----

! Set the radiative lateral boundary conditions.

      else

        call rbcss(idlbcvar,idwbc,idebc,idsbc,idnbc,idnggopt,idlspopt,  &
     &             idvspopt,idlbnews,4,isstp,dts,gtinc,ni,nj,nk,        &
     &             pcpx,pcpy,ppgpv,pptd,ppf)

      end if

! -----

! Solve the pressure to the next time step.

!$omp parallel default(shared) private(k)

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          ppf(i,j,k)=ppf(i,j,k)+dts*fp(i,j,k)/jcb(i,j,k)
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

!!! Exchange the value horizontally.

!! Exchange the value horizontally between sub domain.

! In x direction.

      call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,ppf,1,1,sbuf)

      call s_shiftsx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

      call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,ppf,1,1,rbuf)

! -----

! In y direction.

      call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,ppf,1,1,sbuf)

      call s_shiftsy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

      call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,ppf,1,1,rbuf)

! -----

!! -----

!! Exchange the value horizontally between group domain.

! In x direction.

      call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,ppf,1,1,sbuf)

      call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

      call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,ppf,1,1,rbuf)

! -----

! In y direction.

      call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,ppf,1,1,sbuf)

      call s_shiftgy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

      call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,ppf,1,1,rbuf)

! -----

! In x direction again.

      call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,ppf,1,1,sbuf)

      call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

      call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,ppf,1,1,rbuf)

! -----

!! -----

!!! -----

! Set the periodic boundary conditions.

      call bcycle(idwbc,idebc,idsbc,idnbc,                              &
     &            2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,ppf)

! -----

! Set the boundary conditions at the four corners.

      call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,ni,nj,nk,ppf)

! -----

! Set the lateral boundary conditions.

      if(exbopt.eq.0.or.exbvar(4:4).eq.'x') then

        call lbcs(idwbc,idebc,idsbc,idnbc,idoneopt,ni,nj,nk,ppf)

      end if

! -----

! Set the bottom and the top boundary conditions.

      call vbcp(idbbc,ni,nj,nk,ppf)

! -----

      end subroutine s_steppi

!-----7--------------------------------------------------------------7--

      end module m_steppi
