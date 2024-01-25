!***********************************************************************
      module m_initke
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/03/13
!     Modification: 2003/04/30, 2003/05/19, 2003/10/31, 2003/12/12,
!                   2004/03/05, 2004/03/22, 2004/04/15, 2004/05/07,
!                   2004/05/31, 2004/08/20, 2005/02/10, 2005/10/05,
!                   2006/01/10, 2006/04/03, 2006/09/30, 2006/11/06,
!                   2006/12/04, 2007/01/05, 2007/01/20, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2011/08/09,
!                   2011/11/10, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     initialize the turbulent kinetic energy.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc4news
      use m_bcycle
      use m_combuf
      use m_comindx
      use m_commath
      use m_comphy
      use m_copy3d
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getiname
      use m_getrname
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
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

      public :: initke, s_initke

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface initke

        module procedure s_initke

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_initke(fpwbc,fpebc,fpsbc,fpnbc,fpmpopt,fpmfcopt,     &
     &                    fpadvopt,fpsmtopt,fpisoopt,fpdx,fpdy,fpdz,    &
     &                    ni,nj,nk,jcb,rmf,rbr,rst,rkv,tke,tkep)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: fpisoopt
                       ! Formal parameter of unique index of isoopt

      integer, intent(in) :: fpdx
                       ! Formal parameter of unique index of dx

      integer, intent(in) :: fpdy
                       ! Formal parameter of unique index of dy

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: rkv(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy viscosity

! Output variables

      real, intent(out) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at present

      real, intent(out) :: tkep(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy at past

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor
      integer advopt   ! Option for advection scheme
      integer smtopt   ! Option for numerical smoothing
      integer isoopt   ! Option for grid shape

      integer ni_sub   ! Substitute for ni
      integer nj_sub   ! Substitute for nj

      real dx          ! Grid distance in x direction
      real dy          ! Grid distance in y direction
      real dz          ! Grid distance in z direction

      real ckmds3      ! ckm^3 x dx x dy x dz
      real ckmdz       ! ckm x dz

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real sqrtke      ! Square root of turbulent kinetic energy

!-----7--------------------------------------------------------------7--

!!! Initialize the turbulent kinetic energy.

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getiname(fpadvopt,advopt)
      call getiname(fpsmtopt,smtopt)
      call getiname(fpisoopt,isoopt)
      call getrname(fpdx,dx)
      call getrname(fpdy,dy)
      call getrname(fpdz,dz)

! -----

! Set the common used variables.

      ckmds3=ckm*ckm*ckm*dx*dy*dz
      ckmdz=ckm*dz

! -----

! Set the substituted variables.

      ni_sub=ni
      nj_sub=nj

! -----

!! Set the initial turbulent kinetic energy.

!$omp parallel default(shared) private(k)

! Isotropic case.

      if(isoopt.eq.1) then

        if(mfcopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,sqrtke)

            do j=1,nj-1
            do i=1,ni-1
              sqrtke=rkv(i,j,k)                                         &
     &          /(rbr(i,j,k)*exp(oned3*log(ckmds3*jcb(i,j,k))))

              tkep(i,j,k)=sqrtke*sqrtke

            end do
            end do

!$omp end do

          end do

        else

          if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,sqrtke)

              do j=1,nj-1
              do i=1,ni-1
                sqrtke=rkv(i,j,k)/(rbr(i,j,k)                           &
     &            *exp(oned3*log(ckmds3*rmf(i,j,2)*jcb(i,j,k))))

                tkep(i,j,k)=sqrtke*sqrtke

              end do
              end do

!$omp end do

            end do

          else

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,sqrtke)

              do j=1,nj-1
              do i=1,ni-1
                sqrtke=rkv(i,j,k)/(rbr(i,j,k)                           &
     &            *exp(oned3*log(ckmds3*rmf(i,j,3)*jcb(i,j,k))))

                tkep(i,j,k)=sqrtke*sqrtke

              end do
              end do

!$omp end do

            end do

          end if

        end if

! -----

! Anisotropic case.

      else if(isoopt.eq.2) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,sqrtke)

          do j=1,nj-1
          do i=1,ni-1
            sqrtke=rkv(i,j,k)/(ckmdz*rst(i,j,k))

            tkep(i,j,k)=sqrtke*sqrtke

          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!! -----

!! Exchange the value horizontally or set the periodic boundary
!! conditions.

      if(mod(wbc,10).eq.6.or.mod(ebc,10).eq.6.or.                       &
     &   mod(sbc,10).eq.6.or.mod(nbc,10).eq.6.or.                       &
     &   advopt.ge.2.or.mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

! Exchange the value horizontally between sub domain.

        call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,tkep,1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,tkep,1,1,   &
     &                  rbuf)

        call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,tkep,1,1,sbuf)

        call s_shiftsy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,tkep,1,1,   &
     &                  rbuf)

! -----

! Exchange the value horizontally between group domain.

        call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,tkep,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,tkep,1,1,   &
     &                  rbuf)

        call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,tkep,1,1,sbuf)

        call s_shiftgy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,tkep,1,1,   &
     &                  rbuf)

        call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,tkep,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,tkep,1,1,   &
     &                  rbuf)

! -----

! Set the periodic boundary conditions.

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,nk,tkep)

! -----

! Set the boundary conditions at the four corners.

        call bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub,         &
     &               ni,nj,nk,tkep)

! -----

      end if

!! -----

! Copy the past variables to the present.

      if(advopt.le.3) then

        call copy3d(0,ni+1,0,nj+1,1,nk,tkep,tke)

      end if

! -----

!!! -----

      end subroutine s_initke

!-----7--------------------------------------------------------------7--

      end module m_initke
