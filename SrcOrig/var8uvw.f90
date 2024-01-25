!***********************************************************************
      module m_var8uvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/01/25, 1999/04/06, 1999/07/05,
!                   1999/08/03, 1999/08/18, 1999/08/23, 1999/09/16,
!                   1999/09/30, 1999/10/07, 1999/10/12, 1999/11/01,
!                   1999/11/19, 2000/01/17, 2000/04/18, 2000/12/18,
!                   2001/11/20, 2002/04/02, 2002/06/06, 2003/04/30,
!                   2003/05/19, 2003/11/05, 2005/02/10, 2006/09/30,
!                   2006/12/04, 2007/01/05, 2007/01/20, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     optional variable be averaged to u, v and w points.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc8u
      use m_bc8v
      use m_bc8w
      use m_bcyclex
      use m_bcycley
      use m_combuf
      use m_comindx
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getiname
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

      public :: var8uvw, s_var8uvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface var8uvw

        module procedure s_var8uvw

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
      subroutine s_var8uvw(fpwbc,fpebc,fpexbopt,                        &
     &                     ni,nj,nk,var,var8u,var8v,var8w)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: var(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable

! Output variables

      real, intent(out) :: var8u(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at u points

      real, intent(out) :: var8v(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at v points

      real, intent(out) :: var8w(0:ni+1,0:nj+1,1:nk)
                       ! Optional variable at w points

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer exbopt   ! Option for external boundary forcing

      integer ni_sub   ! Substitute for ni
      integer nj_sub   ! Substitute for nj

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpexbopt,exbopt)

! -----

! Set the substituted variables.

      ni_sub=ni
      nj_sub=nj

! -----

! Be averaged to u, v and w points.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=1,ni
          var8u(i,j,k)=.5e0*(var(i-1,j,k)+var(i,j,k))
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1,nj
        do i=0,ni
          var8v(i,j,k)=.5e0*(var(i,j-1,k)+var(i,j,k))
        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          var8w(i,j,k)=.5e0*(var(i,j,k-1)+var(i,j,k))
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

!! Set the lateral boundary conditions.

      if(exbopt.eq.0) then

! Set the west and east boundary conditions.

        call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,var8u,1,1,    &
     &                  sbuf)

        call s_shiftsx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,var8u,1,1,  &
     &                  rbuf)

        call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,var8u,1,1,    &
     &                  sbuf)

        call s_shiftgx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,var8u,1,1,  &
     &                  rbuf)

        call bcyclex(idwbc,idebc,3,1,ni-2,ni_sub,ni,nj,nk,var8u)

        call bc8u(idwbc,idebc,ni,nj,nk,var8u)

! -----

! Set the south and north boundary conditions.

        call s_putbufsy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,var8v,1,1,    &
     &                  sbuf)

        call s_shiftsy(idsbc,idnbc,'bnd',ni,nk,1,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,var8v,1,1,  &
     &                  rbuf)

        call s_putbufgy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,var8v,1,1,    &
     &                  sbuf)

        call s_shiftgy(idsbc,idnbc,'bnd',ni,nk,1,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,var8v,1,1,  &
     &                  rbuf)

        call bcycley(idsbc,idnbc,3,1,nj-2,nj_sub,ni,nj,nk,var8v)

        call bc8v(idsbc,idnbc,ni,nj,nk,var8v)

! -----

      end if

!! -----

! Set the periodic boundary conditions.

      if(exbopt.ge.1.and.abs(wbc).eq.1.and.abs(ebc).eq.1) then

        call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,var8u,1,1,    &
     &                  sbuf)

        call s_shiftsx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,var8u,1,1,  &
     &                  rbuf)

        call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,var8u,1,1,    &
     &                  sbuf)

        call s_shiftgx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,var8u,1,1,  &
     &                  rbuf)

        call bcyclex(idwbc,idebc,3,1,ni-2,ni_sub,ni,nj,nk,var8u)

        call bc8u(idwbc,idebc,ni,nj,nk,var8u)

      end if

! -----

! Set the bottom and the top boundary conditions.

      call bc8w(idbbc,idtbc,ni,nj,nk,var8w)

! -----

      end subroutine s_var8uvw

!-----7--------------------------------------------------------------7--

      end module m_var8uvw
