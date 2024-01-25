!***********************************************************************
      module m_pblu
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/16
!     Modification: 2001/11/20, 2001/12/10, 2002/01/15, 2002/04/02,
!                   2002/08/15, 2002/12/27, 2003/01/20, 2003/04/30,
!                   2003/05/19, 2003/11/05, 2003/12/12, 2005/02/10,
!                   2006/04/03, 2006/09/30, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/09/22, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the virtical diffusion for the x components of velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc8u
      use m_bcyclex
      use m_combuf
      use m_comindx
      use m_gaussel
      use m_getbufgx
      use m_getbufsx
      use m_getiname
      use m_getrname
      use m_putbufgx
      use m_putbufsx
      use m_shiftgx
      use m_shiftsx
      use m_vbcu

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: pblu, s_pblu

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface pblu

        module procedure s_pblu

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
      subroutine s_pblu(fplevpbl,fpdziv,dtb,ni,nj,nk,jcb8w,rbr,rst8u,   &
     &                  ce,kms,u,rr,ss,tt,tmp1,tmp2)
!***********************************************************************

! Input variables

      integer, intent(in) :: fplevpbl
                       ! Formal parameter of unique index of levpbl

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at u points

      real, intent(in) :: ce(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface momentum flux

      real, intent(in) :: kms(0:ni+1,0:nj+1,1:nk)
                       ! Eddy viscosity in planetaty boundary layer

! Input and output variable

      real, intent(inout) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

! Internal shared variables

      integer levpbl   ! Number of planetary boundary layer

      integer levp1    ! levpbl + 1

      integer ni_sub   ! Substitute for ni

      real dziv        ! Inverse of grid distance in z direction

      real dtdzv       ! 0.5 x dtb x dziv
      real dtdzv2      ! 0.5 x dtb x (dziv^2)

      real, intent(inout) :: rr(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient matrix

      real, intent(inout) :: ss(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient matrix

      real, intent(inout) :: tt(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient matrix

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real a           ! Temporary variable
      real b           ! Temporary variable
      real c           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fplevpbl,levpbl)
      call getrname(fpdziv,dziv)

! -----

! Set the common used variables.

      levp1=levpbl+1

      dtdzv=.5e0*dtb*dziv
      dtdzv2=.5e0*dtb*dziv*dziv

! -----

! Set the substituted variable.

      ni_sub=ni

! -----

!! Calculate the virtical diffusion for the x components of velocity.

! Set the coefficient matrix.

!$omp parallel default(shared) private(k)

      if(levpbl.eq.1) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=2,ni-1
          ss(i,j,2)=1.e0+(ce(i-1,j)+ce(i,j))/rst8u(i,j,2)*dtdzv
        end do
        end do

!$omp end do

      else if(levpbl.ge.2) then

        do k=3,levpbl+1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            tmp1(i,j,k)=(rbr(i,j,k-1)+rbr(i,j,k))*kms(i,j,k)*dtdzv2
          end do
          end do

!$omp end do

        end do

        do k=3,levpbl+1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=2,ni-1
            tmp2(i,j,k)=-(tmp1(i-1,j,k)+tmp1(i,j,k))                    &
     &        /(jcb8w(i-1,j,k)+jcb8w(i,j,k))
          end do
          end do

!$omp end do

        end do

!$omp do schedule(runtime) private(i,j,a)

        do j=1,nj-1
        do i=2,ni-1
          a=1.e0/rst8u(i,j,2)

          rr(i,j,2)=0.e0
          ss(i,j,2)=1.e0+a*((ce(i-1,j)+ce(i,j))*dtdzv-tmp2(i,j,3))
          tt(i,j,2)=a*tmp2(i,j,3)

        end do
        end do

!$omp end do

        if(levpbl.ge.3) then

          do k=3,levpbl

!$omp do schedule(runtime) private(i,j,a,b,c)

            do j=1,nj-1
            do i=2,ni-1
              a=1.e0/rst8u(i,j,k)

              b=a*tmp2(i,j,k)
              c=a*tmp2(i,j,k+1)

              rr(i,j,k)=b
              ss(i,j,k)=1.e0-(b+c)
              tt(i,j,k)=c

            end do
            end do

!$omp end do

          end do

        end if

!$omp do schedule(runtime) private(i,j,a)

        do j=1,nj-1
        do i=2,ni-1
          a=tmp2(i,j,levp1)/rst8u(i,j,levp1)

          rr(i,j,levp1)=a
          ss(i,j,levp1)=1.e0-a
          tt(i,j,levp1)=0.e0

        end do
        end do

!$omp end do

      end if

!$omp end parallel

! -----

! Set the lateral boundary conditions for coefficient matrix.

      if(levpbl.eq.1) then

        call s_putbufsx(idwbc,idebc,'all',3,ni-2,ni,nj,1,               &
     &                  ss(0,0,2),1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'all',nj,1,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'all',1,ni_sub,ni,nj,1,             &
     &                  ss(0,0,2),1,1,rbuf)

        call s_putbufgx(idwbc,idebc,'all',3,ni-2,ni,nj,1,               &
     &                  ss(0,0,2),1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,1,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',1,ni_sub,ni,nj,1,             &
     &                  ss(0,0,2),1,1,rbuf)

        call s_bcyclex(idwbc,idebc,3,1,ni-2,ni_sub,ni,nj,1,ss(0,0,2))

        call s_bc8u(idwbc,idebc,ni,nj,1,ss(0,0,2))

      else if(levpbl.ge.2) then

        call s_putbufsx(idwbc,idebc,'all',3,ni-2,ni,nj,levpbl,          &
     &                  ss(0,0,2),1,3,sbuf)

        call s_putbufsx(idwbc,idebc,'all',3,ni-2,ni,nj,levpbl,          &
     &                  rr(0,0,2),2,3,sbuf)

        call s_putbufsx(idwbc,idebc,'all',3,ni-2,ni,nj,levpbl,          &
     &                  tt(0,0,2),3,3,sbuf)

        call s_shiftsx(idwbc,idebc,'all',nj,levpbl,3,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'all',1,ni_sub,ni,nj,levpbl,        &
     &                  ss(0,0,2),1,3,rbuf)

        call s_getbufsx(idwbc,idebc,'all',1,ni_sub,ni,nj,levpbl,        &
     &                  rr(0,0,2),2,3,rbuf)

        call s_getbufsx(idwbc,idebc,'all',1,ni_sub,ni,nj,levpbl,        &
     &                  tt(0,0,2),3,3,rbuf)

        call s_putbufgx(idwbc,idebc,'all',3,ni-2,ni,nj,levpbl,          &
     &                  ss(0,0,2),1,3,sbuf)

        call s_putbufgx(idwbc,idebc,'all',3,ni-2,ni,nj,levpbl,          &
     &                  rr(0,0,2),2,3,sbuf)

        call s_putbufgx(idwbc,idebc,'all',3,ni-2,ni,nj,levpbl,          &
     &                  tt(0,0,2),3,3,sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,levpbl,3,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',1,ni_sub,ni,nj,levpbl,        &
     &                  ss(0,0,2),1,3,rbuf)

        call s_getbufgx(idwbc,idebc,'all',1,ni_sub,ni,nj,levpbl,        &
     &                  rr(0,0,2),2,3,rbuf)

        call s_getbufgx(idwbc,idebc,'all',1,ni_sub,ni,nj,levpbl,        &
     &                  tt(0,0,2),3,3,rbuf)

        call s_bcyclex(idwbc,idebc,3,1,ni-2,ni_sub,ni,nj,levpbl,        &
     &                 ss(0,0,2))

        call s_bcyclex(idwbc,idebc,3,1,ni-2,ni_sub,ni,nj,levpbl,        &
     &                 rr(0,0,2))

        call s_bcyclex(idwbc,idebc,3,1,ni-2,ni_sub,ni,nj,levpbl,        &
     &                 tt(0,0,2))

        call s_bc8u(idwbc,idebc,ni,nj,levpbl,ss(0,0,2))
        call s_bc8u(idwbc,idebc,ni,nj,levpbl,rr(0,0,2))
        call s_bc8u(idwbc,idebc,ni,nj,levpbl,tt(0,0,2))

      end if

! -----

! Perform the Gauss elimination.

      call gaussel(idoneopt,1,ni_sub,1,nj-1,2,levpbl+1,ni,nj,nk,        &
     &             rr,ss,tt,u,tmp1)

! -----

! Set the bottom boundary conditions.

      call vbcu(ni,nj,nk,u)

! -----

!! -----

      end subroutine s_pblu

!-----7--------------------------------------------------------------7--

      end module m_pblu
