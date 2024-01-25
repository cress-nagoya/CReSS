!***********************************************************************
      module m_tkeflx
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/07/05
!     Modification: 1999/07/21, 1999/08/18, 1999/09/30, 1999/10/12,
!                   1999/11/01, 2000/01/17, 2000/12/19, 2001/02/24,
!                   2001/06/06, 2001/11/20, 2002/01/15, 2002/04/02,
!                   2003/01/04, 2003/03/21, 2003/04/30, 2003/05/19,
!                   2003/12/12, 2004/02/01, 2004/06/10, 2006/02/13,
!                   2006/11/06, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/08/09, 2011/11/10, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the turbulent fluxes for the turbulent kinetic energy.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc8w
      use m_comindx
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: tkeflx, s_tkeflx

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface tkeflx

        module procedure s_tkeflx

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
      subroutine s_tkeflx(fptrnopt,fpmpopt,fpmfcopt,                    &
     &                    fpdxiv,fpdyiv,fpdziv,ni,nj,nk,                &
     &                    j31,j32,jcb,rmf,tke,rkh,rkv,h1,h2,h3,         &
     &                    j31tke,j32tke,jcbtke)
!***********************************************************************

! Input variables

      integer, intent(in) :: fptrnopt
                       ! Formal parameter of unique index of trnopt

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpdxiv
                       ! Formal parameter of unique index of dxiv

      integer, intent(in) :: fpdyiv
                       ! Formal parameter of unique index of dyiv

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! turbulent kinetic energy

      real, intent(in) :: rkh(0:ni+1,0:nj+1,1:nk)
                       ! rbr x horizontal eddy viscosity / Jacobian

      real, intent(in) :: rkv(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy viscosity / Jacobian

! Output variables

      real, intent(out) :: h1(0:ni+1,0:nj+1,1:nk)
                       ! x components of turbulent fluxes

      real, intent(out) :: h2(0:ni+1,0:nj+1,1:nk)
                       ! y components of turbulent fluxes

      real, intent(out) :: h3(0:ni+1,0:nj+1,1:nk)
                       ! z components of turbulent fluxes

! Internal shared variables

      integer trnopt   ! Option for terrain height setting
      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy
      real dziv        ! Inverse of dz

      real dxv5        ! 0.5 x dxiv
      real dyv5        ! 0.5 x dyiv
      real dzv5        ! 0.5 x dziv

      real dzv125      ! 0.125 x dziv

      real, intent(inout) :: j31tke(0:ni+1,0:nj+1,1:nk)
                       ! 4.0 x j31 x tke

      real, intent(inout) :: j32tke(0:ni+1,0:nj+1,1:nk)
                       ! 4.0 x j31 x tke

      real, intent(inout) :: jcbtke(0:ni+1,0:nj+1,1:nk)
                       ! jcb x tke

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

! Remark

!     h3: This variable is also temporary.

!-----7--------------------------------------------------------------7--

! Get the namelist variables.

      call getiname(fptrnopt,trnopt)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)
      call getrname(fpdziv,dziv)

! -----

! Set the common used variables.

      dxv5=.5e0*dxiv
      dyv5=.5e0*dyiv
      dzv5=.5e0*dziv

      dzv125=.125e0*dziv

! -----

! Set the common used array.

      if(trnopt.ge.1) then

!$omp parallel default(shared) private(k)

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=2,ni-1
            j31tke(i,j,k)=j31(i,j,k)                                    &
     &        *((tke(i-1,j,k-1)+tke(i,j,k-1))+(tke(i-1,j,k)+tke(i,j,k)))
          end do
          end do

!$omp end do

        end do

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=1,ni-1
            j32tke(i,j,k)=j32(i,j,k)                                    &
     &        *((tke(i,j-1,k-1)+tke(i,j,k-1))+(tke(i,j-1,k)+tke(i,j,k)))
          end do
          end do

!$omp end do

        end do

!$omp end parallel

        call bc8w(idbbc,idtbc,ni,nj,nk,j31tke)
        call bc8w(idbbc,idtbc,ni,nj,nk,j32tke)

      end if

! -----

!! Calculate the turbulent fluxes.

!$omp parallel default(shared) private(k)

! Set the common used array.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          jcbtke(i,j,k)=tke(i,j,k)*jcb(i,j,k)
        end do
        end do

!$omp end do

      end do

! -----

! Calculate the x components of the turbulent fluxes.

      if(mfcopt.eq.1.and.mpopt.eq.5) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=1,ni-1
            h3(i,j,k)=rmf(i,j,2)*rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

        if(trnopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              h1(i,j,k)=(h3(i-1,j,k)+h3(i,j,k))                         &
     &          *(jcbtke(i,j,k)-jcbtke(i-1,j,k))*dxv5
            end do
            end do

!$omp end do

          end do

        else if(trnopt.ge.1) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              h1(i,j,k)=(h3(i-1,j,k)+h3(i,j,k))                         &
     &          *((jcbtke(i,j,k)-jcbtke(i-1,j,k))*dxv5                  &
     &          +(j31tke(i,j,k+1)-j31tke(i,j,k))*dzv125)
            end do
            end do

!$omp end do

          end do

        end if

      else

        if(trnopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              h1(i,j,k)=(rkh(i-1,j,k)+rkh(i,j,k))                       &
     &          *(jcbtke(i,j,k)-jcbtke(i-1,j,k))*dxv5
            end do
            end do

!$omp end do

          end do

        else if(trnopt.ge.1) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              h1(i,j,k)=(rkh(i-1,j,k)+rkh(i,j,k))                       &
     &          *((jcbtke(i,j,k)-jcbtke(i-1,j,k))*dxv5                  &
     &          +(j31tke(i,j,k+1)-j31tke(i,j,k))*dzv125)
            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Calculate the y components of the turbulent fluxes.

      if(mfcopt.eq.1.and.(mpopt.eq.0.or.mpopt.eq.10)) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=2,ni-2
            h3(i,j,k)=rmf(i,j,2)*rkh(i,j,k)
          end do
          end do

!$omp end do

        end do

        if(trnopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              h2(i,j,k)=(h3(i,j-1,k)+h3(i,j,k))                         &
     &          *(jcbtke(i,j,k)-jcbtke(i,j-1,k))*dyv5
            end do
            end do

!$omp end do

          end do

        else if(trnopt.ge.1) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              h2(i,j,k)=(h3(i,j-1,k)+h3(i,j,k))                         &
     &          *((jcbtke(i,j,k)-jcbtke(i,j-1,k))*dyv5                  &
     &          +(j32tke(i,j,k+1)-j32tke(i,j,k))*dzv125)
            end do
            end do

!$omp end do

          end do

        end if

      else

        if(trnopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              h2(i,j,k)=(rkh(i,j-1,k)+rkh(i,j,k))                       &
     &          *(jcbtke(i,j,k)-jcbtke(i,j-1,k))*dyv5
            end do
            end do

!$omp end do

          end do

        else if(trnopt.ge.1) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              h2(i,j,k)=(rkh(i,j-1,k)+rkh(i,j,k))                       &
     &          *((jcbtke(i,j,k)-jcbtke(i,j-1,k))*dyv5                  &
     &          +(j32tke(i,j,k+1)-j32tke(i,j,k))*dzv125)
            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Get the z components of the turbulent fluxes.

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          h3(i,j,k)                                                     &
     &      =(rkv(i,j,k-1)+rkv(i,j,k))*(tke(i,j,k)-tke(i,j,k-1))*dzv5
        end do
        end do

!$omp end do

      end do

! -----

!$omp end parallel

!! -----

      end subroutine s_tkeflx

!-----7--------------------------------------------------------------7--

      end module m_tkeflx
