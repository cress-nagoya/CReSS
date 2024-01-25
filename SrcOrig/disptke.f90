!***********************************************************************
      module m_disptke
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/10/12
!     Modification: 1999/11/01, 2000/01/17, 2000/07/05, 2001/02/13,
!                   2001/05/29, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2003/11/28, 2003/12/12, 2004/02/01, 2004/03/05,
!                   2004/03/22, 2004/04/15, 2004/05/31, 2004/07/01,
!                   2004/09/01, 2004/12/17, 2006/11/06, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2011/08/09,
!                   2011/11/10, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the dissipation in the turbulent kinetic energy
!     equation.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: disptke, s_disptke

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface disptke

        module procedure s_disptke

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_disptke(fpmpopt,fpmfcopt,fpisoopt,fpdx,fpdy,fpdz,    &
     &                     ni,nj,nk,jcb,rmf,rst,priv,tke,tkefrc)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

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

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: priv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of turbulent Prandtl number

      real, intent(in) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy

! Input and output variable

      real, intent(inout) :: tkefrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in turbulent kinetic energy equation

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor
      integer isoopt   ! Option for grid shape

      real dx          ! Grid distance in x direction
      real dy          ! Grid distance in y direction
      real dz          ! Grid distance in z direction

      real dz05        ! 0.5 x dz

      real ds308       ! 0.125 x dx x dy x dz

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real ln          ! Turbulent length scale

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getiname(fpisoopt,isoopt)
      call getrname(fpdx,dx)
      call getrname(fpdy,dy)
      call getrname(fpdz,dz)

! -----

! Set the common used variables.

      dz05=.5e0*dz

      ds308=.125e0*dx*dy*dz

! -----

!! Calculate the dissipation in the turbulent kinetic energy equation.

!$omp parallel default(shared) private(k)

! Isotropic case.

      if(isoopt.eq.1) then

        if(mfcopt.eq.0) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j,ln)

            do j=2,nj-2
            do i=2,ni-2
              ln=(priv(i,j,k)-1.e0)*exp(oned3*log(ds308*jcb(i,j,k)))+eps

              tkefrc(i,j,k)=tkefrc(i,j,k)-(.37e0*priv(i,j,k)-.18e0)     &
     &          *rst(i,j,k)*tke(i,j,k)*sqrt(tke(i,j,k))/ln

            end do
            end do

!$omp end do

          end do

        else

          if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j,ln)

              do j=2,nj-2
              do i=2,ni-2
                ln=(priv(i,j,k)-1.e0)                                   &
     &            *exp(oned3*log(ds308*rmf(i,j,2)*jcb(i,j,k)))+eps

                tkefrc(i,j,k)=tkefrc(i,j,k)-(.37e0*priv(i,j,k)-.18e0)   &
     &            *rst(i,j,k)*tke(i,j,k)*sqrt(tke(i,j,k))/ln

              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j,ln)

              do j=2,nj-2
              do i=2,ni-2
                ln=(priv(i,j,k)-1.e0)                                   &
     &            *exp(oned3*log(ds308*rmf(i,j,3)*jcb(i,j,k)))+eps

                tkefrc(i,j,k)=tkefrc(i,j,k)-(.37e0*priv(i,j,k)-.18e0)   &
     &            *rst(i,j,k)*tke(i,j,k)*sqrt(tke(i,j,k))/ln

              end do
              end do

!$omp end do

            end do

          end if

        end if

! -----

! Anisotropic case.

      else if(isoopt.eq.2) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j,ln)

          do j=2,nj-2
          do i=2,ni-2
            ln=(priv(i,j,k)-1.e0)*jcb(i,j,k)*dz05+eps

            tkefrc(i,j,k)=tkefrc(i,j,k)-(.37e0*priv(i,j,k)-.18e0)       &
     &        *rst(i,j,k)*tke(i,j,k)*sqrt(tke(i,j,k))/ln

          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_disptke

!-----7--------------------------------------------------------------7--

      end module m_disptke
