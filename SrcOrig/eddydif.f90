!***********************************************************************
      module m_eddydif
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/12/19
!     Modification: 2001/11/20, 2002/04/02, 2003/01/04, 2003/02/13,
!                   2003/03/21, 2003/04/30, 2003/05/19, 2003/11/28,
!                   2004/02/01, 2004/03/05, 2004/06/10, 2006/01/10,
!                   2006/02/13, 2006/11/06, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/01/30, 2009/02/27, 2011/08/09,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the eddy diffusivity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: eddydif, s_eddydif

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface eddydif

        module procedure s_eddydif

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
      subroutine s_eddydif(fpmfcopt,fptubopt,fpisoopt,ni,nj,nk,jcb,mf,  &
     &                     priv,rkh,rkv8w,rkv8s)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: fpisoopt
                       ! Formal parameter of unique index of isoopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: priv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of turbulent Prandtl number

! Input and output variables

      real, intent(inout) :: rkh(0:ni+1,0:nj+1,1:nk)
                       ! rbr x horizontal eddy diffusivity / Jacobian

      real, intent(inout) :: rkv8w(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy diffusivity / Jacobian
                       ! at w points

! Output variable

      real, intent(out) :: rkv8s(0:ni+1,0:nj+1,1:nk)
                       ! Half value of rbr x vertical eddy diffusivity
                       ! at scalar points

! Internal shared variables

      integer mfcopt   ! Option for map scale factor
      integer tubopt   ! Option for turbulent mixing
      integer isoopt   ! Option for grid shape

      real cpriv       ! Inverse of prnum

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real jcbiv       ! Inverse of Jacobian

! Remark

!     rkv8s: This variable is also temporary.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmfcopt,mfcopt)
      call getiname(fptubopt,tubopt)
      call getiname(fpisoopt,isoopt)

! -----

! Set the common used variables.

      cpriv=1.e0/prnum

      if(cpriv.gt.3.e0) then
        cpriv=3.e0
      end if

! -----

!!! Calculate the eddy diffusivity.

!$omp parallel default(shared) private(k)

!! Calculate the eddy diffusivity in the case the Smagorinsky
!! formulation is applied.

      if(tubopt.eq.1) then

! Isotropic case.

        if(isoopt.eq.1) then

          if(mfcopt.eq.0) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                rkv8s(i,j,k)=priv(i,j,k)*rkv8w(i,j,k)/jcb(i,j,k)
                rkh(i,j,k)=rkv8s(i,j,k)
              end do
              end do

!$omp end do

            end do

          else

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1
                rkv8s(i,j,k)=priv(i,j,k)*rkv8w(i,j,k)/jcb(i,j,k)
                rkh(i,j,k)=mf(i,j)*rkv8s(i,j,k)
              end do
              end do

!$omp end do

            end do

          end if

! -----

! Anisotropic case.

        else if(isoopt.eq.2) then

          if(mfcopt.eq.0) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,jcbiv)

              do j=1,nj-1
              do i=1,ni-1
                jcbiv=1.e0/jcb(i,j,k)

                rkh(i,j,k)=jcbiv*cpriv*rkh(i,j,k)
                rkv8s(i,j,k)=jcbiv*priv(i,j,k)*rkv8w(i,j,k)

              end do
              end do

!$omp end do

            end do

          else

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,jcbiv)

              do j=1,nj-1
              do i=1,ni-1
                jcbiv=1.e0/jcb(i,j,k)

                rkh(i,j,k)=jcbiv*cpriv*mf(i,j)*rkh(i,j,k)
                rkv8s(i,j,k)=jcbiv*priv(i,j,k)*rkv8w(i,j,k)

              end do
              end do

!$omp end do

            end do

          end if

        end if

! -----

!! -----

!! Calculate the eddy diffusivity in the case the Deardorff formulation
!! is applied.

      else if(tubopt.ge.2) then

! Isotropic case.

        if(isoopt.eq.1) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              rkv8s(i,j,k)=priv(i,j,k)*rkv8w(i,j,k)
              rkh(i,j,k)=rkv8s(i,j,k)
            end do
            end do

!$omp end do

          end do

! -----

! Anisotropic case.

        else if(isoopt.eq.2) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              rkh(i,j,k)=cpriv*rkh(i,j,k)
              rkv8s(i,j,k)=priv(i,j,k)*rkv8w(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

! Get the virtical eddy diffusivity at the w and scalar points.

      if(tubopt.eq.1) then

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            rkv8w(i,j,k)=.5e0*(rkv8s(i,j,k-1)+rkv8s(i,j,k))
          end do
          end do

!$omp end do

        end do

      else if(tubopt.ge.2) then

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            rkv8w(i,j,k)=.5e0*(rkv8s(i,j,k-1)+rkv8s(i,j,k))
          end do
          end do

!$omp end do

        end do

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            rkv8s(i,j,k)=.5e0*jcb(i,j,k)*rkv8s(i,j,k)
          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!!! -----

      end subroutine s_eddydif

!-----7--------------------------------------------------------------7--

      end module m_eddydif
