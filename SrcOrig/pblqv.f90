!***********************************************************************
      module m_pblqv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/16
!     Modification: 2001/11/20, 2001/12/10, 2002/01/15, 2002/04/02,
!                   2002/08/15, 2002/12/27, 2003/01/20, 2003/04/30,
!                   2003/05/19, 2003/11/05, 2003/12/12, 2006/04/03,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2011/09/22, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the virtical diffusion for the water vapor mixing ratio.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_gaussel
      use m_getiname
      use m_getrname
      use m_vbcs

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: pblqv, s_pblqv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface pblqv

        module procedure s_pblqv

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
      subroutine s_pblqv(fplevpbl,fpdziv,dtb,ni,nj,nk,jcb8w,rbr,rst,    &
     &                   qvsfc,cq,khs,qv,rr,ss,tt,tmp1)
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

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: qvsfc(0:ni+1,0:nj+1)
                       ! Water vapor mixing ratio on surface

      real, intent(in) :: cq(0:ni+1,0:nj+1)
                       ! Exchange coefficient of surface moisture flux

      real, intent(in) :: khs(0:ni+1,0:nj+1,1:nk)
                       ! Eddy diffusivity in planetaty boundary layer

! Input and output variable

      real, intent(inout) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

! Internal shared variables

      integer levpbl   ! Number of planetary boundary layer

      integer levp1    ! levpbl + 1

      real dziv        ! Inverse of grid distance in z direction

      real dtdzv       ! dtb x dziv
      real dtdzv2      ! 0.5 x dtb x (dziv^2)

      real, intent(inout) :: rr(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient matrix

      real, intent(inout) :: ss(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient matrix

      real, intent(inout) :: tt(0:ni+1,0:nj+1,1:nk)
                       ! Coefficient matrix

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
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

      dtdzv=dtb*dziv
      dtdzv2=.5e0*dtb*dziv*dziv

! -----

!! Calculate the virtical diffusion for the water vapor mixing ratio.

! Set the coefficient matrix and the virtical diffusion at lowest level.

!$omp parallel default(shared) private(k)

      if(levpbl.eq.1) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          ss(i,j,2)=1.e0+cq(i,j)/rst(i,j,2)*dtdzv
        end do
        end do

!$omp end do

      else if(levpbl.ge.2) then

        do k=3,levpbl+1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            tmp1(i,j,k)=-(rbr(i,j,k-1)+rbr(i,j,k))*khs(i,j,k)           &
     &        /jcb8w(i,j,k)*dtdzv2
          end do
          end do

!$omp end do

        end do

!$omp do schedule(runtime) private(i,j,a)

        do j=1,nj-1
        do i=1,ni-1
          a=1.e0/rst(i,j,2)

          rr(i,j,2)=0.e0
          ss(i,j,2)=1.e0+a*(cq(i,j)*dtdzv-tmp1(i,j,3))
          tt(i,j,2)=a*tmp1(i,j,3)

        end do
        end do

!$omp end do

        if(levpbl.ge.3) then

          do k=3,levpbl

!$omp do schedule(runtime) private(i,j,a,b,c)

            do j=1,nj-1
            do i=1,ni-1
              a=1.e0/rst(i,j,k)

              b=a*tmp1(i,j,k)
              c=a*tmp1(i,j,k+1)

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
        do i=1,ni-1
          a=tmp1(i,j,levp1)/rst(i,j,levp1)

          rr(i,j,levp1)=a
          ss(i,j,levp1)=1.e0-a
          tt(i,j,levp1)=0.e0

        end do
        end do

!$omp end do

      end if

!$omp do schedule(runtime) private(i,j)

      do j=1,nj-1
      do i=1,ni-1
        qv(i,j,2)=qv(i,j,2)+cq(i,j)*qvsfc(i,j)/rst(i,j,2)*dtdzv
      end do
      end do

!$omp end do

!$omp end parallel

! -----

! Perform the Gauss elimination.

      call gaussel(idoneopt,1,ni-1,1,nj-1,2,levpbl+1,ni,nj,nk,rr,ss,tt, &
     &             qv,tmp1)

! -----

! Set the bottom boundary conditions.

      call vbcs(ni,nj,nk,qv)

! -----

!! -----

      end subroutine s_pblqv

!-----7--------------------------------------------------------------7--

      end module m_pblqv
