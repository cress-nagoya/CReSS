!***********************************************************************
      module m_curveuvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/11/11
!     Modification: 2003/01/04, 2003/04/30, 2003/05/19, 2003/10/31,
!                   2003/11/28, 2003/12/12, 2004/04/15, 2004/06/10,
!                   2006/02/13, 2006/11/06, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2011/08/09, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the curvature of earth.

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

      public :: curveuvw, s_curveuvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface curveuvw

        module procedure s_curveuvw

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
      subroutine s_curveuvw(fpmpopt,fpmfcopt,ni,nj,nk,rmf8u,rmf8v,      &
     &                      rst,u,v,w,ufrc,vfrc,wfrc,tmp1,tmp2,         &
     &                      tmp3,tmp4,tmp5)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

! Input and output variables

      real, intent(inout) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(inout) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(inout) :: wfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in w equation

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      real rev125      ! 0.125 / rearth

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp5(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)

! -----

! Set the common used variable.

      rev125=.125e0/rearth

! -----

!! Calculate the curvature of earth in the x and y components of
!! velocity equation.

!$omp parallel default(shared) private(k)

! Set the common used array.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          tmp1(i,j,k)=u(i,j,k)+u(i+1,j,k)
          tmp2(i,j,k)=v(i,j,k)+v(i,j+1,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          tmp3(i,j,k)=rst(i,j,k)*(w(i,j,k)+w(i,j,k+1))
        end do
        end do

!$omp end do

      end do

! -----

! For x and y components of velocity.

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=1,ni-1
          tmp4(i,j,k)=tmp1(i,j,k)*tmp3(i,j,k)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=2,ni-2
          tmp5(i,j,k)=tmp2(i,j,k)*tmp3(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-1
          ufrc(i,j,k)=ufrc(i,j,k)-rev125*(tmp4(i-1,j,k)+tmp4(i,j,k))
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-2
          vfrc(i,j,k)=vfrc(i,j,k)-rev125*(tmp5(i,j-1,k)+tmp5(i,j,k))
        end do
        end do

!$omp end do

      end do

! -----

! Calculate the curvature of earth in the z components of velocity
! equation.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          tmp3(i,j,k)=rst(i,j,k)                                        &
     &      *(tmp1(i,j,k)*tmp1(i,j,k)+tmp2(i,j,k)*tmp2(i,j,k))
        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          wfrc(i,j,k)=wfrc(i,j,k)+rev125*(tmp3(i,j,k-1)+tmp3(i,j,k))
        end do
        end do

!$omp end do

      end do

! -----

! Add the terms in the case of turning on the option for map scale
! factor.

      if(mfcopt.eq.1) then

        if(mpopt.eq.0.or.mpopt.eq.10) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              tmp3(i,j,k)=rmf8v(i,j,3)*rst(i,j,k)*tmp1(i,j,k)
            end do
            end do

!$omp end do

          end do

        else if(mpopt.eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              tmp3(i,j,k)=rmf8u(i,j,3)*rst(i,j,k)*tmp2(i,j,k)
            end do
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              tmp3(i,j,k)=rst(i,j,k)                                    &
     &          *(rmf8v(i,j,3)*tmp1(i,j,k)-rmf8u(i,j,3)*tmp2(i,j,k))
            end do
            end do

!$omp end do

          end do

        end if

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=2,ni-2
            tmp1(i,j,k)=tmp1(i,j,k)*tmp3(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=1,ni-1
            tmp2(i,j,k)=tmp2(i,j,k)*tmp3(i,j,k)
          end do
          end do

!$omp end do

        end do

        if(mpopt.eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              ufrc(i,j,k)=ufrc(i,j,k)-(tmp2(i-1,j,k)+tmp2(i,j,k))
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              vfrc(i,j,k)=vfrc(i,j,k)+(tmp1(i,j-1,k)+tmp1(i,j,k))
            end do
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              ufrc(i,j,k)=ufrc(i,j,k)+(tmp2(i-1,j,k)+tmp2(i,j,k))
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              vfrc(i,j,k)=vfrc(i,j,k)-(tmp1(i,j-1,k)+tmp1(i,j,k))
            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_curveuvw

!-----7--------------------------------------------------------------7--

      end module m_curveuvw
