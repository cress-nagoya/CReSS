!***********************************************************************
      module m_diverpih
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/12/17
!     Modification: 1999/12/20, 2000/01/17, 2000/03/09, 2000/12/18,
!                   2001/02/24, 2001/06/06, 2001/06/29, 2001/11/20,
!                   2001/12/11, 2002/04/02, 2003/01/04, 2003/03/21,
!                   2003/04/30, 2003/05/19, 2003/11/05, 2003/12/12,
!                   2004/06/10, 2004/07/01, 2006/11/06, 2007/10/19,
!                   2008/05/02, 2008/06/09, 2008/08/25, 2009/02/27,
!                   2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the divergence horizontally in the pressure equation
!     with the horizontally explicit and vertically implicit method.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_diver2d
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: diverpih, s_diverpih

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface diverpih

        module procedure s_diverpih

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
      subroutine s_diverpih(fptrnopt,fpmpopt,fpmfcopt,fpdziv,ni,nj,nk,  &
     &                      j31,j32,jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,     &
     &                      rcsq,u,v,pdiv,tmp1,tmp2,tmp3)
!***********************************************************************

! Input variables

      integer, intent(in) :: fptrnopt
                       ! Formal parameter of unique index of trnopt

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

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

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: rcsq(0:ni+1,0:nj+1,1:nk)
                       ! rbr x sound wave speed squared

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

! Output variable

      real, intent(out) :: pdiv(0:ni+1,0:nj+1,1:nk)
                       ! Horizontal divergence value
                       ! in pressure equation

! Internal shared variables

      integer trnopt   ! Option for terrain height setting
      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      integer nkm1     ! nk - 1

      real dziv        ! Inverse of dz

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fptrnopt,trnopt)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getrname(fpdziv,dziv)

! -----

! Set the common used variable.

      nkm1=nk-1

! -----

! Calculate the divergence horizontally.

      call diver2d(idmpopt,idmfcopt,iddxiv,iddyiv,ni,nj,nk,             &
     &             mf,rmf,rmf8u,rmf8v,jcb8u,jcb8v,u,v,pdiv,tmp1,tmp2)

! -----

!! Calculate the divergence horizontally in the pressure equation with
!! the horizontally explicit and vertically explicit method.

!$omp parallel default(shared) private(k)

! For the flat terrain case.

      if(trnopt.eq.0) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            pdiv(i,j,k)=rcsq(i,j,k)*pdiv(i,j,k)
          end do
          end do

!$omp end do

        end do

! -----

! For the curved grid case.

      else

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-1
            tmp1(i,j,k)=(u(i,j,k)+u(i,j,k-1))*j31(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-2
            tmp2(i,j,k)=(v(i,j,k)+v(i,j,k-1))*j32(i,j,k)
          end do
          end do

!$omp end do

        end do

        if(mfcopt.eq.0) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              tmp3(i,j,k)=.25e0*((tmp1(i,j,k)+tmp1(i+1,j,k))            &
     &          +(tmp2(i,j,k)+tmp2(i,j+1,k)))
            end do
            end do

!$omp end do

          end do

        else

          if(mpopt.eq.0.or.mpopt.eq.10) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                tmp3(i,j,k)=.25e0*(mf(i,j)*(tmp1(i,j,k)+tmp1(i+1,j,k))  &
     &            +(tmp2(i,j,k)+tmp2(i,j+1,k)))
              end do
              end do

!$omp end do

            end do

          else if(mpopt.eq.5) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                tmp3(i,j,k)=.25e0*((tmp1(i,j,k)+tmp1(i+1,j,k))          &
     &            +mf(i,j)*(tmp2(i,j,k)+tmp2(i,j+1,k)))
              end do
              end do

!$omp end do

            end do

          else

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              pdiv(i,j,nk)=.25e0*mf(i,j)
            end do
            end do

!$omp end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                tmp3(i,j,k)=pdiv(i,j,nk)*((tmp1(i,j,k)+tmp1(i+1,j,k))   &
     &            +(tmp2(i,j,k)+tmp2(i,j+1,k)))
              end do
              end do

!$omp end do

            end do

          end if

        end if

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          tmp3(i,j,nkm1)=0.e0
        end do
        end do

!$omp end do

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            pdiv(i,j,k)=rcsq(i,j,k)                                     &
     &        *(pdiv(i,j,k)+(tmp3(i,j,k)-tmp3(i,j,k+1))*dziv)
          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_diverpih

!-----7--------------------------------------------------------------7--

      end module m_diverpih
