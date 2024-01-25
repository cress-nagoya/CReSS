!***********************************************************************
      module m_advs
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/05,
!                   1999/08/18, 1999/09/30, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/04/18, 2000/12/18, 2001/04/15,
!                   2001/05/29, 2001/06/06, 2001/11/20, 2002/04/02,
!                   2003/01/04, 2003/03/21, 2003/04/30, 2003/05/19,
!                   2003/11/28, 2004/03/05, 2004/04/15, 2004/08/01,
!                   2006/02/13, 2006/04/03, 2006/05/12, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2010/02/01,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate optional scalar advection.

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

      public :: advs, s_advs

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface advs

        module procedure s_advs

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
      subroutine s_advs(fpadvopt,fpiwest,fpieast,fpjsouth,fpjnorth,     &
     &                  fpdxiv,fpdyiv,fpdziv,ni,nj,nk,rstxu,rstxv,      &
     &                  rstxwc,s,sfrc,vadv,tmp1,tmp2,tmp3)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpiwest
                       ! Formal parameter of unique index of iwest

      integer, intent(in) :: fpieast
                       ! Formal parameter of unique index of ieast

      integer, intent(in) :: fpjsouth
                       ! Formal parameter of unique index of jsouth

      integer, intent(in) :: fpjnorth
                       ! Formal parameter of unique index of jnorth

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

      real, intent(in) :: rstxu(0:ni+1,0:nj+1,1:nk)
                       ! u x base state density x Jacobian

      real, intent(in) :: rstxv(0:ni+1,0:nj+1,1:nk)
                       ! v x base state density x Jacobian

      real, intent(in) :: rstxwc(0:ni+1,0:nj+1,1:nk)
                       ! wc x base state density x Jacobian

      real, intent(in) :: s(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable

! Input and output variable

      real, intent(inout) :: sfrc(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar forcing term

! Internal shared variables

      integer advopt   ! Option for advection scheme

      integer iwest    ! Added index on west boundary
      integer jsouth   ! Added index on south boundary

      integer ieast    ! Subtracted index on east boundary
      integer jnorth   ! Subtracted index on north boundary

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy
      real dziv        ! Inverse of dz

      real dxv05n      ! - 0.5 x dxiv
      real dyv05n      ! - 0.5 x dyiv
      real dzv05n      ! - 0.5 x dziv

      real dxv24       ! dxiv / 24.0
      real dyv24       ! dyiv / 24.0
      real dzv24       ! dziv / 24.0

      real, intent(inout) :: vadv(0:ni+1,0:nj+1,1:nk)
                       ! Advection value in z direction

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

      call getiname(fpadvopt,advopt)
      call getiname(fpiwest,iwest)
      call getiname(fpieast,ieast)
      call getiname(fpjsouth,jsouth)
      call getiname(fpjnorth,jnorth)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)
      call getrname(fpdziv,dziv)

! -----

! Set the common used variables.

      dxv05n=-.5e0*dxiv
      dyv05n=-.5e0*dyiv
      dzv05n=-.5e0*dziv

      dxv24=oned24*dxiv
      dyv24=oned24*dyiv
      dzv24=oned24*dziv

! -----

!!! Calculate the scalar advection.

!$omp parallel default(shared) private(k)

!! Perform the centered fdm scheme.

      if(advopt.le.3) then

! Calculate the 2nd order scalar advection.

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-1
            tmp1(i,j,k)=rstxu(i,j,k)*(s(i,j,k)-s(i-1,j,k))*dxv05n
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-2
            tmp2(i,j,k)=rstxv(i,j,k)*(s(i,j,k)-s(i,j-1,k))*dyv05n
          end do
          end do

!$omp end do

        end do

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            tmp3(i,j,k)=rstxwc(i,j,k)*(s(i,j,k)-s(i,j,k-1))*dzv05n
          end do
          end do

!$omp end do

        end do

        if(advopt.eq.1) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              sfrc(i,j,k)=(tmp3(i,j,k)+tmp3(i,j,k+1))                   &
     &          +((tmp1(i,j,k)+tmp1(i+1,j,k))                           &
     &          +(tmp2(i,j,k)+tmp2(i,j+1,k)))
            end do
            end do

!$omp end do

          end do

        else

          if(advopt.eq.2) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                sfrc(i,j,k)=(tmp1(i,j,k)+tmp1(i+1,j,k))                 &
     &            +(tmp2(i,j,k)+tmp2(i,j+1,k))
              end do
              end do

!$omp end do

            end do

          else if(advopt.eq.3) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                sfrc(i,j,k)=(tmp1(i,j,k)+tmp1(i+1,j,k))                 &
     &            +(tmp2(i,j,k)+tmp2(i,j+1,k))

                vadv(i,j,k)=tmp3(i,j,k)+tmp3(i,j,k+1)

              end do
              end do

!$omp end do

            end do

          end if

        end if

! -----

! Calculate the 4th order scalar advection.

        if(advopt.eq.2.or.advopt.eq.3) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2+jsouth,nj-2-jnorth
            do i=1+iwest,ni-1-ieast
              tmp1(i,j,k)=(rstxu(i,j,k)+rstxu(i+1,j,k))                 &
     &          *(s(i+1,j,k)-s(i-1,j,k))*dxv24
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-1-jnorth
            do i=2+iwest,ni-2-ieast
              tmp2(i,j,k)=(rstxv(i,j,k)+rstxv(i,j+1,k))                 &
     &          *(s(i,j+1,k)-s(i,j-1,k))*dyv24
            end do
            end do

!$omp end do

          end do

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2+jsouth,nj-2-jnorth
            do i=2+iwest,ni-2-ieast
              sfrc(i,j,k)=fourd3*sfrc(i,j,k)                            &
     &          +((tmp1(i-1,j,k)+tmp1(i+1,j,k))                         &
     &          +(tmp2(i,j-1,k)+tmp2(i,j+1,k)))
            end do
            end do

!$omp end do

          end do

          if(advopt.eq.2) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                sfrc(i,j,k)=sfrc(i,j,k)+(tmp3(i,j,k)+tmp3(i,j,k+1))
              end do
              end do

!$omp end do

            end do

          else if(advopt.eq.3) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2+jsouth,nj-2-jnorth
              do i=2+iwest,ni-2-ieast
                tmp3(i,j,k)=(rstxwc(i,j,k)+rstxwc(i,j,k+1))             &
     &            *(s(i,j,k+1)-s(i,j,k-1))*dzv24
              end do
              end do

!$omp end do

            end do

            do k=3,nk-3

!$omp do schedule(runtime) private(i,j)

              do j=2+jsouth,nj-2-jnorth
              do i=2+iwest,ni-2-ieast
                vadv(i,j,k)                                             &
     &            =fourd3*vadv(i,j,k)+(tmp3(i,j,k-1)+tmp3(i,j,k+1))
              end do
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                sfrc(i,j,k)=sfrc(i,j,k)+vadv(i,j,k)
              end do
              end do

!$omp end do

            end do

          end if

        end if

! -----

!! -----

!! Set the forcing term with 0 in the case the Cubic Lagrange scheme is
!! performed.

      else

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            sfrc(i,j,k)=0.e0
          end do
          end do

!$omp end do

        end do

      end if

!! -----

!$omp end parallel

!!! -----

      end subroutine s_advs

!-----7--------------------------------------------------------------7--

      end module m_advs
