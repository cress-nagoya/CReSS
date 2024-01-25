!***********************************************************************
      module m_vbcw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/05,
!                   1999/08/18, 1999/08/23, 1999/09/06, 1999/09/14,
!                   1999/10/12, 1999/11/01, 1999/12/06, 1999/12/20,
!                   2000/01/17, 2000/03/17, 2001/04/15, 2001/06/06,
!                   2001/07/13, 2001/08/07, 2001/12/11, 2002/04/02,
!                   2002/07/23, 2002/08/15, 2003/04/30, 2003/05/19,
!                   2003/12/12, 2004/03/05, 2004/07/01, 2006/11/06,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the vertical boundary conditions for the z components of
!     velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vbcw, s_vbcw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vbcw

        module procedure s_vbcw

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
      subroutine s_vbcw(fpbbc,fptbc,fpmpopt,fpmfcopt,ni,nj,nk,          &
     &                  j31,j32,jcb8w,mf,uf,vf,wc,wf,mf25,j31u2,j32v2)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpbbc
                       ! Formal parameter of unique index of bbc

      integer, intent(in) :: fptbc
                       ! Formal parameter of unique index of tbc

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

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(in) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(in) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

! Input and output variable

      real, intent(inout) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

! Internal shared variables

      integer bbc      ! Option for bottom boundary conditions
      integer tbc      ! Option for top boundary conditions

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2
      integer nkm3     ! nk - 3

      real, intent(inout) :: mf25(0:ni+1,0:nj+1)
                       ! 0.25 x mf

      real, intent(inout) :: j31u2(0:ni+1,0:nj+1)
                       ! 2.0 x j31 x u

      real, intent(inout) :: j32v2(0:ni+1,0:nj+1)
                       ! 2.0 x j32 x v

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

! Remark

!     wf: This variable is also temporary.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpbbc,bbc)
      call getiname(fptbc,tbc)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)

! -----

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2
      nkm3=nk-3

! -----

!! Set the bottom and top boundary conditions.

!$omp parallel default(shared)

! Set the common used variable.

      if((bbc.eq.2.or.tbc.eq.2).and.(mfcopt.eq.1                        &
     &  .and.(mpopt.ne.0.and.mpopt.ne.5.and.mpopt.ne.10))) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          mf25(i,j)=.25e0*mf(i,j)
        end do
        end do

!$omp end do

      end if

! -----

! Set the bottom boundary conditions.

      if(bbc.eq.2) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni
          j31u2(i,j)=(uf(i,j,2)+uf(i,j,3))*j31(i,j,3)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1,nj
        do i=1,ni-1
          j32v2(i,j)=(vf(i,j,2)+vf(i,j,3))*j32(i,j,3)
        end do
        end do

!$omp end do

        if(mfcopt.eq.0) then

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            wf(i,j,1)=.25e0                                             &
     &        *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
          end do
          end do

!$omp end do

        else

          if(mpopt.eq.0.or.mpopt.eq.10) then

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              wf(i,j,1)=.25e0*(mf(i,j)*(j31u2(i,j)+j31u2(i+1,j))        &
     &          +(j32v2(i,j)+j32v2(i,j+1)))
            end do
            end do

!$omp end do

          else if(mpopt.eq.5) then

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              wf(i,j,1)=.25e0*((j31u2(i,j)+j31u2(i+1,j))                &
     &          +mf(i,j)*(j32v2(i,j)+j32v2(i,j+1)))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              wf(i,j,1)=mf25(i,j)                                       &
     &          *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
            end do
            end do

!$omp end do

          end if

        end if

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          wf(i,j,1)=-jcb8w(i,j,3)*wc(i,j,3)-wf(i,j,1)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni
          j31u2(i,j)=(uf(i,j,1)+uf(i,j,2))*j31(i,j,2)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1,nj
        do i=1,ni-1
          j32v2(i,j)=(vf(i,j,1)+vf(i,j,2))*j32(i,j,2)
        end do
        end do

!$omp end do

        if(mfcopt.eq.0) then

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            wf(i,j,2)=-.25e0                                            &
     &        *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
          end do
          end do

!$omp end do

        else

          if(mpopt.eq.0.or.mpopt.eq.10) then

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              wf(i,j,2)=-.25e0*(mf(i,j)*(j31u2(i,j)+j31u2(i+1,j))       &
     &          +(j32v2(i,j)+j32v2(i,j+1)))
            end do
            end do

!$omp end do

          else if(mpopt.eq.5) then

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              wf(i,j,2)=-.25e0*((j31u2(i,j)+j31u2(i+1,j))               &
     &          +mf(i,j)*(j32v2(i,j)+j32v2(i,j+1)))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              wf(i,j,2)=-mf25(i,j)                                      &
     &          *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
            end do
            end do

!$omp end do

          end if

        end if

      else if(bbc.eq.3) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          wf(i,j,1)=wf(i,j,2)
        end do
        end do

!$omp end do

      end if

! -----

! Set the top boundary conditions.

      if(tbc.eq.2) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni
          j31u2(i,j)=(uf(i,j,nkm3)+uf(i,j,nkm2))*j31(i,j,nkm2)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1,nj
        do i=1,ni-1
          j32v2(i,j)=(vf(i,j,nkm3)+vf(i,j,nkm2))*j32(i,j,nkm2)
        end do
        end do

!$omp end do

        if(mfcopt.eq.0) then

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            wf(i,j,nk)=.25e0                                            &
     &        *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
          end do
          end do

!$omp end do

        else

          if(mpopt.eq.0.or.mpopt.eq.10) then

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              wf(i,j,nk)=.25e0*(mf(i,j)*(j31u2(i,j)+j31u2(i+1,j))       &
     &          +(j32v2(i,j)+j32v2(i,j+1)))
            end do
            end do

!$omp end do

          else if(mpopt.eq.5) then

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              wf(i,j,nk)=.25e0*((j31u2(i,j)+j31u2(i+1,j))               &
     &          +mf(i,j)*(j32v2(i,j)+j32v2(i,j+1)))
            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              wf(i,j,nk)=mf25(i,j)                                      &
     &          *((j31u2(i,j)+j31u2(i+1,j))+(j32v2(i,j)+j32v2(i,j+1)))
            end do
            end do

!$omp end do

          end if

        end if

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          wf(i,j,nk)=-jcb8w(i,j,nkm2)*wc(i,j,nkm2)-wf(i,j,nk)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          wf(i,j,nkm1)=0.e0
        end do
        end do

!$omp end do

      else if(tbc.eq.3) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          wf(i,j,nk)=wf(i,j,nkm1)
        end do
        end do

!$omp end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_vbcw

!-----7--------------------------------------------------------------7--

      end module m_vbcw
