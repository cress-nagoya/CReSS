!***********************************************************************
      module m_smoo4uvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/09/16
!     Modification: 1999/10/12, 1999/11/01, 1999/11/19, 1999/11/24,
!                   2000/01/17, 2000/04/18, 2000/12/18, 2001/06/06,
!                   2001/06/29, 2002/04/02, 2002/07/23, 2002/08/15,
!                   2003/01/04, 2003/04/30, 2003/05/19, 2003/07/15,
!                   2003/11/05, 2003/11/28, 2003/12/12, 2004/01/23,
!                   2004/03/05, 2004/04/15, 2004/08/20, 2004/09/10,
!                   2005/04/04, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the 4th order smoothing for the velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: smoo4uvw, s_smoo4uvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface smoo4uvw

        module procedure s_smoo4uvw

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_smoo4uvw(fpsmtopt,fpiwest,fpieast,fpjsouth,fpjnorth, &
     &                   fpsmhcoe,fpsmvcoe,ni,nj,nk,jcb8u,jcb8v,jcb8w,  &
     &                   ubr,vbr,rst8u,rst8v,rst8w,u,v,w,ufrc,vfrc,wfrc,&
     &                   tmp1,tmp2,tmp3,tmp4,tmp5)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: fpiwest
                       ! Formal parameter of unique index of iwest

      integer, intent(in) :: fpieast
                       ! Formal parameter of unique index of ieast

      integer, intent(in) :: fpjsouth
                       ! Formal parameter of unique index of jsouth

      integer, intent(in) :: fpjnorth
                       ! Formal parameter of unique index of jnorth

      integer, intent(in) :: fpsmhcoe
                       ! Formal parameter of unique index of smhcoe

      integer, intent(in) :: fpsmvcoe
                       ! Formal parameter of unique index of smvcoe

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(in) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(in) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at u points

      real, intent(in) :: rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at v points

      real, intent(in) :: rst8w(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian at w points

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

      integer smtopt   ! Option for numerical smoothing

      integer iwest    ! Added index on west boundary
      integer jsouth   ! Added index on south boundary

      integer ieast    ! Subtracted index on east boundary
      integer jnorth   ! Subtracted index on north boundary

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

      real smhcoe      ! Horizontal smoothig coefficient
      real smvcoe      ! Vertical smoothig coefficient

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

      call getiname(fpsmtopt,smtopt)
      call getiname(fpiwest,iwest)
      call getiname(fpieast,ieast)
      call getiname(fpjsouth,jsouth)
      call getiname(fpjnorth,jnorth)
      call getrname(fpsmhcoe,smhcoe)
      call getrname(fpsmvcoe,smvcoe)

! -----

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2

! -----

!! Calculate the 4th order velocity numerical smoothing.

!$omp parallel default(shared) private(k)

! Calculate the 4th order u smoothing.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=jsouth,nj-jnorth
        do i=iwest,ni+1-ieast
          tmp4(i,j,k)=rst8u(i,j,k)*(u(i,j,k)-ubr(i,j,k))/jcb8u(i,j,k)
          tmp5(i,j,k)=2.e0*tmp4(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=1+iwest,ni-ieast
          tmp1(i,j,k)=(tmp4(i+1,j,k)+tmp4(i-1,j,k))-tmp5(i,j,k)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1+jsouth,nj-1-jnorth
        do i=2,ni-1
          tmp2(i,j,k)=(tmp4(i,j+1,k)+tmp4(i,j-1,k))-tmp5(i,j,k)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-1
          tmp3(i,j,k)=(tmp4(i,j,k+1)+tmp4(i,j,k-1))-tmp5(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-1
          ufrc(i,j,k)=ufrc(i,j,k)                                       &
     &      +smhcoe*(tmp1(i,j,k)+tmp2(i,j,k))+smvcoe*tmp3(i,j,k)
        end do
        end do

!$omp end do

      end do

      if(mod(smtopt,10).eq.2) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2+jsouth,nj-2-jnorth
          do i=2+iwest,ni-1-ieast
            ufrc(i,j,k)=ufrc(i,j,k)                                     &
     &        +smhcoe*((tmp1(i,j,k)-(tmp1(i+1,j,k)+tmp1(i-1,j,k)))      &
     &        +(tmp2(i,j,k)-(tmp2(i,j+1,k)+tmp2(i,j-1,k))))
          end do
          end do

!$omp end do

        end do

      else

!$omp do schedule(runtime) private(i,j)

        do j=2+jsouth,nj-2-jnorth
        do i=2+iwest,ni-1-ieast
          tmp3(i,j,1)=tmp3(i,j,2)
          tmp3(i,j,nkm1)=tmp3(i,j,nkm2)
        end do
        end do

!$omp end do

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2+jsouth,nj-2-jnorth
          do i=2+iwest,ni-1-ieast
            ufrc(i,j,k)=ufrc(i,j,k)                                     &
     &        +(smvcoe*(tmp3(i,j,k)-(tmp3(i,j,k+1)+tmp3(i,j,k-1)))      &
     &        +smhcoe*((tmp1(i,j,k)-(tmp1(i+1,j,k)+tmp1(i-1,j,k)))      &
     &        +(tmp2(i,j,k)-(tmp2(i,j+1,k)+tmp2(i,j-1,k)))))
          end do
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the 4th order v smoothing.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=jsouth,nj+1-jnorth
        do i=iwest,ni-ieast
          tmp4(i,j,k)=rst8v(i,j,k)*(v(i,j,k)-vbr(i,j,k))/jcb8v(i,j,k)
          tmp5(i,j,k)=2.e0*tmp4(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=1+iwest,ni-1-ieast
          tmp1(i,j,k)=(tmp4(i+1,j,k)+tmp4(i-1,j,k))-tmp5(i,j,k)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1+jsouth,nj-jnorth
        do i=2,ni-2
          tmp2(i,j,k)=(tmp4(i,j+1,k)+tmp4(i,j-1,k))-tmp5(i,j,k)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-2
          tmp3(i,j,k)=(tmp4(i,j,k+1)+tmp4(i,j,k-1))-tmp5(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-1
        do i=2,ni-2
          vfrc(i,j,k)=vfrc(i,j,k)                                       &
     &      +smhcoe*(tmp1(i,j,k)+tmp2(i,j,k))+smvcoe*tmp3(i,j,k)
        end do
        end do

!$omp end do

      end do

      if(mod(smtopt,10).eq.2) then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2+jsouth,nj-1-jnorth
          do i=2+iwest,ni-2-ieast
            vfrc(i,j,k)=vfrc(i,j,k)                                     &
     &        +smhcoe*((tmp1(i,j,k)-(tmp1(i+1,j,k)+tmp1(i-1,j,k)))      &
     &        +(tmp2(i,j,k)-(tmp2(i,j+1,k)+tmp2(i,j-1,k))))
          end do
          end do

!$omp end do

        end do

      else

!$omp do schedule(runtime) private(i,j)

        do j=2+jsouth,nj-1-jnorth
        do i=2+iwest,ni-2-ieast
          tmp3(i,j,1)=tmp3(i,j,2)
          tmp3(i,j,nkm1)=tmp3(i,j,nkm2)
        end do
        end do

!$omp end do

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2+jsouth,nj-1-jnorth
          do i=2+iwest,ni-2-ieast
            vfrc(i,j,k)=vfrc(i,j,k)                                     &
     &        +(smvcoe*(tmp3(i,j,k)-(tmp3(i,j,k+1)+tmp3(i,j,k-1)))      &
     &        +smhcoe*((tmp1(i,j,k)-(tmp1(i+1,j,k)+tmp1(i-1,j,k)))      &
     &        +(tmp2(i,j,k)-(tmp2(i,j+1,k)+tmp2(i,j-1,k)))))
          end do
          end do

!$omp end do

        end do

      end if

! -----

! Calculate the 4th order w smoothing.

      do k=1,nk

!$omp do schedule(runtime) private(i,j)

        do j=jsouth,nj-jnorth
        do i=iwest,ni-ieast
          tmp4(i,j,k)=rst8w(i,j,k)*w(i,j,k)/jcb8w(i,j,k)
          tmp5(i,j,k)=2.e0*tmp4(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=1+iwest,ni-1-ieast
          tmp1(i,j,k)=(tmp4(i+1,j,k)+tmp4(i-1,j,k))-tmp5(i,j,k)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1+jsouth,nj-1-jnorth
        do i=2,ni-2
          tmp2(i,j,k)=(tmp4(i,j+1,k)+tmp4(i,j-1,k))-tmp5(i,j,k)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          tmp3(i,j,k)=(tmp4(i,j,k+1)+tmp4(i,j,k-1))-tmp5(i,j,k)
        end do
        end do

!$omp end do

      end do

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=2,nj-2
        do i=2,ni-2
          wfrc(i,j,k)=wfrc(i,j,k)                                       &
     &      +smhcoe*(tmp1(i,j,k)+tmp2(i,j,k))+smvcoe*tmp3(i,j,k)
        end do
        end do

!$omp end do

      end do

      if(mod(smtopt,10).eq.2) then

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2+jsouth,nj-2-jnorth
          do i=2+iwest,ni-2-ieast
            wfrc(i,j,k)=wfrc(i,j,k)                                     &
     &        +smhcoe*((tmp1(i,j,k)-(tmp1(i+1,j,k)+tmp1(i-1,j,k)))      &
     &        +(tmp2(i,j,k)-(tmp2(i,j+1,k)+tmp2(i,j-1,k))))
          end do
          end do

!$omp end do

        end do

      else

!$omp do schedule(runtime) private(i,j)

        do j=2+jsouth,nj-2-jnorth
        do i=2+iwest,ni-2-ieast
          tmp3(i,j,1)=tmp3(i,j,2)
          tmp3(i,j,nk)=tmp3(i,j,nkm1)
        end do
        end do

!$omp end do

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2+jsouth,nj-2-jnorth
          do i=2+iwest,ni-2-ieast
            wfrc(i,j,k)=wfrc(i,j,k)                                     &
     &        +(smvcoe*(tmp3(i,j,k)-(tmp3(i,j,k+1)+tmp3(i,j,k-1)))      &
     &        +smhcoe*((tmp1(i,j,k)-(tmp1(i+1,j,k)+tmp1(i-1,j,k)))      &
     &        +(tmp2(i,j,k)-(tmp2(i,j+1,k)+tmp2(i,j-1,k)))))
          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_smoo4uvw

!-----7--------------------------------------------------------------7--

      end module m_smoo4uvw
