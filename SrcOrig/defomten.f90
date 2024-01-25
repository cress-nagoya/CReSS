!***********************************************************************
      module m_defomten
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/07/05
!     Modification: 1999/07/19, 1999/07/21, 1999/07/28, 1999/08/03,
!                   1999/08/18, 1999/08/23, 1999/09/16, 1999/09/30,
!                   1999/10/07, 1999/10/12, 1999/10/21, 1999/11/01,
!                   1999/11/19, 1999/11/30, 2000/01/17, 2000/04/18,
!                   2000/07/05, 2000/12/19, 2001/02/24, 2001/04/15,
!                   2001/05/29, 2001/06/06, 2001/06/29, 2001/08/07,
!                   2002/04/02, 2002/06/06, 2002/07/23, 2003/01/04,
!                   2003/03/21, 2003/04/30, 2003/05/19, 2003/12/12,
!                   2004/02/01, 2004/06/10, 2004/09/10, 2006/11/06,
!                   2007/10/19, 2008/05/02, 2008/06/09, 2008/08/25,
!                   2009/02/27, 2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the deformation tensor.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc8w
      use m_bc8u
      use m_bc8v
      use m_bcten
      use m_comindx
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: defomten, s_defomten

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface defomten

        module procedure s_defomten

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
      subroutine s_defomten(fptrnopt,fpmpopt,fpmfcopt,                  &
     &                      fpiwest,fpieast,fpjsouth,fpjnorth,          &
     &                      fpdxiv,fpdyiv,fpdziv,ni,nj,nk,j31,j32,      &
     &                      jcb,jcb8u,jcb8v,jcb8w,mf,mf8u,mf8v,u,v,w,   &
     &                      s11,s22,s33,s12,s13,s23,s31,s32,            &
     &                      tmp1,tmp2,tmp3,tmp4)
!***********************************************************************

! Input variables

      integer, intent(in) :: fptrnopt
                       ! Formal parameter of unique index of trnopt

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

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

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: jcb8w(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at w points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(in) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

! Output variables

      real, intent(out) :: s11(0:ni+1,0:nj+1,1:nk)
                       ! x-x components of deformation tensor

      real, intent(out) :: s22(0:ni+1,0:nj+1,1:nk)
                       ! y-y components of deformation tensor

      real, intent(out) :: s33(0:ni+1,0:nj+1,1:nk)
                       ! z-z components of deformation tensor

      real, intent(out) :: s12(0:ni+1,0:nj+1,1:nk)
                       ! x-y components of deformation tensor

      real, intent(out) :: s13(0:ni+1,0:nj+1,1:nk)
                       ! x-z components of deformation tensor

      real, intent(out) :: s23(0:ni+1,0:nj+1,1:nk)
                       ! y-z components of deformation tensor

      real, intent(out) :: s31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of deformation tensor

      real, intent(out) :: s32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of deformation tensor

! Internal shared variables

      integer trnopt   ! Option for terrain height setting
      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      integer iwest    ! Added index on west boundary
      integer jsouth   ! Added index on south boundary

      integer ieast    ! Subtracted index on east boundary
      integer jnorth   ! Subtracted index on north boundary

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy
      real dziv        ! Inverse of dz

      real dxv25       ! 0.25 x dxiv
      real dyv25       ! 0.25 x dyiv

      real dzv25       ! 0.25 x dziv
      real dzv125      ! 0.125 x dziv

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real jcbiv2      ! 2.0 / jcb
      real mfdvj2      ! 2.0 x mf / jcb

! Remark

!     s11,s22,s33,s13,s23: These variables are also temporary.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fptrnopt,trnopt)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getiname(fpiwest,iwest)
      call getiname(fpieast,ieast)
      call getiname(fpjsouth,jsouth)
      call getiname(fpjnorth,jnorth)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)
      call getrname(fpdziv,dziv)

! -----

! Set the common used variables.

      dxv25=.25e0*dxiv
      dyv25=.25e0*dyiv

      dzv25=.25e0*dziv
      dzv125=.125e0*dziv

! -----

! Set the common used array.

      if(trnopt.ge.1) then

!$omp parallel default(shared) private(k)

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=jsouth,nj-jnorth
          do i=1,ni
            s11(i,j,k)=u(i,j,k-1)+u(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=1,nj
          do i=iwest,ni-ieast
            s22(i,j,k)=v(i,j,k-1)+v(i,j,k)
          end do
          end do

!$omp end do

        end do

!$omp end parallel

        call bc8w(idbbc,idtbc,ni,nj,nk,s11)
        call bc8w(idbbc,idtbc,ni,nj,nk,s22)

      end if

! -----

!!! Calculate the deformation tensor.

!$omp parallel default(shared) private(k)

! Calculate the diagonal and the x-y components of the deformation
! tensor.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=jsouth,nj-jnorth
        do i=iwest,ni+1-ieast
          tmp1(i,j,k)=u(i,j,k)*jcb8u(i,j,k)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=jsouth,nj+1-jnorth
        do i=iwest,ni-ieast
          tmp2(i,j,k)=v(i,j,k)*jcb8v(i,j,k)
        end do
        end do

!$omp end do

      end do

      if(mfcopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1+jsouth,nj-jnorth
          do i=1+iwest,ni-ieast
            tmp3(i,j,k)=4.e0/((jcb(i-1,j-1,k)+jcb(i,j,k))               &
     &        +(jcb(i-1,j,k)+jcb(i,j-1,k)))
          end do
          end do

!$omp end do

        end do

      else

        if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

!$omp do schedule(runtime) private(i,j)

          do j=1+jsouth,nj-jnorth
          do i=1+iwest,ni-ieast
            tmp4(i,j,1)=(mf(i-1,j-1)+mf(i,j))+(mf(i-1,j)+mf(i,j-1))
          end do
          end do

!$omp end do

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1+iwest,ni-ieast
              tmp3(i,j,k)=4.e0/((jcb(i-1,j-1,k)+jcb(i,j,k))             &
     &          +(jcb(i-1,j,k)+jcb(i,j-1,k)))
            end do
            end do

!$omp end do

          end do

        else

!$omp do schedule(runtime) private(i,j)

          do j=1+jsouth,nj-jnorth
          do i=1+iwest,ni-ieast
            tmp4(i,j,1)=(mf(i-1,j-1)+mf(i,j))+(mf(i-1,j)+mf(i,j-1))
          end do
          end do

!$omp end do

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1+iwest,ni-ieast
              tmp3(i,j,k)=tmp4(i,j,1)/((jcb(i-1,j-1,k)+jcb(i,j,k))      &
     &          +(jcb(i-1,j,k)+jcb(i,j-1,k)))
            end do
            end do

!$omp end do

          end do

        end if

      end if

      if(trnopt.eq.0) then

        if(mfcopt.eq.1.and.(mpopt.eq.0.or.mpopt.eq.10)) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1+iwest,ni-ieast
              s12(i,j,k)=tmp3(i,j,k)*((tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv  &
     &          +tmp4(i,j,1)*(tmp2(i,j,k)-tmp2(i-1,j,k))*dxv25)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              s11(i,j,k)=(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv
              s22(i,j,k)=(tmp2(i,j+1,k)-tmp2(i,j,k))*dyiv
            end do
            end do

!$omp end do

          end do

        else if(mfcopt.eq.1.and.mpopt.eq.5) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1+iwest,ni-ieast
              s12(i,j,k)=tmp3(i,j,k)*((tmp2(i,j,k)-tmp2(i-1,j,k))*dxiv  &
     &          +tmp4(i,j,1)*(tmp1(i,j,k)-tmp1(i,j-1,k))*dyv25)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              s11(i,j,k)=(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv
              s22(i,j,k)=(tmp2(i,j+1,k)-tmp2(i,j,k))*dyiv
            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1+iwest,ni-ieast
              s12(i,j,k)=tmp3(i,j,k)*((tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv  &
     &          +(tmp2(i,j,k)-tmp2(i-1,j,k))*dxiv)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              s11(i,j,k)=(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv
              s22(i,j,k)=(tmp2(i,j+1,k)-tmp2(i,j,k))*dyiv
            end do
            end do

!$omp end do

          end do

        end if

      else if(trnopt.ge.1) then

        if(mfcopt.eq.1.and.(mpopt.eq.0.or.mpopt.eq.10)) then

          do k=1,nk

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1+iwest,ni-ieast
              s33(i,j,k)=.25e0*tmp4(i,j,1)                              &
     &          *(s22(i-1,j,k)+s22(i,j,k))*(j31(i,j-1,k)+j31(i,j,k))    &
     &          +(s11(i,j-1,k)+s11(i,j,k))*(j32(i-1,j,k)+j32(i,j,k))
            end do
            end do

!$omp end do

          end do

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1+iwest,ni-ieast
              s12(i,j,k)=tmp3(i,j,k)*((s33(i,j,k+1)-s33(i,j,k))*dzv125  &
     &          +(tmp4(i,j,1)*(tmp2(i,j,k)-tmp2(i-1,j,k))*dxv25         &
     &          +(tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv))
            end do
            end do

!$omp end do

          end do

        else if(mfcopt.eq.1.and.mpopt.eq.5) then

          do k=1,nk

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1+iwest,ni-ieast
              s33(i,j,k)=.25e0*tmp4(i,j,1)                              &
     &          *(s11(i,j-1,k)+s11(i,j,k))*(j32(i-1,j,k)+j32(i,j,k))    &
     &          +(s22(i-1,j,k)+s22(i,j,k))*(j31(i,j-1,k)+j31(i,j,k))
            end do
            end do

!$omp end do

          end do

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1+iwest,ni-ieast
              s12(i,j,k)=tmp3(i,j,k)*((s33(i,j,k+1)-s33(i,j,k))*dzv125  &
     &          +(tmp4(i,j,1)*(tmp1(i,j,k)-tmp1(i,j-1,k))*dyv25         &
     &          +(tmp2(i,j,k)-tmp2(i-1,j,k))*dxiv))
            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1+iwest,ni-ieast
              s33(i,j,k)                                                &
     &          =(s22(i-1,j,k)+s22(i,j,k))*(j31(i,j-1,k)+j31(i,j,k))    &
     &          +(s11(i,j-1,k)+s11(i,j,k))*(j32(i-1,j,k)+j32(i,j,k))
            end do
            end do

!$omp end do

          end do

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1+iwest,ni-ieast
              s12(i,j,k)=tmp3(i,j,k)*((s33(i,j,k+1)-s33(i,j,k))*dzv125  &
     &          +((tmp2(i,j,k)-tmp2(i-1,j,k))*dxiv                      &
     &          +(tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv))
            end do
            end do

!$omp end do

          end do

        end if

        do k=1,nk

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni
            s11(i,j,k)=s11(i,j,k)*j31(i,j,k)
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=1,nj
          do i=1,ni-1
            s22(i,j,k)=s22(i,j,k)*j32(i,j,k)
          end do
          end do

!$omp end do

        end do

        do k=1,nk

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            s13(i,j,k)=s11(i,j,k)+s11(i+1,j,k)
            s23(i,j,k)=s22(i,j,k)+s22(i,j+1,k)
          end do
          end do

!$omp end do

        end do

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1
            s11(i,j,k)=(tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv                 &
     &        +(s13(i,j,k+1)-s13(i,j,k))*dzv25

            s22(i,j,k)=(tmp2(i,j+1,k)-tmp2(i,j,k))*dyiv                 &
     &        +(s23(i,j,k+1)-s23(i,j,k))*dzv25

          end do
          end do

!$omp end do

        end do

      end if

      if(mfcopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,jcbiv2)

          do j=1,nj-1
          do i=1,ni-1
            jcbiv2=2.e0/jcb(i,j,k)

            s11(i,j,k)=jcbiv2*s11(i,j,k)
            s22(i,j,k)=jcbiv2*s22(i,j,k)
            s33(i,j,k)=jcbiv2*(w(i,j,k+1)-w(i,j,k))*dziv

          end do
          end do

!$omp end do

        end do

      else

        if(mpopt.eq.0.or.mpopt.eq.10) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,jcbiv2)

            do j=1,nj-1
            do i=1,ni-1
              jcbiv2=2.e0/jcb(i,j,k)

              s11(i,j,k)=jcbiv2*mf(i,j)*s11(i,j,k)
              s22(i,j,k)=jcbiv2*s22(i,j,k)
              s33(i,j,k)=jcbiv2*(w(i,j,k+1)-w(i,j,k))*dziv

            end do
            end do

!$omp end do

          end do

        else if(mpopt.eq.5) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,jcbiv2)

            do j=1,nj-1
            do i=1,ni-1
              jcbiv2=2.e0/jcb(i,j,k)

              s11(i,j,k)=jcbiv2*s11(i,j,k)
              s22(i,j,k)=jcbiv2*mf(i,j)*s22(i,j,k)
              s33(i,j,k)=jcbiv2*(w(i,j,k+1)-w(i,j,k))*dziv

            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,jcbiv2,mfdvj2)

            do j=1,nj-1
            do i=1,ni-1
              jcbiv2=2.e0/jcb(i,j,k)
              mfdvj2=jcbiv2*mf(i,j)

              s11(i,j,k)=mfdvj2*s11(i,j,k)
              s22(i,j,k)=mfdvj2*s22(i,j,k)
              s33(i,j,k)=jcbiv2*(w(i,j,k+1)-w(i,j,k))*dziv

            end do
            end do

!$omp end do

          end do

        end if

      end if

! -----

!! Calculate the x-z, y-z, z-x and z-y components of the deformation
!! tensor.

! Set the common used array.

      do k=1,nk

!$omp do schedule(runtime) private(i,j)

        do j=jsouth,nj-jnorth
        do i=iwest,ni-ieast
          tmp1(i,j,k)=w(i,j,k)*jcb8w(i,j,k)
        end do
        end do

!$omp end do

      end do

! -----

! Calculate the x-z and z-x components of the deformation tensor.

      do k=1,nk

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1+iwest,ni-ieast
          tmp3(i,j,k)=2.e0/(jcb8w(i-1,j,k)+jcb8w(i,j,k))
        end do
        end do

!$omp end do

      end do

      if(trnopt.eq.0) then

        if(mfcopt.eq.1.and.mpopt.ne.5) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1+iwest,ni-ieast
              s13(i,j,k)=(mf8u(i,j)*(tmp1(i,j,k)-tmp1(i-1,j,k))*dxiv    &
     &          +(u(i,j,k)-u(i,j,k-1))*dziv)*tmp3(i,j,k)
            end do
            end do

!$omp end do

          end do

        else

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1+iwest,ni-ieast
              s13(i,j,k)=((tmp1(i,j,k)-tmp1(i-1,j,k))*dxiv              &
     &          +(u(i,j,k)-u(i,j,k-1))*dziv)*tmp3(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

      else if(trnopt.ge.1) then

        do k=1,nk

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1+iwest,ni-ieast
            tmp2(i,j,k)=(w(i-1,j,k)+w(i,j,k))*j31(i,j,k)
          end do
          end do

!$omp end do

        end do

        if(mfcopt.eq.1.and.mpopt.ne.5) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1+iwest,ni-ieast
              tmp4(i,j,k)=tmp2(i,j,k)+tmp2(i,j,k+1)
            end do
            end do

!$omp end do

          end do

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1+iwest,ni-ieast
              s13(i,j,k)=((u(i,j,k)-u(i,j,k-1))*dziv                    &
     &          +mf8u(i,j)*((tmp1(i,j,k)-tmp1(i-1,j,k))*dxiv            &
     &          +(tmp4(i,j,k)-tmp4(i,j,k-1))*dzv25))*tmp3(i,j,k)
            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1+iwest,ni-ieast
              tmp4(i,j,k)=u(i,j,k)+.25e0*(tmp2(i,j,k)+tmp2(i,j,k+1))
            end do
            end do

!$omp end do

          end do

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1+iwest,ni-ieast
              s13(i,j,k)=((tmp1(i,j,k)-tmp1(i-1,j,k))*dxiv              &
     &          +(tmp4(i,j,k)-tmp4(i,j,k-1))*dziv)*tmp3(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

      end if

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1+iwest,ni-ieast
          s31(i,j,k)=s13(i,j,k)
        end do
        end do

!$omp end do

      end do

! -----

! Calculate the y-z and z-y components of the deformation tensor.

      do k=1,nk

!$omp do schedule(runtime) private(i,j)

        do j=1+jsouth,nj-jnorth
        do i=1,ni-1
          tmp3(i,j,k)=2.e0/(jcb8w(i,j-1,k)+jcb8w(i,j,k))
        end do
        end do

!$omp end do

      end do

      if(trnopt.eq.0) then

        if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1,ni-1
              s23(i,j,k)=(mf8v(i,j)*(tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv    &
     &          +(v(i,j,k)-v(i,j,k-1))*dziv)*tmp3(i,j,k)
            end do
            end do

!$omp end do

          end do

        else

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1,ni-1
              s23(i,j,k)=((tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv              &
     &          +(v(i,j,k)-v(i,j,k-1))*dziv)*tmp3(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

      else if(trnopt.ge.1) then

        do k=1,nk

!$omp do schedule(runtime) private(i,j)

          do j=1+jsouth,nj-jnorth
          do i=1,ni-1
            tmp2(i,j,k)=(w(i,j-1,k)+w(i,j,k))*j32(i,j,k)
          end do
          end do

!$omp end do

        end do

        if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1,ni-1
              tmp4(i,j,k)=tmp2(i,j,k)+tmp2(i,j,k+1)
            end do
            end do

!$omp end do

          end do

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1,ni-1
              s23(i,j,k)=((v(i,j,k)-v(i,j,k-1))*dziv                    &
     &          +mf8v(i,j)*((tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv            &
     &          +(tmp4(i,j,k)-tmp4(i,j,k-1))*dzv25))*tmp3(i,j,k)
            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1,ni-1
              tmp4(i,j,k)=v(i,j,k)+.25e0*(tmp2(i,j,k)+tmp2(i,j,k+1))
            end do
            end do

!$omp end do

          end do

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

            do j=1+jsouth,nj-jnorth
            do i=1,ni-1
              s23(i,j,k)=((tmp1(i,j,k)-tmp1(i,j-1,k))*dyiv              &
     &          +(tmp4(i,j,k)-tmp4(i,j,k-1))*dziv)*tmp3(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

      end if

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1+jsouth,nj-jnorth
        do i=1,ni-1
          s32(i,j,k)=s23(i,j,k)
        end do
        end do

!$omp end do

      end do

! -----

!! -----

!$omp end parallel

!!! -----

! Set the boundary conditions.

      call bc8u(idwbc,idebc,ni,nj,nk,s12)
      call bc8v(idsbc,idnbc,ni,nj,nk,s12)

      call bc8u(idwbc,idebc,ni,nj,nk,s31)
      call bc8v(idsbc,idnbc,ni,nj,nk,s32)

      call bcten(idbbc,idtbc,ni,nj,nk,s31)
      call bcten(idbbc,idtbc,ni,nj,nk,s32)

! -----

      end subroutine s_defomten

!-----7--------------------------------------------------------------7--

      end module m_defomten
