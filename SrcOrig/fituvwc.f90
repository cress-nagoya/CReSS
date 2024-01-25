!***********************************************************************
      module m_fituvwc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/06/01
!     Modification: 2000/12/18, 2001/08/07, 2001/11/20, 2002/04/02,
!                   2002/06/06, 2003/04/30, 2003/05/19, 2003/11/05,
!                   2004/09/10, 2005/02/10, 2006/09/30, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     fit the x and the y components of velocity and the zeta components
!     of contravariant velocity to the mass consistent equation by the
!     lamb.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc4news
      use m_bcycle
      use m_combuf
      use m_comindx
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getiname
      use m_getrname
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy
      use m_vbcwc

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: fituvwc, s_fituvwc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface fituvwc

        module procedure s_fituvwc

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
      subroutine s_fituvwc(fpiwest,fpieast,fpjsouth,fpjnorth,           &
     &                     fpdxiv,fpdyiv,fpdziv,fpalpha1,fpalpha2,      &
     &                     ni,nj,nk,jcb8u,jcb8v,jcb8w,lamb,u,v,wc)
!***********************************************************************

! Input variables

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

      integer, intent(in) :: fpalpha1
                       ! Formal parameter of unique index of alpha1

      integer, intent(in) :: fpalpha2
                       ! Formal parameter of unique index of alpha2

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

! Input and output variables

      real, intent(inout) :: lamb(0:ni+1,0:nj+1,1:nk)
                       ! Lagrange multiplier

      real, intent(inout) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(inout) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(inout) :: wc(0:ni+1,0:nj+1,1:nk)
                       ! zeta components of contravariant velocity

! Internal shared variables

      integer iwest    ! Added index on west boundary
      integer jsouth   ! Added index on south boundary

      integer ieast    ! Subtracted index on east boundary
      integer jnorth   ! Subtracted index on north boundary

      integer ni_sub   ! Substitute for ni
      integer nj_sub   ! Substitute for nj

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy
      real dziv        ! Inverse of dz

      real alpha1      ! Weighting coeffient
      real alpha2      ! Weighting coeffient

      real asdxiv      ! 0.5 x dxiv / (alpha1 x alpha1)
      real asdyiv      ! 0.5 x dyiv / (alpha1 x alpha1)
      real asdziv      ! 0.5 x dziv / (alpha2 x alpha2)

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpiwest,iwest)
      call getiname(fpieast,ieast)
      call getiname(fpjsouth,jsouth)
      call getiname(fpjnorth,jnorth)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)
      call getrname(fpdziv,dziv)
      call getrname(fpalpha1,alpha1)
      call getrname(fpalpha2,alpha2)

! -----

! Set the common used variables.

      asdxiv=.5e0*dxiv/(alpha1*alpha1)
      asdyiv=.5e0*dyiv/(alpha1*alpha1)
      asdziv=.5e0*dziv/(alpha2*alpha2)

! -----

! Set the substituted variables.

      ni_sub=ni
      nj_sub=nj

! -----

!! Exchange the value horizontally.

! Exchange the value horizontally between sub domain.

      call s_putbufsx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,lamb,1,1,sbuf)

      call s_shiftsx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

      call s_getbufsx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,lamb,1,1,rbuf)

      call s_putbufsy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,lamb,1,1,sbuf)

      call s_shiftsy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

      call s_getbufsy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,lamb,1,1,rbuf)

! -----

! Exchange the value horizontally between group domain.

      call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,lamb,1,1,sbuf)

      call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

      call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,lamb,1,1,rbuf)

      call s_putbufgy(idsbc,idnbc,'all',3,nj-3,ni,nj,nk,lamb,1,1,sbuf)

      call s_shiftgy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

      call s_getbufgy(idsbc,idnbc,'all',0,nj_sub,ni,nj,nk,lamb,1,1,rbuf)

      call s_putbufgx(idwbc,idebc,'all',3,ni-3,ni,nj,nk,lamb,1,1,sbuf)

      call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

      call s_getbufgx(idwbc,idebc,'all',0,ni_sub,ni,nj,nk,lamb,1,1,rbuf)

! -----

!! -----

! Set the periodic boundary conditions.

      call bcycle(idwbc,idebc,idsbc,idnbc,                              &
     &            3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,nk,lamb)

! -----

! Set the boundary conditions at the four corners.

      call bc4news(idwbc,idebc,idsbc,idnbc,0,ni_sub,0,nj_sub,           &
     &             ni,nj,nk,lamb)

! -----

!! Fit the x and the y components of velocity and the zeta components of
!! contravariant velocity to the mass consistent equation by the lamb.

!$omp parallel default(shared) private(k)

! Fit the x and y components of velocity.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=jsouth,nj-jnorth
        do i=1+iwest,ni-ieast
          u(i,j,k)                                                      &
     &      =u(i,j,k)+(lamb(i,j,k)-lamb(i-1,j,k))/jcb8u(i,j,k)*asdxiv
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

        do j=1+jsouth,nj-jnorth
        do i=iwest,ni-ieast
          v(i,j,k)                                                      &
     &      =v(i,j,k)+(lamb(i,j,k)-lamb(i,j-1,k))/jcb8v(i,j,k)*asdyiv
        end do
        end do

!$omp end do

      end do

! -----

! Fit the zeta components of contravariant velocity.

      do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1
          wc(i,j,k)                                                     &
     &      =wc(i,j,k)+(lamb(i,j,k)-lamb(i,j,k-1))/jcb8w(i,j,k)*asdziv
        end do
        end do

!$omp end do

      end do

! -----

!$omp end parallel

!! -----

! Set the bottom and the top boundary conditions.

      call vbcwc(idbbc,idtbc,ni,nj,nk,wc)

! -----

      end subroutine s_fituvwc

!-----7--------------------------------------------------------------7--

      end module m_fituvwc
