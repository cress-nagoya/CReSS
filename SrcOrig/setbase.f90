!***********************************************************************
      module m_setbase
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/05/10,
!                   1999/05/20, 1999/06/28, 1999/07/05, 1999/08/03,
!                   1999/08/18, 1999/08/23, 1999/09/30, 1999/10/12,
!                   1999/11/01, 2000/01/17, 2000/03/08, 2000/07/05,
!                   2001/01/15, 2001/04/15, 2001/05/29, 2001/06/29,
!                   2001/10/17, 2002/04/02, 2002/06/18, 2002/08/15,
!                   2002/09/09, 2003/04/30, 2003/05/19, 2003/12/12,
!                   2004/01/09, 2004/04/15, 2004/09/10, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the base state variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bcbase
      use m_comindx
      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setbase, s_setbase

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setbase

        module procedure s_setbase

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_setbase(ni,nj,nk,zph,ubr,vbr,pbr,ptbr,qvbr,rbr,      &
     &                     zph8s,pibr,ptvbr)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

! Input and output variables

      real, intent(inout) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(inout) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(inout) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(inout) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(inout) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

! Output variable

      real, intent(out) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

! Internal shared variables

      real rddvcp      ! rd / cp

      real p0iv        ! 1.0 / p0

      real, intent(inout) :: zph8s(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at scalar points

      real, intent(inout) :: pibr(0:ni+1,0:nj+1,1:nk)
                       ! Base state Exner function

      real, intent(inout) :: ptvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state virtual potential temperature

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      rddvcp=rd/cp
      p0iv=1.e0/p0

! -----

!! Set the base state variables.

!$omp parallel default(shared) private(k)

! Calculate the z physical coordinates at the scalar, u and v points.

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          zph8s(i,j,k)=.5e0*(zph(i,j,k+1)+zph(i,j,k))
        end do
        end do

!$omp end do

      end do

! -----

! Get the base state Exner function and the density.

      do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

        do j=0,nj
        do i=0,ni
          ptvbr(i,j,k)=ptbr(i,j,k)                                      &
     &      *(1.e0+epsav*qvbr(i,j,k))/(1.e0+qvbr(i,j,k))

          pibr(i,j,k)=exp(rddvcp*log(p0iv*pbr(i,j,k)))

          rbr(i,j,k)=pbr(i,j,k)/(rd*ptvbr(i,j,k)*pibr(i,j,k))

        end do
        end do

!$omp end do

      end do

! -----

!$omp end parallel

!! -----

! Set the boundary conditions.

      call bcbase(idsmtopt,ni,nj,nk,                                    &
     &            zph8s,ubr,vbr,pbr,ptbr,qvbr,rbr,pibr,ptvbr)

! -----

      end subroutine s_setbase

!-----7--------------------------------------------------------------7--

      end module m_setbase
