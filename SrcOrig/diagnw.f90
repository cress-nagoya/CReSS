!***********************************************************************
      module m_diagnw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/11/26
!     Modification: 2008/05/02, 2008/08/25, 2009/01/30, 2009/02/27,
!                   2009/11/05, 2011/03/18, 2011/09/22, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the diagnostic concentrations of the water hydrometeor.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: diagnw, s_diagnw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface diagnw

        module procedure s_diagnw

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max
      intrinsic min
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_diagnw(ni,nj,nk,nqw,nnw,rbr,qwtr,nwdia)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentrations

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor

! Output variable

      real, intent(out) :: nwdia(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Diagnostic concentrations of water hydrometeor

! Internal shared variables

      real mrmiv2      ! 0.01 / mrmax
      real mr0iv2      ! 100.0 / mr0

      real cdiaqr      ! Coefficient of mean diameter of rain water

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real rbv         ! Inverse of base state density

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      mrmiv2=1.e-2/mrmax
      mr0iv2=1.e2/mr0

      cdiaqr=nr0*nr0*nr0/(cc*rhow)

! -----

!! Get the diagnostic concentrations of the water hydrometeor.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j,rbv)

        do j=1,nj-1
        do i=1,ni-1

! Calculate the inverse of base state density.

          rbv=1.e0/rbr(i,j,k)

! -----

! Get the diagnostic concentrations of cloud water.

          nwdia(i,j,k,1)=nclcst*rbv

! -----

! Get the diagnostic concentrations of rain water.

          nwdia(i,j,k,2)=sqrt(                                          &
     &      sqrt(cdiaqr*rbr(i,j,k)*qwtr(i,j,k,2)))*rbv

          nwdia(i,j,k,2)=min(max(                                       &
     &      nwdia(i,j,k,2),mrmiv2*qwtr(i,j,k,2)),mr0iv2*qwtr(i,j,k,2))

! -----

        end do
        end do

!$omp end do

      end do

!$omp end parallel

!! -----

      end subroutine s_diagnw

!-----7--------------------------------------------------------------7--

      end module m_diagnw
