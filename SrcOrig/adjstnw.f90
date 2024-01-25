!***********************************************************************
      module m_adjstnw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/11/26
!     Modification: 2008/05/02, 2008/08/25, 2009/02/27, 2009/11/05,
!                   2011/03/18, 2011/04/06, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     adjust the concentrations of the cloud water and rain water for
!     the given mixing ratio.

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

      public :: adjstnw, s_adjstnw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface adjstnw

        module procedure s_adjstnw

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
      subroutine s_adjstnw(ni,nj,nk,rbr,rbv,qc,qr,ncc,ncr)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(in) :: qc(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio

      real, intent(in) :: qr(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio

! Input and output variables

      real, intent(inout) :: ncc(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of cloud water

      real, intent(inout) :: ncr(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of rain water

! Internal shared variables

      real mcmiv5      ! 0.5 / mcmax
      real mc0iv2      ! 100.0 / mc0

      real mrmiv2      ! 0.01 / mrmax
      real mr0iv2      ! 100.0 / mr0

      real cdiaqc      ! Coefficient of mean diameter of cloud water
      real cdiaqr      ! Coefficient of mean diameter of rain water

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real ndia        ! Diagnostic concentrations
                       ! of cloud water or rain water

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      mcmiv5=.5e0/mcmax
      mc0iv2=1.e2/mc0

      mrmiv2=1.e-2/mrmax
      mr0iv2=1.e2/mr0

      cdiaqc=nc0*nc0*nc0/(cc*rhow)
      cdiaqr=nr0*nr0*nr0/(cc*rhow)

! -----

!! Adjust the concentrations of the cloud water and rain water.

!$omp parallel default(shared) private(k)

      do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ndia)

        do j=1,nj-1
        do i=1,ni-1

! Adjust the concentrations of the cloud water.

!ORIG     ncc(i,j,k)                                                    &
!ORIG&      =min(max(ncc(i,j,k),mcmiv5*qc(i,j,k)),mc0iv2*qc(i,j,k))

          ndia=sqrt(sqrt(cdiaqc*rbr(i,j,k)*qc(i,j,k)))*rbv(i,j,k)

          ncc(i,j,k)=max(ncc(i,j,k),mcmiv5*qc(i,j,k),5.62341e-3*ndia)
          ncc(i,j,k)=min(ncc(i,j,k),mc0iv2*qc(i,j,k),1.77838e2*ndia)

! -----

! Adjust the concentrations of the rain water.

          ndia=sqrt(sqrt(cdiaqr*rbr(i,j,k)*qr(i,j,k)))*rbv(i,j,k)

          ncr(i,j,k)=max(ncr(i,j,k),mrmiv2*qr(i,j,k),5.62341e-3*ndia)
          ncr(i,j,k)=min(ncr(i,j,k),mr0iv2*qr(i,j,k),1.77838e2*ndia)

! -----

        end do
        end do

!$omp end do

      end do

!$omp end parallel

!! -----

      end subroutine s_adjstnw

!-----7--------------------------------------------------------------7--

      end module m_adjstnw
