!***********************************************************************
      module m_adjstnp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/11/26
!     Modification: 2008/05/02, 2008/08/25, 2009/02/27, 2009/11/05,
!                   2011/03/18, 2011/09/22, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     adjust the concentrations of all of the precipitation categories
!     for the given mixing ratio.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: adjstnp, s_adjstnp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface adjstnp

        module procedure s_adjstnp

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
      subroutine s_adjstnp(fphaiopt,ni,nj,nk,                           &
     &                     rbr,rbv,qr,qs,qg,qh,ncr,ncs,ncg,nch)
!***********************************************************************

! Input variables

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

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

      real, intent(in) :: qr(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio

      real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio

      real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio

      real, intent(in) :: qh(0:ni+1,0:nj+1,1:nk)
                       ! Hail mixing ratio

! Input and output variables

      real, intent(inout) :: ncr(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of rain water

      real, intent(inout) :: ncs(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow

      real, intent(inout) :: ncg(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of graupel

      real, intent(inout) :: nch(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of hail

! Internal shared variables

      integer haiopt   ! Option for additional hail processes

      real mrmiv2      ! 0.01 / mrmax
      real mr0iv2      ! 100.0 / mr0

      real msmiv2      ! 0.01 / msmax
      real ms0iv2      ! 100.0 / ms0

      real mgmiv2      ! 0.01 / mgmax
      real mg0iv2      ! 100.0 / mg0

      real mhmiv2      ! 0.01 / mhmax
      real mh0iv2      ! 100.0 / mh0

      real cdiaqr      ! Coefficient of mean diameter of rain water
      real cdiaqs      ! Coefficient of mean diameter of snow
      real cdiaqg      ! Coefficient of mean diameter of graupel
      real cdiaqh      ! Coefficient of mean diameter of hail

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real ndia        ! Diagnostic concentrations of
                       ! rain water, snow, graupel or hail

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fphaiopt,haiopt)

! -----

! Set the common used variables.

      mrmiv2=1.e-2/mrmax
      mr0iv2=1.e2/mr0

      msmiv2=1.e-2/msmax
      ms0iv2=1.e2/ms0

      mgmiv2=1.e-2/mgmax
      mg0iv2=1.e2/mg0

      mhmiv2=1.e-2/mhmax
      mh0iv2=1.e2/mh0

      cdiaqr=nr0*nr0*nr0/(cc*rhow)
      cdiaqs=ns0*ns0*ns0/(cc*rhos)
      cdiaqg=ng0*ng0*ng0/(cc*rhog)
      cdiaqh=nh0*nh0*nh0/(cc*rhoh)

! -----

!!! Adjust the concentrations of all of the precipitation categories.

!$omp parallel default(shared) private(k)

!! Adjust the concentrations of the rain water, snow and graupel.

      if(haiopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ndia)

          do j=1,nj-1
          do i=1,ni-1

! Adjust the concentrations of the rain water.

            ndia=sqrt(sqrt(cdiaqr*rbr(i,j,k)*qr(i,j,k)))*rbv(i,j,k)

            ncr(i,j,k)=max(ncr(i,j,k),mrmiv2*qr(i,j,k),5.62341e-3*ndia)
            ncr(i,j,k)=min(ncr(i,j,k),mr0iv2*qr(i,j,k),1.77838e2*ndia)

! -----

! Adjust the concentrations of the snow.

            ndia=sqrt(sqrt(cdiaqs*rbr(i,j,k)*qs(i,j,k)))*rbv(i,j,k)

            ncs(i,j,k)=max(ncs(i,j,k),msmiv2*qs(i,j,k),5.62341e-3*ndia)
            ncs(i,j,k)=min(ncs(i,j,k),ms0iv2*qs(i,j,k),1.77838e2*ndia)

! -----

! Adjust the concentrations of the graupel.

            ndia=sqrt(sqrt(cdiaqg*rbr(i,j,k)*qg(i,j,k)))*rbv(i,j,k)

            ncg(i,j,k)=max(ncg(i,j,k),mgmiv2*qg(i,j,k),5.62341e-3*ndia)
            ncg(i,j,k)=min(ncg(i,j,k),mg0iv2*qg(i,j,k),1.77838e2*ndia)

! -----

          end do
          end do

!$omp end do

        end do

!! -----

!! Adjust the concentrations of the rain water, snow, graupel and hail.

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ndia)

          do j=1,nj-1
          do i=1,ni-1

! Adjust the concentrations of the rain water.

            ndia=sqrt(sqrt(cdiaqr*rbr(i,j,k)*qr(i,j,k)))*rbv(i,j,k)

            ncr(i,j,k)=max(ncr(i,j,k),mrmiv2*qr(i,j,k),5.62341e-3*ndia)
            ncr(i,j,k)=min(ncr(i,j,k),mr0iv2*qr(i,j,k),1.77838e2*ndia)

! -----

! Adjust the concentrations of the snow.

            ndia=sqrt(sqrt(cdiaqs*rbr(i,j,k)*qs(i,j,k)))*rbv(i,j,k)

            ncs(i,j,k)=max(ncs(i,j,k),msmiv2*qs(i,j,k),5.62341e-3*ndia)
            ncs(i,j,k)=min(ncs(i,j,k),ms0iv2*qs(i,j,k),1.77838e2*ndia)

! -----

! Adjust the concentrations of the graupel.

            ndia=sqrt(sqrt(cdiaqg*rbr(i,j,k)*qg(i,j,k)))*rbv(i,j,k)

            ncg(i,j,k)=max(ncg(i,j,k),mgmiv2*qg(i,j,k),5.62341e-3*ndia)
            ncg(i,j,k)=min(ncg(i,j,k),mg0iv2*qg(i,j,k),1.77838e2*ndia)

! -----

! Adjust the concentrations of the hail.

            ndia=sqrt(sqrt(cdiaqh*rbr(i,j,k)*qh(i,j,k)))*rbv(i,j,k)

            nch(i,j,k)=max(nch(i,j,k),mhmiv2*qh(i,j,k),5.62341e-3*ndia)
            nch(i,j,k)=min(nch(i,j,k),mh0iv2*qh(i,j,k),1.77838e2*ndia)

! -----

          end do
          end do

!$omp end do

        end do

      end if

!! -----

!$omp end parallel

!!! -----

      end subroutine s_adjstnp

!-----7--------------------------------------------------------------7--

      end module m_adjstnp
