!***********************************************************************
      module m_adjstnsg
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2001/10/18, 2001/11/20, 2002/01/07,
!                   2002/01/15, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2003/10/31, 2003/12/12, 2004/04/01, 2004/05/31,
!                   2004/06/10, 2004/09/01, 2004/09/25, 2004/10/12,
!                   2004/12/17, 2005/04/04, 2005/09/30, 2006/04/03,
!                   2007/10/19, 2007/11/26, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2009/11/05, 2011/03/18, 2011/09/22,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     adjust the concentrations of the precipitation categories of the
!     ice hydrometeor for the given mixing ratio.

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

      public :: adjstnsg, s_adjstnsg

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface adjstnsg

        module procedure s_adjstnsg

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
      subroutine s_adjstnsg(fphaiopt,ni,nj,nk,                          &
     &                      rbr,rbv,qs,qg,qh,ncs,ncg,nch)
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

      real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio

      real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio

      real, intent(in) :: qh(0:ni+1,0:nj+1,1:nk)
                       ! Hail mixing ratio

! Input and output variables

      real, intent(inout) :: ncs(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of snow

      real, intent(inout) :: ncg(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of graupel

      real, intent(inout) :: nch(0:ni+1,0:nj+1,1:nk)
                       ! Concentrations of hail

! Internal shared variables

      integer haiopt   ! Option for additional hail processes

      real msmiv2      ! 0.01 / msmax
      real ms0iv2      ! 100.0 / ms0

      real mgmiv2      ! 0.01 / mgmax
      real mg0iv2      ! 100.0 / mg0

      real mhmiv2      ! 0.01 / mhmax
      real mh0iv2      ! 100.0 / mh0

      real cdiaqs      ! Coefficient of mean diameter of snow
      real cdiaqg      ! Coefficient of mean diameter of graupel
      real cdiaqh      ! Coefficient of mean diameter of hail

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real ndia        ! Diagnostic concentrations of
                       ! snow, graupel or hail

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fphaiopt,haiopt)

! -----

! Set the common used variables.

      msmiv2=1.e-2/msmax
      ms0iv2=1.e2/ms0

      mgmiv2=1.e-2/mgmax
      mg0iv2=1.e2/mg0

      mhmiv2=1.e-2/mhmax
      mh0iv2=1.e2/mh0

      cdiaqs=ns0*ns0*ns0/(cc*rhos)
      cdiaqg=ng0*ng0*ng0/(cc*rhog)
      cdiaqh=nh0*nh0*nh0/(cc*rhoh)

! -----

!!! Adjust the concentrations of the precipitation categories of the ice
!!! hydrometeor.

!$omp parallel default(shared) private(k)

!! Adjust the concentrations of the snow and graupel.

      if(haiopt.eq.0) then

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ndia)

          do j=1,nj-1
          do i=1,ni-1

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

!! Adjust the concentrations of the snow, graupel and hail.

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ndia)

          do j=1,nj-1
          do i=1,ni-1

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

      end subroutine s_adjstnsg

!-----7--------------------------------------------------------------7--

      end module m_adjstnsg
