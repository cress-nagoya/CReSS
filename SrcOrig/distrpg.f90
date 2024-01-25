!***********************************************************************
      module m_distrpg
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2001/05/29, 2001/06/29, 2001/11/20,
!                   2002/04/02, 2002/09/09, 2003/04/30, 2003/05/19,
!                   2003/10/31, 2004/06/10, 2004/09/01, 2004/09/25,
!                   2004/10/12, 2004/12/17, 2005/04/04, 2005/09/30,
!                   2007/10/19, 2007/11/26, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2011/03/18, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the distribution ratio at which the collisions between
!     rain water and snow and reset the collection rate.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: distrpg, s_distrpg

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface distrpg

        module procedure s_distrpg

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_distrpg(cphopt,thresq,ni,nj,nk,qs,t,diaqr,diaqs,     &
     &                     clrs,clsr,clrsn,clsrn,clrsg)
!***********************************************************************

! Input variables

      integer, intent(in) :: cphopt
                       ! Option for cloud micro physics

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: thresq
                       ! Minimum threshold value of mixing ratio

      real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: diaqr(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of rain water

      real, intent(in) :: diaqs(0:ni+1,0:nj+1,1:nk)
                       ! Mean diameter of snow

! Input and output variables

      real, intent(inout) :: clrs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate from rain water to snow

      real, intent(inout) :: clsr(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate from snow to rain water

      real, intent(inout) :: clrsn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between rain water and snow

      real, intent(inout) :: clsrn(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate for concentrations
                       ! between rain water and snow

! Output variable

      real, intent(out) :: clrsg(0:ni+1,0:nj+1,1:nk)
                       ! Production rate of graupel
                       ! from collection rate form rain to snow

! Internal shared variables

      real rhos2       ! rhos x rhos
      real rhow2       ! rhow x rhow

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real alpha       ! Distribution ratio between snow and graupel

      real alpha1      ! 1.0 - alpha

      real a           ! Temporary variable
      real b           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      rhos2=rhos*rhos
      rhow2=rhow*rhow

! -----

!!! Calculate the distribution ratio at which the collisions between
!!! rain water and snow and reset the collection rate.

!$omp parallel default(shared) private(k)

!! In the case nk = 1.

      if(nk.eq.1) then

! Perform calculating in the case the option abs(cphopt) is equal to 2.

        if(abs(cphopt).eq.2) then

!$omp do schedule(runtime) private(i,j,alpha,alpha1,a,b)

          do j=1,nj-1
          do i=1,ni-1

            if(qs(i,j,1).gt.thresq) then

              if(t(i,j,1).lt.t0) then

                a=diaqs(i,j,1)*diaqs(i,j,1)
                b=diaqr(i,j,1)*diaqr(i,j,1)

                a=rhos2*a*a*a
                b=rhow2*b*b*b

                alpha=a/(a+b)
                alpha1=1.e0-alpha

                clrsg(i,j,1)=alpha1*clrs(i,j,1)

                clrs(i,j,1)=alpha*clrs(i,j,1)

                clsr(i,j,1)=alpha1*clsr(i,j,1)

              else

                clrsg(i,j,1)=clrs(i,j,1)

                clrs(i,j,1)=0.e0

              end if

            else

              clrsg(i,j,1)=0.e0

            end if

          end do
          end do

!$omp end do

! -----

! Perform calculating in the case the option abs(cphopt) is greater
! than 2.

        else if(abs(cphopt).ge.3) then

!$omp do schedule(runtime) private(i,j,alpha,alpha1,a,b)

          do j=1,nj-1
          do i=1,ni-1

            if(qs(i,j,1).gt.thresq) then

              if(t(i,j,1).lt.t0) then

                a=diaqs(i,j,1)*diaqs(i,j,1)
                b=diaqr(i,j,1)*diaqr(i,j,1)

                a=rhos2*a*a*a
                b=rhow2*b*b*b

                alpha=a/(a+b)
                alpha1=1.e0-alpha

                clrsg(i,j,1)=alpha1*clrs(i,j,1)

                clrs(i,j,1)=alpha*clrs(i,j,1)

                clsr(i,j,1)=alpha1*clsr(i,j,1)

                clrsn(i,j,1)=alpha1*clrsn(i,j,1)
                clsrn(i,j,1)=alpha1*clsrn(i,j,1)

              else

                clrsg(i,j,1)=clrs(i,j,1)

                clrs(i,j,1)=0.e0

              end if

            else

              clrsg(i,j,1)=0.e0

            end if

          end do
          end do

!$omp end do

        end if

! -----

!! -----

!! In the case nk > 1.

      else

! Perform calculating in the case the option abs(cphopt) is equal to 2.

        if(abs(cphopt).eq.2) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,alpha,alpha1,a,b)

            do j=1,nj-1
            do i=1,ni-1

              if(qs(i,j,k).gt.thresq) then

                if(t(i,j,k).lt.t0) then

                  a=diaqs(i,j,k)*diaqs(i,j,k)
                  b=diaqr(i,j,k)*diaqr(i,j,k)

                  a=rhos2*a*a*a
                  b=rhow2*b*b*b

                  alpha=a/(a+b)
                  alpha1=1.e0-alpha

                  clrsg(i,j,k)=alpha1*clrs(i,j,k)

                  clrs(i,j,k)=alpha*clrs(i,j,k)

                  clsr(i,j,k)=alpha1*clsr(i,j,k)

                else

                  clrsg(i,j,k)=clrs(i,j,k)

                  clrs(i,j,k)=0.e0

                end if

              else

                clrsg(i,j,k)=0.e0

              end if

            end do
            end do

!$omp end do

          end do

! -----

! Perform calculating in the case the option abs(cphopt) is greater
! than 2.

        else if(abs(cphopt).ge.3) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,alpha,alpha1,a,b)

            do j=1,nj-1
            do i=1,ni-1

              if(qs(i,j,k).gt.thresq) then

                if(t(i,j,k).lt.t0) then

                  a=diaqs(i,j,k)*diaqs(i,j,k)
                  b=diaqr(i,j,k)*diaqr(i,j,k)

                  a=rhos2*a*a*a
                  b=rhow2*b*b*b

                  alpha=a/(a+b)
                  alpha1=1.e0-alpha

                  clrsg(i,j,k)=alpha1*clrs(i,j,k)

                  clrs(i,j,k)=alpha*clrs(i,j,k)

                  clsr(i,j,k)=alpha1*clsr(i,j,k)

                  clrsn(i,j,k)=alpha1*clrsn(i,j,k)
                  clsrn(i,j,k)=alpha1*clsrn(i,j,k)

                else

                  clrsg(i,j,k)=clrs(i,j,k)

                  clrs(i,j,k)=0.e0

                end if

              else

                clrsg(i,j,k)=0.e0

              end if

            end do
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

!$omp end parallel

!!! -----

      end subroutine s_distrpg

!-----7--------------------------------------------------------------7--

      end module m_distrpg
