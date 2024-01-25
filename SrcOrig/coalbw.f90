!***********************************************************************
      module m_coalbw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/08/08
!     Modification: 2006/09/30, 2006/11/27, 2007/10/19, 2008/01/11,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2009/03/12,
!                   2011/08/18, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the coalescence processes between water bins.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_remapbw

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: coalbw, s_coalbw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface coalbw

        module procedure s_coalbw

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
      subroutine s_coalbw(k,dtcoal,ni,nj,nk,nqw,nnw,                    &
     &                    bmw,rbmw,dbmw,ewbw,ubw,mwbin,nwbin,mwbrs,     &
     &                    bmwsc,bmwss,mwsc,nwsc,mwss,nwss,pct,tmp1)
!***********************************************************************

! Input variable

      integer, intent(in) :: k
                       ! Array index in z direction

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

      real, intent(in) :: dtcoal
                       ! Time interval of coalescence processes

      real, intent(in) :: bmw(1:nqw+1,1:3)
                       ! Mass at water bin boundary [g]

      real, intent(in) :: rbmw(1:nqw,1:2)
                       ! Related parameters of bmw

      real, intent(in) :: dbmw(1:nqw)
                       ! Differential between adjacent water bins [g]

      real, intent(in) :: ewbw(1:nqw,1:nqw)
                       ! Radius weighted coalescence efficiency
                       ! between water bins [cm^2]

      real, intent(in) :: ubw(0:ni+1,0:nj+1,1:nqw)
                       ! Terminal velocity of water bin [cm/s]

! Input and output variables

      real, intent(inout) :: mwbin(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Total water mass [g/cm^3]

      real, intent(inout) :: nwbin(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations [1/cm^3]

! Internal shared variables

      integer n        ! Array index in water bin categories
      integer ns       ! Array index in current coalescing water bin

      real, intent(inout) :: mwbrs(0:ni+1,0:nj+1)
                       ! Shifted mean water mass

      real, intent(inout) :: bmwsc(0:ni+1,0:nj+1,1:2)
                       ! Shifted mass at water bin boundary
                       ! in continuous coalescence

      real, intent(inout) :: bmwss(0:ni+1,0:nj+1,1:2)
                       ! Shifted mass at water bin boundary
                       ! in stochastic coalescence

      real, intent(inout) :: mwsc(0:ni+1,0:nj+1)
                       ! Shifted water mass
                       ! in continuous coalescence

      real, intent(inout) :: nwsc(0:ni+1,0:nj+1)
                       ! Shifted water concentrations
                       ! in continuous coalescence

      real, intent(inout) :: mwss(0:ni+1,0:nj+1,1:nqw)
                       ! Shifted water mass
                       ! in stochastic coalescence

      real, intent(inout) :: nwss(0:ni+1,0:nj+1,1:nnw)
                       ! Shifted water concentrations
                       ! in stochastic coalescence

      real, intent(inout) :: pct(0:ni+1,0:nj+1,1:nqw)
                       ! Total coalescence probability

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:4)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer n_sub    ! Substitute for n

      real pc          ! Collection probability

      real mwbr        ! Mean water mass

!-----7--------------------------------------------------------------7--

!!! Perform the coalescence processes.

      do ns=nqw,2,-1

!! Get current new value.

!$omp parallel default(shared) private(n_sub)

! Set the common used variables.

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

          if(mwbin(i,j,k,ns).gt.0.e0) then

            mwbrs(i,j)=mwbin(i,j,k,ns)/nwbin(i,j,k,ns)

            bmwsc(i,j,1)=bmw(ns,1)
            bmwsc(i,j,2)=bmw(ns+1,1)

            mwsc(i,j)=mwbrs(i,j)

            mwss(i,j,ns)=0.e0
            nwss(i,j,ns)=0.e0

            pct(i,j,ns)=0.e0

          end if

        end do
        end do

!$omp end do

! -----

! Calculate the coalescence parameters and perform continuous
! coalescence.

        do n_sub=ns-1,1,-1

!$omp do schedule(runtime) private(i,j,pc,mwbr)

          do j=1,nj-1
          do i=1,ni-1

            if(mwbin(i,j,k,ns).gt.0.e0) then

              if(mwbin(i,j,k,n_sub).gt.0.e0) then

                pc=nwbin(i,j,k,n_sub)                                   &
     &            *ewbw(ns,n_sub)*abs(ubw(i,j,ns)-ubw(i,j,n_sub))*dtcoal

                nwss(i,j,n_sub)=pc*nwbin(i,j,k,ns)

                if(nwss(i,j,n_sub).gt.nwbin(i,j,k,n_sub)) then

                  pc=nwbin(i,j,k,n_sub)/nwbin(i,j,k,ns)

                  nwss(i,j,n_sub)=nwbin(i,j,k,n_sub)

                end if

                pct(i,j,n_sub)=pct(i,j,n_sub+1)+pc

                if(pct(i,j,n_sub).gt.1.e0) then

                  mwbr=mwbin(i,j,k,n_sub)/nwbin(i,j,k,n_sub)

                  mwsc(i,j)=mwsc(i,j)+mwbr*pc

                  bmwsc(i,j,1)=bmwsc(i,j,1)+bmw(n_sub,1)*pc
                  bmwsc(i,j,2)=bmwsc(i,j,2)+bmw(n_sub+1,1)*pc

                  mwbin(i,j,k,n_sub)                                    &
     &              =mwbin(i,j,k,n_sub)-nwss(i,j,n_sub)*mwbr

                  nwbin(i,j,k,n_sub)=nwbin(i,j,k,n_sub)-nwss(i,j,n_sub)

                  if(mwbin(i,j,k,n_sub).le.0.e0                         &
     &              .or.nwbin(i,j,k,n_sub).le.0.e0) then

                    mwbin(i,j,k,n_sub)=0.e0
                    nwbin(i,j,k,n_sub)=0.e0

                  end if

                  if(pct(i,j,n_sub+1).le.1.e0) then
                    pct(i,j,ns)=pct(i,j,n_sub+1)
                  end if

                end if

              else

                pct(i,j,n_sub)=pct(i,j,n_sub+1)

              end if

            end if

          end do
          end do

!$omp end do

        end do

! -----

! Perform stochastic coalescence.

        do n_sub=ns-1,1,-1

!$omp do schedule(runtime) private(i,j,mwbr)

          do j=1,nj-1
          do i=1,ni-1

            if(mwbin(i,j,k,ns).gt.0.e0) then

              if(mwbin(i,j,k,n_sub).gt.0.e0) then

                if(pct(i,j,n_sub).le.1.e0) then

                  mwbr=nwss(i,j,n_sub)                                  &
     &              *mwbin(i,j,k,n_sub)/nwbin(i,j,k,n_sub)

                  mwss(i,j,n_sub)=nwss(i,j,n_sub)*mwsc(i,j)+mwbr

                  mwbin(i,j,k,n_sub)=mwbin(i,j,k,n_sub)-mwbr
                  nwbin(i,j,k,n_sub)=nwbin(i,j,k,n_sub)-nwss(i,j,n_sub)

                  if(mwbin(i,j,k,n_sub).le.0.e0                         &
     &              .or.nwbin(i,j,k,n_sub).le.0.e0) then

                    mwbin(i,j,k,n_sub)=0.e0
                    nwbin(i,j,k,n_sub)=0.e0

                  end if

                  mwss(i,j,ns)=mwss(i,j,ns)+nwss(i,j,n_sub)*mwbrs(i,j)
                  nwss(i,j,ns)=nwss(i,j,ns)+nwss(i,j,n_sub)

                else

                  mwss(i,j,n_sub)=0.e0

                end if

              else

                mwss(i,j,n_sub)=0.e0

              end if

            else

              mwss(i,j,n_sub)=0.e0

            end if

          end do
          end do

!$omp end do

        end do

! -----

! Distribute the remnants mass and concentrations after stochastic
! coalescence to the shifted bin for continuous coalescence and
! calculate the removed mass and concentrations for current
! coalescing bin.

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

          if(mwbin(i,j,k,ns).gt.0.e0) then

            if(pct(i,j,1).gt.1.e0) then

              nwsc(i,j)=(1.e0-pct(i,j,ns))*nwbin(i,j,k,ns)
              mwsc(i,j)=nwsc(i,j)*mwsc(i,j)

              mwbin(i,j,k,ns)=0.e0
              nwbin(i,j,k,ns)=0.e0

            else

              mwsc(i,j)=0.e0

              mwbin(i,j,k,ns)=mwbin(i,j,k,ns)-mwss(i,j,ns)
              nwbin(i,j,k,ns)=nwbin(i,j,k,ns)-nwss(i,j,ns)

              if(mwbin(i,j,k,ns).le.0.e0                                &
     &          .or.nwbin(i,j,k,ns).le.0.e0) then

                mwbin(i,j,k,ns)=0.e0
                nwbin(i,j,k,ns)=0.e0

              end if

            end if

          else

            mwsc(i,j)=0.e0

          end if

        end do
        end do

!$omp end do

! -----

!$omp end parallel

!! -----

!! Perform remapping for stochastic coalescense.

        do n=1,ns-1

! Get the shifted bin boundaries.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

            if(mwss(i,j,n).gt.0.e0) then

              if(pct(i,j,n).le.1.e0) then

                bmwss(i,j,1)=bmwsc(i,j,1)+bmw(n,1)
                bmwss(i,j,2)=bmwsc(i,j,2)+bmw(n+1,1)

              end if

            end if

          end do
          end do

!$omp end do

!$omp end parallel

! -----

! Perform remapping.

          call s_remapbw(k,ns,ni,nj,nk,nqw,nnw,1,1,bmw,rbmw,dbmw,       &
     &                  bmwss,mwss(0,0,n),nwss(0,0,n),mwbin,nwbin,      &
     &                  tmp1(0,0,1),tmp1(0,0,2),tmp1(0,0,3),tmp1(0,0,4))

! -----

        end do

!! -----

! Perform remapping for continuous coalescense.

        call s_remapbw(k,ns,ni,nj,nk,nqw,nnw,1,1,                       &
     &                 bmw,rbmw,dbmw,bmwsc,mwsc,nwsc,mwbin,nwbin,       &
     &                 tmp1(0,0,1),tmp1(0,0,2),tmp1(0,0,3),tmp1(0,0,4))

! -----

      end do

!!! -----

      end subroutine s_coalbw

!-----7--------------------------------------------------------------7--

      end module m_coalbw
