!***********************************************************************
      module m_shedding
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2001/06/29, 2001/10/18, 2001/11/20,
!                   2002/01/15, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2003/12/12, 2004/09/01, 2004/09/25, 2004/10/12,
!                   2004/12/17, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the shedding rate from the snow and graupel to the rain
!     water.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: shedding, s_shedding

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface shedding

        module procedure s_shedding

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
      subroutine s_shedding(thresq,ni,nj,nk,qs,qg,t,clcs,clcg,clrs,clrg,&
     &                      clig,clsg,pgwet,shsr,shgr)
!***********************************************************************

! Input variables

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

      real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: clcs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between cloud water and snow

      real, intent(in) :: clcg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between cloud water and graupel

      real, intent(in) :: clrs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between rain water and snow

      real, intent(in) :: clrg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between rain water and graupel

      real, intent(in) :: clig(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between cloud ice and graupel

      real, intent(in) :: clsg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between snow and graupel

      real, intent(in) :: pgwet(0:ni+1,0:nj+1,1:nk)
                       ! Graupel production rate for moist process

! Output variables

      real, intent(out) :: shsr(0:ni+1,0:nj+1,1:nk)
                       ! Shedding rate of liquid water from snow

      real, intent(out) :: shgr(0:ni+1,0:nj+1,1:nk)
                       ! Shedding rate of liquid water from graupel

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

!!! Calculate the shedding rate.

!$omp parallel default(shared) private(k)

!! In the case nk = 1.

      if(nk.eq.1) then

!$omp do schedule(runtime) private(i,j)

        do j=1,nj-1
        do i=1,ni-1

! Calculate the shedding rate from the snow to the rain water.

          if(qs(i,j,1).gt.thresq) then

            if(t(i,j,1).ge.t0) then

              shsr(i,j,1)=clcs(i,j,1)+clrs(i,j,1)

            else

              shsr(i,j,1)=0.e0

            end if

          else

            shsr(i,j,1)=0.e0

          end if

! -----

! Calculate the shedding rate from the graupel to the rain water.

          if(qg(i,j,1).gt.thresq) then

            if(t(i,j,1).ge.t0) then

              shgr(i,j,1)=clcg(i,j,1)+clrg(i,j,1)

            else

              if(pgwet(i,j,1).gt.0.e0) then

                shgr(i,j,1)=(clcg(i,j,1)+clrg(i,j,1)                    &
     &            +clig(i,j,1)+clsg(i,j,1))-pgwet(i,j,1)

              else

                shgr(i,j,1)=0.e0

              end if

            end if

          else

            shgr(i,j,1)=0.e0

          end if

! -----

        end do
        end do

!$omp end do

!! -----

!! In the case nk > 1.

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

! Calculate the shedding rate from the snow to the rain water.

            if(qs(i,j,k).gt.thresq) then

              if(t(i,j,k).ge.t0) then

                shsr(i,j,k)=clcs(i,j,k)+clrs(i,j,k)

              else

                shsr(i,j,k)=0.e0

              end if

            else

              shsr(i,j,k)=0.e0

            end if

! -----

! Calculate the shedding rate from the graupel to the rain water.

            if(qg(i,j,k).gt.thresq) then

              if(t(i,j,k).ge.t0) then

                shgr(i,j,k)=clcg(i,j,k)+clrg(i,j,k)

              else

                if(pgwet(i,j,k).gt.0.e0) then

                  shgr(i,j,k)=(clcg(i,j,k)+clrg(i,j,k)                  &
     &              +clig(i,j,k)+clsg(i,j,k))-pgwet(i,j,k)

                else

                  shgr(i,j,k)=0.e0

                end if

              end if

            else

              shgr(i,j,k)=0.e0

            end if

! -----

          end do
          end do

!$omp end do

        end do

      end if

!! -----

!$omp end parallel

!!! -----

      end subroutine s_shedding

!-----7--------------------------------------------------------------7--

      end module m_shedding
