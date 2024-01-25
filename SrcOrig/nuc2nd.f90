!***********************************************************************
      module m_nuc2nd
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2001/04/15, 2001/05/29,
!                   2001/10/18, 2001/11/20, 2001/12/11, 2002/01/15,
!                   2002/04/02, 2002/12/02, 2003/03/21, 2003/04/30,
!                   2003/05/19, 2003/12/12, 2004/04/01, 2004/04/15,
!                   2004/06/10, 2004/08/01, 2004/09/01, 2004/09/25,
!                   2004/10/12, 2004/12/17, 2005/04/04, 2006/09/30,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the secondary nucleation rate of the ice crystals.

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

      public :: nuc2nd, s_nuc2nd

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface nuc2nd

        module procedure s_nuc2nd

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
      subroutine s_nuc2nd(ni,nj,nk,rbv,t,clcs,clcg,pgwet,spsi,spgi)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: clcs(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between cloud water and snow

      real, intent(in) :: clcg(0:ni+1,0:nj+1,1:nk)
                       ! Collection rate between cloud water and graupel

      real, intent(in) :: pgwet(0:ni+1,0:nj+1,1:nk)
                       ! Graupel production rate for moist process

! Output variables

      real, intent(out) :: spsi(0:ni+1,0:nj+1,1:nk)
                       ! Secondary nucleation rate from snow

      real, intent(out) :: spgi(0:ni+1,0:nj+1,1:nk)
                       ! Secondary nucleation rate from graupel

! Internal shared variables

      real mi0352      ! 3.5 x 1.0 x 10^8 x mi0 / 2.0
      real mi0353      ! 3.5 x 1.0 x 10^8 x mi0 / 3.0

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      mi0352=.5e0*3.5e8*mi0
      mi0353=oned3*3.5e8*mi0

! -----

!! Calculate the secondary nucleation rate.

!$omp parallel default(shared) private(k)

! In the case nk = 1.

      if(nk.eq.1) then

!$omp do schedule(runtime) private(i,j,a)

        do j=1,nj-1
        do i=1,ni-1

          if(pgwet(i,j,1).gt.0.e0) then

            if(t(i,j,1).gt.270.16e0) then

              spsi(i,j,1)=0.e0

            else if(t(i,j,1).gt.268.16e0.and.t(i,j,1).le.270.16e0) then

              spsi(i,j,1)                                               &
     &          =(270.16e0-t(i,j,1))*mi0352*rbv(i,j,1)*clcs(i,j,1)

            else if(t(i,j,1).gt.265.16e0.and.t(i,j,1).le.268.16e0) then

              spsi(i,j,1)                                               &
     &          =(t(i,j,1)-265.16e0)*mi0353*rbv(i,j,1)*clcs(i,j,1)

            else

              spsi(i,j,1)=0.e0

            end if

            spgi(i,j,1)=0.e0

          else

            if(t(i,j,1).gt.270.16e0) then

              spsi(i,j,1)=0.e0
              spgi(i,j,1)=0.e0

            else if(t(i,j,1).gt.268.16e0.and.t(i,j,1).le.270.16e0) then

              a=(270.16e0-t(i,j,1))*mi0352*rbv(i,j,1)

              spsi(i,j,1)=a*clcs(i,j,1)
              spgi(i,j,1)=a*clcg(i,j,1)

            else if(t(i,j,1).gt.265.16e0.and.t(i,j,1).le.268.16e0) then

              a=(t(i,j,1)-265.16e0)*mi0353*rbv(i,j,1)

              spsi(i,j,1)=a*clcs(i,j,1)
              spgi(i,j,1)=a*clcg(i,j,1)

            else

              spsi(i,j,1)=0.e0
              spgi(i,j,1)=0.e0

            end if

          end if

        end do
        end do

!$omp end do

! -----

! In the case nk > 1.

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,a)

          do j=1,nj-1
          do i=1,ni-1

            if(pgwet(i,j,k).gt.0.e0) then

              if(t(i,j,k).gt.270.16e0) then

                spsi(i,j,k)=0.e0

              else if(t(i,j,k).gt.268.16e0                              &
     &           .and.t(i,j,k).le.270.16e0) then

                spsi(i,j,k)                                             &
     &            =(270.16e0-t(i,j,k))*mi0352*rbv(i,j,k)*clcs(i,j,k)

              else if(t(i,j,k).gt.265.16e0                              &
     &           .and.t(i,j,k).le.268.16e0) then

                spsi(i,j,k)                                             &
     &            =(t(i,j,k)-265.16e0)*mi0353*rbv(i,j,k)*clcs(i,j,k)

              else

                spsi(i,j,k)=0.e0

              end if

              spgi(i,j,k)=0.e0

            else

              if(t(i,j,k).gt.270.16e0) then

                spsi(i,j,k)=0.e0
                spgi(i,j,k)=0.e0

              else if(t(i,j,k).gt.268.16e0                              &
     &           .and.t(i,j,k).le.270.16e0) then

                a=(270.16e0-t(i,j,k))*mi0352*rbv(i,j,k)

                spsi(i,j,k)=a*clcs(i,j,k)
                spgi(i,j,k)=a*clcg(i,j,k)

              else if(t(i,j,k).gt.265.16e0                              &
     &           .and.t(i,j,k).le.268.16e0) then

                a=(t(i,j,k)-265.16e0)*mi0353*rbv(i,j,k)

                spsi(i,j,k)=a*clcs(i,j,k)
                spgi(i,j,k)=a*clcg(i,j,k)

              else

                spsi(i,j,k)=0.e0
                spgi(i,j,k)=0.e0

              end if

            end if

          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_nuc2nd

!-----7--------------------------------------------------------------7--

      end module m_nuc2nd
