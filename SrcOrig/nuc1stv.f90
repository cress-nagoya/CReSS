!***********************************************************************
      module m_nuc1stv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/07/05
!     Modification: 2000/08/21, 2000/10/18, 2001/06/29, 2001/10/18,
!                   2001/11/20, 2001/12/11, 2002/01/07, 2002/01/15,
!                   2002/04/02, 2002/12/06, 2003/01/04, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/11/05, 2003/12/12,
!                   2004/03/22, 2004/08/01, 2004/09/01, 2004/09/10,
!                   2004/10/12, 2004/12/17, 2005/04/04, 2006/02/13,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/11/13, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the nucleation rate of the deposition or sorption.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: nuc1stv, s_nuc1stv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface nuc1stv

        module procedure s_nuc1stv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic max
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_nuc1stv(thresq,ni,nj,nk,rbv,qv,qi,t,qvsi,nuvi)
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

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio

      real, intent(in) :: t(0:ni+1,0:nj+1,1:nk)
                       ! Air temperature

      real, intent(in) :: qvsi(0:ni+1,0:nj+1,1:nk)
                       ! Saturation mixing ratio for ice

! Output variable

      real, intent(out) :: nuvi(0:ni+1,0:nj+1,1:nk)
                       ! Nucleation rate of deposition or sorption

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real qvssi       ! Super saturation mixing ratio for ice

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

!! Calculate the nucleation rate of the deposition or sorption.

!$omp parallel default(shared) private(k)

! In the case nk = 1.

      if(nk.eq.1) then

!$omp do schedule(runtime) private(i,j,qvssi,a)

        do j=1,nj-1
        do i=1,ni-1

          if(qv(i,j,1).gt.thresq) then

            if(qv(i,j,1).gt.qvsi(i,j,1)                                 &
     &        .and.t(i,j,1).gt.tlow.and.t(i,j,1).lt.t0) then

              qvssi=qv(i,j,1)-qvsi(i,j,1)

              a=15.25e0*qv(i,j,1)/qvsi(i,j,1)-10.08e0

              if(a.lt.40.e0) then

                nuvi(i,j,1)                                             &
     &            =max(min(mi0*exp(a)*rbv(i,j,1)-qi(i,j,1),qvssi),0.e0)

              else

                nuvi(i,j,1)=qvssi

              end if

            else

              nuvi(i,j,1)=0.e0

            end if

          else

            nuvi(i,j,1)=0.e0

          end if

        end do
        end do

!$omp end do

! -----

! In the case nk > 1.

      else

        do k=1,nk-1

!$omp do schedule(runtime) private(i,j,qvssi,a)

          do j=1,nj-1
          do i=1,ni-1

            if(qv(i,j,k).gt.thresq) then

              if(qv(i,j,k).gt.qvsi(i,j,k)                               &
     &          .and.t(i,j,k).gt.tlow.and.t(i,j,k).lt.t0) then

                qvssi=qv(i,j,k)-qvsi(i,j,k)

                a=15.25e0*qv(i,j,k)/qvsi(i,j,k)-10.08e0

                if(a.lt.40.e0) then

                  nuvi(i,j,k)=max(0.e0,                                 &
     &              min(mi0*exp(a)*rbv(i,j,k)-qi(i,j,k),qvssi))

                else

                  nuvi(i,j,k)=qvssi

                end if

              else

                nuvi(i,j,k)=0.e0

              end if

            else

              nuvi(i,j,k)=0.e0

            end if

          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_nuc1stv

!-----7--------------------------------------------------------------7--

      end module m_nuc1stv
