!***********************************************************************
      module m_vspqv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/20
!     Modification: 1999/07/05, 1999/08/03, 1999/09/30, 1999/10/12,
!                   1999/11/01, 2000/01/17, 2001/04/15, 2002/04/02,
!                   2002/08/15, 2003/04/30, 2003/05/19, 2003/06/27,
!                   2003/12/12, 2004/04/15, 2006/01/10, 2006/09/21,
!                   2007/05/07, 2007/07/30, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2008/12/11, 2009/02/27, 2009/03/23,
!                   2011/09/22, 2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the vertical sponge damping for the water vapor mixing
!     ratio.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getcname
      use m_getiname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vspqv, s_vspqv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vspqv

        module procedure s_vspqv

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
      subroutine s_vspqv(fpgpvvar,fpvspopt,ksp0,gtinc,ni,nj,nk,rst,     &
     &                   qvbr,qvp,rbct,qvgpv,qvtd,qvfrc)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fpvspopt
                       ! Formal parameter of unique index of vspopt

      integer, intent(in) :: ksp0(1:2)
                       ! Index of lowest vertical sponge level

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jabobian

      real, intent(in) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: rbct(1:ni,1:nj,1:nk,1:2)
                       ! Relaxed top sponge damping coefficients

      real, intent(in) :: qvgpv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, intent(in) :: qvtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! water vapor mixing ratio of GPV data

! Input and output variable

      real, intent(inout) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in water vapor mixing ratio equation

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      integer vspopt   ! Option for vertical sponge damping

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(gpvvar)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getiname(fpvspopt,vspopt)

! -----

!! Calculate the vertical sponge damping.

!$omp parallel default(shared) private(k)

! Damp to the GPV data.

      if(vspopt.eq.1.and.gpvvar(2:2).eq.'o') then

        do k=ksp0(1)-1,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            qvfrc(i,j,k)=qvfrc(i,j,k)                                   &
     &        -.5e0*(rbct(i,j,k,1)+rbct(i,j,k+1,1))                     &
     &        *rst(i,j,k)*(qvp(i,j,k)-(qvgpv(i,j,k)+qvtd(i,j,k)*gtinc))
          end do
          end do

!$omp end do

        end do

! -----

! Damp to the 0.

      else

        do k=ksp0(2)-1,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            qvfrc(i,j,k)=qvfrc(i,j,k)                                   &
     &        -.5e0*(rbct(i,j,k,2)+rbct(i,j,k+1,2))                     &
     &        *rst(i,j,k)*(qvp(i,j,k)-qvbr(i,j,k))
          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_vspqv

!-----7--------------------------------------------------------------7--

      end module m_vspqv