!***********************************************************************
      module m_s2gpv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/03/13
!     Modification: 2001/04/15, 2001/05/29, 2001/06/29, 2002/04/02,
!                   2002/08/15, 2002/09/09, 2002/10/31, 2003/04/30,
!                   2003/05/19, 2003/11/28, 2003/12/12, 2006/01/10,
!                   2006/09/21, 2007/05/07, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2009/03/23, 2011/09/22,
!                   2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     perform the analysis nudging to GPV data of optional scalar
!     variable.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getcname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: s2gpv, s_s2gpv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface s2gpv

        module procedure s_s2gpv

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
      subroutine s_s2gpv(fpgpvvar,apg,nggdmp,gtinc,ni,nj,nk,rst,sp,     &
     &                   sgpv,std,sfrc)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: apg
                       ! Pointer of gpvvar

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: nggdmp
                       ! Analysis nudging damping coefficient

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jabobian

      real, intent(in) :: sp(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at past

      real, intent(in) :: sgpv(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable of GPV data
                       ! at marked time

      real, intent(in) :: std(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! optional scalar variable of GPV data

! Input and output variable

      real, intent(inout) :: sfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term of optional scalar

! Internal shared variable

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(gpvvar)

! -----

! Get the required namelist variable.

      call getcname(fpgpvvar,gpvvar)

! -----

! Calculate the analysis nudging terms for scalar variables.

      if(gpvvar(apg:apg).eq.'o') then

!$omp parallel default(shared) private(k)

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            sfrc(i,j,k)=sfrc(i,j,k)+nggdmp                              &
     &        *rst(i,j,k)*((sgpv(i,j,k)+std(i,j,k)*gtinc)-sp(i,j,k))
          end do
          end do

!$omp end do

        end do

!$omp end parallel

      end if

! -----

      end subroutine s_s2gpv

!-----7--------------------------------------------------------------7--

      end module m_s2gpv
