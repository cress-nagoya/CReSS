!***********************************************************************
      module m_bc8u
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/07/05,
!                   1999/07/28, 1999/08/03, 1999/08/18, 1999/08/23,
!                   1999/09/30, 1999/10/07, 1999/11/01, 2000/01/17,
!                   2001/12/11, 2002/04/02, 2003/04/30, 2003/05/19,
!                   2004/08/20, 2006/12/04, 2007/01/05, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the boundary conditions for optional variable at the u points.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: bc8u, s_bc8u

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bc8u

        module procedure s_bc8u

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
      subroutine s_bc8u(fpwbc,fpebc,ni,nj,kmax,var8u)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

! Input and output variable

      real, intent(inout) :: var8u(0:ni+1,0:nj+1,1:kmax)
                       ! Optional variable at u points

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer nim1     ! ni - 1
      integer nim2     ! ni - 2

! Internal private variables

      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)

! -----

! Set the common used variables.

      nim1=ni-1
      nim2=ni-2

! -----

!! Set the west and east boundary conditions.

!$omp parallel default(shared) private(k)

! Set the west boundary conditions.

      if(ebw.eq.1.and.isub.eq.0) then

        if(wbc.eq.2) then

          do k=1,kmax

!$omp do schedule(runtime) private(j)

            do j=0,nj+1
              var8u(1,j,k)=var8u(3,j,k)
            end do

!$omp end do

          end do

        else if(wbc.ge.3) then

          do k=1,kmax

!$omp do schedule(runtime) private(j)

            do j=0,nj+1
              var8u(1,j,k)=var8u(2,j,k)
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Set the east boundary conditions.

      if(ebe.eq.1.and.isub.eq.nisub-1) then

        if(ebc.eq.2) then

          do k=1,kmax

!$omp do schedule(runtime) private(j)

            do j=0,nj+1
              var8u(ni,j,k)=var8u(nim2,j,k)
            end do

!$omp end do

          end do

        else if(ebc.ge.3) then

          do k=1,kmax

!$omp do schedule(runtime) private(j)

            do j=0,nj+1
              var8u(ni,j,k)=var8u(nim1,j,k)
            end do

!$omp end do

          end do

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_bc8u

!-----7--------------------------------------------------------------7--

      end module m_bc8u
