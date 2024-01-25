!***********************************************************************
      module m_bc8v
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
!     set the boundary conditions for optional variable at the v points.

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

      public :: bc8v, s_bc8v

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bc8v

        module procedure s_bc8v

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
      subroutine s_bc8v(fpsbc,fpnbc,ni,nj,kmax,var8v)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

! Input and output variable

      real, intent(inout) :: var8v(0:ni+1,0:nj+1,1:kmax)
                       ! Optional variable at v points

! Internal shared variables

      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer njm1     ! nj - 1
      integer njm2     ! nj - 2

! Internal private variables

      integer i        ! Array index in x direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)

! -----

! Set the common used variables.

      njm1=nj-1
      njm2=nj-2

! -----

!! Set the south and north boundary conditions.

!$omp parallel default(shared) private(k)

! Set the south boundary conditions.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(sbc.eq.2) then

          do k=1,kmax

!$omp do schedule(runtime) private(i)

            do i=0,ni+1
              var8v(i,1,k)=var8v(i,3,k)
            end do

!$omp end do

          end do

        else if(sbc.ge.3) then

          do k=1,kmax

!$omp do schedule(runtime) private(i)

            do i=0,ni+1
              var8v(i,1,k)=var8v(i,2,k)
            end do

!$omp end do

          end do

        end if

      end if

! -----

! Set the north boundary conditions.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(nbc.eq.2) then

          do k=1,kmax

!$omp do schedule(runtime) private(i)

            do i=0,ni+1
              var8v(i,nj,k)=var8v(i,njm2,k)
            end do

!$omp end do

          end do

        else if(nbc.ge.3) then

          do k=1,kmax

!$omp do schedule(runtime) private(i)

            do i=0,ni+1
              var8v(i,nj,k)=var8v(i,njm1,k)
            end do

!$omp end do

          end do

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_bc8v

!-----7--------------------------------------------------------------7--

      end module m_bc8v
