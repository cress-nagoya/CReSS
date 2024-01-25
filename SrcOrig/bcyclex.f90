!***********************************************************************
      module m_bcyclex
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/02/10
!     Modification: 2006/12/04, 2007/01/05, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the periodic boundary conditions in x direction.

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

      public :: bcyclex, s_bcyclex

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcyclex

        module procedure s_bcyclex

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
      subroutine s_bcyclex(fpwbc,fpebc,iwsnd,iwrcv,iesnd,iercv,         &
     &                     ni,nj,kmax,var)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: iwsnd
                       ! Array index of west sended value

      integer, intent(in) :: iwrcv
                       ! Array index of west received value

      integer, intent(in) :: iesnd
                       ! Array index of east sended value

      integer, intent(in) :: iercv
                       ! Array index of east received value

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

! Input and output variable

      real, intent(inout) :: var(0:ni+1,0:nj+1,1:kmax)
                       ! Optional variable

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

! Internal private variables

      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)

! -----

! Set the periodic boundary conditions in x direction.

!$omp parallel default(shared) private(k)

      if(nisub.eq.1) then

        if(wbc.eq.-1.and.ebc.eq.-1) then

          do k=1,kmax

!$omp do schedule(runtime) private(j)

            do j=0,nj+1
              var(iwrcv,j,k)=var(iesnd,j,k)
            end do

!$omp end do

          end do

          do k=1,kmax

!$omp do schedule(runtime) private(j)

            do j=0,nj+1
              var(iercv,j,k)=var(iwsnd,j,k)
            end do

!$omp end do

          end do

        end if

      end if

!$omp end parallel

! -----

      end subroutine s_bcyclex

!-----7--------------------------------------------------------------7--

      end module m_bcyclex
