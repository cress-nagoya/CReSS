!***********************************************************************
      module m_bcycley
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/02/10
!     Modification: 2006/12/04, 2007/01/05, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the periodic boundary conditions in y direction.

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

      public :: bcycley, s_bcycley

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcycley

        module procedure s_bcycley

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
      subroutine s_bcycley(fpsbc,fpnbc,jssnd,jsrcv,jnsnd,jnrcv,         &
     &                     ni,nj,kmax,var)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: jssnd
                       ! Array index of south sended value

      integer, intent(in) :: jsrcv
                       ! Array index of south received value

      integer, intent(in) :: jnsnd
                       ! Array index of north sended value

      integer, intent(in) :: jnrcv
                       ! Array index of north received value

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

      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

! Internal private variables

      integer i        ! Array index in x direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)

! -----

! Set the periodic boundary conditions in y direction.

!$omp parallel default(shared) private(k)

      if(njsub.eq.1) then

        if(sbc.eq.-1.and.nbc.eq.-1) then

          do k=1,kmax

!$omp do schedule(runtime) private(i)

            do i=0,ni+1
              var(i,jsrcv,k)=var(i,jnsnd,k)
            end do

!$omp end do

          end do

          do k=1,kmax

!$omp do schedule(runtime) private(i)

            do i=0,ni+1
              var(i,jnrcv,k)=var(i,jssnd,k)
            end do

!$omp end do

          end do

        end if

      end if

!$omp end parallel

! -----

      end subroutine s_bcycley

!-----7--------------------------------------------------------------7--

      end module m_bcycley
