!***********************************************************************
      module m_smooasl
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     smooth the interpolated aerosol data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_getiname
      use m_gsmoos

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: smooasl, s_smooasl

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface smooasl

        module procedure s_smooasl

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
      subroutine s_smooasl(fpgsmcnt,ni,nj,nk,nqa,qagpv,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgsmcnt
                       ! Formal parameter of unique index of gsmcnt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

! Input and output variable

      real, intent(inout) :: qagpv(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Mixing ratio of aerosol data at marked time

! Internal shared variables

      integer gsmcnt   ! Iteration count

      integer iit      ! Index of iteration

      integer n        ! Array index in 4th direction

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpgsmcnt,gsmcnt)

! -----

! Smooth the interpolated aerosol data.

      do iit=1,gsmcnt

        do n=1,nqa(0)

          call s_gsmoos(idgsmcoe,ni,nj,nk,qagpv(0,0,1,n),tmp1)

        end do

      end do

! -----

      end subroutine s_smooasl

!-----7--------------------------------------------------------------7--

      end module m_smooasl
