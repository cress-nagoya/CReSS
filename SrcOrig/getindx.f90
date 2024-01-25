!***********************************************************************
      module m_getindx
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/04/06
!     Modification: 1999/05/10, 1999/07/05, 1999/09/30, 2000/01/17,
!                   2001/01/15, 2002/07/03, 2003/04/30, 2003/05/19,
!                   2006/09/21, 2008/05/02, 2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the maximum and minimum indices of do loops for optional
!     variable or data grid.

!-----7--------------------------------------------------------------7--

! Module reference

!     none

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getindx, s_getindx

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getindx

        module procedure s_getindx

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
      subroutine s_getindx(xo,imin,imax,jmin,jmax,istr,iend,jstr,jend)
!***********************************************************************

! Input variables

      character(len=2), intent(in) :: xo
                       ! Control flag of variable arrangement

      integer, intent(in) :: imin
                       ! Minimum array index in x direction

      integer, intent(in) :: imax
                       ! Maximum array index in x direction

      integer, intent(in) :: jmin
                       ! Minimum array index in y direction

      integer, intent(in) :: jmax
                       ! Maximum array index in y direction

! Output variables

      integer, intent(out) :: istr
                       ! Minimum do loops index in x direction

      integer, intent(out) :: iend
                       ! Maximum do loops index in x direction

      integer, intent(out) :: jstr
                       ! Minimum do loops index in y direction

      integer, intent(out) :: jend
                       ! Maximum do loops index in y direction

!-----7--------------------------------------------------------------7--

! Get the maximum and minimim indices of do loops for the data grid.

      if(xo(1:2).eq.'oo') then
        istr=imin
        iend=imax
        jstr=jmin
        jend=jmax

! -----

! Get the maximum and minimim indices of do loops for optional variable
! at the u points.

      else if(xo(1:2).eq.'ox') then
        istr=imin+1
        iend=imax-1
        jstr=jmin
        jend=jmax-1

! -----

! Get the maximum and minimim indices of do loops for optional variable
! at the v points.

      else if(xo(1:2).eq.'xo') then
        istr=imin
        iend=imax-1
        jstr=jmin+1
        jend=jmax-1

! -----

! Get the maximum and minimim indices of do loops for optional variable
! at the scalar points.

      else if(xo(1:2).eq.'xx') then
        istr=imin
        iend=imax-1
        jstr=jmin
        jend=jmax-1

! -----

! The other case.

      else
        istr=imin+1
        iend=imax-2
        jstr=jmin+1
        jend=jmax-2

      end if

! -----

      end subroutine s_getindx

!-----7--------------------------------------------------------------7--

      end module m_getindx
