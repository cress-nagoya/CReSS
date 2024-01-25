!***********************************************************************
      module m_rstindx
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/07/01
!     Modification: 2005/08/05, 2006/12/04, 2007/01/05, 2008/05/02,
!                   2008/08/25, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the minimum and maximum and differential do loops
!     indices to reposition.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rstindx, s_rstindx

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rstindx

        module procedure s_rstindx

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
      subroutine s_rstindx(ni,nj,ies,jes,ies_rst,jes_rst,               &
     &                     imin,imax,jmin,jmax,imin_rst,imax_rst,       &
     &                     jmin_rst,jmax_rst,istr,iend,jstr,jend,di,dj, &
     &                     istrb,iendb,jstrb,jendb,dib,djb)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: ies
                       ! Start index in entire domain in x direction

      integer, intent(in) :: jes
                       ! Start index in entire domain in y direction

      integer, intent(in) :: ies_rst
                       ! Start index in entire domain in x direction
                       ! for restructed files

      integer, intent(in) :: jes_rst
                       ! Start index in entire domain in y direction
                       ! for restructed files

      integer, intent(in) :: imin
                       ! Minimum index in group domain in x direction

      integer, intent(in) :: imax
                       ! Maximum index in group domain in x direction

      integer, intent(in) :: jmin
                       ! Minimum index in group domain in y direction

      integer, intent(in) :: jmax
                       ! Maximum index in group domain in y direction

      integer, intent(in) :: imin_rst
                       ! Minimum index in group domain in x direction
                       ! for restructed files

      integer, intent(in) :: imax_rst
                       ! Maximum index in group domain in x direction
                       ! for restructed files

      integer, intent(in) :: jmin_rst
                       ! Minimum index in group domain in y direction
                       ! for restructed files

      integer, intent(in) :: jmax_rst
                       ! Maximum index in group domain in y direction
                       ! for restructed files

! Output variables

      integer, intent(out) :: istr
                       ! Minimum do loops index in x direction

      integer, intent(out) :: iend
                       ! Maximum do loops index in x direction

      integer, intent(out) :: jstr
                       ! Minimum do loops index in y direction

      integer, intent(out) :: jend
                       ! Maximum do loops index in y direction

      integer, intent(out) :: di
                       ! Differential index to istr

      integer, intent(out) :: dj
                       ! Differential index to jstr

      integer, intent(out) :: istrb
                       ! Minimum do loops index in x direction
                       ! of lateral boundary

      integer, intent(out) :: iendb
                       ! Maximum do loops index in x direction
                       ! of lateral boundary

      integer, intent(out) :: jstrb
                       ! Minimum do loops index in y direction
                       ! of lateral boundary

      integer, intent(out) :: jendb
                       ! Maximum do loops index in y direction
                       ! of lateral boundary

      integer, intent(out) :: dib
                       ! Differential index to istrb

      integer, intent(out) :: djb
                       ! Differential index to jstrb

!-----7--------------------------------------------------------------7--

! Calculate the minimum and maximum and differential do loops indices to
! reposition.

      if(imin.le.imin_rst) then

        istr=imin_rst-ies-2
        istrb=imin_rst-ies-1

        di=-istr
        dib=1-istrb

      else

        istr=2
        istrb=2

        di=imin-ies_rst-istr
        dib=imin-ies_rst-istrb

      end if

      if(imax.ge.imax_rst) then

        iend=imax_rst-ies+2
        iendb=imax_rst-ies+1

      else

        iend=ni-2
        iendb=ni-2

      end if

      if(jmin.le.jmin_rst) then

        jstr=jmin_rst-jes-2
        jstrb=jmin_rst-jes-1

        dj=-jstr
        djb=1-jstrb

      else

        jstr=2
        jstrb=2

        dj=jmin-jes_rst-jstr
        djb=jmin-jes_rst-jstrb

      end if

      if(jmax.ge.jmax_rst) then

        jend=jmax_rst-jes+2
        jendb=jmax_rst-jes+1

      else

        jend=nj-2
        jendb=nj-2

      end if

! -----

      end subroutine s_rstindx

!-----7--------------------------------------------------------------7--

      end module m_rstindx
