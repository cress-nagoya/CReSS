!***********************************************************************
      module m_getzph
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/20
!     Modification: 1999/06/07, 1999/06/28, 1999/07/05, 1999/08/03,
!                   1999/08/09, 1999/09/01, 1999/09/30, 1999/10/07,
!                   1999/10/12, 1999/11/01, 2000/01/17, 2000/04/18,
!                   2000/12/18, 2001/01/15, 2001/03/13, 2001/04/15,
!                   2001/05/29, 2001/11/20, 2002/04/02, 2002/06/18,
!                   2002/07/15, 2002/09/09, 2002/10/31, 2003/04/30,
!                   2003/05/19, 2003/09/01, 2003/11/28, 2004/08/20,
!                   2005/02/10, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the z physical coordinates.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc2d
      use m_comindx
      use m_getiname
      use m_gettrn
      use m_getxy
      use m_phycood

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: getzph, s_getzph

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getzph

        module procedure s_getzph

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getzph(fptrnopt,fpexbopt,ni,nj,nk,z,zph,ht,tmp1,tmp2)
!***********************************************************************

! Input variables

      integer, intent(in) :: fptrnopt
                       ! Formal parameter of unique index of trnopt

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: z(1:nk)
                       ! zeta coordinates

! Output variable

      real, intent(out) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

! Internal shared variables

      integer trnopt   ! Option for terrain height setting
      integer exbopt   ! Option for external boundary forcing

      real, intent(inout) :: ht(0:ni+1,0:nj+1)
                       ! Terrain height

      real, intent(inout) :: tmp1(1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(1:nk)
                       ! Temporary array

! Remark

!     zph: This variable is also temporary.

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fptrnopt,trnopt)
      call getiname(fpexbopt,exbopt)

! -----

!! Set the terrain height.

! Get the x and the y coordinates at the model grid points.

      call s_getxy(iddx,iddy,'xx',0,ni+1,0,nj+1,zph(0,0,1),zph(0,0,2))

! -----

! Get the terrain height.

      if(trnopt.eq.2.and.mod(exbopt,10).eq.2) then

        call s_gettrn(idtrnopt,idzsfc,idmnthgh,idmntwx,idmntwy,         &
     &                idmntcx,idmntcy,'terrain.damp',12,1,ni,nj,        &
     &                zph(0,0,1),zph(0,0,2),ht)

      else

        call s_gettrn(idtrnopt,idzsfc,idmnthgh,idmntwx,idmntwy,         &
     &                idmntcx,idmntcy,'terrain     ',7,1,ni,nj,         &
     &                zph(0,0,1),zph(0,0,2),ht)

      end if

      if(exbopt.eq.0) then

        call bc2d(idwbc,idebc,idsbc,idnbc,ni,nj,ht)

      end if

! -----

!! -----

! Calculate the z physical coordinates.

      call phycood(idsthopt,idzsfc,idzflat,ni,nj,nk,z,ht,zph,tmp1,tmp2)

! -----

      end subroutine s_getzph

!-----7--------------------------------------------------------------7--

      end module m_getzph
