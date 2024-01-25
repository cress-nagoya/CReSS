!***********************************************************************
      module m_alloctrn
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/08/08, 2003/11/05, 2005/02/10, 2006/01/10,
!                   2007/01/20, 2008/05/02, 2008/08/25, 2009/01/05,
!                   2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the array for terrain.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commpi
      use m_comtrn
      use m_cpondpe
      use m_destroy
      use m_setcst2d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: alloctrn, s_alloctrn

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface alloctrn

        module procedure s_alloctrn

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_alloctrn(ni,nj,nid_trn,njd_trn)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nid_trn
                       ! Terrain data dimension in x direction

      integer, intent(in) :: njd_trn
                       ! Terrain data dimension in y direction

! Internal shared variables

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

!-----7--------------------------------------------------------------7--

!! Allocate the array for terrain.

! Perform allocate.

      stat=0

      allocate(ri(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(rj(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ht(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp1(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp2(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp3(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(htdat(1:nid_trn,1:njd_trn),stat=cstat)

      stat=stat+abs(cstat)

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('alloctrn',8,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('alloctrn',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

! Fill in all array for the program terrain with 0.

      call setcst2d(0,ni+1,0,nj+1,0.e0,ri)
      call setcst2d(0,ni+1,0,nj+1,0.e0,rj)

      call setcst2d(0,ni+1,0,nj+1,0.e0,ht)

      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp1)
      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp2)
      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp3)

      call setcst2d(1,nid_trn,1,njd_trn,0.e0,htdat)

! -----

      end subroutine s_alloctrn

!-----7--------------------------------------------------------------7--

      end module m_alloctrn
