!***********************************************************************
      module m_allocrdr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/08/08, 2003/11/05, 2003/12/12, 2005/02/10,
!                   2006/01/10, 2007/01/20, 2008/05/02, 2008/08/25,
!                   2009/01/05, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the array for radata.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commpi
      use m_comrdr
      use m_cpondpe
      use m_destroy
      use m_setcst1d
      use m_setcst2d
      use m_setcst3d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: allocrdr, s_allocrdr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allocrdr

        module procedure s_allocrdr

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
      subroutine s_allocrdr(ni,nj,nk,nid_rdr,njd_rdr,nkd_rdr,km_rdr)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nid_rdr
                       ! Radar data dimension in x direction

      integer, intent(in) :: njd_rdr
                       ! Radar data dimension in y direction

      integer, intent(in) :: nkd_rdr
                       ! Radar data dimension in z direction

      integer, intent(in) :: km_rdr
                       ! Dimension of max(nk, nkd_rdr)

! Internal shared variables

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

!-----7--------------------------------------------------------------7--

!! Allocate the array for radata.

! Perform allocate.

      stat=0

      allocate(z(1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(zph(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(u(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(v(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(w(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qp(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp1(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp2(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp3(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp4(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp5(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp6(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(londat(1:nid_rdr,1:njd_rdr),stat=cstat)

      stat=stat+abs(cstat)

      allocate(zdat(1:nid_rdr,1:njd_rdr,1:nkd_rdr),stat=cstat)

      stat=stat+abs(cstat)

      allocate(udat(1:nid_rdr,1:njd_rdr,1:km_rdr),stat=cstat)

      stat=stat+abs(cstat)

      allocate(vdat(1:nid_rdr,1:njd_rdr,1:km_rdr),stat=cstat)

      stat=stat+abs(cstat)

      allocate(wdat(1:nid_rdr,1:njd_rdr,1:km_rdr),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qpdat(1:nid_rdr,1:njd_rdr,1:km_rdr),stat=cstat)

      stat=stat+abs(cstat)

      allocate(dtmp1(1:nid_rdr,1:njd_rdr,1:nk),stat=cstat)

      stat=stat+abs(cstat)

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('allocrdr',8,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('allocrdr',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

! Fill in all array for the program radata with 0.

      call setcst1d(1,nk,0.e0,z)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,zph)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,u)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,v)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,w)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qp)

      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp1)
      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp2)
      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp3)
      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp4)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp5)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp6)

      call setcst2d(1,nid_rdr,1,njd_rdr,0.e0,londat)

      call setcst3d(1,nid_rdr,1,njd_rdr,1,nkd_rdr,0.e0,zdat)

      call setcst3d(1,nid_rdr,1,njd_rdr,1,km_rdr,0.e0,udat)
      call setcst3d(1,nid_rdr,1,njd_rdr,1,km_rdr,0.e0,vdat)
      call setcst3d(1,nid_rdr,1,njd_rdr,1,km_rdr,0.e0,wdat)

      call setcst3d(1,nid_rdr,1,njd_rdr,1,km_rdr,0.e0,qpdat)

      call setcst3d(1,nid_rdr,1,njd_rdr,1,nk,0.e0,dtmp1)

! -----

      end subroutine s_allocrdr

!-----7--------------------------------------------------------------7--

      end module m_allocrdr
