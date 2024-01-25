!***********************************************************************
      module m_allocasl
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2011/11/10

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the array for asldata.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comasl
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_setcst1d
      use m_setcst2d
      use m_setcst3d
      use m_setcst4d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: allocasl, s_allocasl

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allocasl

        module procedure s_allocasl

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
      subroutine s_allocasl(ni,nj,nk,nqa,nid_asl,njd_asl,nkd_asl,km_asl)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

      integer, intent(in) :: nid_asl
                       ! Aerosol data dimension in x direction

      integer, intent(in) :: njd_asl
                       ! Aerosol data dimension in y direction

      integer, intent(in) :: nkd_asl
                       ! Aerosol data dimension in z direction

      integer, intent(in) :: km_asl
                       ! Dimension of max(nk, nkd_asl)

! Internal shared variables

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

!-----7--------------------------------------------------------------7--

!! Allocate the array for asldata.

! Perform allocate.

      stat=0

      allocate(z(1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(zph(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0)),stat=cstat)

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

      allocate(zdat(1:nid_asl,1:njd_asl,1:nkd_asl),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qadat(1:nid_asl,1:njd_asl,1:km_asl,1:nqa(0)),stat=cstat)

      stat=stat+abs(cstat)

      allocate(dtmp1(1:nid_asl,1:njd_asl,1:nk),stat=cstat)

      stat=stat+abs(cstat)

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('allocasl',8,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('allocasl',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

! Fill in all array for the program asldata with 0.

      call setcst1d(1,nk,0.e0,z)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,zph)

      call setcst4d(0,ni+1,0,nj+1,1,nk,1,nqa(0),0.e0,qasl)

      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp1)
      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp2)
      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp3)
      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp4)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp5)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp6)

      call setcst3d(1,nid_asl,1,njd_asl,1,nkd_asl,0.e0,zdat)

      call setcst4d(1,nid_asl,1,njd_asl,1,km_asl,1,nqa(0),0.e0,qadat)

      call setcst3d(1,nid_asl,1,njd_asl,1,nk,0.e0,dtmp1)

! -----

      end subroutine s_allocasl

!-----7--------------------------------------------------------------7--

      end module m_allocasl
