!***********************************************************************
      module m_allocsfc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/07/15, 2003/08/08, 2003/11/05, 2003/12/12,
!                   2004/08/01, 2004/09/01, 2005/01/14, 2005/02/10,
!                   2006/01/10, 2007/01/20, 2007/01/31, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/01/05, 2009/02/27,
!                   2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the array for surface.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commpi
      use m_comsfc
      use m_cpondpe
      use m_destroy
      use m_setcst2d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: allocsfc, s_allocsfc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allocsfc

        module procedure s_allocsfc

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
      subroutine s_allocsfc(ni,nj,nid_lnd,njd_lnd,nid_sst,njd_sst,      &
     &                      nid_ice,njd_ice)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nid_lnd
                       ! Land use data dimension in x direction

      integer, intent(in) :: njd_lnd
                       ! Land use data dimension in y direction

      integer, intent(in) :: nid_sst
                       ! Sea surface temperature data dimension
                       ! in x direction

      integer, intent(in) :: njd_sst
                       ! Sea surface temperature data dimension
                       ! in y direction

      integer, intent(in) :: nid_ice
                       ! Sea ice distribution data dimension
                       ! in x direction

      integer, intent(in) :: njd_ice
                       ! Sea ice distribution data dimension
                       ! in y direction

! Internal shared variables

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer id       ! Data array index in x direction
      integer jd       ! Data array index in y direction

!-----7--------------------------------------------------------------7--

!! Allocate the array for surface.

! Perform allocate.

      stat=0

      allocate(land(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(landat(1:nid_lnd,1:njd_lnd),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ri(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(rj(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(albe(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(beta(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(z0m(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(z0h(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(cap(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(nuu(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(sst(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(kai(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(sstdat(1:nid_sst,1:njd_sst),stat=cstat)

      stat=stat+abs(cstat)

      allocate(icedat(1:nid_ice,1:njd_ice),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp1(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp2(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp3(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ltmp1(1:nid_lnd,1:njd_lnd),stat=cstat)

      stat=stat+abs(cstat)

      allocate(stmp1(1:nid_sst,1:njd_sst),stat=cstat)

      stat=stat+abs(cstat)

      allocate(itmp1(1:nid_ice,1:njd_ice),stat=cstat)

      stat=stat+abs(cstat)

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('allocsfc',8,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('allocsfc',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

!! Fill in all array for the program terrain with 0.

! For the integer variables.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i,j)

      do j=0,nj+1
      do i=0,ni+1
        land(i,j)=0
      end do
      end do

!$omp end do

!$omp do schedule(runtime) private(id,jd)

      do jd=1,njd_lnd
      do id=1,nid_lnd
        landat(id,jd)=0
      end do
      end do

!$omp end do

!$omp end parallel

! -----

! For the real variables.

      call setcst2d(0,ni+1,0,nj+1,0.e0,ri)
      call setcst2d(0,ni+1,0,nj+1,0.e0,rj)

      call setcst2d(0,ni+1,0,nj+1,0.e0,albe)

      call setcst2d(0,ni+1,0,nj+1,0.e0,beta)

      call setcst2d(0,ni+1,0,nj+1,0.e0,z0m)
      call setcst2d(0,ni+1,0,nj+1,0.e0,z0h)

      call setcst2d(0,ni+1,0,nj+1,0.e0,cap)
      call setcst2d(0,ni+1,0,nj+1,0.e0,nuu)

      call setcst2d(0,ni+1,0,nj+1,0.e0,sst)

      call setcst2d(0,ni+1,0,nj+1,0.e0,kai)

      call setcst2d(1,nid_sst,1,njd_sst,0.e0,sstdat)

      call setcst2d(1,nid_ice,1,njd_ice,0.e0,icedat)

      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp1)
      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp2)
      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp3)

      call setcst2d(1,nid_lnd,1,njd_lnd,0.e0,ltmp1)

      call setcst2d(1,nid_sst,1,njd_sst,0.e0,stmp1)

      call setcst2d(1,nid_ice,1,njd_ice,0.e0,itmp1)

! -----

!! -----

      end subroutine s_allocsfc

!-----7--------------------------------------------------------------7--

      end module m_allocsfc
