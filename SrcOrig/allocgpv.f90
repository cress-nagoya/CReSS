!***********************************************************************
      module m_allocgpv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/05/19
!     Modification: 2003/08/08, 2003/11/05, 2003/12/12, 2004/01/09,
!                   2004/07/01, 2005/02/10, 2006/01/10, 2006/02/03,
!                   2007/01/20, 2008/05/02, 2008/08/19, 2008/08/25,
!                   2008/12/11, 2009/01/05, 2009/02/27, 2011/09/22

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the array for gridata.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comgpv
      use m_commpi
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

      public :: allocgpv, s_allocgpv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allocgpv

        module procedure s_allocgpv

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
      subroutine s_allocgpv(ni,nj,nk,nid_gpv,njd_gpv,nkd_gpv,km_gpv)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nid_gpv
                       ! GPV data dimension in x direction

      integer, intent(in) :: njd_gpv
                       ! GPV data dimension in y direction

      integer, intent(in) :: nkd_gpv
                       ! GPV data dimension in z direction

      integer, intent(in) :: km_gpv
                       ! Dimension of max(nk, nkd_gpv)

! Internal shared variables

      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

!-----7--------------------------------------------------------------7--

!! Allocate the array for gridata.

! Perform allocate.

      stat=0

      allocate(land(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(z(1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(zph(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ubr(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(vbr(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(pbr(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ptbr(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qvbr(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(u(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(v(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(w(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(pp(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ptp(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qv(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qc(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qr(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qi(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qs(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qg(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qh(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(z1d(0:4*nk-3),stat=cstat)

      stat=stat+abs(cstat)

      allocate(u1d(0:4*nk-3),stat=cstat)

      stat=stat+abs(cstat)

      allocate(v1d(0:4*nk-3),stat=cstat)

      stat=stat+abs(cstat)

      allocate(p1d(0:4*nk-3),stat=cstat)

      stat=stat+abs(cstat)

      allocate(pt1d(0:4*nk-3),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qv1d(0:4*nk-3),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp1(1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp2(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp3(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp4(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp5(0:ni+1,0:nj+1),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp6(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(tmp7(0:ni+1,0:nj+1,1:nk),stat=cstat)

      stat=stat+abs(cstat)

      allocate(londat(1:nid_gpv,1:njd_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(htdat(1:nid_gpv,1:njd_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(zdat(1:nid_gpv,1:njd_gpv,1:nkd_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(udat(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(vdat(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(wdat(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(pdat(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ptdat(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qvdat(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qcdat(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qrdat(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qidat(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qsdat(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qgdat(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(qhdat(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(pbdat(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(ptbdat(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(dtmp1(1:nid_gpv,1:njd_gpv),stat=cstat)

      stat=stat+abs(cstat)

      allocate(dtmp2(1:nid_gpv,1:njd_gpv,1:km_gpv),stat=cstat)

      stat=stat+abs(cstat)

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('allocgpv',8,'cont',5,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('allocgpv',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! -----

! Fill in all array for the program gridata with 0.

      call setcst1d(1,nk,0.e0,z)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,zph)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ubr)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vbr)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,pbr)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ptbr)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qvbr)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,u)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,v)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,w)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,pp)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ptp)

      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qv)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qc)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qr)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qi)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qs)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qg)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,qh)

      call setcst1d(0,4*nk-3,0.e0,z1d)
      call setcst1d(0,4*nk-3,0.e0,u1d)
      call setcst1d(0,4*nk-3,0.e0,v1d)
      call setcst1d(0,4*nk-3,0.e0,p1d)
      call setcst1d(0,4*nk-3,0.e0,pt1d)
      call setcst1d(0,4*nk-3,0.e0,qv1d)

      call setcst1d(1,nk,0.e0,tmp1)
      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp2)
      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp3)
      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp4)
      call setcst2d(0,ni+1,0,nj+1,0.e0,tmp5)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp6)
      call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,tmp7)

      call setcst2d(1,nid_gpv,1,njd_gpv,0.e0,londat)

      call setcst2d(1,nid_gpv,1,njd_gpv,0.e0,htdat)

      call setcst3d(1,nid_gpv,1,njd_gpv,1,nkd_gpv,0.e0,zdat)

      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,udat)
      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,vdat)
      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,wdat)

      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,pdat)
      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,ptdat)

      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,qvdat)
      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,qcdat)
      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,qrdat)
      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,qidat)
      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,qsdat)
      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,qgdat)
      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,qhdat)

      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,pbdat)
      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,ptbdat)

      call setcst2d(1,nid_gpv,1,njd_gpv,0.e0,dtmp1)
      call setcst3d(1,nid_gpv,1,njd_gpv,1,km_gpv,0.e0,dtmp2)

! -----

      end subroutine s_allocgpv

!-----7--------------------------------------------------------------7--

      end module m_allocgpv
