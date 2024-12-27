!***********************************************************************
      module m_allocflt
!***********************************************************************

!     Author      : Satoki Tsujino
!     Date        : 2017/05/30
!     Modification: 2017/06/11

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the array for Asselin Filter

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commpi
      use m_comflt
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_setcst3d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: allocflt, s_allocflt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface allocflt

        module procedure s_allocflt

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
      subroutine s_allocflt(fpdmpvar,ni,nj,nk)
!***********************************************************************

      use m_comflt

! Input Valriables
      integer, intent(in) :: fpdmpvar
                       ! Formal parameter of unique index of dmpvar

      integer, intent(in) :: ni
                       ! Model dimension in sub x direction

      integer, intent(in) :: nj
                       ! Model dimension in sub y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Internal shared variables
      character(len=108) dmpvar
                       ! Control flag of dump variables

! Internal Valriables
      integer stat     ! Runtime status

      integer cstat    ! Runtime status at current allocate statement

!-----7--------------------------------------------------------------7--

      call getcname(fpdmpvar,dmpvar)

      stat=0

      if(dmpvar(1:1).eq.'+')then

        allocate(dmpfltu(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(netdmpu(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

! Initialize allocated arrays

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,dmpfltu)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,netdmpu)

      end if

      if(dmpvar(2:2).eq.'+')then

        allocate(dmpfltv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(netdmpv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

! Initialize allocated arrays

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,dmpfltv)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,netdmpv)

      end if

      if(dmpvar(3:3).eq.'+')then

        allocate(dmpfltw(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(netdmpw(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

! Initialize allocated arrays

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,dmpfltw)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,netdmpw)

      end if

      if(dmpvar(4:4).eq.'+')then

        allocate(dmpfltp(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(netdmpp(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

! Initialize allocated arrays

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,dmpfltp)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,netdmpp)

      end if

      if(dmpvar(5:5).eq.'+')then

        allocate(dmpfltpt(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(netdmppt(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

! Initialize allocated arrays

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,dmpfltpt)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,netdmppt)

      end if

      if(dmpvar(6:6).eq.'+')then

        allocate(dmpfltqv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(netdmpqv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

! Initialize allocated arrays

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,dmpfltqv)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,netdmpqv)

      end if

! -----

!! -----

      end subroutine s_allocflt

!-----7--------------------------------------------------------------7--

      end module m_allocflt
