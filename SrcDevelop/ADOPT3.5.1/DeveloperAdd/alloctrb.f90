!***********************************************************************
      module m_alloctrb
!***********************************************************************

!     Author      : Satoki Tsujino
!     Date        : 2014/01/28
!     Modification: 2017/06/02   Adding turbulence terms for pt and qv
!                   2017/06/18   Adding numerical diffusion terms
!                   2018/01/05   Adding eddy dif and vis coefficients

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     allocate the array for turbulence physics.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comtub
      use m_commpi
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

      public :: alloctrb, s_alloctrb

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface alloctrb

        module procedure s_alloctrb

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
      subroutine s_alloctrb(fpdmpvar,ni,nj,nk)
!***********************************************************************

! Input Valriables
      integer, intent(in) :: fpdmpvar
                       ! Formal parameter of unique index of dmpvar

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

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

      if(dmpvar(16:16).eq.'o'.or.dmpvar(16:16).eq.'+')then

        stat=0

        allocate(turbu(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(turbv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(turbpt(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(turbqv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,turbu)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,turbv)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,turbpt)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,turbqv)

        if(dmpvar(16:16).eq.'+')then

          allocate(numdu(0:ni+1,0:nj+1,1:nk),stat=cstat)

          stat=stat+abs(cstat)

          allocate(numdv(0:ni+1,0:nj+1,1:nk),stat=cstat)

          stat=stat+abs(cstat)

          allocate(numdw(0:ni+1,0:nj+1,1:nk),stat=cstat)

          stat=stat+abs(cstat)

          allocate(numdpt(0:ni+1,0:nj+1,1:nk),stat=cstat)

          stat=stat+abs(cstat)

          allocate(numdqv(0:ni+1,0:nj+1,1:nk),stat=cstat)

          stat=stat+abs(cstat)

          call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,numdu)
          call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,numdv)
          call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,numdw)
          call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,numdpt)
          call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,numdqv)

        end if

! If error occured, call the procedure destroy.

        call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('alloctrb',8,'cont',5,'              ',14,101, &
     &                   stat)

          end if

          call cpondpe

          call destroy('alloctrb',8,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

      else if(dmpvar(16:16).eq.'v')then

        stat=0

        allocate(difh(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(difv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(vish(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        allocate(visv(0:ni+1,0:nj+1,1:nk),stat=cstat)

        stat=stat+abs(cstat)

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,difh)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,difv)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,vish)
        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,visv)

! If error occured, call the procedure destroy.

        call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('alloctrb',8,'cont',5,'              ',14,101, &
     &                   stat)

          end if

          call cpondpe

          call destroy('alloctrb',8,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

      end if

! -----

!! -----

      end subroutine s_alloctrb

!-----7--------------------------------------------------------------7--

      end module m_alloctrb
