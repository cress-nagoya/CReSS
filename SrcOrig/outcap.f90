!***********************************************************************
      module m_outcap
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2008/01/11
!     Modification: 2008/04/17, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/02/27, 2011/09/22

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the messages to the dumped data checking file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commpi
      use m_cpondpe
      use m_destroy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outcap, s_outcap

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outcap

        module procedure s_outcap

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
      subroutine s_outcap(fproc,vname,vcap,ncvc,iochk,nk)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: fproc
                       ! Control flag of processing type

      character(len=6), intent(in) :: vname
                       ! Optional variable name

      character(len=60), intent(in) :: vcap
                       ! Caption for dumped variable

      integer, intent(in) :: ncvc
                       ! Number of character of vcap

      integer, intent(in) :: iochk
                       ! Unit number of dumped data checking file

      integer, intent(in) :: nk
                       ! Model dimension in z direction

! Internal shared variable

      integer stat     ! Runtime status

!-----7--------------------------------------------------------------7--

! Read in the messages to the dumped data checking file.

      if(mype.eq.root) then

        write(iochk,'(a6,i6,a3,a60,i3)',iostat=stat,err=100)            &
     &                vname(1:6),nk-3,' 0 ',vcap(1:60),ncvc+15

      else

        stat=0

      end if

  100 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          if(fproc(1:3).eq.'geo') then

            call destroy('outgeo  ',6,'cont',3,'              ',14,     &
     &                   iochk,stat)

          else if(fproc(1:3).eq.'dmp') then

            if(nk.eq.3) then

              call destroy('outdmp2d',8,'cont',3,'              ',14,   &
     &                     iochk,stat)

            else

              call destroy('outdmp3d',8,'cont',3,'              ',14,   &
     &                     iochk,stat)

            end if

          end if

        end if

        call cpondpe

        call destroy('outcap  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

      end subroutine s_outcap

!-----7--------------------------------------------------------------7--

      end module m_outcap
