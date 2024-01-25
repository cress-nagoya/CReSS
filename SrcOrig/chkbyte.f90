!***********************************************************************
      module m_chkbyte
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2008/04/17
!     Modification: 2008/05/02, 2008/08/25, 2008/10/10, 2009/01/05,
!                   2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the endian of running machine.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comdmp
      use m_comfile
      use m_comkind
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getunit
      use m_putunit

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chkbyte, s_chkbyte

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chkbyte

        module procedure s_chkbyte

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
      subroutine s_chkbyte
!***********************************************************************

! Internal shared variables

      integer(kind=i8) i8scl
                       ! Optional double precision integer variable

      integer i4scl0   ! Optional single precision integer variable
      integer i4scl1   ! Optional single precision integer variable

      integer iobyte   ! Unit number of endian checking file

      integer wli8     ! Word length of direct access file
                       ! for double precision integer variables

      integer stat     ! Runtime status

!-----7--------------------------------------------------------------7--

!! Check the endian of running machine.

! Initialize the double precision integer variable.

      if(mype.eq.root) then
        i8scl=1_i8
      else
        i8scl=0_i8
      end if

! -----

! Get the unit number.

      if(mype.eq.root) then

        call getunit(iobyte)

      end if

! -----

! Open the byte order checking file.

      if(mype.eq.root) then

        inquire(iolength=wli8) i8scl

        open(iobyte,iostat=stat,err=100,                                &
     &       file=curdir(1:2)//bytefl(1:16),                            &
     &       status='new',access='direct',form='unformatted',           &
     &       recl=wli8,action='readwrite')

      else

        stat=0

      end if

  100 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('chkbyte ',7,'warn',10,'              ',14,      &
     &                 iobyte,stat)

        end if

        call cpondpe

        go to 110

      end if

! -----

! Read in the test data to the temporary binary file.

      if(mype.eq.root) then

        write(iobyte,rec=1,iostat=stat,err=120) i8scl

      else

        stat=0

      end if

  120 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('chkbyte ',7,'warn',10,'              ',14,      &
     &                 iobyte,stat)

        end if

        call cpondpe

        go to 110

      end if

! -----

! Read out the test data from the temporary binary file.

      if(mype.eq.root) then

        read(iobyte,rec=1,iostat=stat,err=130) i4scl0,i4scl1

      else

        i4scl0=0
        i4scl1=0

        stat=0

      end if

  130 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('chkbyte ',7,'warn',10,'              ',14,      &
     &                 iobyte,stat)

        end if

        call cpondpe

        go to 110

      end if

! -----

! Close the byte order checking file.

      if(mype.eq.root) then

        close(iobyte,iostat=stat,err=140,status='delete')

      else

        stat=0

      end if

  140 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('chkbyte ',7,'warn',10,'              ',14,      &
     &                 iobyte,stat)

        end if

        call cpondpe

        go to 110

      end if

! -----

! Check the endian of running machine.

      if(mype.eq.root) then

        if(i4scl0.eq.0) then

          write(border(1:7),'(a7)') 'big    '

        else if(i4scl1.eq.0) then

          write(border(1:7),'(a7)') 'little '

        end if

      end if

! -----

! Return the unit number.

  110 if(mype.eq.root) then

        call putunit(iobyte)

      end if

! -----

!! -----

      end subroutine s_chkbyte

!-----7--------------------------------------------------------------7--

      end module m_chkbyte
