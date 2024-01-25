!***********************************************************************
      module m_outcheck
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/04/11
!     Modification: 2007/06/18, 2007/07/30, 2007/08/24, 2008/01/11,
!                   2008/04/17, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2008/10/10, 2009/01/30, 2009/02/27, 2011/09/22,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the data to the dumped data checking file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chkstd
      use m_comdmp
      use m_comkind
      use m_commpi
      use m_comname
      use m_cpondpe
      use m_destroy
      use m_getiname
      use m_outstd03

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outcheck, s_outcheck

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outcheck

        module procedure s_outcheck

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
      subroutine s_outcheck(fpdmpmon,ctime)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdmpmon
                       ! Formal parameter of unique index of dmpmon

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

! Internal shared variables

      integer dmpmon   ! Option for monitor variables output

      integer in       ! Namelist table index

      integer ncend    ! Number of output character for endian option

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

!-----7--------------------------------------------------------------7--

! Read in the data to the dumped data checking file when the current
! forecast time reaches marked time.

      if(fdmp(1:3).eq.'act') then

        if(mype.eq.root) then

          write(io3c,'(a)',iostat=stat,err=100)                         &
     &         (cname(in),in=1,ncn)

          write(io3c,*,iostat=stat,err=100)                             &
     &         (iname(in),in=1,nin)

          write(io3c,*,iostat=stat,err=100)                             &
     &         (rname(in),in=1,nrn)

          if(border(1:7).eq.'unknown') then

            ncend=16

            write(io3c,'(a30,a30,a15,i3)',iostat=stat,err=100)          &
     &                 'options template              ',                &
     &                 '                              ',                &
     &                 '               ',ncend

          else if(border(1:3).eq.'big') then

            ncend=27

            write(io3c,'(a30,a30,a15,i3)',iostat=stat,err=100)          &
     &                 'options template big_endian   ',                &
     &                 '                              ',                &
     &                 '               ',ncend

          else if(border(1:6).eq.'little') then

            ncend=30

            write(io3c,'(a30,a30,a15,i3)',iostat=stat,err=100)          &
     &                 'options template little_endian',                &
     &                 '                              ',                &
     &                 '               ',ncend

          end if

        else

          stat=0

        end if

  100   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('outcheck',8,'cont',3,'              ',14,     &
     &                   io3c,stat)

          end if

          call cpondpe

          call destroy('outcheck',8,'stop',1001,'              ',14,    &
     &                 101,stat)

        end if

        if(mype.eq.root) then

          call outstd03('outcheck',8,fl3c,108,io3c,4,1,ctime)

        end if

        broot=stat-1

        call chkstd(broot)

      end if

! -----

!! Read in the data to the dumped data checking file for monitor
!! variables when the current forecast time reaches marked time.

      if(fmon(1:3).eq.'act') then

! Get the required namelist variable.

        call getiname(fpdmpmon,dmpmon)

! -----

! Read in the data to the dumped data checking file for monitor
! variables.

        if(dmpmon.eq.1) then

          if(mype.eq.root) then

            write(io2c,'(a)',iostat=stat,err=110)                       &
     &           (cname(in),in=1,ncn)

            write(io2c,*,iostat=stat,err=110)                           &
     &           (iname(in),in=1,nin)

            write(io2c,*,iostat=stat,err=110)                           &
     &           (rname(in),in=1,nrn)

            if(border(1:7).eq.'unknown') then

              ncend=16

              write(io2c,'(a30,a30,a15,i3)',iostat=stat,err=110)        &
     &                   'options template              ',              &
     &                   '                              ',              &
     &                   '               ',ncend

            else if(border(1:3).eq.'big') then

              ncend=27

              write(io2c,'(a30,a30,a15,i3)',iostat=stat,err=110)        &
     &                   'options template big_endian   ',              &
     &                   '                              ',              &
     &                   '               ',ncend

            else if(border(1:6).eq.'little') then

              ncend=30

              write(io2c,'(a30,a30,a15,i3)',iostat=stat,err=110)        &
     &                   'options template little_endian',              &
     &                   '                              ',              &
     &                   '               ',ncend

            end if

          else

            stat=0

          end if

  110     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('outcheck',8,'cont',3,'              ',14,   &
     &                     io2c,stat)

            end if

            call cpondpe

            call destroy('outcheck',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.root) then

            call outstd03('outcheck',8,fl2c,108,io2c,4,1,ctime)

          end if

          broot=stat-1

          call chkstd(broot)

        end if

! -----

      end if

!! -----

      end subroutine s_outcheck

!-----7--------------------------------------------------------------7--

      end module m_outcheck
