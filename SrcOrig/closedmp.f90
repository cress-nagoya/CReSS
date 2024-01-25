!***********************************************************************
      module m_closedmp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/25, 1999/06/14, 1999/06/21, 1999/08/23,
!                   1999/11/19, 2000/01/05, 2000/01/17, 2000/04/18,
!                   2001/04/15, 2001/06/18, 2001/07/15, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2004/05/31,
!                   2004/06/10, 2005/02/10, 2006/09/21, 2007/01/20,
!                   2007/06/18, 2007/07/30, 2008/01/11, 2008/04/17,
!                   2008/05/02, 2008/07/01, 2008/08/25, 2009/02/27,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     close the dumped file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chkstd
      use m_comdmp
      use m_comkind
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getiname
      use m_outstd03
      use m_outstd08
      use m_putunit

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: closedmp, s_closedmp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface closedmp

        module procedure s_closedmp

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
      subroutine s_closedmp(fpdmpmon,ctime)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdmpmon
                       ! Formal parameter of unique index of dmpmon

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

! Internal shared variables

      integer dmpmon   ! Option for monitor variables output

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpdmpmon,dmpmon)

! -----

!! Close the dumped file.

      if(fdmp(1:3).eq.'act'.or.fmon(1:3).eq.'act') then

! Read in the messages to the standard i/o.

        if(dmpmon.eq.0) then

          if(cnt3d.eq.0) then

            if(mype.eq.root) then

              call outstd08(1)

            end if

            call cpondpe

          end if

        else

          if(fdmp(1:3).eq.'act'.and.fmon(1:3).eq.'act') then

            if(cnt3d.eq.0.and.cnt2d.eq.0) then

              if(mype.eq.root) then

                call outstd08(1)

              end if

              call cpondpe

            else if(cnt3d.eq.0.and.cnt2d.gt.0) then

              if(mype.eq.root) then

                call outstd08(2)

              end if

              call cpondpe

            else if(cnt3d.gt.0.and.cnt2d.eq.0) then

              if(mype.eq.root) then

                call outstd08(3)

              end if

              call cpondpe

            end if

          else if(fdmp(1:3).eq.'act'.and.fmon(1:3).eq.'off') then

            if(cnt3d.eq.0) then

              if(mype.eq.root) then

                call outstd08(2)

              end if

              call cpondpe

            end if

          else if(fdmp(1:3).eq.'off'.and.fmon(1:3).eq.'act') then

            if(cnt2d.eq.0) then

              if(mype.eq.root) then

                call outstd08(3)

              end if

              call cpondpe

            end if

          end if

        end if

! -----

! Close the dumped data checking file.

        if(fdmp(1:3).eq.'act') then

          if(mype.eq.root) then

            if(cnt3d.eq.0) then

              close(io3c,iostat=stat,err=100,status='delete')

            else

              close(io3c,iostat=stat,err=100,status='keep')

            end if

          else

            stat=0

          end if

  100     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('closedmp',8,'cont',2,'              ',14,   &
     &                     io3c,stat)

            end if

            call cpondpe

            call destroy('closedmp',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.root) then

            call outstd03('closedmp',8,fl3c,108,io3c,2,1,ctime)

          end if

          broot=stat-1

          call chkstd(broot)

        end if

! -----

! Close the dumped data checking file for monitor variables.

        if(fmon(1:3).eq.'act') then

          if(dmpmon.eq.1) then

            if(mype.eq.root) then

              if(cnt2d.eq.0) then

                close(io2c,iostat=stat,err=110,status='delete')

              else

                close(io2c,iostat=stat,err=110,status='keep')

              end if

            else

              stat=0

            end if

  110       call chkerr(stat)

            if(stat.lt.0) then

              if(mype.eq.-stat-1) then

                call destroy('closedmp',8,'cont',2,'              ',14, &
     &                       io2c,stat)

              end if

              call cpondpe

              call destroy('closedmp',8,'stop',1001,'              ',14,&
     &                     101,stat)

            end if

            if(mype.eq.root) then

              call outstd03('closedmp',8,fl2c,108,io2c,2,1,ctime)

            end if

            broot=stat-1

            call chkstd(broot)

          end if

        end if

! -----

! Close the dumped file.

        if(fdmp(1:3).eq.'act') then

          if(cnt3d.eq.0) then

            close(io3d,iostat=stat,err=120,status='delete')

          else

            close(io3d,iostat=stat,err=120,status='keep')

          end if

  120     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('closedmp',8,'cont',2,'              ',14,   &
     &                     io3d,stat)

            end if

            call cpondpe

            call destroy('closedmp',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.stat-1) then

            call outstd03('closedmp',8,fl3d,108,io3d,2,1,ctime)

          end if

          broot=stat-1

          call chkstd(broot)

        end if

! -----

! Close the dumped file for monitor variables.

        if(fmon(1:3).eq.'act') then

          if(dmpmon.eq.1) then

            if(cnt2d.eq.0) then

              close(io2d,iostat=stat,err=130,status='delete')

            else

              close(io2d,iostat=stat,err=130,status='keep')

            end if

  130       call chkerr(stat)

            if(stat.lt.0) then

              if(mype.eq.-stat-1) then

                call destroy('closedmp',8,'cont',2,'              ',14, &
     &                       io2d,stat)

              end if

              call cpondpe

              call destroy('closedmp',8,'stop',1001,'              ',14,&
     &                     101,stat)

            end if

            if(mype.eq.stat-1) then

              call outstd03('closedmp',8,fl2d,108,io2d,2,1,ctime)

            end if

            broot=stat-1

            call chkstd(broot)

          end if

        end if

! -----

! Return the unit numbers.

        if(fmon(1:3).eq.'act') then

          if(dmpmon.eq.1) then

            call putunit(io2d)

          end if

        end if

        if(fdmp(1:3).eq.'act') then

          call putunit(io3d)

        end if

        if(mype.eq.root) then

          if(fmon(1:3).eq.'act') then

            if(dmpmon.eq.1) then

              call putunit(io2c)

            end if

          end if

          if(fdmp(1:3).eq.'act') then

            call putunit(io3c)

          end if

        end if

! -----

! Set the control flag fdmp and fmon to close the dumped file to the
! module m_comdmp.

        if(fdmp(1:3).eq.'act') then

          write(fdmp(1:3),'(a3)') 'off'

        end if

        if(fmon(1:3).eq.'act') then

          write(fmon(1:3),'(a3)') 'off'

        end if

! -----

      end if

!! -----

      end subroutine s_closedmp

!-----7--------------------------------------------------------------7--

      end module m_closedmp
