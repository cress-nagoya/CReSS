!***********************************************************************
      module m_destroy
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/01/25, 1999/04/06, 1999/06/14, 1999/06/21,
!                   1999/08/23, 2000/01/17, 2000/04/18, 2000/07/05,
!                   2000/12/25, 2001/02/13, 2001/04/15, 2001/05/29,
!                   2001/11/20, 2002/06/18, 2002/12/17, 2003/04/30,
!                   2003/05/19, 2003/10/06, 2004/05/31, 2004/08/20,
!                   2005/02/10, 2006/09/21, 2007/01/20, 2007/11/13,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/01/05,
!                   2009/01/30, 2009/02/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the error messages to the standard i/o and stop the
!     program.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comerr
      use m_endmpi

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: destroy, s_destroy

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface destroy

        module procedure s_destroy

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
      subroutine s_destroy(sname,ncsn,fproc,fmsg,nname,ncnn,ionum,stat)
!***********************************************************************

! Input variables

      character(len=8), intent(in) :: sname
                       ! Procedure name occured error

      character(len=4), intent(in) :: fproc
                       ! Control flag of processing type

      character(len=14), intent(in) :: nname
                       ! Erroneous namalist variable

      integer, intent(in) :: ncsn
                       ! Number of character of sname

      integer, intent(in) :: ncnn
                       ! Number of character of nname

      integer, intent(in) :: fmsg
                       ! Control flag of message type

      integer, intent(in) :: ionum
                       ! File unit number occured error

! Input and output variable

      integer, intent(inout) :: stat
                       ! Runtime status

! Internal shared variables

      integer istat    ! Do loops index

      integer cstat    ! Runtime status in do loops

!-----7--------------------------------------------------------------7--

! Read in the error messages to the standard i/o.

      if(fmsg.lt.200) then

        write(6,*)

      end if

      if(fmsg.lt.100) then

        if(fproc(1:4).eq.'warn') then

          write(6,'(a,a,a)') '  warning: procedure, ',sname(1:ncsn),';'

        else

          write(6,'(a,a,a)') '  error: procedure, ',sname(1:ncsn),';'

        end if

      end if

      if(fmsg.eq.1) then

        write(6,'(a,i5.5,a)')                                           &
     &    '    Can not open the file of unit number ',ionum,'.'

      else if(fmsg.eq.2) then

        write(6,'(a,i5.5,a)')                                           &
     &    '    Can not close the file of unit number ',ionum,'.'

      else if(fmsg.eq.3) then

        write(6,'(a,i5.5,a)')                                           &
     &    '    Can not read out the data from the file of unit number ',&
     &    ionum,'.'

      else if(fmsg.eq.4) then

        write(6,'(a,i5.5,a)')                                           &
     &    '    Can not read in the data to the file of unit number ',   &
     &    ionum,'.'

      else if(fmsg.eq.5) then

        write(6,'(a)') '    Can not allocate the array.'

      else if(fmsg.eq.6) then

        write(6,'(a)') '    Can not deallocate the array.'

      else if(fmsg.eq.7) then

        write(6,'(a)')                                                  &
     &    '    Can not go on with this calculation.'

      else if(fmsg.eq.8) then

        write(6,'(a)') '    This program was aborted.'

      else if(fmsg.eq.9) then

        write(6,'(a)')                                                  &
     &    '    The program, solver was aborted by user handling.'

      else if(fmsg.eq.10) then

        write(6,'(a)')                                                  &
     &    '    Can not perform endian checking.'

      else if(fmsg.eq.11) then

        write(6,'(a)')                                                  &
     &    '    Wrong number of processor elements is specified.'

      else if(fmsg.eq.101) then

        write(6,'(a)') '  bye-bye!!'

      else if(fmsg.eq.102) then

        write(6,'(a,a,a)')                                              &
     &    '  Why do you know new option, ',nname(1:ncnn),'?'

      else if(fmsg.eq.201) then

        if(stat.eq.0) then

          stat=stat+1

          write(errlst(1)(1:14),'(a14)') nname(1:14)

          write(6,*)

          write(6,'(a,a,a)') '  error: procedure, ',sname(1:ncsn),';'

          write(6,'(a,a,a)')                                            &
     &            '    The namelist, "',nname(1:ncnn),'" is wrong.'

        else

          cstat=0

          do_istat_1: do istat=1,stat

            if(nname(1:ncnn).eq.errlst(istat)(1:ncnn)) then

              cstat=1

              exit do_istat_1

            end if

          end do do_istat_1

          if(cstat.eq.0) then

            stat=stat+1

            write(errlst(stat)(1:14),'(a14)') nname(1:14)

            write(6,*)

            write(6,'(a,a,a)') '  error: procedure, ',sname(1:ncsn),';'

            write(6,'(a,a,a)')                                          &
     &              '    The namelist, "',nname(1:ncnn),'" is wrong.'

          end if

        end if

      else if(fmsg.eq.202) then

        if(stat.eq.0) then

          stat=stat+1

          write(errlst(1)(1:14),'(a14)') nname(1:14)

          write(6,*)

          write(6,'(a,a,a)') '  error: procedure, ',sname(1:ncsn),';'

          write(6,'(a,a,a)') '    The namelist, "',nname(1:ncnn),       &
     &            '" is different from previous run.'

        else

          cstat=0

          do_istat_2: do istat=1,stat

            if(nname(1:ncnn).eq.errlst(istat)(1:ncnn)) then

              cstat=1

              exit do_istat_2

            end if

          end do do_istat_2

          if(cstat.eq.0) then

            stat=stat+1

            write(errlst(stat)(1:14),'(a14)') nname(1:14)

            write(6,*)

            write(6,'(a,a,a)') '  error: procedure, ',sname(1:ncsn),';'

            write(6,'(a,a,a)') '    The namelist, "',nname(1:ncnn),     &
     &              '" is different from previous run.'

          end if

        end if

      end if

      if(fproc(1:4).ne.'stop') then

        return

      end if

! -----

! Kill the forked processes.

      if(fproc(1:4).eq.'stop') then

        call endmpi(1)

      end if

      stop

! -----

      end subroutine s_destroy

!-----7--------------------------------------------------------------7--

      end module m_destroy
