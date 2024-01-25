!***********************************************************************
      module m_outstd07
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/01/15
!     Modification: 2001/04/15, 2001/05/29, 2001/12/10, 2002/04/02,
!                   2002/07/03, 2003/03/28, 2003/04/30, 2003/05/19,
!                   2003/09/01, 2003/12/12, 2004/03/05, 2004/06/10,
!                   2004/08/01, 2007/09/19, 2008/01/11, 2008/04/17,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the messages to standard i/o.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outstd07, s_outstd07

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outstd07

        module procedure s_outstd07

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
      subroutine s_outstd07(fpdmplev,nk,z1d)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdmplev
                       ! Formal parameter of unique index of dmplev

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: z1d(1:nk)
                       ! Output z coordinates at scalar point

! Internal shared variables

      integer dmplev   ! Option for z coordinates of dumped variables

      integer nkm2     ! nk - 2

      integer k        ! Array index in z drection

      integer icnt     ! Counter of writing column

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpdmplev,dmplev)

! -----

! Set the common used variable.

      nkm2=nk-2

! -----

! Read in the messages to standard i/o.

      write(6,*)

      write(6,'(a)') '  messages: procedure, opendmp;'

      if(dmplev.lt.10) then

        write(6,'(a,i1,a)',advance='no')                                &
     &          '    Option, dmplev = ',dmplev,','

      else

        write(6,'(a,i2,a)',advance='no')                                &
     &          '    Option, dmplev = ',dmplev,','

      end if

      write(6,'(a)')                                                    &
     &        ' grid points under the ground and above the model top'

      write(6,'(a)')                                                    &
     &        '    have the value -1.0 x 10^35 for the 3d variable(s).'

      write(6,*)

      icnt=0

      do k=2,nk-2

        icnt=icnt+1

        if(icnt.eq.1) then

          if(k.eq.2) then

            write(6,'(a)',advance='no') '    output layers [m] = '

          else

            write(6,'(a)',advance='no') '                        '

          end if

        end if

        if(k.eq.nkm2.or.icnt.eq.3) then

          if(k.eq.nkm2) then

            write(6,'(e13.6e2)') z1d(k)

          else

            write(6,'(e13.6e2,a)') z1d(k),', '

          end if

          icnt=0

        else

          write(6,'(e13.6e2,a)',advance='no') z1d(k),', '

        end if

      end do

! -----

      end subroutine s_outstd07

!-----7--------------------------------------------------------------7--

      end module m_outstd07
