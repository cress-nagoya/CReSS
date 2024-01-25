!***********************************************************************
      module m_rdsnd
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/06/14,
!                   1999/06/21, 1999/11/01, 1999/11/19, 2000/01/17,
!                   2000/04/18, 2001/01/15, 2001/05/29, 2001/10/17,
!                   2002/04/02, 2002/06/18, 2002/07/15, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2004/05/31,
!                   2005/02/10, 2005/04/04, 2006/09/21, 2006/12/04,
!                   2007/01/20, 2007/08/24, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/02/27, 2009/03/31, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the sounding data file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_destroy
      use m_getcname
      use m_getiname
      use m_getunit
      use m_inichar
      use m_outstd03
      use m_putunit

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rdsnd, s_rdsnd

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdsnd

        module procedure s_rdsnd

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
      subroutine s_rdsnd(fpexprim,fpcrsdir,fpncexp,fpnccrs,stat,        &
     &                   nsnd,zsnd,usnd,vsnd,ptsnd,qvsnd)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexprim
                       ! Formal parameter of unique index of exprim

      integer, intent(in) :: fpcrsdir
                       ! Formal parameter of unique index of crsdir

      integer, intent(in) :: fpncexp
                       ! Formal parameter of unique index of ncexp

      integer, intent(in) :: fpnccrs
                       ! Formal parameter of unique index of nccrs

      integer, intent(in) :: nsnd
                       ! Sounding data dimension

! Output variables

      integer, intent(out) :: stat
                       ! Runtime status

      real, intent(out) :: zsnd(1:nsnd)
                       ! z physical coordinates in sounding data

      real, intent(out) :: usnd(1:nsnd)
                       ! x components of velocity in sounding data

      real, intent(out) :: vsnd(1:nsnd)
                       ! y components of velocity in sounding data

      real, intent(out) :: ptsnd(1:nsnd)
                       ! Potential temrerature in sounding data

      real, intent(out) :: qvsnd(1:nsnd)
                       ! Water vapor mixing ratio in sounding data

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) sndfl
                       ! Opened file name

      character(len=1) comflg
                       ! Control flag to find comment line

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer isnd     ! Array index in z direction

      integer iosnd    ! Unit number of sounding file

!-----7--------------------------------------------------------------7--

!! Open and read out the data from the sounding data file.

! Initialize the character variables.

      call inichar(exprim)
      call inichar(crsdir)

! -----

! Get the required namelist variables.

      call getcname(fpexprim,exprim)
      call getcname(fpcrsdir,crsdir)
      call getiname(fpncexp,ncexp)
      call getiname(fpnccrs,nccrs)

! -----

! Initialize the character variable.

      call inichar(sndfl)

! -----

! Get the unit number.

      call getunit(iosnd)

! -----

! Initialize the runtime status.

      stat=0

! -----

! Open the sounding data file.

      sndfl(1:ncexp)=exprim(1:ncexp)

      write(sndfl(ncexp+1:ncexp+12),'(a12)') 'sounding.txt'

      open(iosnd,iostat=stat,err=100,                                   &
     &     file=crsdir(1:nccrs)//sndfl(1:ncexp+12),                     &
     &     status='old',access='sequential',form='formatted',           &
     &     blank='null',position='rewind',action='read')

  100 if(stat.ne.0) then

        call destroy('rdsnd   ',5,'cont',1,'              ',14,iosnd,   &
     &               stat)

        return

      end if

      call outstd03('rdsnd   ',5,sndfl,ncexp+12,iosnd,1,0,0_i8)

! -----

! Read out the data from the sounding data file.

      do_head: do

        read(iosnd,'(a1)',iostat=stat,end=110,err=110) comflg(1:1)

        if(comflg(1:1).ne.'#') then

          backspace(iosnd,iostat=stat,err=110)

          exit do_head

        end if

      end do do_head

      do isnd=nsnd,1,-1

        read(iosnd,*,iostat=stat,end=110,err=110)                       &
     &       zsnd(isnd),ptsnd(isnd),usnd(isnd),vsnd(isnd),qvsnd(isnd)

      end do

  110 if(stat.ne.0) then

        call destroy('rdsnd   ',5,'cont',3,'              ',14,iosnd,   &
     &               stat)

        return

      end if

      call outstd03('rdsnd   ',5,sndfl,108,iosnd,3,0,0_i8)

! -----

! Close the sounding data file.

      close(iosnd,iostat=stat,err=120,status='keep')

  120 if(stat.ne.0) then

        call destroy('rdsnd   ',5,'cont',2,'              ',14,iosnd,   &
     &               stat)

        return

      end if

      call outstd03('rdsnd   ',5,sndfl,108,iosnd,2,0,0_i8)

! -----

! Return the unit number.

      call putunit(iosnd)

! -----

!! -----

      end subroutine s_rdsnd

!-----7--------------------------------------------------------------7--

      end module m_rdsnd
