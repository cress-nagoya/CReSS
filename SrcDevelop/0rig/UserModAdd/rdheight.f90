!***********************************************************************
      module m_rdheight
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/04/06
!     Modification: 1999/05/10, 1999/06/14, 1999/06/21, 1999/11/01,
!                   2000/01/17, 2001/02/13, 2001/04/15, 2001/05/29,
!                   2002/04/02, 2002/06/18, 2002/07/15, 2002/09/09,
!                   2003/03/28, 2003/04/30, 2003/05/19, 2003/06/27,
!                   2004/05/31, 2005/02/10, 2007/01/20, 2007/08/24,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2009/03/31, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the 2 dimensional terrain data
!     file.

!-----7--------------------------------------------------------------7--

! ####################################################################
! ##                                                                ##
! ## You have to modify reading section to adjust your data set.    ##
! ##                                                                ##
! ####################################################################

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_destroy
      use m_getname
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

      public :: rdheight, s_rdheight

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdheight

        module procedure s_rdheight

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
      subroutine s_rdheight(fpdatdir,fpncdat,fpwlngth,stat,             &
     &                      nid,njd,htdat)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdatdir
                       ! Formal parameter of unique index of datdir

      integer, intent(in) :: fpncdat
                       ! Formal parameter of unique index of ncdat

      integer, intent(in) :: fpwlngth
                       ! Formal parameter of unique index of wlngth

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

! Output variables

      integer, intent(out) :: stat
                       ! Runtime status

      real, intent(out) :: htdat(1:nid,1:njd)
                       ! Terrain height in data

! Internal shared variables

      character(len=108) datdir
                       ! User specified directory for external data

      character(len=108) trnfl
                       ! Opened file name

      integer ncdat    ! Number of character of datdir

      integer wlngth   ! Word length of direct access file

      integer id       ! Array index in x direction
      integer jd       ! Array index in y direction

      integer ncfl     ! Number of character of trnfl

      integer iotrn    ! Unit number of terrain file

      integer rectrn   ! Record number of terrain file

      integer siz      ! Record length of terrain file

!-----7--------------------------------------------------------------7--

!! Open and read out the data from the 2 dimensional terrain data file.

! Initialize the character variable.

      call inichar(datdir)

! -----

! Get the required namelist variables.

      call getname(fpdatdir,datdir)
      call getname(fpncdat,ncdat)
      call getname(fpwlngth,wlngth)

! -----

! Initialize the character variable.

      call inichar(trnfl)

! -----

! Get the unit number.

      call getunit(iotrn)

! -----

! Initialize the runtime status.

      stat=0

! -----

! Open the 2 dimensional terrain data file.

! ##### You will have to modify the following lines. #####

      siz=nid*njd*wlngth

      ncfl=16

      write(trnfl(1:ncfl),'(a16)') 'data.terrain.bin'

      open(iotrn,iostat=stat,err=100,                                   &
     &     file=datdir(1:ncdat)//trnfl(1:ncfl),                         &
     &     status='old',access='direct',form='unformatted',             &
     &     recl=siz,action='read')

! ####################

  100 if(stat.ne.0) then

        call destroy('rdheight',8,'cont',1,'              ',14,iotrn,   &
     &               stat)

        return

      end if

      call outstd03('rdheight',8,trnfl,ncfl,iotrn,1,0,0_i8)

! -----

! Read out the data from the 2 dimensional terrain data file.

! ##### You will have to modify the following lines. #####

      rectrn=1

      read(iotrn,rec=rectrn,iostat=stat,err=110)                        &
     &    ((htdat(id,jd),id=1,nid),jd=1,njd)

! ####################

  110 if(stat.ne.0) then

        call destroy('rdheight',8,'cont',3,'              ',14,iotrn,   &
     &               stat)

        return

      end if

      call outstd03('rdheight',8,trnfl,108,iotrn,3,0,0_i8)

! -----

! Close the 2 dimensional terrain data file.

      close(iotrn,iostat=stat,err=120,status='keep')

  120 if(stat.ne.0) then

        call destroy('rdheight',8,'cont',2,'              ',14,iotrn,   &
     &               stat)

        return

      end if

      call outstd03('rdheight',8,trnfl,108,iotrn,2,0,0_i8)

! -----

! Return the unit number.

      call putunit(iotrn)

! -----

!! -----

      end subroutine s_rdheight

!-----7--------------------------------------------------------------7--

      end module m_rdheight
