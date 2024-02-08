!***********************************************************************
      module m_rdice
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/07/15
!     Modification: 2004/05/31, 2005/02/10, 2007/01/20, 2007/08/24,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2009/03/31, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the 2 dimensional sea ice
!     distribution data file.

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

      public :: rdice, s_rdice

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdice

        module procedure s_rdice

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
      subroutine s_rdice(fpdatdir,fpncdat,fpwlngth,stat,nid,njd,icedat)
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

      real, intent(out) :: icedat(1:nid,1:njd)
                       ! Sea ice distribution in data

! Internal shared variables

      character(len=108) datdir
                       ! User specified directory for external data

      character(len=108) icefl
                       ! Opened file name

      integer ncdat    ! Number of character of datdir

      integer wlngth   ! Word length of direct access file

      integer id       ! Array index in x direction
      integer jd       ! Array index in y direction

      integer ncfl     ! Number of character of icefl

      integer ioice    ! Unit number of sea ice distribution file

      integer recice   ! Record number of sea ice distribution file

      integer siz      ! Record length of sea ice distribution file

!-----7--------------------------------------------------------------7--

!! Open and read out the data from the 2 dimensional sea ice
!! distribution data file.

! Initialize the character variable.

      call inichar(datdir)

! -----

! Get the required namelist variables.

      call getname(fpdatdir,datdir)
      call getname(fpncdat,ncdat)
      call getname(fpwlngth,wlngth)

! -----

! Initialize the character variable.

      call inichar(icefl)

! -----

! Get the unit number.

      call getunit(ioice)

! -----

! Initialize the runtime status.

      stat=0

! -----

! Open the 2 dimensional sea ice distribution data file.

! ##### You will have to modify the following lines. #####

      siz=nid*njd*wlngth

      ncfl=12

      write(icefl(1:ncfl),'(a12)') 'data.ice.bin'

      open(ioice,iostat=stat,err=100,                                   &
     &     file=datdir(1:ncdat)//icefl(1:ncfl),                         &
     &     status='old',access='direct',form='unformatted',             &
     &     recl=siz,action='read')

! ####################

  100 if(stat.ne.0) then

        call destroy('rdice   ',5,'cont',1,'              ',14,ioice,   &
     &               stat)

        return

      end if

      call outstd03('rdice   ',5,icefl,ncfl,ioice,1,0,0_i8)

! -----

! Read out the data from the 2 dimensional sea ice distribution data
! file.

! ##### You will have to modify the following lines. #####

      recice=1

      read(ioice,rec=recice,iostat=stat,err=110)                        &
     &    ((icedat(id,jd),id=1,nid),jd=1,njd)

! ####################

  110 if(stat.ne.0) then

        call destroy('rdice   ',5,'cont',3,'              ',14,ioice,   &
     &               stat)

        return

      end if

      call outstd03('rdice   ',5,icefl,108,ioice,3,0,0_i8)

! -----

! Close the 2 dimensional sea ice distribution data file.

      close(ioice,iostat=stat,err=120,status='keep')

  120 if(stat.ne.0) then

        call destroy('rdice   ',5,'cont',2,'              ',14,ioice,   &
     &               stat)

        return

      end if

      call outstd03('rdice   ',5,icefl,108,ioice,2,0,0_i8)

! -----

! Return the unit number.

      call putunit(ioice)

! -----

!! -----

      end subroutine s_rdice

!-----7--------------------------------------------------------------7--

      end module m_rdice
