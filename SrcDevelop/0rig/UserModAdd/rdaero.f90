!***********************************************************************
      module m_rdaero
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the 3 dimensional aerosol data
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

      public :: rdaero, s_rdaero

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdaero

        module procedure s_rdaero

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
      subroutine s_rdaero(fpdatdir,fpncdat,fpwlngth,cdate,ctime,stat,   &
     &                    nid,njd,nkd,nqa,zdat,qadat)
!***********************************************************************

! Input variables

      character(len=12), intent(in) :: cdate
                       ! Current forecast date
                       ! with Gregorian calendar, yyyymmddhhmm

      integer, intent(in) :: fpdatdir
                       ! Formal parameter of unique index of datdir

      integer, intent(in) :: fpncdat
                       ! Formal parameter of unique index of ncdat

      integer, intent(in) :: fpwlngth
                       ! Formal parameter of unique index of wlngth

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

! Output variables

      integer, intent(out) :: stat
                       ! Runtime status

      real, intent(out) :: zdat(1:nid,1:njd,1:nkd)
                       ! z physical coordinates in data

      real, intent(out) :: qadat(1:nid,1:njd,1:nkd,1:nqa(0))
                       ! Aerosol mixing ratio in data

! Internal shared variables

      character(len=108) datdir
                       ! User specified directory for external data

      character(len=108) aslfl
                       ! Opened file name

      integer ncdat    ! Number of character of datdir

      integer wlngth   ! Word length of direct access file

      integer id       ! Array index in x direction
      integer jd       ! Array index in y direction
      integer kd       ! Array index in z direction

      integer n        ! Array index in 4th direction

      integer ncfl     ! Number of character of aslfl

      integer ioasl    ! Unit number of aerosol data file

      integer recasl   ! Record number of aerosol data file

      integer siz      ! Record length of aerosol data file

!-----7--------------------------------------------------------------7--

!! Open and read out the data from the 3 dimensional aerosol data file.

! Initialize the character variable.

      call inichar(datdir)

! -----

! Get the required namelist variables.

      call getname(fpdatdir,datdir)
      call getname(fpncdat,ncdat)
      call getname(fpwlngth,wlngth)

! -----

! Initialize the character variable.

      call inichar(aslfl)

! -----

! Get the unit number.

      call getunit(ioasl)

! -----

! Initialize the runtime status.

      stat=0

! -----

! Open the 3 dimensional aerosol data file.

! ##### You will have to modify the following lines. #####

      siz=nid*njd*wlngth

      ncfl=24

      write(aslfl(1:ncfl),'(a8,a12,a4)') 'data.asl',cdate(1:12),'.bin'

      open(ioasl,iostat=stat,err=100,                                   &
     &     file=datdir(1:ncdat)//aslfl(1:ncfl),                         &
     &     status='old',access='direct',form='unformatted',             &
     &     recl=siz,action='read')

! ####################

  100 if(stat.ne.0) then

        call destroy('rdaero  ',6,'cont',1,'              ',14,ioasl,   &
     &               stat)

        return

      end if

      call outstd03('rdaero  ',6,aslfl,ncfl,ioasl,1,1,ctime)

! -----

! Read out the data from the 3 dimensional aerosol data file.

! ##### You will have to modify the following lines. #####

      recasl=0

      do kd=1,nkd

        recasl=recasl+1

        read(ioasl,rec=recasl,iostat=stat,err=110)                      &
     &      ((zdat(id,jd,kd),id=1,nid),jd=1,njd)

      end do

      do n=1,nqa(0)

        do kd=1,nkd

          recasl=recasl+1

          read(ioasl,rec=recasl,iostat=stat,err=110)                    &
     &        ((qadat(id,jd,kd,n),id=1,nid),jd=1,njd)

        end do

      end do

! ####################

  110 if(stat.ne.0) then

        call destroy('rdaero  ',6,'cont',3,'              ',14,ioasl,   &
     &               stat)

        return

      end if

      call outstd03('rdaero  ',6,aslfl,108,ioasl,3,1,ctime)

! -----

! Close the 3 dimensional aerosol data file.

      close(ioasl,iostat=stat,err=120,status='keep')

  120 if(stat.ne.0) then

        call destroy('rdaero  ',6,'cont',2,'              ',14,ioasl,   &
     &               stat)

        return

      end if

      call outstd03('rdaero  ',6,aslfl,108,ioasl,2,1,ctime)

! -----

! Return the unit number.

      call putunit(ioasl)

! -----

!! -----

      end subroutine s_rdaero

!-----7--------------------------------------------------------------7--

      end module m_rdaero
