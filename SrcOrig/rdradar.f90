!***********************************************************************
      module m_rdradar
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2002/10/31, 2003/03/28, 2003/04/30, 2003/05/19,
!                   2003/06/27, 2004/05/31, 2005/02/10, 2007/01/20,
!                   2007/07/30, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/02/27, 2009/03/31, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the radar data file.

!-----7--------------------------------------------------------------7--

! ####################################################################
! ##                                                                ##
! ## You have to modify reading section to adjust your data set,    ##
! ## and set value lower than -1.0 x 10e34 at undefined points.     ##
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

      public :: rdradar, s_rdradar

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdradar

        module procedure s_rdradar

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
      subroutine s_rdradar(fpdatdir,fprdrvar,fpncdat,fpwlngth,          &
     &                     cdate,ctime,csec,stat,nid,njd,nkd,           &
     &                     zdat,udat,vdat,wdat,qpdat)
!***********************************************************************

! Input variables

      character(len=12), intent(in) :: cdate
                       ! Current forecast date
                       ! with Gregorian calendar, yyyymmddhhmm

      integer, intent(in) :: fpdatdir
                       ! Formal parameter of unique index of datdir

      integer, intent(in) :: fprdrvar
                       ! Formal parameter of unique index of rdrvar

      integer, intent(in) :: fpncdat
                       ! Formal parameter of unique index of ncdat

      integer, intent(in) :: fpwlngth
                       ! Formal parameter of unique index of wlngth

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer(kind=i8), intent(in) :: csec
                       ! mod(ctime, 60000)

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

! Output variables

      integer, intent(out) :: stat
                       ! Runtime status

      real, intent(out) :: zdat(1:nid,1:njd,1:nkd)
                       ! z physical coordinates in data

      real, intent(out) :: udat(1:nid,1:njd,1:nkd)
                       ! x components of velocity in data

      real, intent(out) :: vdat(1:nid,1:njd,1:nkd)
                       ! y components of velocity in data

      real, intent(out) :: wdat(1:nid,1:njd,1:nkd)
                       ! z components of velocity in data

      real, intent(out) :: qpdat(1:nid,1:njd,1:nkd)
                       ! Precipitation mixing ratio in data

! Internal shared variables

      character(len=108) datdir
                       ! User specified directory for external data

      character(len=108) rdrvar
                       ! Control flag of input radar data variables

      character(len=108) rdrfl
                       ! Opened file name

      integer ncdat    ! Number of character of datdir

      integer wlngth   ! Word length of direct access file

      integer id       ! Array index in x direction
      integer jd       ! Array index in y direction
      integer kd       ! Array index in z direction

      integer ncfl     ! Number of character of rdrfl

      integer iordr    ! Unit number of radar data file

      integer recrdr   ! Record number of radar data file

      integer siz      ! Record length of radar data file

!-----7--------------------------------------------------------------7--

!! Open and read out the data from the radar data file.

! Initialize the character variables.

      call inichar(datdir)
      call inichar(rdrvar)

! -----

! Get the required namelist variables.

      call getname(fpdatdir,datdir)
      call getname(fprdrvar,rdrvar)
      call getname(fpncdat,ncdat)
      call getname(fpwlngth,wlngth)

! -----

! Initialize the character variable.

      call inichar(rdrfl)

! -----

! Get the unit number.

      call getunit(iordr)

! -----

! Initialize the runtime status.

      stat=0

! -----

! Open the radar data file.

! ##### You will have to modify the following lines. #####

      siz=nid*njd*wlngth

      ncfl=28

      write(rdrfl(1:ncfl),'(a10,a12,i2.2,a4)')                          &
     &                    'data.radar',cdate(1:12),csec,'.bin'

      open(iordr,iostat=stat,err=100,                                   &
     &     file=datdir(1:ncdat)//rdrfl(1:ncfl),                         &
     &     status='old',access='direct',form='unformatted',             &
     &     recl=siz,action='read')

! ####################

  100 if(stat.ne.0) then

        call destroy('rdradar ',7,'cont',1,'              ',14,iordr,   &
     &               stat)

        return

      end if

      call outstd03('rdradar ',7,rdrfl,ncfl,iordr,1,1,ctime)

! -----

! Read out the data from the radar data file.

! ##### You will have to modify the following lines. #####

      recrdr=0

      do kd=1,nkd

        recrdr=recrdr+1

        read(iordr,rec=recrdr,iostat=stat,err=110)                      &
     &      ((zdat(id,jd,kd),id=1,nid),jd=1,njd)

      end do

      if(rdrvar(1:1).eq.'o') then

        do kd=1,nkd

          recrdr=recrdr+1

          read(iordr,rec=recrdr,iostat=stat,err=110)                    &
     &        ((udat(id,jd,kd),id=1,nid),jd=1,njd)

        end do

      end if

      if(rdrvar(2:2).eq.'o') then

        do kd=1,nkd

          recrdr=recrdr+1

          read(iordr,rec=recrdr,iostat=stat,err=110)                    &
     &        ((vdat(id,jd,kd),id=1,nid),jd=1,njd)

        end do

      end if

      if(rdrvar(3:3).eq.'o') then

        do kd=1,nkd

          recrdr=recrdr+1

          read(iordr,rec=recrdr,iostat=stat,err=110)                    &
     &        ((wdat(id,jd,kd),id=1,nid),jd=1,njd)

        end do

      end if

      if(rdrvar(4:4).eq.'o') then

        do kd=1,nkd

          recrdr=recrdr+1

          read(iordr,rec=recrdr,iostat=stat,err=110)                    &
     &        ((qpdat(id,jd,kd),id=1,nid),jd=1,njd)

        end do

      end if

! ####################

  110 if(stat.ne.0) then

        call destroy('rdradar ',7,'cont',3,'              ',14,iordr,   &
     &               stat)

        return

      end if

      call outstd03('rdradar ',7,rdrfl,108,iordr,3,1,ctime)

! -----

! Close the radar data file.

      close(iordr,iostat=stat,err=120,status='keep')

  120 if(stat.ne.0) then

        call destroy('rdradar ',7,'cont',2,'              ',14,iordr,   &
     &               stat)

        return

      end if

      call outstd03('rdradar ',7,rdrfl,108,iordr,2,1,ctime)

! -----

! Return the unit number.

      call putunit(iordr)

! -----

!! -----

      end subroutine s_rdradar

!-----7--------------------------------------------------------------7--

      end module m_rdradar
