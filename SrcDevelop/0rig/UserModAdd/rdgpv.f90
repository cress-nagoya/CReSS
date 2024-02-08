!***********************************************************************
      module m_rdgpv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/20
!     Modification: 1999/06/14, 1999/06/21, 2000/01/05, 2000/01/17,
!                   2001/01/15, 2001/02/13, 2001/03/13, 2001/04/15,
!                   2001/05/29, 2002/04/02, 2002/04/09, 2002/06/18,
!                   2002/07/15, 2002/09/02, 2002/09/09, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2004/03/05,
!                   2004/04/15, 2004/05/31, 2005/02/10, 2007/01/20,
!                   2007/07/30, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/02/27, 2009/03/31, 2011/09/22, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the 3 dimensional GPV data file.

!-----7--------------------------------------------------------------7--

! ####################################################################
! ##                                                                ##
! ## You have to modify reading section to adjust your data set.    ##
! ##                                                                ##
! ## The data array at kd = 1 are on the surface in the case the    ##
! ## namelist variable refsfc is 1 in the user configuration file.  ##
! ##                                                                ##
! ## The htdat do not have to set in the case the namelist variable ##
! ## exbopt is neither 2 nor 12 in the user configuration file.     ##
! ##                                                                ##
! ## The zdat at kd = 1 may be equal to the htdat in the case the   ##
! ## namelist variable exbopt is 2 or 12 in the user configuration  ##
! ## file.                                                          ##
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

      public :: rdgpv, s_rdgpv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdgpv

        module procedure s_rdgpv

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
      subroutine s_rdgpv(fpdatdir,fpgpvvar,fpncdat,fpwlngth,fpexbopt,   &
     &                   cdate,ctime,stat,nid,njd,nkd,htdat,zdat,       &
     &                   udat,vdat,wdat,pdat,ptdat,qvdat,qcdat,qrdat,   &
     &                   qidat,qsdat,qgdat,qhdat)
!***********************************************************************

! Input variables

      character(len=12), intent(in) :: cdate
                       ! Current forecast date
                       ! with Gregorian calendar, yyyymmddhhmm

      integer, intent(in) :: fpdatdir
                       ! Formal parameter of unique index of datdir

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fpncdat
                       ! Formal parameter of unique index of ncdat

      integer, intent(in) :: fpwlngth
                       ! Formal parameter of unique index of wlngth

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

! Output variables

      integer, intent(out) :: stat
                       ! Runtime status

      real, intent(out) :: htdat(1:nid,1:njd)
                       ! Terrain height in data

      real, intent(out) :: zdat(1:nid,1:njd,1:nkd)
                       ! z physical coordinates in data

      real, intent(out) :: udat(1:nid,1:njd,1:nkd)
                       ! x components of velocity in data

      real, intent(out) :: vdat(1:nid,1:njd,1:nkd)
                       ! y components of velocity in data

      real, intent(out) :: wdat(1:nid,1:njd,1:nkd)
                       ! z components of velocity in data

      real, intent(out) :: pdat(1:nid,1:njd,1:nkd)
                       ! Pressure in data

      real, intent(out) :: ptdat(1:nid,1:njd,1:nkd)
                       ! Potential temperature in data

      real, intent(out) :: qvdat(1:nid,1:njd,1:nkd)
                       ! Water vapor mixing ratio in data

      real, intent(out) :: qcdat(1:nid,1:njd,1:nkd)
                       ! Cloud water mixing ratio in data

      real, intent(out) :: qrdat(1:nid,1:njd,1:nkd)
                       ! Rain water mixing ratio in data

      real, intent(out) :: qidat(1:nid,1:njd,1:nkd)
                       ! Cloud ice mixing ratio in data

      real, intent(out) :: qsdat(1:nid,1:njd,1:nkd)
                       ! Snow mixing ratio in data

      real, intent(out) :: qgdat(1:nid,1:njd,1:nkd)
                       ! Graupel mixing ratio in data

      real, intent(out) :: qhdat(1:nid,1:njd,1:nkd)
                       ! Hail mixing ratio in data

! Internal shared variables

      character(len=108) datdir
                       ! User specified directory for external data

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      character(len=108) gpvfl
                       ! Opened file name

      integer ncdat    ! Number of character of datdir

      integer wlngth   ! Word length of direct access file

      integer exbopt   ! Option for external boundary forcing

      integer id       ! Array index in x direction
      integer jd       ! Array index in y direction
      integer kd       ! Array index in z direction

      integer ncfl     ! Number of character of gpvfl

      integer iogpv    ! Unit number of GPV data file

      integer recgpv   ! Record number of GPV data file

      integer siz      ! Record length of GPV data file

!-----7--------------------------------------------------------------7--

!! Open and read out the data from the 3 dimensional GPV data file.

! Initialize the character variables.

      call inichar(datdir)
      call inichar(gpvvar)

! -----

! Get the required namelist variables.

      call getname(fpdatdir,datdir)
      call getname(fpgpvvar,gpvvar)
      call getname(fpncdat,ncdat)
      call getname(fpwlngth,wlngth)
      call getname(fpexbopt,exbopt)

! -----

! Initialize the character variable.

      call inichar(gpvfl)

! -----

! Get the unit number.

      call getunit(iogpv)

! -----

! Initialize the runtime status.

      stat=0

! -----

! Open the 3 dimensional GPV data file.

! ##### You will have to modify the following lines. #####

      siz=nid*njd*wlngth

      ncfl=24

      write(gpvfl(1:ncfl),'(a8,a12,a4)') 'data.gpv',cdate(1:12),'.bin'

      open(iogpv,iostat=stat,err=100,                                   &
     &     file=datdir(1:ncdat)//gpvfl(1:ncfl),                         &
     &     status='old',access='direct',form='unformatted',             &
     &     recl=siz,action='read')

! ####################

  100 if(stat.ne.0) then

        call destroy('rdgpv   ',5,'cont',1,'              ',14,iogpv,   &
     &               stat)

        return

      end if

      call outstd03('rdgpv   ',5,gpvfl,ncfl,iogpv,1,1,ctime)

! -----

! Read out the data from the 3 dimensional GPV data file.

! ##### You will have to modify the following lines. #####

      recgpv=0

      do kd=1,nkd

        recgpv=recgpv+1

        read(iogpv,rec=recgpv,iostat=stat,err=110)                      &
     &      ((zdat(id,jd,kd),id=1,nid),jd=1,njd)

      end do

      do kd=1,nkd

        recgpv=recgpv+1

        read(iogpv,rec=recgpv,iostat=stat,err=110)                      &
     &      ((udat(id,jd,kd),id=1,nid),jd=1,njd)

      end do

      do kd=1,nkd

        recgpv=recgpv+1

        read(iogpv,rec=recgpv,iostat=stat,err=110)                      &
     &      ((vdat(id,jd,kd),id=1,nid),jd=1,njd)

      end do

      do kd=1,nkd

        recgpv=recgpv+1

        read(iogpv,rec=recgpv,iostat=stat,err=110)                      &
     &      ((pdat(id,jd,kd),id=1,nid),jd=1,njd)

      end do

      do kd=1,nkd

        recgpv=recgpv+1

        read(iogpv,rec=recgpv,iostat=stat,err=110)                      &
     &      ((ptdat(id,jd,kd),id=1,nid),jd=1,njd)

      end do

      if(gpvvar(1:1).eq.'o') then

        do kd=1,nkd

          recgpv=recgpv+1

          read(iogpv,rec=recgpv,iostat=stat,err=110)                    &
     &        ((wdat(id,jd,kd),id=1,nid),jd=1,njd)

        end do

      end if

      if(gpvvar(2:2).eq.'o') then

        do kd=1,nkd

          recgpv=recgpv+1

          read(iogpv,rec=recgpv,iostat=stat,err=110)                    &
     &        ((qvdat(id,jd,kd),id=1,nid),jd=1,njd)

        end do

      end if

      if(gpvvar(3:3).eq.'o') then

        do kd=1,nkd

          recgpv=recgpv+1

          read(iogpv,rec=recgpv,iostat=stat,err=110)                    &
     &        ((qcdat(id,jd,kd),id=1,nid),jd=1,njd)

        end do

      end if

      if(gpvvar(4:4).eq.'o') then

        do kd=1,nkd

          recgpv=recgpv+1

          read(iogpv,rec=recgpv,iostat=stat,err=110)                    &
     &        ((qrdat(id,jd,kd),id=1,nid),jd=1,njd)

        end do

      end if

      if(gpvvar(5:5).eq.'o') then

        do kd=1,nkd

          recgpv=recgpv+1

          read(iogpv,rec=recgpv,iostat=stat,err=110)                    &
     &        ((qidat(id,jd,kd),id=1,nid),jd=1,njd)

        end do

      end if

      if(gpvvar(6:6).eq.'o') then

        do kd=1,nkd

          recgpv=recgpv+1

          read(iogpv,rec=recgpv,iostat=stat,err=110)                    &
     &        ((qsdat(id,jd,kd),id=1,nid),jd=1,njd)

        end do

      end if

      if(gpvvar(7:7).eq.'o') then

        do kd=1,nkd

          recgpv=recgpv+1

          read(iogpv,rec=recgpv,iostat=stat,err=110)                    &
     &        ((qgdat(id,jd,kd),id=1,nid),jd=1,njd)

        end do

      end if

      if(gpvvar(8:8).eq.'o') then

        do kd=1,nkd

          recgpv=recgpv+1

          read(iogpv,rec=recgpv,iostat=stat,err=110)                    &
     &        ((qhdat(id,jd,kd),id=1,nid),jd=1,njd)

        end do

      end if

      if(exbopt.eq.2.or.exbopt.eq.12) then

        do jd=1,njd
        do id=1,nid
          htdat(id,jd)=zdat(id,jd,1)
        end do
        end do

      end if

! ####################

  110 if(stat.ne.0) then

        call destroy('rdgpv   ',5,'cont',3,'              ',14,iogpv,   &
     &               stat)

        return

      end if

      call outstd03('rdgpv   ',5,gpvfl,108,iogpv,3,1,ctime)

! -----

! Close the 3 dimensional GPV data file.

      close(iogpv,iostat=stat,err=120,status='keep')

  120 if(stat.ne.0) then

        call destroy('rdgpv   ',5,'cont',2,'              ',14,iogpv,   &
     &               stat)

        return

      end if

      call outstd03('rdgpv   ',5,gpvfl,108,iogpv,2,1,ctime)

! -----

! Return the unit number.

      call putunit(iogpv)

! -----

!! -----

      end subroutine s_rdgpv

!-----7--------------------------------------------------------------7--

      end module m_rdgpv
