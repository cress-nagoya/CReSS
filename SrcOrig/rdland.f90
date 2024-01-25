!***********************************************************************
      module m_rdland
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/07/03
!     Modification: 2002/07/15, 2002/09/09, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2003/06/27, 2004/05/31, 2005/02/10,
!                   2007/01/20, 2007/08/24, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/02/27, 2009/03/31, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the 2 dimensional land use data
!     file.

!-----7--------------------------------------------------------------7--

! ####################################################################
! ##                                                                ##
! ## You have to modify reading section to adjust your data set,    ##
! ## and reset the your data categories after reading stetment.     ##
! ##                                                                ##
! ##   water(including river and lake)  --->      ~ -1              ##
! ##   ice                              --->   0  ~  4              ##
! ##   snow                             --->   5  ~  9              ##
! ##   the others                       --->  10  ~                 ##
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

      public :: rdland, s_rdland

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdland

        module procedure s_rdland

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
      subroutine s_rdland(fpdatdir,fpncdat,fpwlngth,stat,nid,njd,landat)
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

      integer, intent(out) :: landat(1:nid,1:njd)
                       ! Land use in data

! Internal shared variables

      character(len=108) datdir
                       ! User specified directory for external data

      character(len=108) lndfl
                       ! Opened file name

      integer ncdat    ! Number of character of datdir

      integer wlngth   ! Word length of direct access file

      integer id       ! Array index in x direction
      integer jd       ! Array index in y direction

      integer ncfl     ! Number of character of lndfl

      integer iolnd    ! Unit number of land use file

      integer reclnd   ! Record number of land use file

      integer siz      ! Record length of land use file

!-----7--------------------------------------------------------------7--

!! Open and read out the data from the 2 dimensional land use data file.

! Initialize the character variable.

      call inichar(datdir)

! -----

! Get the required namelist variables.

      call getname(fpdatdir,datdir)
      call getname(fpncdat,ncdat)
      call getname(fpwlngth,wlngth)

! -----

! Initialize the character variable.

      call inichar(lndfl)

! -----

! Get the unit number.

      call getunit(iolnd)

! -----

! Initialize the runtime status.

      stat=0

! -----

! Open the 2 dimensional land use data file.

! ##### You will have to modify the following lines. #####

      siz=nid*njd*wlngth

      ncfl=13

      write(lndfl(1:ncfl),'(a13)') 'data.land.bin'

      open(iolnd,iostat=stat,err=100,                                   &
     &     file=datdir(1:ncdat)//lndfl(1:ncfl),                         &
     &     status='old',access='direct',form='unformatted',             &
     &     recl=siz,action='read')

! ####################

  100 if(stat.ne.0) then

        call destroy('rdland  ',6,'cont',1,'              ',14,iolnd,   &
     &               stat)

        return

      end if

      call outstd03('rdland  ',6,lndfl,ncfl,iolnd,1,0,0_i8)

! -----

! Read out the data from the 2 dimensional land use data file.

! ##### You will have to modify the following lines. #####

      reclnd=1

      read(iolnd,rec=reclnd,iostat=stat,err=110)                        &
     &    ((landat(id,jd),id=1,nid),jd=1,njd)

! ####################

  110 if(stat.ne.0) then

        call destroy('rdland  ',6,'cont',3,'              ',14,iolnd,   &
     &               stat)

        return

      end if

      call outstd03('rdland  ',6,lndfl,108,iolnd,3,0,0_i8)

! -----

! Close the 2 dimensional land use data file.

      close(iolnd,iostat=stat,err=120,status='keep')

  120 if(stat.ne.0) then

        call destroy('rdland  ',6,'cont',2,'              ',14,iolnd,   &
     &               stat)

        return

      end if

      call outstd03('rdland  ',6,lndfl,108,iolnd,2,0,0_i8)

! -----

! Return the unit number.

      call putunit(iolnd)

! -----

!! -----

      end subroutine s_rdland

!-----7--------------------------------------------------------------7--

      end module m_rdland
