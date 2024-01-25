!***********************************************************************
      module m_rmres
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/08/05
!     Modification: 2006/09/21, 2006/12/04, 2007/01/05, 2007/01/20,
!                   2007/05/14, 2007/07/30, 2007/08/24, 2008/03/12,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     remove the original restart files.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chkstd
      use m_comkind
      use m_commpi
      use m_comname
      use m_cpondpe
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

      public :: rmres, s_rmres

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rmres

        module procedure s_rmres

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
      subroutine s_rmres(fpexprim,fpcrsdir,fpncexp,fpnccrs,ctime)
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

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) resfl
                       ! Opened file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer ncfl     ! Number of character of restart file name

      integer iores    ! Unit number of restart file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

!-----7--------------------------------------------------------------7--

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

!! Remove the original restart checking file.

! Initialize the character variable.

      if(mype.eq.root) then

        call inichar(resfl)

      end if

! -----

! Get the unit number.

      if(mype.eq.root) then

        call getunit(iores)

      end if

! -----

! Open the restart checking file.

      if(mype.eq.root) then

        resfl(1:ncexp)=exprim(1:ncexp)

        ncfl=ncexp+21

        write(resfl(ncexp+1:ncexp+21),'(a3,i8.8,a10)')                  &
     &                                'res',ctime/1000_i8,'.check.txt'

        open(iores,iostat=stat,err=100,                                 &
     &       file=crsdir(1:nccrs)//resfl(1:ncfl),                       &
     &       status='old',access='sequential',form='formatted',         &
     &       blank='null',action='readwrite')

      else

        stat=0

      end if

  100 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rmres   ',5,'cont',1,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rmres   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rmres   ',5,resfl,ncfl,iores,1,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the restart checking file.

      if(mype.eq.root) then

        close(iores,iostat=stat,err=110,status='delete')

      else

        stat=0

      end if

  110 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rmres   ',5,'cont',2,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rmres   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rmres   ',5,resfl,108,iores,2,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      if(mype.eq.root) then

        call putunit(iores)

      end if

! -----

!! -----

!! Remove the original restart files.

! Initialize the character variable.

      call inichar(resfl)

! -----

! Get the unit number.

      call getunit(iores)

! -----

! Open the restart file.

      resfl(1:ncexp)=exprim(1:ncexp)

      if(ngrp.eq.1) then

        ncfl=ncexp+22

        write(resfl(ncexp+1:ncexp+22),'(a3,i8.8,a3,i4.4,a4)')           &
     &            'res',ctime/1000_i8,'.pe',mysub,'.bin'

      else

        ncfl=ncexp+31

        write(resfl(ncexp+1:ncexp+31),'(a3,i8.8,2(a4,i4.4),a4)')        &
     &            'res',ctime/1000_i8,'.grp',mygrp,'-sub',mysub,'.bin'

      end if

      open(iores,iostat=stat,err=120,                                   &
     &     file=crsdir(1:nccrs)//resfl(1:ncfl),                         &
     &     status='old',access='sequential',form='unformatted',         &
     &     position='rewind',action='readwrite')

  120 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rmres   ',5,'cont',1,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rmres   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('rmres   ',5,resfl,ncfl,iores,1,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the restart file.

      close(iores,iostat=stat,err=130,status='delete')

  130 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rmres   ',5,'cont',2,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rmres   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('rmres   ',5,resfl,108,iores,2,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      call putunit(iores)

! -----

!! -----

      end subroutine s_rmres

!-----7--------------------------------------------------------------7--

      end module m_rmres
