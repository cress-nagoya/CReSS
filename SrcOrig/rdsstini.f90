!***********************************************************************
      module m_rdsstini
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/11/10
!     Modification: 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the interpolated sea surface
!     temperature data file at forecast start time.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bcsst
      use m_chkerr
      use m_chkfile
      use m_chkstd
      use m_comindx
      use m_comkind
      use m_commpi
      use m_comname
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_getiname
      use m_getrname
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

      public :: rdsstini, s_rdsstini

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdsstini

        module procedure s_rdsstini

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic int

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rdsstini(fpexprim,fpcrsdir,fpncexp,fpnccrs,fpwlngth, &
     &                      fpstime,ni,nj,sst)
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

      integer, intent(in) :: fpwlngth
                       ! Formal parameter of unique index of wlngth

      integer, intent(in) :: fpstime
                       ! Formal parameter of unique index of stime

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

! Output variable

      real, intent(out) :: sst(0:ni+1,0:nj+1)
                       ! Sea surface temperature

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) sstfl
                       ! Opened file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer wlngth   ! Word length of direct access file

      integer in       ! Namelist table index

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer ncfl     ! Number of character of
                       ! sea surface temperature data file

      integer iosst    ! Unit number of
                       ! sea surface temperature data file

      integer siz      ! Record length of
                       ! sea surface temperature data file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

      real stime       ! Forecast start time

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
      call getiname(fpwlngth,wlngth)
      call getrname(fpstime,stime)

! -----

!! Open and read out the data from the interpolated sea surface
!! temperature data checking file at forecast start time.

! Initialize the character variable.

      if(mype.eq.root) then

        call inichar(sstfl)

      end if

! -----

! Get the unit number.

      if(mype.eq.root) then

        call getunit(iosst)

      end if

! -----

! Open the interpolated sea surface temperature data checking file at
! forecast start time.

      if(mype.eq.root) then

        sstfl(1:ncexp)=exprim(1:ncexp)

        write(sstfl(ncexp+1:ncexp+13),'(a13)') 'sst.check.txt'

        open(iosst,iostat=stat,err=100,                                 &
     &       file=crsdir(1:nccrs)//sstfl(1:ncexp+13),                   &
     &       status='old',access='sequential',form='formatted',         &
     &       blank='null',position='rewind',action='read')

      else

        stat=0

      end if

  100 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdsstini',8,'cont',1,'              ',14,iosst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdsstini',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdsstini',8,sstfl,ncexp+13,iosst,1,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the interpolated sea surface temperature data
! checking file at forecast start time.

      if(mype.eq.root) then

        read(iosst,'(a)',iostat=stat,end=110,err=110)                   &
     &      (rcname(in),in=1,ncn)

        read(iosst,*,iostat=stat,end=110,err=110)                       &
     &      (riname(in),in=1,nin)

        read(iosst,*,iostat=stat,end=110,err=110)                       &
     &      (rrname(in),in=1,nrn)

      else

        stat=0

      end if

  110 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdsstini',8,'cont',3,'              ',14,iosst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdsstini',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdsstini',8,sstfl,108,iosst,3,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Check the interpolated sea surface temperature data checking file at
! forecast start time.

      if(mype.eq.root) then

        call chkfile('sst',stat,ncn,nin,nrn,                            &
     &               cname,iname,rname,rcname,riname,rrname)

      else

        stat=0

      end if

      call chkerr(stat)

      if(stat.lt.0) then

        call destroy('chkfile ',7,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

! Close the interpolated sea surface temperature data checking file at
! forecast start time.

      if(mype.eq.root) then

        close(iosst,iostat=stat,err=120,status='keep')

      else

        stat=0

      end if

  120 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdsstini',8,'cont',2,'              ',14,iosst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdsstini',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdsstini',8,sstfl,108,iosst,2,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      if(mype.eq.root) then

        call putunit(iosst)

      end if

! -----

!! -----

!! Open and read out the data from the interpolated sea surface
!! temperature data file at forecast start time.

! Initialize the character variable.

      call inichar(sstfl)

! -----

! Get the unit number.

      call getunit(iosst)

! -----

! Open the interpolated sea surface temperature data file at forecast
! start time.

      siz=(ni+2)*(nj+2)*wlngth

      sstfl(1:ncexp)=exprim(1:ncexp)

      if(ngrp.eq.1) then

        ncfl=ncexp+22

        write(sstfl(ncexp+1:ncexp+22),'(a3,i8.8,a3,i4.4,a4)')           &
     &          'sst',int(stime+.1e0),'.pe',mysub,'.bin'

      else

        ncfl=ncexp+31

        write(sstfl(ncexp+1:ncexp+31),'(a3,i8.8,2(a4,i4.4),a4)')        &
     &     'sst',int(stime+.1e0),'.grp',mygrp,'-sub',mysub,'.bin'

      end if

      open(iosst,iostat=stat,err=130,                                   &
     &     file=crsdir(1:nccrs)//sstfl(1:ncfl),                         &
     &     status='old',access='direct',form='unformatted',             &
     &     recl=siz,action='read')

  130 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdsstini',8,'cont',1,'              ',14,iosst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdsstini',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        if(fpara(1:5).eq.'multi') then

          if(ngrp.eq.1) then

            write(sstfl(ncexp+15:ncexp+18),'(a4)') 'XXXX'

          else

            write(sstfl(ncexp+16:ncexp+19),'(a4)') 'XXXX'
            write(sstfl(ncexp+24:ncexp+27),'(a4)') 'YYYY'

          end if

        end if

        call outstd03('rdsstini',8,sstfl,ncfl,iosst,1,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the interpolated sea surface temperature data
! file at forecast start time.

      read(iosst,rec=1,iostat=stat,err=140)                             &
     &    ((sst(i,j),i=0,ni+1),j=0,nj+1)

  140 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdsstini',8,'cont',3,'              ',14,iosst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdsstini',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('rdsstini',8,sstfl,108,iosst,3,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the interpolated sea surface temperature data file at forecast
! start time.

      close(iosst,iostat=stat,err=150,status='keep')

  150 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdsstini',8,'cont',2,'              ',14,iosst, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdsstini',8,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('rdsstini',8,sstfl,108,iosst,2,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      call putunit(iosst)

! -----

!! -----

! Set the boundary conditions.

      call bcsst(idwbc,idebc,idexbopt,ni,nj,sst)

! -----

      end subroutine s_rdsstini

!-----7--------------------------------------------------------------7--

      end module m_rdsstini
