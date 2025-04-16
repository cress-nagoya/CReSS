!***********************************************************************
      module m_rdsstnxt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/11/10
!     Modification: 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read out the data from the interpolated sea surface temperature
!     data file.

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
      use m_comsave
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_getiname
      use m_getrname
      use m_getunit
      use m_inichar
      use m_outstd03
      use m_putunit
      use m_setsst

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rdsstnxt, s_rdsstnxt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdsstnxt

        module procedure s_rdsstnxt

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic int
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rdsstnxt(fpexprim,fpcrsdir,fpncexp,fpnccrs,fpwlngth, &
     &                      fpetime,fpsstitv,fsst,ctime,ftime,stinc,    &
     &                      ni,nj,sst,sstd)
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

      integer, intent(in) :: fpetime
                       ! Formal parameter of unique index of etime

      integer, intent(in) :: fpsstitv
                       ! Formal parameter of unique index of sstitv

      integer, intent(in) :: fsst
                       ! Descriptor to put into motion
                       ! for sea surface temperature data reading

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer(kind=i8), intent(in) :: ftime
                       ! Model forecast time at 1 step future

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

! Input and output variables

      real, intent(inout) :: sst(0:ni+1,0:nj+1)
                       ! Sea surface temperature of external data
                       ! at marked time

      real, intent(inout) :: sstd(0:ni+1,0:nj+1)
                       ! Time tendency of
                       ! sea surface temperature of external data

! Output variable

      real, intent(out) :: stinc
                       ! Lapse of forecast time
                       ! from sea surface temperature data reading

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

      integer ird      ! Index of count to read out

      integer ncfl     ! Number of character of
                       ! sea surface temperature data file

      integer iosst    ! Unit number of
                       ! sea surface temperature data file

      integer siz      ! Record length of
                       ! sea surface temperature data file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

      integer datest   ! Sea surface temperature data file date

      integer sst01    ! int(sstitv + 0.1)

      integer(kind=i8) sst103
                       ! 1000 x int(sstitv + 0.1)

      integer(kind=i8) crtime
                       ! Current forecast time to read out

      real etime       ! Forecast stop time

      real sstitv      ! Time interval of sea surface temperature data

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpetime,etime)
      call getrname(fpsstitv,sstitv)

! -----

! Set the common used variables.

      sst01=int(sstitv+.1e0)

      sst103=1000_i8*int(sstitv+.1e0,i8)

      crtime=sst103*(ctime/sst103)

      stinc=.001e0*real(ctime-crtime)

      if(fsst.eq.0) then

        nxtsst=crtime/1000_i8

      end if

! -----

!!! Open and read out the data from the interpolated sea surface
!!! temperature data file when the current forecast time reaches
!!! marked time.

      if(ctime/1000_i8.ge.nxtsst                                        &
     &  .and.ctime.lt.1000_i8*int(etime+.1e0,i8)) then

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

! -----

!! Open and read out the data from the interpolated sea surface
!! temperature data checking file.

        if(fsst.eq.0) then

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

! Open the interpolated sea surface temperature data checking file.

          if(mype.eq.root) then

            sstfl(1:ncexp)=exprim(1:ncexp)

            write(sstfl(ncexp+1:ncexp+13),'(a13)') 'sst.check.txt'

            open(iosst,iostat=stat,err=100,                             &
     &           file=crsdir(1:nccrs)//sstfl(1:ncexp+13),               &
     &           status='old',access='sequential',form='formatted',     &
     &           blank='null',position='rewind',action='read')

          else

            stat=0

          end if

  100     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdsstnxt',8,'cont',1,'              ',14,   &
     &                     iosst,stat)

            end if

            call cpondpe

            call destroy('rdsstnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.root) then

            call outstd03('rdsstnxt',8,sstfl,ncexp+13,iosst,1,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Read out the data from the interpolated sea surface temperature data
! checking file.

          if(mype.eq.root) then

            read(iosst,'(a)',iostat=stat,end=110,err=110)               &
     &          (rcname(in),in=1,ncn)

            read(iosst,*,iostat=stat,end=110,err=110)                   &
     &          (riname(in),in=1,nin)

            read(iosst,*,iostat=stat,end=110,err=110)                   &
     &          (rrname(in),in=1,nrn)

          else

            stat=0

          end if

  110     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdsstnxt',8,'cont',3,'              ',14,   &
     &                     iosst,stat)

            end if

            call cpondpe

            call destroy('rdsstnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.root) then

            call outstd03('rdsstnxt',8,sstfl,108,iosst,3,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Check the interpolated sea surface temperature data checking file.

          if(mype.eq.root) then

            call chkfile('sst',stat,ncn,nin,nrn,                        &
     &                   cname,iname,rname,rcname,riname,rrname)

          else

            stat=0

          end if

          call chkerr(stat)

          if(stat.lt.0) then

            call destroy('chkfile ',7,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

! -----

! Close the interpolated sea surface temperature data checking file.

          if(mype.eq.root) then

            close(iosst,iostat=stat,err=120,status='keep')

          else

            stat=0

          end if

  120     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdsstnxt',8,'cont',2,'              ',14,   &
     &                     iosst,stat)

            end if

            call cpondpe

            call destroy('rdsstnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.root) then

            call outstd03('rdsstnxt',8,sstfl,108,iosst,2,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Return the unit number.

          if(mype.eq.root) then

            call putunit(iosst)

          end if

! -----

        end if

!! -----

!! Open and read out the data from the interpolated sea surface
!! temperature data file.

        siz=(ni+2)*(nj+2)*wlngth

        do ird=2,1,-1

! Initialize the character variable.

          call inichar(sstfl)

! -----

! Get the unit number.

          call getunit(iosst)

! -----

! Open the interpolated sea surface temperature data file.

          datest=crtime/1000_i8
          datest=datest+abs(ird-2)*sst01

          sstfl(1:ncexp)=exprim(1:ncexp)

          if(ngrp.eq.1) then

           ncfl=ncexp+22

           write(sstfl(ncexp+1:ncexp+22),'(a3,i8.8,a3,i4.4,a4)')        &
     &                      'sst',datest,'.pe',mysub,'.bin'

          else

           ncfl=ncexp+31

           write(sstfl(ncexp+1:ncexp+31),'(a3,i8.8,2(a4,i4.4),a4)')     &
     &                 'sst',datest,'.grp',mygrp,'-sub',mysub,'.bin'

          end if

          open(iosst,iostat=stat,err=130,                               &
     &         file=crsdir(1:nccrs)//sstfl(1:ncfl),                     &
     &         status='old',access='direct',form='unformatted',         &
     &         recl=siz,action='read')

  130     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdsstnxt',8,'cont',1,'              ',14,   &
     &                     iosst,stat)

            end if

            call cpondpe

            call destroy('rdsstnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

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

            call outstd03('rdsstnxt',8,sstfl,ncfl,iosst,1,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Read out the data from the interpolated sea surface temperature
! data file.

          read(iosst,rec=1,iostat=stat,err=140)                         &
     &        ((sstd(i,j),i=0,ni+1),j=0,nj+1)

  140     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdsstnxt',8,'cont',3,'              ',14,   &
     &                     iosst,stat)

            end if

            call cpondpe

            call destroy('rdsstnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.stat-1) then

            call outstd03('rdsstnxt',8,sstfl,108,iosst,3,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Close the interpolated sea surface temperature data file.

          close(iosst,iostat=stat,err=150,status='keep')

  150     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdsstnxt',8,'cont',2,'              ',14,   &
     &                     iosst,stat)

            end if

            call cpondpe

            call destroy('rdsstnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.stat-1) then

            call outstd03('rdsstnxt',8,sstfl,108,iosst,2,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Return the unit number.

          call putunit(iosst)

! -----

! Set the boundary conditions.

          call bcsst(idwbc,idebc,idexbopt,ni,nj,sstd)

! -----

! Set the time tendency of sea surface temperature at current marked
! time and at initial or restart time.

          call setsst(idsstitv,ird,ni,nj,sst,sstd)

! -----

        end do

!! -----

! Calculate the next time to read out.

        nxtsst=nxtsst+sst01

! -----

      end if

!!! -----

      end subroutine s_rdsstnxt

!-----7--------------------------------------------------------------7--

      end module m_rdsstnxt
