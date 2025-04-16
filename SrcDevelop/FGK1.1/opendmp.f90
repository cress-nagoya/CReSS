!***********************************************************************
      module m_opendmp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/03/25, 1999/04/06, 1999/05/10, 1999/06/14,
!                   1999/06/21, 1999/08/23, 1999/11/19, 2000/01/05,
!                   2000/01/17, 2000/04/18, 2000/07/05, 2001/02/13,
!                   2001/04/15, 2001/05/29, 2001/11/20, 2002/04/02,
!                   2002/06/18, 2002/07/15, 2003/04/30, 2003/05/19,
!                   2003/06/27, 2004/05/31, 2004/06/10, 2004/08/01,
!                   2004/08/20, 2004/09/01, 2004/09/25, 2005/01/14,
!                   2005/02/10, 2006/09/21, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/01/31, 2007/07/30, 2007/08/24,
!                   2007/10/19, 2008/01/11, 2008/04/17, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2008/10/10, 2008/12/11,
!                   2009/02/27, 2011/09/22, 2013/01/28, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open the dumped file when the current forecast time reaches marked
!     time.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chkopen
      use m_chkstd
      use m_comdmp
      use m_comindx
      use m_comkind
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_getiname
      use m_getrname
      use m_getunit
      use m_inichar
      use m_outstd03
      use m_outstd07

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: opendmp, s_opendmp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface opendmp

        module procedure s_opendmp

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic int
      intrinsic mod
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_opendmp(fpexprim,fpcrsdir,fpncexp,fpnccrs,           &
     &                     fpwlngth,fpdmpfmt,fpdmplev,fpdmpmon,         &
     &                     fpetime,fpdmpitv,fpmonitv,fpdz,ctime,        &
     &                     ni,nj,nk,zsth,z1d)
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

      integer, intent(in) :: fpdmpfmt
                       ! Formal parameter of unique index of dmpfmt

      integer, intent(in) :: fpdmplev
                       ! Formal parameter of unique index of dmplev

      integer, intent(in) :: fpdmpmon
                       ! Formal parameter of unique index of dmpmon

      integer, intent(in) :: fpetime
                       ! Formal parameter of unique index of etime

      integer, intent(in) :: fpdmpitv
                       ! Formal parameter of unique index of dmpitv

      integer, intent(in) :: fpmonitv
                       ! Formal parameter of unique index of monitv

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zsth(1:nk)
                       ! 1 dimensional stretched z coordinates

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer wlngth   ! Word length of direct access file

      integer dmpfmt   ! Option for dumped file format
      integer dmplev   ! Option for z coordinates of dumped variables
      integer dmpmon   ! Option for monitor variables output

      integer extpe    ! Control flag of file name extension

      integer siz      ! Record length of dumped file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

      real etime       ! Forecast stop time

      real dmpitv      ! Time interval of dumped file
      real monitv      ! Time interval of dumped file
                       ! for monitor variables

      real dz          ! Grid distance in z direction

      real, intent(inout) :: z1d(1:nk)
                       ! Output z coordinates at scalar point

! Internal private variable

      integer k        ! Array index in z drection

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpdmpmon,dmpmon)
      call getrname(fpetime,etime)
      call getrname(fpdmpitv,dmpitv)
      call getrname(fpmonitv,monitv)

! -----

! Set the control flag fdmp and fmon to open the dumped file to the
! module m_comdmp.

      if(mod(ctime,1000_i8*int(dmpitv+.1e0,i8)).eq.0_i8) then

        write(fdmp(1:3),'(a3)') 'act'

      else

        write(fdmp(1:3),'(a3)') 'off'

      end if

      if(dmpmon.eq.0) then

        if(mod(ctime,1000_i8*int(dmpitv+.1e0,i8)).eq.0_i8) then

          write(fmon(1:3),'(a3)') 'act'

        else

          write(fmon(1:3),'(a3)') 'off'

        end if

      else

        if(mod(ctime,1000_i8*int(monitv+.1e0,i8)).eq.0_i8) then

          write(fmon(1:3),'(a3)') 'act'

        else

          write(fmon(1:3),'(a3)') 'off'

        end if

      end if

      if(ctime.eq.1000_i8*int(etime+.1e0,i8)) then

        write(fdmp(1:3),'(a3)') 'act'
        write(fmon(1:3),'(a3)') 'act'

      end if

! -----

!! Open the dumped file when the current forecast time reaches marked
!! time.

      if(fdmp(1:3).eq.'act'.or.fmon(1:3).eq.'act') then

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
        call getiname(fpdmpfmt,dmpfmt)
        call getiname(fpdmplev,dmplev)
        call getrname(fpdz,dz)

! -----

! Calculate the constant output level.

!$omp parallel default(shared)

        if(fdmp(1:3).eq.'act') then

          if(mod(dmplev,10).eq.2) then

!$omp do schedule(runtime) private(k)

            do k=2,nk-2
              z1d(k)=dz*(real(k)-1.5e0)
            end do

!$omp end do

          else if(mod(dmplev,10).eq.3) then

!$omp do schedule(runtime) private(k)

            do k=2,nk-2
              z1d(k)=.5e0*(zsth(k)+zsth(k+1))
            end do

!$omp end do

          end if

        end if

!$omp end parallel

! -----

! Initialize the character variables.

        if(mype.eq.root) then

          if(fdmp(1:3).eq.'act') then

            call inichar(fl3c)

          end if

          if(fmon(1:3).eq.'act') then

            if(dmpmon.eq.1) then

              call inichar(fl2c)

            end if

          end if

        end if

        if(fdmp(1:3).eq.'act') then

          call inichar(fl3d)

        end if

        if(fmon(1:3).eq.'act') then

          if(dmpmon.eq.1) then

            call inichar(fl2d)

          end if

        end if

! -----

! Get the unit numbers.

        if(mype.eq.root) then

          if(fdmp(1:3).eq.'act') then

            call getunit(io3c)

          end if

          if(fmon(1:3).eq.'act') then

            if(dmpmon.eq.1) then

              call getunit(io2c)

            end if

          end if

        end if

        if(fdmp(1:3).eq.'act') then

          call getunit(io3d)

        end if

        if(fmon(1:3).eq.'act') then

          if(dmpmon.eq.1) then

            call getunit(io2d)

          end if

        end if

! -----

! Open the dumped data checking file.

        if(fdmp(1:3).eq.'act') then

          if(mype.eq.root) then

            fl3c(1:ncexp)=exprim(1:ncexp)

            write(fl3c(ncexp+1:ncexp+13),'(a13)') 'dmp.check.txt'

            if(ctime.eq.0_i8) then

              open(io3c,iostat=stat,err=100,                            &
     &             file=crsdir(1:nccrs)//fl3c(1:ncexp+13),              &
     &             status='new',access='sequential',form='formatted',   &
     &             blank='null',action='write')

  100         if(stat.eq.0) then

                nc3c=ncexp+13

              else

                nc3c=ncexp+18

                write(fl3c(ncexp+14:ncexp+18),'(a5)') '.swap'

                open(io3c,iostat=stat,err=110,                          &
     &               file=crsdir(1:nccrs)//fl3c(1:ncexp+18),            &
     &               status='new',access='sequential',form='formatted', &
     &               blank='null',action='write')

              end if

            else

              nc3c=ncexp+13

              open(io3c,iostat=stat,err=110,                            &
     &            file=crsdir(1:nccrs)//fl3c(1:nc3c),                   &
     &            status='replace',access='sequential',form='formatted',&
     &            blank='null',action='write')

            end if

          else

            stat=0

          end if

  110     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('opendmp ',7,'cont',1,'              ',14,   &
     &                     io3c,stat)

            end if

            call cpondpe

            call destroy('opendmp ',7,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.root) then

            call outstd03('opendmp ',7,fl3c,nc3c,io3c,1,1,ctime)

          end if

          broot=stat-1

          call chkstd(broot)

        end if

! -----

! Open the dumped data checking file for monitor variables.

        if(fmon(1:3).eq.'act') then

          if(dmpmon.eq.1) then

            if(mype.eq.root) then

              fl2c(1:ncexp)=exprim(1:ncexp)

              write(fl2c(ncexp+1:ncexp+13),'(a13)') 'mon.check.txt'

              if(ctime.eq.0_i8) then

                open(io2c,iostat=stat,err=120,                          &
     &               file=crsdir(1:nccrs)//fl2c(1:ncexp+13),            &
     &               status='new',access='sequential',form='formatted', &
     &               blank='null',action='write')

  120           if(stat.eq.0) then

                  nc2c=ncexp+13

                else

                  nc2c=ncexp+18

                  write(fl2c(ncexp+14:ncexp+18),'(a5)') '.swap'

                  open(io2c,iostat=stat,err=130,                        &
     &                file=crsdir(1:nccrs)//fl2c(1:ncexp+18),           &
     &                status='new',access='sequential',form='formatted',&
     &                blank='null',action='write')

                end if

              else

                nc2c=ncexp+13

                open(io2c,iostat=stat,err=130,                          &
     &            file=crsdir(1:nccrs)//fl2c(1:nc2c),                   &
     &            status='replace',access='sequential',form='formatted',&
     &            blank='null',action='write')

              end if

            else

              stat=0

            end if

  130       call chkerr(stat)

            if(stat.lt.0) then

              if(mype.eq.-stat-1) then

                call destroy('opendmp ',7,'cont',1,'              ',14, &
     &                       io2c,stat)

              end if

              call cpondpe

              call destroy('opendmp ',7,'stop',1001,'              ',14,&
     &                     101,stat)

            end if

            if(mype.eq.root) then

              call outstd03('opendmp ',7,fl2c,nc2c,io2c,1,1,ctime)

            end if

            broot=stat-1

            call chkstd(broot)

          end if

        end if

! -----

! Open the dumped file.

        if(fdmp(1:3).eq.'act') then

          fl3d(1:ncexp)=exprim(1:ncexp)

          if(ngrp.eq.1) then

            nc3d=ncexp+22

            write(fl3d(ncexp+1:ncexp+18),'(a3,i8.8,a3,i4.4)')           &
     &               'dmp',ctime/1000_i8,'.pe',mysub

          else

            nc3d=ncexp+31

            write(fl3d(ncexp+1:ncexp+27),'(a3,i8.8,2(a4,i4.4))')        &
     &               'dmp',ctime/1000_i8,'.grp',mygrp,'-sub',mysub

          end if

          if(dmpfmt.eq.1) then

            write(fl3d(nc3d-3:nc3d),'(a4)') '.txt'

            open(io3d,iostat=stat,err=140,                              &
     &           file=crsdir(1:nccrs)//fl3d(1:nc3d),                    &
     &           status='new',access='sequential',form='formatted',     &
     &           blank='null',action='write')

  140       if(stat.eq.0) then

              extpe=0

            else

              extpe=1

              nc3d=nc3d+5

              write(fl3d(nc3d-4:nc3d),'(a5)') '.swap'

              open(io3d,iostat=stat,err=150,                            &
     &             file=crsdir(1:nccrs)//fl3d(1:nc3d),                  &
     &             status='new',access='sequential',form='formatted',   &
     &             blank='null',action='write')

            end if

          else if(dmpfmt.eq.2) then

            siz=(ni-3)*(nj-3)*wlngth

            write(fl3d(nc3d-3:nc3d),'(a4)') '.bin'

            open(io3d,iostat=stat,err=160,                              &
     &           file=crsdir(1:nccrs)//fl3d(1:nc3d),                    &
     &           status='new',access='direct',form='unformatted',       &
     &           recl=siz,action='write')

  160       if(stat.eq.0) then

              extpe=0

            else

              extpe=1

              nc3d=nc3d+5

              write(fl3d(nc3d-4:nc3d),'(a5)') '.swap'

              open(io3d,iostat=stat,err=150,                            &
     &             file=crsdir(1:nccrs)//fl3d(1:nc3d),                  &
     &             status='new',access='direct',form='unformatted',     &
     &             recl=siz,action='write')

            end if

          end if

  150     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('opendmp ',7,'cont',1,'              ',14,   &
     &                     io3d,stat)

            end if

            call cpondpe

            call destroy('opendmp ',7,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          call chkopen(extpe,io3d)

          if(mype.eq.stat-1) then

            if(fpara(1:5).eq.'multi') then

              if(ngrp.eq.1) then

                write(fl3d(ncexp+15:ncexp+18),'(a4)') 'XXXX'

              else

                write(fl3d(ncexp+16:ncexp+19),'(a4)') 'XXXX'
                write(fl3d(ncexp+24:ncexp+27),'(a4)') 'YYYY'

              end if

            end if

            call outstd03('opendmp ',7,fl3d,nc3d,io3d,1,1,ctime)

          end if

          broot=stat-1

          call chkstd(broot)

        end if

! -----

! Open the dumped file for monitor variables.

        if(fmon(1:3).eq.'act') then

          if(dmpmon.eq.1) then

            fl2d(1:ncexp)=exprim(1:ncexp)

            if(ngrp.eq.1) then

              nc2d=ncexp+22

              write(fl2d(ncexp+1:ncexp+18),'(a3,i8.8,a3,i4.4)')         &
     &                 'mon',ctime/1000_i8,'.pe',mysub

            else

              nc2d=ncexp+31

              write(fl2d(ncexp+1:ncexp+27),'(a3,i8.8,2(a4,i4.4))')      &
     &                 'mon',ctime/1000_i8,'.grp',mygrp,'-sub',mysub

            end if

            if(dmpfmt.eq.1) then

              write(fl2d(nc2d-3:nc2d),'(a4)') '.txt'

              open(io2d,iostat=stat,err=170,                            &
     &             file=crsdir(1:nccrs)//fl2d(1:nc2d),                  &
     &             status='new',access='sequential',form='formatted',   &
     &             blank='null',action='write')

  170         if(stat.eq.0) then

                extpe=0

              else

                extpe=1

                nc2d=nc2d+5

                write(fl2d(nc2d-4:nc2d),'(a5)') '.swap'

                open(io2d,iostat=stat,err=180,                          &
     &               file=crsdir(1:nccrs)//fl2d(1:nc2d),                &
     &               status='new',access='sequential',form='formatted', &
     &               blank='null',action='write')

              end if

            else if(dmpfmt.eq.2) then

              siz=(ni-3)*(nj-3)*wlngth

              write(fl2d(nc2d-3:nc2d),'(a4)') '.bin'

              open(io2d,iostat=stat,err=190,                            &
     &             file=crsdir(1:nccrs)//fl2d(1:nc2d),                  &
     &             status='new',access='direct',form='unformatted',     &
     &             recl=siz,action='write')

  190         if(stat.eq.0) then

                extpe=0

              else

                extpe=1

                nc2d=nc2d+5

                write(fl2d(nc2d-4:nc2d),'(a5)') '.swap'

                open(io2d,iostat=stat,err=180,                          &
     &               file=crsdir(1:nccrs)//fl2d(1:nc2d),                &
     &               status='new',access='direct',form='unformatted',   &
     &               recl=siz,action='write')

              end if

            end if

  180       call chkerr(stat)

            if(stat.lt.0) then

              if(mype.eq.-stat-1) then

                call destroy('opendmp ',7,'cont',1,'              ',14, &
     &                       io2d,stat)

              end if

              call cpondpe

              call destroy('opendmp ',7,'stop',1001,'              ',14,&
     &                     101,stat)

            end if

            call chkopen(extpe,io2d)

            if(mype.eq.stat-1) then

              if(fpara(1:5).eq.'multi') then

                if(ngrp.eq.1) then

                  write(fl2d(ncexp+15:ncexp+18),'(a4)') 'XXXX'

                else

                  write(fl2d(ncexp+16:ncexp+19),'(a4)') 'XXXX'
                  write(fl2d(ncexp+24:ncexp+27),'(a4)') 'YYYY'

                end if

              end if

              call outstd03('opendmp ',7,fl2d,nc2d,io2d,1,1,ctime)

            end if

            broot=stat-1

            call chkstd(broot)

          end if

        end if

! -----

! Read in the messages to standard i/o.

        if(fdmp(1:3).eq.'act') then

          if(mype.eq.root) then

            if(mod(dmplev,10).eq.2.or.mod(dmplev,10).eq.3) then

              call outstd07(iddmplev,nk,z1d)

            end if

          end if

          broot=stat-1

          call chkstd(broot)

        end if

! -----

! Reset the counter and the record number of the dumped file to 0.

        if(fdmp(1:3).eq.'act') then
          cnt3d=0
          rec3d=0
        end if

        if(fmon(1:3).eq.'act') then

          if(dmpmon.eq.1) then
            cnt2d=0
            rec2d=0
          end if

        end if

! -----

      end if

!! -----

      end subroutine s_opendmp

!-----7--------------------------------------------------------------7--

      end module m_opendmp
