!***********************************************************************
      module m_rdaslnxt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read out the data from the interpolated aerosol data file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_adjstqa
      use m_bcasl
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
      use m_setasl
      use m_smooasl

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rdaslnxt, s_rdaslnxt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdaslnxt

        module procedure s_rdaslnxt

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
      subroutine s_rdaslnxt(fpexprim,fpcrsdir,fpncexp,fpnccrs,fpwlngth, &
     &                      fpgsmopt,fpetime,fpaslitv,fasl,ctime,ftime, &
     &                      atinc,ni,nj,nk,nqa,qagpv,qatd,tmp1)
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

      integer, intent(in) :: fpgsmopt
                       ! Formal parameter of unique index of gsmopt

      integer, intent(in) :: fpetime
                       ! Formal parameter of unique index of etime

      integer, intent(in) :: fpaslitv
                       ! Formal parameter of unique index of aslitv

      integer, intent(in) :: fasl
                       ! Descriptor to put into motion
                       ! for aerosol nudging

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer(kind=i8), intent(in) :: ftime
                       ! Model forecast time at 1 step future

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

! Input and output variables

      real, intent(inout) :: qagpv(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Mixing ratio of aerosol data at marked time

      real, intent(inout) :: qatd(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Time tendency of mixing ratio of aerosol data

! Output variable

      real, intent(out) :: atinc
                       ! Lapse of forecast time
                       ! from aerosol data reading

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) aslfl
                       ! Opened file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer wlngth   ! Word length of direct access file

      integer gsmopt   ! Option for GPV data smoothing

      integer in       ! Namelist table index

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in 4th direction

      integer ird      ! Index of count to read out

      integer ncfl     ! Number of character of aerosol data file

      integer ioasl    ! Unit number of aerosol data file

      integer siz      ! Record length of aerosol data file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

      integer dateas   ! Aerosol data file date

      integer asl01    ! int(aslitv + 0.1)

      integer(kind=i8) asl103
                       ! 1000 x int(aslitv + 0.1)

      integer(kind=i8) crtime
                       ! Current forecast time to read out

      real etime       ! Forecast stop time

      real aslitv      ! Time interval of aerosol data

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpetime,etime)
      call getrname(fpaslitv,aslitv)

! -----

! Set the common used variables.

      asl01=int(aslitv+.1e0)

      asl103=1000_i8*int(aslitv+.1e0,i8)

      crtime=asl103*(ctime/asl103)

      atinc=.001e0*real(ctime-crtime)

      if(fasl.eq.0) then

        nxtasl=crtime/1000_i8

      end if

! -----

!!! Open and read out the data from the interpolated aerosol data file
!!! when the current forecast time reaches marked time.

      if(ctime/1000_i8.ge.nxtasl                                        &
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
        call getiname(fpgsmopt,gsmopt)

! -----

!! Open and read out the data from the interpolated aerosol data
!! checking file.

        if(fasl.eq.0) then

! Initialize the character variable.

          if(mype.eq.root) then

            call inichar(aslfl)

          end if

! -----

! Get the unit number.

          if(mype.eq.root) then

            call getunit(ioasl)

          end if

! -----

! Open the interpolated aerosol data checking file.

          if(mype.eq.root) then

            aslfl(1:ncexp)=exprim(1:ncexp)

            write(aslfl(ncexp+1:ncexp+13),'(a13)') 'asl.check.txt'

            open(ioasl,iostat=stat,err=100,                             &
     &           file=crsdir(1:nccrs)//aslfl(1:ncexp+13),               &
     &           status='old',access='sequential',form='formatted',     &
     &           blank='null',position='rewind',action='read')

          else

            stat=0

          end if

  100     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdaslnxt',8,'cont',1,'              ',14,   &
     &                     ioasl,stat)

            end if

            call cpondpe

            call destroy('rdaslnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.root) then

            call outstd03('rdaslnxt',8,aslfl,ncexp+13,ioasl,1,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Read out the data from the interpolated aerosol data checking file.

          if(mype.eq.root) then

            read(ioasl,'(a)',iostat=stat,end=110,err=110)               &
     &          (rcname(in),in=1,ncn)

            read(ioasl,*,iostat=stat,end=110,err=110)                   &
     &          (riname(in),in=1,nin)

            read(ioasl,*,iostat=stat,end=110,err=110)                   &
     &          (rrname(in),in=1,nrn)

          else

            stat=0

          end if

  110     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdaslnxt',8,'cont',3,'              ',14,   &
     &                     ioasl,stat)

            end if

            call cpondpe

            call destroy('rdaslnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.root) then

            call outstd03('rdaslnxt',8,aslfl,108,ioasl,3,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Check the interpolated aerosol data checking file.

          if(mype.eq.root) then

            call chkfile('asl',stat,ncn,nin,nrn,                        &
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

! Close the interpolated aerosol data checking file.

          if(mype.eq.root) then

            close(ioasl,iostat=stat,err=120,status='keep')

          else

            stat=0

          end if

  120     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdaslnxt',8,'cont',2,'              ',14,   &
     &                     ioasl,stat)

            end if

            call cpondpe

            call destroy('rdaslnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.root) then

            call outstd03('rdaslnxt',8,aslfl,108,ioasl,2,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Return the unit number.

          if(mype.eq.root) then

            call putunit(ioasl)

          end if

! -----

        end if

!! -----

!! Open and read out the data from the interpolated aerosol data file.

        siz=(ni+2)*(nj+2)*nk*wlngth

        do ird=2,1,-1

! Initialize the character variable.

          call inichar(aslfl)

! -----

! Get the unit number.

          call getunit(ioasl)

! -----

! Open the interpolated aerosol data file.

          dateas=crtime/1000_i8
          dateas=dateas+abs(ird-2)*asl01

          aslfl(1:ncexp)=exprim(1:ncexp)

          if(ngrp.eq.1) then

            ncfl=ncexp+22

            write(aslfl(ncexp+1:ncexp+22),'(a3,i8.8,a3,i4.4,a4)')       &
     &                       'asl',dateas,'.pe',mysub,'.bin'

          else

            ncfl=ncexp+31

            write(aslfl(ncexp+1:ncexp+31),'(a3,i8.8,2(a4,i4.4),a4)')    &
     &                  'asl',dateas,'.grp',mygrp,'-sub',mysub,'.bin'

          end if

          open(ioasl,iostat=stat,err=130,                               &
     &         file=crsdir(1:nccrs)//aslfl(1:ncfl),                     &
     &         status='old',access='direct',form='unformatted',         &
     &         recl=siz,action='read')

  130     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdaslnxt',8,'cont',1,'              ',14,   &
     &                     ioasl,stat)

            end if

            call cpondpe

            call destroy('rdaslnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.stat-1) then

            if(fpara(1:5).eq.'multi') then

              if(ngrp.eq.1) then

                write(aslfl(ncexp+15:ncexp+18),'(a4)') 'XXXX'

              else

                write(aslfl(ncexp+16:ncexp+19),'(a4)') 'XXXX'
                write(aslfl(ncexp+24:ncexp+27),'(a4)') 'YYYY'

              end if

            end if

            call outstd03('rdaslnxt',8,aslfl,ncfl,ioasl,1,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Read out the data from the interpolated aerosol data file.

          do n=1,nqa(0)

            read(ioasl,rec=n,iostat=stat,err=140)                       &
     &          (((qatd(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk)

          end do

  140     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdaslnxt',8,'cont',3,'              ',14,   &
     &                     ioasl,stat)

            end if

            call cpondpe

            call destroy('rdaslnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.stat-1) then

            call outstd03('rdaslnxt',8,aslfl,108,ioasl,3,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Close the interpolated aerosol data file.

          close(ioasl,iostat=stat,err=150,status='keep')

  150     call chkerr(stat)

          if(stat.lt.0) then

            if(mype.eq.-stat-1) then

              call destroy('rdaslnxt',8,'cont',2,'              ',14,   &
     &                     ioasl,stat)

            end if

            call cpondpe

            call destroy('rdaslnxt',8,'stop',1001,'              ',14,  &
     &                   101,stat)

          end if

          if(mype.eq.stat-1) then

            call outstd03('rdaslnxt',8,aslfl,108,ioasl,2,1,ftime)

          end if

          broot=stat-1

          call chkstd(broot)

! -----

! Return the unit number.

          call putunit(ioasl)

! -----

! Set the boundary conditions.

          call bcasl(idwbc,idebc,idexbopt,ni,nj,nk,nqa,qatd)

! -----

! Smooth the interpolated aerosol data.

          if(gsmopt.eq.1) then

            call smooasl(idgsmcnt,ni,nj,nk,nqa,qatd,tmp1)

          end if

! -----

! Force mixing ratio more than user specified value.

          call adjstqa(ni,nj,nk,nqa,qatd)

! -----

! Set the time tendency of variables at current marked time and the
! variables at initial or restart time.

          call setasl(idaslitv,ird,ni,nj,nk,nqa,qagpv,qatd)

! -----

        end do

!! -----

! Calculate the next time to read out.

        nxtasl=nxtasl+asl01

! -----

      end if

!!! -----

      end subroutine s_rdaslnxt

!-----7--------------------------------------------------------------7--

      end module m_rdaslnxt
