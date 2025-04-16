!***********************************************************************
      module m_outasl
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/09/22
!     Modification: 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read in the data to the interpolated aerosol data file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chkopen
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

      public :: outasl, s_outasl

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outasl

        module procedure s_outasl

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
      subroutine s_outasl(fpexprim,fpcrsdir,fpncexp,fpnccrs,fpwlngth,   &
     &                    it,nstp0,ctime,ni,nj,nk,nqa,qasl)
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

      integer(kind=i8), intent(in) :: it
                       ! Index of main do loop in upper procedure

      integer(kind=i8), intent(in) :: nstp0
                       ! Start index of main do loop in upper procedure

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqa(0:4)
                       ! Number of types of aerosol

      real, intent(in) :: qasl(0:ni+1,0:nj+1,1:nk,1:nqa(0))
                       ! Aerosol mixing ratio

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

      integer in       ! Namelist table index

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer n        ! Array index in 4th direction

      integer ncfl     ! Number of character of aerosol data file name

      integer extpe    ! Control flag of file name extension

      integer ioasl    ! Unit number of aerosol data file

      integer siz      ! Record length of aerosol data file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

!-----7--------------------------------------------------------------7--

!!! Open and read in the data to the interpolated aerosol data file.

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

!! Open and read in the data to the interpolated aerosol data checking
!! file.

      if(it.eq.nstp0) then

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

          open(ioasl,iostat=stat,err=100,                               &
     &         file=crsdir(1:nccrs)//aslfl(1:ncexp+13),                 &
     &         status='new',access='sequential',form='formatted',       &
     &         blank='null',action='write')

  100     if(stat.eq.0) then

            ncfl=ncexp+13

          else

            ncfl=ncexp+18

            write(aslfl(ncexp+14:ncexp+18),'(a5)') '.swap'

            open(ioasl,iostat=stat,err=110,                             &
     &           file=crsdir(1:nccrs)//aslfl(1:ncexp+18),               &
     &           status='new',access='sequential',form='formatted',     &
     &           blank='null',action='write')

          end if

        else

          stat=0

        end if

  110   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('outasl  ',6,'cont',1,'              ',14,     &
     &                   ioasl,stat)

          end if

          call cpondpe

          call destroy('outasl  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('outasl  ',6,aslfl,ncfl,ioasl,1,1,ctime)

        end if

        broot=stat-1

        call chkstd(broot)

! -----

! Read in the data to the interpolated aerosol data checking file.

        if(mype.eq.root) then

          write(ioasl,'(a)',iostat=stat,err=120)                        &
     &         (cname(in),in=1,ncn)

          write(ioasl,*,iostat=stat,err=120)                            &
     &         (iname(in),in=1,nin)

          write(ioasl,*,iostat=stat,err=120)                            &
     &         (rname(in),in=1,nrn)

        else

          stat=0

        end if

  120   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('outasl  ',6,'cont',3,'              ',14,     &
     &                   ioasl,stat)

          end if

          call cpondpe

          call destroy('outasl  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('outasl  ',6,aslfl,108,ioasl,4,1,ctime)

        end if

        broot=stat-1

        call chkstd(broot)

! -----

! Close the interpolated aerosol data checking file.

        if(mype.eq.root) then

          close(ioasl,iostat=stat,err=130,status='keep')

        else

          stat=0

        end if

  130   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('outasl  ',6,'cont',2,'              ',14,     &
     &                   ioasl,stat)

          end if

          call cpondpe

          call destroy('outasl  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('outasl  ',6,aslfl,108,ioasl,2,1,ctime)

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

!! Open and read in the data to the interpolated aerosol data file.

! Initialize the character variable.

      call inichar(aslfl)

! -----

! Get the unit number.

      call getunit(ioasl)

! -----

! Open the interpolated aerosol data file.

      siz=(ni+2)*(nj+2)*nk*wlngth

      aslfl(1:ncexp)=exprim(1:ncexp)

      if(ngrp.eq.1) then

        ncfl=ncexp+22

        write(aslfl(ncexp+1:ncexp+22),'(a3,i8.8,a3,i4.4,a4)')           &
     &            'asl',ctime/1000_i8,'.pe',mysub,'.bin'

      else

        ncfl=ncexp+31

        write(aslfl(ncexp+1:ncexp+31),'(a3,i8.8,2(a4,i4.4),a4)')        &
     &            'asl',ctime/1000_i8,'.grp',mygrp,'-sub',mysub,'.bin'

      end if

      open(ioasl,iostat=stat,err=140,                                   &
     &     file=crsdir(1:nccrs)//aslfl(1:ncfl),                         &
     &     status='new',access='direct',form='unformatted',             &
     &     recl=siz,action='write')

  140 if(stat.eq.0) then

        extpe=0

      else

        extpe=1

        ncfl=ncfl+5

        write(aslfl(ncfl-4:ncfl),'(a5)') '.swap'

        open(ioasl,iostat=stat,err=150,                                 &
     &       file=crsdir(1:nccrs)//aslfl(1:ncfl),                       &
     &       status='new',access='direct',form='unformatted',           &
     &       recl=siz,action='write')

      end if

  150 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outasl  ',6,'cont',1,'              ',14,ioasl, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outasl  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      call chkopen(extpe,ioasl)

      if(mype.eq.stat-1) then

        if(fpara(1:5).eq.'multi') then

          if(ngrp.eq.1) then

            write(aslfl(ncexp+15:ncexp+18),'(a4)') 'XXXX'

          else

            write(aslfl(ncexp+16:ncexp+19),'(a4)') 'XXXX'
            write(aslfl(ncexp+24:ncexp+27),'(a4)') 'YYYY'

          end if

        end if

        call outstd03('outasl  ',6,aslfl,ncfl,ioasl,1,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read in the data to the interpolated aerosol data file.

      do n=1,nqa(0)

        write(ioasl,rec=n,iostat=stat,err=160)                          &
     &       (((qasl(i,j,k,n),i=0,ni+1),j=0,nj+1),k=1,nk)

      end do

  160 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outasl  ',6,'cont',3,'              ',14,ioasl, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outasl  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('outasl  ',6,aslfl,108,ioasl,4,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the interpolated aerosol data file.

      close(ioasl,iostat=stat,err=170,status='keep')

  170 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outasl  ',6,'cont',2,'              ',14,ioasl, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outasl  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('outasl  ',6,aslfl,108,ioasl,2,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      call putunit(ioasl)

! -----

!! -----

!!! -----

      end subroutine s_outasl

!-----7--------------------------------------------------------------7--

      end module m_outasl
