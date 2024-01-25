!***********************************************************************
      module m_outgpv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/10
!     Modification: 1999/06/07, 1999/06/14, 1999/06/21, 1999/07/05,
!                   1999/08/03, 1999/09/30, 1999/11/01, 1999/11/19,
!                   2000/01/05, 2000/01/17, 2000/04/18, 2001/01/15,
!                   2001/02/13, 2001/03/13, 2001/04/15, 2001/05/29,
!                   2001/06/29, 2002/01/07, 2002/04/02, 2002/04/09,
!                   2002/06/18, 2002/07/15, 2002/08/15, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2003/12/12,
!                   2004/01/09, 2004/05/31, 2004/08/20, 2004/08/31,
!                   2004/09/25, 2005/01/14, 2005/02/10, 2006/02/03,
!                   2006/09/21, 2006/12/04, 2007/01/05, 2007/01/20,
!                   2007/07/30, 2007/08/24, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2008/12/11, 2009/01/30, 2009/02/27,
!                   2011/09/22, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read in the data to the interpolated GPV data file.

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

      public :: outgpv, s_outgpv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outgpv

        module procedure s_outgpv

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
      subroutine s_outgpv(fpexprim,fpcrsdir,fpgpvvar,fpncexp,fpnccrs,   &
     &                    fpwlngth,it,nstp0,ctime,ni,nj,nk,ubr,vbr,pbr, &
     &                    ptbr,qvbr,u,v,w,pp,ptp,qv,qc,qr,qi,qs,qg,qh)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexprim
                       ! Formal parameter of unique index of exprim

      integer, intent(in) :: fpcrsdir
                       ! Formal parameter of unique index of crsdir

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

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

      real, intent(in) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(in) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: qc(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio

      real, intent(in) :: qr(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio

      real, intent(in) :: qi(0:ni+1,0:nj+1,1:nk)
                       ! Cloud ice mixing ratio

      real, intent(in) :: qs(0:ni+1,0:nj+1,1:nk)
                       ! Snow mixing ratio

      real, intent(in) :: qg(0:ni+1,0:nj+1,1:nk)
                       ! Graupel mixing ratio

      real, intent(in) :: qh(0:ni+1,0:nj+1,1:nk)
                       ! Hail mixing ratio

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      character(len=108) gpvfl
                       ! Opened file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer wlngth   ! Word length of direct access file

      integer in       ! Namelist table index

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer ncfl     ! Number of character of GPV data file name

      integer extpe    ! Control flag of file name extension

      integer recgpv   ! Current record number of GPV data file

      integer iogpv    ! Unit number of GPV data file

      integer siz      ! Record length of GPV data file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

!-----7--------------------------------------------------------------7--

!!! Open and read in the data to the interpolated GPV data file.

! Initialize the character variables.

      call inichar(exprim)
      call inichar(crsdir)
      call inichar(gpvvar)

! -----

! Get the required namelist variables.

      call getcname(fpexprim,exprim)
      call getcname(fpcrsdir,crsdir)
      call getcname(fpgpvvar,gpvvar)
      call getiname(fpncexp,ncexp)
      call getiname(fpnccrs,nccrs)
      call getiname(fpwlngth,wlngth)

! -----

!! Open and read in the data to the interpolated GPV data checking file.

      if(it.eq.nstp0) then

! Initialize the character variable.

        if(mype.eq.root) then

          call inichar(gpvfl)

        end if

! -----

! Get the unit number.

        if(mype.eq.root) then

          call getunit(iogpv)

        end if

! -----

! Open the interpolated GPV data checking file.

        if(mype.eq.root) then

          gpvfl(1:ncexp)=exprim(1:ncexp)

          write(gpvfl(ncexp+1:ncexp+13),'(a13)') 'gpv.check.txt'

          open(iogpv,iostat=stat,err=100,                               &
     &         file=crsdir(1:nccrs)//gpvfl(1:ncexp+13),                 &
     &         status='new',access='sequential',form='formatted',       &
     &         blank='null',action='write')

  100     if(stat.eq.0) then

            ncfl=ncexp+13

          else

            ncfl=ncexp+18

            write(gpvfl(ncexp+14:ncexp+18),'(a5)') '.swap'

            open(iogpv,iostat=stat,err=110,                             &
     &           file=crsdir(1:nccrs)//gpvfl(1:ncexp+18),               &
     &           status='new',access='sequential',form='formatted',     &
     &           blank='null',action='write')

          end if

        else

          stat=0

        end if

  110   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('outgpv  ',6,'cont',1,'              ',14,     &
     &                   iogpv,stat)

          end if

          call cpondpe

          call destroy('outgpv  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('outgpv  ',6,gpvfl,ncfl,iogpv,1,1,ctime)

        end if

        broot=stat-1

        call chkstd(broot)

! -----

! Read in the data to the interpolated GPV data checking file.

        if(mype.eq.root) then

          write(iogpv,'(a)',iostat=stat,err=120)                        &
     &         (cname(in),in=1,ncn)

          write(iogpv,*,iostat=stat,err=120)                            &
     &         (iname(in),in=1,nin)

          write(iogpv,*,iostat=stat,err=120)                            &
     &         (rname(in),in=1,nrn)

        else

          stat=0

        end if

  120   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('outgpv  ',6,'cont',3,'              ',14,     &
     &                   iogpv,stat)

          end if

          call cpondpe

          call destroy('outgpv  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('outgpv  ',6,gpvfl,108,iogpv,4,1,ctime)

        end if

        broot=stat-1

        call chkstd(broot)

! -----

! Close the interpolated GPV data checking file.

        if(mype.eq.root) then

          close(iogpv,iostat=stat,err=130,status='keep')

        else

          stat=0

        end if

  130   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('outgpv  ',6,'cont',2,'              ',14,     &
     &                   iogpv,stat)

          end if

          call cpondpe

          call destroy('outgpv  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('outgpv  ',6,gpvfl,108,iogpv,2,1,ctime)

        end if

        broot=stat-1

        call chkstd(broot)

! -----

! Return the unit number.

        if(mype.eq.root) then

          call putunit(iogpv)

        end if

! -----

      end if

!! -----

!! Open and read in the data to the interpolated GPV data file.

! Initialize the character variable.

      call inichar(gpvfl)

! -----

! Get the unit number.

      call getunit(iogpv)

! -----

! Open the interpolated GPV data file.

      siz=(ni+2)*(nj+2)*nk*wlngth

      gpvfl(1:ncexp)=exprim(1:ncexp)

      if(ngrp.eq.1) then

        ncfl=ncexp+22

        write(gpvfl(ncexp+1:ncexp+22),'(a3,i8.8,a3,i4.4,a4)')           &
     &            'gpv',ctime/1000_i8,'.pe',mysub,'.bin'

      else

        ncfl=ncexp+31

        write(gpvfl(ncexp+1:ncexp+31),'(a3,i8.8,2(a4,i4.4),a4)')        &
     &            'gpv',ctime/1000_i8,'.grp',mygrp,'-sub',mysub,'.bin'

      end if

      open(iogpv,iostat=stat,err=140,                                   &
     &     file=crsdir(1:nccrs)//gpvfl(1:ncfl),                         &
     &     status='new',access='direct',form='unformatted',             &
     &     recl=siz,action='write')

  140 if(stat.eq.0) then

        extpe=0

      else

        extpe=1

        ncfl=ncfl+5

        write(gpvfl(ncfl-4:ncfl),'(a5)') '.swap'

        open(iogpv,iostat=stat,err=150,                                 &
     &       file=crsdir(1:nccrs)//gpvfl(1:ncfl),                       &
     &       status='new',access='direct',form='unformatted',           &
     &       recl=siz,action='write')

      end if

  150 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outgpv  ',6,'cont',1,'              ',14,iogpv, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outgpv  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      call chkopen(extpe,iogpv)

      if(mype.eq.stat-1) then

        if(fpara(1:5).eq.'multi') then

          if(ngrp.eq.1) then

            write(gpvfl(ncexp+15:ncexp+18),'(a4)') 'XXXX'

          else

            write(gpvfl(ncexp+16:ncexp+19),'(a4)') 'XXXX'
            write(gpvfl(ncexp+24:ncexp+27),'(a4)') 'YYYY'

          end if

        end if

        call outstd03('outgpv  ',6,gpvfl,ncfl,iogpv,1,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read in the data to the interpolated GPV data file.

      if(ctime.eq.0_i8) then

        if(gpvvar(2:2).eq.'o') then

          recgpv=9

          write(iogpv,rec=1,iostat=stat,err=160)                        &
     &         (((ubr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=2,iostat=stat,err=160)                        &
     &         (((vbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=3,iostat=stat,err=160)                        &
     &         (((pbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=4,iostat=stat,err=160)                        &
     &         (((ptbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=5,iostat=stat,err=160)                        &
     &         (((qvbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=6,iostat=stat,err=160)                        &
     &         (((u(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=7,iostat=stat,err=160)                        &
     &         (((v(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=8,iostat=stat,err=160)                        &
     &         (((pp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=9,iostat=stat,err=160)                        &
     &         (((ptp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        else

          recgpv=8

          write(iogpv,rec=1,iostat=stat,err=160)                        &
     &         (((ubr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=2,iostat=stat,err=160)                        &
     &         (((vbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=3,iostat=stat,err=160)                        &
     &         (((pbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=4,iostat=stat,err=160)                        &
     &         (((ptbr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=5,iostat=stat,err=160)                        &
     &         (((u(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=6,iostat=stat,err=160)                        &
     &         (((v(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=7,iostat=stat,err=160)                        &
     &         (((pp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

          write(iogpv,rec=8,iostat=stat,err=160)                        &
     &         (((ptp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        end if

      else

        recgpv=4

        write(iogpv,rec=1,iostat=stat,err=160)                          &
     &       (((u(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        write(iogpv,rec=2,iostat=stat,err=160)                          &
     &       (((v(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        write(iogpv,rec=3,iostat=stat,err=160)                          &
     &       (((pp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

        write(iogpv,rec=4,iostat=stat,err=160)                          &
     &       (((ptp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(gpvvar(1:1).eq.'o') then
        recgpv=recgpv+1

        write(iogpv,rec=recgpv,iostat=stat,err=160)                     &
     &       (((w(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(gpvvar(2:2).eq.'o') then
        recgpv=recgpv+1

        write(iogpv,rec=recgpv,iostat=stat,err=160)                     &
     &       (((qv(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(gpvvar(3:3).eq.'o') then
        recgpv=recgpv+1

        write(iogpv,rec=recgpv,iostat=stat,err=160)                     &
     &       (((qc(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(gpvvar(4:4).eq.'o') then
        recgpv=recgpv+1

        write(iogpv,rec=recgpv,iostat=stat,err=160)                     &
     &       (((qr(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(gpvvar(5:5).eq.'o') then
        recgpv=recgpv+1

        write(iogpv,rec=recgpv,iostat=stat,err=160)                     &
     &       (((qi(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(gpvvar(6:6).eq.'o') then
        recgpv=recgpv+1

        write(iogpv,rec=recgpv,iostat=stat,err=160)                     &
     &       (((qs(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(gpvvar(7:7).eq.'o') then
        recgpv=recgpv+1

        write(iogpv,rec=recgpv,iostat=stat,err=160)                     &
     &       (((qg(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(gpvvar(8:8).eq.'o') then
        recgpv=recgpv+1

        write(iogpv,rec=recgpv,iostat=stat,err=160)                     &
     &       (((qh(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

  160 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outgpv  ',6,'cont',3,'              ',14,iogpv, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outgpv  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('outgpv  ',6,gpvfl,108,iogpv,4,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the interpolated GPV data file.

      close(iogpv,iostat=stat,err=170,status='keep')

  170 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outgpv  ',6,'cont',2,'              ',14,iogpv, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outgpv  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('outgpv  ',6,gpvfl,108,iogpv,2,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      call putunit(iogpv)

! -----

!! -----

!!! -----

      end subroutine s_outgpv

!-----7--------------------------------------------------------------7--

      end module m_outgpv
