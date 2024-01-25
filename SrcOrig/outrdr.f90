!***********************************************************************
      module m_outrdr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2002/12/02, 2003/03/28, 2003/04/30, 2003/05/19,
!                   2003/06/27, 2004/01/09, 2004/05/31, 2004/08/20,
!                   2004/09/01, 2004/09/25, 2005/01/14, 2005/02/10,
!                   2006/09/21, 2006/12/04, 2007/01/05, 2007/01/20,
!                   2007/07/30, 2007/08/24, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/01/30, 2009/02/27, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read in the data to the interpolated radar data file.

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

      public :: outrdr, s_outrdr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outrdr

        module procedure s_outrdr

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
      subroutine s_outrdr(fpexprim,fpcrsdir,fprdrvar,fpncexp,fpnccrs,   &
     &                    fpwlngth,it,nstp0,ctime,ni,nj,nk,u,v,w,qp)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexprim
                       ! Formal parameter of unique index of exprim

      integer, intent(in) :: fpcrsdir
                       ! Formal parameter of unique index of crsdir

      integer, intent(in) :: fprdrvar
                       ! Formal parameter of unique index of rdrvar

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

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity

      real, intent(in) :: qp(0:ni+1,0:nj+1,1:nk)
                       ! Precipitation mixing ratio

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) rdrvar
                       ! Control flag of input radar data variables

      character(len=108) rdrfl
                       ! Opened file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer wlngth   ! Word length of direct access file

      integer in       ! Namelist table index

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer ncfl     ! Number of character of radar data file name

      integer extpe    ! Control flag of file name extension

      integer recrdr   ! Current record number of radar data file

      integer iordr    ! Unit number of radar data file

      integer siz      ! Record length of radar data file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

!-----7--------------------------------------------------------------7--

!!! Open and read in the data to the interpolated radar data file.

! Initialize the character variables.

      call inichar(exprim)
      call inichar(crsdir)
      call inichar(rdrvar)

! -----

! Get the required namelist variables.

      call getcname(fpexprim,exprim)
      call getcname(fpcrsdir,crsdir)
      call getcname(fprdrvar,rdrvar)
      call getiname(fpncexp,ncexp)
      call getiname(fpnccrs,nccrs)
      call getiname(fpwlngth,wlngth)

! -----

!! Open and read in the data to the interpolated radar data checking
!! file.

      if(it.eq.nstp0) then

! Initialize the character variable.

        if(mype.eq.root) then

          call inichar(rdrfl)

        end if

! -----

! Get the unit number.

        if(mype.eq.root) then

          call getunit(iordr)

        end if

! -----

! Open the interpolated radar data checking file.

        if(mype.eq.root) then

          rdrfl(1:ncexp)=exprim(1:ncexp)

          write(rdrfl(ncexp+1:ncexp+13),'(a13)') 'rdr.check.txt'

          open(iordr,iostat=stat,err=100,                               &
     &         file=crsdir(1:nccrs)//rdrfl(1:ncexp+13),                 &
     &         status='new',access='sequential',form='formatted',       &
     &         blank='null',action='write')

  100     if(stat.eq.0) then

            ncfl=ncexp+13

          else

            ncfl=ncexp+18

            write(rdrfl(ncexp+14:ncexp+18),'(a5)') '.swap'

            open(iordr,iostat=stat,err=110,                             &
     &           file=crsdir(1:nccrs)//rdrfl(1:ncexp+18),               &
     &           status='new',access='sequential',form='formatted',     &
     &           blank='null',action='write')

          end if

        else

          stat=0

        end if

  110   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('outrdr  ',6,'cont',1,'              ',14,     &
     &                   iordr,stat)

          end if

          call cpondpe

          call destroy('outrdr  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('outrdr  ',6,rdrfl,ncfl,iordr,1,1,ctime)

        end if

        broot=stat-1

        call chkstd(broot)

! -----

! Read in the data to the interpolated radar data checking file.

        if(mype.eq.root) then

          write(iordr,'(a)',iostat=stat,err=120)                        &
     &         (cname(in),in=1,ncn)

          write(iordr,*,iostat=stat,err=120)                            &
     &         (iname(in),in=1,nin)

          write(iordr,*,iostat=stat,err=120)                            &
     &         (rname(in),in=1,nrn)

        else

          stat=0

        end if

  120   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('outrdr  ',6,'cont',3,'              ',14,     &
     &                   iordr,stat)

          end if

          call cpondpe

          call destroy('outrdr  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('outrdr  ',6,rdrfl,108,iordr,4,1,ctime)

        end if

        broot=stat-1

        call chkstd(broot)

! -----

! Close the interpolated radar data checking file.

        if(mype.eq.root) then

          close(iordr,iostat=stat,err=130,status='keep')

        else

          stat=0

        end if

  130   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('outrdr  ',6,'cont',2,'              ',14,     &
     &                   iordr,stat)

          end if

          call cpondpe

          call destroy('outrdr  ',6,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('outrdr  ',6,rdrfl,108,iordr,2,1,ctime)

        end if

        broot=stat-1

        call chkstd(broot)

! -----

! Return the unit number.

        if(mype.eq.root) then

          call putunit(iordr)

        end if

! -----

      end if

!! -----

!! Open and read in the data to the interpolated radar data file.

! Initialize the character variable.

      call inichar(rdrfl)

! -----

! Get the unit number.

      call getunit(iordr)

! -----

! Open the interpolated radar data file.

      siz=(ni+2)*(nj+2)*nk*wlngth

      rdrfl(1:ncexp)=exprim(1:ncexp)

      if(ngrp.eq.1) then

        ncfl=ncexp+22

        write(rdrfl(ncexp+1:ncexp+22),'(a3,i8.8,a3,i4.4,a4)')           &
     &            'rdr',ctime/1000_i8,'.pe',mysub,'.bin'

      else

        ncfl=ncexp+31

        write(rdrfl(ncexp+1:ncexp+31),'(a3,i8.8,2(a4,i4.4),a4)')        &
     &            'rdr',ctime/1000_i8,'.grp',mygrp,'-sub',mysub,'.bin'

      end if

      open(iordr,iostat=stat,err=140,                                   &
     &     file=crsdir(1:nccrs)//rdrfl(1:ncfl),                         &
     &     status='new',access='direct',form='unformatted',             &
     &     recl=siz,action='write')

  140 if(stat.eq.0) then

        extpe=0

      else

        extpe=1

        ncfl=ncfl+5

        write(rdrfl(ncfl-4:ncfl),'(a5)') '.swap'

        open(iordr,iostat=stat,err=150,                                 &
     &       file=crsdir(1:nccrs)//rdrfl(1:ncfl),                       &
     &       status='new',access='direct',form='unformatted',           &
     &       recl=siz,action='write')

      end if

  150 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outrdr  ',6,'cont',1,'              ',14,iordr, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outrdr  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      call chkopen(extpe,iordr)

      if(mype.eq.stat-1) then

        if(fpara(1:5).eq.'multi') then

          if(ngrp.eq.1) then

            write(rdrfl(ncexp+15:ncexp+18),'(a4)') 'XXXX'

          else

            write(rdrfl(ncexp+16:ncexp+19),'(a4)') 'XXXX'
            write(rdrfl(ncexp+24:ncexp+27),'(a4)') 'YYYY'

          end if

        end if

        call outstd03('outrdr  ',6,rdrfl,ncfl,iordr,1,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read in the data to the interpolated radar data file.

      recrdr=0

      if(rdrvar(1:1).eq.'o') then
        recrdr=recrdr+1

        write(iordr,rec=recrdr,iostat=stat,err=160)                     &
     &       (((u(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(rdrvar(2:2).eq.'o') then
        recrdr=recrdr+1

        write(iordr,rec=recrdr,iostat=stat,err=160)                     &
     &       (((v(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(rdrvar(3:3).eq.'o') then
        recrdr=recrdr+1

        write(iordr,rec=recrdr,iostat=stat,err=160)                     &
     &       (((w(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

      if(rdrvar(4:4).eq.'o') then
        recrdr=recrdr+1

        write(iordr,rec=recrdr,iostat=stat,err=160)                     &
     &       (((qp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nk)

      end if

  160 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outrdr  ',6,'cont',3,'              ',14,iordr, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outrdr  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('outrdr  ',6,rdrfl,108,iordr,4,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the interpolated radar data file.

      close(iordr,iostat=stat,err=170,status='keep')

  170 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outrdr  ',6,'cont',2,'              ',14,iordr, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outrdr  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('outrdr  ',6,rdrfl,108,iordr,2,1,ctime)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      call putunit(iordr)

! -----

!! -----

!!! -----

      end subroutine s_outrdr

!-----7--------------------------------------------------------------7--

      end module m_outrdr
