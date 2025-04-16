!***********************************************************************
      module m_outtrn
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/04/06
!     Modification: 1999/06/14, 1999/06/21, 1999/09/30, 1999/11/01,
!                   1999/11/19, 2000/01/05, 2000/01/17, 2000/04/18,
!                   2001/02/13, 2001/04/15, 2001/05/29, 2002/04/02,
!                   2002/06/18, 2002/07/15, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2003/06/27, 2004/01/09, 2004/05/31,
!                   2004/08/20, 2004/08/31, 2004/09/25, 2005/01/14,
!                   2005/02/10, 2006/09/21, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/08/24, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/01/30, 2009/02/27, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read in the data to the interpolated terrain file.

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

      public :: outtrn, s_outtrn

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outtrn

        module procedure s_outtrn

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
      subroutine s_outtrn(fpexprim,fpcrsdir,fpncexp,fpnccrs,fpwlngth,   &
     &                    dvname,ncdvn,fmsg,ni,nj,ht)
!***********************************************************************

! Input variables

      character(len=12), intent(in) :: dvname
                       ! Optional data variable name

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

      integer, intent(in) :: ncdvn
                       ! Number of character of dvname

      integer, intent(in) :: fmsg
                       ! Control flag of message type in outstd03

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      real, intent(in) :: ht(0:ni+1,0:nj+1)
                       ! Terrain height

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) trnfl
                       ! Opened file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer wlngth   ! Word length of direct access file

      integer in       ! Namelist table index

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer ncfl     ! Number of character of terrain file name

      integer extpe    ! Control flag of file name extension

      integer iotrn    ! Unit number of terrain file

      integer siz      ! Record length of terrain file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

!-----7--------------------------------------------------------------7--

!!! Open and read in the data to the interpolated terrain file.

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

!! Open and read in the data to the interpolated terrain checking file.

! Initialize the character variable.

      if(mype.eq.root) then

        call inichar(trnfl)

      end if

! -----

! Get the unit number.

      if(mype.eq.root) then

        call getunit(iotrn)

      end if

! -----

! Open the interpolated terrain checking file.

      if(mype.eq.root) then

        trnfl(1:ncexp)=exprim(1:ncexp)
        trnfl(ncexp+1:ncexp+ncdvn)=dvname(1:ncdvn)

        write(trnfl(ncexp+ncdvn+1:ncexp+ncdvn+10),'(a10)') '.check.txt'

        open(iotrn,iostat=stat,err=100,                                 &
     &       file=crsdir(1:nccrs)//trnfl(1:ncexp+ncdvn+10),             &
     &       status='new',access='sequential',form='formatted',         &
     &       blank='null',action='write')

  100   if(stat.eq.0) then

          ncfl=ncexp+ncdvn+10

        else

          ncfl=ncexp+ncdvn+15

          write(trnfl(ncexp+ncdvn+11:ncexp+ncdvn+15),'(a5)') '.swap'

          open(iotrn,iostat=stat,err=110,                               &
     &         file=crsdir(1:nccrs)//trnfl(1:ncexp+ncdvn+15),           &
     &         status='new',access='sequential',form='formatted',       &
     &         blank='null',action='write')

        end if

      else

        stat=0

      end if

  110 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outtrn  ',6,'cont',1,'              ',14,iotrn, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outtrn  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('outtrn  ',6,trnfl,ncfl,iotrn,1,fmsg,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read in the data to the interpolated terrain checking file.

      if(mype.eq.root) then

        write(iotrn,'(a)',iostat=stat,err=120)                          &
     &       (cname(in),in=1,ncn)

        write(iotrn,*,iostat=stat,err=120)                              &
     &       (iname(in),in=1,nin)

        write(iotrn,*,iostat=stat,err=120)                              &
     &       (rname(in),in=1,nrn)

      else

        stat=0

      end if

  120 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outtrn  ',6,'cont',4,'              ',14,iotrn, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outtrn  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('outtrn  ',6,trnfl,108,iotrn,4,fmsg,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the interpolated terrain checking file.

      if(mype.eq.root) then

        close(iotrn,iostat=stat,err=130,status='keep')

      else

        stat=0

      end if

  130 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outtrn  ',6,'cont',2,'              ',14,iotrn, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outtrn  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('outtrn  ',6,trnfl,108,iotrn,2,fmsg,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      if(mype.eq.root) then

        call putunit(iotrn)

      end if

! -----

!! -----

!! Open and read in the data to the interpolated terrain file.

! Initialize the character variable.

      call inichar(trnfl)

! -----

! Get the unit number.

      call getunit(iotrn)

! -----

! Open the interpolated terrain file.

      siz=(ni+2)*(nj+2)*wlngth

      trnfl(1:ncexp)=exprim(1:ncexp)
      trnfl(ncexp+1:ncexp+ncdvn)=dvname(1:ncdvn)

      if(ngrp.eq.1) then

        ncfl=ncexp+ncdvn+11

        write(trnfl(ncexp+ncdvn+1:ncexp+ncdvn+11),'(a3,i4.4,a4)')       &
     &                                '.pe',mysub,'.bin'

      else

        ncfl=ncexp+ncdvn+20

        write(trnfl(ncexp+ncdvn+1:ncexp+ncdvn+20),'(2(a4,i4.4),a4)')    &
     &                               '.grp',mygrp,'-sub',mysub,'.bin'

      end if

      open(iotrn,iostat=stat,err=140,                                   &
     &     file=crsdir(1:nccrs)//trnfl(1:ncfl),                         &
     &     status='new',access='direct',form='unformatted',             &
     &     recl=siz,action='write')

  140 if(stat.eq.0) then

        extpe=0

      else

        extpe=1

        ncfl=ncfl+5

        write(trnfl(ncfl-4:ncfl),'(a5)') '.swap'

        open(iotrn,iostat=stat,err=150,                                 &
     &       file=crsdir(1:nccrs)//trnfl(1:ncfl),                       &
     &       status='new',access='direct',form='unformatted',           &
     &       recl=siz,action='write')

      end if

  150 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outtrn  ',6,'cont',1,'              ',14,iotrn, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outtrn  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      call chkopen(extpe,iotrn)

      if(mype.eq.stat-1) then

        if(fpara(1:5).eq.'multi') then

          if(ngrp.eq.1) then

            write(trnfl(ncexp+ncdvn+4:ncexp+ncdvn+7),'(a4)') 'XXXX'

          else

            write(trnfl(ncexp+ncdvn+5:ncexp+ncdvn+8),'(a4)') 'XXXX'
            write(trnfl(ncexp+ncdvn+13:ncexp+ncdvn+16),'(a4)') 'YYYY'

          end if

        end if

        call outstd03('outtrn  ',6,trnfl,ncfl,iotrn,1,fmsg,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read in the data to the interpolated terrain file.

      write(iotrn,rec=1,iostat=stat,err=160)                            &
     &     ((ht(i,j),i=0,ni+1),j=0,nj+1)

  160 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outtrn  ',6,'cont',4,'              ',14,iotrn, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outtrn  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('outtrn  ',6,trnfl,108,iotrn,4,fmsg,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the interpolated terrain file.

      close(iotrn,iostat=stat,err=170,status='keep')

  170 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('outtrn  ',6,'cont',2,'              ',14,iotrn, &
     &                 stat)

        end if

        call cpondpe

        call destroy('outtrn  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('outtrn  ',6,trnfl,108,iotrn,2,fmsg,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      call putunit(iotrn)

! -----

!! -----

!!! -----

      end subroutine s_outtrn

!-----7--------------------------------------------------------------7--

      end module m_outtrn
