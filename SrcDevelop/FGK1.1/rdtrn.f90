!***********************************************************************
      module m_rdtrn
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/03/16
!     Modification: 1999/03/25, 1999/04/06, 1999/05/20, 1999/06/14,
!                   1999/06/21, 1999/08/23, 1999/09/30, 1999/11/01,
!                   1999/11/19, 2000/01/05, 2000/01/17, 2000/04/18,
!                   2000/12/06, 2001/02/13, 2001/04/15, 2001/05/29,
!                   2002/04/02, 2002/06/18, 2002/07/15, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2004/01/09,
!                   2004/05/31, 2004/08/20, 2004/09/25, 2005/01/14,
!                   2005/02/10, 2006/09/21, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/08/24, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/01/30, 2009/02/27, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the interpolated terrain file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_chkfile
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

      public :: rdtrn, s_rdtrn

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdtrn

        module procedure s_rdtrn

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
      subroutine s_rdtrn(fpexprim,fpcrsdir,fpncexp,fpnccrs,fpwlngth,    &
     &                   dvname,ncdvn,fmsg,ni,nj,ht)
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

! Output variable

      real, intent(out) :: ht(0:ni+1,0:nj+1)
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

      integer ncfl     ! Number of character of terrain file

      integer iotrn    ! Unit number of terrain file

      integer siz      ! Record length of terrain file

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
      call getiname(fpwlngth,wlngth)

! -----

!! Open and read out the data from the interpolated terrain checking
!! file.

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
     &       status='old',access='sequential',form='formatted',         &
     &       blank='null',position='rewind',action='read')

      else

        stat=0

      end if

  100 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdtrn   ',5,'cont',1,'              ',14,iotrn, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdtrn   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdtrn   ',5,trnfl,ncexp+ncdvn+10,iotrn,1,fmsg,   &
     &                0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the interpolated terrain checking file.

      if(mype.eq.root) then

        read(iotrn,'(a)',iostat=stat,end=110,err=110)                   &
     &      (rcname(in),in=1,ncn)

        read(iotrn,*,iostat=stat,end=110,err=110)                       &
     &      (riname(in),in=1,nin)

        read(iotrn,*,iostat=stat,end=110,err=110)                       &
     &      (rrname(in),in=1,nrn)

      else

        stat=0

      end if

  110 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdtrn   ',5,'cont',3,'              ',14,iotrn, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdtrn   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdtrn   ',5,trnfl,108,iotrn,3,fmsg,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Check the interpolated terrain checking file.

      if(mype.eq.root) then

        call chkfile('trn',stat,ncn,nin,nrn,                            &
     &               cname,iname,rname,rcname,riname,rrname)

      else

        stat=0

      end if

      call chkerr(stat)

      if(stat.lt.0) then

        call destroy('chkfile ',7,'stop',1001,'              ',14,10,   &
     &               stat)

      end if

! -----

! Close the interpolated terrain checking file.

      if(mype.eq.root) then

        close(iotrn,iostat=stat,err=120,status='keep')

      else

        stat=0

      end if

  120 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdtrn   ',5,'cont',2,'              ',14,iotrn, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdtrn   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdtrn   ',5,trnfl,108,iotrn,2,fmsg,0_i8)

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

!! Open and read out the data from the interpolated terrain file.

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

      open(iotrn,iostat=stat,err=130,                                   &
     &     file=crsdir(1:nccrs)//trnfl(1:ncfl),                         &
     &     status='old',access='direct',form='unformatted',             &
     &     recl=siz,action='read')

  130 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdtrn   ',5,'cont',1,'              ',14,iotrn, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdtrn   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        if(fpara(1:5).eq.'multi') then

          if(ngrp.eq.1) then

            write(trnfl(ncexp+ncdvn+4:ncexp+ncdvn+7),'(a4)') 'XXXX'

          else

            write(trnfl(ncexp+ncdvn+5:ncexp+ncdvn+8),'(a4)') 'XXXX'
            write(trnfl(ncexp+ncdvn+13:ncexp+ncdvn+16),'(a4)') 'YYYY'

          end if

        end if

        call outstd03('rdtrn   ',5,trnfl,ncfl,iotrn,1,fmsg,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the interpolated terrain file.

      read(iotrn,rec=1,iostat=stat,err=140)                             &
     &    ((ht(i,j),i=0,ni+1),j=0,nj+1)

  140 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdtrn   ',5,'cont',3,'              ',14,iotrn, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdtrn   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('rdtrn   ',5,trnfl,108,iotrn,3,fmsg,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the interpolated terrain file.

      close(iotrn,iostat=stat,err=150,status='keep')

  150 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdtrn   ',5,'cont',2,'              ',14,iotrn, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdtrn   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('rdtrn   ',5,trnfl,108,iotrn,2,fmsg,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      call putunit(iotrn)

! -----

!! -----

      end subroutine s_rdtrn

!-----7--------------------------------------------------------------7--

      end module m_rdtrn
