!***********************************************************************
      module m_rdcheck
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/04/11
!     Modification: 2007/07/30, 2007/08/24, 2008/01/11, 2008/04/17,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/01/30,
!                   2009/02/27, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the dumped data checking file.

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

      public :: rdcheck, s_rdcheck

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdcheck

        module procedure s_rdcheck

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
      subroutine s_rdcheck(fpexprim,fpcrsdir,fpncexp,fpnccrs,           &
     &                     fproc,ctime)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: fproc
                       ! Control flag of processing type

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

      character(len=108) chkfl
                       ! Opened file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer in       ! Namelist table index

      integer iochk    ! Unit number of checking file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

!-----7--------------------------------------------------------------7--

!! Open and read out the data from the dumped data checking file.

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

! Initialize the character variable.

      if(mype.eq.root) then

        call inichar(chkfl)

      end if

! -----

! Get the unit number.

      if(mype.eq.root) then

        call getunit(iochk)

      end if

! -----

! Open the dumped data checking file.

      if(mype.eq.root) then

        chkfl(1:ncexp)=exprim(1:ncexp)

        if(fproc(1:3).eq.'dmp') then

          write(chkfl(ncexp+1:ncexp+13),'(a13)') 'dmp.check.txt'

          open(iochk,iostat=stat,err=100,                               &
     &         file=crsdir(1:nccrs)//chkfl(1:ncexp+13),                 &
     &         status='old',access='sequential',form='formatted',       &
     &         blank='null',position='rewind',action='read')

        else if(fproc(1:3).eq.'mon') then

          write(chkfl(ncexp+1:ncexp+13),'(a13)') 'mon.check.txt'

          open(iochk,iostat=stat,err=100,                               &
     &         file=crsdir(1:nccrs)//chkfl(1:ncexp+13),                 &
     &         status='old',access='sequential',form='formatted',       &
     &         blank='null',position='rewind',action='read')

        else if(fproc(1:3).eq.'geo') then

          write(chkfl(ncexp+1:ncexp+19),'(a19)') 'geography.check.txt'

          open(iochk,iostat=stat,err=100,                               &
     &         file=crsdir(1:nccrs)//chkfl(1:ncexp+19),                 &
     &         status='old',access='sequential',form='formatted',       &
     &         blank='null',position='rewind',action='read')

        end if

      else

        stat=0

      end if

  100 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdcheck ',7,'cont',1,'              ',14,       &
     &                 iochk,stat)

        end if

        call cpondpe

        call destroy('rdcheck ',7,'stop',1001,'              ',14,      &
     &               101,stat)

      end if

      if(mype.eq.root) then

        if(fproc(1:3).eq.'dmp') then

          call outstd03('rdcheck ',7,chkfl,ncexp+13,iochk,1,1,ctime)

        else if(fproc(1:3).eq.'mon') then

          call outstd03('rdcheck ',7,chkfl,ncexp+13,iochk,1,1,ctime)

        else if(fproc(1:3).eq.'geo') then

          call outstd03('rdcheck ',7,chkfl,ncexp+19,iochk,1,0,0_i8)

        end if

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the dumped data checking file.

      if(mype.eq.root) then

        read(iochk,'(a)',iostat=stat,end=110,err=110)                   &
     &      (rcname(in),in=1,ncn)

        read(iochk,*,iostat=stat,end=110,err=110)                       &
     &      (riname(in),in=1,nin)

        read(iochk,*,iostat=stat,end=110,err=110)                       &
     &      (rrname(in),in=1,nrn)

      else

        stat=0

      end if

  110 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdcheck ',7,'cont',3,'              ',14,       &
     &                 iochk,stat)

        end if

        call cpondpe

        call destroy('rdcheck ',7,'stop',1001,'              ',14,      &
     &               101,stat)

      end if

      if(mype.eq.root) then

        if(fproc(1:3).eq.'dmp') then

          call outstd03('rdcheck ',7,chkfl,108,iochk,3,1,ctime)

        else if(fproc(1:3).eq.'mon') then

          call outstd03('rdcheck ',7,chkfl,108,iochk,3,1,ctime)

        else if(fproc(1:3).eq.'geo') then

          call outstd03('rdcheck ',7,chkfl,108,iochk,3,0,0_i8)

        end if

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Check the dumped data checking file.

      if(mype.eq.root) then

        call chkfile(fproc,stat,ncn,nin,nrn,                            &
     &               cname,iname,rname,rcname,riname,rrname)

      else

        stat=0

      end if

      call chkerr(stat)

      if(stat.lt.0) then

        call destroy('chkfile ',7,'stop',1001,'              ',14,      &
     &               101,stat)

      end if

! -----

! Close the dumped data checking file.

      if(mype.eq.root) then

        close(iochk,iostat=stat,err=120,status='keep')

      else

        stat=0

      end if

  120 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdcheck ',7,'cont',2,'              ',14,       &
     &                 iochk,stat)

        end if

        call cpondpe

        call destroy('rdcheck ',7,'stop',1001,'              ',14,      &
     &               101,stat)

      end if

      if(mype.eq.root) then

        if(fproc(1:3).eq.'dmp') then

          call outstd03('rdcheck ',7,chkfl,108,iochk,2,1,ctime)

        else if(fproc(1:3).eq.'mon') then

          call outstd03('rdcheck ',7,chkfl,108,iochk,2,1,ctime)

        else if(fproc(1:3).eq.'geo') then

          call outstd03('rdcheck ',7,chkfl,108,iochk,2,0,0_i8)

        end if

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      if(mype.eq.root) then

        call putunit(iochk)

      end if

! -----

!! -----

      end subroutine s_rdcheck

!-----7--------------------------------------------------------------7--

      end module m_rdcheck
