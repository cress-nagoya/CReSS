!***********************************************************************
      module m_rdtund
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/08/27
!     Modification: 2003/03/28, 2003/04/30, 2003/05/19, 2003/06/27,
!                   2003/12/12, 2004/01/09, 2004/05/31, 2004/08/01,
!                   2004/08/20, 2005/01/14, 2005/02/10, 2005/04/04,
!                   2006/04/03, 2006/09/21, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/08/24, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2008/10/10, 2009/01/30, 2009/02/27,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the soil and sea temperature data from the
!     restart file.

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

      public :: rdtund, s_rdtund

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdtund

        module procedure s_rdtund

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
      subroutine s_rdtund(fpcrsdir,fpprvres,fpnccrs,fpncprv,fpadvopt,   &
     &                    ni,nj,nund,tund,tundp)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcrsdir
                       ! Formal parameter of unique index of crsdir

      integer, intent(in) :: fpprvres
                       ! Formal parameter of unique index of prvres

      integer, intent(in) :: fpnccrs
                       ! Formal parameter of unique index of nccrs

      integer, intent(in) :: fpncprv
                       ! Formal parameter of unique index of ncprv

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

! Output variables

      real, intent(out) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at present

      real, intent(out) :: tundp(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at past

! Internal shared variables

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) prvres
                       ! Restart file name
                       ! without extension of previous running

      character(len=108) resfl
                       ! Opened file name

      integer nccrs    ! Number of character of crsdir
      integer ncprv    ! Number of character of prvres

      integer advopt   ! Option for advection scheme

      integer in       ! Namelist table index

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer ncfl     ! Number of character of restart file

      integer iores    ! Unit number of restart file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(crsdir)
      call inichar(prvres)

! -----

! Get the required namelist variables.

      call getcname(fpcrsdir,crsdir)
      call getcname(fpprvres,prvres)
      call getiname(fpnccrs,nccrs)
      call getiname(fpncprv,ncprv)
      call getiname(fpadvopt,advopt)

! -----

!! Open and read out the data from the restart checking file.

! Initialize the character variable.

      if(mype.eq.root) then

        call inichar(resfl)

      end if

! -----

! Get the unit number.

      if(mype.eq.root) then

        call getunit(iores)

      end if

! -----

! Open the restart checking file.

      if(mype.eq.root) then

        resfl(1:ncprv)=prvres(1:ncprv)

        write(resfl(ncprv+1:ncprv+10),'(a10)') '.check.txt'

        open(iores,iostat=stat,err=100,                                 &
     &       file=crsdir(1:nccrs)//resfl(1:ncprv+10),                   &
     &       status='old',access='sequential',form='formatted',         &
     &       blank='null',position='rewind',action='read')

      else

        stat=0

      end if

  100 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdtund  ',6,'cont',1,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdtund  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdtund  ',6,resfl,ncprv+10,iores,1,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the restart checking file.

      if(mype.eq.root) then

        read(iores,'(a)',iostat=stat,end=110,err=110)                   &
     &      (rcname(in),in=1,ncn)

        read(iores,*,iostat=stat,end=110,err=110)                       &
     &      (riname(in),in=1,nin)

        read(iores,*,iostat=stat,end=110,err=110)                       &
     &      (rrname(in),in=1,nrn)

      else

        stat=0

      end if

  110 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdtund  ',6,'cont',3,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdtund  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdtund  ',6,resfl,108,iores,3,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Check the restart checking file.

      if(mype.eq.root) then

        call chkfile('und',stat,ncn,nin,nrn,                            &
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

! Close the restart checking file.

      if(mype.eq.root) then

        close(iores,iostat=stat,err=120,status='keep')

      else

        stat=0

      end if

  120 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdtund  ',6,'cont',2,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdtund  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.root) then

        call outstd03('rdtund  ',6,resfl,108,iores,2,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      if(mype.eq.root) then

        call putunit(iores)

      end if

! -----

!! -----

!! Open and read out the data from the restart file.

! Initialize the character variable.

      call inichar(resfl)

! -----

! Get the unit number.

      call getunit(iores)

! -----

! Open the restart file.

      resfl(1:ncprv)=prvres(1:ncprv)

      if(ngrp.eq.1) then

        ncfl=ncprv+11

        write(resfl(ncprv+1:ncprv+11),'(a3,i4.4,a4)') '.pe',mysub,'.bin'

      else

        ncfl=ncprv+20

        write(resfl(ncprv+1:ncprv+20),'(2(a4,i4.4),a4)')                &
     &                                '.grp',mygrp,'-sub',mysub,'.bin'

      end if

      open(iores,iostat=stat,err=130,                                   &
     &     file=crsdir(1:nccrs)//resfl(1:ncfl),                         &
     &     status='old',access='sequential',form='unformatted',         &
     &     position='rewind',action='read')

  130 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdtund  ',6,'cont',1,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdtund  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        if(fpara(1:5).eq.'multi') then

          if(ngrp.eq.1) then

            write(resfl(ncprv+4:ncprv+7),'(a4)') 'XXXX'

          else

            write(resfl(ncprv+5:ncprv+8),'(a4)') 'XXXX'
            write(resfl(ncprv+13:ncprv+16),'(a4)') 'YYYY'

          end if

        end if

        call outstd03('rdtund  ',6,resfl,ncfl,iores,1,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Read out the data from the restart file.

      if(advopt.le.3) then

        read(iores,iostat=stat,end=140,err=140)                         &
     &      (((tund(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nund)

      end if

      read(iores,iostat=stat,end=140,err=140)                           &
     &    (((tundp(i,j,k),i=0,ni+1),j=0,nj+1),k=1,nund)

  140 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdtund  ',6,'cont',3,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdtund  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('rdtund  ',6,resfl,108,iores,3,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Close the restart file.

      close(iores,iostat=stat,err=150,status='keep')

  150 call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdtund  ',6,'cont',2,'              ',14,iores, &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdtund  ',6,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

      if(mype.eq.stat-1) then

        call outstd03('rdtund  ',6,resfl,108,iores,2,0,0_i8)

      end if

      broot=stat-1

      call chkstd(broot)

! -----

! Return the unit number.

      call putunit(iores)

! -----

!! -----

      end subroutine s_rdtund

!-----7--------------------------------------------------------------7--

      end module m_rdtund
