!***********************************************************************
      module m_rdgrp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/12/04
!     Modification: 2007/01/05, 2007/01/20, 2007/04/11, 2007/08/24,
!                   2007/10/19, 2008/05/02, 2008/07/25, 2008/08/25,
!                   2008/10/10, 2009/01/05, 2009/02/27, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     open and read out the data from the group domain arrangement file
!     and check them.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_castgrp
      use m_chkerr
      use m_chkstd
      use m_comgrp
      use m_comkind
      use m_commpi
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

      public :: rdgrp, s_rdgrp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rdgrp

        module procedure s_rdgrp

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
      subroutine s_rdgrp(fpexprim,fpcrsdir,fpncexp,fpnccrs)
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

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) grpfl
                       ! Opened file name

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer igc      ! Array index in x direction
      integer jgc      ! Array index in y direction

      integer iered    ! Index of group domain at east boundary
                       ! in reductional entire domain
      integer jnred    ! Index of group domain at north boundary
                       ! in reductional entire domain

      integer nsrtmp   ! Temporary used maximum serial number
                       ! of group domain in entire domain

      integer iogrp    ! Unit number of group domain arrangement file

      integer stat     ! Runtime status

      integer broot    ! Broadcasting root

! Internal private variables

      integer igc_sub  ! Substitute for igc
      integer jgc_sub  ! Substitute for jgc

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

! -----

!! Open and read out the data from the group domain arrangement file.

      if(ngrp.ge.3.and.ngrp.ne.nsrl) then

! Initialize the character variable.

        if(mype.eq.root) then

          call inichar(grpfl)

        end if

! -----

! Get the unit number.

        if(mype.eq.root) then

          call getunit(iogrp)

        end if

! -----

! Open the group domain arrangement file.

        if(mype.eq.root) then

          grpfl(1:ncexp)=exprim(1:ncexp)

          write(grpfl(ncexp+1:ncexp+11),'(a11)') 'arrange.txt'

          open(iogrp,iostat=stat,err=100,                               &
     &         file=crsdir(1:nccrs)//grpfl(1:ncexp+11),                 &
     &         status='old',access='sequential',form='formatted',       &
     &         blank='null',position='rewind',action='read')

        else

          stat=0

        end if

  100   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('rdgrp   ',5,'cont',1,'              ',14,     &
     &                   iogrp,stat)

          end if

          call cpondpe

          call destroy('rdgrp   ',5,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('rdgrp   ',5,grpfl,ncexp+11,iogrp,1,0,0_i8)

        end if

        broot=stat-1

        call chkstd(broot)

! -----

! Read out the data from the group domain arrangement file.

        if(mype.eq.root) then

          do_head: do

            read(iogrp,'(a1)',iostat=stat,end=110,err=110)              &
     &        chrgrp(1,1)(1:1)

            if(chrgrp(1,1)(1:1).ne.'#') then

              backspace(iogrp,iostat=stat,err=110)

              exit do_head

            end if

          end do do_head

          do jgc=njgrp,1,-1

            read(iogrp,'(4999a1)',iostat=stat,end=110,err=110)          &
     &          (chrgrp(igc,jgc)(1:1),igc=1,nigrp)

          end do

        else

          stat=0

        end if

  110   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('rdgrp   ',5,'cont',3,'              ',14,     &
     &                   iogrp,stat)

          end if

          call cpondpe

          call destroy('rdgrp   ',5,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('rdgrp   ',5,grpfl,108,iogrp,3,0,0_i8)

        end if

        broot=stat-1

        call chkstd(broot)

! -----

! Close the group domain arrangement file.

        if(mype.eq.root) then

          close(iogrp,iostat=stat,err=120,status='keep')

        else

          stat=0

        end if

  120   call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('rdgrp   ',5,'cont',2,'              ',14,     &
     &                   iogrp,stat)

          end if

          call cpondpe

          call destroy('rdgrp   ',5,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(mype.eq.root) then

          call outstd03('rdgrp   ',5,grpfl,108,iogrp,2,0,0_i8)

        end if

        broot=stat-1

        call chkstd(broot)

! -----

! Return the unit number.

        if(mype.eq.root) then

          call putunit(iogrp)

        end if

! -----

! Broadcast the group domain arrangemnet table.

        call castgrp

! -----

!! -----

! Initialize the group domain arrangement.

      else

        do jgc=1,njgrp
        do igc=1,nigrp

          write(chrgrp(igc,jgc)(1:1),'(a1)') 'o'

        end do
        end do

      end if

! -----

! Check the group domain arrangement.

      nsrtmp=-1

      do_jgc: do jgc=1,njgrp
      do_igc: do igc=1,nigrp

        if(chrgrp(igc,jgc)(1:1).eq.'x') then

          grpxy(igc,jgc)=-1

        else if(chrgrp(igc,jgc)(1:1).eq.'o') then

          nsrtmp=nsrtmp+1

          if(nsrtmp.ge.nsrl) then

            exit do_jgc

          end if

          grpxy(igc,jgc)=nsrtmp

          xgrp(nsrtmp)=igc
          ygrp(nsrtmp)=jgc

        else

          nsrtmp=-1

          exit do_jgc

        end if

      end do do_igc
      end do do_jgc

      if(nsrtmp.eq.nsrl-1) then
        stat=0
      else
        stat=1
      end if

! -----

! If error occured, call the procedure destroy.

      call chkerr(stat)

      if(stat.lt.0) then

        if(mype.eq.-stat-1) then

          call destroy('rdgrp   ',5,'cont',11,'              ',14,101,  &
     &                 stat)

        end if

        call cpondpe

        call destroy('rdgrp   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

!! Set the parameters of reductional entire domain.

! Initialize the parameters of reductional entire domain.

      iwred=nigrp+1
      jsred=njgrp+1

      iered=0
      jnred=0

! -----

! Get the parameters of reductional entire domain.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(igc_sub,jgc_sub)                     &
!$omp&   reduction(min: iwred,jsred) reduction(max: iered,jnred)

      do jgc_sub=1,njgrp
      do igc_sub=1,nigrp

        if(grpxy(igc_sub,jgc_sub).ge.0) then

          if(igc_sub.lt.iwred) then
            iwred=igc_sub
          end if

          if(igc_sub.gt.iered) then
            iered=igc_sub
          end if

          if(jgc_sub.lt.jsred) then
            jsred=jgc_sub
          end if

          if(jgc_sub.gt.jnred) then
            jnred=jgc_sub
          end if

        end if

      end do
      end do

!$omp end do

!$omp end parallel

! -----

! Finally reset the parameters of reductional entire domain.

      iwred=iwred-1
      jsred=jsred-1

      iered=iered-1
      jnred=jnred-1

      nired=iered-iwred+1
      njred=jnred-jsred+1

      nred=nired*njred

! -----

!! -----

      end subroutine s_rdgrp

!-----7--------------------------------------------------------------7--

      end module m_rdgrp
