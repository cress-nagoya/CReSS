!***********************************************************************
      module m_unimongr
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/08/10
!     Modification: 2001/01/15, 2001/02/13, 2001/03/13, 2001/05/29,
!                   2001/06/29, 2001/08/07, 2001/11/20, 2002/01/15,
!                   2002/04/02, 2002/06/18, 2002/07/03, 2002/09/02,
!                   2002/09/09, 2003/02/05, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2003/06/27, 2004/01/09, 2004/05/31,
!                   2004/06/10, 2004/07/01, 2004/08/01, 2004/08/20,
!                   2004/09/01, 2005/01/14, 2005/02/10, 2006/09/21,
!                   2006/09/30, 2006/12/04, 2007/01/05, 2007/01/20,
!                   2007/04/11, 2007/05/14, 2007/05/21, 2007/07/30,
!                   2007/08/24, 2008/04/17, 2008/05/02, 2008/08/25,
!                   2008/10/10, 2009/02/27, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     generate the united history file from the dumped files in group
!     domain for monitor variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
      use m_commpi
      use m_currpe
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

      public :: unimongr, s_unimongr

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface unimongr

        module procedure s_unimongr

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
      subroutine s_unimongr(fpexprim,fpcrsdir,fpncexp,fpnccrs,          &
     &                      fpwlngth,fprmopt_uni,ctime,nx,ny,           &
     &                      ni,nj,ni_uni,nj_uni,var,nio_uni,iodmp)
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

      integer, intent(in) :: fprmopt_uni
                       ! Formal parameter of unique index of rmopt_uni

      integer(kind=i8), intent(in) :: ctime
                       ! Model current forecast time

      integer, intent(in) :: nx
                       ! Composite model dimension in x direction

      integer, intent(in) :: ny
                       ! Composite model dimension in y direction

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: ni_uni
                       ! Model dimension of unite in x direction

      integer, intent(in) :: nj_uni
                       ! Model dimension of unite in y direction

      integer, intent(in) :: nio_uni
                       ! Maximum unit number of unite

! Internal shared variables

      character(len=108) exprim
                       ! Optional run name

      character(len=108) crsdir
                       ! User specified directory for CReSS files

      character(len=108) unifl
                       ! Name of united file

      character(len=108) dmpfl
                       ! Name of dumped file

      integer ncexp    ! Number of character of exprim
      integer nccrs    ! Number of character of crsdir

      integer wlngth   ! Word length of direct access file

      integer rmopt_uni
                       ! Option for original dumped files removing

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer istr     ! Start index in x direction
      integer iend     ! End index in x direction

      integer jstr     ! Start index in y direction
      integer jend     ! End index in y direction

      integer njrec    ! Number of record for each i and j plane

      integer recmst   ! Tail of record number at each process
                       ! in outermost loops

      integer recj     ! Tail of record number at each process
                       ! in y direction

      integer ncuni    ! Number of character of united file
      integer ncdmp    ! Number of character of dumped file

      integer recuni   ! Current record number of united file
      integer recdmp   ! Current record number of dumped file

      integer sizuni   ! Record length of united file
      integer sizdmp   ! Record length of dumped file

      integer stat     ! Runtime status

      integer iouni    ! Unit number of united file

      integer, intent(inout) :: iodmp(0:nio_uni)
                       ! Unit number of dumped file

      real, intent(inout) :: var(2:ni_uni-2,2:nj_uni-2)
                       ! Optional variable in dumped file

!-----7--------------------------------------------------------------7--

! Initialize character variables.

      call inichar(exprim)
      call inichar(crsdir)

! -----

! Get the required namelist variables.

      call getcname(fpexprim,exprim)
      call getcname(fpcrsdir,crsdir)
      call getiname(fpncexp,ncexp)
      call getiname(fpnccrs,nccrs)
      call getiname(fpwlngth,wlngth)
      call getiname(fprmopt_uni,rmopt_uni)

! -----

! Set the common used variable.

      njrec=nsub*(nj-3)

! -----

!!!! Generate the united file from the dumped files for monitor
!!!! variables.

! Initialize the character variable.

      call inichar(unifl)

! -----

! Get the unit number.

      call getunit(iouni)

! -----

! Open the united file.

      if(nio_uni.le.nsub-1) then

        if(nio_uni.le.nisub-1) then

          sizuni=(ni-3)*wlngth

        else

          sizuni=(nx-3)*(nj-3)*wlngth

        end if

      else

        sizuni=(nx-3)*(ny-3)*wlngth

      end if

      unifl(1:ncexp)=exprim(1:ncexp)

      if(ngrp.eq.1) then

        ncuni=ncexp+22

        write(unifl(ncexp+1:ncexp+22),'(a3,i8.8,a11)')                  &
     &            'mon',ctime/1000_i8,'.united.bin'

      else

        ncuni=ncexp+30

        write(unifl(ncexp+1:ncexp+30),'(a3,i8.8,a4,i4.4,a11)')          &
     &            'mon',ctime/1000_i8,'.grp',mygrp,'.united.bin'

      end if

      open(iouni,iostat=stat,err=100,                                   &
     &     file=crsdir(1:nccrs)//unifl(1:ncuni),                        &
     &     status='new',access='direct',form='unformatted',             &
     &     recl=sizuni,action='write')

  100 if(stat.ne.0) then

        call destroy('unimongr',8,'stop',1,'              ',14,iouni,   &
     &               stat)

      end if

      call outstd03('unimongr',8,unifl,ncuni,iouni,1,1,ctime)

! -----

!!! Perform uniting.

!! For single processing case.

      if(nio_uni.le.nsub-1) then

        if(nio_uni.le.nisub-1) then

          do mysub=0,nsub-1

! Set the parameters of parallelizing.

            call currpe('unite   ',5,'ijsub')

! -----

! Get the unit number.

            call getunit(iodmp(0))

! -----

! Open the dumped file.

            sizdmp=(ni-3)*(nj-3)*wlngth

            dmpfl(1:ncexp)=exprim(1:ncexp)

            if(ngrp.eq.1) then

              ncdmp=ncexp+22

              write(dmpfl(ncexp+1:ncexp+22),'(a3,i8.8,a3,i4.4,a4)')     &
     &                  'mon',ctime/1000_i8,'.pe',mysub,'.bin'

            else

              ncdmp=ncexp+31

              write(dmpfl(ncexp+1:ncexp+31),'(a3,i8.8,2(a4,i4.4),a4)')  &
     &             'mon',ctime/1000_i8,'.grp',mygrp,'-sub',mysub,'.bin'

            end if

            if(rmopt_uni.eq.1) then

              open(iodmp(0),iostat=stat,err=110,                        &
     &             file=crsdir(1:nccrs)//dmpfl(1:ncdmp),                &
     &             status='old',access='direct',form='unformatted',     &
     &             recl=sizdmp,action='readwrite')

            else

              open(iodmp(0),iostat=stat,err=110,                        &
     &             file=crsdir(1:nccrs)//dmpfl(1:ncdmp),                &
     &             status='old',access='direct',form='unformatted',     &
     &             recl=sizdmp,action='read')

            end if

  110       if(stat.ne.0) then

              call destroy('unimongr',8,'cont',1,'              ',14,   &
     &                     iodmp(0),stat)

              go to 220

            end if

            call outstd03('unimongr',8,dmpfl,ncdmp,iodmp(0),1,1,ctime)

! -----

! Read out the data from the dumped file and read in the data to the
! united file.

            recdmp=0

            recj=nisub*jsub*(nj-3)+isub+1

            do_inf_1: do

              recdmp=recdmp+1

              recmst=njrec*(recdmp-1)+recj

              read(iodmp(0),rec=recdmp,iostat=stat,err=120)             &
     &            ((var(i,j),i=2,ni-2),j=2,nj-2)

              do j=2,nj-2
                recuni=recmst+(j-2)*nisub

                write(iouni,rec=recuni,iostat=stat,err=210)             &
     &               (var(i,j),i=2,ni-2)

              end do

            end do do_inf_1

  120       if(recdmp.eq.1) then

              call destroy('unimongr',8,'cont',3,'              ',14,   &
     &                     iodmp(0),stat)

              go to 220

            else

              call outstd03('unimongr',8,dmpfl,108,iodmp(0),3,1,ctime)

            end if

            call outstd03('unimongr',8,unifl,108,iouni,4,1,ctime)

! -----

! Close the dumped file.

            if(rmopt_uni.eq.1) then

              close(iodmp(0),iostat=stat,err=130,status='delete')

            else

              close(iodmp(0),iostat=stat,err=130,status='keep')

            end if

  130       if(stat.ne.0) then

              call destroy('unimongr',8,'cont',2,'              ',14,   &
     &                     iodmp(0),stat)

              go to 220

            end if

            call outstd03('unimongr',8,unifl,108,iodmp(0),2,1,ctime)

! -----

! Return the unit number.

            call putunit(iodmp(0))

! -----

          end do

!! -----

!! For batch processing case in x direction.

        else

          do jsub=0,njsub-1

! Get the unit number.

            do isub=0,nisub-1

              call getunit(iodmp(isub))

            end do

! -----

! Open the dumped file.

            sizdmp=(ni-3)*(nj-3)*wlngth

            dmpfl(1:ncexp)=exprim(1:ncexp)

            if(ngrp.eq.1) then

              ncdmp=ncexp+22

              write(dmpfl(ncexp+1:ncexp+14),'(a3,i8.8,a3)')             &
     &                  'mon',ctime/1000_i8,'.pe'

            else

              ncdmp=ncexp+31

              write(dmpfl(ncexp+1:ncexp+23),'(a3,i8.8,a4,i4.4,a4)')     &
     &                  'mon',ctime/1000_i8,'.grp',mygrp,'-sub'

            end if

            do isub=0,nisub-1

              call currpe('unite   ',5,'mysub')

              write(dmpfl(ncdmp-7:ncdmp),'(i4.4,a4)') mysub,'.bin'

              if(rmopt_uni.eq.1) then

               open(iodmp(isub),iostat=stat,err=140,                    &
     &              file=crsdir(1:nccrs)//dmpfl(1:ncdmp),               &
     &              status='old',access='direct',form='unformatted',    &
     &              recl=sizdmp,action='readwrite')

              else

               open(iodmp(isub),iostat=stat,err=140,                    &
     &              file=crsdir(1:nccrs)//dmpfl(1:ncdmp),               &
     &              status='old',access='direct',form='unformatted',    &
     &              recl=sizdmp,action='read')

              end if

  140         if(stat.ne.0) then

                call destroy('unimongr',8,'cont',1,'              ',14, &
     &                       iodmp(isub),stat)

                go to 220

              end if

              call outstd03('unimongr',8,dmpfl,ncdmp,iodmp(isub),1,1,   &
     &                      ctime)

            end do

! -----

! Read out the data from the dumped file and read in the data to the
! united file.

            recdmp=0

            do_inf_2: do

              recdmp=recdmp+1

              do isub=0,nisub-1

                istr=isub*(ni-3)+2
                iend=(isub+1)*(ni-3)+1

                read(iodmp(isub),rec=recdmp,iostat=stat,err=150)        &
     &              ((var(i,j),i=istr,iend),j=2,nj-2)

              end do

              recuni=njsub*(recdmp-1)+jsub+1

              write(iouni,rec=recuni,iostat=stat,err=210)               &
     &             ((var(i,j),i=2,nx-2),j=2,nj-2)

            end do do_inf_2

  150       if(recdmp.eq.1) then

              call destroy('unimongr',8,'cont',3,'              ',14,   &
     &                     iodmp(isub),stat)

              go to 220

            else

              do isub=0,nisub-1

                call outstd03('unimongr',8,dmpfl,108,iodmp(isub),3,1,   &
     &                        ctime)

              end do

            end if

            call outstd03('unimongr',8,unifl,108,iouni,4,1,ctime)

! -----

! Close the dumped file.

            do isub=nisub-1,0,-1

              if(rmopt_uni.eq.1) then

               close(iodmp(isub),iostat=stat,err=160,status='delete')

              else

               close(iodmp(isub),iostat=stat,err=160,status='keep')

              end if

  160         if(stat.ne.0) then

                call destroy('unimongr',8,'cont',2,'              ',14, &
     &                       iodmp(isub),stat)

                go to 220

              end if

              call outstd03('unimongr',8,dmpfl,108,iodmp(isub),2,1,     &
     &                      ctime)

            end do

! -----

! Return the unit number.

            do isub=nisub-1,0,-1

              call putunit(iodmp(isub))

            end do

! -----

          end do

        end if

!! -----

!! For batch processing case in x and y directions.

      else

! Get the unit number.

        do mysub=0,nsub-1

          call getunit(iodmp(mysub))

        end do

! -----

! Open the dumped file.

        sizdmp=(ni-3)*(nj-3)*wlngth

        dmpfl(1:ncexp)=exprim(1:ncexp)

        if(ngrp.eq.1) then

          ncdmp=ncexp+22

          write(dmpfl(ncexp+1:ncexp+14),'(a3,i8.8,a3)')                 &
     &              'mon',ctime/1000_i8,'.pe'

        else

          ncdmp=ncexp+31

          write(dmpfl(ncexp+1:ncexp+23),'(a3,i8.8,a4,i4.4,a4)')         &
     &              'mon',ctime/1000_i8,'.grp',mygrp,'-sub'

        end if

        do mysub=0,nsub-1

          write(dmpfl(ncdmp-7:ncdmp),'(i4.4,a4)') mysub,'.bin'

          if(rmopt_uni.eq.1) then

            open(iodmp(mysub),iostat=stat,err=170,                      &
     &           file=crsdir(1:nccrs)//dmpfl(1:ncdmp),                  &
     &           status='old',access='direct',form='unformatted',       &
     &           recl=sizdmp,action='readwrite')

          else

            open(iodmp(mysub),iostat=stat,err=170,                      &
     &           file=crsdir(1:nccrs)//dmpfl(1:ncdmp),                  &
     &           status='old',access='direct',form='unformatted',       &
     &           recl=sizdmp,action='read')

          end if

  170     if(stat.ne.0) then

            call destroy('unimongr',8,'cont',1,'              ',14,     &
     &                   iodmp(mysub),stat)

            go to 220

          end if

          call outstd03('unimongr',8,dmpfl,ncdmp,iodmp(mysub),1,1,ctime)

        end do

! -----

! Read out the data from the dumped file and read in the data to the
! united file.

        recdmp=0

        do_inf_3: do

          recdmp=recdmp+1

          do mysub=0,nsub-1

            call currpe('unite   ',5,'ijsub')

            istr=isub*(ni-3)+2
            iend=(isub+1)*(ni-3)+1

            jstr=jsub*(nj-3)+2
            jend=(jsub+1)*(nj-3)+1

            read(iodmp(mysub),rec=recdmp,iostat=stat,err=180)           &
     &          ((var(i,j),i=istr,iend),j=jstr,jend)

          end do

          recuni=recdmp

          write(iouni,rec=recuni,iostat=stat,err=210)                   &
     &         ((var(i,j),i=2,nx-2),j=2,ny-2)

        end do do_inf_3

  180   if(recdmp.eq.1) then

          call destroy('unimongr',8,'cont',3,'              ',14,       &
     &                 iodmp(mysub),stat)

          go to 220

        else

          do mysub=0,nsub-1

            call outstd03('unimongr',8,dmpfl,108,iodmp(mysub),3,1,ctime)

          end do

        end if

        call outstd03('unimongr',8,unifl,108,iouni,4,1,ctime)

! -----

! Close the dumped file.

        do mysub=nsub-1,0,-1

          if(rmopt_uni.eq.1) then

            close(iodmp(mysub),iostat=stat,err=190,status='delete')

          else

            close(iodmp(mysub),iostat=stat,err=190,status='keep')

          end if

  190     if(stat.ne.0) then

            call destroy('unimongr',8,'cont',2,'              ',14,     &
     &                   iodmp(mysub),stat)

            go to 220

          end if

          call outstd03('unimongr',8,dmpfl,108,iodmp(mysub),2,1,ctime)

        end do

! -----

! Return the unit number.

        do mysub=nsub-1,0,-1

          call putunit(iodmp(mysub))

        end do

! -----

      end if

!! -----

!!! -----

! Close the united file.

      close(iouni,iostat=stat,err=200,status='keep')

  200 if(stat.ne.0) then

        call destroy('unimongr',8,'stop',2,'              ',14,iouni,   &
     &               stat)

      end if

      call outstd03('unimongr',8,unifl,108,iouni,2,1,ctime)

! -----

! Return the unit number.

      call putunit(iouni)

! -----

!!!! -----

      return

! If error occured, call the procedure destroy.

  210 call destroy('unimongr',8,'cont',4,'              ',14,iouni,     &
     &             stat)

  220 close(iouni,iostat=stat,err=230,status='delete')

  230 if(stat.ne.0) then

        call destroy('unimongr',8,'stop',2,'              ',14,iouni,   &
     &               stat)

      end if

      call outstd03('unimongr',8,unifl,108,iouni,2,1,ctime)

      call destroy('unimongr',8,'stop',1001,'              ',14,101,    &
     &             stat)

! -----

      end subroutine s_unimongr

!-----7--------------------------------------------------------------7--

      end module m_unimongr
