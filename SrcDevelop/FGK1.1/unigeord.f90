!***********************************************************************
      module m_unigeord
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
!                   2007/04/11, 2007/05/14, 2007/05/21, 2007/08/24,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     generate the united geography file from the dumped files in
!     reductional entire domain.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comgrp
      use m_comkind
      use m_commath
      use m_commpi
      use m_currpe
      use m_destroy
      use m_getcname
      use m_getiname
      use m_getunit
      use m_inichar
      use m_outstd03
      use m_putunit
      use m_setcst2d

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: unigeord, s_unigeord

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface unigeord

        module procedure s_unigeord

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
      subroutine s_unigeord(fpexprim,fpcrsdir,fpncexp,fpnccrs,fpwlngth, &
     &                      fprmopt_uni,nx,ny,var,nio_uni,iodmp)
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

      integer, intent(in) :: nx
                       ! Composite model dimension in x direction

      integer, intent(in) :: ny
                       ! Composite model dimension in y direction

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

      integer njrec    ! Number of record for each i and j plane

      integer recmst   ! Tail of record number at each process
                       ! in outermost loops

      integer recj     ! Tail of record number at each process
                       ! in y direction

      integer recuni   ! Current record number of united file
      integer recdmp   ! Current record number of dumped file
      integer recmax   ! Maximum record number of dumped file

      integer sizuni   ! Record length of united file
      integer sizdmp   ! Record length of dumped file

      integer stat     ! Runtime status

      integer iouni    ! Unit number of united file

      integer, intent(inout) :: iodmp(0:nio_uni)
                       ! Unit number of dumped file

      real, intent(inout) :: var(2:nx-2,2:ny-2)
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

      njrec=nred*(ny-3)

! -----

!!!! Generate the united file from the dumped files.

! Initialize the character variable.

      call inichar(unifl)

! -----

! Get the unit number.

      call getunit(iouni)

! -----

! Open the united file.

      sizuni=(nx-3)*wlngth

      unifl(1:ncexp)=exprim(1:ncexp)

      write(unifl(ncexp+1:ncexp+20),'(a20)') 'geography.united.bin'

      open(iouni,iostat=stat,err=100,                                   &
     &     file=crsdir(1:nccrs)//unifl(1:ncexp+20),                     &
     &     status='new',access='direct',form='unformatted',             &
     &     recl=sizuni,action='write')

  100 if(stat.ne.0) then

        call destroy('unigeord',8,'stop',1,'              ',14,iouni,   &
     &               stat)

      end if

      call outstd03('unigeord',8,unifl,ncexp+20,iouni,1,0,0_i8)

! -----

!!! Perform uniting.

!! For the defined domain.

      recmax=0

      do myred=0,nred-1

! Set the parameters of parallelizing.

        call currpe('unite   ',5,'ijred')

! -----

! Get the unit number.

        if(grpxy(igrp+1,jgrp+1).ge.0) then

          call getunit(iodmp(0))

        end if

! -----

! Open the dumped file.

        if(grpxy(igrp+1,jgrp+1).ge.0) then

          sizdmp=(nx-3)*(ny-3)*wlngth

          dmpfl(1:ncexp)=exprim(1:ncexp)

          if(nsub.eq.1) then

            write(dmpfl(ncexp+1:ncexp+29),'(a13,i4.4,a12)')             &
     &                              'geography.grp',mygrp,'-sub0000.bin'

            if(rmopt_uni.eq.1) then

              open(iodmp(0),iostat=stat,err=110,                        &
     &             file=crsdir(1:nccrs)//dmpfl(1:ncexp+29),             &
     &             status='old',access='direct',form='unformatted',     &
     &             recl=sizdmp,action='readwrite')

            else

              open(iodmp(0),iostat=stat,err=110,                        &
     &             file=crsdir(1:nccrs)//dmpfl(1:ncexp+29),             &
     &             status='old',access='direct',form='unformatted',     &
     &             recl=sizdmp,action='read')

            end if

  110       if(stat.ne.0) then

              call destroy('unigeord',8,'cont',1,'              ',14,   &
     &                     iodmp(0),stat)

              go to 170

            end if

            call outstd03('unigeord',8,dmpfl,ncexp+29,iodmp(0),1,0,0_i8)

          else

            write(dmpfl(ncexp+1:ncexp+28),'(a13,i4.4,a11)')             &
     &                               'geography.grp',mygrp,'.united.bin'

            if(rmopt_uni.eq.1) then

              open(iodmp(0),iostat=stat,err=120,                        &
     &             file=crsdir(1:nccrs)//dmpfl(1:ncexp+28),             &
     &             status='old',access='direct',form='unformatted',     &
     &             recl=sizdmp,action='readwrite')

            else

              open(iodmp(0),iostat=stat,err=120,                        &
     &             file=crsdir(1:nccrs)//dmpfl(1:ncexp+28),             &
     &             status='old',access='direct',form='unformatted',     &
     &             recl=sizdmp,action='read')

            end if

  120       if(stat.ne.0) then

              call destroy('unigeord',8,'cont',1,'              ',14,   &
     &                     iodmp(0),stat)

              go to 170

            end if

            call outstd03('unigeord',8,dmpfl,ncexp+28,iodmp(0),1,0,0_i8)

          end if

        end if

! -----

! Read out the data from the dumped file and read in the data to the
! united file.

        if(grpxy(igrp+1,jgrp+1).ge.0) then

          recdmp=0

          recj=nired*jred*(ny-3)+ired+1

          do_inf_1: do

            recdmp=recdmp+1

            recmst=njrec*(recdmp-1)+recj

            read(iodmp(0),rec=recdmp,iostat=stat,err=130)               &
     &          ((var(i,j),i=2,nx-2),j=2,ny-2)

            do j=2,ny-2
              recuni=recmst+(j-2)*nired

              write(iouni,rec=recuni,iostat=stat,err=160)               &
     &             (var(i,j),i=2,nx-2)

            end do

          end do do_inf_1

  130     if(recdmp.eq.1) then

            call destroy('unigeord',8,'cont',3,'              ',14,     &
     &                   iodmp(0),stat)

            go to 170

          else

            if(recmax.lt.recdmp) then
              recmax=recdmp
            end if

            call outstd03('unigeord',8,dmpfl,108,iodmp(0),3,0,0_i8)

          end if

          call outstd03('unigeord',8,unifl,108,iouni,4,0,0_i8)

        end if

! -----

! Close the dumped file.

        if(grpxy(igrp+1,jgrp+1).ge.0) then

          if(rmopt_uni.eq.1) then

            close(iodmp(0),iostat=stat,err=140,status='delete')

          else

            close(iodmp(0),iostat=stat,err=140,status='keep')

          end if

  140     if(stat.ne.0) then

            call destroy('unigeord',8,'cont',2,'              ',14,     &
     &                   iodmp(0),stat)

            go to 170

          end if

          call outstd03('unigeord',8,unifl,108,iodmp(0),2,0,0_i8)

        end if

! -----

! Return the unit number.

        if(grpxy(igrp+1,jgrp+1).ge.0) then

          call putunit(iodmp(0))

        end if

! -----

      end do

!! -----

!! For the undefined domain.

      do myred=0,nred-1

! Set the parameters of parallelizing.

        call currpe('unite   ',5,'ijred')

! -----

! Read in the undefined data to the united file.

        if(grpxy(igrp+1,jgrp+1).lt.0) then

          recdmp=0

          recj=nired*jred*(ny-3)+ired+1

          do_inf_2: do

            recdmp=recdmp+1

            if(recdmp.ge.recmax) then

              exit do_inf_2

            end if

            recmst=njrec*(recdmp-1)+recj

            call setcst2d(2,nx-2,2,ny-2,lim35n,var)

            do j=2,ny-2
              recuni=recmst+(j-2)*nired

              write(iouni,rec=recuni,iostat=stat,err=160)               &
     &             (var(i,j),i=2,nx-2)

            end do

          end do do_inf_2

          call outstd03('unigeord',8,unifl,108,iouni,4,0,0_i8)

        end if

! -----

      end do

!! -----

!!! -----

! Close the united file.

      close(iouni,iostat=stat,err=150,status='keep')

  150 if(stat.ne.0) then

        call destroy('unigeord',8,'stop',2,'              ',14,iouni,   &
     &               stat)

      end if

      call outstd03('unigeord',8,unifl,108,iouni,2,0,0_i8)

! -----

! Return the unit number.

      call putunit(iouni)

! -----

!!!! -----

      return

! If error occured, call the procedure destroy.

  160 call destroy('unigeord',8,'cont',4,'              ',14,iouni,     &
     &             stat)

  170 close(iouni,iostat=stat,err=180,status='delete')

  180 if(stat.ne.0) then

        call destroy('unigeord',8,'stop',2,'              ',14,iouni,   &
     &               stat)

      end if

      call outstd03('unigeord',8,unifl,108,iouni,2,0,0_i8)

      call destroy('unigeord',8,'stop',1001,'              ',14,101,    &
     &             stat)

! -----

      end subroutine s_unigeord

!-----7--------------------------------------------------------------7--

      end module m_unigeord
