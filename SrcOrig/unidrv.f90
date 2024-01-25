!***********************************************************************
      module m_unidrv
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
!                   2007/04/11, 2007/05/21, 2007/06/27, 2007/07/30,
!                   2007/08/24, 2007/09/28, 2008/01/11, 2008/04/17,
!                   2008/05/02, 2008/08/25, 2008/10/10, 2008/12/11,
!                   2009/01/30, 2009/02/27, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     control the inferior procedures for uniting.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_comkind
      use m_commpi
      use m_currpe
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar
      use m_outllsw
      use m_outstd05
      use m_rdcheck
      use m_unidmpen
      use m_unidmpgr
      use m_unidmprd
      use m_unigeoen
      use m_unigeogr
      use m_unigeord
      use m_unimonen
      use m_unimongr
      use m_unimonrd
      use m_unistep

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: unidrv, s_unidrv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface unidrv

        module procedure s_unidrv

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic int
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_unidrv(fpfltyp_uni,fpiniopt,fpdmpmon,                &
     &                    fpuniopt_uni,fpugroup_uni,fpflitv_uni,        &
     &                    ni,nj,nk,tmp1,tmp2,tmp3,tmp4,                 &
     &                    ni_uni,nj_uni,var,nio_uni,iodmp)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpfltyp_uni
                       ! Formal parameter of unique index of fltyp_uni

      integer, intent(in) :: fpiniopt
                       ! Formal parameter of unique index of iniopt

      integer, intent(in) :: fpdmpmon
                       ! Formal parameter of unique index of dmpmon

      integer, intent(in) :: fpuniopt_uni
                       ! Formal parameter of unique index of uniopt_uni

      integer, intent(in) :: fpugroup_uni
                       ! Formal parameter of unique index of ugroup_uni

      integer, intent(in) :: fpflitv_uni
                       ! Formal parameter of unique index of flitv_uni

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: ni_uni
                       ! Model dimension of unite in x direction

      integer, intent(in) :: nj_uni
                       ! Model dimension of unite in y direction

      integer, intent(in) :: nio_uni
                       ! Maximum unit number of unite

! Internal shared variables

      character(len=108) fltyp_uni
                       ! Control flag of processed file type

      integer iniopt   ! Option for model initialization
      integer dmpmon   ! Option for monitor variables output

      integer uniopt_uni
                       ! Option for uniting process

      integer ugroup_uni
                       ! User specified and united group number
                       ! in entire domain

      integer nx       ! Composite model dimension in x direction
      integer ny       ! Composite model dimension in y direction

      integer(kind=i8) it
                       ! Index of main do loop

      integer(kind=i8) nstp0
                       ! Start index of main do loop

      integer(kind=i8) nstp1
                       ! End index of main do loop

      integer(kind=i8) ctime
                       ! Model current forecast time

      integer(kind=i8) fl103
                       ! 1000 x int(flitv_uni + 0.1)

      integer, intent(inout) :: iodmp(0:nio_uni)
                       ! Unit number of dumped file

      real flitv_uni   ! Time interval of processed file

      real, intent(inout) :: tmp1(1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp4(1:(nj-3)*njgrp*njsub)
                       ! Temporary array

      real, intent(inout) :: var(2:ni_uni-2,2:nj_uni-2)
                       ! Optional variable in dumped file

!-----7--------------------------------------------------------------7--

! Initialize character variable.

      call inichar(fltyp_uni)

! -----

! Get the required namelist variables.

      call getcname(fpfltyp_uni,fltyp_uni)
      call getiname(fpiniopt,iniopt)
      call getiname(fpdmpmon,dmpmon)
      call getiname(fpuniopt_uni,uniopt_uni)
      call getiname(fpugroup_uni,ugroup_uni)
      call getrname(fpflitv_uni,flitv_uni)

! -----

! Set the common used variables.

      nx=(ni-3)*nisub+3
      ny=(nj-3)*njsub+3

      fl103=1000_i8*int(flitv_uni+.1e0,i8)

! -----

!!!! Generate the united geography file from the dumped files.

      if(fltyp_uni(1:3).eq.'all'.or.fltyp_uni(1:3).eq.'geo') then

!!! Generate the united file from the dumped files.

! Check the geography data checking file.

        if(iniopt.ne.12) then

          if(abs(mod(uniopt_uni,10)).le.4.or.                           &
     &      (ngrp.eq.1.and.abs(mod(uniopt_uni,10)).ge.5)) then

            call rdcheck(idexprim,idcrsdir,idncexp,idnccrs,'geo',0_i8)

          end if

        end if

! -----

!! Generate the united file from the dumped files in group domain.

! In the case of uniting all files to group domain files.

        if(abs(mod(uniopt_uni,10)).eq.1.or.abs(mod(uniopt_uni,10)).eq.2 &
     &    .or.abs(mod(uniopt_uni,10)).eq.5) then

          if(nsub.gt.1) then

            do mysrl=0,nsrl-1

              call currpe('unite   ',5,'mygrp')

              call unigeogr(idexprim,idcrsdir,idncexp,idnccrs,idwlngth, &
     &                      idrmopt_uni,nx,ny,ni,nj,ni_uni,nj_uni,var,  &
     &                      nio_uni,iodmp)

            end do

          end if

! -----

! In the case of uniting specified files to group domain file.

        else if(abs(mod(uniopt_uni,10)).eq.6) then

          if(nsub.gt.1) then

            mygrp=ugroup_uni

            call unigeogr(idexprim,idcrsdir,idncexp,idnccrs,idwlngth,   &
     &                    idrmopt_uni,nx,ny,ni,nj,ni_uni,nj_uni,var,    &
     &                    nio_uni,iodmp)

          end if

        end if

! -----

!! -----

! Generate the united file from the dumped files in entire domain.

        if(abs(mod(uniopt_uni,10)).eq.1                                 &
     &    .or.abs(mod(uniopt_uni,10)).eq.3) then

          call unigeoen(idexprim,idcrsdir,idncexp,idnccrs,idwlngth,     &
     &                  idrmopt_uni,nx,ny,var,nio_uni,iodmp)

        end if

! -----

! Generate the united file from the dumped files in reductional entire
! domain.

        if(abs(mod(uniopt_uni,10)).eq.2                                 &
     &    .or.abs(mod(uniopt_uni,10)).eq.4) then

          call unigeord(idexprim,idcrsdir,idncexp,idnccrs,idwlngth,     &
     &                  idrmopt_uni,nx,ny,var,nio_uni,iodmp)

        end if

! -----

!!! -----

! Read in the latitude and the longitude at south-west corner to
! standard i/o.

        if(abs(mod(uniopt_uni,10)).le.4.or.                             &
     &    (ngrp.eq.1.and.abs(mod(uniopt_uni,10)).ge.5)) then

          call outllsw(idfltyp_uni,iddmpmon,iduniopt_uni,               &
     &                 'geo',1_i8,1_i8,ni,nj,nk,tmp1,tmp2,tmp3,tmp4)

        end if

! -----

      end if

!!!! -----

!!!! Generate the united history file from the dumped files.

      if(fltyp_uni(1:3).eq.'all'                                        &
     &  .or.fltyp_uni(1:3).eq.'dmp'.or.fltyp_uni(1:3).eq.'mon') then

! Read in the message to standard i/o.

        call outstd05(0)

! -----

! Calculate the number of steps of the main do loop.

        call unistep(idfltyp_uni,idflitv_uni,idstime,idetime,           &
     &               nstp0,nstp1)

! -----

!!! Generate the united file from the dumped files.

        do it=nstp0,nstp1

! Calculate the current forecast time.

          ctime=fl103*(it-1_i8)

! -----

! Check the history data checking file.

          if(iniopt.ne.12) then

            if(abs(mod(uniopt_uni,10)).le.4.or.                         &
     &        (ngrp.eq.1.and.abs(mod(uniopt_uni,10)).ge.5)) then

              if(dmpmon.eq.0) then

                call rdcheck(idexprim,idcrsdir,idncexp,idnccrs,'dmp',   &
     &                       ctime)

              else

                if(fltyp_uni(1:3).ne.'mon') then

                  call rdcheck(idexprim,idcrsdir,idncexp,idnccrs,'dmp', &
     &                         ctime)

                end if

                if(fltyp_uni(1:3).ne.'dmp') then

                  call rdcheck(idexprim,idcrsdir,idncexp,idnccrs,'mon', &
     &                         ctime)

                end if

              end if

            end if

          end if

! -----

!! Generate the united file from the dumped files in group domain.

! In the case of uniting all files to group domain files.

          if(abs(mod(uniopt_uni,10)).eq.1                               &
     &      .or.abs(mod(uniopt_uni,10)).eq.2                            &
     &      .or.abs(mod(uniopt_uni,10)).eq.5) then

            if(nsub.gt.1) then

              if(dmpmon.eq.0) then

                do mysrl=0,nsrl-1

                  call currpe('unite   ',5,'mygrp')

                  call unidmpgr(idexprim,idcrsdir,idncexp,idnccrs,      &
     &                          idwlngth,idrmopt_uni,ctime,nx,ny,       &
     &                          ni,nj,ni_uni,nj_uni,var,nio_uni,iodmp)

                end do

              else

                if(fltyp_uni(1:3).ne.'mon') then

                  do mysrl=0,nsrl-1

                    call currpe('unite   ',5,'mygrp')

                    call unidmpgr(idexprim,idcrsdir,idncexp,idnccrs,    &
     &                            idwlngth,idrmopt_uni,ctime,nx,ny,     &
     &                            ni,nj,ni_uni,nj_uni,var,nio_uni,iodmp)

                  end do

                end if

                if(fltyp_uni(1:3).ne.'dmp') then

                  do mysrl=0,nsrl-1

                    call currpe('unite   ',5,'mygrp')

                    call unimongr(idexprim,idcrsdir,idncexp,idnccrs,    &
     &                            idwlngth,idrmopt_uni,ctime,nx,ny,     &
     &                            ni,nj,ni_uni,nj_uni,var,nio_uni,iodmp)

                  end do

                end if

              end if

            end if

! -----

! In the case of uniting specified files to group domain file.

          else if(abs(mod(uniopt_uni,10)).eq.6) then

            if(nsub.gt.1) then

              mygrp=ugroup_uni

              if(dmpmon.eq.0) then

                call unidmpgr(idexprim,idcrsdir,idncexp,idnccrs,        &
     &                        idwlngth,idrmopt_uni,ctime,nx,ny,         &
     &                        ni,nj,ni_uni,nj_uni,var,nio_uni,iodmp)

              else

                if(fltyp_uni(1:3).ne.'mon') then

                  call unidmpgr(idexprim,idcrsdir,idncexp,idnccrs,      &
     &                          idwlngth,idrmopt_uni,ctime,nx,ny,       &
     &                          ni,nj,ni_uni,nj_uni,var,nio_uni,iodmp)

                end if

                if(fltyp_uni(1:3).ne.'dmp') then

                  call unimongr(idexprim,idcrsdir,idncexp,idnccrs,      &
     &                          idwlngth,idrmopt_uni,ctime,nx,ny,       &
     &                          ni,nj,ni_uni,nj_uni,var,nio_uni,iodmp)

                end if

              end if

            end if

          end if

! -----

!! -----

! Generate the united file from the dumped files in entire domain.

          if(abs(mod(uniopt_uni,10)).eq.1                               &
     &      .or.abs(mod(uniopt_uni,10)).eq.3) then

            if(dmpmon.eq.0) then

              call unidmpen(idexprim,idcrsdir,idncexp,idnccrs,idwlngth, &
     &                      idrmopt_uni,ctime,nx,ny,var,nio_uni,iodmp)

            else

              if(fltyp_uni(1:3).ne.'mon') then

                call unidmpen(idexprim,idcrsdir,idncexp,idnccrs,        &
     &                        idwlngth,idrmopt_uni,ctime,nx,ny,         &
     &                        var,nio_uni,iodmp)

              end if

              if(fltyp_uni(1:3).ne.'dmp') then

                call unimonen(idexprim,idcrsdir,idncexp,idnccrs,        &
     &                        idwlngth,idrmopt_uni,ctime,nx,ny,         &
     &                        var,nio_uni,iodmp)

              end if

            end if

          end if

! -----

! Generate the united file from the dumped files in reductional entire
! domain.

          if(abs(mod(uniopt_uni,10)).eq.2                               &
     &      .or.abs(mod(uniopt_uni,10)).eq.4) then

            if(dmpmon.eq.0) then

              call unidmprd(idexprim,idcrsdir,idncexp,idnccrs,idwlngth, &
     &                      idrmopt_uni,ctime,nx,ny,var,nio_uni,iodmp)

            else

              if(fltyp_uni(1:3).ne.'mon') then

                call unidmprd(idexprim,idcrsdir,idncexp,idnccrs,        &
     &                        idwlngth,idrmopt_uni,ctime,nx,ny,         &
     &                        var,nio_uni,iodmp)

              end if

              if(fltyp_uni(1:3).ne.'dmp') then

                call unimonrd(idexprim,idcrsdir,idncexp,idnccrs,        &
     &                        idwlngth,idrmopt_uni,ctime,nx,ny,         &
     &                        var,nio_uni,iodmp)

              end if

            end if

          end if

! -----

! Read in the message to standard i/o.

          call outstd05(0)

! -----

        end do

!!! -----

! Read in the latitude and the longitude at south-west corner to
! standard i/o.

        if(abs(mod(uniopt_uni,10)).le.4.or.                             &
     &    (ngrp.eq.1.and.abs(mod(uniopt_uni,10)).ge.5)) then

          if(dmpmon.eq.0) then

            call outllsw(idfltyp_uni,iddmpmon,iduniopt_uni,             &
     &                   'dmp',nstp0,nstp1,ni,nj,nk,tmp1,tmp2,tmp3,tmp4)

          else

            if(fltyp_uni(1:3).ne.'mon') then

              call outllsw(idfltyp_uni,iddmpmon,iduniopt_uni,           &
     &                   'dmp',nstp0,nstp1,ni,nj,nk,tmp1,tmp2,tmp3,tmp4)

            end if

            if(fltyp_uni(1:3).ne.'dmp') then

              call outllsw(idfltyp_uni,iddmpmon,iduniopt_uni,           &
     &                   'mon',nstp0,nstp1,ni,nj,nk,tmp1,tmp2,tmp3,tmp4)

            end if

          end if

        end if

! -----

      end if

!!!! -----

      end subroutine s_unidrv

!-----7--------------------------------------------------------------7--

      end module m_unidrv
