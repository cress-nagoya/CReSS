!***********************************************************************
      module m_bcini
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2005/02/10
!     Modification: 2006/01/10, 2006/09/21, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/03/23, 2011/08/18, 2011/09/22, 2011/11/10,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the boundary conditions for GPV data at forecast start time.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bcycle
      use m_bcyclex
      use m_combuf
      use m_comindx
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getcname
      use m_getiname
      use m_inichar
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: bcini, s_bcini

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcini

        module procedure s_bcini

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_bcini(fpgpvvar,fpwbc,fpebc,fpexbopt,                 &
     &                   fpcphopt,fphaiopt,ni,nj,nk,nqw,nqi,            &
     &                   ubr,vbr,pbr,ptbr,qvbr,ugpv,vgpv,wgpv,          &
     &                   ppgpv,ptpgpv,qvgpv,qwgpv,qigpv)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

! Input and output variables

      real, intent(inout) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(inout) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(inout) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(inout) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(inout) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(inout) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(inout) :: ppgpv(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation of GPV data
                       ! at marked time

      real, intent(inout) :: ptpgpv(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation of GPV data
                       ! at marked time

      real, intent(inout) :: qvgpv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio of GPV data
                       ! at marked time

      real, intent(inout) :: qwgpv(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of GPV data at marked time

      real, intent(inout) :: qigpv(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of GPV data at marked time

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer exbopt   ! Option for external boundary forcing
      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes

      integer ib       ! Exchanging variables number
      integer nb       ! Number of exchanging variables
      integer nbb      ! Number of exchanging base state variables

      integer ni_sub   ! Substitute for ni
      integer nj_sub   ! Substitute for nj

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(gpvvar)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpexbopt,exbopt)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)

! -----

! Set the substituted variables.

      ni_sub=ni
      nj_sub=nj

! -----

! Count the number of variables.

      if(exbopt.eq.0                                                    &
     &  .or.(exbopt.ge.1.and.abs(wbc).eq.1.and.abs(ebc).eq.1)) then

        if(gpvvar(2:2).eq.'o') then
          nb=9
          nbb=5
        else
          nb=8
          nbb=4
        end if

        if(gpvvar(1:1).eq.'o') then
          nb=nb+1
        end if

        if(gpvvar(2:2).eq.'o') then
          nb=nb+1
        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then
              nb=nb+1
            end if

            if(gpvvar(4:4).eq.'o') then
              nb=nb+1
            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then
              nb=nb+1
            end if

            if(gpvvar(6:6).eq.'o') then
              nb=nb+1
            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then
              nb=nb+1
            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then
                nb=nb+1
              end if

            end if

          end if

        end if

      end if

! -----

!!!! Set the lateral boundary conditions.

      if(exbopt.eq.0) then

!!! Exchange the value between first halo regions.

!! Exchange the value between sub domain.

! In x direction.

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ubr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,pbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qvbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ugpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vgpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ppgpv,      &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,sbuf)

        else

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ubr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,pbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ugpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vgpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ppgpv,      &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,sbuf)

        end if

        if(gpvvar(1:1).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,wgpv,       &
     &                    ib,nb,sbuf)

        end if

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qvgpv,      &
     &                    ib,nb,sbuf)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qwgpv(0,0,1,1),ib,nb,sbuf)

            end if

            if(gpvvar(4:4).eq.'o') then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qwgpv(0,0,1,2),ib,nb,sbuf)

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,1),ib,nb,sbuf)

            end if

            if(gpvvar(6:6).eq.'o') then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,2),ib,nb,sbuf)

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,3),ib,nb,sbuf)

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,      &
     &                          qigpv(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          end if

        end if

        call s_shiftsx(idwbc,idebc,'bnd',nj,nk,nb,sbuf,rbuf)

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ubr,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,pbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qvbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ugpv,     &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vgpv,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ppgpv,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,rbuf)

        else

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ubr,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,pbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ugpv,     &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vgpv,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ppgpv,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,rbuf)

        end if

        if(gpvvar(1:1).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,wgpv,       &
     &                    ib,nb,rbuf)

        end if

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qvgpv,      &
     &                    ib,nb,rbuf)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qwgpv(0,0,1,1),ib,nb,rbuf)

            end if

            if(gpvvar(4:4).eq.'o') then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qwgpv(0,0,1,2),ib,nb,rbuf)

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,1),ib,nb,rbuf)

            end if

            if(gpvvar(6:6).eq.'o') then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,2),ib,nb,rbuf)

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,3),ib,nb,rbuf)

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,      &
     &                          qigpv(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          end if

        end if

! -----

! In y direction.

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ubr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,vbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,pbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ptbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,qvbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ugpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,vgpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ppgpv,      &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,sbuf)

        else

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ubr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,vbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,pbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ptbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ugpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,vgpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ppgpv,      &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,sbuf)

        end if

        if(gpvvar(1:1).eq.'o') then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,wgpv,       &
     &                    ib,nb,sbuf)

        end if

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,qvgpv,      &
     &                    ib,nb,sbuf)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,        &
     &                        qwgpv(0,0,1,1),ib,nb,sbuf)

            end if

            if(gpvvar(4:4).eq.'o') then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,        &
     &                        qwgpv(0,0,1,2),ib,nb,sbuf)

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,1),ib,nb,sbuf)

            end if

            if(gpvvar(6:6).eq.'o') then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,2),ib,nb,sbuf)

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              ib=ib+1

              call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,3),ib,nb,sbuf)

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,      &
     &                          qigpv(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          end if

        end if

        call s_shiftsy(idsbc,idnbc,'bnd',ni,nk,nb,sbuf,rbuf)

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ubr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,vbr,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,pbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ptbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,qvbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ugpv,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,vgpv,     &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ppgpv,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,rbuf)

        else

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ubr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,vbr,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,pbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ptbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ugpv,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,vgpv,     &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ppgpv,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,rbuf)

        end if

        if(gpvvar(1:1).eq.'o') then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,wgpv,       &
     &                    ib,nb,rbuf)

        end if

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,qvgpv,      &
     &                    ib,nb,rbuf)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,        &
     &                        qwgpv(0,0,1,1),ib,nb,rbuf)

            end if

            if(gpvvar(4:4).eq.'o') then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,        &
     &                        qwgpv(0,0,1,2),ib,nb,rbuf)

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,1),ib,nb,rbuf)

            end if

            if(gpvvar(6:6).eq.'o') then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,2),ib,nb,rbuf)

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              ib=ib+1

              call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,3),ib,nb,rbuf)

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,      &
     &                          qigpv(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          end if

        end if

! -----

!! -----

!! Exchange the value betweein group domain.

! In x direction.

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ubr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,pbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qvbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ugpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vgpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ppgpv,      &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,sbuf)

        else

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ubr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,pbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ugpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vgpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ppgpv,      &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,sbuf)

        end if

        if(gpvvar(1:1).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,wgpv,       &
     &                    ib,nb,sbuf)

        end if

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qvgpv,      &
     &                    ib,nb,sbuf)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qwgpv(0,0,1,1),ib,nb,sbuf)

            end if

            if(gpvvar(4:4).eq.'o') then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qwgpv(0,0,1,2),ib,nb,sbuf)

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,1),ib,nb,sbuf)

            end if

            if(gpvvar(6:6).eq.'o') then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,2),ib,nb,sbuf)

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,3),ib,nb,sbuf)

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,      &
     &                          qigpv(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          end if

        end if

        call s_shiftgx(idwbc,idebc,'bnd',nj,nk,nb,sbuf,rbuf)

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ubr,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,pbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qvbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ugpv,     &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vgpv,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ppgpv,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,rbuf)

        else

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ubr,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,pbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ugpv,     &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vgpv,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ppgpv,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,rbuf)

        end if

        if(gpvvar(1:1).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,wgpv,       &
     &                    ib,nb,rbuf)

        end if

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qvgpv,      &
     &                    ib,nb,rbuf)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qwgpv(0,0,1,1),ib,nb,rbuf)

            end if

            if(gpvvar(4:4).eq.'o') then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qwgpv(0,0,1,2),ib,nb,rbuf)

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,1),ib,nb,rbuf)

            end if

            if(gpvvar(6:6).eq.'o') then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,2),ib,nb,rbuf)

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,3),ib,nb,rbuf)

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,      &
     &                          qigpv(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          end if

        end if

! -----

! In y direction.

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ubr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,vbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,pbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ptbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,qvbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ugpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,vgpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ppgpv,      &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,sbuf)

        else

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ubr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,vbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,pbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ptbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ugpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',3,nj-2,ni,nj,nk,vgpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ppgpv,      &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,sbuf)

        end if

        if(gpvvar(1:1).eq.'o') then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,wgpv,       &
     &                    ib,nb,sbuf)

        end if

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,qvgpv,      &
     &                    ib,nb,sbuf)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,        &
     &                        qwgpv(0,0,1,1),ib,nb,sbuf)

            end if

            if(gpvvar(4:4).eq.'o') then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,        &
     &                        qwgpv(0,0,1,2),ib,nb,sbuf)

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,1),ib,nb,sbuf)

            end if

            if(gpvvar(6:6).eq.'o') then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,2),ib,nb,sbuf)

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              ib=ib+1

              call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,3),ib,nb,sbuf)

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,      &
     &                          qigpv(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          end if

        end if

        call s_shiftgy(idsbc,idnbc,'bnd',ni,nk,nb,sbuf,rbuf)

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ubr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,vbr,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,pbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ptbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,qvbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ugpv,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,vgpv,     &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ppgpv,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,rbuf)

        else

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ubr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,vbr,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,pbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ptbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ugpv,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj_sub,ni,nj,nk,vgpv,     &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ppgpv,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,rbuf)

        end if

        if(gpvvar(1:1).eq.'o') then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,wgpv,       &
     &                    ib,nb,rbuf)

        end if

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,qvgpv,      &
     &                    ib,nb,rbuf)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,        &
     &                        qwgpv(0,0,1,1),ib,nb,rbuf)

            end if

            if(gpvvar(4:4).eq.'o') then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,        &
     &                        qwgpv(0,0,1,2),ib,nb,rbuf)

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,1),ib,nb,rbuf)

            end if

            if(gpvvar(6:6).eq.'o') then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,2),ib,nb,rbuf)

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              ib=ib+1

              call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,3),ib,nb,rbuf)

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,      &
     &                          qigpv(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          end if

        end if

! -----

!! -----

! Set the periodic boundary conditions.

        if(gpvvar(2:2).eq.'o') then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,1,ni-2,ni_sub,2,1,nj-2,nj-1,ni,nj,nk,ubr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,3,1,nj-2,nj_sub,ni,nj,nk,vbr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,pbr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,ptbr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,qvbr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,1,ni-2,ni_sub,2,1,nj-2,nj-1,ni,nj,nk,ugpv)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,3,1,nj-2,nj_sub,ni,nj,nk,vgpv)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,ppgpv)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,ptpgpv)

        else

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,1,ni-2,ni_sub,2,1,nj-2,nj-1,ni,nj,nk,ubr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,3,1,nj-2,nj_sub,ni,nj,nk,vbr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,pbr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,ptbr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,1,ni-2,ni_sub,2,1,nj-2,nj-1,ni,nj,nk,ugpv)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,3,1,nj-2,nj_sub,ni,nj,nk,vgpv)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,ppgpv)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,ptpgpv)

        end if

        if(gpvvar(1:1).eq.'o') then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,wgpv)

        end if

        if(gpvvar(2:2).eq.'o') then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,qvgpv)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qwgpv(0,0,1,1))

            end if

            if(gpvvar(4:4).eq.'o') then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qwgpv(0,0,1,2))

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qigpv(0,0,1,1))

            end if

            if(gpvvar(6:6).eq.'o') then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qigpv(0,0,1,2))

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              call s_bcycle(idwbc,idebc,idsbc,idnbc,                    &
     &                      2,1,ni-2,ni-1,2,1,nj-2,nj-1,                &
     &                      ni,nj,nk,qigpv(0,0,1,3))

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                call s_bcycle(idwbc,idebc,idsbc,idnbc,                  &
     &                        2,1,ni-2,ni-1,2,1,nj-2,nj-1,              &
     &                        ni,nj,nk,qigpv(0,0,1,4))

              end if

            end if

          end if

        end if

! -----

!!! -----

!!! Exchange the value between second halo regions.

!! Exchange the value between sub domain.

! In x direction.

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',4,ni-3,ni,nj,nk,ubr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,vbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,pbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,ptbr,       &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,qvbr,       &
     &                    ib,nbb,sbuf)

        else

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',4,ni-3,ni,nj,nk,ubr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,vbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,pbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,ptbr,       &
     &                    ib,nbb,sbuf)

        end if

        call s_shiftsx(idwbc,idebc,'bnd',nj,nk,nbb,sbuf,rbuf)

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni+1,ni,nj,nk,ubr,        &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,vbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,pbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,ptbr,     &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,qvbr,     &
     &                    ib,nbb,rbuf)

        else

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni+1,ni,nj,nk,ubr,        &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,vbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,pbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,ptbr,     &
     &                    ib,nbb,rbuf)

        end if

! -----

! In y direction.

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',3,nj-3,ni,nj,nk,ubr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',4,nj-3,ni,nj,nk,vbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',3,nj-3,ni,nj,nk,pbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',3,nj-3,ni,nj,nk,ptbr,       &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',3,nj-3,ni,nj,nk,qvbr,       &
     &                    ib,nbb,sbuf)

        else

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',3,nj-3,ni,nj,nk,ubr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',4,nj-3,ni,nj,nk,vbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',3,nj-3,ni,nj,nk,pbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsy(idsbc,idnbc,'bnd',3,nj-3,ni,nj,nk,ptbr,       &
     &                    ib,nbb,sbuf)

        end if

        call s_shiftsy(idsbc,idnbc,'bnd',ni,nk,nbb,sbuf,rbuf)

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',0,nj_sub,ni,nj,nk,ubr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',0,nj+1,ni,nj,nk,vbr,        &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',0,nj_sub,ni,nj,nk,pbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',0,nj_sub,ni,nj,nk,ptbr,     &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',0,nj_sub,ni,nj,nk,qvbr,     &
     &                    ib,nbb,rbuf)

        else

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',0,nj_sub,ni,nj,nk,ubr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',0,nj+1,ni,nj,nk,vbr,        &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',0,nj_sub,ni,nj,nk,pbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsy(idsbc,idnbc,'bnd',0,nj_sub,ni,nj,nk,ptbr,     &
     &                    ib,nbb,rbuf)

        end if

! -----

!! -----

!! Exchange the value betweein group domain.

! In x direction.

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',4,ni-3,ni,nj,nk,ubr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,vbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,pbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,ptbr,       &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,qvbr,       &
     &                    ib,nbb,sbuf)

        else

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',4,ni-3,ni,nj,nk,ubr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,vbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,pbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,ptbr,       &
     &                    ib,nbb,sbuf)

        end if

        call s_shiftgx(idwbc,idebc,'bnd',nj,nk,nbb,sbuf,rbuf)

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni+1,ni,nj,nk,ubr,        &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,vbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,pbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,ptbr,     &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,qvbr,     &
     &                    ib,nbb,rbuf)

        else

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni+1,ni,nj,nk,ubr,        &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,vbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,pbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,ptbr,     &
     &                    ib,nbb,rbuf)

        end if

! -----

! In y direction.

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',3,nj-3,ni,nj,nk,ubr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',4,nj-3,ni,nj,nk,vbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',3,nj-3,ni,nj,nk,pbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',3,nj-3,ni,nj,nk,ptbr,       &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',3,nj-3,ni,nj,nk,qvbr,       &
     &                    ib,nbb,sbuf)

        else

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',3,nj-3,ni,nj,nk,ubr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',4,nj-3,ni,nj,nk,vbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',3,nj-3,ni,nj,nk,pbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgy(idsbc,idnbc,'bnd',3,nj-3,ni,nj,nk,ptbr,       &
     &                    ib,nbb,sbuf)

        end if

        call s_shiftgy(idsbc,idnbc,'bnd',ni,nk,nbb,sbuf,rbuf)

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',0,nj_sub,ni,nj,nk,ubr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',0,nj+1,ni,nj,nk,vbr,        &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',0,nj_sub,ni,nj,nk,pbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',0,nj_sub,ni,nj,nk,ptbr,     &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',0,nj_sub,ni,nj,nk,qvbr,     &
     &                    ib,nbb,rbuf)

        else

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',0,nj_sub,ni,nj,nk,ubr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',0,nj+1,ni,nj,nk,vbr,        &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',0,nj_sub,ni,nj,nk,pbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgy(idsbc,idnbc,'bnd',0,nj_sub,ni,nj,nk,ptbr,     &
     &                    ib,nbb,rbuf)

        end if

! -----

!! -----

! Set the periodic boundary conditions.

        if(gpvvar(2:2).eq.'o') then

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                4,0,ni-3,ni+1,3,0,nj-3,nj_sub,ni,nj,nk,ubr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,0,ni-3,ni_sub,4,0,nj-3,nj+1,ni,nj,nk,vbr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,nk,pbr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,nk,ptbr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,nk,qvbr)

        else

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                4,0,ni-3,ni+1,3,0,nj-3,nj_sub,ni,nj,nk,ubr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,0,ni-3,ni_sub,4,0,nj-3,nj+1,ni,nj,nk,vbr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,nk,pbr)

          call bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                3,0,ni-3,ni_sub,3,0,nj-3,nj_sub,ni,nj,nk,ptbr)

        end if

! -----

!!! -----

      end if

!!!! -----

!!! Set the periodic boundary conditions in x direction.

      if(exbopt.ge.1.and.abs(wbc).eq.1.and.abs(ebc).eq.1) then

!! Exchange the value between first halo regions.

! Exchange the value between sub domain.

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ubr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,pbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qvbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ugpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vgpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ppgpv,      &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,sbuf)

        else

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ubr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,pbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ugpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vgpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ppgpv,      &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,sbuf)

        end if

        if(gpvvar(1:1).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,wgpv,       &
     &                    ib,nb,sbuf)

        end if

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qvgpv,      &
     &                    ib,nb,sbuf)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qwgpv(0,0,1,1),ib,nb,sbuf)

            end if

            if(gpvvar(4:4).eq.'o') then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qwgpv(0,0,1,2),ib,nb,sbuf)

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,1),ib,nb,sbuf)

            end if

            if(gpvvar(6:6).eq.'o') then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,2),ib,nb,sbuf)

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              ib=ib+1

              call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,3),ib,nb,sbuf)

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,      &
     &                          qigpv(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          end if

        end if

        call s_shiftsx(idwbc,idebc,'bnd',nj,nk,nb,sbuf,rbuf)

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ubr,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,pbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qvbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ugpv,     &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vgpv,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ppgpv,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,rbuf)

        else

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ubr,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,pbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ugpv,     &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vgpv,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ppgpv,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,rbuf)

        end if

        if(gpvvar(1:1).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,wgpv,       &
     &                    ib,nb,rbuf)

        end if

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qvgpv,      &
     &                    ib,nb,rbuf)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qwgpv(0,0,1,1),ib,nb,rbuf)

            end if

            if(gpvvar(4:4).eq.'o') then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qwgpv(0,0,1,2),ib,nb,rbuf)

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,1),ib,nb,rbuf)

            end if

            if(gpvvar(6:6).eq.'o') then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,2),ib,nb,rbuf)

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              ib=ib+1

              call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,3),ib,nb,rbuf)

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,      &
     &                          qigpv(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          end if

        end if

! -----

! Exchange the value betweein group domain.

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ubr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,pbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qvbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ugpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vgpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ppgpv,      &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,sbuf)

        else

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ubr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,pbr,        &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptbr,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-2,ni,nj,nk,ugpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,vgpv,       &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ppgpv,      &
     &                    ib,nb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,sbuf)

        end if

        if(gpvvar(1:1).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,wgpv,       &
     &                    ib,nb,sbuf)

        end if

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qvgpv,      &
     &                    ib,nb,sbuf)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qwgpv(0,0,1,1),ib,nb,sbuf)

            end if

            if(gpvvar(4:4).eq.'o') then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qwgpv(0,0,1,2),ib,nb,sbuf)

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,1),ib,nb,sbuf)

            end if

            if(gpvvar(6:6).eq.'o') then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,2),ib,nb,sbuf)

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              ib=ib+1

              call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,        &
     &                        qigpv(0,0,1,3),ib,nb,sbuf)

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                ib=ib+1

                call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,      &
     &                          qigpv(0,0,1,4),ib,nb,sbuf)

              end if

            end if

          end if

        end if

        call s_shiftgx(idwbc,idebc,'bnd',nj,nk,nb,sbuf,rbuf)

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ubr,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,pbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qvbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ugpv,     &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vgpv,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ppgpv,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,rbuf)

        else

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ubr,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,pbr,        &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptbr,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni_sub,ni,nj,nk,ugpv,     &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,vgpv,       &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ppgpv,      &
     &                    ib,nb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptpgpv,     &
     &                    ib,nb,rbuf)

        end if

        if(gpvvar(1:1).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,wgpv,       &
     &                    ib,nb,rbuf)

        end if

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qvgpv,      &
     &                    ib,nb,rbuf)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qwgpv(0,0,1,1),ib,nb,rbuf)

            end if

            if(gpvvar(4:4).eq.'o') then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qwgpv(0,0,1,2),ib,nb,rbuf)

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,1),ib,nb,rbuf)

            end if

            if(gpvvar(6:6).eq.'o') then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,2),ib,nb,rbuf)

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              ib=ib+1

              call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,        &
     &                        qigpv(0,0,1,3),ib,nb,rbuf)

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                ib=ib+1

                call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,      &
     &                          qigpv(0,0,1,4),ib,nb,rbuf)

              end if

            end if

          end if

        end if

! -----

! Set the periodic boundary conditions.

        if(gpvvar(2:2).eq.'o') then

          call bcyclex(idwbc,idebc,3,1,ni-2,ni_sub,ni,nj,nk,ubr)
          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,vbr)
          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,pbr)
          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,ptbr)
          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,qvbr)

          call bcyclex(idwbc,idebc,3,1,ni-2,ni_sub,ni,nj,nk,ugpv)
          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,vgpv)
          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,ppgpv)
          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,ptpgpv)

        else

          call bcyclex(idwbc,idebc,3,1,ni-2,ni_sub,ni,nj,nk,ubr)
          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,vbr)
          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,pbr)
          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,ptbr)

          call bcyclex(idwbc,idebc,3,1,ni-2,ni_sub,ni,nj,nk,ugpv)
          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,vgpv)
          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,ppgpv)
          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,ptpgpv)

        end if

        if(gpvvar(1:1).eq.'o') then

          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,wgpv)

        end if

        if(gpvvar(2:2).eq.'o') then

          call bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,nk,qvgpv)

        end if

        if(abs(cphopt).lt.10) then

          if(abs(cphopt).ge.1) then

            if(gpvvar(3:3).eq.'o') then

              call s_bcyclex(idwbc,idebc,                               &
     &                       2,1,ni-2,ni-1,ni,nj,nk,qwgpv(0,0,1,1))

            end if

            if(gpvvar(4:4).eq.'o') then

              call s_bcyclex(idwbc,idebc,                               &
     &                       2,1,ni-2,ni-1,ni,nj,nk,qwgpv(0,0,1,2))

            end if

          end if

          if(abs(cphopt).ge.2) then

            if(gpvvar(5:5).eq.'o') then

              call s_bcyclex(idwbc,idebc,                               &
     &                       2,1,ni-2,ni-1,ni,nj,nk,qigpv(0,0,1,1))

            end if

            if(gpvvar(6:6).eq.'o') then

              call s_bcyclex(idwbc,idebc,                               &
     &                       2,1,ni-2,ni-1,ni,nj,nk,qigpv(0,0,1,2))

            end if

            if(gpvvar(7:7).eq.'o'.or.gpvvar(8:8).eq.'o') then

              call s_bcyclex(idwbc,idebc,                               &
     &                       2,1,ni-2,ni-1,ni,nj,nk,qigpv(0,0,1,3))

            end if

            if(gpvvar(8:8).eq.'o') then

              if(haiopt.eq.1) then

                call s_bcyclex(idwbc,idebc,                             &
     &                         2,1,ni-2,ni-1,ni,nj,nk,qigpv(0,0,1,4))

              end if

            end if

          end if

        end if

! -----

!! -----

!! Exchange the value between second halo regions.

! Exchange the value between sub domain.

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',4,ni-3,ni,nj,nk,ubr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,vbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,pbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,ptbr,       &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,qvbr,       &
     &                    ib,nbb,sbuf)

        else

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',4,ni-3,ni,nj,nk,ubr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,vbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,pbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufsx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,ptbr,       &
     &                    ib,nbb,sbuf)

        end if

        call s_shiftsx(idwbc,idebc,'bnd',nj,nk,nbb,sbuf,rbuf)

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni+1,ni,nj,nk,ubr,        &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,vbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,pbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,ptbr,     &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,qvbr,     &
     &                    ib,nbb,rbuf)

        else

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni+1,ni,nj,nk,ubr,        &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,vbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,pbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufsx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,ptbr,     &
     &                    ib,nbb,rbuf)

        end if

! -----

! Exchange the value betweein group domain.

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',4,ni-3,ni,nj,nk,ubr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,vbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,pbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,ptbr,       &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,qvbr,       &
     &                    ib,nbb,sbuf)

        else

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',4,ni-3,ni,nj,nk,ubr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,vbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,pbr,        &
     &                    ib,nbb,sbuf)

          ib=ib+1

          call s_putbufgx(idwbc,idebc,'bnd',3,ni-3,ni,nj,nk,ptbr,       &
     &                    ib,nbb,sbuf)

        end if

        call s_shiftgx(idwbc,idebc,'bnd',nj,nk,nbb,sbuf,rbuf)

        ib=0

        if(gpvvar(2:2).eq.'o') then

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni+1,ni,nj,nk,ubr,        &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,vbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,pbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,ptbr,     &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,qvbr,     &
     &                    ib,nbb,rbuf)

        else

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni+1,ni,nj,nk,ubr,        &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,vbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,pbr,      &
     &                    ib,nbb,rbuf)

          ib=ib+1

          call s_getbufgx(idwbc,idebc,'bnd',0,ni_sub,ni,nj,nk,ptbr,     &
     &                    ib,nbb,rbuf)

        end if

! -----

! Set the periodic boundary conditions.

        if(gpvvar(2:2).eq.'o') then

          call bcyclex(idwbc,idebc,4,0,ni-3,ni+1,ni,nj,nk,ubr)
          call bcyclex(idwbc,idebc,3,0,ni-3,ni_sub,ni,nj,nk,vbr)
          call bcyclex(idwbc,idebc,3,0,ni-3,ni_sub,ni,nj,nk,pbr)
          call bcyclex(idwbc,idebc,3,0,ni-3,ni_sub,ni,nj,nk,ptbr)
          call bcyclex(idwbc,idebc,3,0,ni-3,ni_sub,ni,nj,nk,qvbr)

        else

          call bcyclex(idwbc,idebc,4,0,ni-3,ni+1,ni,nj,nk,ubr)
          call bcyclex(idwbc,idebc,3,0,ni-3,ni_sub,ni,nj,nk,vbr)
          call bcyclex(idwbc,idebc,3,0,ni-3,ni_sub,ni,nj,nk,pbr)
          call bcyclex(idwbc,idebc,3,0,ni-3,ni_sub,ni,nj,nk,ptbr)

        end if

! -----

!! -----

      end if

!!! -----

      end subroutine s_bcini

!-----7--------------------------------------------------------------7--

      end module m_bcini
