!***********************************************************************
      module m_chknltrn
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/01/20
!     Modification: 2008/01/11, 2008/03/12, 2008/04/17, 2008/05/02,
!                   2008/08/25, 2008/10/10, 2009/01/05, 2009/02/27,
!                   2011/08/09, 2011/09/22, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the namelist variables for terrain.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy
      use m_defname
      use m_destroy
      use m_numchar
      use m_outstd14
      use m_outstd15

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: chknltrn, s_chknltrn

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chknltrn

        module procedure s_chknltrn

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic mod
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_chknltrn(pname,ncpn,stat)
!***********************************************************************

! Input variables

      character(len=8), intent(in) :: pname
                       ! Running program name

      integer, intent(in) :: ncpn
                       ! Number of character of pname

! Output variable

      integer, intent(out) :: stat
                       ! Runtime status

! Internal shared variables

      integer fmsg     ! Control flag of message type for standard i/o

      integer ncspc    ! Number of space character of exprim and prvres

      real r4scl       ! Optional real scalar variable

!-----7--------------------------------------------------------------7--

! Set the word length for direct access file.

      r4scl=0.e0

      inquire(iolength=wlngth) r4scl

! -----

!! Check the namelist variables.

! Initialize the controler.

      stat=0

      fmsg=0

! -----

! For the section, runame.

      call numchar(exprim,1,ncexp,ncspc)

      if(ncexp.gt.64.or.ncspc.ne.0) then

        call destroy('chknltrn',8,'cont',201,'exprim        ',6,101,    &
     &               stat)

      end if

      if(pname(1:ncpn).ne.'check') then

        write(exprim(ncexp+1:ncexp+1),'(a1)') '.'

        ncexp=ncexp+1

      end if

! -----

! For the section, drname.

      call numchar(crsdir,1,nccrs,ncspc)

      if(crsdir(nccrs:nccrs).ne.'/') then

        nccrs=nccrs+1

        if(nccrs.le.108) then

          write(crsdir(nccrs:nccrs),'(a1)') '/'

        end if

      end if

      if(nccrs.gt.108.or.ncspc.ne.0) then

        call destroy('chknltrn',8,'cont',201,'crsdir        ',6,101,    &
     &               stat)

      end if

      call numchar(datdir,1,ncdat,ncspc)

      if(datdir(ncdat:ncdat).ne.'/') then

        ncdat=ncdat+1

        if(ncdat.le.108) then

          write(datdir(ncdat:ncdat),'(a1)') '/'

        end if

      end if

      if(ncdat.gt.108.or.ncspc.ne.0) then

        call destroy('chknltrn',8,'cont',201,'datdir        ',6,101,    &
     &               stat)

      end if

      if(pname(1:ncpn).ne.'check') then

        call outstd14('chknltrn',8,crsdir,datdir,nccrs,ncdat)

      end if

! -----

! For the section, dimset.

      if(xdim.le.4) then

        call destroy('chknltrn',8,'cont',201,'xdim          ',4,101,    &
     &               stat)

      end if

      if(mpopt.eq.5) then

        if(ydim.ne.4) then

          call destroy('chknltrn',8,'cont',201,'ydim          ',4,101,  &
     &                 stat)

        end if

      else

        if(ydim.le.4) then

          call destroy('chknltrn',8,'cont',201,'ydim          ',4,101,  &
     &                 stat)

        end if

      end if

      if(zdim.lt.7) then

        call destroy('chknltrn',8,'cont',201,'zdim          ',4,101,    &
     &               stat)

      end if

! -----

! For the section, project.

      if(mpopt.ne.0.and.mpopt.ne.1.and.mpopt.ne.2.and.mpopt.ne.3.and.   &
     &   mpopt.ne.4.and.mpopt.ne.5.and.mpopt.ne.10.and.mpopt.ne.13) then

        call destroy('chknltrn',8,'cont',201,'mpopt         ',5,101,    &
     &               stat)

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.3.or.mpopt.eq.13) then

        if(nspol.ne.1.and.nspol.ne.-1) then

          call destroy('chknltrn',8,'cont',201,'nspol         ',5,101,  &
     &                 stat)

        end if

        if(tlat1.lt.-90.e0.or.tlat1.gt.90.e0) then

          call destroy('chknltrn',8,'cont',201,'tlat1         ',5,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.2) then

        if(tlat2.lt.-90.e0.or.tlat2.gt.90.e0) then

          call destroy('chknltrn',8,'cont',201,'tlat2         ',5,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.4) then

        if(tlon.lt.-180.e0.or.tlon.gt.180.e0) then

          call destroy('chknltrn',8,'cont',201,'tlon          ',4,101,  &
     &                 stat)

        end if

      end if

      if(ulat.lt.-90.e0.or.ulat.gt.90.e0) then

        call destroy('chknltrn',8,'cont',201,'ulat          ',4,101,    &
     &               stat)

      end if

      if(ulon.lt.-180.e0.or.ulon.gt.180.e0) then

        call destroy('chknltrn',8,'cont',201,'ulon          ',4,101,    &
     &               stat)

      end if

      if(mpopt.eq.5) then

        if(disr.lt.eps) then

          call destroy('chknltrn',8,'cont',201,'disr          ',4,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, gridset.

      if(mpopt.lt.10) then

        if(dx.lt.eps) then

          call destroy('chknltrn',8,'cont',201,'dx            ',2,101,  &
     &                 stat)

        else

          if(pname(1:ncpn).ne.'check') then

            if(mpopt.eq.0) then
              dx=rearth*d2r*dx
            end if

            dxiv=1.e0/dx

          end if

        end if

      else

        if(pname(1:ncpn).ne.'check') then

          if(xdim.ge.4) then

            if(mpopt.eq.10) then

              fmsg=fmsg+1

              dx=360.e0/real(xdim)

              call outstd15('dx    ',2,fmsg,2,dx)

              dx=rearth*d2r*dx

            end if

            if(mpopt.eq.13) then

              fmsg=fmsg+1

              dx=2.e0*cc*rearth/real(xdim-3)

              call outstd15('dx    ',2,fmsg,1,dx)

            end if

            dxiv=1.e0/dx

          end if

        end if

      end if

      if(mpopt.eq.5) then

        if(2.e0*dx.lt.disr) then

          call destroy('chknltrn',8,'cont',201,'dx            ',2,101,  &
     &                 stat)

        end if

      end if

      if(dy.lt.eps) then

        call destroy('chknltrn',8,'cont',201,'dy            ',2,101,    &
     &               stat)

      else

        if(pname(1:ncpn).ne.'check') then

          if(mpopt.eq.0.or.mpopt.eq.10) then
            dy=rearth*d2r*dy
          end if

          dyiv=1.e0/dy

        end if

      end if

      if(dz.lt.eps) then

        call destroy('chknltrn',8,'cont',201,'dz            ',2,101,    &
     &               stat)

      else

        if(pname(1:ncpn).ne.'check') then
          dziv=1.e0/dz
        end if

      end if

! -----

! For the section, terrain.

      if(trnopt.ne.2) then

        call destroy('chknltrn',8,'cont',201,'trnopt        ',6,101,    &
     &               stat)

      end if

! -----

! For the sections, project_trn, gridset_trn and datconf_trn.

      if(mpopt_trn.ne.0.and.                                            &
     &   mpopt_trn.ne.1.and.mpopt_trn.ne.2.and.                         &
     &   mpopt_trn.ne.3.and.mpopt_trn.ne.4.and.                         &
     &   mpopt_trn.ne.10.and.mpopt_trn.ne.13) then

        call destroy('chknltrn',8,'cont',201,'mpopt_trn     ',9,101,    &
     &               stat)

      end if

      if(mpopt.ge.10) then

        if(mpopt_trn.lt.10) then

          call destroy('chknltrn',8,'cont',201,'mpopt_trn     ',9,101,  &
     &                 stat)

        end if

      end if

      if(mpopt_trn.eq.1.or.mpopt_trn.eq.2                               &
     &  .or.mpopt_trn.eq.3.or.mpopt_trn.eq.13) then

        if(nspol_trn.ne.1.and.nspol_trn.ne.-1) then

          call destroy('chknltrn',8,'cont',201,'nspol_trn     ',9,101,  &
     &                 stat)

        end if

        if(tlat1_trn.lt.-90.e0.or.tlat1_trn.gt.90.e0) then

          call destroy('chknltrn',8,'cont',201,'tlat1_trn     ',9,101,  &
     &                 stat)

        end if

      end if

      if(mpopt_trn.eq.2) then

        if(tlat2_trn.lt.-90.e0.or.tlat2_trn.gt.90.e0) then

          call destroy('chknltrn',8,'cont',201,'tlat2_trn     ',9,101,  &
     &                 stat)

        end if

      end if

      if(mpopt_trn.eq.1.or.mpopt_trn.eq.2.or.mpopt_trn.eq.4) then

        if(tlon_trn.lt.-180.e0.or.tlon_trn.gt.180.e0) then

          call destroy('chknltrn',8,'cont',201,'tlon_trn      ',8,101,  &
     &                 stat)

        end if

      end if

      if(ulat_trn.lt.-90.e0.or.ulat_trn.gt.90.e0) then

        call destroy('chknltrn',8,'cont',201,'ulat_trn      ',8,101,    &
     &               stat)

      end if

      if(ulon_trn.lt.-180.e0.or.ulon_trn.gt.180.e0) then

        call destroy('chknltrn',8,'cont',201,'ulon_trn      ',8,101,    &
     &               stat)

      end if

      if(xdim_trn.lt.2) then

        call destroy('chknltrn',8,'cont',201,'xdim_trn      ',8,101,    &
     &               stat)

      end if

      if(ydim_trn.lt.2) then

        call destroy('chknltrn',8,'cont',201,'ydim_trn      ',8,101,    &
     &               stat)

      end if

      if(mpopt_trn.lt.10) then

        if(dx_trn.lt.eps) then

          call destroy('chknltrn',8,'cont',201,'dx_trn        ',6,101,  &
     &                 stat)

        else

          if(pname(1:ncpn).ne.'check') then
            dxiv_trn=1.e0/dx_trn
          end if

        end if

      else

        if(pname(1:ncpn).ne.'check') then

          if(xdim_trn.ge.4) then

            if(mpopt_trn.eq.10) then

              fmsg=fmsg+1

              dx_trn=360.e0/real(xdim_trn)

              call outstd15('dx_trn',6,fmsg,2,dx_trn)

            end if

            if(mpopt_trn.eq.13) then

              fmsg=fmsg+1

              dx_trn=2.e0*cc*rearth/real(xdim_trn)

              call outstd15('dx_trn',6,fmsg,1,dx_trn)

            end if

            dxiv_trn=1.e0/dx_trn

          end if

        end if

      end if

      if(dy_trn.lt.eps) then

        call destroy('chknltrn',8,'cont',201,'dy_trn        ',6,101,    &
     &               stat)

      else

        if(pname(1:ncpn).ne.'check') then
          dyiv_trn=1.e0/dy_trn
        end if

      end if

      if(intopt_trn.ne.1.and.intopt_trn.ne.2) then

        call destroy('chknltrn',8,'cont',201,'intopt_trn    ',10,101,   &
     &               stat)

      end if

! -----

!! -----

! Set the constant value.

      if(pname(1:ncpn).ne.'check') then

        write(gpvvar(9:9),'(a1)') 'o'
        write(exbvar(9:9),'(a1)') 'x'

      end if

! -----

! Check the parameters of parallelizing.

      if(xgroup.eq.0.or.xsub.eq.0.or.ygroup.eq.0.or.ysub.eq.0) then

        stat=stat+1

        call destroy('chknltrn',8,'cont',11,'              ',14,101,    &
     &               stat)

      else

        if(mod((xdim-3),xgroup*xsub).ne.0                               &
     &    .or.mod((ydim-3),ygroup*ysub).ne.0) then

          stat=stat+1

          call destroy('chknltrn',8,'cont',11,'              ',14,101,  &
     &                 stat)

        end if

      end if

      if(numpe.lt.xsub*ysub.or.numpe.gt.xsub*ysub*xgroup*ygroup) then

        stat=stat+1

        call destroy('chknltrn',8,'cont',11,'              ',14,101,    &
     &               stat)

      end if

! -----

      end subroutine s_chknltrn

!-----7--------------------------------------------------------------7--

      end module m_chknltrn
