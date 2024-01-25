!***********************************************************************
      module m_chknlsfc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/01/20
!     Modification: 2008/01/11, 2008/03/12, 2008/04/17, 2008/05/02,
!                   2008/08/25, 2008/10/10, 2009/01/05, 2009/02/27,
!                   2011/08/09, 2011/09/22, 2011/11/10, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     check the namelist variables for surface.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comkind
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

      public :: chknlsfc, s_chknlsfc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface chknlsfc

        module procedure s_chknlsfc

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic aint
      intrinsic int
      intrinsic mod
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_chknlsfc(pname,ncpn,stat)
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

      integer iid      ! Index of do loops

      integer ierr     ! Error descriptor

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

      ierr=0

      fmsg=0

! -----

! For the section, runame.

      call numchar(exprim,1,ncexp,ncspc)

      if(ncexp.gt.64.or.ncspc.ne.0) then

        call destroy('chknlsfc',8,'cont',201,'exprim        ',6,101,    &
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

        call destroy('chknlsfc',8,'cont',201,'crsdir        ',6,101,    &
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

        call destroy('chknlsfc',8,'cont',201,'datdir        ',6,101,    &
     &               stat)

      end if

      if(pname(1:ncpn).ne.'check') then

        call outstd14('chknlsfc',8,crsdir,datdir,nccrs,ncdat)

      end if

! -----

! For the section, dimset.

      if(xdim.le.4) then

        call destroy('chknlsfc',8,'cont',201,'xdim          ',4,101,    &
     &               stat)

      end if

      if(mpopt.eq.5) then

        if(ydim.ne.4) then

          call destroy('chknlsfc',8,'cont',201,'ydim          ',4,101,  &
     &                 stat)

        end if

      else

        if(ydim.le.4) then

          call destroy('chknlsfc',8,'cont',201,'ydim          ',4,101,  &
     &                 stat)

        end if

      end if

      if(zdim.lt.7) then

        call destroy('chknlsfc',8,'cont',201,'zdim          ',4,101,    &
     &               stat)

      end if

! -----

! For the section, project.

      if(mpopt.ne.0.and.mpopt.ne.1.and.mpopt.ne.2.and.mpopt.ne.3.and.   &
     &   mpopt.ne.4.and.mpopt.ne.5.and.mpopt.ne.10.and.mpopt.ne.13) then

        call destroy('chknlsfc',8,'cont',201,'mpopt         ',5,101,    &
     &               stat)

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.3.or.mpopt.eq.13) then

        if(nspol.ne.1.and.nspol.ne.-1) then

          call destroy('chknlsfc',8,'cont',201,'nspol         ',5,101,  &
     &                 stat)

        end if

        if(tlat1.lt.-90.e0.or.tlat1.gt.90.e0) then

          call destroy('chknlsfc',8,'cont',201,'tlat1         ',5,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.2) then

        if(tlat2.lt.-90.e0.or.tlat2.gt.90.e0) then

          call destroy('chknlsfc',8,'cont',201,'tlat2         ',5,101,  &
     &                 stat)

        end if

      end if

      if(mpopt.eq.1.or.mpopt.eq.2.or.mpopt.eq.4) then

        if(tlon.lt.-180.e0.or.tlon.gt.180.e0) then

          call destroy('chknlsfc',8,'cont',201,'tlon          ',4,101,  &
     &                 stat)

        end if

      end if

      if(ulat.lt.-90.e0.or.ulat.gt.90.e0) then

        call destroy('chknlsfc',8,'cont',201,'ulat          ',4,101,    &
     &               stat)

      end if

      if(ulon.lt.-180.e0.or.ulon.gt.180.e0) then

        call destroy('chknlsfc',8,'cont',201,'ulon          ',4,101,    &
     &               stat)

      end if

      if(mpopt.eq.5) then

        if(disr.lt.eps) then

          call destroy('chknlsfc',8,'cont',201,'disr          ',4,101,  &
     &                 stat)

        end if

      end if

! -----

! For the section, gridset.

      if(mpopt.lt.10) then

        if(dx.lt.eps) then

          call destroy('chknlsfc',8,'cont',201,'dx            ',2,101,  &
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

          call destroy('chknlsfc',8,'cont',201,'dx            ',2,101,  &
     &                 stat)

        end if

      end if

      if(dy.lt.eps) then

        call destroy('chknlsfc',8,'cont',201,'dy            ',2,101,    &
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

        call destroy('chknlsfc',8,'cont',201,'dz            ',2,101,    &
     &               stat)

      else

        if(pname(1:ncpn).ne.'check') then
          dziv=1.e0/dz
        end if

      end if

! -----

! For the section, flength.

      if(pname(1:ncpn).ne.'check') then
        idate(1:4)=sfcast(1:4)
        idate(5:6)=sfcast(6:7)
        idate(7:8)=sfcast(9:10)
        idate(9:10)=sfcast(12:13)
        idate(11:12)=sfcast(15:16)
      end if

      if(sfcopt.eq.3.or.sfcopt.eq.13) then

       if(stime.lt.0.e0) then

         call destroy('chknlsfc',8,'cont',201,'stime         ',5,101,   &
     &                stat)

       end if

       if(etime.lt.0.e0) then

         call destroy('chknlsfc',8,'cont',201,'etime         ',5,101,   &
     &                stat)

       end if

       if(etime.lt.stime) then

         call destroy('chknlsfc',8,'cont',201,'etime         ',5,101,   &
     &                stat)

       end if

       if(1000_i8*int(stime+.1e0,i8).gt.1000_i8*int(etime+.1e0,i8)) then

         call destroy('chknlsfc',8,'cont',201,'stime         ',5,101,   &
     &                stat)

       end if

      end if

! -----

! For the section, sfcphys.

      ierr=0

      if(sfcdat(1:1).ne.'o'.and.sfcdat(1:1).ne.'x') then
        ierr=ierr+1
      end if

      if(sfcdat(2:2).ne.'o'.and.sfcdat(2:2).ne.'x') then
        ierr=ierr+1
      end if

      if(sfcdat(3:3).ne.'o'.and.sfcdat(3:3).ne.'x') then
        ierr=ierr+1
      end if

      if(ierr.ne.0) then

        call destroy('chknlsfc',8,'cont',201,'sfcdat        ',6,101,    &
     &               stat)

      end if

      if(sfcopt.ne.1.and.sfcopt.ne.2.and.sfcopt.ne.3                    &
     &  .and.sfcopt.ne.11.and.sfcopt.ne.12.and.sfcopt.ne.13) then

        call destroy('chknlsfc',8,'cont',201,'sfcopt        ',6,101,    &
     &               stat)

      end if

      if(sfcdat(1:3).eq.'xxx') then

        call destroy('chknlsfc',8,'cont',201,'sfcdat        ',6,101,    &
     &               stat)

      end if

      if(sfcdat(2:2).eq.'x') then

        if(sfcopt.eq.3.or.sfcopt.eq.13) then

          call destroy('chknlsfc',8,'cont',201,'sfcopt        ',6,101,  &
     &                 stat)

        end if

      else if(sfcdat(2:2).eq.'o') then

        if(sfcopt.eq.3.or.sfcopt.eq.13) then

          if(int(sstitv+.1e0).lt.60) then

            call destroy('chknlsfc',8,'cont',201,'sstitv        ',6,101,&
     &                   stat)

          end if

          if(pname(1:ncpn).ne.'check') then
            sstitv=aint(sstitv+.1e0)
          end if

        end if

      end if

! -----

! For the sections, project_lnd, gridset_lnd and datconf_lnd.

      if(sfcdat(1:1).eq.'o') then

        if(mpopt_lnd.ne.0.and.                                          &
     &     mpopt_lnd.ne.1.and.mpopt_lnd.ne.2.and.                       &
     &     mpopt_lnd.ne.3.and.mpopt_lnd.ne.4.and.                       &
     &     mpopt_lnd.ne.10.and.mpopt_lnd.ne.13) then

          call destroy('chknlsfc',8,'cont',201,'mpopt_lnd     ',9,101,  &
     &                 stat)

        end if

        if(mpopt.ge.10) then

          if(mpopt_lnd.lt.10) then

            call destroy('chknlsfc',8,'cont',201,'mpopt_lnd     ',9,101,&
     &                   stat)

          end if

        end if

        if(mpopt_lnd.eq.1.or.mpopt_lnd.eq.2                             &
     &    .or.mpopt_lnd.eq.3.or.mpopt_lnd.eq.13) then

          if(nspol_lnd.ne.1.and.nspol_lnd.ne.-1) then

            call destroy('chknlsfc',8,'cont',201,'nspol_lnd     ',9,101,&
     &                   stat)

          end if

          if(tlat1_lnd.lt.-90.e0.or.tlat1_lnd.gt.90.e0) then

            call destroy('chknlsfc',8,'cont',201,'tlat1_lnd     ',9,101,&
     &                   stat)

          end if

        end if

        if(mpopt_lnd.eq.2) then

          if(tlat2_lnd.lt.-90.e0.or.tlat2_lnd.gt.90.e0) then

            call destroy('chknlsfc',8,'cont',201,'tlat2_lnd     ',9,101,&
     &                   stat)

          end if

        end if

        if(mpopt_lnd.eq.1.or.mpopt_lnd.eq.2.or.mpopt_lnd.eq.4) then

          if(tlon_lnd.lt.-180.e0.or.tlon_lnd.gt.180.e0) then

            call destroy('chknlsfc',8,'cont',201,'tlon_lnd      ',8,101,&
     &                   stat)

          end if

        end if

        if(ulat_lnd.lt.-90.e0.or.ulat_lnd.gt.90.e0) then

          call destroy('chknlsfc',8,'cont',201,'ulat_lnd      ',8,101,  &
     &                 stat)

        end if

        if(ulon_lnd.lt.-180.e0.or.ulon_lnd.gt.180.e0) then

          call destroy('chknlsfc',8,'cont',201,'ulon_lnd      ',8,101,  &
     &                 stat)

        end if

        if(xdim_lnd.lt.2) then

          call destroy('chknlsfc',8,'cont',201,'xdim_lnd      ',8,101,  &
     &                 stat)

        end if

        if(ydim_lnd.lt.2) then

          call destroy('chknlsfc',8,'cont',201,'ydim_lnd      ',8,101,  &
     &                 stat)

        end if

        if(mpopt_lnd.lt.10) then

          if(dx_lnd.lt.eps) then

            call destroy('chknlsfc',8,'cont',201,'dx_lnd        ',6,101,&
     &                   stat)

          else

            if(pname(1:ncpn).ne.'check') then
              dxiv_lnd=1.e0/dx_lnd
            end if

          end if

        else

          if(pname(1:ncpn).ne.'check') then

            if(xdim_lnd.ge.4) then

              if(mpopt_lnd.eq.10) then

                fmsg=fmsg+1

                dx_lnd=360.e0/real(xdim_lnd)

                call outstd15('dx_lnd',6,fmsg,2,dx_lnd)

              end if

              if(mpopt_lnd.eq.13) then

                fmsg=fmsg+1

                dx_lnd=2.e0*cc*rearth/real(xdim_lnd)

                call outstd15('dx_lnd',6,fmsg,1,dx_lnd)

              end if

              dxiv_lnd=1.e0/dx_lnd

            end if

          end if

        end if

        if(dy_lnd.lt.eps) then

          call destroy('chknlsfc',8,'cont',201,'dy_lnd        ',6,101,  &
     &                 stat)

        else

          if(pname(1:ncpn).ne.'check') then
            dyiv_lnd=1.e0/dy_lnd
          end if

        end if

        if(intopt_lnd.ne.0.and.intopt_lnd.ne.1) then

          call destroy('chknlsfc',8,'cont',201,'intopt_lnd    ',10,101, &
     &                 stat)

        end if

        if(numctg_lnd.lt.1.or.numctg_lnd.gt.100) then

          call destroy('chknlsfc',8,'cont',201,'numctg_lnd    ',10,101, &
     &                 stat)

        end if

        if(numctg_lnd.ge.1.and.numctg_lnd.le.100) then

          do iid=1,numctg_lnd

            if(albe_lnd(iid).lt.0.e0.or.albe_lnd(iid).gt.1.e0) then

              call destroy('chknlsfc',8,'cont',201,'albe_lnd      ',8,  &
     &                     101,stat)

            end if

            if(beta_lnd(iid).lt.0.e0.or.beta_lnd(iid).gt.1.e0) then

              call destroy('chknlsfc',8,'cont',201,'beta_lnd      ',8,  &
     &                     101,stat)

            end if

            if(z0m_lnd(iid).lt.0.e0) then

              call destroy('chknlsfc',8,'cont',201,'z0m_lnd       ',7,  &
     &                     101,stat)

            end if

            if(z0h_lnd(iid).lt.0.e0) then

              call destroy('chknlsfc',8,'cont',201,'z0h_lnd       ',7,  &
     &                     101,stat)

            end if

            if(cap_lnd(iid).lt.0.e0) then

              call destroy('chknlsfc',8,'cont',201,'cap_lnd       ',7,  &
     &                     101,stat)

            end if

            if(nuu_lnd(iid).lt.0.e0) then

              call destroy('chknlsfc',8,'cont',201,'nuu_lnd       ',7,  &
     &                     101,stat)

            end if

          end do

        end if

      end if

! -----

! For the sections, project_sst and gridset_sst.

      if(sfcdat(2:2).eq.'o') then

        if(mpopt_sst.ne.0.and.                                          &
     &     mpopt_sst.ne.1.and.mpopt_sst.ne.2.and.                       &
     &     mpopt_sst.ne.3.and.mpopt_sst.ne.4.and.                       &
     &     mpopt_sst.ne.10.and.mpopt_sst.ne.13) then

          call destroy('chknlsfc',8,'cont',201,'mpopt_sst     ',9,101,  &
     &                 stat)

        end if

        if(mpopt.ge.10) then

          if(mpopt_sst.lt.10) then

            call destroy('chknlsfc',8,'cont',201,'mpopt_sst     ',9,101,&
     &                   stat)

          end if

        end if

        if(mpopt_sst.eq.1.or.mpopt_sst.eq.2                             &
     &    .or.mpopt_sst.eq.3.or.mpopt_sst.eq.13) then

          if(nspol_sst.ne.1.and.nspol_sst.ne.-1) then

            call destroy('chknlsfc',8,'cont',201,'nspol_sst     ',9,101,&
     &                   stat)

          end if

          if(tlat1_sst.lt.-90.e0.or.tlat1_sst.gt.90.e0) then

            call destroy('chknlsfc',8,'cont',201,'tlat1_sst     ',9,101,&
     &                   stat)

          end if

        end if

        if(mpopt_sst.eq.2) then

          if(tlat2_sst.lt.-90.e0.or.tlat2_sst.gt.90.e0) then

            call destroy('chknlsfc',8,'cont',201,'tlat2_sst     ',9,101,&
     &                   stat)

          end if

        end if

        if(mpopt_sst.eq.1.or.mpopt_sst.eq.2.or.mpopt_sst.eq.4) then

          if(tlon_sst.lt.-180.e0.or.tlon_sst.gt.180.e0) then

            call destroy('chknlsfc',8,'cont',201,'tlon_sst      ',8,101,&
     &                   stat)

          end if

        end if

        if(ulat_sst.lt.-90.e0.or.ulat_sst.gt.90.e0) then

          call destroy('chknlsfc',8,'cont',201,'ulat_sst      ',8,101,  &
     &                 stat)

        end if

        if(ulon_sst.lt.-180.e0.or.ulon_sst.gt.180.e0) then

          call destroy('chknlsfc',8,'cont',201,'ulon_sst      ',8,101,  &
     &                 stat)

        end if

        if(xdim_sst.lt.2) then

          call destroy('chknlsfc',8,'cont',201,'xdim_sst      ',8,101,  &
     &                 stat)

        end if

        if(ydim_sst.lt.2) then

          call destroy('chknlsfc',8,'cont',201,'ydim_sst      ',8,101,  &
     &                 stat)

        end if

        if(mpopt_sst.lt.10) then

          if(dx_sst.lt.eps) then

            call destroy('chknlsfc',8,'cont',201,'dx_sst        ',6,101,&
     &                   stat)

          else

            if(pname(1:ncpn).ne.'check') then
              dxiv_sst=1.e0/dx_sst
            end if

          end if

        else

          if(pname(1:ncpn).ne.'check') then

            if(xdim_sst.ge.4) then

              if(mpopt_sst.eq.10) then

                fmsg=fmsg+1

                dx_sst=360.e0/real(xdim_sst)

                call outstd15('dx_sst',6,fmsg,2,dx_sst)

              end if

              if(mpopt_sst.eq.13) then

                fmsg=fmsg+1

                dx_sst=2.e0*cc*rearth/real(xdim_sst)

                call outstd15('dx_sst',6,fmsg,1,dx_sst)

              end if

              dxiv_sst=1.e0/dx_sst

            end if

          end if

        end if

        if(dy_sst.lt.eps) then

          call destroy('chknlsfc',8,'cont',201,'dy_sst        ',6,101,  &
     &                 stat)

        else

          if(pname(1:ncpn).ne.'check') then
            dyiv_sst=1.e0/dy_sst
          end if

        end if

      end if

! -----

! For the sections, project_ice and gridset_ice.

      if(sfcdat(3:3).eq.'o') then

        if(mpopt_ice.ne.0.and.                                          &
     &     mpopt_ice.ne.1.and.mpopt_ice.ne.2.and.                       &
     &     mpopt_ice.ne.3.and.mpopt_ice.ne.4.and.                       &
     &     mpopt_ice.ne.10.and.mpopt_ice.ne.13) then

          call destroy('chknlsfc',8,'cont',201,'mpopt_ice     ',9,101,  &
     &                 stat)

        end if

        if(mpopt.ge.10) then

          if(mpopt_ice.lt.10) then

            call destroy('chknlsfc',8,'cont',201,'mpopt_ice     ',9,101,&
     &                   stat)

          end if

        end if

        if(mpopt_ice.eq.1.or.mpopt_ice.eq.2                             &
     &    .or.mpopt_ice.eq.3.or.mpopt_ice.eq.13) then

          if(nspol_ice.ne.1.and.nspol_ice.ne.-1) then

            call destroy('chknlsfc',8,'cont',201,'nspol_ice     ',9,101,&
     &                   stat)

          end if

          if(tlat1_ice.lt.-90.e0.or.tlat1_ice.gt.90.e0) then

            call destroy('chknlsfc',8,'cont',201,'tlat1_ice     ',9,101,&
     &                   stat)

          end if

        end if

        if(mpopt_ice.eq.2) then

          if(tlat2_ice.lt.-90.e0.or.tlat2_ice.gt.90.e0) then

            call destroy('chknlsfc',8,'cont',201,'tlat2_ice     ',9,101,&
     &                   stat)

          end if

        end if

        if(mpopt_ice.eq.1.or.mpopt_ice.eq.2.or.mpopt_ice.eq.4) then

          if(tlon_ice.lt.-180.e0.or.tlon_ice.gt.180.e0) then

            call destroy('chknlsfc',8,'cont',201,'tlon_ice      ',8,101,&
     &                   stat)

          end if

        end if

        if(ulat_ice.lt.-90.e0.or.ulat_ice.gt.90.e0) then

          call destroy('chknlsfc',8,'cont',201,'ulat_ice      ',8,101,  &
     &                 stat)

        end if

        if(ulon_ice.lt.-180.e0.or.ulon_ice.gt.180.e0) then

          call destroy('chknlsfc',8,'cont',201,'ulon_ice      ',8,101,  &
     &                 stat)

        end if

        if(xdim_ice.lt.2) then

          call destroy('chknlsfc',8,'cont',201,'xdim_ice      ',8,101,  &
     &                 stat)

        end if

        if(ydim_ice.lt.2) then

          call destroy('chknlsfc',8,'cont',201,'ydim_ice      ',8,101,  &
     &                 stat)

        end if

        if(mpopt_ice.lt.10) then

          if(dx_ice.lt.eps) then

            call destroy('chknlsfc',8,'cont',201,'dx_ice        ',6,101,&
     &                   stat)

          else

            if(pname(1:ncpn).ne.'check') then
              dxiv_ice=1.e0/dx_ice
            end if

          end if

        else

          if(pname(1:ncpn).ne.'check') then

            if(xdim_ice.ge.4) then

              if(mpopt_ice.eq.10) then

                fmsg=fmsg+1

                dx_ice=360.e0/real(xdim_ice)

                call outstd15('dx_ice',6,fmsg,2,dx_ice)

              end if

              if(mpopt_ice.eq.13) then

                fmsg=fmsg+1

                dx_ice=2.e0*cc*rearth/real(xdim_ice)

                call outstd15('dx_ice',6,fmsg,1,dx_ice)

              end if

              dxiv_ice=1.e0/dx_ice

            end if

          end if

        end if

        if(dy_ice.lt.eps) then

          call destroy('chknlsfc',8,'cont',201,'dy_ice        ',6,101,  &
     &                 stat)

        else

          if(pname(1:ncpn).ne.'check') then
            dyiv_ice=1.e0/dy_ice
          end if

        end if

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

        call destroy('chknlsfc',8,'cont',11,'              ',14,101,    &
     &               stat)

      else

        if(mod((xdim-3),xgroup*xsub).ne.0                               &
     &    .or.mod((ydim-3),ygroup*ysub).ne.0) then

          stat=stat+1

          call destroy('chknlsfc',8,'cont',11,'              ',14,101,  &
     &                 stat)

        end if

      end if

      if(numpe.lt.xsub*ysub.or.numpe.gt.xsub*ysub*xgroup*ygroup) then

        stat=stat+1

        call destroy('chknlsfc',8,'cont',11,'              ',14,101,    &
     &               stat)

      end if

! -----

      end subroutine s_chknlsfc

!-----7--------------------------------------------------------------7--

      end module m_chknlsfc
