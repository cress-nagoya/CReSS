!***********************************************************************
      module m_setdim
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/01/20
!     Modification: 2007/12/25, 2008/04/17, 2008/05/02, 2008/08/25,
!                   2008/12/11, 2009/01/30, 2009/02/27, 2009/11/13,
!                   2011/08/18, 2011/09/22, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the model dimension variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_commpi
      use m_defdim
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setdim, s_setdim

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setdim

        module procedure s_setdim

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_setdim(fpcphopt,fphaiopt,fpaslopt)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fpaslopt
                       ! Formal parameter of unique index of aslopt

! Internal shared variables

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes
      integer aslopt   ! Option for aerosol processes

      integer numpe    ! Total number of processor elements

      integer xgroup   ! Number of group domain in entire domain
                       ! in x direction

      integer ygroup   ! Number of group domain in entire domain
                       ! in y direction

      integer xsub     ! Number of sub domain in group domain
                       ! in x direction

      integer ysub     ! Number of sub domain in group domain
                       ! in y direction

      integer xdim     ! Model dimension in x direction
      integer ydim     ! Model dimension in y direction
      integer zdim     ! Model dimension in z direction

      integer levund   ! Number of soil and sea layers

      integer snddim   ! Sounding data dimension

      integer ncbinw   ! Number of categories for warm bin

      integer xdim_gpv ! GPV data dimension in x direction
      integer ydim_gpv ! GPV data dimension in y direction
      integer zdim_gpv ! GPV data dimension in z direction

      integer xdim_asl ! Aerosol data dimension in x direction
      integer ydim_asl ! Aerosol data dimension in y direction
      integer zdim_asl ! Aerosol data dimension in z direction

      integer xdim_rdr ! Radar data dimension in x direction
      integer ydim_rdr ! Radar data dimension in y direction
      integer zdim_rdr ! Radar data dimension in z direction

      integer xdim_trn ! Terrain data dimension in x direction
      integer ydim_trn ! Terrain data dimension in y direction

      integer xdim_lnd ! Land use data dimension in x direction
      integer ydim_lnd ! Land use data dimension in y direction

      integer xdim_sst ! Sea surface temperature data dimension
                       ! in x direction

      integer ydim_sst ! Sea surface temperature data dimension
                       ! in y direction

      integer xdim_ice ! Sea ice distribution data dimension
                       ! in x direction

      integer ydim_ice ! Sea ice distribution data dimension
                       ! in y direction

      integer xsub_rst ! Number of sub domain in group domain
                       ! in x direction for restructed files

      integer ysub_rst ! Number of sub domain in group domain
                       ! in y direction for restructed files

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getiname(fpaslopt,aslopt)

      call getiname(idnumpe,numpe)
      call getiname(idxgroup,xgroup)
      call getiname(idygroup,ygroup)
      call getiname(idxsub,xsub)
      call getiname(idysub,ysub)
      call getiname(idxdim,xdim)
      call getiname(idydim,ydim)
      call getiname(idzdim,zdim)
      call getiname(idlevund,levund)
      call getiname(idsnddim,snddim)
      call getiname(idncbinw,ncbinw)
      call getiname(idxdim_gpv,xdim_gpv)
      call getiname(idydim_gpv,ydim_gpv)
      call getiname(idzdim_gpv,zdim_gpv)
      call getiname(idxdim_asl,xdim_asl)
      call getiname(idydim_asl,ydim_asl)
      call getiname(idzdim_asl,zdim_asl)
      call getiname(idxdim_rdr,xdim_rdr)
      call getiname(idydim_rdr,ydim_rdr)
      call getiname(idzdim_rdr,zdim_rdr)
      call getiname(idxdim_trn,xdim_trn)
      call getiname(idydim_trn,ydim_trn)
      call getiname(idxdim_lnd,xdim_lnd)
      call getiname(idydim_lnd,ydim_lnd)
      call getiname(idxdim_sst,xdim_sst)
      call getiname(idydim_sst,ydim_sst)
      call getiname(idxdim_ice,xdim_ice)
      call getiname(idydim_ice,ydim_ice)
      call getiname(idxsub_rst,xsub_rst)
      call getiname(idysub_rst,ysub_rst)

! -----

!! Set the model dimension variables.

! Set the dimension variables of CReSS grid points.

      ni=max((xdim-3)/(xsub*xgroup)+3,4)
      nj=max((ydim-3)/(ysub*ygroup)+3,4)

      nk=max(zdim,4)

! -----

! Set the dimension variables of categories of cloud physics.

      if(abs(cphopt).lt.10) then

        if(abs(cphopt).eq.0) then
          nqw=1
        else
          nqw=2
        end if

        if(abs(cphopt).le.1) then
          nnw=1
        else
          nnw=nqw
        end if

        if(haiopt.eq.0) then

          if(abs(cphopt).le.1) then
            nqi=1
            nni=1

          else
            nqi=3
            nni=nqi

          end if

        else

          if(abs(cphopt).le.1) then
            nqi=1
            nni=1

          else
            nqi=4
            nni=nqi

          end if

        end if

      else if(abs(cphopt).gt.10.and.abs(cphopt).lt.20) then

        nqw=ncbinw
        nnw=nqw

        if(abs(cphopt).le.11) then
          nqi=1
          nni=1

        else
          nqi=ncbinw
          nni=nqi

        end if

      end if

      nqw=max(nqw,1)
      nnw=max(nnw,1)

      nqi=max(nqi,1)
      nni=max(nni,1)

! -----

! Set the subordinate dimension variable.

      if(abs(cphopt).lt.10) then

        km=nk

      else

        if(nqw.gt.nqi) then
          km=nqw+1
        else
          km=nqi+1
        end if

        km=max(nk,km)

      end if

! -----

! Set the dimension variable of categories of aerosol.

      if(aslopt.eq.0) then
        nqa(1)=0
        nqa(2)=0
        nqa(3)=0
        nqa(4)=0

        nqa(0)=1

      else
        nqa(1)=6+nqw+nqi
        nqa(2)=4+nqw+nqi
        nqa(3)=3+nqw+nqi
        nqa(4)=4+nqw+nqi

        nqa(0)=nqa(1)+nqa(2)+nqa(3)+nqa(4)

      end if

! -----

! Set the dimension variable of soil temperature calculation.

      nund=max(levund,1)

! -----

! Set the dimension variable of sounding data.

      nlev=max(4*snddim,4)

! -----

! Set the dimension variables of GPV data.

      nid_gpv=max(xdim_gpv,1)
      njd_gpv=max(ydim_gpv,1)
      nkd_gpv=max(zdim_gpv,1)

      km_gpv=max(nk,nkd_gpv)

! -----

! Set the dimension variables of aerosol data.

      nid_asl=max(xdim_asl,1)
      njd_asl=max(ydim_asl,1)
      nkd_asl=max(zdim_asl,1)

      km_asl=max(nk,nkd_asl)

! -----

! Set the dimension variables of radar data.

      nid_rdr=max(xdim_rdr,1)
      njd_rdr=max(ydim_rdr,1)
      nkd_rdr=max(zdim_rdr,1)

      km_rdr=max(nk,nkd_rdr)

! -----

! Set the dimension variables of terrain data.

      nid_trn=max(xdim_trn,1)
      njd_trn=max(ydim_trn,1)

! -----

! Set the dimension variables of surface data.

      nid_lnd=max(xdim_lnd,1)
      njd_lnd=max(ydim_lnd,1)

      nid_sst=max(xdim_sst,1)
      njd_sst=max(ydim_sst,1)

      nid_ice=max(xdim_ice,1)
      njd_ice=max(ydim_ice,1)

! -----

!! -----

! Set the parameters of parallelizing again.

      npe=max(numpe,1)

      nsub=max(xsub*ysub,1)

      ngrp=max(xgroup*ygroup,1)

      nsrl=max(numpe/(xsub*ysub),1)

      nsub_rst=max(xsub_rst*ysub_rst,1)

      nisub=max(xsub,1)
      njsub=max(ysub,1)

      nigrp=max(xgroup,1)
      njgrp=max(ygroup,1)

      nisub_rst=max(xsub_rst,1)
      njsub_rst=max(ysub_rst,1)

! -----

      end subroutine s_setdim

!-----7--------------------------------------------------------------7--

      end module m_setdim
