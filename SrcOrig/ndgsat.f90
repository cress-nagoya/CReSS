!***********************************************************************
      module m_ndgsat
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/05/07
!     Modification: 2007/07/30, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/03/23, 2009/03/31, 2010/05/17, 2011/09/22,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     Control the inferior procedures for analysis nudging to radar data
!     of potential temperature and water wapor mixing ratio.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_getcname
      use m_getetop
      use m_getiname
      use m_getqvs
      use m_gettlcl
      use m_inichar
      use m_qv2rdr

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: ndgsat, s_ndgsat

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface ndgsat

        module procedure s_ndgsat

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
      subroutine s_ndgsat(fpngrvar,fpngropt,ngrdmp,rtinc,               &
     &                    ni,nj,nk,nqw,nqi,zph,pbr,ptbr,rst,            &
     &                    ppp,ptpp,qvp,qwrdr,qwrtd,qirdr,qirtd,         &
     &                    qvfrc,qvs,qprdr,zph8s,etop,tlcl,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpngrvar
                       ! Formal parameter of unique index of ngrvar

      integer, intent(in) :: fpngropt
                       ! Formal parameter of unique index of ngropt

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

      real, intent(in) :: ngrdmp(1:2)
                       ! Analysis nudging damping coefficient for radar

      real, intent(in) :: rtinc(1:2)
                       ! Lapse of forecast time from radar data reading

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qwrdr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of radar data
                       ! at marked time

      real, intent(in) :: qwrtd(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Time tendency of
                       ! water hydrometeor of radar data

      real, intent(in) :: qirdr(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of radar data
                       ! at marked time

      real, intent(in) :: qirtd(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Time tendency of
                       ! ice hydrometeor of radar data

! Input and output variable

      real, intent(inout) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in water vapor mixing ratio equation

! Internal shared variables

      character(len=108) ngrvar
                       ! Control flag of
                       ! analysis nudged variables to radar data

      integer ngropt   ! Option for analysis nudging to radar

      real, intent(inout) :: qvs(0:ni+1,0:nj+1,1:nk)
                       ! Saturation mixing ratio

      real, intent(inout) :: qprdr(0:ni+1,0:nj+1,1:nk)
                       ! Total precipitation mixing ratio of radar data
                       ! at marked time

      real, intent(inout) :: zph8s(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at scalar points

      real, intent(inout) :: etop(0:ni+1,0:nj+1)
                       ! z physical coordinates at radar echo top

      real, intent(inout) :: tlcl(0:ni+1,0:nj+1)
                       ! Air temperature at lifting condensation level

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1)
                       ! Temporary array

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(ngrvar)

! -----

! Get the required namelist variables.

      call getcname(fpngrvar,ngrvar)
      call getiname(fpngropt,ngropt)

! -----

!! Perform analysis nudging to radar data of potential temperature and
!! water wapor mixing ratio.

      if(ngropt.ge.1.and.(ngrdmp(1).gt.0.e0.or.ngrdmp(2).gt.0.e0)) then

        if(ngrvar(4:4).eq.'o') then

! Calculate the saturation mixing ratio.

          call getqvs(ni,nj,nk,pbr,ptbr,ppp,ptpp,qvs)

! -----

! Get the z physical coordinates at scalar points and radar echo top and
! total precipitation mixing ratio of radar data

          call getetop(idngropt,idhaiopt,ngrdmp,rtinc,ni,nj,nk,nqw,nqi, &
     &                 zph,qwrdr,qwrtd,qirdr,qirtd,zph8s,etop,qprdr)

! -----

! Get the air temperature at Lifting Condensation Level.

          call gettlcl('m',2,ni,nj,nk,pbr,ptbr,ppp,ptpp,qvp,tlcl)

! -----

! Perform the analysis nudging to radar data of water vapor mixing
! ratio.

          call qv2rdr(idngropt,ngrdmp,ni,nj,nk,zph,zph8s,etop,tlcl,     &
     &                pbr,ptbr,rst,ppp,ptpp,qvp,qvs,qvfrc,tmp1)

! -----

        end if

      end if

!! -----

      end subroutine s_ndgsat

!-----7--------------------------------------------------------------7--

      end module m_ndgsat
