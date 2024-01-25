!***********************************************************************
      module m_outpbl
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/12/10
!     Modification: 2002/04/02, 2002/07/03, 2002/12/02, 2003/01/20,
!                   2003/03/28, 2003/04/30, 2003/05/19, 2003/07/15,
!                   2003/09/01, 2003/10/31, 2003/12/12, 2004/02/01,
!                   2004/03/05, 2004/04/01, 2004/05/07, 2004/05/31,
!                   2004/06/10, 2004/08/20, 2004/09/01, 2004/09/10,
!                   2004/09/25, 2005/01/31, 2006/01/10, 2006/02/13,
!                   2006/05/12, 2006/09/21, 2006/11/06, 2007/01/05,
!                   2007/01/20, 2007/05/14, 2007/05/21, 2007/06/27,
!                   2007/07/30, 2007/09/04, 2007/10/19, 2008/01/11,
!                   2008/03/12, 2008/04/17, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2008/10/10, 2009/02/27, 2009/08/20,
!                   2011/08/18, 2011/09/22, 2013/01/28, 2013/02/13,
!                   2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     read in the surface variables to the dumped file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bulksfc
      use m_comcapt
      use m_comdmp
      use m_comindx
      use m_comphy
      use m_comtable
      use m_getcname
      use m_getiname
      use m_inichar
      use m_outdmp2d
      use m_rotuvm2s
      use m_setcst2d
      use m_setproj

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: outpbl, s_outpbl

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface outpbl

        module procedure s_outpbl

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic log

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_outpbl(fpdmpvar,fpdmplev,fmois,ni,nj,nk,nund,        &
     &                    za,lon,p,u,v,qv,ufrc,vfrc,ptfrc,qvfrc,        &
     &                    land,kai,z0m,z0h,ptv,qvsfc,rch,cm,ch,         &
     &                    tund,tice,hs,le,rgd,rsd,rld,rlu,cdl,cdm,cdh,  &
     &                    z10,z15,cm10,ch15,u10,v10,p15,pt15,qv15,      &
     &                    tsfc,cdave,usflx,vsflx,ptsflx,qvsflx)
!***********************************************************************

! Input variables

      character(len=5), intent(in) :: fmois
                       ! Control flag of air moisture

      integer, intent(in) :: fpdmpvar
                       ! Formal parameter of unique index of dmpvar

      integer, intent(in) :: fpdmplev
                       ! Formal parameter of unique index of dmplev

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nund
                       ! Number of soil and sea layers

      integer, intent(in) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(in) :: za(0:ni+1,0:nj+1)
                       ! z physical coordinates at lowest plane

      real, intent(in) :: lon(0:ni+1,0:nj+1)
                       ! Longitude

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity

      real, intent(in) :: qv(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio

      real, intent(in) :: ufrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in u equation

      real, intent(in) :: vfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in v equation

      real, intent(in) :: ptfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in potential temperature equation

      real, intent(in) :: qvfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term
                       ! in water vapor mixing ratio equation

      real, intent(in) :: kai(0:ni+1,0:nj+1)
                       ! Sea ice distribution

      real, intent(in) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(in) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(in) :: ptv(0:ni+1,0:nj+1,1:nk)
                       ! Virtual potential temperature

      real, intent(in) :: qvsfc(0:ni+1,0:nj+1)
                       ! Water vapor mixing ratio on surface

      real, intent(in) :: rch(0:ni+1,0:nj+1)
                       ! Bulk Richardson number

      real, intent(in) :: cm(0:ni+1,0:nj+1)
                       ! Bulk coefficient for velocity

      real, intent(in) :: ch(0:ni+1,0:nj+1)
                       ! Bulk coefficient for scalar

      real, intent(in) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature

      real, intent(in) :: tice(0:ni+1,0:nj+1)
                       ! Mixed ice surface temperature

      real, intent(in) :: hs(0:ni+1,0:nj+1)
                       ! Sensible heat

      real, intent(in) :: le(0:ni+1,0:nj+1)
                       ! Latent heat

      real, intent(in) :: rgd(0:ni+1,0:nj+1)
                       ! Global solar radiation

      real, intent(in) :: rsd(0:ni+1,0:nj+1)
                       ! Net downward short wave radiation

      real, intent(in) :: rld(0:ni+1,0:nj+1)
                       ! Downward long wave radiation

      real, intent(in) :: rlu(0:ni+1,0:nj+1)
                       ! Upward long wave radiation

      real, intent(in) :: cdl(0:ni+1,0:nj+1)
                       ! Cloud cover in lower layer

      real, intent(in) :: cdm(0:ni+1,0:nj+1)
                       ! Cloud cover in middle layer

      real, intent(in) :: cdh(0:ni+1,0:nj+1)
                       ! Cloud cover in upper layer

! Internal shared variables

      character(len=108) dmpvar
                       ! Control flag of dumped variables

      integer dmplev   ! Option for z coordinates of dumped variables

      real rddwkp      ! ln(da0 / dv0) / wkappa

      real cpj(1:7)    ! Map projection parameters

      real, intent(inout) :: z10(0:ni+1,0:nj+1)
                       ! Constant height of 10m

      real, intent(inout) :: z15(0:ni+1,0:nj+1)
                       ! Constant height of 1.5m

      real, intent(inout) :: cm10(0:ni+1,0:nj+1)
                       ! Bulk coefficient for velocity
                       ! at an altitude of 10m

      real, intent(inout) :: ch15(0:ni+1,0:nj+1)
                       ! Bulk coefficient for scalar
                       ! at an altitude of 1.5m

      real, intent(inout) :: u10(0:ni+1,0:nj+1)
                       ! x components of velocity at an altitude of 10m

      real, intent(inout) :: v10(0:ni+1,0:nj+1)
                       ! y components of velocity at an altitude of 10m

      real, intent(inout) :: p15(0:ni+1,0:nj+1)
                       ! Pressure at an altitude of 1.5m

      real, intent(inout) :: pt15(0:ni+1,0:nj+1)
                       ! Potential temperature at an altitude of 1.5m

      real, intent(inout) :: qv15(0:ni+1,0:nj+1)
                       ! Water vapor mixing ratio at an altitude of 1.5m

      real, intent(inout) :: tsfc(0:ni+1,0:nj+1)
                       ! Soil and sea surface temperature

      real, intent(inout) :: cdave(0:ni+1,0:nj+1)
                       ! Averaged cloud cover

      real, intent(inout) :: usflx(0:ni+1,0:nj+1)
                       ! Surface stress in x-z coordinates

      real, intent(inout) :: vsflx(0:ni+1,0:nj+1)
                       ! Surface stress in y-z coordinates

      real, intent(inout) :: ptsflx(0:ni+1,0:nj+1)
                       ! Surface heat flux

      real, intent(inout) :: qvsflx(0:ni+1,0:nj+1)
                       ! Surface moisture flux

! Internal private variables

      integer i        ! Array index in x drection
      integer j        ! Array index in y drection

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(dmpvar)

! -----

! Get the required namelist variable.

      call getcname(fpdmpvar,dmpvar)

! -----

!! Read in the surface predicted variables to the dumped file.

      if(fmon(1:3).eq.'act'.and.dmpvar(13:13).eq.'o') then

! Get the required namelist variable.

        call getiname(fpdmplev,dmplev)

! -----

! Set the common used variable.

        rddwkp=log(da0/dv0)/wkappa

! -----

! Set the constant height of 10.0 and 1.5 m.

        call setcst2d(0,ni+1,0,nj+1,10.e0,z10)
        call setcst2d(0,ni+1,0,nj+1,1.5e0,z15)

! -----

! Calculate the bulk coefficients.

        call bulksfc(ni,nj,z10,land,kai,z0m,z0h,rch,cm10,u10)
        call bulksfc(ni,nj,z15,land,kai,z0m,z0h,rch,p15,ch15)

! -----

! Calculate the x and y components of velocity at an altitude of 10 m,
! the pressure, potential temperature and water vapor mixing ratio at
! an altitude of 1.5 m, the surface temperature, total cloud cover and
! surface fluxes.

!$omp parallel default(shared)

        if(fmois(1:3).eq.'dry') then

!$omp do schedule(runtime) private(i,j,a)

          do j=1,nj-1
          do i=1,ni-1
            a=.5e0*cm(i,j)/cm10(i,j)

            u10(i,j)=a*(u(i,j,2)+u(i+1,j,2))
            v10(i,j)=a*(v(i,j,2)+v(i,j+1,2))

            a=ch(i,j)/ch15(i,j)

            pt15(i,j)=ptv(i,j,1)+a*(ptv(i,j,2)-ptv(i,j,1))
            qv15(i,j)=0.e0

            p15(i,j)=.5e0*(p(i,j,1)*(za(i,j)-1.5e0)                     &
     &        +p(i,j,2)*(za(i,j)+1.5e0))/za(i,j)

            if(land(i,j).eq.1) then
              tsfc(i,j)=kai(i,j)*tice(i,j)+(1.e0-kai(i,j))*tund(i,j,1)
            else
              tsfc(i,j)=tund(i,j,1)
            end if

            cdave(i,j)=0.e0

            usflx(i,j)=-.5e0*(ufrc(i,j,1)+ufrc(i+1,j,1))
            vsflx(i,j)=-.5e0*(vfrc(i,j,1)+vfrc(i,j+1,1))

            ptsflx(i,j)=-ptfrc(i,j,1)

            qvsflx(i,j)=0.e0

          end do
          end do

!$omp end do

        else if(fmois(1:5).eq.'moist') then

!$omp do schedule(runtime) private(i,j,a)

          do j=1,nj-1
          do i=1,ni-1
            a=.5e0*cm(i,j)/cm10(i,j)

            u10(i,j)=a*(u(i,j,2)+u(i+1,j,2))
            v10(i,j)=a*(v(i,j,2)+v(i,j+1,2))

            a=ch(i,j)/ch15(i,j)

            pt15(i,j)=ptv(i,j,1)+a*(ptv(i,j,2)-ptv(i,j,1))

            if(land(i,j).lt.0) then

              qv15(i,j)=qvsfc(i,j)+a*(qv(i,j,2)-qvsfc(i,j))             &
     &          *((1.e0+rddwkp*ch15(i,j))/(1.e0+rddwkp*ch(i,j)))

            else

              qv15(i,j)=qvsfc(i,j)+a*(qv(i,j,2)-qvsfc(i,j))

            end if

            pt15(i,j)=pt15(i,j)*(1.e0+qv15(i,j))/(1.e0+epsav*qv15(i,j))

            p15(i,j)=.5e0*(p(i,j,1)*(za(i,j)-1.5e0)                     &
     &        +p(i,j,2)*(za(i,j)+1.5e0))/za(i,j)

            if(land(i,j).eq.1) then
              tsfc(i,j)=kai(i,j)*tice(i,j)+(1.e0-kai(i,j))*tund(i,j,1)
            else
              tsfc(i,j)=tund(i,j,1)
            end if

            cdave(i,j)                                                  &
     &        =cdrat(1)*cdl(i,j)+cdrat(2)*cdm(i,j)+cdrat(3)*cdh(i,j)

            usflx(i,j)=-.5e0*(ufrc(i,j,1)+ufrc(i+1,j,1))
            vsflx(i,j)=-.5e0*(vfrc(i,j,1)+vfrc(i,j+1,1))

            ptsflx(i,j)=-ptfrc(i,j,1)
            qvsflx(i,j)=-qvfrc(i,j,1)

          end do
          end do

!$omp end do

        end if

!$omp end parallel

! -----

! Dump the surface monitor variables.

        if(dmplev.lt.10) then

          call outdmp2d('us    ',2,capt(48),ncpt(48),'xx',ni,nj,u10)
          call outdmp2d('vs    ',2,capt(49),ncpt(49),'xx',ni,nj,v10)

        else

          call setproj(idmpopt,idnspol,idtlat1,idtlat2,cpj)

          call s_rotuvm2s(idmpopt,idnspol,idtlon,1,ni-1,1,nj-1,1,1,cpj, &
     &                    0,ni+1,0,nj+1,1,1,lon,u10,v10)

          call outdmp2d('us    ',2,capt(53),ncpt(53),'xx',ni,nj,u10)
          call outdmp2d('vs    ',2,capt(54),ncpt(54),'xx',ni,nj,v10)

        end if

        call outdmp2d('ps    ',2,capt(50),ncpt(50),'xx',ni,nj,p15)
        call outdmp2d('pts   ',3,capt(51),ncpt(51),'xx',ni,nj,pt15)
        call outdmp2d('qvs   ',3,capt(52),ncpt(52),'xx',ni,nj,qv15)

        call outdmp2d('tgs   ',3,capt(55),ncpt(55),'xx',ni,nj,tsfc)

        call outdmp2d('hs    ',2,capt(56),ncpt(56),'xx',ni,nj,hs)
        call outdmp2d('le    ',2,capt(57),ncpt(57),'xx',ni,nj,le)

        call outdmp2d('rgd   ',3,capt(78),ncpt(78),'xx',ni,nj,rgd)

        call outdmp2d('rsd   ',3,capt(58),ncpt(58),'xx',ni,nj,rsd)
        call outdmp2d('rld   ',3,capt(59),ncpt(59),'xx',ni,nj,rld)
        call outdmp2d('rlu   ',3,capt(60),ncpt(60),'xx',ni,nj,rlu)

        call outdmp2d('cdl   ',3,capt(79),ncpt(79),'xx',ni,nj,cdl)
        call outdmp2d('cdm   ',3,capt(80),ncpt(80),'xx',ni,nj,cdm)
        call outdmp2d('cdh   ',3,capt(81),ncpt(81),'xx',ni,nj,cdh)

        call outdmp2d('cdave ',5,capt(61),ncpt(61),'xx',ni,nj,cdave)

        call outdmp2d('usflx ',5,capt(62),ncpt(62),'xx',ni,nj,usflx)
        call outdmp2d('vsflx ',5,capt(63),ncpt(63),'xx',ni,nj,vsflx)
        call outdmp2d('ptsflx',6,capt(64),ncpt(64),'xx',ni,nj,ptsflx)
        call outdmp2d('qvsflx',6,capt(65),ncpt(65),'xx',ni,nj,qvsflx)

! -----

      end if

!! -----

      end subroutine s_outpbl

!-----7--------------------------------------------------------------7--

      end module m_outpbl
