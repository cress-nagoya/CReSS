!***********************************************************************
      module m_adjstuv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/10/28
!     Modification: 2003/04/30, 2003/05/19, 2003/12/12, 2004/09/10,
!                   2005/04/04, 2006/04/03, 2006/11/06, 2006/12/04,
!                   2007/01/05, 2007/01/31, 2007/05/07, 2007/05/14,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2011/08/09, 2011/09/22, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     adjust the x and y components of velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_comphy
      use m_getiname
      use m_getrname
      use m_reducelb
      use m_reducevb

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: adjstuv, s_adjstuv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface adjstuv

        module procedure s_adjstuv

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
      subroutine s_adjstuv(fpwbc,fpebc,fpadvopt,fpmpopt,fpmfcopt,       &
     &                     fpdx,fpdy,fpdz,dtb,gtinc,area,ni,nj,nk,      &
     &                     rmf,rmf8u,rmf8v,rst8u,rst8v,ppp,ppf,         &
     &                     ugpv,utd,vgpv,vtd,uf,vf)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpdx
                       ! Formal parameter of unique index of dx

      integer, intent(in) :: fpdy
                       ! Formal parameter of unique index of dy

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: area(0:4)
                       ! Area of each boundary plane

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: rst8u(0:ni+1,0:nj+1,1:nk)
                       ! Bese state density x Jacobian at u points

      real, intent(in) :: rst8v(0:ni+1,0:nj+1,1:nk)
                       ! Bese state density x Jacobian at v points

      real, intent(in) :: ppp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at past

      real, intent(in) :: ppf(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation at future

      real, intent(in) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: utd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of x components of velocity
                       ! of GPV data

      real, intent(in) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: vtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of y components of velocity
                       ! of GPV data

! Input and output variables

      real, intent(inout) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(inout) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer advopt   ! Option for advection scheme
      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction

      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction

      real dx          ! Grid distance in x direction
      real dy          ! Grid distance in y direction
      real dz          ! Grid distance in z direction

      real dxdy        ! dx x dy
      real dxdz        ! dx x dz
      real dydz        ! dy x dz

      real tpdt        ! gtinc + dtb or gtinc + 2.0 x dtb

      real dpsp2       ! 2.0 x total difference
                       ! between top and bottom pressure at past

      real dpsf2       ! 2.0 x total difference
                       ! between top and bottom pressure at future

      real dflw        ! Flux on west boundary
      real dfle        ! Flux on east boundary
      real dfls        ! Flux on south boundary
      real dfln        ! Flux on north boundary

      real adj         ! Adjustment value

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpadvopt,advopt)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getrname(fpdx,dx)
      call getrname(fpdy,dy)
      call getrname(fpdz,dz)

! -----

! Set the common used variables.

      nkm1=nk-1
      nkm2=nk-2

      dxdy=dx*dy
      dxdz=dx*dz
      dydz=dy*dz

      if(advopt.le.3) then
        tpdt=gtinc+2.e0*dtb
      else
        tpdt=gtinc+dtb
      end if

! -----

! Get the maximum and minimum do loops index.

      if(ebw.eq.1.and.isub.eq.0) then
        istr=1
      else
        istr=2
      end if

      if(ebe.eq.1.and.isub.eq.nisub-1) then
        iend=ni-1
      else
        iend=ni-2
      end if

      if(ebs.eq.1.and.jsub.eq.0) then
        jstr=1
      else
        jstr=2
      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then
        jend=nj-1
      else
        jend=nj-2
      end if

! -----

! Initialize the sumed variables.

      dpsp2=0.e0
      dpsf2=0.e0

      dflw=0.e0
      dfle=0.e0
      dfls=0.e0
      dfln=0.e0

! -----

!! Adjust the x and y components of velocity.

! Calculate the 2.0 x total difference between bottom and top pressure
! and the flux on lateral boundary.

!$omp parallel default(shared)

      if(mfcopt.eq.0) then

!$omp do schedule(runtime) private(i,j) reduction(+: dpsp2,dpsf2)

        do j=jstr,jend
        do i=istr,iend

         dpsp2=dpsp2                                                    &
     &     +dxdy*((ppp(i,j,1)+ppp(i,j,2))-(ppp(i,j,nkm1)+ppp(i,j,nkm2)))

         dpsf2=dpsf2                                                    &
     &     +dxdy*((ppf(i,j,1)+ppf(i,j,2))-(ppf(i,j,nkm1)+ppf(i,j,nkm2)))

        end do
        end do

!$omp end do

      else

        if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

!$omp do schedule(runtime) private(i,j,a) reduction(+: dpsp2,dpsf2)

          do j=jstr,jend
          do i=istr,iend

            a=dxdy*rmf(i,j,2)

            dpsp2=dpsp2                                                 &
     &        +a*((ppp(i,j,1)+ppp(i,j,2))-(ppp(i,j,nkm1)+ppp(i,j,nkm2)))

            dpsf2=dpsf2                                                 &
     &        +a*((ppf(i,j,1)+ppf(i,j,2))-(ppf(i,j,nkm1)+ppf(i,j,nkm2)))

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,j,a) reduction(+: dpsp2,dpsf2)

          do j=jstr,jend
          do i=istr,iend

            a=dxdy*rmf(i,j,3)

            dpsp2=dpsp2                                                 &
     &        +a*((ppp(i,j,1)+ppp(i,j,2))-(ppp(i,j,nkm1)+ppp(i,j,nkm2)))

            dpsf2=dpsf2                                                 &
     &        +a*((ppf(i,j,1)+ppf(i,j,2))-(ppf(i,j,nkm1)+ppf(i,j,nkm2)))

          end do
          end do

!$omp end do

        end if

      end if

      if(ebw.eq.1.and.isub.eq.0.and.abs(wbc).ne.1) then

        if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

!$omp do schedule(runtime) private(j,k) reduction(+: dflw)

          do k=2,nk-2
          do j=jstr,jend
            dflw=dflw+dydz*rmf8u(1,j,2)*rst8u(1,j,k)                    &
     &        *(uf(1,j,k)-(ugpv(1,j,k)+utd(1,j,k)*tpdt))
          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(j,k) reduction(+: dflw)

          do k=2,nk-2
          do j=jstr,jend
            dflw=dflw+dydz*rst8u(1,j,k)                                 &
     &        *(uf(1,j,k)-(ugpv(1,j,k)+utd(1,j,k)*tpdt))
          end do
          end do

!$omp end do

        end if

      end if

      if(ebe.eq.1.and.isub.eq.nisub-1.and.abs(ebc).ne.1) then

        if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

!$omp do schedule(runtime) private(j,k) reduction(+: dfle)

          do k=2,nk-2
          do j=jstr,jend
            dfle=dfle+dydz*rmf8u(ni,j,2)*rst8u(ni,j,k)                  &
     &        *(uf(ni,j,k)-(ugpv(ni,j,k)+utd(ni,j,k)*tpdt))
          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(j,k) reduction(+: dfle)

          do k=2,nk-2
          do j=jstr,jend
            dfle=dfle+dydz*rst8u(ni,j,k)                                &
     &        *(uf(ni,j,k)-(ugpv(ni,j,k)+utd(ni,j,k)*tpdt))
          end do
          end do

!$omp end do

        end if

      end if

      if(ebs.eq.1.and.jsub.eq.0) then

        if(mfcopt.eq.1.and.mpopt.ne.5) then

!$omp do schedule(runtime) private(i,k) reduction(+: dfls)

          do k=2,nk-2
          do i=istr,iend
            dfls=dfls+dxdz*rmf8v(i,1,2)*rst8v(i,1,k)                    &
     &        *(vf(i,1,k)-(vgpv(i,1,k)+vtd(i,1,k)*tpdt))
          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,k) reduction(+: dfls)

          do k=2,nk-2
          do i=istr,iend
            dfls=dfls+dxdz*rst8v(i,1,k)                                 &
     &        *(vf(i,1,k)-(vgpv(i,1,k)+vtd(i,1,k)*tpdt))
          end do
          end do

!$omp end do

        end if

      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(mfcopt.eq.1.and.mpopt.ne.5) then

!$omp do schedule(runtime) private(i,k) reduction(+: dfln)

          do k=2,nk-2
          do i=istr,iend
            dfln=dfln+dxdz*rmf8v(i,nj,2)*rst8v(i,nj,k)                  &
     &        *(vf(i,nj,k)-(vgpv(i,nj,k)+vtd(i,nj,k)*tpdt))
          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,k) reduction(+: dfln)

          do k=2,nk-2
          do i=istr,iend
            dfln=dfln+dxdz*rst8v(i,nj,k)                                &
     &        *(vf(i,nj,k)-(vgpv(i,nj,k)+vtd(i,nj,k)*tpdt))
          end do
          end do

!$omp end do

        end if

      end if

!$omp end parallel

! -----

! Get the adjstment value.

      call reducevb(dpsp2)
      call reducevb(dpsf2)

      call reducelb(dflw,dfle,dfls,dfln)

      if(advopt.le.3) then

        adj=(.25e0*(dpsf2-dpsp2)/(g*dtb*area(0))                        &
     &    -((dflw-dfle)+(dfls-dfln)))/(area(1)+area(2)+area(3)+area(4))

      else

        adj=(.5e0*(dpsf2-dpsp2)/(g*dtb*area(0))                         &
     &    -((dflw-dfle)+(dfls-dfln)))/(area(1)+area(2)+area(3)+area(4))

      end if

! -----

! Finally adjust the x and y components of velocity.

!$omp parallel default(shared)

      if(ebw.eq.1.and.isub.eq.0.and.abs(wbc).ne.1) then

!$omp do schedule(runtime) private(j,k)

        do k=1,nk-1
        do j=jstr,jend
          uf(1,j,k)=uf(1,j,k)+adj/rst8u(1,j,k)
        end do
        end do

!$omp end do

      end if

      if(ebe.eq.1.and.isub.eq.nisub-1.and.abs(ebc).ne.1) then

!$omp do schedule(runtime) private(j,k)

        do k=1,nk-1
        do j=jstr,jend
          uf(ni,j,k)=uf(ni,j,k)-adj/rst8u(ni,j,k)
        end do
        end do

!$omp end do

      end if

      if(ebs.eq.1.and.jsub.eq.0) then

!$omp do schedule(runtime) private(i,k)

        do k=1,nk-1
        do i=istr,iend
          vf(i,1,k)=vf(i,1,k)+adj/rst8v(i,1,k)
        end do
        end do

!$omp end do

      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

!$omp do schedule(runtime) private(i,k)

        do k=1,nk-1
        do i=istr,iend
          vf(i,nj,k)=vf(i,nj,k)-adj/rst8v(i,nj,k)
        end do
        end do

!$omp end do

      end if

!$omp end parallel

! -----

!! -----

      end subroutine s_adjstuv

!-----7--------------------------------------------------------------7--

      end module m_adjstuv
