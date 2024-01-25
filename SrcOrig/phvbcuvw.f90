!***********************************************************************
      module m_phvbcuvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/07/13
!     Modification: 2001/08/07, 2001/09/13, 2001/12/11, 2002/04/02,
!                   2002/06/06, 2002/07/23, 2002/08/15, 2002/10/31,
!                   2003/01/04, 2003/02/13, 2003/03/13, 2003/04/30,
!                   2003/05/19, 2003/06/27, 2003/07/28, 2003/10/31,
!                   2003/11/05, 2003/11/28, 2003/12/12, 2004/04/01,
!                   2004/04/15, 2004/08/01, 2004/08/20, 2004/09/01,
!                   2006/04/03, 2006/09/21, 2006/11/06, 2006/12/04,
!                   2007/01/05, 2007/01/31, 2007/05/07, 2007/05/21,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/02/27, 2009/03/23, 2009/11/13, 2010/12/01,
!                   2011/01/19, 2011/08/09, 2011/09/22, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the diffrential phase speed term between the external
!     boundary and model grid for velocity variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_commpi
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: phvbcuvw, s_phvbcuvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface phvbcuvw

        module procedure s_phvbcuvw

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic max
      intrinsic min
      intrinsic mod
      intrinsic real
      intrinsic sign

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_phvbcuvw(fpexbvar,fpwbc,fpebc,fpsbc,fpnbc,           &
     &                      fpmpopt,fpmfcopt,fpdxiv,fpdyiv,fpgwave,     &
     &                      dtb,dts,gtinc,ni,nj,nk,rmf,rmf8u,rmf8v,     &
     &                      u,up,uf,v,vp,vf,w,wp,wf,ugpv,utd,vgpv,vtd,  &
     &                      wgpv,wtd,ucpx,ucpy,vcpx,vcpy,wcpx,wcpy,     &
     &                      u8v,v8u,cpavex,cpavey)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpexbvar
                       ! Formal parameter of unique index of exbvar

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpdxiv
                       ! Formal parameter of unique index of dxiv

      integer, intent(in) :: fpdyiv
                       ! Formal parameter of unique index of dyiv

      integer, intent(in) :: fpgwave
                       ! Formal parameter of unique index of gwave

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: dts
                       ! Small time steps interval

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(in) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(in) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(in) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(in) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(in) :: w(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at present

      real, intent(in) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

      real, intent(in) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

      real, intent(in) :: ugpv(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: utd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! x components of velocity of GPV data

      real, intent(in) :: vgpv(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: vtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! y components of velocity of GPV data

      real, intent(in) :: wgpv(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity of GPV data
                       ! at marked time

      real, intent(in) :: wtd(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! z components of velocity of GPV data

! Output variables

      real, intent(out) :: ucpx(1:nj,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on west and east boundary

      real, intent(out) :: ucpy(1:ni,1:nk,1:2)
                       ! Phase speed of x components of velocity
                       ! on south and north boundary

      real, intent(out) :: vcpx(1:nj,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on west and east boundary

      real, intent(out) :: vcpy(1:ni,1:nk,1:2)
                       ! Phase speed of y components of velocity
                       ! on south and north boundary

      real, intent(out) :: wcpx(1:nj,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on west and east boundary

      real, intent(out) :: wcpy(1:ni,1:nk,1:2)
                       ! Phase speed of z components of velocity
                       ! on south and north boundary

! Internal shared variables

      character(len=108) exbvar
                       ! Control flag of
                       ! extrenal boundary forced variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      integer istr     ! Minimum do loops index in x direction
      integer iend     ! Maximum do loops index in x direction
      integer jstr     ! Minimum do loops index in y direction
      integer jend     ! Maximum do loops index in y direction

      integer nim1     ! ni - 1
      integer nim2     ! ni - 2
      integer nim3     ! ni - 3

      integer njm1     ! nj - 1
      integer njm2     ! nj - 2
      integer njm3     ! nj - 3

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dx

      real gwave       ! Fastest gravity wave speed

      real gdxdt       ! gwave x dxiv x dtb
      real gdydt       ! gwave x dyiv x dtb

      real gdxdtn      ! - gwave x dxiv x dtb
      real gdydtn      ! - gwave x dyiv x dtb

      real dxdt        ! dxiv x dtb
      real dydt        ! dyiv x dtb

      real dxdt5       ! 0.5 x dxiv x dtb
      real dydt5       ! 0.5 x dyiv x dtb

      real gtinc0      ! gtinc
      real gtinc1      ! gtinc + dtb
      real gtinc2      ! gtinc + 2.0 x dtb

      real dtsdb       ! dts / dtb

      real nkm2v       ! 1.0 / real(nk - 2)
      real nkm3v       ! 1.0 / real(nk - 3)

      real, intent(inout) :: u8v(0:nj+1,1:nk)
                       ! x components of velocity at v points

      real, intent(inout) :: v8u(0:ni+1,1:nk)
                       ! y components of velocity at u points

      real, intent(inout) :: cpavex(0:nj+1)
                       ! Vertically averaged phase speed in x direction

      real, intent(inout) :: cpavey(0:ni+1)
                       ! Vertically averaged phase speed in y direction

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real bc0         ! Temporary variable
      real bc1         ! Temporary variable
      real bc2         ! Temporary variable

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(exbvar)

! -----

! Get the required namelist variables.

      call getcname(fpexbvar,exbvar)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)
      call getrname(fpgwave,gwave)

! -----

! Set the common used variables.

      nim1=ni-1
      nim2=ni-2
      nim3=ni-3

      njm1=nj-1
      njm2=nj-2
      njm3=nj-3

      if(ebw.eq.1.and.isub.eq.0) then
        istr=2
      else
        istr=1
      end if

      if(ebe.eq.1.and.isub.eq.nisub-1) then
        iend=ni-1
      else
        iend=ni
      end if

      if(ebs.eq.1.and.jsub.eq.0) then
        jstr=2
      else
        jstr=1
      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then
        jend=nj-1
      else
        jend=nj
      end if

      gdxdt=gwave*dtb*dxiv
      gdydt=gwave*dtb*dyiv

      gdxdtn=-gwave*dtb*dxiv
      gdydtn=-gwave*dtb*dyiv

      dxdt=dtb*dxiv
      dydt=dtb*dyiv

      dxdt5=.5e0*dtb*dxiv
      dydt5=.5e0*dtb*dyiv

      gtinc0=gtinc
      gtinc1=gtinc+dtb
      gtinc2=gtinc+2.e0*dtb

      dtsdb=dts/dtb

      nkm2v=1.e0/real(nk-2)
      nkm3v=1.e0/real(nk-3)

! -----

!!!! Calculate the diffrential phase speed term between the external
!!!! boundary and model grid velocity variables.

!$omp parallel default(shared) private(k)

!!! Calculate the differential phase speed term for x components of
!!! velocity.

!! Calculate the differential phase speed term for x components of
!! velocity on west and east boundary.

      if(exbvar(1:1).eq.'-') then

! Calculate the differential phase speed term for x components of
! velocity on west boundary.

        if(ebw.eq.1.and.isub.eq.0) then

          if(mod(wbc,10).eq.4.or.mod(wbc,10).eq.5) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j,bc0,bc1,bc2)

              do j=1,nj-1
                bc0=up(2,j,k)-(ugpv(2,j,k)+utd(2,j,k)*gtinc0)
                bc1=u(3,j,k)-(ugpv(3,j,k)+utd(3,j,k)*gtinc1)
                bc2=uf(2,j,k)-(ugpv(2,j,k)+utd(2,j,k)*gtinc2)

                ucpx(j,k,1)=bc2+bc0-2.e0*bc1

                if(abs(ucpx(j,k,1)).lt.eps) then

                  ucpx(j,k,1)=sign(eps,ucpx(j,k,1))

                end if

                ucpx(j,k,1)=min((bc2-bc0)/ucpx(j,k,1),gdxdtn)

              end do

!$omp end do

            end do

            if(mod(wbc,10).eq.5) then

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                cpavex(j)=0.e0
              end do

!$omp end do

              do k=2,nk-2

!$omp do schedule(runtime) private(j)

                do j=1,nj-1
                  cpavex(j)=cpavex(j)+ucpx(j,k,1)*nkm3v
                end do

!$omp end do

              end do

              do k=2,nk-2

!$omp do schedule(runtime) private(j)

                do j=1,nj-1
                  ucpx(j,k,1)=cpavex(j)
                end do

!$omp end do

              end do

            end if

          else if(mod(wbc,10).eq.6) then

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                ucpx(j,k,1)=min(u(2,j,k)*dxdt,gdxdtn)
              end do

!$omp end do

            end do

          else if(wbc.eq.7) then

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                ucpx(j,k,1)=gdxdtn
              end do

!$omp end do

            end do

          end if

          if(mfcopt.eq.1.and.mpopt.ne.5) then

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                ucpx(j,k,1)=max(ucpx(j,k,1),-rmf8u(2,j,2))
              end do

!$omp end do

            end do

          else

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                ucpx(j,k,1)=max(ucpx(j,k,1),-1.e0)
              end do

!$omp end do

            end do

          end if

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              ucpx(j,k,1)=ucpx(j,k,1)*dtsdb
            end do

!$omp end do

          end do

        end if

! -----

! Calculate the differential phase speed term for x components of
! velocity on east boundary.

        if(ebe.eq.1.and.isub.eq.nisub-1) then

          if(mod(ebc,10).eq.4.or.mod(ebc,10).eq.5) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j,bc0,bc1,bc2)

              do j=1,nj-1
                bc0=up(nim1,j,k)-(ugpv(nim1,j,k)+utd(nim1,j,k)*gtinc0)
                bc1=u(nim2,j,k)-(ugpv(nim2,j,k)+utd(nim2,j,k)*gtinc1)
                bc2=uf(nim1,j,k)-(ugpv(nim1,j,k)+utd(nim1,j,k)*gtinc2)

                ucpx(j,k,2)=2.e0*bc1-bc2-bc0

                if(abs(ucpx(j,k,2)).lt.eps) then

                  ucpx(j,k,2)=sign(eps,ucpx(j,k,2))

                end if

                ucpx(j,k,2)=max((bc2-bc0)/ucpx(j,k,2),gdxdt)

              end do

!$omp end do

            end do

            if(mod(ebc,10).eq.5) then

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                cpavex(j)=0.e0
              end do

!$omp end do

              do k=2,nk-2

!$omp do schedule(runtime) private(j)

                do j=1,nj-1
                  cpavex(j)=cpavex(j)+ucpx(j,k,2)*nkm3v
                end do

!$omp end do

              end do

              do k=2,nk-2

!$omp do schedule(runtime) private(j)

                do j=1,nj-1
                  ucpx(j,k,2)=cpavex(j)
                end do

!$omp end do

              end do

            end if

          else if(mod(ebc,10).eq.6) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                ucpx(j,k,2)=max(u(nim1,j,k)*dxdt,gdxdt)
              end do

!$omp end do

            end do

          else if(ebc.eq.7) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                ucpx(j,k,2)=gdxdt
              end do

!$omp end do

            end do

          end if

          if(mfcopt.eq.1.and.mpopt.ne.5) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                ucpx(j,k,2)=min(ucpx(j,k,2),rmf8u(nim1,j,2))
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                ucpx(j,k,2)=min(ucpx(j,k,2),1.e0)
              end do

!$omp end do

            end do

          end if

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              ucpx(j,k,2)=ucpx(j,k,2)*dtsdb
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

!! Calculate the differential phase speed term for x components of
!! velocity on south and north boundary.

      if(exbvar(1:1).eq.'-'.or.exbvar(1:1).eq.'+') then

! Calculate the differential phase speed term for x components of
! velocity on south boundary.

        if(ebs.eq.1.and.jsub.eq.0) then

          if(sbc.eq.4.or.sbc.eq.5) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,bc0,bc1,bc2)

              do i=1,ni
                bc0=up(i,2,k)-(ugpv(i,2,k)+utd(i,2,k)*gtinc0)
                bc1=u(i,3,k)-(ugpv(i,3,k)+utd(i,3,k)*gtinc1)
                bc2=uf(i,2,k)-(ugpv(i,2,k)+utd(i,2,k)*gtinc2)

                ucpy(i,k,1)=bc2+bc0-2.e0*bc1

                if(abs(ucpy(i,k,1)).lt.eps) then

                  ucpy(i,k,1)=sign(eps,ucpy(i,k,1))

                end if

                ucpy(i,k,1)=min((bc2-bc0)/ucpy(i,k,1),gdydtn)

              end do

!$omp end do

            end do

            if(sbc.eq.5) then

!$omp do schedule(runtime) private(i)

              do i=1,ni
                cpavey(i)=0.e0
              end do

!$omp end do

              do k=2,nk-2

!$omp do schedule(runtime) private(i)

                do i=1,ni
                  cpavey(i)=cpavey(i)+ucpy(i,k,1)*nkm3v
                end do

!$omp end do

              end do

              do k=2,nk-2

!$omp do schedule(runtime) private(i)

                do i=1,ni
                  ucpy(i,k,1)=cpavey(i)
                end do

!$omp end do

              end do

            end if

          else if(sbc.eq.6) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=istr,iend
                v8u(i,k)=v(i-1,2,k)+v(i,2,k)
              end do

!$omp end do

            end do

            if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime)

              do k=2,nk-2
                v8u(1,k)=v8u(2,k)
              end do

!$omp end do

            else if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime)

              do k=2,nk-2
                v8u(ni,k)=v8u(ni-1,k)
              end do

!$omp end do

            end if

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                ucpy(i,k,1)=min(v8u(i,k)*dydt5,gdydtn)
              end do

!$omp end do

            end do

          else if(sbc.ge.7) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                ucpy(i,k,1)=gdydtn
              end do

!$omp end do

            end do

          end if

          if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                ucpy(i,k,1)=max(ucpy(i,k,1),-rmf8u(i,2,2))
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                ucpy(i,k,1)=max(ucpy(i,k,1),-1.e0)
              end do

!$omp end do

            end do

          end if

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni
              ucpy(i,k,1)=ucpy(i,k,1)*dtsdb
            end do

!$omp end do

          end do

        end if

! -----

! Calculate the differential phase speed term for x components of
! velocity on north boundary.

        if(ebn.eq.1.and.jsub.eq.njsub-1) then

          if(nbc.eq.4.or.nbc.eq.5) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,bc0,bc1,bc2)

              do i=1,ni
                bc0=up(i,njm2,k)-(ugpv(i,njm2,k)+utd(i,njm2,k)*gtinc0)
                bc1=u(i,njm3,k)-(ugpv(i,njm3,k)+utd(i,njm3,k)*gtinc1)
                bc2=uf(i,njm2,k)-(ugpv(i,njm2,k)+utd(i,njm2,k)*gtinc2)

                ucpy(i,k,2)=2.e0*bc1-bc2-bc0

                if(abs(ucpy(i,k,2)).lt.eps) then

                  ucpy(i,k,2)=sign(eps,ucpy(i,k,2))

                end if

                ucpy(i,k,2)=max((bc2-bc0)/ucpy(i,k,2),gdydt)

              end do

!$omp end do

            end do

            if(nbc.eq.5) then

!$omp do schedule(runtime) private(i)

              do i=1,ni
                cpavey(i)=0.e0
              end do

!$omp end do

              do k=2,nk-2

!$omp do schedule(runtime) private(i)

                do i=1,ni
                  cpavey(i)=cpavey(i)+ucpy(i,k,2)*nkm3v
                end do

!$omp end do

              end do

              do k=2,nk-2

!$omp do schedule(runtime) private(i)

                do i=1,ni
                  ucpy(i,k,2)=cpavey(i)
                end do

!$omp end do

              end do

            end if

          else if(nbc.eq.6) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=istr,iend
                v8u(i,k)=v(i-1,njm1,k)+v(i,njm1,k)
              end do

!$omp end do

            end do

            if(ebw.eq.1.and.isub.eq.0) then

!$omp do schedule(runtime)

              do k=2,nk-2
                v8u(1,k)=v8u(2,k)
              end do

!$omp end do

            else if(ebe.eq.1.and.isub.eq.nisub-1) then

!$omp do schedule(runtime)

              do k=2,nk-2
                v8u(ni,k)=v8u(ni-1,k)
              end do

!$omp end do

            end if

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                ucpy(i,k,2)=max(v8u(i,k)*dydt5,gdydt)
              end do

!$omp end do

            end do

          else if(nbc.ge.7) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                ucpy(i,k,2)=gdydt
              end do

!$omp end do

            end do

          end if

          if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                ucpy(i,k,2)=min(ucpy(i,k,2),rmf8u(i,njm2,2))
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni
                ucpy(i,k,2)=min(ucpy(i,k,2),1.e0)
              end do

!$omp end do

            end do

          end if

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni
              ucpy(i,k,2)=ucpy(i,k,2)*dtsdb
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

!!! -----

!!! Calculate the differential phase speed term for y components of
!!! velocity.

!! Calculate the differential phase speed term for y components of
!! velocity on west and east boundary.

      if(exbvar(2:2).eq.'-'.or.exbvar(2:2).eq.'+') then

! Calculate the differential phase speed term for y components of
! velocity on west boundary.

        if(ebw.eq.1.and.isub.eq.0) then

          if(wbc.eq.4.or.wbc.eq.5) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j,bc0,bc1,bc2)

              do j=1,nj
                bc0=vp(2,j,k)-(vgpv(2,j,k)+vtd(2,j,k)*gtinc0)
                bc1=v(3,j,k)-(vgpv(3,j,k)+vtd(3,j,k)*gtinc1)
                bc2=vf(2,j,k)-(vgpv(2,j,k)+vtd(2,j,k)*gtinc2)

                vcpx(j,k,1)=bc2+bc0-2.e0*bc1

                if(abs(vcpx(j,k,1)).lt.eps) then

                  vcpx(j,k,1)=sign(eps,vcpx(j,k,1))

                end if

                vcpx(j,k,1)=min((bc2-bc0)/vcpx(j,k,1),gdxdtn)

              end do

!$omp end do

            end do

            if(wbc.eq.5) then

!$omp do schedule(runtime) private(j)

              do j=1,nj
                cpavex(j)=0.e0
              end do

!$omp end do

              do k=2,nk-2

!$omp do schedule(runtime) private(j)

                do j=1,nj
                  cpavex(j)=cpavex(j)+vcpx(j,k,1)*nkm3v
                end do

!$omp end do

              end do

              do k=2,nk-2

!$omp do schedule(runtime) private(j)

                do j=1,nj
                  vcpx(j,k,1)=cpavex(j)
                end do

!$omp end do

              end do

            end if

          else if(wbc.eq.6) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=jstr,jend
                u8v(j,k)=u(2,j-1,k)+u(2,j,k)
              end do

!$omp end do

            end do

            if(ebs.eq.1.and.jsub.eq.0) then

!$omp do schedule(runtime)

              do k=2,nk-2
                u8v(1,k)=u8v(2,k)
              end do

!$omp end do

            else if(ebn.eq.1.and.jsub.eq.njsub-1) then

!$omp do schedule(runtime)

              do k=2,nk-2
                u8v(nj,k)=u8v(nj-1,k)
              end do

!$omp end do

            end if

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj
                vcpx(j,k,1)=min(u8v(j,k)*dxdt5,gdxdtn)
              end do

!$omp end do

            end do

          else if(wbc.ge.7) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj
                vcpx(j,k,1)=gdxdtn
              end do

!$omp end do

            end do

          end if

          if(mfcopt.eq.1.and.mpopt.ne.5) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj
                vcpx(j,k,1)=max(vcpx(j,k,1),-rmf8v(2,j,2))
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj
                vcpx(j,k,1)=max(vcpx(j,k,1),-1.e0)
              end do

!$omp end do

            end do

          end if

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj
              vcpx(j,k,1)=vcpx(j,k,1)*dtsdb
            end do

!$omp end do

          end do

        end if

! -----

! Calculate the differential phase speed term for y components of
! velocity on east boundary.

        if(ebe.eq.1.and.isub.eq.nisub-1) then

          if(ebc.eq.4.or.ebc.eq.5) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j,bc0,bc1,bc2)

              do j=1,nj
                bc0=vp(nim2,j,k)-(vgpv(nim2,j,k)+vtd(nim2,j,k)*gtinc0)
                bc1=v(nim3,j,k)-(vgpv(nim3,j,k)+vtd(nim3,j,k)*gtinc1)
                bc2=vf(nim2,j,k)-(vgpv(nim2,j,k)+vtd(nim2,j,k)*gtinc2)

                vcpx(j,k,2)=2.e0*bc1-bc2-bc0

                if(abs(vcpx(j,k,2)).lt.eps) then

                  vcpx(j,k,2)=sign(eps,vcpx(j,k,2))

                end if

                vcpx(j,k,2)=max((bc2-bc0)/vcpx(j,k,2),gdxdt)

              end do

!$omp end do

            end do

            if(ebc.eq.5) then

!$omp do schedule(runtime) private(j)

              do j=1,nj
                cpavex(j)=0.e0
              end do

!$omp end do

              do k=2,nk-2

!$omp do schedule(runtime) private(j)

                do j=1,nj
                  cpavex(j)=cpavex(j)+vcpx(j,k,2)*nkm3v
                end do

!$omp end do

              end do

              do k=2,nk-2

!$omp do schedule(runtime) private(j)

                do j=1,nj
                  vcpx(j,k,2)=cpavex(j)
                end do

!$omp end do

              end do

            end if

          else if(ebc.eq.6) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=jstr,jend
                u8v(j,k)=u(nim1,j-1,k)+u(nim1,j,k)
              end do

!$omp end do

            end do

            if(ebs.eq.1.and.jsub.eq.0) then

!$omp do schedule(runtime)

              do k=2,nk-2
                u8v(1,k)=u8v(2,k)
              end do

!$omp end do

            else if(ebn.eq.1.and.jsub.eq.njsub-1) then

!$omp do schedule(runtime)

              do k=2,nk-2
                u8v(nj,k)=u8v(nj-1,k)
              end do

!$omp end do

            end if

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj
                vcpx(j,k,2)=max(u8v(j,k)*dxdt5,gdxdt)
              end do

!$omp end do

            end do

          else if(ebc.ge.7) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj
                vcpx(j,k,2)=gdxdt
              end do

!$omp end do

            end do

          end if

          if(mfcopt.eq.1.and.mpopt.ne.5) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj
                vcpx(j,k,2)=min(vcpx(j,k,2),rmf8v(nim2,j,2))
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj
                vcpx(j,k,2)=min(vcpx(j,k,2),1.e0)
              end do

!$omp end do

            end do

          end if

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj
              vcpx(j,k,2)=vcpx(j,k,2)*dtsdb
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

!! Calculate the differential phase speed term for y components of
!! velocity on south and north boundary.

      if(exbvar(2:2).eq.'-') then

! Calculate the differential phase speed term for y components of
! velocity on south boundary.

        if(ebs.eq.1.and.jsub.eq.0) then

          if(mod(sbc,10).eq.4.or.mod(sbc,10).eq.5) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,bc0,bc1,bc2)

              do i=1,ni-1
                bc0=vp(i,2,k)-(vgpv(i,2,k)+vtd(i,2,k)*gtinc0)
                bc1=v(i,3,k)-(vgpv(i,3,k)+vtd(i,3,k)*gtinc1)
                bc2=vf(i,2,k)-(vgpv(i,2,k)+vtd(i,2,k)*gtinc2)

                vcpy(i,k,1)=bc2+bc0-2.e0*bc1

                if(abs(vcpy(i,k,1)).lt.eps) then

                  vcpy(i,k,1)=sign(eps,vcpy(i,k,1))

                end if

                vcpy(i,k,1)=min((bc2-bc0)/vcpy(i,k,1),gdydtn)

              end do

!$omp end do

            end do

            if(mod(sbc,10).eq.5) then

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                cpavey(i)=0.e0
              end do

!$omp end do

              do k=2,nk-2

!$omp do schedule(runtime) private(i)

                do i=1,ni-1
                  cpavey(i)=cpavey(i)+vcpy(i,k,1)*nkm3v
                end do

!$omp end do

              end do

              do k=2,nk-2

!$omp do schedule(runtime) private(i)

                do i=1,ni-1
                  vcpy(i,k,1)=cpavey(i)
                end do

!$omp end do

              end do

            end if

          else if(mod(sbc,10).eq.6) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                vcpy(i,k,1)=min(v(i,2,k)*dydt,gdydtn)
              end do

!$omp end do

            end do

          else if(sbc.eq.7) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                vcpy(i,k,1)=gdydtn
              end do

!$omp end do

            end do

          end if

          if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                vcpy(i,k,1)=max(vcpy(i,k,1),-rmf8v(i,2,2))
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                vcpy(i,k,1)=max(vcpy(i,k,1),-1.e0)
              end do

!$omp end do

            end do

          end if

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              vcpy(i,k,1)=vcpy(i,k,1)*dtsdb
            end do

!$omp end do

          end do

        end if

! -----

! Calculate the differential phase speed term for y components of
! velocity on north boundary.

        if(ebn.eq.1.and.jsub.eq.njsub-1) then

          if(mod(nbc,10).eq.4.or.mod(nbc,10).eq.5) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,bc0,bc1,bc2)

              do i=1,ni-1
                bc0=vp(i,njm1,k)-(vgpv(i,njm1,k)+vtd(i,njm1,k)*gtinc0)
                bc1=v(i,njm2,k)-(vgpv(i,njm2,k)+vtd(i,njm2,k)*gtinc1)
                bc2=vf(i,njm1,k)-(vgpv(i,njm1,k)+vtd(i,njm1,k)*gtinc2)

                vcpy(i,k,2)=2.e0*bc1-bc2-bc0

                if(abs(vcpy(i,k,2)).lt.eps) then

                  vcpy(i,k,2)=sign(eps,vcpy(i,k,2))

                end if

                vcpy(i,k,2)=max((bc2-bc0)/vcpy(i,k,2),gdydt)

              end do

!$omp end do

            end do

            if(mod(nbc,10).eq.5) then

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                cpavey(i)=0.e0
              end do

!$omp end do

              do k=2,nk-2

!$omp do schedule(runtime) private(i)

                do i=1,ni-1
                  cpavey(i)=cpavey(i)+vcpy(i,k,2)*nkm3v
                end do

!$omp end do

              end do

              do k=2,nk-2

!$omp do schedule(runtime) private(i)

                do i=1,ni-1
                  vcpy(i,k,2)=cpavey(i)
                end do

!$omp end do

              end do

            end if

          else if(mod(nbc,10).eq.6) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                vcpy(i,k,2)=max(v(i,njm1,k)*dydt,gdydt)
              end do

!$omp end do

            end do

          else if(nbc.eq.7) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                vcpy(i,k,2)=gdydt
              end do

!$omp end do

            end do

          end if

          if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                vcpy(i,k,2)=min(vcpy(i,k,2),rmf8v(i,njm1,2))
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                vcpy(i,k,2)=min(vcpy(i,k,2),1.e0)
              end do

!$omp end do

            end do

          end if

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              vcpy(i,k,2)=vcpy(i,k,2)*dtsdb
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

!!! -----

!! Calculate the differential phase speed term for z components of
!! velocity.

      if(exbvar(3:3).eq.'-') then

! Calculate the differential phase speed term for z components of
! velocity on west boundary.

        if(ebw.eq.1.and.isub.eq.0) then

          if(wbc.eq.4.or.wbc.eq.5) then

            do k=2,nk-1

!$omp do schedule(runtime) private(j,bc0,bc1,bc2)

              do j=1,nj-1
                bc0=wp(2,j,k)-(wgpv(2,j,k)+wtd(2,j,k)*gtinc0)
                bc1=w(3,j,k)-(wgpv(3,j,k)+wtd(3,j,k)*gtinc1)
                bc2=wf(2,j,k)-(wgpv(2,j,k)+wtd(2,j,k)*gtinc2)

                wcpx(j,k,1)=bc2+bc0-2.e0*bc1

                if(abs(wcpx(j,k,1)).lt.eps) then

                  wcpx(j,k,1)=sign(eps,wcpx(j,k,1))

                end if

                wcpx(j,k,1)=min((bc2-bc0)/wcpx(j,k,1),gdxdtn)

              end do

!$omp end do

            end do

            if(wbc.eq.5) then

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                cpavex(j)=0.e0
              end do

!$omp end do

              do k=2,nk-1

!$omp do schedule(runtime) private(j)

                do j=1,nj-1
                  cpavex(j)=cpavex(j)+wcpx(j,k,1)*nkm2v
                end do

!$omp end do

              end do

              do k=2,nk-1

!$omp do schedule(runtime) private(j)

                do j=1,nj-1
                  wcpx(j,k,1)=cpavex(j)
                end do

!$omp end do

              end do

            end if

          else if(wbc.eq.6) then

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                wcpx(j,k,1)=min((u(2,j,k-1)+u(2,j,k))*dxdt5,gdxdtn)
              end do

!$omp end do

            end do

          else if(wbc.ge.7) then

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                wcpx(j,k,1)=gdxdtn
              end do

!$omp end do

            end do

          end if

          if(mfcopt.eq.1.and.mpopt.ne.5) then

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                wcpx(j,k,1)=max(wcpx(j,k,1),-rmf(2,j,2))
              end do

!$omp end do

            end do

          else

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                wcpx(j,k,1)=max(wcpx(j,k,1),-1.e0)
              end do

!$omp end do

            end do

          end if

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wcpx(j,k,1)=wcpx(j,k,1)*dtsdb
            end do

!$omp end do

          end do

        end if

! -----

! Calculate the differential phase speed term for z components of
! velocity on east boundary.

        if(ebe.eq.1.and.isub.eq.nisub-1) then

          if(ebc.eq.4.or.ebc.eq.5) then

            do k=2,nk-1

!$omp do schedule(runtime) private(j,bc0,bc1,bc2)

              do j=1,nj-1
                bc0=wp(nim2,j,k)-(wgpv(nim2,j,k)+wtd(nim2,j,k)*gtinc0)
                bc1=w(nim3,j,k)-(wgpv(nim3,j,k)+wtd(nim3,j,k)*gtinc1)
                bc2=wf(nim2,j,k)-(wgpv(nim2,j,k)+wtd(nim2,j,k)*gtinc2)

                wcpx(j,k,2)=2.e0*bc1-bc2-bc0

                if(abs(wcpx(j,k,2)).lt.eps) then

                  wcpx(j,k,2)=sign(eps,wcpx(j,k,2))

                end if

                wcpx(j,k,2)=max((bc2-bc0)/wcpx(j,k,2),gdxdt)

              end do

!$omp end do

            end do

            if(ebc.eq.5) then

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                cpavex(j)=0.e0
              end do

!$omp end do

              do k=2,nk-1

!$omp do schedule(runtime) private(j)

                do j=1,nj-1
                  cpavex(j)=cpavex(j)+wcpx(j,k,2)*nkm2v
                end do

!$omp end do

              end do

              do k=2,nk-1

!$omp do schedule(runtime) private(j)

                do j=1,nj-1
                  wcpx(j,k,2)=cpavex(j)
                end do

!$omp end do

              end do

            end if

          else if(ebc.eq.6) then

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                wcpx(j,k,2)=max((u(nim1,j,k-1)+u(nim1,j,k))*dxdt5,gdxdt)
              end do

!$omp end do

            end do

          else if(ebc.ge.7) then

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                wcpx(j,k,2)=gdxdt
              end do

!$omp end do

            end do

          end if

          if(mfcopt.eq.1.and.mpopt.ne.5) then

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                wcpx(j,k,2)=min(wcpx(j,k,2),rmf(nim2,j,2))
              end do

!$omp end do

            end do

          else

            do k=2,nk-1

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                wcpx(j,k,2)=min(wcpx(j,k,2),1.e0)
              end do

!$omp end do

            end do

          end if

          do k=2,nk-1

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              wcpx(j,k,2)=wcpx(j,k,2)*dtsdb
            end do

!$omp end do

          end do

        end if

! -----

! Calculate the differential phase speed term for z components of
! velocity on south boundary.

        if(ebs.eq.1.and.jsub.eq.0) then

          if(sbc.eq.4.or.sbc.eq.5) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i,bc0,bc1,bc2)

              do i=1,ni-1
                bc0=wp(i,2,k)-(wgpv(i,2,k)+wtd(i,2,k)*gtinc0)
                bc1=w(i,3,k)-(wgpv(i,3,k)+wtd(i,3,k)*gtinc1)
                bc2=wf(i,2,k)-(wgpv(i,2,k)+wtd(i,2,k)*gtinc2)

                wcpy(i,k,1)=bc2+bc0-2.e0*bc1

                if(abs(wcpy(i,k,1)).lt.eps) then

                  wcpy(i,k,1)=sign(eps,wcpy(i,k,1))

                end if

                wcpy(i,k,1)=min((bc2-bc0)/wcpy(i,k,1),gdydtn)

              end do

!$omp end do

            end do

            if(sbc.eq.5) then

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                cpavey(i)=0.e0
              end do

!$omp end do

              do k=2,nk-1

!$omp do schedule(runtime) private(i)

                do i=1,ni-1
                  cpavey(i)=cpavey(i)+wcpy(i,k,1)*nkm2v
                end do

!$omp end do

              end do

              do k=2,nk-1

!$omp do schedule(runtime) private(i)

                do i=1,ni-1
                  wcpy(i,k,1)=cpavey(i)
                end do

!$omp end do

              end do

            end if

          else if(sbc.eq.6) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                wcpy(i,k,1)=min((v(i,2,k-1)+v(i,2,k))*dydt5,gdydtn)
              end do

!$omp end do

            end do

          else if(sbc.ge.7) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                wcpy(i,k,1)=gdydtn
              end do

!$omp end do

            end do

          end if

          if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                wcpy(i,k,1)=max(wcpy(i,k,1),-rmf(i,2,2))
              end do

!$omp end do

            end do

          else

            do k=2,nk-1

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                wcpy(i,k,1)=max(wcpy(i,k,1),-1.e0)
              end do

!$omp end do

            end do

          end if

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wcpy(i,k,1)=wcpy(i,k,1)*dtsdb
            end do

!$omp end do

          end do

        end if

! -----

! Calculate the differential phase speed term for z components of
! velocity on north boundary.

        if(ebn.eq.1.and.jsub.eq.njsub-1) then

          if(nbc.eq.4.or.nbc.eq.5) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i,bc0,bc1,bc2)

              do i=1,ni-1
                bc0=wp(i,njm2,k)-(wgpv(i,njm2,k)+wtd(i,njm2,k)*gtinc0)
                bc1=w(i,njm3,k)-(wgpv(i,njm3,k)+wtd(i,njm3,k)*gtinc1)
                bc2=wf(i,njm2,k)-(wgpv(i,njm2,k)+wtd(i,njm2,k)*gtinc2)

                wcpy(i,k,2)=2.e0*bc1-bc2-bc0

                if(abs(wcpy(i,k,2)).lt.eps) then

                  wcpy(i,k,2)=sign(eps,wcpy(i,k,2))

                end if

                wcpy(i,k,2)=max((bc2-bc0)/wcpy(i,k,2),gdydt)

              end do

!$omp end do

            end do

            if(nbc.eq.5) then

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                cpavey(i)=0.e0
              end do

!$omp end do

              do k=2,nk-1

!$omp do schedule(runtime) private(i)

                do i=1,ni-1
                  cpavey(i)=cpavey(i)+wcpy(i,k,2)*nkm2v
                end do

!$omp end do

              end do

              do k=2,nk-1

!$omp do schedule(runtime) private(i)

                do i=1,ni-1
                  wcpy(i,k,2)=cpavey(i)
                end do

!$omp end do

              end do

            end if

          else if(nbc.eq.6) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                wcpy(i,k,2)=max((v(i,njm1,k-1)+v(i,njm1,k))*dydt5,gdydt)
              end do

!$omp end do

            end do

          else if(nbc.ge.7) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                wcpy(i,k,2)=gdydt
              end do

!$omp end do

            end do

          end if

          if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                wcpy(i,k,2)=min(wcpy(i,k,2),rmf(i,njm2,2))
              end do

!$omp end do

            end do

          else

            do k=2,nk-1

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                wcpy(i,k,2)=min(wcpy(i,k,2),1.e0)
              end do

!$omp end do

            end do

          end if

          do k=2,nk-1

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              wcpy(i,k,2)=wcpy(i,k,2)*dtsdb
            end do

!$omp end do

          end do

        end if

! -----

      end if

!! -----

!$omp end parallel

!!!! -----

      end subroutine s_phvbcuvw

!-----7--------------------------------------------------------------7--

      end module m_phvbcuvw
