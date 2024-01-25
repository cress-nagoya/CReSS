!***********************************************************************
      module m_phvbcs
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/07/13
!     Modification: 2001/08/07, 2001/09/13, 2001/12/11, 2002/04/02,
!                   2002/06/06, 2002/07/23, 2002/08/15, 2002/10/31,
!                   2003/01/04, 2003/02/13, 2003/03/13, 2003/04/30,
!                   2003/05/19, 2003/10/31, 2003/11/05, 2003/11/28,
!                   2003/12/12, 2004/04/15, 2004/08/20, 2004/09/01,
!                   2006/04/03, 2006/11/06, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/01/31, 2007/05/07, 2007/05/21,
!                   2007/10/19, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2008/10/10, 2009/02/27, 2009/03/23, 2009/11/13,
!                   2010/12/01, 2011/01/19, 2011/07/15, 2011/08/09,
!                   2011/09/22, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the diffrential phase speed term between the external
!     boundary and model grid for scalar variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_commpi
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: phvbcs, s_phvbcs

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface phvbcs

        module procedure s_phvbcs

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic max
      intrinsic min
      intrinsic real
      intrinsic sign

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_phvbcs(fpwbc,fpebc,fpsbc,fpnbc,fpadvopt,             &
     &                    fpmpopt,fpmfcopt,fpdxiv,fpdyiv,fpgwave,       &
     &                    fproc,dtb,dts,dtsep,gtinc,ni,nj,nk,rmf,       &
     &                    u,v,s,sp,sf,sgpv,std,scpx,scpy,cpavex,cpavey)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

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

      real, intent(in) :: dtsep
                       ! Time steps interval
                       ! of vertical Cubic Lagrange advection

      real, intent(in) :: gtinc
                       ! Lapse of forecast time from GPV data reading

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: u(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at present

      real, intent(in) :: v(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at present

      real, intent(in) :: s(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at present

      real, intent(in) :: sp(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at past

      real, intent(in) :: sf(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at future

      real, intent(in) :: sgpv(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable of GPV data
                       ! at marked time

      real, intent(in) :: std(0:ni+1,0:nj+1,1:nk)
                       ! Time tendency of
                       ! optional scalar variable of GPV data

! Output variables

      real, intent(out) :: scpx(1:nj,1:nk,1:2)
                       ! Phase speed of optional scalar variable
                       ! on west and east boundary

      real, intent(out) :: scpy(1:ni,1:nk,1:2)
                       ! Phase speed of optional scalar variable
                       ! on south and north boundary

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer advopt   ! Option for advection scheme
      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

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

      real gtinc0      ! gtinc
      real gtinc1      ! gtinc + dtb
      real gtinc2      ! gtinc + 2.0 x dtb

      real dtdvb       ! dts / dtb or dtsep / dtb

      real nkm3v       ! 1.0 / real(nk - 3)

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

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)
      call getiname(fpadvopt,advopt)
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

      gdxdt=gwave*dtb*dxiv
      gdydt=gwave*dtb*dyiv

      gdxdtn=-gwave*dtb*dxiv
      gdydtn=-gwave*dtb*dyiv

      dxdt=dtb*dxiv
      dydt=dtb*dyiv

      gtinc0=gtinc
      gtinc1=gtinc+dtb
      gtinc2=gtinc+2.e0*dtb

      if(fproc(1:3).eq.'sml') then
        dtdvb=dts/dtb
      else

        if(advopt.ge.4) then
          dtdvb=dtsep/dtb
        end if

      end if

      nkm3v=1.e0/real(nk-3)

! -----

!! Calculate the diffrential phase speed term between the external
!! boundary and model grid scalar variables.

!$omp parallel default(shared) private(k)

! Calculate the differential phase speed term for optional scalar
! variable on west boundary.

      if(ebw.eq.1.and.isub.eq.0) then

        if(wbc.eq.4.or.wbc.eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j,bc0,bc1,bc2)

            do j=1,nj-1
              bc0=sp(2,j,k)-(sgpv(2,j,k)+std(2,j,k)*gtinc0)
              bc1=s(3,j,k)-(sgpv(3,j,k)+std(3,j,k)*gtinc1)
              bc2=sf(2,j,k)-(sgpv(2,j,k)+std(2,j,k)*gtinc2)

              scpx(j,k,1)=bc2+bc0-2.e0*bc1

              if(abs(scpx(j,k,1)).lt.eps) then

                scpx(j,k,1)=sign(eps,scpx(j,k,1))

              end if

              scpx(j,k,1)=min((bc2-bc0)/scpx(j,k,1),gdxdtn)

            end do

!$omp end do

          end do

          if(wbc.eq.5) then

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              cpavex(j)=0.e0
            end do

!$omp end do

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                cpavex(j)=cpavex(j)+scpx(j,k,1)*nkm3v
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                scpx(j,k,1)=cpavex(j)
              end do

!$omp end do

            end do

          end if

        else if(wbc.eq.6) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              scpx(j,k,1)=min(u(2,j,k)*dxdt,gdxdtn)
            end do

!$omp end do

          end do

        else if(wbc.ge.7) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              scpx(j,k,1)=gdxdtn
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.mpopt.ne.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              scpx(j,k,1)=max(scpx(j,k,1),-rmf(2,j,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              scpx(j,k,1)=max(scpx(j,k,1),-1.e0)
            end do

!$omp end do

          end do

        end if

        if(fproc(1:3).eq.'sml') then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              scpx(j,k,1)=scpx(j,k,1)*dtdvb
            end do

!$omp end do

          end do

        else

          if(advopt.ge.4) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                scpx(j,k,1)=scpx(j,k,1)*dtdvb
              end do

!$omp end do

            end do

          end if

        end if

      end if

! -----

! Calculate the differential phase speed term for optional scalar
! variable on east boundary.

      if(ebe.eq.1.and.isub.eq.nisub-1) then

        if(ebc.eq.4.or.ebc.eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j,bc0,bc1,bc2)

            do j=1,nj-1
              bc0=sp(nim2,j,k)-(sgpv(nim2,j,k)+std(nim2,j,k)*gtinc0)
              bc1=s(nim3,j,k)-(sgpv(nim3,j,k)+std(nim3,j,k)*gtinc1)
              bc2=sf(nim2,j,k)-(sgpv(nim2,j,k)+std(nim2,j,k)*gtinc2)

              scpx(j,k,2)=2.e0*bc1-bc2-bc0

              if(abs(scpx(j,k,2)).lt.eps) then

                scpx(j,k,2)=sign(eps,scpx(j,k,2))

              end if

              scpx(j,k,2)=max((bc2-bc0)/scpx(j,k,2),gdxdt)

            end do

!$omp end do

          end do

          if(ebc.eq.5) then

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              cpavex(j)=0.e0
            end do

!$omp end do

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                cpavex(j)=cpavex(j)+scpx(j,k,2)*nkm3v
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                scpx(j,k,2)=cpavex(j)
              end do

!$omp end do

            end do

          end if

        else if(ebc.eq.6) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              scpx(j,k,2)=max(u(nim1,j,k)*dxdt,gdxdt)
            end do

!$omp end do

          end do

        else if(ebc.ge.7) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              scpx(j,k,2)=gdxdt
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.mpopt.ne.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              scpx(j,k,2)=min(scpx(j,k,2),rmf(nim2,j,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              scpx(j,k,2)=min(scpx(j,k,2),1.e0)
            end do

!$omp end do

          end do

        end if

        if(fproc(1:3).eq.'sml') then

          do k=2,nk-2

!$omp do schedule(runtime) private(j)

            do j=1,nj-1
              scpx(j,k,2)=scpx(j,k,2)*dtdvb
            end do

!$omp end do

          end do

        else

          if(advopt.ge.4) then

            do k=2,nk-2

!$omp do schedule(runtime) private(j)

              do j=1,nj-1
                scpx(j,k,2)=scpx(j,k,2)*dtdvb
              end do

!$omp end do

            end do

          end if

        end if

      end if

! -----

! Calculate the differential phase speed term for optional scalar
! variable on south boundary.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(sbc.eq.4.or.sbc.eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,bc0,bc1,bc2)

            do i=1,ni-1
              bc0=sp(i,2,k)-(sgpv(i,2,k)+std(i,2,k)*gtinc0)
              bc1=s(i,3,k)-(sgpv(i,3,k)+std(i,3,k)*gtinc1)
              bc2=sf(i,2,k)-(sgpv(i,2,k)+std(i,2,k)*gtinc2)

              scpy(i,k,1)=bc2+bc0-2.e0*bc1

              if(abs(scpy(i,k,1)).lt.eps) then

                scpy(i,k,1)=sign(eps,scpy(i,k,1))

              end if

              scpy(i,k,1)=min((bc2-bc0)/scpy(i,k,1),gdydtn)

            end do

!$omp end do

          end do

          if(sbc.eq.5) then

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              cpavey(i)=0.e0
            end do

!$omp end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                cpavey(i)=cpavey(i)+scpy(i,k,1)*nkm3v
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                scpy(i,k,1)=cpavey(i)
              end do

!$omp end do

            end do

          end if

        else if(sbc.eq.6) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              scpy(i,k,1)=min(v(i,2,k)*dydt,gdydtn)
            end do

!$omp end do

          end do

        else if(sbc.ge.7) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              scpy(i,k,1)=gdydtn
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              scpy(i,k,1)=max(scpy(i,k,1),-rmf(i,2,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              scpy(i,k,1)=max(scpy(i,k,1),-1.e0)
            end do

!$omp end do

          end do

        end if

        if(fproc(1:3).eq.'sml') then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              scpy(i,k,1)=scpy(i,k,1)*dtdvb
            end do

!$omp end do

          end do

        else

          if(advopt.ge.4) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                scpy(i,k,1)=scpy(i,k,1)*dtdvb
              end do

!$omp end do

            end do

          end if

        end if

      end if

! -----

! Calculate the differential phase speed term for optional scalar
! variable on north boundary.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(nbc.eq.4.or.nbc.eq.5) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,bc0,bc1,bc2)

            do i=1,ni-1
              bc0=sp(i,njm2,k)-(sgpv(i,njm2,k)+std(i,njm2,k)*gtinc0)
              bc1=s(i,njm3,k)-(sgpv(i,njm3,k)+std(i,njm3,k)*gtinc1)
              bc2=sf(i,njm2,k)-(sgpv(i,njm2,k)+std(i,njm2,k)*gtinc2)

              scpy(i,k,2)=2.e0*bc1-bc2-bc0

              if(abs(scpy(i,k,2)).lt.eps) then

                scpy(i,k,2)=sign(eps,scpy(i,k,2))

              end if

              scpy(i,k,2)=max((bc2-bc0)/scpy(i,k,2),gdydt)

            end do

!$omp end do

          end do

          if(nbc.eq.5) then

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              cpavey(i)=0.e0
            end do

!$omp end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                cpavey(i)=cpavey(i)+scpy(i,k,2)*nkm3v
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                scpy(i,k,2)=cpavey(i)
              end do

!$omp end do

            end do

          end if

        else if(nbc.eq.6) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              scpy(i,k,2)=max(v(i,njm1,k)*dydt,gdydt)
            end do

!$omp end do

          end do

        else if(nbc.ge.7) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              scpy(i,k,2)=gdydt
            end do

!$omp end do

          end do

        end if

        if(mfcopt.eq.1.and.(mpopt.ne.0.and.mpopt.ne.10)) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              scpy(i,k,2)=min(scpy(i,k,2),rmf(i,njm2,2))
            end do

!$omp end do

          end do

        else

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              scpy(i,k,2)=min(scpy(i,k,2),1.e0)
            end do

!$omp end do

          end do

        end if

        if(fproc(1:3).eq.'sml') then

          do k=2,nk-2

!$omp do schedule(runtime) private(i)

            do i=1,ni-1
              scpy(i,k,2)=scpy(i,k,2)*dtdvb
            end do

!$omp end do

          end do

        else

          if(advopt.ge.4) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i)

              do i=1,ni-1
                scpy(i,k,2)=scpy(i,k,2)*dtdvb
              end do

!$omp end do

            end do

          end if

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_phvbcs

!-----7--------------------------------------------------------------7--

      end module m_phvbcs
