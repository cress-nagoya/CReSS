!***********************************************************************
      module m_stretch
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/06/07, 1999/06/21,
!                   1999/11/01, 1999/11/19, 1999/11/24, 2000/01/17,
!                   2000/07/05, 2001/04/15, 2001/05/29, 2001/06/06,
!                   2001/06/29, 2001/10/18, 2002/04/02, 2002/06/18,
!                   2003/01/04, 2003/04/30, 2003/05/19, 2003/10/31,
!                   2004/09/01, 2004/09/10, 2005/02/10, 2005/08/05,
!                   2007/01/20, 2007/01/31, 2007/10/19, 2008/05/02,
!                   2008/08/19, 2008/08/25, 2009/01/05, 2009/02/27,
!                   2009/11/13, 2012/01/25, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the 1 dimensional physical coordinates.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commath
      use m_commpi
      use m_cpondpe
      use m_destroy
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: stretch, s_stretch

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface stretch

        module procedure s_stretch

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic int
      intrinsic max
      intrinsic min
      intrinsic real
      intrinsic sqrt
      intrinsic tanh

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_stretch(fpsthopt,fpzsfc,fpdzmin,fplayer1,fplayer2,   &
     &                     ztop1,ztop2,nk,zsth,dzsth)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsthopt
                       ! Formal parameter of unique index of sthopt

      integer, intent(in) :: fpzsfc
                       ! Formal parameter of unique index of zsfc

      integer, intent(in) :: fpdzmin
                       ! Formal parameter of unique index of dzmin

      integer, intent(in) :: fplayer1
                       ! Formal parameter of unique index of layer1

      integer, intent(in) :: fplayer2
                       ! Formal parameter of unique index of layer2

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: ztop1
                       ! zeta coordinate at nk - 2

      real, intent(in) :: ztop2
                       ! zeta coordinate at nk - 1

! Output variables

      real, intent(out) :: zsth(1:nk)
                       ! 1 dimensional stretched z coordinates

! Internal shared variables

      integer sthopt   ! Option for vertical grid stretching

      integer stat     ! Runtime status

      integer klow     ! Maximum k index for lower levels
      integer kmid     ! Maximum k index for middle levels

      real zsfc        ! Sea surface terrain height

      real dzmin       ! Minimum dz in lowest layer

      real layer1      ! Lowest stretching level
      real layer2      ! Highest stretching level

      real z1          ! zsfc + layer1
      real z2          ! z1 + layer2

      real dzall       ! Distance of unstretched z coordinates

      real dzmid       ! Distance of
                       ! stretched z coordinates for middle levels

      real dzhigh      ! Distance of
                       ! stretched z coordinates for higher levels

      real, intent(inout) :: dzsth(1:nk)
                       ! Distance of
                       ! 1 dimensional stretched z coordinates

      real a           ! Temporary variable
      real b           ! Temporary variable
      real c           ! Temporary variable

! Internal private variables

      integer k        ! Array index in z direction

      real tmp         ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsthopt,sthopt)
      call getrname(fpzsfc,zsfc)
      call getrname(fpdzmin,dzmin)
      call getrname(fplayer1,layer1)
      call getrname(fplayer2,layer2)

! -----

! Get the 3 levels physical height.

      layer1=max(min(ztop1-zsfc,layer1),0.e0)

      z1=zsfc+layer1

      layer2=max(min(ztop2-z1,layer2),0.e0)

      z2=z1+layer2

! -----

!! Calculate the z physical coordinates in the case vertical stretching
!! is not applied.

      if(int((z1+eps)/(ztop2+eps)+.1e0).eq.1                            &
     &  .and.int((z1-zsfc+eps)/(real(nk-3)*dzmin+eps)+.1e0).eq.1) then

! Set the common used variable, dzall.

        dzall=(ztop2-zsfc)/real(nk-3)

! -----

! Get the constant stretching function.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(k)

        do k=1,nk
          zsth(k)=zsfc+real(k-2)*dzall
        end do

!$omp end do

!$omp end parallel

! -----

!! -----

!!! Calculate the z physical coordinates in the case vertical stretching
!!! is applied.

      else

! Get the 3 levels dz and k index.

        z1=max(zsfc,z1)
        z2=min(ztop2,z2)

        klow=int((z1-zsfc)/dzmin+.1e0)

        stat=0

        if((z1-zsfc).ge.real(nk-3)*dzmin) then
          stat=1
        end if

        call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('stretch ',7,'cont',7,'              ',14,101, &
     &                   stat)

          end if

          call cpondpe

          call destroy('stretch ',7,'stop',1001,'              ',14,101,&
     &                 stat)

        end if

        if(z2.ge.ztop2) then

          dzmid=(ztop2-zsfc-real(klow)*dzmin)/real(nk-3-klow)
          dzhigh=2.e0*dzmid-dzmin

          kmid=0

        else if(z2.lt.ztop2) then

          a=2.e0*real(nk-3-klow)
          b=2.e0*zsfc-ztop2-z2-real(nk-3-3*klow)*dzmin
          c=dzmin*(z2-zsfc-real(klow)*dzmin)

          dzmid=.5e0*(sqrt(b*b-4.e0*a*c)-b)/a

          kmid=int((ztop2-z2)/(2.e0*dzmid-dzmin)+.1e0)

          a=zsfc+real(klow-kmid)*dzmin+real(nk-3-klow+kmid)*dzmid

          if(kmid.ne.0) then

            dzhigh=(2.e0*dzmid-dzmin)+(ztop2-a)/real(kmid)

          else if(nk-3-klow-kmid.ne.0) then

            dzmid=dzmid+(ztop2-a)/real(nk-3-klow-kmid)
            dzhigh=2.e0*dzmid-dzmin

          end if

        end if

! -----

! Set the constants for stretching function.

        if(sthopt.eq.1) then

          a=dzmid-dzmin

          b=1.e0/real(nk-4-kmid-klow)

        else

          a=.5e0*real(nk-kmid+klow)

          if(2*(klow+2)-(nk-kmid+klow).eq.0) then
            b=0.e0
          else
            b=2.e0/(real(klow+2)-a)
          end if

          c=(dzmin-dzmid)/tanh2

        end if

! -----

! Set the bottom boundary condition.

        zsth(2)=zsfc

! -----

!! Apply the stretching function.

!$omp parallel default(shared)

! Calculate the dz of stretching for the low level.

!$omp do schedule(runtime) private(k)

        do k=1,klow+1
          dzsth(k)=dzmin
        end do

!$omp end do

! -----

! Calculate the dz of stretching with the cubic function.

        if(sthopt.eq.1) then

!$omp do schedule(runtime) private(k,tmp)

          do k=klow+2,nk-2-kmid
            tmp=2.e0*b*real(k-2-klow)-1.e0

            dzsth(k)=dzmid+a*tmp*tmp*tmp

          end do

!$omp end do

! -----

! Calculate the dz of stretching with the hyperbolic tangent function.

        else if(sthopt.eq.2) then

!$omp do schedule(runtime) private(k)

          do k=klow+2,nk-2-kmid
            dzsth(k)=dzmid+c*tanh(b*(real(k)-a))
          end do

!$omp end do

        end if

! -----

! Calculate the dz of stretching for the higher levels.

!$omp do schedule(runtime) private(k)

        do k=nk-1-kmid,nk-2
          dzsth(k)=dzhigh
        end do

!$omp end do

! -----

! Finally get the z physical coordinates.

!$omp single private(k)

        do k=3,nk-2
          zsth(k)=zsth(k-1)+dzsth(k-1)
        end do

!$omp end single

! -----

!$omp end parallel

!! -----

! Set the bottom and top boundary conditions.

        zsth(1)=zsth(2)-dzsth(1)

        zsth(nk-1)=ztop2
        zsth(nk)=zsth(nk-1)+dzsth(nk-2)

! -----

      end if

!!! -----

      end subroutine s_stretch

!-----7--------------------------------------------------------------7--

      end module m_stretch
