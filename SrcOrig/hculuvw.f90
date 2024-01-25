!***********************************************************************
      module m_hculuvw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/04/03
!     Modification: 2006/05/12, 2006/06/21, 2006/11/06, 2007/01/31,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2011/07/15, 2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the horizontal velocity advection by the Cubic Lagrange
!     scheme.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_commath
      use m_getiname
      use m_getrname
      use m_lbculu
      use m_lbculv
      use m_lbculw

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: hculuvw, s_hculuvw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface hculuvw

        module procedure s_hculuvw

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
      subroutine s_hculuvw(fpmpopt,fpmfcopt,fpadvopt,fpdxiv,fpdyiv,dtb, &
     &                     ni,nj,nk,mf8u,mf8v,up,vp,wp,uf,vf,wf,advd,   &
     &                     dxt22,dxt16,dxt36w,dxt36e,                   &
     &                     dyt22,dyt16,dyt36s,dyt36n)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpdxiv
                       ! Formal parameter of unique index of dxiv

      integer, intent(in) :: fpdyiv
                       ! Formal parameter of unique index of dyiv

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: mf8u(0:ni+1,0:nj+1)
                       ! Map scale factors at u points

      real, intent(in) :: mf8v(0:ni+1,0:nj+1)
                       ! Map scale factors at v points

! Input and output variables

      real, intent(inout) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(inout) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

      real, intent(inout) :: wp(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at past

! Output variables

      real, intent(out) :: uf(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at future

      real, intent(out) :: vf(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at future

      real, intent(out) :: wf(0:ni+1,0:nj+1,1:nk)
                       ! z components of velocity at future

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor
      integer advopt   ! Option for advection scheme

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy

      real cdxt2       ! dxiv^2 x dtb^2 / 2.0
      real cdxt3       ! - dxiv^3 x dtb^3 / 6.0
      real cdxt1       ! - dxiv x dtb / 6.0

      real cdyt2       ! dyiv^2 x dtb^2 / 2.0
      real cdyt3       ! - dyiv^3 x dtb^3 / 6.0
      real cdyt1       ! - dyiv x dtb / 6.0

      real, intent(inout) :: advd(0:ni+1,0:nj+1,1:nk)
                       ! Advected value in x or y direction

      real, intent(inout) :: dxt22(0:ni+1)
                       ! Array of cdxt2

      real, intent(inout) :: dxt16(0:ni+1)
                       ! Array of cdxt1

      real, intent(inout) :: dxt36w(0:ni+1)
                       ! Array of cdxt3 and 0 on west boundary

      real, intent(inout) :: dxt36e(0:ni+1)
                       ! Array of cdxt3 and 0 on east boundary

      real, intent(inout) :: dyt22(0:nj+1)
                       ! Array of cdyt2

      real, intent(inout) :: dyt16(0:nj+1)
                       ! Array of cdyt1

      real, intent(inout) :: dyt36s(0:nj+1)
                       ! Array of cdyt3 and 0 on south boundary

      real, intent(inout) :: dyt36n(0:nj+1)
                       ! Array of cdyt3 and 0 on north boundary

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real um          ! mf8u x u
      real vm          ! mf8v x v

      real u8v         ! u at v points
      real u8w         ! u at w points

      real v8u         ! v at u points
      real v8w         ! v at w points

      real advij       ! Temporary variable

      real advim2      ! Temporary variable
      real advim1      ! Temporary variable

      real advip2      ! Temporary variable
      real advip1      ! Temporary variable

      real advjm2      ! Temporary variable
      real advjm1      ! Temporary variable

      real advjp2      ! Temporary variable
      real advjp1      ! Temporary variable

      real a           ! Temporary variable
      real b           ! Temporary variable
      real c           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getiname(fpadvopt,advopt)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)

! -----

! Set the common used variables.

      cdxt2=.5e0*dxiv*dxiv*dtb*dtb
      cdxt3=-oned6*dxiv*dxiv*dxiv*dtb*dtb*dtb
      cdxt1=-oned6*dxiv*dtb

      cdyt2=.5e0*dyiv*dyiv*dtb*dtb
      cdyt3=-oned6*dyiv*dyiv*dyiv*dtb*dtb*dtb
      cdyt1=-oned6*dyiv*dtb

! -----

!!!!! Calculate the u advection horizontally.

! Set the common used variables.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i)

      do i=0,ni+1
        dxt22(i)=cdxt2
        dxt16(i)=cdxt1

        dxt36w(i)=cdxt3
        dxt36e(i)=cdxt3

      end do

!$omp end do

!$omp do schedule(runtime) private(j)

      do j=0,nj
        dyt22(j)=cdyt2
        dyt16(j)=cdyt1

        dyt36s(j)=cdyt3
        dyt36n(j)=cdyt3

      end do

!$omp end do

!$omp end parallel

! -----

! Set the lateral boundary conditions.

      call lbculu(idwbc,idebc,idsbc,idnbc,ni,nj,nk,up,                  &
     &            dxt36w,dxt36e,dyt36s,dyt36n)

! -----

!!!! Perform Cubic Lagrange scheme.

!$omp parallel default(shared) private(k)

!!! Perform horizontal-vertical seperated Cubic Lagrange scheme.

      if(advopt.eq.4) then

! Perform calculation without the map scale factor.

        if(mfcopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,v8u,advij,advjm2,advjm1,advjp1,advjp2,a,b,c)

            do j=2,nj-2
            do i=2,ni-1

              v8u=.25e0                                                 &
     &          *((vp(i-1,j,k)+vp(i,j,k))+(vp(i-1,j+1,k)+vp(i,j+1,k)))

              if(up(i,j,k).gt.0.e0) then

                if(v8u.gt.0.e0) then

                  a=dxt36w(i)*(up(i+1,j-2,k)-3.e0*up(i,j-2,k)           &
     &              +3.e0*up(i-1,j-2,k)-up(i-2,j-2,k))

                  b=dxt22(i)*(up(i+1,j-2,k)                             &
     &              -2.e0*up(i,j-2,k)+up(i-1,j-2,k))

                  c=dxt16(i)*(3.e0*up(i,j-2,k)+2.e0*up(i+1,j-2,k)       &
     &              +up(i-2,j-2,k)-6.e0*up(i-1,j-2,k))

                  advjm2=up(i,j-2,k)                                    &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dxt36w(i)*(up(i+1,j-1,k)-3.e0*up(i,j-1,k)           &
     &              +3.e0*up(i-1,j-1,k)-up(i-2,j-1,k))

                  b=dxt22(i)*(up(i+1,j-1,k)                             &
     &              -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                  c=dxt16(i)*(3.e0*up(i,j-1,k)+2.e0*up(i+1,j-1,k)       &
     &              +up(i-2,j-1,k)-6.e0*up(i-1,j-1,k))

                  advjm1=up(i,j-1,k)                                    &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dxt36w(i)*(up(i+1,j,k)-3.e0*up(i,j,k)               &
     &              +3.e0*up(i-1,j,k)-up(i-2,j,k))

                  b=dxt22(i)*(up(i+1,j,k)                               &
     &              -2.e0*up(i,j,k)+up(i-1,j,k))

                  c=dxt16(i)*(3.e0*up(i,j,k)+2.e0*up(i+1,j,k)           &
     &              +up(i-2,j,k)-6.e0*up(i-1,j,k))

                  advij=up(i,j,k)                                       &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dxt36w(i)*(up(i+1,j+1,k)-3.e0*up(i,j+1,k)           &
     &              +3.e0*up(i-1,j+1,k)-up(i-2,j+1,k))

                  b=dxt22(i)*(up(i+1,j+1,k)                             &
     &              -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                  c=dxt16(i)*(3.e0*up(i,j+1,k)+2.e0*up(i+1,j+1,k)       &
     &              +up(i-2,j+1,k)-6.e0*up(i-1,j+1,k))

                  advjp1=up(i,j+1,k)                                    &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                  b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                  c=dyt16(j)                                            &
     &              *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                  uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                else

                  a=dxt36w(i)*(up(i+1,j-1,k)-3.e0*up(i,j-1,k)           &
     &              +3.e0*up(i-1,j-1,k)-up(i-2,j-1,k))

                  b=dxt22(i)*(up(i+1,j-1,k)                             &
     &              -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                  c=dxt16(i)*(3.e0*up(i,j-1,k)+2.e0*up(i+1,j-1,k)       &
     &              +up(i-2,j-1,k)-6.e0*up(i-1,j-1,k))

                  advjm1=up(i,j-1,k)                                    &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dxt36w(i)*(up(i+1,j,k)-3.e0*up(i,j,k)               &
     &              +3.e0*up(i-1,j,k)-up(i-2,j,k))

                  b=dxt22(i)*(up(i+1,j,k)                               &
     &              -2.e0*up(i,j,k)+up(i-1,j,k))

                  c=dxt16(i)*(3.e0*up(i,j,k)+2.e0*up(i+1,j,k)           &
     &              +up(i-2,j,k)-6.e0*up(i-1,j,k))

                  advij=up(i,j,k)                                       &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dxt36w(i)*(up(i+1,j+1,k)-3.e0*up(i,j+1,k)           &
     &              +3.e0*up(i-1,j+1,k)-up(i-2,j+1,k))

                  b=dxt22(i)*(up(i+1,j+1,k)                             &
     &              -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                  c=dxt16(i)*(3.e0*up(i,j+1,k)+2.e0*up(i+1,j+1,k)       &
     &              +up(i-2,j+1,k)-6.e0*up(i-1,j+1,k))

                  advjp1=up(i,j+1,k)                                    &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dxt36w(i)*(up(i+1,j+2,k)-3.e0*up(i,j+2,k)           &
     &              +3.e0*up(i-1,j+2,k)-up(i-2,j+2,k))

                  b=dxt22(i)*(up(i+1,j+2,k)                             &
     &              -2.e0*up(i,j+2,k)+up(i-1,j+2,k))

                  c=dxt16(i)*(3.e0*up(i,j+2,k)+2.e0*up(i+1,j+2,k)       &
     &              +up(i-2,j+2,k)-6.e0*up(i-1,j+2,k))

                  advjp2=up(i,j+2,k)                                    &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                  b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                  c=dyt16(j)                                            &
     &              *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                  uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                end if

              else

                if(v8u.gt.0.e0) then

                  a=dxt36e(i)*(up(i+2,j-2,k)-3.e0*up(i+1,j-2,k)         &
     &              +3.e0*up(i,j-2,k)-up(i-1,j-2,k))

                  b=dxt22(i)*(up(i+1,j-2,k)                             &
     &              -2.e0*up(i,j-2,k)+up(i-1,j-2,k))

                  c=dxt16(i)*(6.e0*up(i+1,j-2,k)-up(i+2,j-2,k)          &
     &              -2.e0*up(i-1,j-2,k)-3.e0*up(i,j-2,k))

                  advjm2=up(i,j-2,k)                                    &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dxt36e(i)*(up(i+2,j-1,k)-3.e0*up(i+1,j-1,k)         &
     &              +3.e0*up(i,j-1,k)-up(i-1,j-1,k))

                  b=dxt22(i)*(up(i+1,j-1,k)                             &
     &              -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                  c=dxt16(i)*(6.e0*up(i+1,j-1,k)-up(i+2,j-1,k)          &
     &              -2.e0*up(i-1,j-1,k)-3.e0*up(i,j-1,k))

                  advjm1=up(i,j-1,k)                                    &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dxt36e(i)*(up(i+2,j,k)-3.e0*up(i+1,j,k)             &
     &              +3.e0*up(i,j,k)-up(i-1,j,k))

                  b=dxt22(i)*(up(i+1,j,k)                               &
     &              -2.e0*up(i,j,k)+up(i-1,j,k))

                  c=dxt16(i)*(6.e0*up(i+1,j,k)-up(i+2,j,k)              &
     &              -2.e0*up(i-1,j,k)-3.e0*up(i,j,k))

                  advij=up(i,j,k)                                       &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dxt36e(i)*(up(i+2,j+1,k)-3.e0*up(i+1,j+1,k)         &
     &              +3.e0*up(i,j+1,k)-up(i-1,j+1,k))

                  b=dxt22(i)*(up(i+1,j+1,k)                             &
     &              -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                  c=dxt16(i)*(6.e0*up(i+1,j+1,k)-up(i+2,j+1,k)          &
     &              -2.e0*up(i-1,j+1,k)-3.e0*up(i,j+1,k))

                  advjp1=up(i,j+1,k)                                    &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                  b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                  c=dyt16(j)                                            &
     &              *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                  uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                else

                  a=dxt36e(i)*(up(i+2,j-1,k)-3.e0*up(i+1,j-1,k)         &
     &              +3.e0*up(i,j-1,k)-up(i-1,j-1,k))

                  b=dxt22(i)*(up(i+1,j-1,k)                             &
     &              -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                  c=dxt16(i)*(6.e0*up(i+1,j-1,k)-up(i+2,j-1,k)          &
     &              -2.e0*up(i-1,j-1,k)-3.e0*up(i,j-1,k))

                  advjm1=up(i,j-1,k)                                    &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dxt36e(i)*(up(i+2,j,k)-3.e0*up(i+1,j,k)             &
     &              +3.e0*up(i,j,k)-up(i-1,j,k))

                  b=dxt22(i)*(up(i+1,j,k)                               &
     &              -2.e0*up(i,j,k)+up(i-1,j,k))

                  c=dxt16(i)*(6.e0*up(i+1,j,k)-up(i+2,j,k)              &
     &              -2.e0*up(i-1,j,k)-3.e0*up(i,j,k))

                  advij=up(i,j,k)                                       &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dxt36e(i)*(up(i+2,j+1,k)-3.e0*up(i+1,j+1,k)         &
     &              +3.e0*up(i,j+1,k)-up(i-1,j+1,k))

                  b=dxt22(i)*(up(i+1,j+1,k)                             &
     &              -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                  c=dxt16(i)*(6.e0*up(i+1,j+1,k)-up(i+2,j+1,k)          &
     &              -2.e0*up(i-1,j+1,k)-3.e0*up(i,j+1,k))

                  advjp1=up(i,j+1,k)                                    &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dxt36e(i)*(up(i+2,j+2,k)-3.e0*up(i+1,j+2,k)         &
     &              +3.e0*up(i,j+2,k)-up(i-1,j+2,k))

                  b=dxt22(i)*(up(i+1,j+2,k)                             &
     &              -2.e0*up(i,j+2,k)+up(i-1,j+2,k))

                  c=dxt16(i)*(6.e0*up(i+1,j+2,k)-up(i+2,j+2,k)          &
     &              -2.e0*up(i-1,j+2,k)-3.e0*up(i,j+2,k))

                  advjp2=up(i,j+2,k)                                    &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                  a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                  b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                  c=dyt16(j)                                            &
     &              *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                  uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                end if

              end if

            end do
            end do

!$omp end do

          end do

! -----

!! Perform calculation with the map scale factor.

        else

! Applied the map scale factor for the x direction.

          if(mpopt.eq.0.or.mpopt.eq.10) then

            do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,um,v8u,advij,advjm2,advjm1,advjp1,advjp2,a,b,c)

              do j=2,nj-2
              do i=2,ni-1

                um=mf8u(i,j)*up(i,j,k)

                v8u=.25e0                                               &
     &            *((vp(i-1,j,k)+vp(i,j,k))+(vp(i-1,j+1,k)+vp(i,j+1,k)))

                if(um.gt.0.e0) then

                  if(v8u.gt.0.e0) then

                    a=dxt36w(i)*(up(i+1,j-2,k)-3.e0*up(i,j-2,k)         &
     &                +3.e0*up(i-1,j-2,k)-up(i-2,j-2,k))

                    b=dxt22(i)*(up(i+1,j-2,k)                           &
     &                -2.e0*up(i,j-2,k)+up(i-1,j-2,k))

                    c=dxt16(i)*(3.e0*up(i,j-2,k)+2.e0*up(i+1,j-2,k)     &
     &                +up(i-2,j-2,k)-6.e0*up(i-1,j-2,k))

                    advjm2=up(i,j-2,k)+((a*um+b)*um+c)*um

                    a=dxt36w(i)*(up(i+1,j-1,k)-3.e0*up(i,j-1,k)         &
     &                +3.e0*up(i-1,j-1,k)-up(i-2,j-1,k))

                    b=dxt22(i)*(up(i+1,j-1,k)                           &
     &                -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*up(i,j-1,k)+2.e0*up(i+1,j-1,k)     &
     &                +up(i-2,j-1,k)-6.e0*up(i-1,j-1,k))

                    advjm1=up(i,j-1,k)+((a*um+b)*um+c)*um

                    a=dxt36w(i)*(up(i+1,j,k)-3.e0*up(i,j,k)             &
     &                +3.e0*up(i-1,j,k)-up(i-2,j,k))

                    b=dxt22(i)*(up(i+1,j,k)                             &
     &                -2.e0*up(i,j,k)+up(i-1,j,k))

                    c=dxt16(i)*(3.e0*up(i,j,k)+2.e0*up(i+1,j,k)         &
     &                +up(i-2,j,k)-6.e0*up(i-1,j,k))

                    advij=up(i,j,k)+((a*um+b)*um+c)*um

                    a=dxt36w(i)*(up(i+1,j+1,k)-3.e0*up(i,j+1,k)         &
     &                +3.e0*up(i-1,j+1,k)-up(i-2,j+1,k))

                    b=dxt22(i)*(up(i+1,j+1,k)                           &
     &                -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*up(i,j+1,k)+2.e0*up(i+1,j+1,k)     &
     &                +up(i-2,j+1,k)-6.e0*up(i-1,j+1,k))

                    advjp1=up(i,j+1,k)+((a*um+b)*um+c)*um

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                  else

                    a=dxt36w(i)*(up(i+1,j-1,k)-3.e0*up(i,j-1,k)         &
     &                +3.e0*up(i-1,j-1,k)-up(i-2,j-1,k))

                    b=dxt22(i)*(up(i+1,j-1,k)                           &
     &                -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*up(i,j-1,k)+2.e0*up(i+1,j-1,k)     &
     &                +up(i-2,j-1,k)-6.e0*up(i-1,j-1,k))

                    advjm1=up(i,j-1,k)+((a*um+b)*um+c)*um

                    a=dxt36w(i)*(up(i+1,j,k)-3.e0*up(i,j,k)             &
     &                +3.e0*up(i-1,j,k)-up(i-2,j,k))

                    b=dxt22(i)*(up(i+1,j,k)                             &
     &                -2.e0*up(i,j,k)+up(i-1,j,k))

                    c=dxt16(i)*(3.e0*up(i,j,k)+2.e0*up(i+1,j,k)         &
     &                +up(i-2,j,k)-6.e0*up(i-1,j,k))

                    advij=up(i,j,k)+((a*um+b)*um+c)*um

                    a=dxt36w(i)*(up(i+1,j+1,k)-3.e0*up(i,j+1,k)         &
     &                +3.e0*up(i-1,j+1,k)-up(i-2,j+1,k))

                    b=dxt22(i)*(up(i+1,j+1,k)                           &
     &                -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*up(i,j+1,k)+2.e0*up(i+1,j+1,k)     &
     &                +up(i-2,j+1,k)-6.e0*up(i-1,j+1,k))

                    advjp1=up(i,j+1,k)+((a*um+b)*um+c)*um

                    a=dxt36w(i)*(up(i+1,j+2,k)-3.e0*up(i,j+2,k)         &
     &                +3.e0*up(i-1,j+2,k)-up(i-2,j+2,k))

                    b=dxt22(i)*(up(i+1,j+2,k)                           &
     &                -2.e0*up(i,j+2,k)+up(i-1,j+2,k))

                    c=dxt16(i)*(3.e0*up(i,j+2,k)+2.e0*up(i+1,j+2,k)     &
     &                +up(i-2,j+2,k)-6.e0*up(i-1,j+2,k))

                    advjp2=up(i,j+2,k)+((a*um+b)*um+c)*um

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                  end if

                else

                  if(v8u.gt.0.e0) then

                    a=dxt36e(i)*(up(i+2,j-2,k)-3.e0*up(i+1,j-2,k)       &
     &                +3.e0*up(i,j-2,k)-up(i-1,j-2,k))

                    b=dxt22(i)*(up(i+1,j-2,k)                           &
     &                -2.e0*up(i,j-2,k)+up(i-1,j-2,k))

                    c=dxt16(i)*(6.e0*up(i+1,j-2,k)-up(i+2,j-2,k)        &
     &                -2.e0*up(i-1,j-2,k)-3.e0*up(i,j-2,k))

                    advjm2=up(i,j-2,k)+((a*um+b)*um+c)*um

                    a=dxt36e(i)*(up(i+2,j-1,k)-3.e0*up(i+1,j-1,k)       &
     &                +3.e0*up(i,j-1,k)-up(i-1,j-1,k))

                    b=dxt22(i)*(up(i+1,j-1,k)                           &
     &                -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*up(i+1,j-1,k)-up(i+2,j-1,k)        &
     &                -2.e0*up(i-1,j-1,k)-3.e0*up(i,j-1,k))

                    advjm1=up(i,j-1,k)+((a*um+b)*um+c)*um

                    a=dxt36e(i)*(up(i+2,j,k)-3.e0*up(i+1,j,k)           &
     &                +3.e0*up(i,j,k)-up(i-1,j,k))

                    b=dxt22(i)*(up(i+1,j,k)                             &
     &                -2.e0*up(i,j,k)+up(i-1,j,k))

                    c=dxt16(i)*(6.e0*up(i+1,j,k)-up(i+2,j,k)            &
     &                -2.e0*up(i-1,j,k)-3.e0*up(i,j,k))

                    advij=up(i,j,k)+((a*um+b)*um+c)*um

                    a=dxt36e(i)*(up(i+2,j+1,k)-3.e0*up(i+1,j+1,k)       &
     &                +3.e0*up(i,j+1,k)-up(i-1,j+1,k))

                    b=dxt22(i)*(up(i+1,j+1,k)                           &
     &                -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*up(i+1,j+1,k)-up(i+2,j+1,k)        &
     &                -2.e0*up(i-1,j+1,k)-3.e0*up(i,j+1,k))

                    advjp1=up(i,j+1,k)+((a*um+b)*um+c)*um

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                  else

                    a=dxt36e(i)*(up(i+2,j-1,k)-3.e0*up(i+1,j-1,k)       &
     &                +3.e0*up(i,j-1,k)-up(i-1,j-1,k))

                    b=dxt22(i)*(up(i+1,j-1,k)                           &
     &                -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*up(i+1,j-1,k)-up(i+2,j-1,k)        &
     &                -2.e0*up(i-1,j-1,k)-3.e0*up(i,j-1,k))

                    advjm1=up(i,j-1,k)+((a*um+b)*um+c)*um

                    a=dxt36e(i)*(up(i+2,j,k)-3.e0*up(i+1,j,k)           &
     &                +3.e0*up(i,j,k)-up(i-1,j,k))

                    b=dxt22(i)*(up(i+1,j,k)                             &
     &                -2.e0*up(i,j,k)+up(i-1,j,k))

                    c=dxt16(i)*(6.e0*up(i+1,j,k)-up(i+2,j,k)            &
     &                -2.e0*up(i-1,j,k)-3.e0*up(i,j,k))

                    advij=up(i,j,k)+((a*um+b)*um+c)*um

                    a=dxt36e(i)*(up(i+2,j+1,k)-3.e0*up(i+1,j+1,k)       &
     &                +3.e0*up(i,j+1,k)-up(i-1,j+1,k))

                    b=dxt22(i)*(up(i+1,j+1,k)                           &
     &                -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*up(i+1,j+1,k)-up(i+2,j+1,k)        &
     &                -2.e0*up(i-1,j+1,k)-3.e0*up(i,j+1,k))

                    advjp1=up(i,j+1,k)+((a*um+b)*um+c)*um

                    a=dxt36e(i)*(up(i+2,j+2,k)-3.e0*up(i+1,j+2,k)       &
     &                +3.e0*up(i,j+2,k)-up(i-1,j+2,k))

                    b=dxt22(i)*(up(i+1,j+2,k)                           &
     &                -2.e0*up(i,j+2,k)+up(i-1,j+2,k))

                    c=dxt16(i)*(6.e0*up(i+1,j+2,k)-up(i+2,j+2,k)        &
     &                -2.e0*up(i-1,j+2,k)-3.e0*up(i,j+2,k))

                    advjp2=up(i,j+2,k)+((a*um+b)*um+c)*um

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                  end if

                end if

              end do
              end do

!$omp end do

            end do

! -----

! Applied the map scale factor for the y direction.

          else if(mpopt.eq.5) then

            do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,v8u,advij,advjm2,advjm1,advjp1,advjp2,a,b,c)

              do j=2,nj-2
              do i=2,ni-1

                v8u=.25e0*((mf8v(i-1,j)*vp(i-1,j,k)+mf8v(i,j)*vp(i,j,k))&
     &           +(mf8v(i-1,j+1)*vp(i-1,j+1,k)+mf8v(i,j+1)*vp(i,j+1,k)))

                if(up(i,j,k).gt.0.e0) then

                  if(v8u.gt.0.e0) then

                    a=dxt36w(i)*(up(i+1,j-2,k)-3.e0*up(i,j-2,k)         &
     &                +3.e0*up(i-1,j-2,k)-up(i-2,j-2,k))

                    b=dxt22(i)*(up(i+1,j-2,k)                           &
     &                -2.e0*up(i,j-2,k)+up(i-1,j-2,k))

                    c=dxt16(i)*(3.e0*up(i,j-2,k)+2.e0*up(i+1,j-2,k)     &
     &                +up(i-2,j-2,k)-6.e0*up(i-1,j-2,k))

                    advjm2=up(i,j-2,k)                                  &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dxt36w(i)*(up(i+1,j-1,k)-3.e0*up(i,j-1,k)         &
     &                +3.e0*up(i-1,j-1,k)-up(i-2,j-1,k))

                    b=dxt22(i)*(up(i+1,j-1,k)                           &
     &                -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*up(i,j-1,k)+2.e0*up(i+1,j-1,k)     &
     &                +up(i-2,j-1,k)-6.e0*up(i-1,j-1,k))

                    advjm1=up(i,j-1,k)                                  &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dxt36w(i)*(up(i+1,j,k)-3.e0*up(i,j,k)             &
     &                +3.e0*up(i-1,j,k)-up(i-2,j,k))

                    b=dxt22(i)*(up(i+1,j,k)                             &
     &                -2.e0*up(i,j,k)+up(i-1,j,k))

                    c=dxt16(i)*(3.e0*up(i,j,k)+2.e0*up(i+1,j,k)         &
     &                +up(i-2,j,k)-6.e0*up(i-1,j,k))

                    advij=up(i,j,k)                                     &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dxt36w(i)*(up(i+1,j+1,k)-3.e0*up(i,j+1,k)         &
     &                +3.e0*up(i-1,j+1,k)-up(i-2,j+1,k))

                    b=dxt22(i)*(up(i+1,j+1,k)                           &
     &                -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*up(i,j+1,k)+2.e0*up(i+1,j+1,k)     &
     &                +up(i-2,j+1,k)-6.e0*up(i-1,j+1,k))

                    advjp1=up(i,j+1,k)                                  &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                  else

                    a=dxt36w(i)*(up(i+1,j-1,k)-3.e0*up(i,j-1,k)         &
     &                +3.e0*up(i-1,j-1,k)-up(i-2,j-1,k))

                    b=dxt22(i)*(up(i+1,j-1,k)                           &
     &                -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*up(i,j-1,k)+2.e0*up(i+1,j-1,k)     &
     &                +up(i-2,j-1,k)-6.e0*up(i-1,j-1,k))

                    advjm1=up(i,j-1,k)                                  &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dxt36w(i)*(up(i+1,j,k)-3.e0*up(i,j,k)             &
     &                +3.e0*up(i-1,j,k)-up(i-2,j,k))

                    b=dxt22(i)*(up(i+1,j,k)                             &
     &                -2.e0*up(i,j,k)+up(i-1,j,k))

                    c=dxt16(i)*(3.e0*up(i,j,k)+2.e0*up(i+1,j,k)         &
     &                +up(i-2,j,k)-6.e0*up(i-1,j,k))

                    advij=up(i,j,k)                                     &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dxt36w(i)*(up(i+1,j+1,k)-3.e0*up(i,j+1,k)         &
     &                +3.e0*up(i-1,j+1,k)-up(i-2,j+1,k))

                    b=dxt22(i)*(up(i+1,j+1,k)                           &
     &                -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*up(i,j+1,k)+2.e0*up(i+1,j+1,k)     &
     &                +up(i-2,j+1,k)-6.e0*up(i-1,j+1,k))

                    advjp1=up(i,j+1,k)                                  &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dxt36w(i)*(up(i+1,j+2,k)-3.e0*up(i,j+2,k)         &
     &                +3.e0*up(i-1,j+2,k)-up(i-2,j+2,k))

                    b=dxt22(i)*(up(i+1,j+2,k)                           &
     &                -2.e0*up(i,j+2,k)+up(i-1,j+2,k))

                    c=dxt16(i)*(3.e0*up(i,j+2,k)+2.e0*up(i+1,j+2,k)     &
     &                +up(i-2,j+2,k)-6.e0*up(i-1,j+2,k))

                    advjp2=up(i,j+2,k)                                  &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                  end if

                else

                  if(v8u.gt.0.e0) then

                    a=dxt36e(i)*(up(i+2,j-2,k)-3.e0*up(i+1,j-2,k)       &
     &                +3.e0*up(i,j-2,k)-up(i-1,j-2,k))

                    b=dxt22(i)*(up(i+1,j-2,k)                           &
     &                -2.e0*up(i,j-2,k)+up(i-1,j-2,k))

                    c=dxt16(i)*(6.e0*up(i+1,j-2,k)-up(i+2,j-2,k)        &
     &                -2.e0*up(i-1,j-2,k)-3.e0*up(i,j-2,k))

                    advjm2=up(i,j-2,k)                                  &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dxt36e(i)*(up(i+2,j-1,k)-3.e0*up(i+1,j-1,k)       &
     &                +3.e0*up(i,j-1,k)-up(i-1,j-1,k))

                    b=dxt22(i)*(up(i+1,j-1,k)                           &
     &                -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*up(i+1,j-1,k)-up(i+2,j-1,k)        &
     &                -2.e0*up(i-1,j-1,k)-3.e0*up(i,j-1,k))

                    advjm1=up(i,j-1,k)                                  &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dxt36e(i)*(up(i+2,j,k)-3.e0*up(i+1,j,k)           &
     &                +3.e0*up(i,j,k)-up(i-1,j,k))

                    b=dxt22(i)*(up(i+1,j,k)                             &
     &                -2.e0*up(i,j,k)+up(i-1,j,k))

                    c=dxt16(i)*(6.e0*up(i+1,j,k)-up(i+2,j,k)            &
     &                -2.e0*up(i-1,j,k)-3.e0*up(i,j,k))

                    advij=up(i,j,k)                                     &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dxt36e(i)*(up(i+2,j+1,k)-3.e0*up(i+1,j+1,k)       &
     &                +3.e0*up(i,j+1,k)-up(i-1,j+1,k))

                    b=dxt22(i)*(up(i+1,j+1,k)                           &
     &                -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*up(i+1,j+1,k)-up(i+2,j+1,k)        &
     &                -2.e0*up(i-1,j+1,k)-3.e0*up(i,j+1,k))

                    advjp1=up(i,j+1,k)                                  &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                  else

                    a=dxt36e(i)*(up(i+2,j-1,k)-3.e0*up(i+1,j-1,k)       &
     &                +3.e0*up(i,j-1,k)-up(i-1,j-1,k))

                    b=dxt22(i)*(up(i+1,j-1,k)                           &
     &                -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*up(i+1,j-1,k)-up(i+2,j-1,k)        &
     &                -2.e0*up(i-1,j-1,k)-3.e0*up(i,j-1,k))

                    advjm1=up(i,j-1,k)                                  &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dxt36e(i)*(up(i+2,j,k)-3.e0*up(i+1,j,k)           &
     &                +3.e0*up(i,j,k)-up(i-1,j,k))

                    b=dxt22(i)*(up(i+1,j,k)                             &
     &                -2.e0*up(i,j,k)+up(i-1,j,k))

                    c=dxt16(i)*(6.e0*up(i+1,j,k)-up(i+2,j,k)            &
     &                -2.e0*up(i-1,j,k)-3.e0*up(i,j,k))

                    advij=up(i,j,k)                                     &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dxt36e(i)*(up(i+2,j+1,k)-3.e0*up(i+1,j+1,k)       &
     &                +3.e0*up(i,j+1,k)-up(i-1,j+1,k))

                    b=dxt22(i)*(up(i+1,j+1,k)                           &
     &                -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*up(i+1,j+1,k)-up(i+2,j+1,k)        &
     &                -2.e0*up(i-1,j+1,k)-3.e0*up(i,j+1,k))

                    advjp1=up(i,j+1,k)                                  &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dxt36e(i)*(up(i+2,j+2,k)-3.e0*up(i+1,j+2,k)       &
     &                +3.e0*up(i,j+2,k)-up(i-1,j+2,k))

                    b=dxt22(i)*(up(i+1,j+2,k)                           &
     &                -2.e0*up(i,j+2,k)+up(i-1,j+2,k))

                    c=dxt16(i)*(6.e0*up(i+1,j+2,k)-up(i+2,j+2,k)        &
     &                -2.e0*up(i-1,j+2,k)-3.e0*up(i,j+2,k))

                    advjp2=up(i,j+2,k)                                  &
     &                +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                  end if

                end if

              end do
              end do

!$omp end do

            end do

! -----

! Applied the map scale factor for the x and y direction.

          else

            do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,um,v8u,advij,advjm2,advjm1,advjp1,advjp2,a,b,c)

              do j=2,nj-2
              do i=2,ni-1

                um=mf8u(i,j)*up(i,j,k)

                v8u=.25e0*((mf8v(i-1,j)*vp(i-1,j,k)+mf8v(i,j)*vp(i,j,k))&
     &           +(mf8v(i-1,j+1)*vp(i-1,j+1,k)+mf8v(i,j+1)*vp(i,j+1,k)))

                if(um.gt.0.e0) then

                  if(v8u.gt.0.e0) then

                    a=dxt36w(i)*(up(i+1,j-2,k)-3.e0*up(i,j-2,k)         &
     &                +3.e0*up(i-1,j-2,k)-up(i-2,j-2,k))

                    b=dxt22(i)*(up(i+1,j-2,k)                           &
     &                -2.e0*up(i,j-2,k)+up(i-1,j-2,k))

                    c=dxt16(i)*(3.e0*up(i,j-2,k)+2.e0*up(i+1,j-2,k)     &
     &                +up(i-2,j-2,k)-6.e0*up(i-1,j-2,k))

                    advjm2=up(i,j-2,k)+((a*um+b)*um+c)*um

                    a=dxt36w(i)*(up(i+1,j-1,k)-3.e0*up(i,j-1,k)         &
     &                +3.e0*up(i-1,j-1,k)-up(i-2,j-1,k))

                    b=dxt22(i)*(up(i+1,j-1,k)                           &
     &                -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*up(i,j-1,k)+2.e0*up(i+1,j-1,k)     &
     &                +up(i-2,j-1,k)-6.e0*up(i-1,j-1,k))

                    advjm1=up(i,j-1,k)+((a*um+b)*um+c)*um

                    a=dxt36w(i)*(up(i+1,j,k)-3.e0*up(i,j,k)             &
     &                +3.e0*up(i-1,j,k)-up(i-2,j,k))

                    b=dxt22(i)*(up(i+1,j,k)                             &
     &                -2.e0*up(i,j,k)+up(i-1,j,k))

                    c=dxt16(i)*(3.e0*up(i,j,k)+2.e0*up(i+1,j,k)         &
     &                +up(i-2,j,k)-6.e0*up(i-1,j,k))

                    advij=up(i,j,k)+((a*um+b)*um+c)*um

                    a=dxt36w(i)*(up(i+1,j+1,k)-3.e0*up(i,j+1,k)         &
     &                +3.e0*up(i-1,j+1,k)-up(i-2,j+1,k))

                    b=dxt22(i)*(up(i+1,j+1,k)                           &
     &                -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*up(i,j+1,k)+2.e0*up(i+1,j+1,k)     &
     &                +up(i-2,j+1,k)-6.e0*up(i-1,j+1,k))

                    advjp1=up(i,j+1,k)+((a*um+b)*um+c)*um

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                  else

                    a=dxt36w(i)*(up(i+1,j-1,k)-3.e0*up(i,j-1,k)         &
     &                +3.e0*up(i-1,j-1,k)-up(i-2,j-1,k))

                    b=dxt22(i)*(up(i+1,j-1,k)                           &
     &                -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*up(i,j-1,k)+2.e0*up(i+1,j-1,k)     &
     &                +up(i-2,j-1,k)-6.e0*up(i-1,j-1,k))

                    advjm1=up(i,j-1,k)+((a*um+b)*um+c)*um

                    a=dxt36w(i)*(up(i+1,j,k)-3.e0*up(i,j,k)             &
     &                +3.e0*up(i-1,j,k)-up(i-2,j,k))

                    b=dxt22(i)*(up(i+1,j,k)                             &
     &                -2.e0*up(i,j,k)+up(i-1,j,k))

                    c=dxt16(i)*(3.e0*up(i,j,k)+2.e0*up(i+1,j,k)         &
     &                +up(i-2,j,k)-6.e0*up(i-1,j,k))

                    advij=up(i,j,k)+((a*um+b)*um+c)*um

                    a=dxt36w(i)*(up(i+1,j+1,k)-3.e0*up(i,j+1,k)         &
     &                +3.e0*up(i-1,j+1,k)-up(i-2,j+1,k))

                    b=dxt22(i)*(up(i+1,j+1,k)                           &
     &                -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*up(i,j+1,k)+2.e0*up(i+1,j+1,k)     &
     &                +up(i-2,j+1,k)-6.e0*up(i-1,j+1,k))

                    advjp1=up(i,j+1,k)+((a*um+b)*um+c)*um

                    a=dxt36w(i)*(up(i+1,j+2,k)-3.e0*up(i,j+2,k)         &
     &                +3.e0*up(i-1,j+2,k)-up(i-2,j+2,k))

                    b=dxt22(i)*(up(i+1,j+2,k)                           &
     &                -2.e0*up(i,j+2,k)+up(i-1,j+2,k))

                    c=dxt16(i)*(3.e0*up(i,j+2,k)+2.e0*up(i+1,j+2,k)     &
     &                +up(i-2,j+2,k)-6.e0*up(i-1,j+2,k))

                    advjp2=up(i,j+2,k)+((a*um+b)*um+c)*um

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                  end if

                else

                  if(v8u.gt.0.e0) then

                    a=dxt36e(i)*(up(i+2,j-2,k)-3.e0*up(i+1,j-2,k)       &
     &                +3.e0*up(i,j-2,k)-up(i-1,j-2,k))

                    b=dxt22(i)*(up(i+1,j-2,k)                           &
     &                -2.e0*up(i,j-2,k)+up(i-1,j-2,k))

                    c=dxt16(i)*(6.e0*up(i+1,j-2,k)-up(i+2,j-2,k)        &
     &                -2.e0*up(i-1,j-2,k)-3.e0*up(i,j-2,k))

                    advjm2=up(i,j-2,k)+((a*um+b)*um+c)*um

                    a=dxt36e(i)*(up(i+2,j-1,k)-3.e0*up(i+1,j-1,k)       &
     &                +3.e0*up(i,j-1,k)-up(i-1,j-1,k))

                    b=dxt22(i)*(up(i+1,j-1,k)                           &
     &                -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*up(i+1,j-1,k)-up(i+2,j-1,k)        &
     &                -2.e0*up(i-1,j-1,k)-3.e0*up(i,j-1,k))

                    advjm1=up(i,j-1,k)+((a*um+b)*um+c)*um

                    a=dxt36e(i)*(up(i+2,j,k)-3.e0*up(i+1,j,k)           &
     &                +3.e0*up(i,j,k)-up(i-1,j,k))

                    b=dxt22(i)*(up(i+1,j,k)                             &
     &                -2.e0*up(i,j,k)+up(i-1,j,k))

                    c=dxt16(i)*(6.e0*up(i+1,j,k)-up(i+2,j,k)            &
     &                -2.e0*up(i-1,j,k)-3.e0*up(i,j,k))

                    advij=up(i,j,k)+((a*um+b)*um+c)*um

                    a=dxt36e(i)*(up(i+2,j+1,k)-3.e0*up(i+1,j+1,k)       &
     &                +3.e0*up(i,j+1,k)-up(i-1,j+1,k))

                    b=dxt22(i)*(up(i+1,j+1,k)                           &
     &                -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*up(i+1,j+1,k)-up(i+2,j+1,k)        &
     &                -2.e0*up(i-1,j+1,k)-3.e0*up(i,j+1,k))

                    advjp1=up(i,j+1,k)+((a*um+b)*um+c)*um

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                  else

                    a=dxt36e(i)*(up(i+2,j-1,k)-3.e0*up(i+1,j-1,k)       &
     &                +3.e0*up(i,j-1,k)-up(i-1,j-1,k))

                    b=dxt22(i)*(up(i+1,j-1,k)                           &
     &                -2.e0*up(i,j-1,k)+up(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*up(i+1,j-1,k)-up(i+2,j-1,k)        &
     &                -2.e0*up(i-1,j-1,k)-3.e0*up(i,j-1,k))

                    advjm1=up(i,j-1,k)+((a*um+b)*um+c)*um

                    a=dxt36e(i)*(up(i+2,j,k)-3.e0*up(i+1,j,k)           &
     &                +3.e0*up(i,j,k)-up(i-1,j,k))

                    b=dxt22(i)*(up(i+1,j,k)                             &
     &                -2.e0*up(i,j,k)+up(i-1,j,k))

                    c=dxt16(i)*(6.e0*up(i+1,j,k)-up(i+2,j,k)            &
     &                -2.e0*up(i-1,j,k)-3.e0*up(i,j,k))

                    advij=up(i,j,k)+((a*um+b)*um+c)*um

                    a=dxt36e(i)*(up(i+2,j+1,k)-3.e0*up(i+1,j+1,k)       &
     &                +3.e0*up(i,j+1,k)-up(i-1,j+1,k))

                    b=dxt22(i)*(up(i+1,j+1,k)                           &
     &                -2.e0*up(i,j+1,k)+up(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*up(i+1,j+1,k)-up(i+2,j+1,k)        &
     &                -2.e0*up(i-1,j+1,k)-3.e0*up(i,j+1,k))

                    advjp1=up(i,j+1,k)+((a*um+b)*um+c)*um

                    a=dxt36e(i)*(up(i+2,j+2,k)-3.e0*up(i+1,j+2,k)       &
     &                +3.e0*up(i,j+2,k)-up(i-1,j+2,k))

                    b=dxt22(i)*(up(i+1,j+2,k)                           &
     &                -2.e0*up(i,j+2,k)+up(i-1,j+2,k))

                    c=dxt16(i)*(6.e0*up(i+1,j+2,k)-up(i+2,j+2,k)        &
     &                -2.e0*up(i-1,j+2,k)-3.e0*up(i,j+2,k))

                    advjp2=up(i,j+2,k)+((a*um+b)*um+c)*um

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    uf(i,j,k)=advij+((a*v8u+b)*v8u+c)*v8u

                  end if

                end if

              end do
              end do

!$omp end do

            end do

          end if

! -----

        end if

!! -----

!!! -----

!!! Perform x-y-z seperated Cubic Lagrange scheme.

      else if(advopt.eq.5) then

! Perform calculation without the map scale factor.

        if(mfcopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,a,b,c)

            do j=0,nj
            do i=2,ni-1

              if(up(i,j,k).gt.0.e0) then

                a=dxt36w(i)*(up(i+1,j,k)-3.e0*up(i,j,k)                 &
     &            +3.e0*up(i-1,j,k)-up(i-2,j,k))

                b=dxt22(i)*(up(i+1,j,k)                                 &
     &            -2.e0*up(i,j,k)+up(i-1,j,k))

                c=dxt16(i)*(3.e0*up(i,j,k)+2.e0*up(i+1,j,k)             &
     &            +up(i-2,j,k)-6.e0*up(i-1,j,k))

                advd(i,j,k)=up(i,j,k)                                   &
     &            +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

              else

                a=dxt36e(i)*(up(i+2,j,k)-3.e0*up(i+1,j,k)               &
     &            +3.e0*up(i,j,k)-up(i-1,j,k))

                b=dxt22(i)*(up(i+1,j,k)                                 &
     &            -2.e0*up(i,j,k)+up(i-1,j,k))

                c=dxt16(i)*(6.e0*up(i+1,j,k)-up(i+2,j,k)                &
     &            -2.e0*up(i-1,j,k)-3.e0*up(i,j,k))

                advd(i,j,k)=up(i,j,k)                                   &
     &            +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

              end if

            end do
            end do

!$omp end do

          end do

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,v8u,a,b,c)

            do j=2,nj-2
            do i=2,ni-1

              v8u=.25e0                                                 &
     &          *((vp(i-1,j,k)+vp(i,j,k))+(vp(i-1,j+1,k)+vp(i,j+1,k)))

              if(v8u.gt.0.e0) then

                a=dyt36s(j)*(advd(i,j+1,k)-3.e0*advd(i,j,k)             &
     &            +3.e0*advd(i,j-1,k)-advd(i,j-2,k))

                b=dyt22(j)*(advd(i,j+1,k)                               &
     &            -2.e0*advd(i,j,k)+advd(i,j-1,k))

                c=dyt16(j)*(3.e0*advd(i,j,k)+2.e0*advd(i,j+1,k)         &
     &            +advd(i,j-2,k)-6.e0*advd(i,j-1,k))

                uf(i,j,k)=advd(i,j,k)+((a*v8u+b)*v8u+c)*v8u

              else

                a=dyt36n(j)*(advd(i,j+2,k)-3.e0*advd(i,j+1,k)           &
     &            +3.e0*advd(i,j,k)-advd(i,j-1,k))

                b=dyt22(j)*(advd(i,j+1,k)                               &
     &            -2.e0*advd(i,j,k)+advd(i,j-1,k))

                c=dyt16(j)*(6.e0*advd(i,j+1,k)-advd(i,j+2,k)            &
     &            -2.e0*advd(i,j-1,k)-3.e0*advd(i,j,k))

                uf(i,j,k)=advd(i,j,k)+((a*v8u+b)*v8u+c)*v8u

              end if

            end do
            end do

!$omp end do

          end do

! -----

!! Perform calculation with the map scale factor.

        else

! Applied the map scale factor for the x direction.

          if(mpopt.eq.0.or.mpopt.eq.10) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,um,a,b,c)

              do j=0,nj
              do i=2,ni-1

                um=mf8u(i,j)*up(i,j,k)

                if(um.gt.0.e0) then

                  a=dxt36w(i)*(up(i+1,j,k)-3.e0*up(i,j,k)               &
     &              +3.e0*up(i-1,j,k)-up(i-2,j,k))

                  b=dxt22(i)*(up(i+1,j,k)                               &
     &              -2.e0*up(i,j,k)+up(i-1,j,k))

                  c=dxt16(i)*(3.e0*up(i,j,k)+2.e0*up(i+1,j,k)           &
     &              +up(i-2,j,k)-6.e0*up(i-1,j,k))

                  advd(i,j,k)=up(i,j,k)+((a*um+b)*um+c)*um

                else

                  a=dxt36e(i)*(up(i+2,j,k)-3.e0*up(i+1,j,k)             &
     &              +3.e0*up(i,j,k)-up(i-1,j,k))

                  b=dxt22(i)*(up(i+1,j,k)                               &
     &              -2.e0*up(i,j,k)+up(i-1,j,k))

                  c=dxt16(i)*(6.e0*up(i+1,j,k)-up(i+2,j,k)              &
     &              -2.e0*up(i-1,j,k)-3.e0*up(i,j,k))

                  advd(i,j,k)=up(i,j,k)+((a*um+b)*um+c)*um

                end if

              end do
              end do

!$omp end do

            end do

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,v8u,a,b,c)

              do j=2,nj-2
              do i=2,ni-1

                v8u=.25e0                                               &
     &            *((vp(i-1,j,k)+vp(i,j,k))+(vp(i-1,j+1,k)+vp(i,j+1,k)))

                if(v8u.gt.0.e0) then

                  a=dyt36s(j)*(advd(i,j+1,k)-3.e0*advd(i,j,k)           &
     &              +3.e0*advd(i,j-1,k)-advd(i,j-2,k))

                  b=dyt22(j)*(advd(i,j+1,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i,j-1,k))

                  c=dyt16(j)*(3.e0*advd(i,j,k)+2.e0*advd(i,j+1,k)       &
     &              +advd(i,j-2,k)-6.e0*advd(i,j-1,k))

                  uf(i,j,k)=advd(i,j,k)+((a*v8u+b)*v8u+c)*v8u

                else

                  a=dyt36n(j)*(advd(i,j+2,k)-3.e0*advd(i,j+1,k)         &
     &              +3.e0*advd(i,j,k)-advd(i,j-1,k))

                  b=dyt22(j)*(advd(i,j+1,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i,j-1,k))

                  c=dyt16(j)*(6.e0*advd(i,j+1,k)-advd(i,j+2,k)          &
     &              -2.e0*advd(i,j-1,k)-3.e0*advd(i,j,k))

                  uf(i,j,k)=advd(i,j,k)+((a*v8u+b)*v8u+c)*v8u

                end if

              end do
              end do

!$omp end do

            end do

! -----

! Applied the map scale factor for the y direction.

          else if(mpopt.eq.5) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,a,b,c)

              do j=0,nj
              do i=2,ni-1

                if(up(i,j,k).gt.0.e0) then

                  a=dxt36w(i)*(up(i+1,j,k)-3.e0*up(i,j,k)               &
     &              +3.e0*up(i-1,j,k)-up(i-2,j,k))

                  b=dxt22(i)*(up(i+1,j,k)                               &
     &              -2.e0*up(i,j,k)+up(i-1,j,k))

                  c=dxt16(i)*(3.e0*up(i,j,k)+2.e0*up(i+1,j,k)           &
     &              +up(i-2,j,k)-6.e0*up(i-1,j,k))

                  advd(i,j,k)=up(i,j,k)                                 &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                else

                  a=dxt36e(i)*(up(i+2,j,k)-3.e0*up(i+1,j,k)             &
     &              +3.e0*up(i,j,k)-up(i-1,j,k))

                  b=dxt22(i)*(up(i+1,j,k)                               &
     &              -2.e0*up(i,j,k)+up(i-1,j,k))

                  c=dxt16(i)*(6.e0*up(i+1,j,k)-up(i+2,j,k)              &
     &              -2.e0*up(i-1,j,k)-3.e0*up(i,j,k))

                  advd(i,j,k)=up(i,j,k)                                 &
     &              +((a*up(i,j,k)+b)*up(i,j,k)+c)*up(i,j,k)

                end if

              end do
              end do

!$omp end do

            end do

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,v8u,a,b,c)

              do j=2,nj-2
              do i=2,ni-1

                v8u=.25e0*((mf8v(i-1,j)*vp(i-1,j,k)+mf8v(i,j)*vp(i,j,k))&
     &           +(mf8v(i-1,j+1)*vp(i-1,j+1,k)+mf8v(i,j+1)*vp(i,j+1,k)))

                if(v8u.gt.0.e0) then

                  a=dyt36s(j)*(advd(i,j+1,k)-3.e0*advd(i,j,k)           &
     &              +3.e0*advd(i,j-1,k)-advd(i,j-2,k))

                  b=dyt22(j)*(advd(i,j+1,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i,j-1,k))

                  c=dyt16(j)*(3.e0*advd(i,j,k)+2.e0*advd(i,j+1,k)       &
     &              +advd(i,j-2,k)-6.e0*advd(i,j-1,k))

                  uf(i,j,k)=advd(i,j,k)+((a*v8u+b)*v8u+c)*v8u

                else

                  a=dyt36n(j)*(advd(i,j+2,k)-3.e0*advd(i,j+1,k)         &
     &              +3.e0*advd(i,j,k)-advd(i,j-1,k))

                  b=dyt22(j)*(advd(i,j+1,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i,j-1,k))

                  c=dyt16(j)*(6.e0*advd(i,j+1,k)-advd(i,j+2,k)          &
     &              -2.e0*advd(i,j-1,k)-3.e0*advd(i,j,k))

                  uf(i,j,k)=advd(i,j,k)+((a*v8u+b)*v8u+c)*v8u

                end if

              end do
              end do

!$omp end do

            end do

! -----

! Applied the map scale factor for the x and y direction.

          else

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,um,a,b,c)

              do j=0,nj
              do i=2,ni-1

                um=mf8u(i,j)*up(i,j,k)

                if(um.gt.0.e0) then

                  a=dxt36w(i)*(up(i+1,j,k)-3.e0*up(i,j,k)               &
     &              +3.e0*up(i-1,j,k)-up(i-2,j,k))

                  b=dxt22(i)*(up(i+1,j,k)                               &
     &              -2.e0*up(i,j,k)+up(i-1,j,k))

                  c=dxt16(i)*(3.e0*up(i,j,k)+2.e0*up(i+1,j,k)           &
     &              +up(i-2,j,k)-6.e0*up(i-1,j,k))

                  advd(i,j,k)=up(i,j,k)+((a*um+b)*um+c)*um

                else

                  a=dxt36e(i)*(up(i+2,j,k)-3.e0*up(i+1,j,k)             &
     &              +3.e0*up(i,j,k)-up(i-1,j,k))

                  b=dxt22(i)*(up(i+1,j,k)                               &
     &              -2.e0*up(i,j,k)+up(i-1,j,k))

                  c=dxt16(i)*(6.e0*up(i+1,j,k)-up(i+2,j,k)              &
     &              -2.e0*up(i-1,j,k)-3.e0*up(i,j,k))

                  advd(i,j,k)=up(i,j,k)+((a*um+b)*um+c)*um

                end if

              end do
              end do

!$omp end do

            end do

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,v8u,a,b,c)

              do j=2,nj-2
              do i=2,ni-1

                v8u=.25e0*((mf8v(i-1,j)*vp(i-1,j,k)+mf8v(i,j)*vp(i,j,k))&
     &           +(mf8v(i-1,j+1)*vp(i-1,j+1,k)+mf8v(i,j+1)*vp(i,j+1,k)))

                if(v8u.gt.0.e0) then

                  a=dyt36s(j)*(advd(i,j+1,k)-3.e0*advd(i,j,k)           &
     &              +3.e0*advd(i,j-1,k)-advd(i,j-2,k))

                  b=dyt22(j)*(advd(i,j+1,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i,j-1,k))

                  c=dyt16(j)*(3.e0*advd(i,j,k)+2.e0*advd(i,j+1,k)       &
     &              +advd(i,j-2,k)-6.e0*advd(i,j-1,k))

                  uf(i,j,k)=advd(i,j,k)+((a*v8u+b)*v8u+c)*v8u

                else

                  a=dyt36n(j)*(advd(i,j+2,k)-3.e0*advd(i,j+1,k)         &
     &              +3.e0*advd(i,j,k)-advd(i,j-1,k))

                  b=dyt22(j)*(advd(i,j+1,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i,j-1,k))

                  c=dyt16(j)*(6.e0*advd(i,j+1,k)-advd(i,j+2,k)          &
     &              -2.e0*advd(i,j-1,k)-3.e0*advd(i,j,k))

                  uf(i,j,k)=advd(i,j,k)+((a*v8u+b)*v8u+c)*v8u

                end if

              end do
              end do

!$omp end do

            end do

          end if

! -----

        end if

!! -----

      end if

!!! -----

!$omp end parallel

!!!! -----

!!!!! -----

!!!!! Calculate the v advection horizontally.

! Set the common used variables.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i)

      do i=0,ni
        dxt22(i)=cdxt2
        dxt16(i)=cdxt1

        dxt36w(i)=cdxt3
        dxt36e(i)=cdxt3

      end do

!$omp end do

!$omp do schedule(runtime) private(j)

      do j=0,nj+1
        dyt22(j)=cdyt2
        dyt16(j)=cdyt1

        dyt36s(j)=cdyt3
        dyt36n(j)=cdyt3

      end do

!$omp end do

!$omp end parallel

! -----

! Set the lateral boundary conditions.

      call lbculv(idwbc,idebc,idsbc,idnbc,ni,nj,nk,vp,                  &
     &            dxt36w,dxt36e,dyt36s,dyt36n)

! -----

!!!! Perform Cubic Lagrange scheme.

!$omp parallel default(shared) private(k)

!!! Perform horizontal-vertical seperated Cubic Lagrange scheme.

      if(advopt.eq.4) then

! Perform calculation without the map scale factor.

        if(mfcopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,u8v,advij,advim2,advim1,advip1,advip2,a,b,c)

            do j=2,nj-1
            do i=2,ni-2

              u8v=.25e0                                                 &
     &          *((up(i,j-1,k)+up(i,j,k))+(up(i+1,j-1,k)+up(i+1,j,k)))

              if(vp(i,j,k).gt.0.e0) then

                if(u8v.gt.0.e0) then

                  a=dyt36s(j)*(vp(i-2,j+1,k)-3.e0*vp(i-2,j,k)           &
     &              +3.e0*vp(i-2,j-1,k)-vp(i-2,j-2,k))

                  b=dyt22(j)*(vp(i-2,j+1,k)                             &
     &              -2.e0*vp(i-2,j,k)+vp(i-2,j-1,k))

                  c=dyt16(j)*(3.e0*vp(i-2,j,k)+2.e0*vp(i-2,j+1,k)       &
     &              +vp(i-2,j-2,k)-6.e0*vp(i-2,j-1,k))

                  advim2=vp(i-2,j,k)                                    &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dyt36s(j)*(vp(i-1,j+1,k)-3.e0*vp(i-1,j,k)           &
     &              +3.e0*vp(i-1,j-1,k)-vp(i-1,j-2,k))

                  b=dyt22(j)*(vp(i-1,j+1,k)                             &
     &              -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                  c=dyt16(j)*(3.e0*vp(i-1,j,k)+2.e0*vp(i-1,j+1,k)       &
     &              +vp(i-1,j-2,k)-6.e0*vp(i-1,j-1,k))

                  advim1=vp(i-1,j,k)                                    &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dyt36s(j)*(vp(i,j+1,k)-3.e0*vp(i,j,k)               &
     &              +3.e0*vp(i,j-1,k)-vp(i,j-2,k))

                  b=dyt22(j)*(vp(i,j+1,k)                               &
     &              -2.e0*vp(i,j,k)+vp(i,j-1,k))

                  c=dyt16(j)*(3.e0*vp(i,j,k)+2.e0*vp(i,j+1,k)           &
     &              +vp(i,j-2,k)-6.e0*vp(i,j-1,k))

                  advij=vp(i,j,k)                                       &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dyt36s(j)*(vp(i+1,j+1,k)-3.e0*vp(i+1,j,k)           &
     &              +3.e0*vp(i+1,j-1,k)-vp(i+1,j-2,k))

                  b=dyt22(j)*(vp(i+1,j+1,k)                             &
     &              -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                  c=dyt16(j)*(3.e0*vp(i+1,j,k)+2.e0*vp(i+1,j+1,k)       &
     &              +vp(i+1,j-2,k)-6.e0*vp(i+1,j-1,k))

                  advip1=vp(i+1,j,k)                                    &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dxt36w(i)*(advip1-3.e0*advij+3.e0*advim1-advim2)

                  b=dxt22(i)*(advip1-2.e0*advij+advim1)

                  c=dxt16(i)                                            &
     &              *(3.e0*advij+2.e0*advip1+advim2-6.e0*advim1)

                  vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                else

                  a=dyt36s(j)*(vp(i-1,j+1,k)-3.e0*vp(i-1,j,k)           &
     &              +3.e0*vp(i-1,j-1,k)-vp(i-1,j-2,k))

                  b=dyt22(j)*(vp(i-1,j+1,k)                             &
     &              -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                  c=dyt16(j)*(3.e0*vp(i-1,j,k)+2.e0*vp(i-1,j+1,k)       &
     &              +vp(i-1,j-2,k)-6.e0*vp(i-1,j-1,k))

                  advim1=vp(i-1,j,k)                                    &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dyt36s(j)*(vp(i,j+1,k)-3.e0*vp(i,j,k)               &
     &              +3.e0*vp(i,j-1,k)-vp(i,j-2,k))

                  b=dyt22(j)*(vp(i,j+1,k)                               &
     &              -2.e0*vp(i,j,k)+vp(i,j-1,k))

                  c=dyt16(j)*(3.e0*vp(i,j,k)+2.e0*vp(i,j+1,k)           &
     &              +vp(i,j-2,k)-6.e0*vp(i,j-1,k))

                  advij=vp(i,j,k)                                       &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dyt36s(j)*(vp(i+1,j+1,k)-3.e0*vp(i+1,j,k)           &
     &              +3.e0*vp(i+1,j-1,k)-vp(i+1,j-2,k))

                  b=dyt22(j)*(vp(i+1,j+1,k)                             &
     &              -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                  c=dyt16(j)*(3.e0*vp(i+1,j,k)+2.e0*vp(i+1,j+1,k)       &
     &              +vp(i+1,j-2,k)-6.e0*vp(i+1,j-1,k))

                  advip1=vp(i+1,j,k)                                    &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dyt36s(j)*(vp(i+2,j+1,k)-3.e0*vp(i+2,j,k)           &
     &              +3.e0*vp(i+2,j-1,k)-vp(i+2,j-2,k))

                  b=dyt22(j)*(vp(i+2,j+1,k)                             &
     &              -2.e0*vp(i+2,j,k)+vp(i+2,j-1,k))

                  c=dyt16(j)*(3.e0*vp(i+2,j,k)+2.e0*vp(i+2,j+1,k)       &
     &              +vp(i+2,j-2,k)-6.e0*vp(i+2,j-1,k))

                  advip2=vp(i+2,j,k)                                    &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dxt36e(i)*(advip2-3.e0*advip1+3.e0*advij-advim1)

                  b=dxt22(i)*(advip1-2.e0*advij+advim1)

                  c=dxt16(i)                                            &
     &              *(6.e0*advip1-advip2-2.e0*advim1-3.e0*advij)

                  vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                end if

              else

                if(u8v.gt.0.e0) then

                  a=dyt36n(j)*(vp(i-2,j+2,k)-3.e0*vp(i-2,j+1,k)         &
     &              +3.e0*vp(i-2,j,k)-vp(i-2,j-1,k))

                  b=dyt22(j)*(vp(i-2,j+1,k)                             &
     &              -2.e0*vp(i-2,j,k)+vp(i-2,j-1,k))

                  c=dyt16(j)*(6.e0*vp(i-2,j+1,k)-vp(i-2,j+2,k)          &
     &              -2.e0*vp(i-2,j-1,k)-3.e0*vp(i-2,j,k))

                  advim2=vp(i-2,j,k)                                    &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dyt36n(j)*(vp(i-1,j+2,k)-3.e0*vp(i-1,j+1,k)         &
     &              +3.e0*vp(i-1,j,k)-vp(i-1,j-1,k))

                  b=dyt22(j)*(vp(i-1,j+1,k)                             &
     &              -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                  c=dyt16(j)*(6.e0*vp(i-1,j+1,k)-vp(i-1,j+2,k)          &
     &              -2.e0*vp(i-1,j-1,k)-3.e0*vp(i-1,j,k))

                  advim1=vp(i-1,j,k)                                    &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dyt36n(j)*(vp(i,j+2,k)-3.e0*vp(i,j+1,k)             &
     &              +3.e0*vp(i,j,k)-vp(i,j-1,k))

                  b=dyt22(j)*(vp(i,j+1,k)                               &
     &              -2.e0*vp(i,j,k)+vp(i,j-1,k))

                  c=dyt16(j)*(6.e0*vp(i,j+1,k)-vp(i,j+2,k)              &
     &              -2.e0*vp(i,j-1,k)-3.e0*vp(i,j,k))

                  advij=vp(i,j,k)                                       &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dyt36n(j)*(vp(i+1,j+2,k)-3.e0*vp(i+1,j+1,k)         &
     &              +3.e0*vp(i+1,j,k)-vp(i+1,j-1,k))

                  b=dyt22(j)*(vp(i+1,j+1,k)                             &
     &              -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                  c=dyt16(j)*(6.e0*vp(i+1,j+1,k)-vp(i+1,j+2,k)          &
     &              -2.e0*vp(i+1,j-1,k)-3.e0*vp(i+1,j,k))

                  advip1=vp(i+1,j,k)                                    &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dxt36w(i)*(advip1-3.e0*advij+3.e0*advim1-advim2)

                  b=dxt22(i)*(advip1-2.e0*advij+advim1)

                  c=dxt16(i)                                            &
     &              *(3.e0*advij+2.e0*advip1+advim2-6.e0*advim1)

                  vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                else

                  a=dyt36n(j)*(vp(i-1,j+2,k)-3.e0*vp(i-1,j+1,k)         &
     &              +3.e0*vp(i-1,j,k)-vp(i-1,j-1,k))

                  b=dyt22(j)*(vp(i-1,j+1,k)                             &
     &              -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                  c=dyt16(j)*(6.e0*vp(i-1,j+1,k)-vp(i-1,j+2,k)          &
     &              -2.e0*vp(i-1,j-1,k)-3.e0*vp(i-1,j,k))

                  advim1=vp(i-1,j,k)                                    &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dyt36n(j)*(vp(i,j+2,k)-3.e0*vp(i,j+1,k)             &
     &              +3.e0*vp(i,j,k)-vp(i,j-1,k))

                  b=dyt22(j)*(vp(i,j+1,k)                               &
     &              -2.e0*vp(i,j,k)+vp(i,j-1,k))

                  c=dyt16(j)*(6.e0*vp(i,j+1,k)-vp(i,j+2,k)              &
     &              -2.e0*vp(i,j-1,k)-3.e0*vp(i,j,k))

                  advij=vp(i,j,k)                                       &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dyt36n(j)*(vp(i+1,j+2,k)-3.e0*vp(i+1,j+1,k)         &
     &              +3.e0*vp(i+1,j,k)-vp(i+1,j-1,k))

                  b=dyt22(j)*(vp(i+1,j+1,k)                             &
     &              -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                  c=dyt16(j)*(6.e0*vp(i+1,j+1,k)-vp(i+1,j+2,k)          &
     &              -2.e0*vp(i+1,j-1,k)-3.e0*vp(i+1,j,k))

                  advip1=vp(i+1,j,k)                                    &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dyt36n(j)*(vp(i+2,j+2,k)-3.e0*vp(i+2,j+1,k)         &
     &              +3.e0*vp(i+2,j,k)-vp(i+2,j-1,k))

                  b=dyt22(j)*(vp(i+2,j+1,k)                             &
     &              -2.e0*vp(i+2,j,k)+vp(i+2,j-1,k))

                  c=dyt16(j)*(6.e0*vp(i+2,j+1,k)-vp(i+2,j+2,k)          &
     &              -2.e0*vp(i+2,j-1,k)-3.e0*vp(i+2,j,k))

                  advip2=vp(i+2,j,k)                                    &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                  a=dxt36e(i)*(advip2-3.e0*advip1+3.e0*advij-advim1)

                  b=dxt22(i)*(advip1-2.e0*advij+advim1)

                  c=dxt16(i)                                            &
     &              *(6.e0*advip1-advip2-2.e0*advim1-3.e0*advij)

                  vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                end if

              end if

            end do
            end do

!$omp end do

          end do

! -----

!! Perform calculation with the map scale factor.

        else

! Applied the map scale factor for the x direction.

          if(mpopt.eq.0.or.mpopt.eq.10) then

            do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,u8v,advij,advim2,advim1,advip1,advip2,a,b,c)

              do j=2,nj-1
              do i=2,ni-2

                u8v=.25e0*((mf8u(i,j-1)*up(i,j-1,k)+mf8u(i,j)*up(i,j,k))&
     &           +(mf8u(i+1,j-1)*up(i+1,j-1,k)+mf8u(i+1,j)*up(i+1,j,k)))

                if(vp(i,j,k).gt.0.e0) then

                  if(u8v.gt.0.e0) then

                    a=dyt36s(j)*(vp(i-2,j+1,k)-3.e0*vp(i-2,j,k)         &
     &                +3.e0*vp(i-2,j-1,k)-vp(i-2,j-2,k))

                    b=dyt22(j)*(vp(i-2,j+1,k)                           &
     &                -2.e0*vp(i-2,j,k)+vp(i-2,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i-2,j,k)+2.e0*vp(i-2,j+1,k)     &
     &                +vp(i-2,j-2,k)-6.e0*vp(i-2,j-1,k))

                    advim2=vp(i-2,j,k)                                  &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dyt36s(j)*(vp(i-1,j+1,k)-3.e0*vp(i-1,j,k)         &
     &                +3.e0*vp(i-1,j-1,k)-vp(i-1,j-2,k))

                    b=dyt22(j)*(vp(i-1,j+1,k)                           &
     &                -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i-1,j,k)+2.e0*vp(i-1,j+1,k)     &
     &                +vp(i-1,j-2,k)-6.e0*vp(i-1,j-1,k))

                    advim1=vp(i-1,j,k)                                  &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dyt36s(j)*(vp(i,j+1,k)-3.e0*vp(i,j,k)             &
     &                +3.e0*vp(i,j-1,k)-vp(i,j-2,k))

                    b=dyt22(j)*(vp(i,j+1,k)                             &
     &                -2.e0*vp(i,j,k)+vp(i,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i,j,k)+2.e0*vp(i,j+1,k)         &
     &                +vp(i,j-2,k)-6.e0*vp(i,j-1,k))

                    advij=vp(i,j,k)                                     &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dyt36s(j)*(vp(i+1,j+1,k)-3.e0*vp(i+1,j,k)         &
     &                +3.e0*vp(i+1,j-1,k)-vp(i+1,j-2,k))

                    b=dyt22(j)*(vp(i+1,j+1,k)                           &
     &                -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i+1,j,k)+2.e0*vp(i+1,j+1,k)     &
     &                +vp(i+1,j-2,k)-6.e0*vp(i+1,j-1,k))

                    advip1=vp(i+1,j,k)                                  &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dxt36w(i)*(advip1-3.e0*advij+3.e0*advim1-advim2)

                    b=dxt22(i)*(advip1-2.e0*advij+advim1)

                    c=dxt16(i)                                          &
     &                *(3.e0*advij+2.e0*advip1+advim2-6.e0*advim1)

                    vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                  else

                    a=dyt36s(j)*(vp(i-1,j+1,k)-3.e0*vp(i-1,j,k)         &
     &                +3.e0*vp(i-1,j-1,k)-vp(i-1,j-2,k))

                    b=dyt22(j)*(vp(i-1,j+1,k)                           &
     &                -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i-1,j,k)+2.e0*vp(i-1,j+1,k)     &
     &                +vp(i-1,j-2,k)-6.e0*vp(i-1,j-1,k))

                    advim1=vp(i-1,j,k)                                  &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dyt36s(j)*(vp(i,j+1,k)-3.e0*vp(i,j,k)             &
     &                +3.e0*vp(i,j-1,k)-vp(i,j-2,k))

                    b=dyt22(j)*(vp(i,j+1,k)                             &
     &                -2.e0*vp(i,j,k)+vp(i,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i,j,k)+2.e0*vp(i,j+1,k)         &
     &                +vp(i,j-2,k)-6.e0*vp(i,j-1,k))

                    advij=vp(i,j,k)                                     &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dyt36s(j)*(vp(i+1,j+1,k)-3.e0*vp(i+1,j,k)         &
     &                +3.e0*vp(i+1,j-1,k)-vp(i+1,j-2,k))

                    b=dyt22(j)*(vp(i+1,j+1,k)                           &
     &                -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i+1,j,k)+2.e0*vp(i+1,j+1,k)     &
     &                +vp(i+1,j-2,k)-6.e0*vp(i+1,j-1,k))

                    advip1=vp(i+1,j,k)                                  &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dyt36s(j)*(vp(i+2,j+1,k)-3.e0*vp(i+2,j,k)         &
     &                +3.e0*vp(i+2,j-1,k)-vp(i+2,j-2,k))

                    b=dyt22(j)*(vp(i+2,j+1,k)                           &
     &                -2.e0*vp(i+2,j,k)+vp(i+2,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i+2,j,k)+2.e0*vp(i+2,j+1,k)     &
     &                +vp(i+2,j-2,k)-6.e0*vp(i+2,j-1,k))

                    advip2=vp(i+2,j,k)                                  &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dxt36e(i)*(advip2-3.e0*advip1+3.e0*advij-advim1)

                    b=dxt22(i)*(advip1-2.e0*advij+advim1)

                    c=dxt16(i)                                          &
     &                *(6.e0*advip1-advip2-2.e0*advim1-3.e0*advij)

                    vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                  end if

                else

                  if(u8v.gt.0.e0) then

                    a=dyt36n(j)*(vp(i-2,j+2,k)-3.e0*vp(i-2,j+1,k)       &
     &                +3.e0*vp(i-2,j,k)-vp(i-2,j-1,k))

                    b=dyt22(j)*(vp(i-2,j+1,k)                           &
     &                -2.e0*vp(i-2,j,k)+vp(i-2,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i-2,j+1,k)-vp(i-2,j+2,k)        &
     &                -2.e0*vp(i-2,j-1,k)-3.e0*vp(i-2,j,k))

                    advim2=vp(i-2,j,k)                                  &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dyt36n(j)*(vp(i-1,j+2,k)-3.e0*vp(i-1,j+1,k)       &
     &                +3.e0*vp(i-1,j,k)-vp(i-1,j-1,k))

                    b=dyt22(j)*(vp(i-1,j+1,k)                           &
     &                -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i-1,j+1,k)-vp(i-1,j+2,k)        &
     &                -2.e0*vp(i-1,j-1,k)-3.e0*vp(i-1,j,k))

                    advim1=vp(i-1,j,k)                                  &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dyt36n(j)*(vp(i,j+2,k)-3.e0*vp(i,j+1,k)           &
     &                +3.e0*vp(i,j,k)-vp(i,j-1,k))

                    b=dyt22(j)*(vp(i,j+1,k)                             &
     &                -2.e0*vp(i,j,k)+vp(i,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i,j+1,k)-vp(i,j+2,k)            &
     &                -2.e0*vp(i,j-1,k)-3.e0*vp(i,j,k))

                    advij=vp(i,j,k)                                     &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dyt36n(j)*(vp(i+1,j+2,k)-3.e0*vp(i+1,j+1,k)       &
     &                +3.e0*vp(i+1,j,k)-vp(i+1,j-1,k))

                    b=dyt22(j)*(vp(i+1,j+1,k)                           &
     &                -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i+1,j+1,k)-vp(i+1,j+2,k)        &
     &                -2.e0*vp(i+1,j-1,k)-3.e0*vp(i+1,j,k))

                    advip1=vp(i+1,j,k)                                  &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dxt36w(i)*(advip1-3.e0*advij+3.e0*advim1-advim2)

                    b=dxt22(i)*(advip1-2.e0*advij+advim1)

                    c=dxt16(i)                                          &
     &                *(3.e0*advij+2.e0*advip1+advim2-6.e0*advim1)

                    vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                  else

                    a=dyt36n(j)*(vp(i-1,j+2,k)-3.e0*vp(i-1,j+1,k)       &
     &                +3.e0*vp(i-1,j,k)-vp(i-1,j-1,k))

                    b=dyt22(j)*(vp(i-1,j+1,k)                           &
     &                -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i-1,j+1,k)-vp(i-1,j+2,k)        &
     &                -2.e0*vp(i-1,j-1,k)-3.e0*vp(i-1,j,k))

                    advim1=vp(i-1,j,k)                                  &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dyt36n(j)*(vp(i,j+2,k)-3.e0*vp(i,j+1,k)           &
     &                +3.e0*vp(i,j,k)-vp(i,j-1,k))

                    b=dyt22(j)*(vp(i,j+1,k)                             &
     &                -2.e0*vp(i,j,k)+vp(i,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i,j+1,k)-vp(i,j+2,k)            &
     &                -2.e0*vp(i,j-1,k)-3.e0*vp(i,j,k))

                    advij=vp(i,j,k)                                     &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dyt36n(j)*(vp(i+1,j+2,k)-3.e0*vp(i+1,j+1,k)       &
     &                +3.e0*vp(i+1,j,k)-vp(i+1,j-1,k))

                    b=dyt22(j)*(vp(i+1,j+1,k)                           &
     &                -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i+1,j+1,k)-vp(i+1,j+2,k)        &
     &                -2.e0*vp(i+1,j-1,k)-3.e0*vp(i+1,j,k))

                    advip1=vp(i+1,j,k)                                  &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dyt36n(j)*(vp(i+2,j+2,k)-3.e0*vp(i+2,j+1,k)       &
     &                +3.e0*vp(i+2,j,k)-vp(i+2,j-1,k))

                    b=dyt22(j)*(vp(i+2,j+1,k)                           &
     &                -2.e0*vp(i+2,j,k)+vp(i+2,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i+2,j+1,k)-vp(i+2,j+2,k)        &
     &                -2.e0*vp(i+2,j-1,k)-3.e0*vp(i+2,j,k))

                    advip2=vp(i+2,j,k)                                  &
     &                +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                    a=dxt36e(i)*(advip2-3.e0*advip1+3.e0*advij-advim1)

                    b=dxt22(i)*(advip1-2.e0*advij+advim1)

                    c=dxt16(i)                                          &
     &                *(6.e0*advip1-advip2-2.e0*advim1-3.e0*advij)

                    vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                  end if

                end if

              end do
              end do

!$omp end do

            end do

! -----

! Applied the map scale factor for the y direction.

          else if(mpopt.eq.5) then

            do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,vm,u8v,advij,advim2,advim1,advip1,advip2,a,b,c)

              do j=2,nj-1
              do i=2,ni-2

                vm=mf8v(i,j)*vp(i,j,k)

                u8v=.25e0                                               &
     &            *((up(i,j-1,k)+up(i,j,k))+(up(i+1,j-1,k)+up(i+1,j,k)))

                if(vm.gt.0.e0) then

                  if(u8v.gt.0.e0) then

                    a=dyt36s(j)*(vp(i-2,j+1,k)-3.e0*vp(i-2,j,k)         &
     &                +3.e0*vp(i-2,j-1,k)-vp(i-2,j-2,k))

                    b=dyt22(j)*(vp(i-2,j+1,k)                           &
     &                -2.e0*vp(i-2,j,k)+vp(i-2,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i-2,j,k)+2.e0*vp(i-2,j+1,k)     &
     &                +vp(i-2,j-2,k)-6.e0*vp(i-2,j-1,k))

                    advim2=vp(i-2,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36s(j)*(vp(i-1,j+1,k)-3.e0*vp(i-1,j,k)         &
     &                +3.e0*vp(i-1,j-1,k)-vp(i-1,j-2,k))

                    b=dyt22(j)*(vp(i-1,j+1,k)                           &
     &                -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i-1,j,k)+2.e0*vp(i-1,j+1,k)     &
     &                +vp(i-1,j-2,k)-6.e0*vp(i-1,j-1,k))

                    advim1=vp(i-1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36s(j)*(vp(i,j+1,k)-3.e0*vp(i,j,k)             &
     &                +3.e0*vp(i,j-1,k)-vp(i,j-2,k))

                    b=dyt22(j)*(vp(i,j+1,k)                             &
     &                -2.e0*vp(i,j,k)+vp(i,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i,j,k)+2.e0*vp(i,j+1,k)         &
     &                +vp(i,j-2,k)-6.e0*vp(i,j-1,k))

                    advij=vp(i,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36s(j)*(vp(i+1,j+1,k)-3.e0*vp(i+1,j,k)         &
     &                +3.e0*vp(i+1,j-1,k)-vp(i+1,j-2,k))

                    b=dyt22(j)*(vp(i+1,j+1,k)                           &
     &                -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i+1,j,k)+2.e0*vp(i+1,j+1,k)     &
     &                +vp(i+1,j-2,k)-6.e0*vp(i+1,j-1,k))

                    advip1=vp(i+1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dxt36w(i)*(advip1-3.e0*advij+3.e0*advim1-advim2)

                    b=dxt22(i)*(advip1-2.e0*advij+advim1)

                    c=dxt16(i)                                          &
     &                *(3.e0*advij+2.e0*advip1+advim2-6.e0*advim1)

                    vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                  else

                    a=dyt36s(j)*(vp(i-1,j+1,k)-3.e0*vp(i-1,j,k)         &
     &                +3.e0*vp(i-1,j-1,k)-vp(i-1,j-2,k))

                    b=dyt22(j)*(vp(i-1,j+1,k)                           &
     &                -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i-1,j,k)+2.e0*vp(i-1,j+1,k)     &
     &                +vp(i-1,j-2,k)-6.e0*vp(i-1,j-1,k))

                    advim1=vp(i-1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36s(j)*(vp(i,j+1,k)-3.e0*vp(i,j,k)             &
     &                +3.e0*vp(i,j-1,k)-vp(i,j-2,k))

                    b=dyt22(j)*(vp(i,j+1,k)                             &
     &                -2.e0*vp(i,j,k)+vp(i,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i,j,k)+2.e0*vp(i,j+1,k)         &
     &                +vp(i,j-2,k)-6.e0*vp(i,j-1,k))

                    advij=vp(i,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36s(j)*(vp(i+1,j+1,k)-3.e0*vp(i+1,j,k)         &
     &                +3.e0*vp(i+1,j-1,k)-vp(i+1,j-2,k))

                    b=dyt22(j)*(vp(i+1,j+1,k)                           &
     &                -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i+1,j,k)+2.e0*vp(i+1,j+1,k)     &
     &                +vp(i+1,j-2,k)-6.e0*vp(i+1,j-1,k))

                    advip1=vp(i+1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36s(j)*(vp(i+2,j+1,k)-3.e0*vp(i+2,j,k)         &
     &                +3.e0*vp(i+2,j-1,k)-vp(i+2,j-2,k))

                    b=dyt22(j)*(vp(i+2,j+1,k)                           &
     &                -2.e0*vp(i+2,j,k)+vp(i+2,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i+2,j,k)+2.e0*vp(i+2,j+1,k)     &
     &                +vp(i+2,j-2,k)-6.e0*vp(i+2,j-1,k))

                    advip2=vp(i+2,j,k)+((a*vm+b)*vm+c)*vm

                    a=dxt36e(i)*(advip2-3.e0*advip1+3.e0*advij-advim1)

                    b=dxt22(i)*(advip1-2.e0*advij+advim1)

                    c=dxt16(i)                                          &
     &                *(6.e0*advip1-advip2-2.e0*advim1-3.e0*advij)

                    vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                  end if

                else

                  if(u8v.gt.0.e0) then

                    a=dyt36n(j)*(vp(i-2,j+2,k)-3.e0*vp(i-2,j+1,k)       &
     &                +3.e0*vp(i-2,j,k)-vp(i-2,j-1,k))

                    b=dyt22(j)*(vp(i-2,j+1,k)                           &
     &                -2.e0*vp(i-2,j,k)+vp(i-2,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i-2,j+1,k)-vp(i-2,j+2,k)        &
     &                -2.e0*vp(i-2,j-1,k)-3.e0*vp(i-2,j,k))

                    advim2=vp(i-2,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36n(j)*(vp(i-1,j+2,k)-3.e0*vp(i-1,j+1,k)       &
     &                +3.e0*vp(i-1,j,k)-vp(i-1,j-1,k))

                    b=dyt22(j)*(vp(i-1,j+1,k)                           &
     &                -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i-1,j+1,k)-vp(i-1,j+2,k)        &
     &                -2.e0*vp(i-1,j-1,k)-3.e0*vp(i-1,j,k))

                    advim1=vp(i-1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36n(j)*(vp(i,j+2,k)-3.e0*vp(i,j+1,k)           &
     &                +3.e0*vp(i,j,k)-vp(i,j-1,k))

                    b=dyt22(j)*(vp(i,j+1,k)                             &
     &                -2.e0*vp(i,j,k)+vp(i,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i,j+1,k)-vp(i,j+2,k)            &
     &                -2.e0*vp(i,j-1,k)-3.e0*vp(i,j,k))

                    advij=vp(i,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36n(j)*(vp(i+1,j+2,k)-3.e0*vp(i+1,j+1,k)       &
     &                +3.e0*vp(i+1,j,k)-vp(i+1,j-1,k))

                    b=dyt22(j)*(vp(i+1,j+1,k)                           &
     &                -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i+1,j+1,k)-vp(i+1,j+2,k)        &
     &                -2.e0*vp(i+1,j-1,k)-3.e0*vp(i+1,j,k))

                    advip1=vp(i+1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dxt36w(i)*(advip1-3.e0*advij+3.e0*advim1-advim2)

                    b=dxt22(i)*(advip1-2.e0*advij+advim1)

                    c=dxt16(i)                                          &
     &                *(3.e0*advij+2.e0*advip1+advim2-6.e0*advim1)

                    vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                  else

                    a=dyt36n(j)*(vp(i-1,j+2,k)-3.e0*vp(i-1,j+1,k)       &
     &                +3.e0*vp(i-1,j,k)-vp(i-1,j-1,k))

                    b=dyt22(j)*(vp(i-1,j+1,k)                           &
     &                -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i-1,j+1,k)-vp(i-1,j+2,k)        &
     &                -2.e0*vp(i-1,j-1,k)-3.e0*vp(i-1,j,k))

                    advim1=vp(i-1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36n(j)*(vp(i,j+2,k)-3.e0*vp(i,j+1,k)           &
     &                +3.e0*vp(i,j,k)-vp(i,j-1,k))

                    b=dyt22(j)*(vp(i,j+1,k)                             &
     &                -2.e0*vp(i,j,k)+vp(i,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i,j+1,k)-vp(i,j+2,k)            &
     &                -2.e0*vp(i,j-1,k)-3.e0*vp(i,j,k))

                    advij=vp(i,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36n(j)*(vp(i+1,j+2,k)-3.e0*vp(i+1,j+1,k)       &
     &                +3.e0*vp(i+1,j,k)-vp(i+1,j-1,k))

                    b=dyt22(j)*(vp(i+1,j+1,k)                           &
     &                -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i+1,j+1,k)-vp(i+1,j+2,k)        &
     &                -2.e0*vp(i+1,j-1,k)-3.e0*vp(i+1,j,k))

                    advip1=vp(i+1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36n(j)*(vp(i+2,j+2,k)-3.e0*vp(i+2,j+1,k)       &
     &                +3.e0*vp(i+2,j,k)-vp(i+2,j-1,k))

                    b=dyt22(j)*(vp(i+2,j+1,k)                           &
     &                -2.e0*vp(i+2,j,k)+vp(i+2,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i+2,j+1,k)-vp(i+2,j+2,k)        &
     &                -2.e0*vp(i+2,j-1,k)-3.e0*vp(i+2,j,k))

                    advip2=vp(i+2,j,k)+((a*vm+b)*vm+c)*vm

                    a=dxt36e(i)*(advip2-3.e0*advip1+3.e0*advij-advim1)

                    b=dxt22(i)*(advip1-2.e0*advij+advim1)

                    c=dxt16(i)                                          &
     &                *(6.e0*advip1-advip2-2.e0*advim1-3.e0*advij)

                    vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                  end if

                end if

              end do
              end do

!$omp end do

            end do

! -----

! Applied the map scale factor for the x and y direction.

          else

            do k=1,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,vm,u8v,advij,advim2,advim1,advip1,advip2,a,b,c)

              do j=2,nj-1
              do i=2,ni-2

                vm=mf8v(i,j)*vp(i,j,k)

                u8v=.25e0*((mf8u(i,j-1)*up(i,j-1,k)+mf8u(i,j)*up(i,j,k))&
     &           +(mf8u(i+1,j-1)*up(i+1,j-1,k)+mf8u(i+1,j)*up(i+1,j,k)))

                if(vm.gt.0.e0) then

                  if(u8v.gt.0.e0) then

                    a=dyt36s(j)*(vp(i-2,j+1,k)-3.e0*vp(i-2,j,k)         &
     &                +3.e0*vp(i-2,j-1,k)-vp(i-2,j-2,k))

                    b=dyt22(j)*(vp(i-2,j+1,k)                           &
     &                -2.e0*vp(i-2,j,k)+vp(i-2,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i-2,j,k)+2.e0*vp(i-2,j+1,k)     &
     &                +vp(i-2,j-2,k)-6.e0*vp(i-2,j-1,k))

                    advim2=vp(i-2,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36s(j)*(vp(i-1,j+1,k)-3.e0*vp(i-1,j,k)         &
     &                +3.e0*vp(i-1,j-1,k)-vp(i-1,j-2,k))

                    b=dyt22(j)*(vp(i-1,j+1,k)                           &
     &                -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i-1,j,k)+2.e0*vp(i-1,j+1,k)     &
     &                +vp(i-1,j-2,k)-6.e0*vp(i-1,j-1,k))

                    advim1=vp(i-1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36s(j)*(vp(i,j+1,k)-3.e0*vp(i,j,k)             &
     &                +3.e0*vp(i,j-1,k)-vp(i,j-2,k))

                    b=dyt22(j)*(vp(i,j+1,k)                             &
     &                -2.e0*vp(i,j,k)+vp(i,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i,j,k)+2.e0*vp(i,j+1,k)         &
     &                +vp(i,j-2,k)-6.e0*vp(i,j-1,k))

                    advij=vp(i,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36s(j)*(vp(i+1,j+1,k)-3.e0*vp(i+1,j,k)         &
     &                +3.e0*vp(i+1,j-1,k)-vp(i+1,j-2,k))

                    b=dyt22(j)*(vp(i+1,j+1,k)                           &
     &                -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i+1,j,k)+2.e0*vp(i+1,j+1,k)     &
     &                +vp(i+1,j-2,k)-6.e0*vp(i+1,j-1,k))

                    advip1=vp(i+1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dxt36w(i)*(advip1-3.e0*advij+3.e0*advim1-advim2)

                    b=dxt22(i)*(advip1-2.e0*advij+advim1)

                    c=dxt16(i)                                          &
     &                *(3.e0*advij+2.e0*advip1+advim2-6.e0*advim1)

                    vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                  else

                    a=dyt36s(j)*(vp(i-1,j+1,k)-3.e0*vp(i-1,j,k)         &
     &                +3.e0*vp(i-1,j-1,k)-vp(i-1,j-2,k))

                    b=dyt22(j)*(vp(i-1,j+1,k)                           &
     &                -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i-1,j,k)+2.e0*vp(i-1,j+1,k)     &
     &                +vp(i-1,j-2,k)-6.e0*vp(i-1,j-1,k))

                    advim1=vp(i-1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36s(j)*(vp(i,j+1,k)-3.e0*vp(i,j,k)             &
     &                +3.e0*vp(i,j-1,k)-vp(i,j-2,k))

                    b=dyt22(j)*(vp(i,j+1,k)                             &
     &                -2.e0*vp(i,j,k)+vp(i,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i,j,k)+2.e0*vp(i,j+1,k)         &
     &                +vp(i,j-2,k)-6.e0*vp(i,j-1,k))

                    advij=vp(i,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36s(j)*(vp(i+1,j+1,k)-3.e0*vp(i+1,j,k)         &
     &                +3.e0*vp(i+1,j-1,k)-vp(i+1,j-2,k))

                    b=dyt22(j)*(vp(i+1,j+1,k)                           &
     &                -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i+1,j,k)+2.e0*vp(i+1,j+1,k)     &
     &                +vp(i+1,j-2,k)-6.e0*vp(i+1,j-1,k))

                    advip1=vp(i+1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36s(j)*(vp(i+2,j+1,k)-3.e0*vp(i+2,j,k)         &
     &                +3.e0*vp(i+2,j-1,k)-vp(i+2,j-2,k))

                    b=dyt22(j)*(vp(i+2,j+1,k)                           &
     &                -2.e0*vp(i+2,j,k)+vp(i+2,j-1,k))

                    c=dyt16(j)*(3.e0*vp(i+2,j,k)+2.e0*vp(i+2,j+1,k)     &
     &                +vp(i+2,j-2,k)-6.e0*vp(i+2,j-1,k))

                    advip2=vp(i+2,j,k)+((a*vm+b)*vm+c)*vm

                    a=dxt36e(i)*(advip2-3.e0*advip1+3.e0*advij-advim1)

                    b=dxt22(i)*(advip1-2.e0*advij+advim1)

                    c=dxt16(i)                                          &
     &                *(6.e0*advip1-advip2-2.e0*advim1-3.e0*advij)

                    vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                  end if

                else

                  if(u8v.gt.0.e0) then

                    a=dyt36n(j)*(vp(i-2,j+2,k)-3.e0*vp(i-2,j+1,k)       &
     &                +3.e0*vp(i-2,j,k)-vp(i-2,j-1,k))

                    b=dyt22(j)*(vp(i-2,j+1,k)                           &
     &                -2.e0*vp(i-2,j,k)+vp(i-2,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i-2,j+1,k)-vp(i-2,j+2,k)        &
     &                -2.e0*vp(i-2,j-1,k)-3.e0*vp(i-2,j,k))

                    advim2=vp(i-2,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36n(j)*(vp(i-1,j+2,k)-3.e0*vp(i-1,j+1,k)       &
     &                +3.e0*vp(i-1,j,k)-vp(i-1,j-1,k))

                    b=dyt22(j)*(vp(i-1,j+1,k)                           &
     &                -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i-1,j+1,k)-vp(i-1,j+2,k)        &
     &                -2.e0*vp(i-1,j-1,k)-3.e0*vp(i-1,j,k))

                    advim1=vp(i-1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36n(j)*(vp(i,j+2,k)-3.e0*vp(i,j+1,k)           &
     &                +3.e0*vp(i,j,k)-vp(i,j-1,k))

                    b=dyt22(j)*(vp(i,j+1,k)                             &
     &                -2.e0*vp(i,j,k)+vp(i,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i,j+1,k)-vp(i,j+2,k)            &
     &                -2.e0*vp(i,j-1,k)-3.e0*vp(i,j,k))

                    advij=vp(i,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36n(j)*(vp(i+1,j+2,k)-3.e0*vp(i+1,j+1,k)       &
     &                +3.e0*vp(i+1,j,k)-vp(i+1,j-1,k))

                    b=dyt22(j)*(vp(i+1,j+1,k)                           &
     &                -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i+1,j+1,k)-vp(i+1,j+2,k)        &
     &                -2.e0*vp(i+1,j-1,k)-3.e0*vp(i+1,j,k))

                    advip1=vp(i+1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dxt36w(i)*(advip1-3.e0*advij+3.e0*advim1-advim2)

                    b=dxt22(i)*(advip1-2.e0*advij+advim1)

                    c=dxt16(i)                                          &
     &                *(3.e0*advij+2.e0*advip1+advim2-6.e0*advim1)

                    vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                  else

                    a=dyt36n(j)*(vp(i-1,j+2,k)-3.e0*vp(i-1,j+1,k)       &
     &                +3.e0*vp(i-1,j,k)-vp(i-1,j-1,k))

                    b=dyt22(j)*(vp(i-1,j+1,k)                           &
     &                -2.e0*vp(i-1,j,k)+vp(i-1,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i-1,j+1,k)-vp(i-1,j+2,k)        &
     &                -2.e0*vp(i-1,j-1,k)-3.e0*vp(i-1,j,k))

                    advim1=vp(i-1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36n(j)*(vp(i,j+2,k)-3.e0*vp(i,j+1,k)           &
     &                +3.e0*vp(i,j,k)-vp(i,j-1,k))

                    b=dyt22(j)*(vp(i,j+1,k)                             &
     &                -2.e0*vp(i,j,k)+vp(i,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i,j+1,k)-vp(i,j+2,k)            &
     &                -2.e0*vp(i,j-1,k)-3.e0*vp(i,j,k))

                    advij=vp(i,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36n(j)*(vp(i+1,j+2,k)-3.e0*vp(i+1,j+1,k)       &
     &                +3.e0*vp(i+1,j,k)-vp(i+1,j-1,k))

                    b=dyt22(j)*(vp(i+1,j+1,k)                           &
     &                -2.e0*vp(i+1,j,k)+vp(i+1,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i+1,j+1,k)-vp(i+1,j+2,k)        &
     &                -2.e0*vp(i+1,j-1,k)-3.e0*vp(i+1,j,k))

                    advip1=vp(i+1,j,k)+((a*vm+b)*vm+c)*vm

                    a=dyt36n(j)*(vp(i+2,j+2,k)-3.e0*vp(i+2,j+1,k)       &
     &                +3.e0*vp(i+2,j,k)-vp(i+2,j-1,k))

                    b=dyt22(j)*(vp(i+2,j+1,k)                           &
     &                -2.e0*vp(i+2,j,k)+vp(i+2,j-1,k))

                    c=dyt16(j)*(6.e0*vp(i+2,j+1,k)-vp(i+2,j+2,k)        &
     &                -2.e0*vp(i+2,j-1,k)-3.e0*vp(i+2,j,k))

                    advip2=vp(i+2,j,k)+((a*vm+b)*vm+c)*vm

                    a=dxt36e(i)*(advip2-3.e0*advip1+3.e0*advij-advim1)

                    b=dxt22(i)*(advip1-2.e0*advij+advim1)

                    c=dxt16(i)                                          &
     &                *(6.e0*advip1-advip2-2.e0*advim1-3.e0*advij)

                    vf(i,j,k)=advij+((a*u8v+b)*u8v+c)*u8v

                  end if

                end if

              end do
              end do

!$omp end do

            end do

          end if

! -----

        end if

!! -----

!!! -----

!!! Perform x-y-z seperated Cubic Lagrange scheme.

      else if(advopt.eq.5) then

! Perform calculation without the map scale factor.

        if(mfcopt.eq.0) then

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,a,b,c)

            do j=2,nj-1
            do i=0,ni

              if(vp(i,j,k).gt.0.e0) then

                a=dyt36s(j)*(vp(i,j+1,k)-3.e0*vp(i,j,k)                 &
     &            +3.e0*vp(i,j-1,k)-vp(i,j-2,k))

                b=dyt22(j)*(vp(i,j+1,k)                                 &
     &            -2.e0*vp(i,j,k)+vp(i,j-1,k))

                c=dyt16(j)*(3.e0*vp(i,j,k)+2.e0*vp(i,j+1,k)             &
     &            +vp(i,j-2,k)-6.e0*vp(i,j-1,k))

                advd(i,j,k)=vp(i,j,k)                                   &
     &            +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

              else

                a=dyt36n(j)*(vp(i,j+2,k)-3.e0*vp(i,j+1,k)               &
     &            +3.e0*vp(i,j,k)-vp(i,j-1,k))

                b=dyt22(j)*(vp(i,j+1,k)                                 &
     &            -2.e0*vp(i,j,k)+vp(i,j-1,k))

                c=dyt16(j)*(6.e0*vp(i,j+1,k)-vp(i,j+2,k)                &
     &            -2.e0*vp(i,j-1,k)-3.e0*vp(i,j,k))

                advd(i,j,k)=vp(i,j,k)                                   &
     &            +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

              end if

            end do
            end do

!$omp end do

          end do

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,u8v,a,b,c)

            do j=2,nj-1
            do i=2,ni-2

              u8v=.25e0                                                 &
     &          *((up(i,j-1,k)+up(i,j,k))+(up(i+1,j-1,k)+up(i+1,j,k)))

              if(u8v.gt.0.e0) then

                a=dxt36w(i)*(advd(i+1,j,k)-3.e0*advd(i,j,k)             &
     &            +3.e0*advd(i-1,j,k)-advd(i-2,j,k))

                b=dxt22(i)*(advd(i+1,j,k)                               &
     &            -2.e0*advd(i,j,k)+advd(i-1,j,k))

                c=dxt16(i)*(3.e0*advd(i,j,k)+2.e0*advd(i+1,j,k)         &
     &            +advd(i-2,j,k)-6.e0*advd(i-1,j,k))

                vf(i,j,k)=advd(i,j,k)+((a*u8v+b)*u8v+c)*u8v

              else

                a=dxt36e(i)*(advd(i+2,j,k)-3.e0*advd(i+1,j,k)           &
     &            +3.e0*advd(i,j,k)-advd(i-1,j,k))

                b=dxt22(i)*(advd(i+1,j,k)                               &
     &            -2.e0*advd(i,j,k)+advd(i-1,j,k))

                c=dxt16(i)*(6.e0*advd(i+1,j,k)-advd(i+2,j,k)            &
     &            -2.e0*advd(i-1,j,k)-3.e0*advd(i,j,k))

                vf(i,j,k)=advd(i,j,k)+((a*u8v+b)*u8v+c)*u8v

              end if

            end do
            end do

!$omp end do

          end do

! -----

!! Perform calculation with the map scale factor.

        else

! Applied the map scale factor for the x direction.

          if(mpopt.eq.0.or.mpopt.eq.10) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,a,b,c)

              do j=2,nj-1
              do i=0,ni

                if(vp(i,j,k).gt.0.e0) then

                  a=dyt36s(j)*(vp(i,j+1,k)-3.e0*vp(i,j,k)               &
     &              +3.e0*vp(i,j-1,k)-vp(i,j-2,k))

                  b=dyt22(j)*(vp(i,j+1,k)                               &
     &              -2.e0*vp(i,j,k)+vp(i,j-1,k))

                  c=dyt16(j)*(3.e0*vp(i,j,k)+2.e0*vp(i,j+1,k)           &
     &              +vp(i,j-2,k)-6.e0*vp(i,j-1,k))

                  advd(i,j,k)=vp(i,j,k)                                 &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                else

                  a=dyt36n(j)*(vp(i,j+2,k)-3.e0*vp(i,j+1,k)             &
     &              +3.e0*vp(i,j,k)-vp(i,j-1,k))

                  b=dyt22(j)*(vp(i,j+1,k)                               &
     &              -2.e0*vp(i,j,k)+vp(i,j-1,k))

                  c=dyt16(j)*(6.e0*vp(i,j+1,k)-vp(i,j+2,k)              &
     &              -2.e0*vp(i,j-1,k)-3.e0*vp(i,j,k))

                  advd(i,j,k)=vp(i,j,k)                                 &
     &              +((a*vp(i,j,k)+b)*vp(i,j,k)+c)*vp(i,j,k)

                end if

              end do
              end do

!$omp end do

            end do

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,u8v,a,b,c)

              do j=2,nj-1
              do i=2,ni-2

                u8v=.25e0*((mf8u(i,j-1)*up(i,j-1,k)+mf8u(i,j)*up(i,j,k))&
     &           +(mf8u(i+1,j-1)*up(i+1,j-1,k)+mf8u(i+1,j)*up(i+1,j,k)))

                if(u8v.gt.0.e0) then

                  a=dxt36w(i)*(advd(i+1,j,k)-3.e0*advd(i,j,k)           &
     &              +3.e0*advd(i-1,j,k)-advd(i-1,j,k))

                  b=dxt22(i)*(advd(i+1,j,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i-1,j,k))

                  c=dxt16(i)*(3.e0*advd(i,j,k)+2.e0*advd(i+1,j,k)       &
     &              +advd(i-1,j,k)-6.e0*advd(i-1,j,k))

                  vf(i,j,k)=advd(i,j,k)+((a*u8v+b)*u8v+c)*u8v

                else

                  a=dxt36e(i)*(advd(i+2,j,k)-3.e0*advd(i+1,j,k)         &
     &              +3.e0*advd(i,j,k)-advd(i-1,j,k))

                  b=dxt22(i)*(advd(i+1,j,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i-1,j,k))

                  c=dxt16(i)*(6.e0*advd(i+1,j,k)-advd(i+2,j,k)          &
     &              -2.e0*advd(i-1,j,k)-3.e0*advd(i,j,k))

                  vf(i,j,k)=advd(i,j,k)+((a*u8v+b)*u8v+c)*u8v

                end if

              end do
              end do

!$omp end do

            end do

! -----

! Applied the map scale factor for the y direction.

          else if(mpopt.eq.5) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,vm,a,b,c)

              do j=2,nj-1
              do i=0,ni

                vm=mf8v(i,j)*vp(i,j,k)

                if(vm.gt.0.e0) then

                  a=dyt36s(j)*(vp(i,j+1,k)-3.e0*vp(i,j,k)               &
     &              +3.e0*vp(i,j-1,k)-vp(i,j-2,k))

                  b=dyt22(j)*(vp(i,j+1,k)                               &
     &              -2.e0*vp(i,j,k)+vp(i,j-1,k))

                  c=dyt16(j)*(3.e0*vp(i,j,k)+2.e0*vp(i,j+1,k)           &
     &              +vp(i,j-2,k)-6.e0*vp(i,j-1,k))

                  advd(i,j,k)=vp(i,j,k)+((a*vm+b)*vm+c)*vm

                else

                  a=dyt36n(j)*(vp(i,j+2,k)-3.e0*vp(i,j+1,k)             &
     &              +3.e0*vp(i,j,k)-vp(i,j-1,k))

                  b=dyt22(j)*(vp(i,j+1,k)                               &
     &              -2.e0*vp(i,j,k)+vp(i,j-1,k))

                  c=dyt16(j)*(6.e0*vp(i,j+1,k)-vp(i,j+2,k)              &
     &              -2.e0*vp(i,j-1,k)-3.e0*vp(i,j,k))

                  advd(i,j,k)=vp(i,j,k)+((a*vm+b)*vm+c)*vm

                end if

              end do
              end do

!$omp end do

            end do

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,u8v,a,b,c)

              do j=2,nj-1
              do i=2,ni-2

                u8v=.25e0                                               &
     &            *((up(i,j-1,k)+up(i,j,k))+(up(i+1,j-1,k)+up(i+1,j,k)))

                if(u8v.gt.0.e0) then

                  a=dxt36w(i)*(advd(i+1,j,k)-3.e0*advd(i,j,k)           &
     &              +3.e0*advd(i-1,j,k)-advd(i-2,j,k))

                  b=dxt22(i)*(advd(i+1,j,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i-1,j,k))

                  c=dxt16(i)*(3.e0*advd(i,j,k)+2.e0*advd(i+1,j,k)       &
     &              +advd(i-2,j,k)-6.e0*advd(i-1,j,k))

                  vf(i,j,k)=advd(i,j,k)+((a*u8v+b)*u8v+c)*u8v

                else

                  a=dxt36e(i)*(advd(i+2,j,k)-3.e0*advd(i+1,j,k)         &
     &              +3.e0*advd(i,j,k)-advd(i-1,j,k))

                  b=dxt22(i)*(advd(i+1,j,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i-1,j,k))

                  c=dxt16(i)*(6.e0*advd(i+1,j,k)-advd(i+2,j,k)          &
     &              -2.e0*advd(i-1,j,k)-3.e0*advd(i,j,k))

                  vf(i,j,k)=advd(i,j,k)+((a*u8v+b)*u8v+c)*u8v

                end if

              end do
              end do

!$omp end do

            end do

! -----

! Applied the map scale factor for the x and y direction.

          else

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,vm,a,b,c)

              do j=2,nj-1
              do i=0,ni

                vm=mf8v(i,j)*vp(i,j,k)

                if(vm.gt.0.e0) then

                  a=dyt36s(j)*(vp(i,j+1,k)-3.e0*vp(i,j,k)               &
     &              +3.e0*vp(i,j-1,k)-vp(i,j-2,k))

                  b=dyt22(j)*(vp(i,j+1,k)                               &
     &              -2.e0*vp(i,j,k)+vp(i,j-1,k))

                  c=dyt16(j)*(3.e0*vp(i,j,k)+2.e0*vp(i,j+1,k)           &
     &              +vp(i,j-2,k)-6.e0*vp(i,j-1,k))

                  advd(i,j,k)=vp(i,j,k)+((a*vm+b)*vm+c)*vm

                else

                  a=dyt36n(j)*(vp(i,j+2,k)-3.e0*vp(i,j+1,k)             &
     &              +3.e0*vp(i,j,k)-vp(i,j-1,k))

                  b=dyt22(j)*(vp(i,j+1,k)                               &
     &              -2.e0*vp(i,j,k)+vp(i,j-1,k))

                  c=dyt16(j)*(6.e0*vp(i,j+1,k)-vp(i,j+2,k)              &
     &              -2.e0*vp(i,j-1,k)-3.e0*vp(i,j,k))

                  advd(i,j,k)=vp(i,j,k)+((a*vm+b)*vm+c)*vm

                end if

              end do
              end do

!$omp end do

            end do

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,u8v,a,b,c)

              do j=2,nj-1
              do i=2,ni-2

                u8v=.25e0*((mf8u(i,j-1)*up(i,j-1,k)+mf8u(i,j)*up(i,j,k))&
     &           +(mf8u(i+1,j-1)*up(i+1,j-1,k)+mf8u(i+1,j)*up(i+1,j,k)))

                if(u8v.gt.0.e0) then

                  a=dxt36w(i)*(advd(i+1,j,k)-3.e0*advd(i,j,k)           &
     &              +3.e0*advd(i-1,j,k)-advd(i-2,j,k))

                  b=dxt22(i)*(advd(i+1,j,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i-1,j,k))

                  c=dxt16(i)*(3.e0*advd(i,j,k)+2.e0*advd(i+1,j,k)       &
     &              +advd(i-2,j,k)-6.e0*advd(i-1,j,k))

                  vf(i,j,k)=advd(i,j,k)+((a*u8v+b)*u8v+c)*u8v

                else

                  a=dxt36e(i)*(advd(i+2,j,k)-3.e0*advd(i+1,j,k)         &
     &              +3.e0*advd(i,j,k)-advd(i-1,j,k))

                  b=dxt22(i)*(advd(i+1,j,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i-1,j,k))

                  c=dxt16(i)*(6.e0*advd(i+1,j,k)-advd(i+2,j,k)          &
     &              -2.e0*advd(i-1,j,k)-3.e0*advd(i,j,k))

                  vf(i,j,k)=advd(i,j,k)+((a*u8v+b)*u8v+c)*u8v

                end if

              end do
              end do

!$omp end do

            end do

          end if

! -----

        end if

      end if

!!! -----

!$omp end parallel

!!!! -----

!!!!! -----

!!!!! Calculate the w advection horizontally.

! Set the common used variables.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(i)

      do i=0,ni
        dxt22(i)=cdxt2
        dxt16(i)=cdxt1

        dxt36w(i)=cdxt3
        dxt36e(i)=cdxt3

      end do

!$omp end do

!$omp do schedule(runtime) private(j)

      do j=0,nj
        dyt22(j)=cdyt2
        dyt16(j)=cdyt1

        dyt36s(j)=cdyt3
        dyt36n(j)=cdyt3

      end do

!$omp end do

!$omp end parallel

! -----

! Set the lateral boundary conditions.

      call lbculw(idwbc,idebc,idsbc,idnbc,ni,nj,nk,wp,                  &
     &            dxt36w,dxt36e,dyt36s,dyt36n)

! -----

!!!! Perform Cubic Lagrange scheme.

!$omp parallel default(shared) private(k)

!!! Perform horizontal-vertical seperated Cubic Lagrange scheme.

      if(advopt.eq.4) then

! Perform calculation without the map scale factor.

        if(mfcopt.eq.0) then

          do k=2,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,u8w,v8w,advij,advjm2,advjm1,advjp1,advjp2,a,b,c)

            do j=2,nj-2
            do i=2,ni-2

              u8w=.25e0                                                 &
     &          *((up(i,j,k-1)+up(i,j,k))+(up(i+1,j,k-1)+up(i+1,j,k)))

              v8w=.25e0                                                 &
     &          *((vp(i,j,k-1)+vp(i,j,k))+(vp(i,j+1,k-1)+vp(i,j+1,k)))

              if(u8w.gt.0.e0) then

                if(v8w.gt.0.e0) then

                  a=dxt36w(i)*(wp(i+1,j-2,k)-3.e0*wp(i,j-2,k)           &
     &              +3.e0*wp(i-1,j-2,k)-wp(i-2,j-2,k))

                  b=dxt22(i)*(wp(i+1,j-2,k)                             &
     &              -2.e0*wp(i,j-2,k)+wp(i-1,j-2,k))

                  c=dxt16(i)*(3.e0*wp(i,j-2,k)+2.e0*wp(i+1,j-2,k)       &
     &              +wp(i-2,j-2,k)-6.e0*wp(i-1,j-2,k))

                  advjm2=wp(i,j-2,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dxt36w(i)*(wp(i+1,j-1,k)-3.e0*wp(i,j-1,k)           &
     &              +3.e0*wp(i-1,j-1,k)-wp(i-2,j-1,k))

                  b=dxt22(i)*(wp(i+1,j-1,k)                             &
     &              -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                  c=dxt16(i)*(3.e0*wp(i,j-1,k)+2.e0*wp(i+1,j-1,k)       &
     &              +wp(i-2,j-1,k)-6.e0*wp(i-1,j-1,k))

                  advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dxt36w(i)*(wp(i+1,j,k)-3.e0*wp(i,j,k)               &
     &              +3.e0*wp(i-1,j,k)-wp(i-2,j,k))

                  b=dxt22(i)*(wp(i+1,j,k)                               &
     &              -2.e0*wp(i,j,k)+wp(i-1,j,k))

                  c=dxt16(i)*(3.e0*wp(i,j,k)+2.e0*wp(i+1,j,k)           &
     &              +wp(i-2,j,k)-6.e0*wp(i-1,j,k))

                  advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dxt36w(i)*(wp(i+1,j+1,k)-3.e0*wp(i,j+1,k)           &
     &              +3.e0*wp(i-1,j+1,k)-wp(i-2,j+1,k))

                  b=dxt22(i)*(wp(i+1,j+1,k)                             &
     &              -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                  c=dxt16(i)*(3.e0*wp(i,j+1,k)+2.e0*wp(i+1,j+1,k)       &
     &              +wp(i-2,j+1,k)-6.e0*wp(i-1,j+1,k))

                  advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                  b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                  c=dyt16(j)                                            &
     &              *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                  wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                else

                  a=dxt36w(i)*(wp(i+1,j-1,k)-3.e0*wp(i,j-1,k)           &
     &              +3.e0*wp(i-1,j-1,k)-wp(i-2,j-1,k))

                  b=dxt22(i)*(wp(i+1,j-1,k)                             &
     &              -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                  c=dxt16(i)*(3.e0*wp(i,j-1,k)+2.e0*wp(i+1,j-1,k)       &
     &              +wp(i-2,j-1,k)-6.e0*wp(i-1,j-1,k))

                  advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dxt36w(i)*(wp(i+1,j,k)-3.e0*wp(i,j,k)               &
     &              +3.e0*wp(i-1,j,k)-wp(i-2,j,k))

                  b=dxt22(i)*(wp(i+1,j,k)                               &
     &              -2.e0*wp(i,j,k)+wp(i-1,j,k))

                  c=dxt16(i)*(3.e0*wp(i,j,k)+2.e0*wp(i+1,j,k)           &
     &              +wp(i-2,j,k)-6.e0*wp(i-1,j,k))

                  advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dxt36w(i)*(wp(i+1,j+1,k)-3.e0*wp(i,j+1,k)           &
     &              +3.e0*wp(i-1,j+1,k)-wp(i-2,j+1,k))

                  b=dxt22(i)*(wp(i+1,j+1,k)                             &
     &              -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                  c=dxt16(i)*(3.e0*wp(i,j+1,k)+2.e0*wp(i+1,j+1,k)       &
     &              +wp(i-2,j+1,k)-6.e0*wp(i-1,j+1,k))

                  advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dxt36w(i)*(wp(i+1,j+2,k)-3.e0*wp(i,j+2,k)           &
     &              +3.e0*wp(i-1,j+2,k)-wp(i-2,j+2,k))

                  b=dxt22(i)*(wp(i+1,j+2,k)                             &
     &              -2.e0*wp(i,j+2,k)+wp(i-1,j+2,k))

                  c=dxt16(i)*(3.e0*wp(i,j+2,k)+2.e0*wp(i+1,j+2,k)       &
     &              +wp(i-2,j+2,k)-6.e0*wp(i-1,j+2,k))

                  advjp2=wp(i,j+2,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                  b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                  c=dyt16(j)                                            &
     &              *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                  wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                end if

              else

                if(v8w.gt.0.e0) then

                  a=dxt36e(i)*(wp(i+2,j-2,k)-3.e0*wp(i+1,j-2,k)         &
     &              +3.e0*wp(i,j-2,k)-wp(i-1,j-2,k))

                  b=dxt22(i)*(wp(i+1,j-2,k)                             &
     &              -2.e0*wp(i,j-2,k)+wp(i-1,j-2,k))

                  c=dxt16(i)*(6.e0*wp(i+1,j-2,k)-wp(i+2,j-2,k)          &
     &              -2.e0*wp(i-1,j-2,k)-3.e0*wp(i,j-2,k))

                  advjm2=wp(i,j-2,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dxt36e(i)*(wp(i+2,j-1,k)-3.e0*wp(i+1,j-1,k)         &
     &              +3.e0*wp(i,j-1,k)-wp(i-1,j-1,k))

                  b=dxt22(i)*(wp(i+1,j-1,k)                             &
     &              -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                  c=dxt16(i)*(6.e0*wp(i+1,j-1,k)-wp(i+2,j-1,k)          &
     &              -2.e0*wp(i-1,j-1,k)-3.e0*wp(i,j-1,k))

                  advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dxt36e(i)*(wp(i+2,j,k)-3.e0*wp(i+1,j,k)             &
     &              +3.e0*wp(i,j,k)-wp(i-1,j,k))

                  b=dxt22(i)*(wp(i+1,j,k)                               &
     &              -2.e0*wp(i,j,k)+wp(i-1,j,k))

                  c=dxt16(i)*(6.e0*wp(i+1,j,k)-wp(i+2,j,k)              &
     &              -2.e0*wp(i-1,j,k)-3.e0*wp(i,j,k))

                  advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dxt36e(i)*(wp(i+2,j+1,k)-3.e0*wp(i+1,j+1,k)         &
     &              +3.e0*wp(i,j+1,k)-wp(i-1,j+1,k))

                  b=dxt22(i)*(wp(i+1,j+1,k)                             &
     &              -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                  c=dxt16(i)*(6.e0*wp(i+1,j+1,k)-wp(i+2,j+1,k)          &
     &              -2.e0*wp(i-1,j+1,k)-3.e0*wp(i,j+1,k))

                  advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                  b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                  c=dyt16(j)                                            &
     &              *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                  wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                else

                  a=dxt36e(i)*(wp(i+2,j-1,k)-3.e0*wp(i+1,j-1,k)         &
     &              +3.e0*wp(i,j-1,k)-wp(i-1,j-1,k))

                  b=dxt22(i)*(wp(i+1,j-1,k)                             &
     &              -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                  c=dxt16(i)*(6.e0*wp(i+1,j-1,k)-wp(i+2,j-1,k)          &
     &              -2.e0*wp(i-1,j-1,k)-3.e0*wp(i,j-1,k))

                  advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dxt36e(i)*(wp(i+2,j,k)-3.e0*wp(i+1,j,k)             &
     &              +3.e0*wp(i,j,k)-wp(i-1,j,k))

                  b=dxt22(i)*(wp(i+1,j,k)                               &
     &              -2.e0*wp(i,j,k)+wp(i-1,j,k))

                  c=dxt16(i)*(6.e0*wp(i+1,j,k)-wp(i+2,j,k)              &
     &              -2.e0*wp(i-1,j,k)-3.e0*wp(i,j,k))

                  advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dxt36e(i)*(wp(i+2,j+1,k)-3.e0*wp(i+1,j+1,k)         &
     &              +3.e0*wp(i,j+1,k)-wp(i-1,j+1,k))

                  b=dxt22(i)*(wp(i+1,j+1,k)                             &
     &              -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                  c=dxt16(i)*(6.e0*wp(i+1,j+1,k)-wp(i+2,j+1,k)          &
     &              -2.e0*wp(i-1,j+1,k)-3.e0*wp(i,j+1,k))

                  advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dxt36e(i)*(wp(i+2,j+2,k)-3.e0*wp(i+1,j+2,k)         &
     &              +3.e0*wp(i,j+2,k)-wp(i-1,j+2,k))

                  b=dxt22(i)*(wp(i+1,j+2,k)                             &
     &              -2.e0*wp(i,j+2,k)+wp(i-1,j+2,k))

                  c=dxt16(i)*(6.e0*wp(i+1,j+2,k)-wp(i+2,j+2,k)          &
     &              -2.e0*wp(i-1,j+2,k)-3.e0*wp(i,j+2,k))

                  advjp2=wp(i,j+2,k)+((a*u8w+b)*u8w+c)*u8w

                  a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                  b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                  c=dyt16(j)                                            &
     &              *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                  wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                end if

              end if

            end do
            end do

!$omp end do

          end do

! -----

!! Perform calculation with the map scale factor.

        else

! Applied the map scale factor for the x direction.

          if(mpopt.eq.0.or.mpopt.eq.10) then

            do k=2,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,u8w,v8w,advij,advjm2,advjm1,advjp1,advjp2,a,b,c)

              do j=2,nj-2
              do i=2,ni-2

                u8w=.25e0*(mf8u(i,j)*(up(i,j,k-1)+up(i,j,k))            &
     &            +mf8u(i+1,j)*(up(i+1,j,k-1)+up(i+1,j,k)))

                v8w=.25e0                                               &
     &            *((vp(i,j,k-1)+vp(i,j,k))+(vp(i,j+1,k-1)+vp(i,j+1,k)))

                if(u8w.gt.0.e0) then

                  if(v8w.gt.0.e0) then

                    a=dxt36w(i)*(wp(i+1,j-2,k)-3.e0*wp(i,j-2,k)         &
     &                +3.e0*wp(i-1,j-2,k)-wp(i-2,j-2,k))

                    b=dxt22(i)*(wp(i+1,j-2,k)                           &
     &                -2.e0*wp(i,j-2,k)+wp(i-1,j-2,k))

                    c=dxt16(i)*(3.e0*wp(i,j-2,k)+2.e0*wp(i+1,j-2,k)     &
     &                +wp(i-2,j-2,k)-6.e0*wp(i-1,j-2,k))

                    advjm2=wp(i,j-2,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j-1,k)-3.e0*wp(i,j-1,k)         &
     &                +3.e0*wp(i-1,j-1,k)-wp(i-2,j-1,k))

                    b=dxt22(i)*(wp(i+1,j-1,k)                           &
     &                -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*wp(i,j-1,k)+2.e0*wp(i+1,j-1,k)     &
     &                +wp(i-2,j-1,k)-6.e0*wp(i-1,j-1,k))

                    advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j,k)-3.e0*wp(i,j,k)             &
     &                +3.e0*wp(i-1,j,k)-wp(i-2,j,k))

                    b=dxt22(i)*(wp(i+1,j,k)                             &
     &                -2.e0*wp(i,j,k)+wp(i-1,j,k))

                    c=dxt16(i)*(3.e0*wp(i,j,k)+2.e0*wp(i+1,j,k)         &
     &                +wp(i-2,j,k)-6.e0*wp(i-1,j,k))

                    advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j+1,k)-3.e0*wp(i,j+1,k)         &
     &                +3.e0*wp(i-1,j+1,k)-wp(i-2,j+1,k))

                    b=dxt22(i)*(wp(i+1,j+1,k)                           &
     &                -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*wp(i,j+1,k)+2.e0*wp(i+1,j+1,k)     &
     &                +wp(i-2,j+1,k)-6.e0*wp(i-1,j+1,k))

                    advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                  else

                    a=dxt36w(i)*(wp(i+1,j-1,k)-3.e0*wp(i,j-1,k)         &
     &                +3.e0*wp(i-1,j-1,k)-wp(i-2,j-1,k))

                    b=dxt22(i)*(wp(i+1,j-1,k)                           &
     &                -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*wp(i,j-1,k)+2.e0*wp(i+1,j-1,k)     &
     &                +wp(i-2,j-1,k)-6.e0*wp(i-1,j-1,k))

                    advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j,k)-3.e0*wp(i,j,k)             &
     &                +3.e0*wp(i-1,j,k)-wp(i-2,j,k))

                    b=dxt22(i)*(wp(i+1,j,k)                             &
     &                -2.e0*wp(i,j,k)+wp(i-1,j,k))

                    c=dxt16(i)*(3.e0*wp(i,j,k)+2.e0*wp(i+1,j,k)         &
     &                +wp(i-2,j,k)-6.e0*wp(i-1,j,k))

                    advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j+1,k)-3.e0*wp(i,j+1,k)         &
     &                +3.e0*wp(i-1,j+1,k)-wp(i-2,j+1,k))

                    b=dxt22(i)*(wp(i+1,j+1,k)                           &
     &                -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*wp(i,j+1,k)+2.e0*wp(i+1,j+1,k)     &
     &                +wp(i-2,j+1,k)-6.e0*wp(i-1,j+1,k))

                    advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j+2,k)-3.e0*wp(i,j+2,k)         &
     &                +3.e0*wp(i-1,j+2,k)-wp(i-2,j+2,k))

                    b=dxt22(i)*(wp(i+1,j+2,k)                           &
     &                -2.e0*wp(i,j+2,k)+wp(i-1,j+2,k))

                    c=dxt16(i)*(3.e0*wp(i,j+2,k)+2.e0*wp(i+1,j+2,k)     &
     &                +wp(i-2,j+2,k)-6.e0*wp(i-1,j+2,k))

                    advjp2=wp(i,j+2,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                  end if

                else

                  if(v8w.gt.0.e0) then

                    a=dxt36e(i)*(wp(i+2,j-2,k)-3.e0*wp(i+1,j-2,k)       &
     &                +3.e0*wp(i,j-2,k)-wp(i-1,j-2,k))

                    b=dxt22(i)*(wp(i+1,j-2,k)                           &
     &                -2.e0*wp(i,j-2,k)+wp(i-1,j-2,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j-2,k)-wp(i+2,j-2,k)        &
     &                -2.e0*wp(i-1,j-2,k)-3.e0*wp(i,j-2,k))

                    advjm2=wp(i,j-2,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j-1,k)-3.e0*wp(i+1,j-1,k)       &
     &                +3.e0*wp(i,j-1,k)-wp(i-1,j-1,k))

                    b=dxt22(i)*(wp(i+1,j-1,k)                           &
     &                -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j-1,k)-wp(i+2,j-1,k)        &
     &                -2.e0*wp(i-1,j-1,k)-3.e0*wp(i,j-1,k))

                    advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j,k)-3.e0*wp(i+1,j,k)           &
     &                +3.e0*wp(i,j,k)-wp(i-1,j,k))

                    b=dxt22(i)*(wp(i+1,j,k)                             &
     &                -2.e0*wp(i,j,k)+wp(i-1,j,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j,k)-wp(i+2,j,k)            &
     &                -2.e0*wp(i-1,j,k)-3.e0*wp(i,j,k))

                    advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j+1,k)-3.e0*wp(i+1,j+1,k)       &
     &                +3.e0*wp(i,j+1,k)-wp(i-1,j+1,k))

                    b=dxt22(i)*(wp(i+1,j+1,k)                           &
     &                -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j+1,k)-wp(i+2,j+1,k)        &
     &                -2.e0*wp(i-1,j+1,k)-3.e0*wp(i,j+1,k))

                    advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                  else

                    a=dxt36e(i)*(wp(i+2,j-1,k)-3.e0*wp(i+1,j-1,k)       &
     &                +3.e0*wp(i,j-1,k)-wp(i-1,j-1,k))

                    b=dxt22(i)*(wp(i+1,j-1,k)                           &
     &                -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j-1,k)-wp(i+2,j-1,k)        &
     &                -2.e0*wp(i-1,j-1,k)-3.e0*wp(i,j-1,k))

                    advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j,k)-3.e0*wp(i+1,j,k)           &
     &                +3.e0*wp(i,j,k)-wp(i-1,j,k))

                    b=dxt22(i)*(wp(i+1,j,k)                             &
     &                -2.e0*wp(i,j,k)+wp(i-1,j,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j,k)-wp(i+2,j,k)            &
     &                -2.e0*wp(i-1,j,k)-3.e0*wp(i,j,k))

                    advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j+1,k)-3.e0*wp(i+1,j+1,k)       &
     &                +3.e0*wp(i,j+1,k)-wp(i-1,j+1,k))

                    b=dxt22(i)*(wp(i+1,j+1,k)                           &
     &                -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j+1,k)-wp(i+2,j+1,k)        &
     &                -2.e0*wp(i-1,j+1,k)-3.e0*wp(i,j+1,k))

                    advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j+2,k)-3.e0*wp(i+1,j+2,k)       &
     &                +3.e0*wp(i,j+2,k)-wp(i-1,j+2,k))

                    b=dxt22(i)*(wp(i+1,j+2,k)                           &
     &                -2.e0*wp(i,j+2,k)+wp(i-1,j+2,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j+2,k)-wp(i+2,j+2,k)        &
     &                -2.e0*wp(i-1,j+2,k)-3.e0*wp(i,j+2,k))

                    advjp2=wp(i,j+2,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                  end if

                end if

              end do
              end do

!$omp end do

            end do

! -----

! Applied the map scale factor for the y direction.

          else if(mpopt.eq.5) then

            do k=2,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,u8w,v8w,advij,advjm2,advjm1,advjp1,advjp2,a,b,c)

              do j=2,nj-2
              do i=2,ni-2

                u8w=.25e0                                               &
     &            *((up(i,j,k-1)+up(i,j,k))+(up(i+1,j,k-1)+up(i+1,j,k)))

                v8w=.25e0*(mf8v(i,j)*(vp(i,j,k-1)+vp(i,j,k))            &
     &            +mf8v(i,j+1)*(vp(i,j+1,k-1)+vp(i,j+1,k)))

                if(u8w.gt.0.e0) then

                  if(v8w.gt.0.e0) then

                    a=dxt36w(i)*(wp(i+1,j-2,k)-3.e0*wp(i,j-2,k)         &
     &                +3.e0*wp(i-1,j-2,k)-wp(i-2,j-2,k))

                    b=dxt22(i)*(wp(i+1,j-2,k)                           &
     &                -2.e0*wp(i,j-2,k)+wp(i-1,j-2,k))

                    c=dxt16(i)*(3.e0*wp(i,j-2,k)+2.e0*wp(i+1,j-2,k)     &
     &                +wp(i-2,j-2,k)-6.e0*wp(i-1,j-2,k))

                    advjm2=wp(i,j-2,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j-1,k)-3.e0*wp(i,j-1,k)         &
     &                +3.e0*wp(i-1,j-1,k)-wp(i-2,j-1,k))

                    b=dxt22(i)*(wp(i+1,j-1,k)                           &
     &                -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*wp(i,j-1,k)+2.e0*wp(i+1,j-1,k)     &
     &                +wp(i-2,j-1,k)-6.e0*wp(i-1,j-1,k))

                    advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j,k)-3.e0*wp(i,j,k)             &
     &                +3.e0*wp(i-1,j,k)-wp(i-2,j,k))

                    b=dxt22(i)*(wp(i+1,j,k)                             &
     &                -2.e0*wp(i,j,k)+wp(i-1,j,k))

                    c=dxt16(i)*(3.e0*wp(i,j,k)+2.e0*wp(i+1,j,k)         &
     &                +wp(i-2,j,k)-6.e0*wp(i-1,j,k))

                    advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j+1,k)-3.e0*wp(i,j+1,k)         &
     &                +3.e0*wp(i-1,j+1,k)-wp(i-2,j+1,k))

                    b=dxt22(i)*(wp(i+1,j+1,k)                           &
     &                -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*wp(i,j+1,k)+2.e0*wp(i+1,j+1,k)     &
     &                +wp(i-2,j+1,k)-6.e0*wp(i-1,j+1,k))

                    advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                  else

                    a=dxt36w(i)*(wp(i+1,j-1,k)-3.e0*wp(i,j-1,k)         &
     &                +3.e0*wp(i-1,j-1,k)-wp(i-2,j-1,k))

                    b=dxt22(i)*(wp(i+1,j-1,k)                           &
     &                -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*wp(i,j-1,k)+2.e0*wp(i+1,j-1,k)     &
     &                +wp(i-2,j-1,k)-6.e0*wp(i-1,j-1,k))

                    advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j,k)-3.e0*wp(i,j,k)             &
     &                +3.e0*wp(i-1,j,k)-wp(i-2,j,k))

                    b=dxt22(i)*(wp(i+1,j,k)                             &
     &                -2.e0*wp(i,j,k)+wp(i-1,j,k))

                    c=dxt16(i)*(3.e0*wp(i,j,k)+2.e0*wp(i+1,j,k)         &
     &                +wp(i-2,j,k)-6.e0*wp(i-1,j,k))

                    advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j+1,k)-3.e0*wp(i,j+1,k)         &
     &                +3.e0*wp(i-1,j+1,k)-wp(i-2,j+1,k))

                    b=dxt22(i)*(wp(i+1,j+1,k)                           &
     &                -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*wp(i,j+1,k)+2.e0*wp(i+1,j+1,k)     &
     &                +wp(i-2,j+1,k)-6.e0*wp(i-1,j+1,k))

                    advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j+2,k)-3.e0*wp(i,j+2,k)         &
     &                +3.e0*wp(i-1,j+2,k)-wp(i-2,j+2,k))

                    b=dxt22(i)*(wp(i+1,j+2,k)                           &
     &                -2.e0*wp(i,j+2,k)+wp(i-1,j+2,k))

                    c=dxt16(i)*(3.e0*wp(i,j+2,k)+2.e0*wp(i+1,j+2,k)     &
     &                +wp(i-2,j+2,k)-6.e0*wp(i-1,j+2,k))

                    advjp2=wp(i,j+2,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                  end if

                else

                  if(v8w.gt.0.e0) then

                    a=dxt36e(i)*(wp(i+2,j-2,k)-3.e0*wp(i+1,j-2,k)       &
     &                +3.e0*wp(i,j-2,k)-wp(i-1,j-2,k))

                    b=dxt22(i)*(wp(i+1,j-2,k)                           &
     &                -2.e0*wp(i,j-2,k)+wp(i-1,j-2,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j-2,k)-wp(i+2,j-2,k)        &
     &                -2.e0*wp(i-1,j-2,k)-3.e0*wp(i,j-2,k))

                    advjm2=wp(i,j-2,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j-1,k)-3.e0*wp(i+1,j-1,k)       &
     &                +3.e0*wp(i,j-1,k)-wp(i-1,j-1,k))

                    b=dxt22(i)*(wp(i+1,j-1,k)                           &
     &                -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j-1,k)-wp(i+2,j-1,k)        &
     &                -2.e0*wp(i-1,j-1,k)-3.e0*wp(i,j-1,k))

                    advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j,k)-3.e0*wp(i+1,j,k)           &
     &                +3.e0*wp(i,j,k)-wp(i-1,j,k))

                    b=dxt22(i)*(wp(i+1,j,k)                             &
     &                -2.e0*wp(i,j,k)+wp(i-1,j,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j,k)-wp(i+2,j,k)            &
     &                -2.e0*wp(i-1,j,k)-3.e0*wp(i,j,k))

                    advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j+1,k)-3.e0*wp(i+1,j+1,k)       &
     &                +3.e0*wp(i,j+1,k)-wp(i-1,j+1,k))

                    b=dxt22(i)*(wp(i+1,j+1,k)                           &
     &                -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j+1,k)-wp(i+2,j+1,k)        &
     &                -2.e0*wp(i-1,j+1,k)-3.e0*wp(i,j+1,k))

                    advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                  else

                    a=dxt36e(i)*(wp(i+2,j-1,k)-3.e0*wp(i+1,j-1,k)       &
     &                +3.e0*wp(i,j-1,k)-wp(i-1,j-1,k))

                    b=dxt22(i)*(wp(i+1,j-1,k)                           &
     &                -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j-1,k)-wp(i+2,j-1,k)        &
     &                -2.e0*wp(i-1,j-1,k)-3.e0*wp(i,j-1,k))

                    advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j,k)-3.e0*wp(i+1,j,k)           &
     &                +3.e0*wp(i,j,k)-wp(i-1,j,k))

                    b=dxt22(i)*(wp(i+1,j,k)                             &
     &                -2.e0*wp(i,j,k)+wp(i-1,j,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j,k)-wp(i+2,j,k)            &
     &                -2.e0*wp(i-1,j,k)-3.e0*wp(i,j,k))

                    advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j+1,k)-3.e0*wp(i+1,j+1,k)       &
     &                +3.e0*wp(i,j+1,k)-wp(i-1,j+1,k))

                    b=dxt22(i)*(wp(i+1,j+1,k)                           &
     &                -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j+1,k)-wp(i+2,j+1,k)        &
     &                -2.e0*wp(i-1,j+1,k)-3.e0*wp(i,j+1,k))

                    advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j+2,k)-3.e0*wp(i+1,j+2,k)       &
     &                +3.e0*wp(i,j+2,k)-wp(i-1,j+2,k))

                    b=dxt22(i)*(wp(i+1,j+2,k)                           &
     &                -2.e0*wp(i,j+2,k)+wp(i-1,j+2,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j+2,k)-wp(i+2,j+2,k)        &
     &                -2.e0*wp(i-1,j+2,k)-3.e0*wp(i,j+2,k))

                    advjp2=wp(i,j+2,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                  end if

                end if

              end do
              end do

!$omp end do

            end do

! -----

! Applied the map scale factor for the x and y direction.

          else

            do k=2,nk-1

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,u8w,v8w,advij,advjm2,advjm1,advjp1,advjp2,a,b,c)

              do j=2,nj-2
              do i=2,ni-2

                u8w=.25e0*(mf8u(i,j)*(up(i,j,k-1)+up(i,j,k))            &
     &            +mf8u(i+1,j)*(up(i+1,j,k-1)+up(i+1,j,k)))

                v8w=.25e0*(mf8v(i,j)*(vp(i,j,k-1)+vp(i,j,k))            &
     &            +mf8v(i,j+1)*(vp(i,j+1,k-1)+vp(i,j+1,k)))

                if(u8w.gt.0.e0) then

                  if(v8w.gt.0.e0) then

                    a=dxt36w(i)*(wp(i+1,j-2,k)-3.e0*wp(i,j-2,k)         &
     &                +3.e0*wp(i-1,j-2,k)-wp(i-2,j-2,k))

                    b=dxt22(i)*(wp(i+1,j-2,k)                           &
     &                -2.e0*wp(i,j-2,k)+wp(i-1,j-2,k))

                    c=dxt16(i)*(3.e0*wp(i,j-2,k)+2.e0*wp(i+1,j-2,k)     &
     &                +wp(i-2,j-2,k)-6.e0*wp(i-1,j-2,k))

                    advjm2=wp(i,j-2,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j-1,k)-3.e0*wp(i,j-1,k)         &
     &                +3.e0*wp(i-1,j-1,k)-wp(i-2,j-1,k))

                    b=dxt22(i)*(wp(i+1,j-1,k)                           &
     &                -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*wp(i,j-1,k)+2.e0*wp(i+1,j-1,k)     &
     &                +wp(i-2,j-1,k)-6.e0*wp(i-1,j-1,k))

                    advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j,k)-3.e0*wp(i,j,k)             &
     &                +3.e0*wp(i-1,j,k)-wp(i-2,j,k))

                    b=dxt22(i)*(wp(i+1,j,k)                             &
     &                -2.e0*wp(i,j,k)+wp(i-1,j,k))

                    c=dxt16(i)*(3.e0*wp(i,j,k)+2.e0*wp(i+1,j,k)         &
     &                +wp(i-2,j,k)-6.e0*wp(i-1,j,k))

                    advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j+1,k)-3.e0*wp(i,j+1,k)         &
     &                +3.e0*wp(i-1,j+1,k)-wp(i-2,j+1,k))

                    b=dxt22(i)*(wp(i+1,j+1,k)                           &
     &                -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*wp(i,j+1,k)+2.e0*wp(i+1,j+1,k)     &
     &                +wp(i-2,j+1,k)-6.e0*wp(i-1,j+1,k))

                    advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                  else

                    a=dxt36w(i)*(wp(i+1,j-1,k)-3.e0*wp(i,j-1,k)         &
     &                +3.e0*wp(i-1,j-1,k)-wp(i-2,j-1,k))

                    b=dxt22(i)*(wp(i+1,j-1,k)                           &
     &                -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*wp(i,j-1,k)+2.e0*wp(i+1,j-1,k)     &
     &                +wp(i-2,j-1,k)-6.e0*wp(i-1,j-1,k))

                    advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j,k)-3.e0*wp(i,j,k)             &
     &                +3.e0*wp(i-1,j,k)-wp(i-2,j,k))

                    b=dxt22(i)*(wp(i+1,j,k)                             &
     &                -2.e0*wp(i,j,k)+wp(i-1,j,k))

                    c=dxt16(i)*(3.e0*wp(i,j,k)+2.e0*wp(i+1,j,k)         &
     &                +wp(i-2,j,k)-6.e0*wp(i-1,j,k))

                    advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j+1,k)-3.e0*wp(i,j+1,k)         &
     &                +3.e0*wp(i-1,j+1,k)-wp(i-2,j+1,k))

                    b=dxt22(i)*(wp(i+1,j+1,k)                           &
     &                -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*wp(i,j+1,k)+2.e0*wp(i+1,j+1,k)     &
     &                +wp(i-2,j+1,k)-6.e0*wp(i-1,j+1,k))

                    advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36w(i)*(wp(i+1,j+2,k)-3.e0*wp(i,j+2,k)         &
     &                +3.e0*wp(i-1,j+2,k)-wp(i-2,j+2,k))

                    b=dxt22(i)*(wp(i+1,j+2,k)                           &
     &                -2.e0*wp(i,j+2,k)+wp(i-1,j+2,k))

                    c=dxt16(i)*(3.e0*wp(i,j+2,k)+2.e0*wp(i+1,j+2,k)     &
     &                +wp(i-2,j+2,k)-6.e0*wp(i-1,j+2,k))

                    advjp2=wp(i,j+2,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                  end if

                else

                  if(v8w.gt.0.e0) then

                    a=dxt36e(i)*(wp(i+2,j-2,k)-3.e0*wp(i+1,j-2,k)       &
     &                +3.e0*wp(i,j-2,k)-wp(i-1,j-2,k))

                    b=dxt22(i)*(wp(i+1,j-2,k)                           &
     &                -2.e0*wp(i,j-2,k)+wp(i-1,j-2,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j-2,k)-wp(i+2,j-2,k)        &
     &                -2.e0*wp(i-1,j-2,k)-3.e0*wp(i,j-2,k))

                    advjm2=wp(i,j-2,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j-1,k)-3.e0*wp(i+1,j-1,k)       &
     &                +3.e0*wp(i,j-1,k)-wp(i-1,j-1,k))

                    b=dxt22(i)*(wp(i+1,j-1,k)                           &
     &                -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j-1,k)-wp(i+2,j-1,k)        &
     &                -2.e0*wp(i-1,j-1,k)-3.e0*wp(i,j-1,k))

                    advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j,k)-3.e0*wp(i+1,j,k)           &
     &                +3.e0*wp(i,j,k)-wp(i-1,j,k))

                    b=dxt22(i)*(wp(i+1,j,k)                             &
     &                -2.e0*wp(i,j,k)+wp(i-1,j,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j,k)-wp(i+2,j,k)            &
     &                -2.e0*wp(i-1,j,k)-3.e0*wp(i,j,k))

                    advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j+1,k)-3.e0*wp(i+1,j+1,k)       &
     &                +3.e0*wp(i,j+1,k)-wp(i-1,j+1,k))

                    b=dxt22(i)*(wp(i+1,j+1,k)                           &
     &                -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j+1,k)-wp(i+2,j+1,k)        &
     &                -2.e0*wp(i-1,j+1,k)-3.e0*wp(i,j+1,k))

                    advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                  else

                    a=dxt36e(i)*(wp(i+2,j-1,k)-3.e0*wp(i+1,j-1,k)       &
     &                +3.e0*wp(i,j-1,k)-wp(i-1,j-1,k))

                    b=dxt22(i)*(wp(i+1,j-1,k)                           &
     &                -2.e0*wp(i,j-1,k)+wp(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j-1,k)-wp(i+2,j-1,k)        &
     &                -2.e0*wp(i-1,j-1,k)-3.e0*wp(i,j-1,k))

                    advjm1=wp(i,j-1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j,k)-3.e0*wp(i+1,j,k)           &
     &                +3.e0*wp(i,j,k)-wp(i-1,j,k))

                    b=dxt22(i)*(wp(i+1,j,k)                             &
     &                -2.e0*wp(i,j,k)+wp(i-1,j,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j,k)-wp(i+2,j,k)            &
     &                -2.e0*wp(i-1,j,k)-3.e0*wp(i,j,k))

                    advij=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j+1,k)-3.e0*wp(i+1,j+1,k)       &
     &                +3.e0*wp(i,j+1,k)-wp(i-1,j+1,k))

                    b=dxt22(i)*(wp(i+1,j+1,k)                           &
     &                -2.e0*wp(i,j+1,k)+wp(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j+1,k)-wp(i+2,j+1,k)        &
     &                -2.e0*wp(i-1,j+1,k)-3.e0*wp(i,j+1,k))

                    advjp1=wp(i,j+1,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dxt36e(i)*(wp(i+2,j+2,k)-3.e0*wp(i+1,j+2,k)       &
     &                +3.e0*wp(i,j+2,k)-wp(i-1,j+2,k))

                    b=dxt22(i)*(wp(i+1,j+2,k)                           &
     &                -2.e0*wp(i,j+2,k)+wp(i-1,j+2,k))

                    c=dxt16(i)*(6.e0*wp(i+1,j+2,k)-wp(i+2,j+2,k)        &
     &                -2.e0*wp(i-1,j+2,k)-3.e0*wp(i,j+2,k))

                    advjp2=wp(i,j+2,k)+((a*u8w+b)*u8w+c)*u8w

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    wf(i,j,k)=advij+((a*v8w+b)*v8w+c)*v8w

                  end if

                end if

              end do
              end do

!$omp end do

            end do

          end if

! -----

        end if

!! -----

!!! -----

!!! Perform x-y-z seperated Cubic Lagrange scheme.

      else if(advopt.eq.5) then

! Perform calculation without the map scale factor.

        if(mfcopt.eq.0) then

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j,u8w,a,b,c)

            do j=0,nj
            do i=2,ni-2

              u8w=.25e0                                                 &
     &          *((up(i,j,k-1)+up(i,j,k))+(up(i+1,j,k-1)+up(i+1,j,k)))

              if(u8w.gt.0.e0) then

                a=dxt36w(i)*(wp(i+1,j,k)-3.e0*wp(i,j,k)                 &
     &            +3.e0*wp(i-1,j,k)-wp(i-2,j,k))

                b=dxt22(i)*(wp(i+1,j,k)                                 &
     &            -2.e0*wp(i,j,k)+wp(i-1,j,k))

                c=dxt16(i)*(3.e0*wp(i,j,k)+2.e0*wp(i+1,j,k)             &
     &            +wp(i-2,j,k)-6.e0*wp(i-1,j,k))

                advd(i,j,k)=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

              else

                a=dxt36e(i)*(wp(i+2,j,k)-3.e0*wp(i+1,j,k)               &
     &            +3.e0*wp(i,j,k)-wp(i-1,j,k))

                b=dxt22(i)*(wp(i+1,j,k)                                 &
     &            -2.e0*wp(i,j,k)+wp(i-1,j,k))

                c=dxt16(i)*(6.e0*wp(i+1,j,k)-wp(i+2,j,k)                &
     &            -2.e0*wp(i-1,j,k)-3.e0*wp(i,j,k))

                advd(i,j,k)=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

              end if

            end do
            end do

!$omp end do

          end do

          do k=2,nk-1

!$omp do schedule(runtime) private(i,j,v8w,a,b,c)

            do j=2,nj-2
            do i=2,ni-2

              v8w=.25e0                                                 &
     &          *((vp(i,j,k-1)+vp(i,j,k))+(vp(i,j+1,k-1)+vp(i,j+1,k)))

              if(v8w.gt.0.e0) then

                a=dyt36s(j)*(advd(i,j+1,k)-3.e0*advd(i,j,k)             &
     &            +3.e0*advd(i,j-1,k)-advd(i,j-2,k))

                b=dyt22(j)*(advd(i,j+1,k)                               &
     &            -2.e0*advd(i,j,k)+advd(i,j-1,k))

                c=dyt16(j)*(3.e0*advd(i,j,k)+2.e0*advd(i,j+1,k)         &
     &            +advd(i,j-2,k)-6.e0*advd(i,j-1,k))

                wf(i,j,k)=advd(i,j,k)+((a*v8w+b)*v8w+c)*v8w

              else

                a=dyt36n(j)*(advd(i,j+2,k)-3.e0*advd(i,j+1,k)           &
     &            +3.e0*advd(i,j,k)-advd(i,j-1,k))

                b=dyt22(j)*(advd(i,j+1,k)                               &
     &            -2.e0*advd(i,j,k)+advd(i,j-1,k))

                c=dyt16(j)*(6.e0*advd(i,j+1,k)-advd(i,j+2,k)            &
     &            -2.e0*advd(i,j-1,k)-3.e0*advd(i,j,k))

                wf(i,j,k)=advd(i,j,k)+((a*v8w+b)*v8w+c)*v8w

              end if

            end do
            end do

!$omp end do

          end do

! -----

!! Perform calculation with the map scale factor.

        else

! Applied the map scale factor for the x direction.

          if(mpopt.eq.0.or.mpopt.eq.10) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j,u8w,a,b,c)

              do j=0,nj
              do i=2,ni-2

                u8w=.25e0*(mf8u(i,j)*(up(i,j,k-1)+up(i,j,k))            &
     &            +mf8u(i+1,j)*(up(i+1,j,k-1)+up(i+1,j,k)))

                if(u8w.gt.0.e0) then

                  a=dxt36w(i)*(wp(i+1,j,k)-3.e0*wp(i,j,k)               &
     &              +3.e0*wp(i-1,j,k)-wp(i-2,j,k))

                  b=dxt22(i)*(wp(i+1,j,k)                               &
     &              -2.e0*wp(i,j,k)+wp(i-1,j,k))

                  c=dxt16(i)*(3.e0*wp(i,j,k)+2.e0*wp(i+1,j,k)           &
     &              +wp(i-2,j,k)-6.e0*wp(i-1,j,k))

                  advd(i,j,k)=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                else

                  a=dxt36e(i)*(wp(i+2,j,k)-3.e0*wp(i+1,j,k)             &
     &              +3.e0*wp(i,j,k)-wp(i-1,j,k))

                  b=dxt22(i)*(wp(i+1,j,k)                               &
     &              -2.e0*wp(i,j,k)+wp(i-1,j,k))

                  c=dxt16(i)*(6.e0*wp(i+1,j,k)-wp(i+2,j,k)              &
     &              -2.e0*wp(i-1,j,k)-3.e0*wp(i,j,k))

                  advd(i,j,k)=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                end if

              end do
              end do

!$omp end do

            end do

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j,v8w,a,b,c)

              do j=2,nj-2
              do i=2,ni-2

                v8w=.25e0                                               &
     &            *((vp(i,j,k-1)+vp(i,j,k))+(vp(i,j+1,k-1)+vp(i,j+1,k)))

                if(v8w.gt.0.e0) then

                  a=dyt36s(j)*(advd(i,j+1,k)-3.e0*advd(i,j,k)           &
     &              +3.e0*advd(i,j-1,k)-advd(i,j-2,k))

                  b=dyt22(j)*(advd(i,j+1,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i,j-1,k))

                  c=dyt16(j)*(3.e0*advd(i,j,k)+2.e0*advd(i,j+1,k)       &
     &              +advd(i,j-2,k)-6.e0*advd(i,j-1,k))

                  wf(i,j,k)=advd(i,j,k)+((a*v8w+b)*v8w+c)*v8w

                else

                  a=dyt36n(j)*(advd(i,j+2,k)-3.e0*advd(i,j+1,k)         &
     &              +3.e0*advd(i,j,k)-advd(i,j-1,k))

                  b=dyt22(j)*(advd(i,j+1,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i,j-1,k))

                  c=dyt16(j)*(6.e0*advd(i,j+1,k)-advd(i,j+2,k)          &
     &              -2.e0*advd(i,j-1,k)-3.e0*advd(i,j,k))

                  wf(i,j,k)=advd(i,j,k)+((a*v8w+b)*v8w+c)*v8w

                end if

              end do
              end do

!$omp end do

            end do

! -----

! Applied the map scale factor for the y direction.

          else if(mpopt.eq.5) then

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j,u8w,a,b,c)

              do j=0,nj
              do i=2,ni-2

                u8w=.25e0                                               &
     &            *((up(i,j,k-1)+up(i,j,k))+(up(i+1,j,k-1)+up(i+1,j,k)))

                if(u8w.gt.0.e0) then

                  a=dxt36w(i)*(wp(i+1,j,k)-3.e0*wp(i,j,k)               &
     &              +3.e0*wp(i-1,j,k)-wp(i-2,j,k))

                  b=dxt22(i)*(wp(i+1,j,k)                               &
     &              -2.e0*wp(i,j,k)+wp(i-1,j,k))

                  c=dxt16(i)*(3.e0*wp(i,j,k)+2.e0*wp(i+1,j,k)           &
     &              +wp(i-2,j,k)-6.e0*wp(i-1,j,k))

                  advd(i,j,k)=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                else

                  a=dxt36e(i)*(wp(i+2,j,k)-3.e0*wp(i+1,j,k)             &
     &              +3.e0*wp(i,j,k)-wp(i-1,j,k))

                  b=dxt22(i)*(wp(i+1,j,k)                               &
     &              -2.e0*wp(i,j,k)+wp(i-1,j,k))

                  c=dxt16(i)*(6.e0*wp(i+1,j,k)-wp(i+2,j,k)              &
     &              -2.e0*wp(i-1,j,k)-3.e0*wp(i,j,k))

                  advd(i,j,k)=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                end if

              end do
              end do

!$omp end do

            end do

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j,v8w,a,b,c)

              do j=2,nj-2
              do i=2,ni-2

                v8w=.25e0*(mf8v(i,j)*(vp(i,j,k-1)+vp(i,j,k))            &
     &            +mf8v(i,j+1)*(vp(i,j+1,k-1)+vp(i,j+1,k)))

                if(v8w.gt.0.e0) then

                  a=dyt36s(j)*(advd(i,j+1,k)-3.e0*advd(i,j,k)           &
     &              +3.e0*advd(i,j-1,k)-advd(i,j-2,k))

                  b=dyt22(j)*(advd(i,j+1,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i,j-1,k))

                  c=dyt16(j)*(3.e0*advd(i,j,k)+2.e0*advd(i,j+1,k)       &
     &              +advd(i,j-2,k)-6.e0*advd(i,j-1,k))

                  wf(i,j,k)=advd(i,j,k)+((a*v8w+b)*v8w+c)*v8w

                else

                  a=dyt36n(j)*(advd(i,j+2,k)-3.e0*advd(i,j+1,k)         &
     &              +3.e0*advd(i,j,k)-advd(i,j-1,k))

                  b=dyt22(j)*(advd(i,j+1,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i,j-1,k))

                  c=dyt16(j)*(6.e0*advd(i,j+1,k)-advd(i,j+2,k)          &
     &              -2.e0*advd(i,j-1,k)-3.e0*advd(i,j,k))

                  wf(i,j,k)=advd(i,j,k)+((a*v8w+b)*v8w+c)*v8w

                end if

              end do
              end do

!$omp end do

            end do

! -----

! Applied the map scale factor for the x and y direction.

          else

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j,u8w,a,b,c)

              do j=0,nj
              do i=2,ni-2

                u8w=.25e0*(mf8u(i,j)*(up(i,j,k-1)+up(i,j,k))            &
     &            +mf8u(i+1,j)*(up(i+1,j,k-1)+up(i+1,j,k)))

                if(u8w.gt.0.e0) then

                  a=dxt36w(i)*(wp(i+1,j,k)-3.e0*wp(i,j,k)               &
     &              +3.e0*wp(i-1,j,k)-wp(i-2,j,k))

                  b=dxt22(i)*(wp(i+1,j,k)                               &
     &              -2.e0*wp(i,j,k)+wp(i-1,j,k))

                  c=dxt16(i)*(3.e0*wp(i,j,k)+2.e0*wp(i+1,j,k)           &
     &              +wp(i-2,j,k)-6.e0*wp(i-1,j,k))

                  advd(i,j,k)=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                else

                  a=dxt36e(i)*(wp(i+2,j,k)-3.e0*wp(i+1,j,k)             &
     &              +3.e0*wp(i,j,k)-wp(i-1,j,k))

                  b=dxt22(i)*(wp(i+1,j,k)                               &
     &              -2.e0*wp(i,j,k)+wp(i-1,j,k))

                  c=dxt16(i)*(6.e0*wp(i+1,j,k)-wp(i+2,j,k)              &
     &              -2.e0*wp(i-1,j,k)-3.e0*wp(i,j,k))

                  advd(i,j,k)=wp(i,j,k)+((a*u8w+b)*u8w+c)*u8w

                end if

              end do
              end do

!$omp end do

            end do

            do k=2,nk-1

!$omp do schedule(runtime) private(i,j,v8w,a,b,c)

              do j=2,nj-2
              do i=2,ni-2

                v8w=.25e0*(mf8v(i,j)*(vp(i,j,k-1)+vp(i,j,k))            &
     &            +mf8v(i,j+1)*(vp(i,j+1,k-1)+vp(i,j+1,k)))

                if(v8w.gt.0.e0) then

                  a=dyt36s(j)*(advd(i,j+1,k)-3.e0*advd(i,j,k)           &
     &              +3.e0*advd(i,j-1,k)-advd(i,j-2,k))

                  b=dyt22(j)*(advd(i,j+1,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i,j-1,k))

                  c=dyt16(j)*(3.e0*advd(i,j,k)+2.e0*advd(i,j+1,k)       &
     &              +advd(i,j-2,k)-6.e0*advd(i,j-1,k))

                  wf(i,j,k)=advd(i,j,k)+((a*v8w+b)*v8w+c)*v8w

                else

                  a=dyt36n(j)*(advd(i,j+2,k)-3.e0*advd(i,j+1,k)         &
     &              +3.e0*advd(i,j,k)-advd(i,j-1,k))

                  b=dyt22(j)*(advd(i,j+1,k)                             &
     &              -2.e0*advd(i,j,k)+advd(i,j-1,k))

                  c=dyt16(j)*(6.e0*advd(i,j+1,k)-advd(i,j+2,k)          &
     &              -2.e0*advd(i,j-1,k)-3.e0*advd(i,j,k))

                  wf(i,j,k)=advd(i,j,k)+((a*v8w+b)*v8w+c)*v8w

                end if

              end do
              end do

!$omp end do

            end do

          end if

! -----

        end if

!! -----

      end if

!!! -----

!$omp end parallel

!!!! -----

!!!!! -----

      end subroutine s_hculuvw

!-----7--------------------------------------------------------------7--

      end module m_hculuvw
