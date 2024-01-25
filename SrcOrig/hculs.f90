!***********************************************************************
      module m_hculs
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/04/03
!     Modification: 2006/05/12, 2006/06/21, 2006/11/06, 2007/01/31,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2011/07/15, 2011/08/09, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the horizontal scalar advection by the Cubic Lagrange
!     scheme.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_commath
      use m_getiname
      use m_getrname
      use m_lbculs

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: hculs, s_hculs

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface hculs

        module procedure s_hculs

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
      subroutine s_hculs(fpmpopt,fpmfcopt,fpadvopt,fpdxiv,fpdyiv,dtb,   &
     &                   ni,nj,nk,mf8u,mf8v,up,vp,sp,sf,advx,           &
     &                   dxt22,dxt16,dxt36w,dxt36e,                     &
     &                   dyt22,dyt16,dyt36s,dyt36n)
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

      real, intent(in) :: up(0:ni+1,0:nj+1,1:nk)
                       ! x components of velocity at past

      real, intent(in) :: vp(0:ni+1,0:nj+1,1:nk)
                       ! y components of velocity at past

! Input and output variable

      real, intent(inout) :: sp(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at past

! Output variable

      real, intent(out) :: sf(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at future

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

      real, intent(inout) :: advx(0:ni+1,0:nj+1,1:nk)
                       ! Advected value in x direction

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

      real u8s         ! u at scalar points
      real v8s         ! v at scalar points

      real advij       ! Temporary variable

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

!!!!! Calculate the scalar advection horizontally.

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

      call lbculs(idwbc,idebc,idsbc,idnbc,ni,nj,nk,sp,                  &
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
!$omp&   private(i,j,u8s,v8s,advij,advjm2,advjm1,advjp1,advjp2,a,b,c)

            do j=2,nj-2
            do i=2,ni-2

              u8s=.5e0*(up(i,j,k)+up(i+1,j,k))
              v8s=.5e0*(vp(i,j,k)+vp(i,j+1,k))

              if(u8s.gt.0.e0) then

                if(v8s.gt.0.e0) then

                  a=dxt36w(i)*(sp(i+1,j-2,k)-3.e0*sp(i,j-2,k)           &
     &              +3.e0*sp(i-1,j-2,k)-sp(i-2,j-2,k))

                  b=dxt22(i)*(sp(i+1,j-2,k)                             &
     &              -2.e0*sp(i,j-2,k)+sp(i-1,j-2,k))

                  c=dxt16(i)*(3.e0*sp(i,j-2,k)+2.e0*sp(i+1,j-2,k)       &
     &              +sp(i-2,j-2,k)-6.e0*sp(i-1,j-2,k))

                  advjm2=sp(i,j-2,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dxt36w(i)*(sp(i+1,j-1,k)-3.e0*sp(i,j-1,k)           &
     &              +3.e0*sp(i-1,j-1,k)-sp(i-2,j-1,k))

                  b=dxt22(i)*(sp(i+1,j-1,k)                             &
     &              -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                  c=dxt16(i)*(3.e0*sp(i,j-1,k)+2.e0*sp(i+1,j-1,k)       &
     &              +sp(i-2,j-1,k)-6.e0*sp(i-1,j-1,k))

                  advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dxt36w(i)*(sp(i+1,j,k)-3.e0*sp(i,j,k)               &
     &              +3.e0*sp(i-1,j,k)-sp(i-2,j,k))

                  b=dxt22(i)*(sp(i+1,j,k)                               &
     &              -2.e0*sp(i,j,k)+sp(i-1,j,k))

                  c=dxt16(i)*(3.e0*sp(i,j,k)+2.e0*sp(i+1,j,k)           &
     &              +sp(i-2,j,k)-6.e0*sp(i-1,j,k))

                  advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dxt36w(i)*(sp(i+1,j+1,k)-3.e0*sp(i,j+1,k)           &
     &              +3.e0*sp(i-1,j+1,k)-sp(i-2,j+1,k))

                  b=dxt22(i)*(sp(i+1,j+1,k)                             &
     &              -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                  c=dxt16(i)*(3.e0*sp(i,j+1,k)+2.e0*sp(i+1,j+1,k)       &
     &              +sp(i-2,j+1,k)-6.e0*sp(i-1,j+1,k))

                  advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                  b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                  c=dyt16(j)                                            &
     &              *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                  sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

                else

                  a=dxt36w(i)*(sp(i+1,j-1,k)-3.e0*sp(i,j-1,k)           &
     &              +3.e0*sp(i-1,j-1,k)-sp(i-2,j-1,k))

                  b=dxt22(i)*(sp(i+1,j-1,k)                             &
     &              -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                  c=dxt16(i)*(3.e0*sp(i,j-1,k)+2.e0*sp(i+1,j-1,k)       &
     &              +sp(i-2,j-1,k)-6.e0*sp(i-1,j-1,k))

                  advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dxt36w(i)*(sp(i+1,j,k)-3.e0*sp(i,j,k)               &
     &              +3.e0*sp(i-1,j,k)-sp(i-2,j,k))

                  b=dxt22(i)*(sp(i+1,j,k)                               &
     &              -2.e0*sp(i,j,k)+sp(i-1,j,k))

                  c=dxt16(i)*(3.e0*sp(i,j,k)+2.e0*sp(i+1,j,k)           &
     &              +sp(i-2,j,k)-6.e0*sp(i-1,j,k))

                  advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dxt36w(i)*(sp(i+1,j+1,k)-3.e0*sp(i,j+1,k)           &
     &              +3.e0*sp(i-1,j+1,k)-sp(i-2,j+1,k))

                  b=dxt22(i)*(sp(i+1,j+1,k)                             &
     &              -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                  c=dxt16(i)*(3.e0*sp(i,j+1,k)+2.e0*sp(i+1,j+1,k)       &
     &              +sp(i-2,j+1,k)-6.e0*sp(i-1,j+1,k))

                  advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dxt36w(i)*(sp(i+1,j+2,k)-3.e0*sp(i,j+2,k)           &
     &              +3.e0*sp(i-1,j+2,k)-sp(i-2,j+2,k))

                  b=dxt22(i)*(sp(i+1,j+2,k)                             &
     &              -2.e0*sp(i,j+2,k)+sp(i-1,j+2,k))

                  c=dxt16(i)*(3.e0*sp(i,j+2,k)+2.e0*sp(i+1,j+2,k)       &
     &              +sp(i-2,j+2,k)-6.e0*sp(i-1,j+2,k))

                  advjp2=sp(i,j+2,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                  b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                  c=dyt16(j)                                            &
     &              *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                  sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

                end if

              else

                if(v8s.gt.0.e0) then

                  a=dxt36e(i)*(sp(i+2,j-2,k)-3.e0*sp(i+1,j-2,k)         &
     &              +3.e0*sp(i,j-2,k)-sp(i-1,j-2,k))

                  b=dxt22(i)*(sp(i+1,j-2,k)                             &
     &              -2.e0*sp(i,j-2,k)+sp(i-1,j-2,k))

                  c=dxt16(i)*(6.e0*sp(i+1,j-2,k)-sp(i+2,j-2,k)          &
     &              -2.e0*sp(i-1,j-2,k)-3.e0*sp(i,j-2,k))

                  advjm2=sp(i,j-2,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dxt36e(i)*(sp(i+2,j-1,k)-3.e0*sp(i+1,j-1,k)         &
     &              +3.e0*sp(i,j-1,k)-sp(i-1,j-1,k))

                  b=dxt22(i)*(sp(i+1,j-1,k)                             &
     &              -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                  c=dxt16(i)*(6.e0*sp(i+1,j-1,k)-sp(i+2,j-1,k)          &
     &              -2.e0*sp(i-1,j-1,k)-3.e0*sp(i,j-1,k))

                  advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dxt36e(i)*(sp(i+2,j,k)-3.e0*sp(i+1,j,k)             &
     &              +3.e0*sp(i,j,k)-sp(i-1,j,k))

                  b=dxt22(i)*(sp(i+1,j,k)                               &
     &              -2.e0*sp(i,j,k)+sp(i-1,j,k))

                  c=dxt16(i)*(6.e0*sp(i+1,j,k)-sp(i+2,j,k)              &
     &              -2.e0*sp(i-1,j,k)-3.e0*sp(i,j,k))

                  advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dxt36e(i)*(sp(i+2,j+1,k)-3.e0*sp(i+1,j+1,k)         &
     &              +3.e0*sp(i,j+1,k)-sp(i-1,j+1,k))

                  b=dxt22(i)*(sp(i+1,j+1,k)                             &
     &              -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                  c=dxt16(i)*(6.e0*sp(i+1,j+1,k)-sp(i+2,j+1,k)          &
     &              -2.e0*sp(i-1,j+1,k)-3.e0*sp(i,j+1,k))

                  advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                  b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                  c=dyt16(j)                                            &
     &              *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                  sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

                else

                  a=dxt36e(i)*(sp(i+2,j-1,k)-3.e0*sp(i+1,j-1,k)         &
     &              +3.e0*sp(i,j-1,k)-sp(i-1,j-1,k))

                  b=dxt22(i)*(sp(i+1,j-1,k)                             &
     &              -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                  c=dxt16(i)*(6.e0*sp(i+1,j-1,k)-sp(i+2,j-1,k)          &
     &              -2.e0*sp(i-1,j-1,k)-3.e0*sp(i,j-1,k))

                  advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dxt36e(i)*(sp(i+2,j,k)-3.e0*sp(i+1,j,k)             &
     &              +3.e0*sp(i,j,k)-sp(i-1,j,k))

                  b=dxt22(i)*(sp(i+1,j,k)                               &
     &              -2.e0*sp(i,j,k)+sp(i-1,j,k))

                  c=dxt16(i)*(6.e0*sp(i+1,j,k)-sp(i+2,j,k)              &
     &              -2.e0*sp(i-1,j,k)-3.e0*sp(i,j,k))

                  advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dxt36e(i)*(sp(i+2,j+1,k)-3.e0*sp(i+1,j+1,k)         &
     &              +3.e0*sp(i,j+1,k)-sp(i-1,j+1,k))

                  b=dxt22(i)*(sp(i+1,j+1,k)                             &
     &              -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                  c=dxt16(i)*(6.e0*sp(i+1,j+1,k)-sp(i+2,j+1,k)          &
     &              -2.e0*sp(i-1,j+1,k)-3.e0*sp(i,j+1,k))

                  advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dxt36e(i)*(sp(i+2,j+2,k)-3.e0*sp(i+1,j+2,k)         &
     &              +3.e0*sp(i,j+2,k)-sp(i-1,j+2,k))

                  b=dxt22(i)*(sp(i+1,j+2,k)                             &
     &              -2.e0*sp(i,j+2,k)+sp(i-1,j+2,k))

                  c=dxt16(i)*(6.e0*sp(i+1,j+2,k)-sp(i+2,j+2,k)          &
     &              -2.e0*sp(i-1,j+2,k)-3.e0*sp(i,j+2,k))

                  advjp2=sp(i,j+2,k)+((a*u8s+b)*u8s+c)*u8s

                  a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                  b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                  c=dyt16(j)                                            &
     &              *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                  sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

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
!$omp&   private(i,j,u8s,v8s,advij,advjm2,advjm1,advjp1,advjp2,a,b,c)

              do j=2,nj-2
              do i=2,ni-2

                u8s=.5e0*(mf8u(i,j)*up(i,j,k)+mf8u(i+1,j)*up(i+1,j,k))
                v8s=.5e0*(vp(i,j,k)+vp(i,j+1,k))

                if(u8s.gt.0.e0) then

                  if(v8s.gt.0.e0) then

                    a=dxt36w(i)*(sp(i+1,j-2,k)-3.e0*sp(i,j-2,k)         &
     &                +3.e0*sp(i-1,j-2,k)-sp(i-2,j-2,k))

                    b=dxt22(i)*(sp(i+1,j-2,k)                           &
     &                -2.e0*sp(i,j-2,k)+sp(i-1,j-2,k))

                    c=dxt16(i)*(3.e0*sp(i,j-2,k)+2.e0*sp(i+1,j-2,k)     &
     &                +sp(i-2,j-2,k)-6.e0*sp(i-1,j-2,k))

                    advjm2=sp(i,j-2,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j-1,k)-3.e0*sp(i,j-1,k)         &
     &                +3.e0*sp(i-1,j-1,k)-sp(i-2,j-1,k))

                    b=dxt22(i)*(sp(i+1,j-1,k)                           &
     &                -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*sp(i,j-1,k)+2.e0*sp(i+1,j-1,k)     &
     &                +sp(i-2,j-1,k)-6.e0*sp(i-1,j-1,k))

                    advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j,k)-3.e0*sp(i,j,k)             &
     &                +3.e0*sp(i-1,j,k)-sp(i-2,j,k))

                    b=dxt22(i)*(sp(i+1,j,k)                             &
     &                -2.e0*sp(i,j,k)+sp(i-1,j,k))

                    c=dxt16(i)*(3.e0*sp(i,j,k)+2.e0*sp(i+1,j,k)         &
     &                +sp(i-2,j,k)-6.e0*sp(i-1,j,k))

                    advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j+1,k)-3.e0*sp(i,j+1,k)         &
     &                +3.e0*sp(i-1,j+1,k)-sp(i-2,j+1,k))

                    b=dxt22(i)*(sp(i+1,j+1,k)                           &
     &                -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*sp(i,j+1,k)+2.e0*sp(i+1,j+1,k)     &
     &                +sp(i-2,j+1,k)-6.e0*sp(i-1,j+1,k))

                    advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

                  else

                    a=dxt36w(i)*(sp(i+1,j-1,k)-3.e0*sp(i,j-1,k)         &
     &                +3.e0*sp(i-1,j-1,k)-sp(i-2,j-1,k))

                    b=dxt22(i)*(sp(i+1,j-1,k)                           &
     &                -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*sp(i,j-1,k)+2.e0*sp(i+1,j-1,k)     &
     &                +sp(i-2,j-1,k)-6.e0*sp(i-1,j-1,k))

                    advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j,k)-3.e0*sp(i,j,k)             &
     &                +3.e0*sp(i-1,j,k)-sp(i-2,j,k))

                    b=dxt22(i)*(sp(i+1,j,k)                             &
     &                -2.e0*sp(i,j,k)+sp(i-1,j,k))

                    c=dxt16(i)*(3.e0*sp(i,j,k)+2.e0*sp(i+1,j,k)         &
     &                +sp(i-2,j,k)-6.e0*sp(i-1,j,k))

                    advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j+1,k)-3.e0*sp(i,j+1,k)         &
     &                +3.e0*sp(i-1,j+1,k)-sp(i-2,j+1,k))

                    b=dxt22(i)*(sp(i+1,j+1,k)                           &
     &                -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*sp(i,j+1,k)+2.e0*sp(i+1,j+1,k)     &
     &                +sp(i-2,j+1,k)-6.e0*sp(i-1,j+1,k))

                    advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j+2,k)-3.e0*sp(i,j+2,k)         &
     &                +3.e0*sp(i-1,j+2,k)-sp(i-2,j+2,k))

                    b=dxt22(i)*(sp(i+1,j+2,k)                           &
     &                -2.e0*sp(i,j+2,k)+sp(i-1,j+2,k))

                    c=dxt16(i)*(3.e0*sp(i,j+2,k)+2.e0*sp(i+1,j+2,k)     &
     &                +sp(i-2,j+2,k)-6.e0*sp(i-1,j+2,k))

                    advjp2=sp(i,j+2,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

                  end if

                else

                  if(v8s.gt.0.e0) then

                    a=dxt36e(i)*(sp(i+2,j-2,k)-3.e0*sp(i+1,j-2,k)       &
     &                +3.e0*sp(i,j-2,k)-sp(i-1,j-2,k))

                    b=dxt22(i)*(sp(i+1,j-2,k)                           &
     &                -2.e0*sp(i,j-2,k)+sp(i-1,j-2,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j-2,k)-sp(i+2,j-2,k)        &
     &                -2.e0*sp(i-1,j-2,k)-3.e0*sp(i,j-2,k))

                    advjm2=sp(i,j-2,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j-1,k)-3.e0*sp(i+1,j-1,k)       &
     &                +3.e0*sp(i,j-1,k)-sp(i-1,j-1,k))

                    b=dxt22(i)*(sp(i+1,j-1,k)                           &
     &                -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j-1,k)-sp(i+2,j-1,k)        &
     &                -2.e0*sp(i-1,j-1,k)-3.e0*sp(i,j-1,k))

                    advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j,k)-3.e0*sp(i+1,j,k)           &
     &                +3.e0*sp(i,j,k)-sp(i-1,j,k))

                    b=dxt22(i)*(sp(i+1,j,k)                             &
     &                -2.e0*sp(i,j,k)+sp(i-1,j,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j,k)-sp(i+2,j,k)            &
     &                -2.e0*sp(i-1,j,k)-3.e0*sp(i,j,k))

                    advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j+1,k)-3.e0*sp(i+1,j+1,k)       &
     &                +3.e0*sp(i,j+1,k)-sp(i-1,j+1,k))

                    b=dxt22(i)*(sp(i+1,j+1,k)                           &
     &                -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j+1,k)-sp(i+2,j+1,k)        &
     &                -2.e0*sp(i-1,j+1,k)-3.e0*sp(i,j+1,k))

                    advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

                  else

                    a=dxt36e(i)*(sp(i+2,j-1,k)-3.e0*sp(i+1,j-1,k)       &
     &                +3.e0*sp(i,j-1,k)-sp(i-1,j-1,k))

                    b=dxt22(i)*(sp(i+1,j-1,k)                           &
     &                -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j-1,k)-sp(i+2,j-1,k)        &
     &                -2.e0*sp(i-1,j-1,k)-3.e0*sp(i,j-1,k))

                    advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j,k)-3.e0*sp(i+1,j,k)           &
     &                +3.e0*sp(i,j,k)-sp(i-1,j,k))

                    b=dxt22(i)*(sp(i+1,j,k)                             &
     &                -2.e0*sp(i,j,k)+sp(i-1,j,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j,k)-sp(i+2,j,k)            &
     &                -2.e0*sp(i-1,j,k)-3.e0*sp(i,j,k))

                    advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j+1,k)-3.e0*sp(i+1,j+1,k)       &
     &                +3.e0*sp(i,j+1,k)-sp(i-1,j+1,k))

                    b=dxt22(i)*(sp(i+1,j+1,k)                           &
     &                -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j+1,k)-sp(i+2,j+1,k)        &
     &                -2.e0*sp(i-1,j+1,k)-3.e0*sp(i,j+1,k))

                    advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j+2,k)-3.e0*sp(i+1,j+2,k)       &
     &                +3.e0*sp(i,j+2,k)-sp(i-1,j+2,k))

                    b=dxt22(i)*(sp(i+1,j+2,k)                           &
     &                -2.e0*sp(i,j+2,k)+sp(i-1,j+2,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j+2,k)-sp(i+2,j+2,k)        &
     &                -2.e0*sp(i-1,j+2,k)-3.e0*sp(i,j+2,k))

                    advjp2=sp(i,j+2,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

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
!$omp&   private(i,j,u8s,v8s,advij,advjm2,advjm1,advjp1,advjp2,a,b,c)

              do j=2,nj-2
              do i=2,ni-2

                u8s=.5e0*(up(i,j,k)+up(i+1,j,k))
                v8s=.5e0*(mf8v(i,j)*vp(i,j,k)+mf8v(i,j+1)*vp(i,j+1,k))

                if(u8s.gt.0.e0) then

                  if(v8s.gt.0.e0) then

                    a=dxt36w(i)*(sp(i+1,j-2,k)-3.e0*sp(i,j-2,k)         &
     &                +3.e0*sp(i-1,j-2,k)-sp(i-2,j-2,k))

                    b=dxt22(i)*(sp(i+1,j-2,k)                           &
     &                -2.e0*sp(i,j-2,k)+sp(i-1,j-2,k))

                    c=dxt16(i)*(3.e0*sp(i,j-2,k)+2.e0*sp(i+1,j-2,k)     &
     &                +sp(i-2,j-2,k)-6.e0*sp(i-1,j-2,k))

                    advjm2=sp(i,j-2,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j-1,k)-3.e0*sp(i,j-1,k)         &
     &                +3.e0*sp(i-1,j-1,k)-sp(i-2,j-1,k))

                    b=dxt22(i)*(sp(i+1,j-1,k)                           &
     &                -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*sp(i,j-1,k)+2.e0*sp(i+1,j-1,k)     &
     &                +sp(i-2,j-1,k)-6.e0*sp(i-1,j-1,k))

                    advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j,k)-3.e0*sp(i,j,k)             &
     &                +3.e0*sp(i-1,j,k)-sp(i-2,j,k))

                    b=dxt22(i)*(sp(i+1,j,k)                             &
     &                -2.e0*sp(i,j,k)+sp(i-1,j,k))

                    c=dxt16(i)*(3.e0*sp(i,j,k)+2.e0*sp(i+1,j,k)         &
     &                +sp(i-2,j,k)-6.e0*sp(i-1,j,k))

                    advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j+1,k)-3.e0*sp(i,j+1,k)         &
     &                +3.e0*sp(i-1,j+1,k)-sp(i-2,j+1,k))

                    b=dxt22(i)*(sp(i+1,j+1,k)                           &
     &                -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*sp(i,j+1,k)+2.e0*sp(i+1,j+1,k)     &
     &                +sp(i-2,j+1,k)-6.e0*sp(i-1,j+1,k))

                    advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

                  else

                    a=dxt36w(i)*(sp(i+1,j-1,k)-3.e0*sp(i,j-1,k)         &
     &                +3.e0*sp(i-1,j-1,k)-sp(i-2,j-1,k))

                    b=dxt22(i)*(sp(i+1,j-1,k)                           &
     &                -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*sp(i,j-1,k)+2.e0*sp(i+1,j-1,k)     &
     &                +sp(i-2,j-1,k)-6.e0*sp(i-1,j-1,k))

                    advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j,k)-3.e0*sp(i,j,k)             &
     &                +3.e0*sp(i-1,j,k)-sp(i-2,j,k))

                    b=dxt22(i)*(sp(i+1,j,k)                             &
     &                -2.e0*sp(i,j,k)+sp(i-1,j,k))

                    c=dxt16(i)*(3.e0*sp(i,j,k)+2.e0*sp(i+1,j,k)         &
     &                +sp(i-2,j,k)-6.e0*sp(i-1,j,k))

                    advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j+1,k)-3.e0*sp(i,j+1,k)         &
     &                +3.e0*sp(i-1,j+1,k)-sp(i-2,j+1,k))

                    b=dxt22(i)*(sp(i+1,j+1,k)                           &
     &                -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*sp(i,j+1,k)+2.e0*sp(i+1,j+1,k)     &
     &                +sp(i-2,j+1,k)-6.e0*sp(i-1,j+1,k))

                    advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j+2,k)-3.e0*sp(i,j+2,k)         &
     &                +3.e0*sp(i-1,j+2,k)-sp(i-2,j+2,k))

                    b=dxt22(i)*(sp(i+1,j+2,k)                           &
     &                -2.e0*sp(i,j+2,k)+sp(i-1,j+2,k))

                    c=dxt16(i)*(3.e0*sp(i,j+2,k)+2.e0*sp(i+1,j+2,k)     &
     &                +sp(i-2,j+2,k)-6.e0*sp(i-1,j+2,k))

                    advjp2=sp(i,j+2,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

                  end if

                else

                  if(v8s.gt.0.e0) then

                    a=dxt36e(i)*(sp(i+2,j-2,k)-3.e0*sp(i+1,j-2,k)       &
     &                +3.e0*sp(i,j-2,k)-sp(i-1,j-2,k))

                    b=dxt22(i)*(sp(i+1,j-2,k)                           &
     &                -2.e0*sp(i,j-2,k)+sp(i-1,j-2,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j-2,k)-sp(i+2,j-2,k)        &
     &                -2.e0*sp(i-1,j-2,k)-3.e0*sp(i,j-2,k))

                    advjm2=sp(i,j-2,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j-1,k)-3.e0*sp(i+1,j-1,k)       &
     &                +3.e0*sp(i,j-1,k)-sp(i-1,j-1,k))

                    b=dxt22(i)*(sp(i+1,j-1,k)                           &
     &                -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j-1,k)-sp(i+2,j-1,k)        &
     &                -2.e0*sp(i-1,j-1,k)-3.e0*sp(i,j-1,k))

                    advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j,k)-3.e0*sp(i+1,j,k)           &
     &                +3.e0*sp(i,j,k)-sp(i-1,j,k))

                    b=dxt22(i)*(sp(i+1,j,k)                             &
     &                -2.e0*sp(i,j,k)+sp(i-1,j,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j,k)-sp(i+2,j,k)            &
     &                -2.e0*sp(i-1,j,k)-3.e0*sp(i,j,k))

                    advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j+1,k)-3.e0*sp(i+1,j+1,k)       &
     &                +3.e0*sp(i,j+1,k)-sp(i-1,j+1,k))

                    b=dxt22(i)*(sp(i+1,j+1,k)                           &
     &                -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j+1,k)-sp(i+2,j+1,k)        &
     &                -2.e0*sp(i-1,j+1,k)-3.e0*sp(i,j+1,k))

                    advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

                  else

                    a=dxt36e(i)*(sp(i+2,j-1,k)-3.e0*sp(i+1,j-1,k)       &
     &                +3.e0*sp(i,j-1,k)-sp(i-1,j-1,k))

                    b=dxt22(i)*(sp(i+1,j-1,k)                           &
     &                -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j-1,k)-sp(i+2,j-1,k)        &
     &                -2.e0*sp(i-1,j-1,k)-3.e0*sp(i,j-1,k))

                    advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j,k)-3.e0*sp(i+1,j,k)           &
     &                +3.e0*sp(i,j,k)-sp(i-1,j,k))

                    b=dxt22(i)*(sp(i+1,j,k)                             &
     &                -2.e0*sp(i,j,k)+sp(i-1,j,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j,k)-sp(i+2,j,k)            &
     &                -2.e0*sp(i-1,j,k)-3.e0*sp(i,j,k))

                    advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j+1,k)-3.e0*sp(i+1,j+1,k)       &
     &                +3.e0*sp(i,j+1,k)-sp(i-1,j+1,k))

                    b=dxt22(i)*(sp(i+1,j+1,k)                           &
     &                -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j+1,k)-sp(i+2,j+1,k)        &
     &                -2.e0*sp(i-1,j+1,k)-3.e0*sp(i,j+1,k))

                    advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j+2,k)-3.e0*sp(i+1,j+2,k)       &
     &                +3.e0*sp(i,j+2,k)-sp(i-1,j+2,k))

                    b=dxt22(i)*(sp(i+1,j+2,k)                           &
     &                -2.e0*sp(i,j+2,k)+sp(i-1,j+2,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j+2,k)-sp(i+2,j+2,k)        &
     &                -2.e0*sp(i-1,j+2,k)-3.e0*sp(i,j+2,k))

                    advjp2=sp(i,j+2,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

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
!$omp&   private(i,j,u8s,v8s,advij,advjm2,advjm1,advjp1,advjp2,a,b,c)

              do j=2,nj-2
              do i=2,ni-2

                u8s=.5e0*(mf8u(i,j)*up(i,j,k)+mf8u(i+1,j)*up(i+1,j,k))
                v8s=.5e0*(mf8v(i,j)*vp(i,j,k)+mf8v(i,j+1)*vp(i,j+1,k))

                if(u8s.gt.0.e0) then

                  if(v8s.gt.0.e0) then

                    a=dxt36w(i)*(sp(i+1,j-2,k)-3.e0*sp(i,j-2,k)         &
     &                +3.e0*sp(i-1,j-2,k)-sp(i-2,j-2,k))

                    b=dxt22(i)*(sp(i+1,j-2,k)                           &
     &                -2.e0*sp(i,j-2,k)+sp(i-1,j-2,k))

                    c=dxt16(i)*(3.e0*sp(i,j-2,k)+2.e0*sp(i+1,j-2,k)     &
     &                +sp(i-2,j-2,k)-6.e0*sp(i-1,j-2,k))

                    advjm2=sp(i,j-2,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j-1,k)-3.e0*sp(i,j-1,k)         &
     &                +3.e0*sp(i-1,j-1,k)-sp(i-2,j-1,k))

                    b=dxt22(i)*(sp(i+1,j-1,k)                           &
     &                -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*sp(i,j-1,k)+2.e0*sp(i+1,j-1,k)     &
     &                +sp(i-2,j-1,k)-6.e0*sp(i-1,j-1,k))

                    advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j,k)-3.e0*sp(i,j,k)             &
     &                +3.e0*sp(i-1,j,k)-sp(i-2,j,k))

                    b=dxt22(i)*(sp(i+1,j,k)                             &
     &                -2.e0*sp(i,j,k)+sp(i-1,j,k))

                    c=dxt16(i)*(3.e0*sp(i,j,k)+2.e0*sp(i+1,j,k)         &
     &                +sp(i-2,j,k)-6.e0*sp(i-1,j,k))

                    advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j+1,k)-3.e0*sp(i,j+1,k)         &
     &                +3.e0*sp(i-1,j+1,k)-sp(i-2,j+1,k))

                    b=dxt22(i)*(sp(i+1,j+1,k)                           &
     &                -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*sp(i,j+1,k)+2.e0*sp(i+1,j+1,k)     &
     &                +sp(i-2,j+1,k)-6.e0*sp(i-1,j+1,k))

                    advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

                  else

                    a=dxt36w(i)*(sp(i+1,j-1,k)-3.e0*sp(i,j-1,k)         &
     &                +3.e0*sp(i-1,j-1,k)-sp(i-2,j-1,k))

                    b=dxt22(i)*(sp(i+1,j-1,k)                           &
     &                -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                    c=dxt16(i)*(3.e0*sp(i,j-1,k)+2.e0*sp(i+1,j-1,k)     &
     &                +sp(i-2,j-1,k)-6.e0*sp(i-1,j-1,k))

                    advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j,k)-3.e0*sp(i,j,k)             &
     &                +3.e0*sp(i-1,j,k)-sp(i-2,j,k))

                    b=dxt22(i)*(sp(i+1,j,k)                             &
     &                -2.e0*sp(i,j,k)+sp(i-1,j,k))

                    c=dxt16(i)*(3.e0*sp(i,j,k)+2.e0*sp(i+1,j,k)         &
     &                +sp(i-2,j,k)-6.e0*sp(i-1,j,k))

                    advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j+1,k)-3.e0*sp(i,j+1,k)         &
     &                +3.e0*sp(i-1,j+1,k)-sp(i-2,j+1,k))

                    b=dxt22(i)*(sp(i+1,j+1,k)                           &
     &                -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                    c=dxt16(i)*(3.e0*sp(i,j+1,k)+2.e0*sp(i+1,j+1,k)     &
     &                +sp(i-2,j+1,k)-6.e0*sp(i-1,j+1,k))

                    advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36w(i)*(sp(i+1,j+2,k)-3.e0*sp(i,j+2,k)         &
     &                +3.e0*sp(i-1,j+2,k)-sp(i-2,j+2,k))

                    b=dxt22(i)*(sp(i+1,j+2,k)                           &
     &                -2.e0*sp(i,j+2,k)+sp(i-1,j+2,k))

                    c=dxt16(i)*(3.e0*sp(i,j+2,k)+2.e0*sp(i+1,j+2,k)     &
     &                +sp(i-2,j+2,k)-6.e0*sp(i-1,j+2,k))

                    advjp2=sp(i,j+2,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

                  end if

                else

                  if(v8s.gt.0.e0) then

                    a=dxt36e(i)*(sp(i+2,j-2,k)-3.e0*sp(i+1,j-2,k)       &
     &                +3.e0*sp(i,j-2,k)-sp(i-1,j-2,k))

                    b=dxt22(i)*(sp(i+1,j-2,k)                           &
     &                -2.e0*sp(i,j-2,k)+sp(i-1,j-2,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j-2,k)-sp(i+2,j-2,k)        &
     &                -2.e0*sp(i-1,j-2,k)-3.e0*sp(i,j-2,k))

                    advjm2=sp(i,j-2,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j-1,k)-3.e0*sp(i+1,j-1,k)       &
     &                +3.e0*sp(i,j-1,k)-sp(i-1,j-1,k))

                    b=dxt22(i)*(sp(i+1,j-1,k)                           &
     &                -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j-1,k)-sp(i+2,j-1,k)        &
     &                -2.e0*sp(i-1,j-1,k)-3.e0*sp(i,j-1,k))

                    advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j,k)-3.e0*sp(i+1,j,k)           &
     &                +3.e0*sp(i,j,k)-sp(i-1,j,k))

                    b=dxt22(i)*(sp(i+1,j,k)                             &
     &                -2.e0*sp(i,j,k)+sp(i-1,j,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j,k)-sp(i+2,j,k)            &
     &                -2.e0*sp(i-1,j,k)-3.e0*sp(i,j,k))

                    advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j+1,k)-3.e0*sp(i+1,j+1,k)       &
     &                +3.e0*sp(i,j+1,k)-sp(i-1,j+1,k))

                    b=dxt22(i)*(sp(i+1,j+1,k)                           &
     &                -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j+1,k)-sp(i+2,j+1,k)        &
     &                -2.e0*sp(i-1,j+1,k)-3.e0*sp(i,j+1,k))

                    advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dyt36s(j)*(advjp1-3.e0*advij+3.e0*advjm1-advjm2)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(3.e0*advij+2.e0*advjp1+advjm2-6.e0*advjm1)

                    sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

                  else

                    a=dxt36e(i)*(sp(i+2,j-1,k)-3.e0*sp(i+1,j-1,k)       &
     &                +3.e0*sp(i,j-1,k)-sp(i-1,j-1,k))

                    b=dxt22(i)*(sp(i+1,j-1,k)                           &
     &                -2.e0*sp(i,j-1,k)+sp(i-1,j-1,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j-1,k)-sp(i+2,j-1,k)        &
     &                -2.e0*sp(i-1,j-1,k)-3.e0*sp(i,j-1,k))

                    advjm1=sp(i,j-1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j,k)-3.e0*sp(i+1,j,k)           &
     &                +3.e0*sp(i,j,k)-sp(i-1,j,k))

                    b=dxt22(i)*(sp(i+1,j,k)                             &
     &                -2.e0*sp(i,j,k)+sp(i-1,j,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j,k)-sp(i+2,j,k)            &
     &                -2.e0*sp(i-1,j,k)-3.e0*sp(i,j,k))

                    advij=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j+1,k)-3.e0*sp(i+1,j+1,k)       &
     &                +3.e0*sp(i,j+1,k)-sp(i-1,j+1,k))

                    b=dxt22(i)*(sp(i+1,j+1,k)                           &
     &                -2.e0*sp(i,j+1,k)+sp(i-1,j+1,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j+1,k)-sp(i+2,j+1,k)        &
     &                -2.e0*sp(i-1,j+1,k)-3.e0*sp(i,j+1,k))

                    advjp1=sp(i,j+1,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dxt36e(i)*(sp(i+2,j+2,k)-3.e0*sp(i+1,j+2,k)       &
     &                +3.e0*sp(i,j+2,k)-sp(i-1,j+2,k))

                    b=dxt22(i)*(sp(i+1,j+2,k)                           &
     &                -2.e0*sp(i,j+2,k)+sp(i-1,j+2,k))

                    c=dxt16(i)*(6.e0*sp(i+1,j+2,k)-sp(i+2,j+2,k)        &
     &                -2.e0*sp(i-1,j+2,k)-3.e0*sp(i,j+2,k))

                    advjp2=sp(i,j+2,k)+((a*u8s+b)*u8s+c)*u8s

                    a=dyt36n(j)*(advjp2-3.e0*advjp1+3.e0*advij-advjm1)

                    b=dyt22(j)*(advjp1-2.e0*advij+advjm1)

                    c=dyt16(j)                                          &
     &                *(6.e0*advjp1-advjp2-2.e0*advjm1-3.e0*advij)

                    sf(i,j,k)=advij+((a*v8s+b)*v8s+c)*v8s

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

!$omp do schedule(runtime) private(i,j,u8s,a,b,c)

            do j=0,nj
            do i=2,ni-2

              u8s=.5e0*(up(i,j,k)+up(i+1,j,k))

              if(u8s.gt.0.e0) then

                a=dxt36w(i)*(sp(i+1,j,k)-3.e0*sp(i,j,k)                 &
     &            +3.e0*sp(i-1,j,k)-sp(i-2,j,k))

                b=dxt22(i)*(sp(i+1,j,k)                                 &
     &            -2.e0*sp(i,j,k)+sp(i-1,j,k))

                c=dxt16(i)*(3.e0*sp(i,j,k)+2.e0*sp(i+1,j,k)             &
     &            +sp(i-2,j,k)-6.e0*sp(i-1,j,k))

                advx(i,j,k)=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

              else

                a=dxt36e(i)*(sp(i+2,j,k)-3.e0*sp(i+1,j,k)               &
     &            +3.e0*sp(i,j,k)-sp(i-1,j,k))

                b=dxt22(i)*(sp(i+1,j,k)                                 &
     &            -2.e0*sp(i,j,k)+sp(i-1,j,k))

                c=dxt16(i)*(6.e0*sp(i+1,j,k)-sp(i+2,j,k)                &
     &            -2.e0*sp(i-1,j,k)-3.e0*sp(i,j,k))

                advx(i,j,k)=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

              end if

            end do
            end do

!$omp end do

          end do

          do k=1,nk-1

!$omp do schedule(runtime) private(i,j,v8s,a,b,c)

            do j=2,nj-2
            do i=2,ni-2

              v8s=.5e0*(vp(i,j,k)+vp(i,j+1,k))

              if(v8s.gt.0.e0) then

                a=dyt36s(j)*(advx(i,j+1,k)-3.e0*advx(i,j,k)             &
     &            +3.e0*advx(i,j-1,k)-advx(i,j-2,k))

                b=dyt22(j)*(advx(i,j+1,k)                               &
     &            -2.e0*advx(i,j,k)+advx(i,j-1,k))

                c=dyt16(j)*(3.e0*advx(i,j,k)+2.e0*advx(i,j+1,k)         &
     &            +advx(i,j-2,k)-6.e0*advx(i,j-1,k))

                sf(i,j,k)=advx(i,j,k)+((a*v8s+b)*v8s+c)*v8s

              else

                a=dyt36n(j)*(advx(i,j+2,k)-3.e0*advx(i,j+1,k)           &
     &            +3.e0*advx(i,j,k)-advx(i,j-1,k))

                b=dyt22(j)*(advx(i,j+1,k)                               &
     &            -2.e0*advx(i,j,k)+advx(i,j-1,k))

                c=dyt16(j)*(6.e0*advx(i,j+1,k)-advx(i,j+2,k)            &
     &            -2.e0*advx(i,j-1,k)-3.e0*advx(i,j,k))

                sf(i,j,k)=advx(i,j,k)+((a*v8s+b)*v8s+c)*v8s

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

!$omp do schedule(runtime) private(i,j,u8s,a,b,c)

              do j=0,nj
              do i=2,ni-2

                u8s=.5e0*(mf8u(i,j)*up(i,j,k)+mf8u(i+1,j)*up(i+1,j,k))

                if(u8s.gt.0.e0) then

                  a=dxt36w(i)*(sp(i+1,j,k)-3.e0*sp(i,j,k)               &
     &              +3.e0*sp(i-1,j,k)-sp(i-2,j,k))

                  b=dxt22(i)*(sp(i+1,j,k)                               &
     &              -2.e0*sp(i,j,k)+sp(i-1,j,k))

                  c=dxt16(i)*(3.e0*sp(i,j,k)+2.e0*sp(i+1,j,k)           &
     &              +sp(i-2,j,k)-6.e0*sp(i-1,j,k))

                  advx(i,j,k)=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                else

                  a=dxt36e(i)*(sp(i+2,j,k)-3.e0*sp(i+1,j,k)             &
     &              +3.e0*sp(i,j,k)-sp(i-1,j,k))

                  b=dxt22(i)*(sp(i+1,j,k)                               &
     &              -2.e0*sp(i,j,k)+sp(i-1,j,k))

                  c=dxt16(i)*(6.e0*sp(i+1,j,k)-sp(i+2,j,k)              &
     &              -2.e0*sp(i-1,j,k)-3.e0*sp(i,j,k))

                  advx(i,j,k)=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                end if

              end do
              end do

!$omp end do

            end do

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,v8s,a,b,c)

              do j=2,nj-2
              do i=2,ni-2

                v8s=.5e0*(vp(i,j,k)+vp(i,j+1,k))

                if(v8s.gt.0.e0) then

                  a=dyt36s(j)*(advx(i,j+1,k)-3.e0*advx(i,j,k)           &
     &              +3.e0*advx(i,j-1,k)-advx(i,j-2,k))

                  b=dyt22(j)*(advx(i,j+1,k)                             &
     &              -2.e0*advx(i,j,k)+advx(i,j-1,k))

                  c=dyt16(j)*(3.e0*advx(i,j,k)+2.e0*advx(i,j+1,k)       &
     &              +advx(i,j-2,k)-6.e0*advx(i,j-1,k))

                  sf(i,j,k)=advx(i,j,k)+((a*v8s+b)*v8s+c)*v8s

                else

                  a=dyt36n(j)*(advx(i,j+2,k)-3.e0*advx(i,j+1,k)         &
     &              +3.e0*advx(i,j,k)-advx(i,j-1,k))

                  b=dyt22(j)*(advx(i,j+1,k)                             &
     &              -2.e0*advx(i,j,k)+advx(i,j-1,k))

                  c=dyt16(j)*(6.e0*advx(i,j+1,k)-advx(i,j+2,k)          &
     &              -2.e0*advx(i,j-1,k)-3.e0*advx(i,j,k))

                  sf(i,j,k)=advx(i,j,k)+((a*v8s+b)*v8s+c)*v8s

                end if

              end do
              end do

!$omp end do

            end do

! -----

! Applied the map scale factor for the y direction.

          else if(mpopt.eq.5) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,u8s,a,b,c)

              do j=0,nj
              do i=2,ni-2

                u8s=.5e0*(up(i,j,k)+up(i+1,j,k))

                if(u8s.gt.0.e0) then

                  a=dxt36w(i)*(sp(i+1,j,k)-3.e0*sp(i,j,k)               &
     &              +3.e0*sp(i-1,j,k)-sp(i-2,j,k))

                  b=dxt22(i)*(sp(i+1,j,k)                               &
     &              -2.e0*sp(i,j,k)+sp(i-1,j,k))

                  c=dxt16(i)*(3.e0*sp(i,j,k)+2.e0*sp(i+1,j,k)           &
     &              +sp(i-2,j,k)-6.e0*sp(i-1,j,k))

                  advx(i,j,k)=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                else

                  a=dxt36e(i)*(sp(i+2,j,k)-3.e0*sp(i+1,j,k)             &
     &              +3.e0*sp(i,j,k)-sp(i-1,j,k))

                  b=dxt22(i)*(sp(i+1,j,k)                               &
     &              -2.e0*sp(i,j,k)+sp(i-1,j,k))

                  c=dxt16(i)*(6.e0*sp(i+1,j,k)-sp(i+2,j,k)              &
     &              -2.e0*sp(i-1,j,k)-3.e0*sp(i,j,k))

                  advx(i,j,k)=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                end if

              end do
              end do

!$omp end do

            end do

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,v8s,a,b,c)

              do j=2,nj-2
              do i=2,ni-2

                v8s=.5e0*(mf8v(i,j)*vp(i,j,k)+mf8v(i,j+1)*vp(i,j+1,k))

                if(v8s.gt.0.e0) then

                  a=dyt36s(j)*(advx(i,j+1,k)-3.e0*advx(i,j,k)           &
     &              +3.e0*advx(i,j-1,k)-advx(i,j-2,k))

                  b=dyt22(j)*(advx(i,j+1,k)                             &
     &              -2.e0*advx(i,j,k)+advx(i,j-1,k))

                  c=dyt16(j)*(3.e0*advx(i,j,k)+2.e0*advx(i,j+1,k)       &
     &              +advx(i,j-2,k)-6.e0*advx(i,j-1,k))

                  sf(i,j,k)=advx(i,j,k)+((a*v8s+b)*v8s+c)*v8s

                else

                  a=dyt36n(j)*(advx(i,j+2,k)-3.e0*advx(i,j+1,k)         &
     &              +3.e0*advx(i,j,k)-advx(i,j-1,k))

                  b=dyt22(j)*(advx(i,j+1,k)                             &
     &              -2.e0*advx(i,j,k)+advx(i,j-1,k))

                  c=dyt16(j)*(6.e0*advx(i,j+1,k)-advx(i,j+2,k)          &
     &              -2.e0*advx(i,j-1,k)-3.e0*advx(i,j,k))

                  sf(i,j,k)=advx(i,j,k)+((a*v8s+b)*v8s+c)*v8s

                end if

              end do
              end do

!$omp end do

            end do

! -----

! Applied the map scale factor for the x and y direction.

          else

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,u8s,a,b,c)

              do j=0,nj
              do i=2,ni-2

                u8s=.5e0*(mf8u(i,j)*up(i,j,k)+mf8u(i+1,j)*up(i+1,j,k))

                if(u8s.gt.0.e0) then

                  a=dxt36w(i)*(sp(i+1,j,k)-3.e0*sp(i,j,k)               &
     &              +3.e0*sp(i-1,j,k)-sp(i-2,j,k))

                  b=dxt22(i)*(sp(i+1,j,k)                               &
     &              -2.e0*sp(i,j,k)+sp(i-1,j,k))

                  c=dxt16(i)*(3.e0*sp(i,j,k)+2.e0*sp(i+1,j,k)           &
     &              +sp(i-2,j,k)-6.e0*sp(i-1,j,k))

                  advx(i,j,k)=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                else

                  a=dxt36e(i)*(sp(i+2,j,k)-3.e0*sp(i+1,j,k)             &
     &              +3.e0*sp(i,j,k)-sp(i-1,j,k))

                  b=dxt22(i)*(sp(i+1,j,k)                               &
     &              -2.e0*sp(i,j,k)+sp(i-1,j,k))

                  c=dxt16(i)*(6.e0*sp(i+1,j,k)-sp(i+2,j,k)              &
     &              -2.e0*sp(i-1,j,k)-3.e0*sp(i,j,k))

                  advx(i,j,k)=sp(i,j,k)+((a*u8s+b)*u8s+c)*u8s

                end if

              end do
              end do

!$omp end do

            end do

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,v8s,a,b,c)

              do j=2,nj-2
              do i=2,ni-2

                v8s=.5e0*(mf8v(i,j)*vp(i,j,k)+mf8v(i,j+1)*vp(i,j+1,k))

                if(v8s.gt.0.e0) then

                  a=dyt36s(j)*(advx(i,j+1,k)-3.e0*advx(i,j,k)           &
     &              +3.e0*advx(i,j-1,k)-advx(i,j-2,k))

                  b=dyt22(j)*(advx(i,j+1,k)                             &
     &              -2.e0*advx(i,j,k)+advx(i,j-1,k))

                  c=dyt16(j)*(3.e0*advx(i,j,k)+2.e0*advx(i,j+1,k)       &
     &              +advx(i,j-2,k)-6.e0*advx(i,j-1,k))

                  sf(i,j,k)=advx(i,j,k)+((a*v8s+b)*v8s+c)*v8s

                else

                  a=dyt36n(j)*(advx(i,j+2,k)-3.e0*advx(i,j+1,k)         &
     &              +3.e0*advx(i,j,k)-advx(i,j-1,k))

                  b=dyt22(j)*(advx(i,j+1,k)                             &
     &              -2.e0*advx(i,j,k)+advx(i,j-1,k))

                  c=dyt16(j)*(6.e0*advx(i,j+1,k)-advx(i,j+2,k)          &
     &              -2.e0*advx(i,j-1,k)-3.e0*advx(i,j,k))

                  sf(i,j,k)=advx(i,j,k)+((a*v8s+b)*v8s+c)*v8s

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

      end subroutine s_hculs

!-----7--------------------------------------------------------------7--

      end module m_hculs
