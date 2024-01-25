!***********************************************************************
      module m_eddyvis
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/07/05
!     Modification: 1999/07/21, 1999/09/30, 1999/10/12, 1999/11/01,
!                   1999/11/24, 2000/01/17, 2000/07/05, 2000/12/19,
!                   2001/02/13, 2001/04/15, 2001/05/29, 2001/10/18,
!                   2001/11/20, 2002/04/02, 2003/01/04, 2003/02/13,
!                   2003/03/21, 2003/04/30, 2003/05/19, 2003/10/31,
!                   2003/12/12, 2004/02/01, 2004/03/05, 2004/03/22,
!                   2004/04/15, 2004/05/31, 2004/07/01, 2004/09/01,
!                   2006/02/13, 2006/05/12, 2006/11/06, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/01/30, 2009/02/27,
!                   2009/11/13, 2011/08/09, 2011/11/10, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the eddy viscosity and the turbulent length scale.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: eddyvis, s_eddyvis

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface eddyvis

        module procedure s_eddyvis

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic exp
      intrinsic log
      intrinsic max
      intrinsic min
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_eddyvis(fpmpopt,fpmfcopt,fpsfcopt,fptubopt,fpisoopt, &
     &                     fpdx,fpdy,fpdz,dtb,ni,nj,nk,zph,jcb,rmf,rbr, &
     &                     tke,ssq,nsq8w,rkh,rkv,priv)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fptubopt
                       ! Formal parameter of unique index of tubopt

      integer, intent(in) :: fpisoopt
                       ! Formal parameter of unique index of isoopt

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

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: tke(0:ni+1,0:nj+1,1:nk)
                       ! Turbulent kinetic energy

      real, intent(in) :: ssq(0:ni+1,0:nj+1,1:nk)
                       ! Magnitude of deformation squared

      real, intent(in) :: nsq8w(0:ni+1,0:nj+1,1:nk)
                       ! Half value of Brunt Vaisara frequency squared
                       ! at w points

! Output variables

      real, intent(out) :: rkh(0:ni+1,0:nj+1,1:nk)
                       ! rbr x horizontal eddy viscosity

      real, intent(out) :: rkv(0:ni+1,0:nj+1,1:nk)
                       ! rbr x vertical eddy viscosity

      real, intent(out) :: priv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of turbulent Prandtl number

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor
      integer sfcopt   ! Option for surface physics
      integer tubopt   ! Option for turbulent mixing
      integer isoopt   ! Option for grid shape

      real dx          ! Grid distance in x direction
      real dy          ! Grid distance in y direction
      real dz          ! Grid distance in z direction

      real ds3         ! dx x dy x dz

      real cpriv       ! Inverse of prnum

      real lnh         ! Horizontal turbulent length scale

      real ckmax       ! Coefficient of upper limit of eddy viscosity

      real khmin       ! Minimum horizontal eddy viscosity
      real khmax       ! Maximum horizontal eddy viscosity

      real csnum2      ! csnum x csnum
      real cslnh2      ! csnum x csnum x lnh x lnh

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real ln          ! Turbulent length scale

      real ln0         ! Reference turbulent length scale

      real ln2         ! ln x ln

      real ln02        ! ln0 x ln0

      real htskp       ! kappa x terrain height from surface

      real nsq         ! Brunt Vaisara frequency squared

      real stabc       ! Stability cheker

      real a           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getiname(fpsfcopt,sfcopt)
      call getiname(fptubopt,tubopt)
      call getiname(fpisoopt,isoopt)
      call getrname(fpdx,dx)
      call getrname(fpdy,dy)
      call getrname(fpdz,dz)

! -----

! Set the common used variables.

      ds3=dx*dy*dz

      cpriv=1.e0/prnum

      if(cpriv.gt.3.e0) then
        cpriv=3.e0
      end if

      lnh=sqrt(dx*dy)

      ckmax=.125e0/dtb

      khmin=dx*dy*ckmin
      khmax=dx*dy*ckmax

      csnum2=csnum*csnum
      cslnh2=csnum*csnum*lnh*lnh

! -----

!!!! Calculate the eddy viscosity.

!$omp parallel default(shared) private(k)

!! Calculate the eddy viscosity with the Smagorinsky formulation.

      if(tubopt.eq.1) then

! Isotropic case.

        if(isoopt.eq.1) then

          if(mfcopt.eq.0) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln2)

              do j=1,nj-1
              do i=1,ni-1
                ln=exp(oned3*log(ds3*jcb(i,j,k)))

                ln2=ln*ln

                priv(i,j,k)=cpriv

                rkv(i,j,k)=csnum2*ln2*sqrt(max(ssq(i,j,k)               &
     &            -(nsq8w(i,j,k)+nsq8w(i,j,k+1))*cpriv,0.e0))

                rkv(i,j,k)=rbr(i,j,k)*min(rkv(i,j,k),ckmax*ln2)
                rkh(i,j,k)=rkv(i,j,k)

              end do
              end do

!$omp end do

            end do

          else

            if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln2)

                do j=1,nj-1
                do i=1,ni-1
                  ln=exp(oned3*log(ds3*rmf(i,j,2)*jcb(i,j,k)))

                  ln2=ln*ln

                  priv(i,j,k)=cpriv

                  rkv(i,j,k)=csnum2*ln2*sqrt(max(ssq(i,j,k)             &
     &              -(nsq8w(i,j,k)+nsq8w(i,j,k+1))*cpriv,0.e0))

                  rkv(i,j,k)=rbr(i,j,k)*min(rkv(i,j,k),ckmax*ln2)
                  rkh(i,j,k)=rkv(i,j,k)

                end do
                end do

!$omp end do

              end do

            else

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln2)

                do j=1,nj-1
                do i=1,ni-1
                  ln=exp(oned3*log(ds3*rmf(i,j,3)*jcb(i,j,k)))

                  ln2=ln*ln

                  priv(i,j,k)=cpriv

                  rkv(i,j,k)=csnum2*ln2*sqrt(max(ssq(i,j,k)             &
     &              -(nsq8w(i,j,k)+nsq8w(i,j,k+1))*cpriv,0.e0))

                  rkv(i,j,k)=rbr(i,j,k)*min(rkv(i,j,k),ckmax*ln2)
                  rkh(i,j,k)=rkv(i,j,k)

                end do
                end do

!$omp end do

              end do

            end if

          end if

! -----

! Anisotropic case.

        else if(isoopt.eq.2) then

          if(mfcopt.eq.0) then

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln2,stabc)

              do j=1,nj-1
              do i=1,ni-1
                ln=zph(i,j,k+1)-zph(i,j,k)

                ln2=ln*ln

                priv(i,j,k)=cpriv

                stabc=sqrt(max(ssq(i,j,k)                               &
     &            -(nsq8w(i,j,k)+nsq8w(i,j,k+1))*cpriv,0.e0))

                rkv(i,j,k)=rbr(i,j,k)*min(stabc*csnum2,ckmax)*ln2
                rkh(i,j,k)=rbr(i,j,k)*min(stabc*cslnh2,khmax)

              end do
              end do

!$omp end do

            end do

          else

            if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln2,stabc)

                do j=1,nj-1
                do i=1,ni-1
                  ln=zph(i,j,k+1)-zph(i,j,k)

                  ln2=ln*ln

                  priv(i,j,k)=cpriv

                  stabc=sqrt(max(ssq(i,j,k)                             &
     &              -(nsq8w(i,j,k)+nsq8w(i,j,k+1))*cpriv,0.e0))

                  rkv(i,j,k)=rbr(i,j,k)                                 &
     &              *min(stabc*csnum2,ckmax)*ln2

                  rkh(i,j,k)=rbr(i,j,k)                                 &
     &              *min(stabc*cslnh2,khmax)*rmf(i,j,2)

                end do
                end do

!$omp end do

              end do

            else

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln2,stabc)

                do j=1,nj-1
                do i=1,ni-1
                  ln=zph(i,j,k+1)-zph(i,j,k)

                  ln2=ln*ln

                  priv(i,j,k)=cpriv

                  stabc=sqrt(max(ssq(i,j,k)                             &
     &              -(nsq8w(i,j,k)+nsq8w(i,j,k+1))*cpriv,0.e0))

                  rkv(i,j,k)=rbr(i,j,k)                                 &
     &              *min(stabc*csnum2,ckmax)*ln2

                  rkh(i,j,k)=rbr(i,j,k)                                 &
     &              *min(stabc*cslnh2,khmax)*rmf(i,j,3)

                end do
                end do

!$omp end do

              end do

            end if

          end if

        end if

! -----

!! -----

!!! Calculate the eddy viscosity with the Deardorff formulation.

      else

!! In the case of no surface process.

        if(sfcopt.eq.0) then

! Isotropic case.

          if(isoopt.eq.1) then

            if(mfcopt.eq.0) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln0,ln02,nsq)

                do j=1,nj-1
                do i=1,ni-1
                  ln0=exp(oned3*log(ds3*jcb(i,j,k)))

                  ln02=ln0*ln0

                  nsq=nsq8w(i,j,k)+nsq8w(i,j,k+1)

                  if(nsq.lt.0.e0) then

                    ln=ln0

                  else

                    ln=max(.1e0*ln0,                                    &
     &                 min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)),ln0))

                  end if

                  priv(i,j,k)=1.e0+2.e0*ln/ln0

                  rkv(i,j,k)=ckm*ln*sqrt(tke(i,j,k))

                  if(priv(i,j,k)*nsq.lt.ssq(i,j,k)) then

                    rkv(i,j,k)=max(rkv(i,j,k),ckmin*ln02)

                  end if

                  rkv(i,j,k)=rbr(i,j,k)*min(rkv(i,j,k),ckmax*ln02)
                  rkh(i,j,k)=rkv(i,j,k)

                end do
                end do

!$omp end do

              end do

            else

              if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln0,ln02,nsq)

                  do j=1,nj-1
                  do i=1,ni-1
                    ln0=exp(oned3*log(ds3*rmf(i,j,2)*jcb(i,j,k)))

                    ln02=ln0*ln0

                    nsq=nsq8w(i,j,k)+nsq8w(i,j,k+1)

                    if(nsq.lt.0.e0) then

                      ln=ln0

                    else

                      ln=max(.1e0*ln0,                                  &
     &                   min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)),ln0))

                    end if

                    priv(i,j,k)=1.e0+2.e0*ln/ln0

                    rkv(i,j,k)=ckm*ln*sqrt(tke(i,j,k))

                    if(priv(i,j,k)*nsq.lt.ssq(i,j,k)) then

                      rkv(i,j,k)=max(rkv(i,j,k),ckmin*ln02)

                    end if

                    rkv(i,j,k)=rbr(i,j,k)*min(rkv(i,j,k),ckmax*ln02)
                    rkh(i,j,k)=rkv(i,j,k)

                  end do
                  end do

!$omp end do

                end do

              else

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln0,ln02,nsq)

                  do j=1,nj-1
                  do i=1,ni-1
                    ln0=exp(oned3*log(ds3*rmf(i,j,3)*jcb(i,j,k)))

                    ln02=ln0*ln0

                    nsq=nsq8w(i,j,k)+nsq8w(i,j,k+1)

                    if(nsq.lt.0.e0) then

                      ln=ln0

                    else

                      ln=max(.1e0*ln0,                                  &
     &                   min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)),ln0))

                    end if

                    priv(i,j,k)=1.e0+2.e0*ln/ln0

                    rkv(i,j,k)=ckm*ln*sqrt(tke(i,j,k))

                    if(priv(i,j,k)*nsq.lt.ssq(i,j,k)) then

                      rkv(i,j,k)=max(rkv(i,j,k),ckmin*ln02)

                    end if

                    rkv(i,j,k)=rbr(i,j,k)*min(rkv(i,j,k),ckmax*ln02)
                    rkh(i,j,k)=rkv(i,j,k)

                  end do
                  end do

!$omp end do

                end do

              end if

            end if

! -----

! Anisotropic case.

          else if(isoopt.eq.2) then

            if(mfcopt.eq.0) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln0,ln02,nsq,a)

                do j=1,nj-1
                do i=1,ni-1
                  ln0=zph(i,j,k+1)-zph(i,j,k)

                  ln02=ln0*ln0

                  nsq=nsq8w(i,j,k)+nsq8w(i,j,k+1)

                  if(nsq.lt.0.e0) then

                    ln=ln0

                  else

                    ln=max(.1e0*ln0,                                    &
     &                 min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)),ln0))

                  end if

                  priv(i,j,k)=1.e0+2.e0*ln/ln0

                  a=ckm*sqrt(tke(i,j,k))

                  rkv(i,j,k)=a*ln
                  rkh(i,j,k)=a*lnh

                  if(priv(i,j,k)*nsq.lt.ssq(i,j,k)) then

                    rkv(i,j,k)=max(rkv(i,j,k),ckmin*ln02)
                    rkh(i,j,k)=max(rkh(i,j,k),khmin)

                  end if

                  rkv(i,j,k)=rbr(i,j,k)*min(rkv(i,j,k),ckmax*ln02)
                  rkh(i,j,k)=rbr(i,j,k)*min(rkh(i,j,k),khmax)

                end do
                end do

!$omp end do

              end do

            else

              if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln0,ln02,nsq,a)

                  do j=1,nj-1
                  do i=1,ni-1
                    ln0=zph(i,j,k+1)-zph(i,j,k)

                    ln02=ln0*ln0

                    nsq=nsq8w(i,j,k)+nsq8w(i,j,k+1)

                    if(nsq.lt.0.e0) then

                      ln=ln0

                    else

                      ln=max(.1e0*ln0,                                  &
     &                   min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)),ln0))

                    end if

                    priv(i,j,k)=1.e0+2.e0*ln/ln0

                    a=ckm*sqrt(tke(i,j,k))

                    rkv(i,j,k)=a*ln
                    rkh(i,j,k)=a*lnh*rmf(i,j,4)

                    if(priv(i,j,k)*nsq.lt.ssq(i,j,k)) then

                      rkv(i,j,k)=max(rkv(i,j,k),ckmin*ln02)
                      rkh(i,j,k)=max(rkh(i,j,k),khmin*rmf(i,j,2))

                    end if

                    rkv(i,j,k)=rbr(i,j,k)                               &
     &                *min(rkv(i,j,k),ckmax*ln02)

                    rkh(i,j,k)=rbr(i,j,k)                               &
     &                *min(rkh(i,j,k),khmax*rmf(i,j,2))

                  end do
                  end do

!$omp end do

                end do

              else

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln0,ln02,nsq,a)

                  do j=1,nj-1
                  do i=1,ni-1
                    ln0=zph(i,j,k+1)-zph(i,j,k)

                    ln02=ln0*ln0

                    nsq=nsq8w(i,j,k)+nsq8w(i,j,k+1)

                    if(nsq.lt.0.e0) then

                      ln=ln0

                    else

                      ln=max(.1e0*ln0,                                  &
     &                   min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)),ln0))

                    end if

                    priv(i,j,k)=1.e0+2.e0*ln/ln0

                    a=ckm*sqrt(tke(i,j,k))

                    rkv(i,j,k)=a*ln
                    rkh(i,j,k)=a*lnh*rmf(i,j,2)

                    if(priv(i,j,k)*nsq.lt.ssq(i,j,k)) then

                      rkv(i,j,k)=max(rkv(i,j,k),ckmin*ln02)
                      rkh(i,j,k)=max(rkh(i,j,k),khmin*rmf(i,j,3))

                    end if

                    rkv(i,j,k)=rbr(i,j,k)                               &
     &                *min(rkv(i,j,k),ckmax*ln02)

                    rkh(i,j,k)=rbr(i,j,k)                               &
     &                *min(rkh(i,j,k),khmax*rmf(i,j,3))

                  end do
                  end do

!$omp end do

                end do

              end if

            end if

          end if

! -----

!! -----

!! In the case of performing surface process.

        else

! Isotropic case.

          if(isoopt.eq.1) then

            if(mfcopt.eq.0) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln0,ln02,htskp,nsq)

                do j=1,nj-1
                do i=1,ni-1
                  ln0=exp(oned3*log(ds3*jcb(i,j,k)))

                  ln02=ln0*ln0

                  htskp=kappa                                           &
     &              *abs(.5e0*(zph(i,j,k)+zph(i,j,k+1))-zph(i,j,2))

                  nsq=nsq8w(i,j,k)+nsq8w(i,j,k+1)

                  if(nsq.lt.0.e0) then

                    ln=ln0

                  else

                    ln=max(.1e0*ln0,                                    &
     &                 min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)),ln0))

                  end if

                  ln=ln*htskp/(htskp+ln)

                  priv(i,j,k)=1.e0+2.e0*ln/ln0

                  rkv(i,j,k)=ckm*ln*sqrt(tke(i,j,k))

                  if(priv(i,j,k)*nsq.lt.ssq(i,j,k)) then

                    rkv(i,j,k)=max(rkv(i,j,k),ckmin*ln02)

                  end if

                  rkv(i,j,k)=rbr(i,j,k)*min(rkv(i,j,k),ckmax*ln02)
                  rkh(i,j,k)=rkv(i,j,k)

                end do
                end do

!$omp end do

              end do

            else

              if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln0,ln02,htskp,nsq)

                  do j=1,nj-1
                  do i=1,ni-1
                    ln0=exp(oned3*log(ds3*rmf(i,j,2)*jcb(i,j,k)))

                    ln02=ln0*ln0

                    htskp=kappa                                         &
     &                *abs(.5e0*(zph(i,j,k)+zph(i,j,k+1))-zph(i,j,2))

                    nsq=nsq8w(i,j,k)+nsq8w(i,j,k+1)

                    if(nsq.lt.0.e0) then

                      ln=ln0

                    else

                      ln=max(.1e0*ln0,                                  &
     &                   min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)),ln0))

                    end if

                    ln=ln*htskp/(htskp+ln)

                    priv(i,j,k)=1.e0+2.e0*ln/ln0

                    rkv(i,j,k)=ckm*ln*sqrt(tke(i,j,k))

                    if(priv(i,j,k)*nsq.lt.ssq(i,j,k)) then

                      rkv(i,j,k)=max(rkv(i,j,k),ckmin*ln02)

                    end if

                    rkv(i,j,k)=rbr(i,j,k)*min(rkv(i,j,k),ckmax*ln02)
                    rkh(i,j,k)=rkv(i,j,k)

                  end do
                  end do

!$omp end do

                end do

              else

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln0,ln02,htskp,nsq)

                  do j=1,nj-1
                  do i=1,ni-1
                    ln0=exp(oned3*log(ds3*rmf(i,j,3)*jcb(i,j,k)))

                    ln02=ln0*ln0

                    htskp=kappa                                         &
     &                *abs(.5e0*(zph(i,j,k)+zph(i,j,k+1))-zph(i,j,2))

                    nsq=nsq8w(i,j,k)+nsq8w(i,j,k+1)

                    if(nsq.lt.0.e0) then

                      ln=ln0

                    else

                      ln=max(.1e0*ln0,                                  &
     &                   min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)),ln0))

                    end if

                    ln=ln*htskp/(htskp+ln)

                    priv(i,j,k)=1.e0+2.e0*ln/ln0

                    rkv(i,j,k)=ckm*ln*sqrt(tke(i,j,k))

                    if(priv(i,j,k)*nsq.lt.ssq(i,j,k)) then

                      rkv(i,j,k)=max(rkv(i,j,k),ckmin*ln02)

                    end if

                    rkv(i,j,k)=rbr(i,j,k)*min(rkv(i,j,k),ckmax*ln02)
                    rkh(i,j,k)=rkv(i,j,k)

                  end do
                  end do

!$omp end do

                end do

              end if

            end if

! -----

! Anisotropic case.

          else if(isoopt.eq.2) then

            if(mfcopt.eq.0) then

              do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln0,ln02,htskp,nsq,a)

                do j=1,nj-1
                do i=1,ni-1
                  ln0=zph(i,j,k+1)-zph(i,j,k)

                  ln02=ln0*ln0

                  htskp=kappa                                           &
     &              *abs(.5e0*(zph(i,j,k)+zph(i,j,k+1))-zph(i,j,2))

                  nsq=nsq8w(i,j,k)+nsq8w(i,j,k+1)

                  if(nsq.lt.0.e0) then

                    ln=ln0

                  else

                    ln=max(.1e0*ln0,                                    &
     &                 min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)),ln0))

                  end if

                  ln=ln*htskp/(htskp+ln)

                  priv(i,j,k)=1.e0+2.e0*ln/ln0

                  a=ckm*sqrt(tke(i,j,k))

                  rkv(i,j,k)=a*ln
                  rkh(i,j,k)=a*lnh

                  if(priv(i,j,k)*nsq.lt.ssq(i,j,k)) then

                    rkv(i,j,k)=max(rkv(i,j,k),ckmin*ln02)
                    rkh(i,j,k)=max(rkh(i,j,k),khmin)

                  end if

                  rkv(i,j,k)=rbr(i,j,k)*min(rkv(i,j,k),ckmax*ln02)
                  rkh(i,j,k)=rbr(i,j,k)*min(rkh(i,j,k),khmax)

                end do
                end do

!$omp end do

              end do

            else

              if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln0,ln02,htskp,nsq,a)

                  do j=1,nj-1
                  do i=1,ni-1
                    ln0=zph(i,j,k+1)-zph(i,j,k)

                    ln02=ln0*ln0

                    htskp=kappa                                         &
     &                *abs(.5e0*(zph(i,j,k)+zph(i,j,k+1))-zph(i,j,2))

                    nsq=nsq8w(i,j,k)+nsq8w(i,j,k+1)

                    if(nsq.lt.0.e0) then

                      ln=ln0

                    else

                      ln=max(.1e0*ln0,                                  &
     &                   min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)),ln0))

                    end if

                    ln=ln*htskp/(htskp+ln)

                    priv(i,j,k)=1.e0+2.e0*ln/ln0

                    a=ckm*sqrt(tke(i,j,k))

                    rkv(i,j,k)=a*ln
                    rkh(i,j,k)=a*lnh*rmf(i,j,4)

                    if(priv(i,j,k)*nsq.lt.ssq(i,j,k)) then

                      rkv(i,j,k)=max(rkv(i,j,k),ckmin*ln02)
                      rkh(i,j,k)=max(rkh(i,j,k),khmin*rmf(i,j,2))

                    end if

                    rkv(i,j,k)=rbr(i,j,k)                               &
     &                *min(rkv(i,j,k),ckmax*ln02)

                    rkh(i,j,k)=rbr(i,j,k)                               &
     &                *min(rkh(i,j,k),khmax*rmf(i,j,2))

                  end do
                  end do

!$omp end do

                end do

              else

                do k=1,nk-1

!$omp do schedule(runtime) private(i,j,ln,ln0,ln02,htskp,nsq,a)

                  do j=1,nj-1
                  do i=1,ni-1
                    ln0=zph(i,j,k+1)-zph(i,j,k)

                    ln02=ln0*ln0

                    htskp=kappa                                         &
     &                *abs(.5e0*(zph(i,j,k)+zph(i,j,k+1))-zph(i,j,2))

                    nsq=nsq8w(i,j,k)+nsq8w(i,j,k+1)

                    if(nsq.lt.0.e0) then

                      ln=ln0

                    else

                      ln=max(.1e0*ln0,                                  &
     &                   min(.76e0*sqrt(tke(i,j,k)/(nsq+eps)),ln0))

                    end if

                    ln=ln*htskp/(htskp+ln)

                    priv(i,j,k)=1.e0+2.e0*ln/ln0

                    a=ckm*sqrt(tke(i,j,k))

                    rkv(i,j,k)=a*ln
                    rkh(i,j,k)=a*lnh*rmf(i,j,2)

                    if(priv(i,j,k)*nsq.lt.ssq(i,j,k)) then

                      rkv(i,j,k)=max(rkv(i,j,k),ckmin*ln02)
                      rkh(i,j,k)=max(rkh(i,j,k),khmin*rmf(i,j,3))

                    end if

                    rkv(i,j,k)=rbr(i,j,k)                               &
     &                *min(rkv(i,j,k),ckmax*ln02)

                    rkh(i,j,k)=rbr(i,j,k)                               &
     &                *min(rkh(i,j,k),khmax*rmf(i,j,3))

                  end do
                  end do

!$omp end do

                end do

              end if

            end if

          end if

! -----

        end if

!! -----

      end if

!!! -----

!$omp end parallel

!!!! -----

      end subroutine s_eddyvis

!-----7--------------------------------------------------------------7--

      end module m_eddyvis
