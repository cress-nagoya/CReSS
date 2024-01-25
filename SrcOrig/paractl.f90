!***********************************************************************
      module m_paractl
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/04/11
!     Modification: 2007/08/24, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2009/11/13, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the parameters of GrADS control file.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_commath
      use m_commpi
      use m_currpe
      use m_getiname
      use m_getrname
      use m_llnews
      use m_setproj

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: paractl, s_paractl

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface paractl

        module procedure s_paractl

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic atan
      intrinsic exp
      intrinsic int
      intrinsic max
      intrinsic min
      intrinsic mod
      intrinsic real
      intrinsic tan

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_paractl(fpmpopt,fpnspol,fpuniopt_uni,                &
     &                     fpdx,fpdy,fpdxiv,fpdyiv,ni,nj,               &
     &                     latsw,lonsw,latmin,lonmin,dlat,dlon,         &
     &                     ipole,jpole,mlat)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpnspol
                       ! Formal parameter of unique index of nspol

      integer, intent(in) :: fpuniopt_uni
                       ! Formal parameter of unique index of uniopt_uni

      integer, intent(in) :: fpdx
                       ! Formal parameter of unique index of dx

      integer, intent(in) :: fpdy
                       ! Formal parameter of unique index of dy

      integer, intent(in) :: fpdxiv
                       ! Formal parameter of unique index of dxiv

      integer, intent(in) :: fpdyiv
                       ! Formal parameter of unique index of dyiv

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

! Output variables

      real, intent(out) :: latsw
                       ! Latitude at south-west corner

      real, intent(out) :: lonsw
                       ! Longitude at south-west corner

      real, intent(out) :: latmin
                       ! Minimum latitude in model domain

      real, intent(out) :: lonmin
                       ! Minimum longitude in model domain

      real, intent(out) :: dlat
                       ! Grid distance in latitude direction

      real, intent(out) :: dlon
                       ! Grid distance in longitude direction

      real, intent(out) :: ipole
                       ! Reference real index at pole in x direction

      real, intent(out) :: jpole
                       ! Reference real index at pole in y direction

      real, intent(out) :: mlat(1:(nj-3)*njgrp*njsub)
                       ! Latitude with Mercator projection

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer nspol    ! Option for projected region

      integer uniopt_uni
                       ! Option for uniting process

      integer inpole   ! Control flag of position of pole

      real dx          ! Grid distance in x direction
      real dy          ! Grid distance in y direction

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy

      real x0          ! x origin
      real y0          ! y origin

      real cpj(1:7)    ! Map projection parameters

      real latmax      ! Maximum latitude in model domain
      real lonmax      ! Maximum longitude in model domain

      real x4          ! x coordinates at model corners
      real y4          ! y coordinates at model corners

      real lat4(1:5)   ! Latitude at model corners
      real lon4(1:5)   ! Longitude at model corners

! Internal private variables

      integer j        ! Array index in y direction

      integer ic       ! Index of do loop

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpnspol,nspol)
      call getiname(fpuniopt_uni,uniopt_uni)
      call getrname(fpdx,dx)
      call getrname(fpdy,dy)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)

! -----

! Get the map projection parameters.

      call setproj(idmpopt,idnspol,iddx,iddy,idulat,idulon,idriu,idrju, &
     &             idtlat1,idtlat2,idtlon,'solver  ',6,x0,y0,cpj)

! -----

!! Get the latitude and longitude at south-west corner.

! Set the parameters of parallelizing.

      if(abs(uniopt_uni).eq.1.or.abs(uniopt_uni).eq.3) then

        mygrp=0

        call currpe('unite   ',5,'ijgrp')

      else if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

        myred=0

        call currpe('unite   ',5,'ijred')

      else if(ngrp.eq.1.and.abs(uniopt_uni).ge.5) then

        mygrp=0

        call currpe('unite   ',5,'ijgrp')

      end if

      x4=.5e0*real(2*(ni-3)*nisub*igrp+1)*dx
      y4=.5e0*real(2*(nj-3)*njsub*jgrp+1)*dy

! -----

! Calculate the latitude and the longitude at south-west corner.

      call llnews(idmpopt,idnspol,idtlon,x0,y0,cpj,x4,y4,               &
     &            lat4(1),lon4(1))

! -----

!! -----

!! Get the latitude and longitude at south-east corner.

! Set the parameters of parallelizing.

      if(abs(uniopt_uni).eq.1.or.abs(uniopt_uni).eq.3) then

        mygrp=nigrp-1

        call currpe('unite   ',5,'ijgrp')

      else if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

        myred=nired-1

        call currpe('unite   ',5,'ijred')

      else if(ngrp.eq.1.and.abs(uniopt_uni).ge.5) then

        mygrp=nigrp-1

        call currpe('unite   ',5,'ijgrp')

      end if

      x4=.5e0*real(2*(ni-3)*nisub*(igrp+1)-1)*dx
      y4=.5e0*real(2*(nj-3)*njsub*jgrp+1)*dy

! -----

! Calculate the latitude and the longitude at south-west corner.

      call llnews(idmpopt,idnspol,idtlon,x0,y0,cpj,x4,y4,               &
     &            lat4(2),lon4(2))

! -----

!! -----

!! Get the latitude and longitude at north-west corner.

! Set the parameters of parallelizing.

      if(abs(uniopt_uni).eq.1.or.abs(uniopt_uni).eq.3) then

        mygrp=nigrp*(njgrp-1)

        call currpe('unite   ',5,'ijgrp')

      else if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

        myred=nired*(njred-1)

        call currpe('unite   ',5,'ijred')

      else if(ngrp.eq.1.and.abs(uniopt_uni).ge.5) then

        mygrp=nigrp*(njgrp-1)

        call currpe('unite   ',5,'ijgrp')

      end if

      x4=.5e0*real(2*(ni-3)*nisub*igrp+1)*dx
      y4=.5e0*real(2*(nj-3)*njsub*(jgrp+1)-1)*dy

! -----

! Calculate the latitude and the longitude at south-west corner.

      call llnews(idmpopt,idnspol,idtlon,x0,y0,cpj,x4,y4,               &
     &            lat4(3),lon4(3))

! -----

!! -----

!! Get the latitude and longitude at north-east corner.

! Set the parameters of parallelizing.

      if(abs(uniopt_uni).eq.1.or.abs(uniopt_uni).eq.3) then

        mygrp=ngrp-1

        call currpe('unite   ',5,'ijgrp')

      else if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

        myred=nred-1

        call currpe('unite   ',5,'ijred')

      else if(ngrp.eq.1.and.abs(uniopt_uni).ge.5) then

        mygrp=ngrp-1

        call currpe('unite   ',5,'ijgrp')

      end if

      x4=.5e0*real(2*(ni-3)*nisub*(igrp+1)-1)*dx
      y4=.5e0*real(2*(nj-3)*njsub*(jgrp+1)-1)*dy

! -----

! Calculate the latitude and the longitude at south-west corner.

      call llnews(idmpopt,idnspol,idtlon,x0,y0,cpj,x4,y4,               &
     &            lat4(4),lon4(4))

! -----

!! -----

!! Get the latitude and longitude at true longitude point.

      if(mpopt.eq.1.or.mpopt.eq.2) then

! Set the parameters of parallelizing.

        igrp=int((-x0*dxiv+.5e0)/real((ni-3)*nisub))

        if(igrp.ge.0.and.igrp.le.nigrp-1) then

          if(nspol.eq.1) then

            if(abs(uniopt_uni).eq.1.or.abs(uniopt_uni).eq.3) then

              mygrp=nigrp*(njgrp-1)+igrp

              call currpe('unite   ',5,'ijgrp')

            else if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

              ired=igrp-iwred

              myred=nired*(njred-1)+ired

              call currpe('unite   ',5,'ijred')

            else if(ngrp.eq.1.and.abs(uniopt_uni).ge.5) then

              mygrp=nigrp*(njgrp-1)+igrp

              call currpe('unite   ',5,'ijgrp')

            end if

            x4=.5e0*real(2*(ni-3)*nisub*igrp                            &
     &        +2*mod(int(-x0*dxiv+.5e0),((ni-3)*nisub))+1)*dx

            y4=.5e0*real(2*(nj-3)*njsub*(jgrp+1)-1)*dy

          else

            if(abs(uniopt_uni).eq.1.or.abs(uniopt_uni).eq.3) then

              mygrp=igrp

              call currpe('unite   ',5,'ijgrp')

            else if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

              myred=igrp-iwred

              call currpe('unite   ',5,'ijred')

            else if(ngrp.eq.1.and.abs(uniopt_uni).ge.5) then

              mygrp=igrp

              call currpe('unite   ',5,'ijgrp')

            end if

            x4=.5e0*real(2*(ni-3)*nisub*igrp                            &
     &        +2*mod(int(-x0*dxiv+.5e0),((ni-3)*nisub))+1)*dx

            y4=.5e0*real(2*(nj-3)*njsub*jgrp+1)*dy

          end if

        else

          if(abs(uniopt_uni).eq.1.or.abs(uniopt_uni).eq.3) then

            mygrp=ngrp-1

            call currpe('unite   ',5,'ijgrp')

          else if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

            myred=nred-1

            call currpe('unite   ',5,'ijred')

          else if(ngrp.eq.1.and.abs(uniopt_uni).ge.5) then

            mygrp=ngrp-1

            call currpe('unite   ',5,'ijgrp')

          end if

          x4=.5e0*real(2*(ni-3)*nisub*(igrp+1)-1)*dx
          y4=.5e0*real(2*(nj-3)*njsub*(jgrp+1)-1)*dy

        end if

! -----

! Calculate the latitude and the longitude at south-west corner.

        call llnews(idmpopt,idnspol,idtlon,x0,y0,cpj,x4,y4,             &
     &              lat4(5),lon4(5))

! -----

      end if

!! -----

!! Get the minimum latitude and longitude.

! Initialize the processed variables.

      latmin=100.e0
      latmax=-100.e0

      lonmin=200.e0
      lonmax=-200.e0

! -----

! Calculate the minimum latitude and longitude.

!$omp parallel default(shared)

      if(mpopt.eq.1.or.mpopt.eq.2) then

!$omp do schedule(runtime) private(ic)                                  &
!$omp&   reduction(min: latmin,lonmin) reduction(max: latmax,lonmax)

        do ic=1,5

          latmin=min(lat4(ic),latmin)
          latmax=max(lat4(ic),latmax)

          lonmin=min(lon4(ic),lonmin)
          lonmax=max(lon4(ic),lonmax)

        end do

!$omp end do

      else

!$omp do schedule(runtime) private(ic)                                  &
!$omp&   reduction(min: latmin,lonmin) reduction(max: latmax,lonmax)

        do ic=1,4

          latmin=min(lat4(ic),latmin)
          latmax=max(lat4(ic),latmax)

          lonmin=min(lon4(ic),lonmin)
          lonmax=max(lon4(ic),lonmax)

        end do

!$omp end do

      end if

!$omp end parallel

! -----

! Set the latitude and longitude at south-west corner.

      latsw=lat4(1)
      lonsw=lon4(1)

! -----

!! -----

! Get the reference indices at pole.

      if(mpopt.eq.1.or.mpopt.eq.2) then

        if(nspol.eq.1) then
          ipole=-x0*dxiv+.5e0
          jpole=-y0*dyiv+.5e0
        else
          ipole=-x0*dxiv+.5e0
          jpole=y0*dyiv+.5e0
        end if

        inpole=0

        if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

          ipole=ipole-(ni-3)*nisub*iwred
          jpole=jpole-(nj-3)*njsub*jsred

          if((ipole.ge.1.and.ipole.le.(ni-3)*nired*nisub)               &
     &      .and.(jpole.ge.1.and.jpole.le.(nj-3)*njred*njsub)) then

            inpole=1

          end if

        else

          if((ipole.ge.1.and.ipole.le.(ni-3)*nigrp*nisub)               &
     &      .and.(jpole.ge.1.and.jpole.le.(nj-3)*njgrp*njsub)) then

            inpole=1

          end if

        end if

      else

        ipole=1.e0
        jpole=1.e0

        inpole=0

      end if

! -----

! Get the grid distance in latitude and longitude direction.

      dlon=abs(lonmax-lonmin)

      if(inpole.eq.1) then

        if(nspol.eq.1) then
          dlat=90.e0-latmin
        else
          dlat=90.e0-abs(latmax)
        end if

        if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

          if((nj-3)*njred*njsub.gt.1) then
            dlat=dlat/real((nj-3)*njred*njsub+1)

            latmin=latmin-2.e0*dlat

          else
            dlat=1.e0

          end if

          if((ni-3)*nired*nisub.gt.1) then
            dlon=360.e0/real((ni-3)*nired*nisub+2)

          else
            dlon=1.e0

          end if

        else

          if((nj-3)*njgrp*njsub.gt.1) then
            dlat=dlat/real((nj-3)*njgrp*njsub+1)

            latmin=latmin-2.e0*dlat

          else
            dlat=1.e0

          end if

          if((ni-3)*nigrp*nisub.gt.1) then
            dlon=360.e0/real((ni-3)*nigrp*nisub+2)

          else
            dlon=1.e0

          end if

        end if

      else

        dlat=latmax-latmin

        if(lon4(2).lt.lon4(1).or.lon4(4).lt.lon4(3)) then

          if(dlon.gt.180.e0) then
            dlon=abs(dlon-360.e0)
          end if

        end if

        if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

          if((nj-3)*njred*njsub.gt.1) then
            dlat=dlat/real((nj-3)*njred*njsub-1)

            if(mpopt.eq.1.or.mpopt.eq.2) then
              latmin=latmin-dlat
            end if

          else
            dlat=1.e0

          end if

          if((ni-3)*nired*nisub.gt.1) then
            dlon=dlon/real((ni-3)*nired*nisub-1)

            if(mpopt.eq.1.or.mpopt.eq.2) then
              lonmin=lonmin-dlon
            end if

          else
            dlon=1.e0

          end if

        else

          if((nj-3)*njgrp*njsub.gt.1) then
            dlat=dlat/real((nj-3)*njgrp*njsub-1)

            if(mpopt.eq.1.or.mpopt.eq.2) then
              latmin=latmin-dlat
            end if

          else
            dlat=1.e0

          end if

          if((ni-3)*nigrp*nisub.gt.1) then
            dlon=dlon/real((ni-3)*nigrp*nisub-1)

            if(mpopt.eq.1.or.mpopt.eq.2) then
              lonmin=lonmin-dlon
            end if

          else
            dlon=1.e0

          end if

        end if

      end if

! -----

! Calculate the latitude with the Mercator projection method.

!$omp parallel default(shared)

      if(mpopt.eq.3.or.mpopt.eq.13) then

        if(abs(uniopt_uni).eq.2.or.abs(uniopt_uni).eq.4) then

!$omp do schedule(runtime) private(j)

          do j=(nj-3)*njsub*jsred+1,(nj-3)*njsub*(jsred+njred)

            mlat(j)=90.e0                                               &
     &        -2.e0*atan(exp((.5e0*real(1-2*j)*dy-y0)*cpj(3)))*r2d

            mlat(j)=max(min(mlat(j),90.e0),-90.e0)

          end do

!$omp end do

        else

!$omp do schedule(runtime) private(j)

          do j=1,(nj-3)*njgrp*njsub

            mlat(j)=90.e0                                               &
     &        -2.e0*atan(exp((.5e0*real(1-2*j)*dy-y0)*cpj(3)))*r2d

            mlat(j)=max(min(mlat(j),90.e0),-90.e0)

          end do

!$omp end do

        end if

      end if

!$omp end parallel

! -----

      end subroutine s_paractl

!-----7--------------------------------------------------------------7--

      end module m_paractl
