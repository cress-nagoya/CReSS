!***********************************************************************
      module m_inisfc
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/18
!     Modification: 2001/11/20, 2002/02/05, 2002/04/02, 2002/07/03,
!                   2002/08/27, 2002/12/02, 2003/04/30, 2003/05/19,
!                   2003/07/15, 2003/08/08, 2003/11/05, 2003/12/12,
!                   2004/02/01, 2004/03/05, 2004/04/01, 2004/04/15,
!                   2004/07/01, 2004/08/01, 2004/09/01, 2005/01/14,
!                   2005/02/10, 2005/08/05, 2005/12/13, 2006/09/21,
!                   2006/12/04, 2007/01/05, 2007/01/20, 2007/01/31,
!                   2007/05/14, 2007/10/19, 2008/01/11, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2008/10/10, 2009/02/27,
!                   2011/11/10, 2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     initialize the surface physical parameters.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc2d
      use m_bcycle
      use m_bcyclex
      use m_combuf
      use m_comindx
      use m_commpi
      use m_comphy
      use m_convland
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar
      use m_rdsfc
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_shiftgx
      use m_shiftgy
      use m_shiftsx
      use m_shiftsy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: inisfc, s_inisfc

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface inisfc

        module procedure s_inisfc

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic int
      intrinsic mod
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_inisfc(fpsfcdat,fpwbc,fpebc,fpexbopt,                &
     &                    fplnduse,fpdstopt,fpzsfc,fpgralbe,fpgrbeta,   &
     &                    fpgrz0m,fpgrz0h,fpgrcap,fpgrnuu,ni,nj,nk,zph, &
     &                    land,albe,beta,z0m,z0h,cap,nuu,kai,rland)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsfcdat
                       ! Formal parameter of unique index of sfcdat

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpexbopt
                       ! Formal parameter of unique index of exbopt

      integer, intent(in) :: fplnduse
                       ! Formal parameter of unique index of lnduse

      integer, intent(in) :: fpdstopt
                       ! Formal parameter of unique index of dstopt

      integer, intent(in) :: fpzsfc
                       ! Formal parameter of unique index of zsfc

      integer, intent(in) :: fpgralbe
                       ! Formal parameter of unique index of gralbe

      integer, intent(in) :: fpgrbeta
                       ! Formal parameter of unique index of grbeta

      integer, intent(in) :: fpgrz0m
                       ! Formal parameter of unique index of grz0m

      integer, intent(in) :: fpgrz0h
                       ! Formal parameter of unique index of grz0h

      integer, intent(in) :: fpgrcap
                       ! Formal parameter of unique index of grcap

      integer, intent(in) :: fpgrnuu
                       ! Formal parameter of unique index of grnuu

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

! Output variables

      integer, intent(out) :: land(0:ni+1,0:nj+1)
                       ! Land use of surface

      real, intent(out) :: albe(0:ni+1,0:nj+1)
                       ! Albedo

      real, intent(out) :: beta(0:ni+1,0:nj+1)
                       ! Evapotranspiration efficiency

      real, intent(out) :: z0m(0:ni+1,0:nj+1)
                       ! Roughness length for velocity

      real, intent(out) :: z0h(0:ni+1,0:nj+1)
                       ! Roughness length for scalar

      real, intent(out) :: cap(0:ni+1,0:nj+1)
                       ! Thermal capacity

      real, intent(out) :: nuu(0:ni+1,0:nj+1)
                       ! Thermal diffusivity

      real, intent(out) :: kai(0:ni+1,0:nj+1)
                       ! Sea ice distribution

! Internal shared variables

      character(len=108) sfcdat
                       ! Control flag of input surface data type

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions

      integer exbopt   ! Option for external boundary forcing

      integer lnduse   ! User specified land use category

      integer dstopt   ! Option for sea ice distribution

      integer nx       ! Composite model dimension in x direction
      integer ny       ! Composite model dimension in y direction

      integer icmin    ! Minimum array index
                       ! in x direction in composite model dimension

      integer icmax    ! Maximum array index
                       ! in x direction in composite model dimension

      integer jcmin    ! Minimum array index
                       ! in y direction in composite model dimension

      integer jcmax    ! Maximum array index
                       ! in y direction in composite model dimension

      integer ic       ! Array index
                       ! in x direction in composite model dimension

      integer jc       ! Array index
                       ! in y direction in composite model dimension

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer ix       ! Parameter to make random noize

      real zsfc        ! Sea surface terrain height

      real gralbe      ! Albedo on soil surface

      real grbeta      ! Evapotranspiration efficiency on soil surface

      real grz0m       ! Roughness length for velocity on soil surface
      real grz0h       ! Roughness length for scalar on soil surface

      real grcap       ! Thermal capacity of soil
      real grnuu       ! Thermal diffusivity of soil

      real, intent(inout) :: rland(0:ni+1,0:nj+1)
                       ! Real land use of surface

! Internal private variables

      integer i_sub    ! Substitute for i
      integer j_sub    ! Substitute for j

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(sfcdat)

! -----

! Get the required namelist variables.

      call getcname(fpsfcdat,sfcdat)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpexbopt,exbopt)
      call getiname(fplnduse,lnduse)
      call getiname(fpdstopt,dstopt)
      call getrname(fpzsfc,zsfc)
      call getrname(fpgralbe,gralbe)
      call getrname(fpgrbeta,grbeta)
      call getrname(fpgrz0m,grz0m)
      call getrname(fpgrz0h,grz0h)
      call getrname(fpgrcap,grcap)
      call getrname(fpgrnuu,grnuu)

! -----

! Set the common used variables.

      nx=(ni-3)*nigrp*nisub+3
      ny=(nj-3)*njgrp*njsub+3

      icmin=(ni-3)*nisub*igrp+(ni-3)*isub
      icmax=(ni-3)*nisub*igrp+(ni-3)*isub+ni
      jcmin=(nj-3)*njsub*jgrp+(nj-3)*jsub
      jcmax=(nj-3)*njsub*jgrp+(nj-3)*jsub+nj

      ix=0

! -----

! Read out the data from the interpolated surface file.

      if(sfcdat(1:1).eq.'o'.or.sfcdat(3:3).eq.'o') then

        call rdsfc(idexprim,idcrsdir,idsfcdat,idncexp,idnccrs,idwlngth, &
     &             ni,nj,land,albe,beta,z0m,z0h,cap,nuu,kai)

      end if

! -----

!!! Reset the sea surface temperature, the sea ice distribution, the
!!! land use categories, evapotranspiration efficiency, albedo,
!!! roughness length, thermal capacity and thermal diffusivity.

! Reset the sea ice distribution by the all-or-nothing arrangement.

      if(sfcdat(3:3).eq.'o') then

        if(dstopt.eq.2) then

          do jc=1,ny-1
          do ic=1,nx-1

            ix=ix+1
            ix=mod(97*mod(ix,65536),65536)

            if((ic.gt.icmin.and.ic.lt.icmax)                            &
     &        .and.(jc.gt.jcmin.and.jc.lt.jcmax)) then

              i=ic-icmin
              j=jc-jcmin

              if(655.36e0*kai(i,j).gt.real(ix)) then
                kai(i,j)=100.1e0
              else
                kai(i,j)=.1e0
              end if

            end if

          end do
          end do

        end if

      end if

! -----

! Reset the user specified land use category in case of no land use
! data.

      if(sfcdat(1:1).eq.'x') then

        if(lnduse.eq.1) then
          lnduse=10
        else if(lnduse.eq.2) then
          lnduse=5
        end if

      end if

! -----

!! Reset the sea surface temperature, the land use categories,
!! evapotranspiration efficiency, albedo, roughness length,
!! thermal capacity and thermal diffusivity.

!$omp parallel default(shared)

! Reset the land use categories.

      if(sfcdat(1:1).eq.'o'.and.sfcdat(3:3).eq.'o') then

!$omp do schedule(runtime) private(i_sub,j_sub)

        do j_sub=1,nj-1
        do i_sub=1,ni-1

          if(land(i_sub,j_sub).lt.0) then

            if(int(kai(i_sub,j_sub)).gt.0                               &
     &        .and.int(kai(i_sub,j_sub)).lt.100) then

              land(i_sub,j_sub)=1

            else if(int(kai(i_sub,j_sub)).eq.100) then

              land(i_sub,j_sub)=3

            end if

          else if(land(i_sub,j_sub).ge.0                                &
     &      .and.land(i_sub,j_sub).lt.5) then

            land(i_sub,j_sub)=3

          end if

        end do
        end do

!$omp end do

      else if(sfcdat(1:1).eq.'o'.and.sfcdat(3:3).eq.'x') then

!$omp do schedule(runtime) private(i_sub,j_sub)

        do j_sub=1,nj-1
        do i_sub=1,ni-1

          if(land(i_sub,j_sub).ge.0.and.land(i_sub,j_sub).lt.5) then

            land(i_sub,j_sub)=3

          end if

        end do
        end do

!$omp end do

      else if(sfcdat(1:1).eq.'x'.and.sfcdat(3:3).eq.'o') then

!$omp do schedule(runtime) private(i_sub,j_sub)

        do j_sub=1,nj-1
        do i_sub=1,ni-1

          if(zph(i_sub,j_sub,2).le.zsfc) then

            if(int(kai(i_sub,j_sub)).gt.0                               &
     &        .and.int(kai(i_sub,j_sub)).lt.100) then

              land(i_sub,j_sub)=1

            else if(int(kai(i_sub,j_sub)).eq.100) then

              land(i_sub,j_sub)=3

            else

              land(i_sub,j_sub)=-1

            end if

          else

            land(i_sub,j_sub)=lnduse

          end if

        end do
        end do

!$omp end do

      else if(sfcdat(1:1).eq.'x'.and.sfcdat(3:3).eq.'x') then

!$omp do schedule(runtime) private(i_sub,j_sub)

        do j_sub=1,nj-1
        do i_sub=1,ni-1

          if(zph(i_sub,j_sub,2).le.zsfc) then

            land(i_sub,j_sub)=-1

          else

            land(i_sub,j_sub)=lnduse

          end if

        end do
        end do

!$omp end do

      end if

! -----

! Reset the sea ice distribution by the all-or-nothing arrangement.

      if(sfcdat(3:3).eq.'o') then

        if(dstopt.eq.1) then

!$omp do schedule(runtime) private(i_sub,j_sub)

          do j_sub=1,nj-1
          do i_sub=1,ni-1
            kai(i_sub,j_sub)=.01e0*kai(i_sub,j_sub)
          end do
          end do

!$omp end do

        end if

      end if

! -----

! Reset the evapotranspiration efficiency, albedo, roughness length,
! thermal capacity and thermal diffusivity.

      if(sfcdat(1:1).eq.'o') then

!$omp do schedule(runtime) private(i_sub,j_sub)

        do j_sub=1,nj-1
        do i_sub=1,ni-1

          if(land(i_sub,j_sub).lt.3) then

            albe(i_sub,j_sub)=sealbe

            beta(i_sub,j_sub)=sebeta

            z0m(i_sub,j_sub)=sez0m
            z0h(i_sub,j_sub)=sez0h

            cap(i_sub,j_sub)=secap
            nuu(i_sub,j_sub)=senuu

          else if(land(i_sub,j_sub).ge.3                                &
     &      .and.land(i_sub,j_sub).lt.5) then

            albe(i_sub,j_sub)=icalbe

            beta(i_sub,j_sub)=icbeta

            z0m(i_sub,j_sub)=icz0m
            z0h(i_sub,j_sub)=icz0h

            cap(i_sub,j_sub)=0.e0
            nuu(i_sub,j_sub)=0.e0

          else if(land(i_sub,j_sub).ge.5                                &
     &      .and.land(i_sub,j_sub).lt.10) then

            albe(i_sub,j_sub)=snalbe

            beta(i_sub,j_sub)=snbeta

            z0m(i_sub,j_sub)=snz0m
            z0h(i_sub,j_sub)=snz0h

            cap(i_sub,j_sub)=0.e0
            nuu(i_sub,j_sub)=0.e0

          end if

        end do
        end do

!$omp end do

      else

!$omp do schedule(runtime) private(i_sub,j_sub)

        do j_sub=1,nj-1
        do i_sub=1,ni-1

          if(land(i_sub,j_sub).lt.3) then

            albe(i_sub,j_sub)=sealbe

            beta(i_sub,j_sub)=sebeta

            z0m(i_sub,j_sub)=sez0m
            z0h(i_sub,j_sub)=sez0h

            cap(i_sub,j_sub)=secap
            nuu(i_sub,j_sub)=senuu

          else if(land(i_sub,j_sub).ge.3                                &
     &      .and.land(i_sub,j_sub).lt.5) then

            albe(i_sub,j_sub)=icalbe

            beta(i_sub,j_sub)=icbeta

            z0m(i_sub,j_sub)=icz0m
            z0h(i_sub,j_sub)=icz0h

            cap(i_sub,j_sub)=0.e0
            nuu(i_sub,j_sub)=0.e0

          else if(land(i_sub,j_sub).ge.5                                &
     &      .and.land(i_sub,j_sub).lt.10) then

            albe(i_sub,j_sub)=snalbe

            beta(i_sub,j_sub)=snbeta

            z0m(i_sub,j_sub)=snz0m
            z0h(i_sub,j_sub)=snz0h

            cap(i_sub,j_sub)=0.e0
            nuu(i_sub,j_sub)=0.e0

          else

            albe(i_sub,j_sub)=gralbe

            beta(i_sub,j_sub)=grbeta

            z0m(i_sub,j_sub)=grz0m
            z0h(i_sub,j_sub)=grz0h

            cap(i_sub,j_sub)=grcap
            nuu(i_sub,j_sub)=grnuu

          end if

        end do
        end do

!$omp end do

      end if

! ----

!$omp end parallel

!! -----

!!! -----

!! Set the lateral boundary conditions.

      if(exbopt.eq.0) then

! Convert the integer land use categories to real.

        call convland('real   ',ni,nj,land,rland)

! -----

! Exchange the variables horizontally.

        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,rland,1,8,sbuf)
        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,albe,2,8,sbuf)
        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,beta,3,8,sbuf)
        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,z0m,4,8,sbuf)
        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,z0h,5,8,sbuf)
        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,cap,6,8,sbuf)
        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,nuu,7,8,sbuf)
        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,kai,8,8,sbuf)

        call s_shiftsx(idwbc,idebc,'bnd',nj,1,8,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,rland,1,8,rbuf)
        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,albe,2,8,rbuf)
        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,beta,3,8,rbuf)
        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,z0m,4,8,rbuf)
        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,z0h,5,8,rbuf)
        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,cap,6,8,rbuf)
        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,nuu,7,8,rbuf)
        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,kai,8,8,rbuf)

        call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,rland,1,8,sbuf)
        call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,albe,2,8,sbuf)
        call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,beta,3,8,sbuf)
        call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,z0m,4,8,sbuf)
        call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,z0h,5,8,sbuf)
        call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,cap,6,8,sbuf)
        call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,nuu,7,8,sbuf)
        call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,kai,8,8,sbuf)

        call s_shiftsy(idsbc,idnbc,'bnd',ni,1,8,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,rland,1,8,rbuf)
        call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,albe,2,8,rbuf)
        call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,beta,3,8,rbuf)
        call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,z0m,4,8,rbuf)
        call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,z0h,5,8,rbuf)
        call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,cap,6,8,rbuf)
        call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,nuu,7,8,rbuf)
        call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,kai,8,8,rbuf)

        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,rland,1,8,sbuf)
        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,albe,2,8,sbuf)
        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,beta,3,8,sbuf)
        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,z0m,4,8,sbuf)
        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,z0h,5,8,sbuf)
        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,cap,6,8,sbuf)
        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,nuu,7,8,sbuf)
        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,kai,8,8,sbuf)

        call s_shiftgx(idwbc,idebc,'bnd',nj,1,8,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,rland,1,8,rbuf)
        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,albe,2,8,rbuf)
        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,beta,3,8,rbuf)
        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,z0m,4,8,rbuf)
        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,z0h,5,8,rbuf)
        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,cap,6,8,rbuf)
        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,nuu,7,8,rbuf)
        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,kai,8,8,rbuf)

        call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,rland,1,8,sbuf)
        call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,albe,2,8,sbuf)
        call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,beta,3,8,sbuf)
        call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,z0m,4,8,sbuf)
        call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,z0h,5,8,sbuf)
        call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,cap,6,8,sbuf)
        call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,nuu,7,8,sbuf)
        call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,1,kai,8,8,sbuf)

        call s_shiftgy(idsbc,idnbc,'bnd',ni,1,8,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,rland,1,8,rbuf)
        call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,albe,2,8,rbuf)
        call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,beta,3,8,rbuf)
        call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,z0m,4,8,rbuf)
        call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,z0h,5,8,rbuf)
        call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,cap,6,8,rbuf)
        call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,nuu,7,8,rbuf)
        call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,1,kai,8,8,rbuf)

        call s_bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,1,rland)

        call s_bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,1,albe)

        call s_bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,1,beta)

        call s_bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,1,z0m)

        call s_bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,1,z0h)

        call s_bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,1,cap)

        call s_bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,1,nuu)

        call s_bcycle(idwbc,idebc,idsbc,idnbc,                          &
     &                2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,1,kai)

        call bc2d(idwbc,idebc,idsbc,idnbc,ni,nj,rland)
        call bc2d(idwbc,idebc,idsbc,idnbc,ni,nj,albe)
        call bc2d(idwbc,idebc,idsbc,idnbc,ni,nj,beta)
        call bc2d(idwbc,idebc,idsbc,idnbc,ni,nj,z0m)
        call bc2d(idwbc,idebc,idsbc,idnbc,ni,nj,z0h)
        call bc2d(idwbc,idebc,idsbc,idnbc,ni,nj,cap)
        call bc2d(idwbc,idebc,idsbc,idnbc,ni,nj,nuu)
        call bc2d(idwbc,idebc,idsbc,idnbc,ni,nj,kai)

! -----

! Convert the real land use categories to integer.

        call convland('integer',ni,nj,land,rland)

! -----

      end if

!! -----

!! Set the periodic boundary conditions in x direction.

      if(exbopt.ge.1.and.abs(wbc).eq.1.and.abs(ebc).eq.1) then

! Convert the integer land use categories to real.

        call convland('real   ',ni,nj,land,rland)

! -----

! Exchange the variables horizontally.

        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,rland,1,8,sbuf)
        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,albe,2,8,sbuf)
        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,beta,3,8,sbuf)
        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,z0m,4,8,sbuf)
        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,z0h,5,8,sbuf)
        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,cap,6,8,sbuf)
        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,nuu,7,8,sbuf)
        call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,kai,8,8,sbuf)

        call s_shiftsx(idwbc,idebc,'bnd',nj,1,8,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,rland,1,8,rbuf)
        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,albe,2,8,rbuf)
        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,beta,3,8,rbuf)
        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,z0m,4,8,rbuf)
        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,z0h,5,8,rbuf)
        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,cap,6,8,rbuf)
        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,nuu,7,8,rbuf)
        call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,kai,8,8,rbuf)

        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,rland,1,8,sbuf)
        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,albe,2,8,sbuf)
        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,beta,3,8,sbuf)
        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,z0m,4,8,sbuf)
        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,z0h,5,8,sbuf)
        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,cap,6,8,sbuf)
        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,nuu,7,8,sbuf)
        call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,1,kai,8,8,sbuf)

        call s_shiftgx(idwbc,idebc,'bnd',nj,1,8,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,rland,1,8,rbuf)
        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,albe,2,8,rbuf)
        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,beta,3,8,rbuf)
        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,z0m,4,8,rbuf)
        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,z0h,5,8,rbuf)
        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,cap,6,8,rbuf)
        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,nuu,7,8,rbuf)
        call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,1,kai,8,8,rbuf)

        call s_bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,1,rland)
        call s_bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,1,albe)
        call s_bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,1,beta)
        call s_bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,1,z0m)
        call s_bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,1,z0h)
        call s_bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,1,cap)
        call s_bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,1,nuu)
        call s_bcyclex(idwbc,idebc,2,1,ni-2,ni-1,ni,nj,1,kai)

! -----

! Convert the real land use categories to integer.

        call convland('integer',ni,nj,land,rland)

! -----

      end if

!! -----

      end subroutine s_inisfc

!-----7--------------------------------------------------------------7--

      end module m_inisfc
