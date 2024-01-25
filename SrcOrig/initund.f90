!***********************************************************************
      module m_initund
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2001/10/18
!     Modification: 2001/11/20, 2002/02/05, 2002/04/02, 2002/07/03,
!                   2002/08/27, 2002/12/02, 2003/04/30, 2003/05/19,
!                   2003/07/15, 2003/08/08, 2003/11/05, 2003/12/12,
!                   2004/02/01, 2004/03/05, 2004/04/01, 2004/04/15,
!                   2004/07/01, 2004/08/01, 2004/09/01, 2004/09/10,
!                   2005/04/04, 2006/04/03, 2006/11/06, 2007/01/31,
!                   2007/06/27, 2007/09/14, 2007/10/19, 2008/05/02,
!                   2008/07/01, 2008/08/25, 2008/10/10, 2009/01/30,
!                   2009/02/27, 2009/11/13, 2011/11/10, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     initialize the soil and sea temperature.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_comphy
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar
      use m_rdsstini
      use m_rdtund

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: initund, s_initund

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface initund

        module procedure s_initund

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log
      intrinsic min
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_initund(fpsfcdat,fpsfcopt,fpadvopt,fpdzgrd,fptgdeep, &
     &                     fpsstcst,ni,nj,nk,nund,pbr,ptbr,pp,ptp,      &
     &                     land,tund,tundp,sst,ek)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsfcdat
                       ! Formal parameter of unique index of sfcdat

      integer, intent(in) :: fpsfcopt
                       ! Formal parameter of unique index of sfcopt

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fpdzgrd
                       ! Formal parameter of unique index of dzgrd

      integer, intent(in) :: fptgdeep
                       ! Formal parameter of unique index of tgdeep

      integer, intent(in) :: fpsstcst
                       ! Formal parameter of unique index of sstcst

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

      real, intent(in) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: pp(0:ni+1,0:nj+1,1:nk)
                       ! Pressure perturbation

      real, intent(in) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbarion

! Output variables

      real, intent(out) :: tund(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at present

      real, intent(out) :: tundp(0:ni+1,0:nj+1,1:nund)
                       ! Soil and sea temperature at past

! Internal shared variables

      character(len=108) sfcdat
                       ! Control flag of input surface data type

      integer sfcopt   ! Option for surface physics
      integer advopt   ! Option for advection scheme

      real rddvcp      ! rd / cp

      real p0iv        ! 1.0 / p0

      real dzgrd       ! Grid distance in soil layers in z direction

      real tgdeep      ! Constant soil temperature in deepest layer
      real sstcst      ! Constant sea surface temperature

      real enk         ! Temporary variable
      real enkm1v      ! Temporary variable

      real, intent(inout) :: sst(0:ni+1,0:nj+1)
                       ! Sea surface temperature

      real, intent(inout) :: ek(1:nk)
                       ! Temporary variable

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(sfcdat)

! -----

! Get the required namelist variables.

      call getcname(fpsfcdat,sfcdat)
      call getiname(fpsfcopt,sfcopt)
      call getiname(fpadvopt,advopt)
      call getrname(fpdzgrd,dzgrd)
      call getrname(fptgdeep,tgdeep)
      call getrname(fpsstcst,sstcst)

! -----

! Set the common used variables.

      rddvcp=rd/cp
      p0iv=1.e0/p0

      enk=exp(real(1-nund)*dzgrd)
      enkm1v=1.e0/(exp(real(1-nund)*dzgrd)-1.e0)

! -----

! Read the sea surface temperature from external data file.

      if(sfcdat(2:2).eq.'o') then

        call rdsstini(idexprim,idcrsdir,idncexp,idnccrs,                &
     &                idwlngth,idstime,ni,nj,sst)

      end if

! -----

! Read the soil and sea temperature from restart file.

      if(sfcopt.gt.10) then

        call rdtund(idcrsdir,idprvres,idnccrs,idncprv,idadvopt,         &
     &              ni,nj,nund,tund,tundp)

      end if

! -----

!!! Initialize the soil and sea temperature.

!$omp parallel default(shared) private(k)

!! Initialized by diagnostic value.

      if(sfcopt.eq.1.or.sfcopt.eq.2.or.sfcopt.eq.3) then

! Set the surface temperature.

        if(sfcdat(2:2).eq.'o') then

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

            if(land(i,j).lt.3) then

              tundp(i,j,1)=sst(i,j)

            else

              tundp(i,j,1)=(ptbr(i,j,1)+ptp(i,j,1))                     &
     &          *exp(rddvcp*log(p0iv*(pbr(i,j,1)+pp(i,j,1))))

              tundp(i,j,1)=.5e0*(tundp(i,j,1)+(ptbr(i,j,2)+ptp(i,j,2))  &
     &          *exp(rddvcp*log(p0iv*(pbr(i,j,2)+pp(i,j,2)))))

              if(land(i,j).lt.10) then

                tundp(i,j,1)=min(tundp(i,j,1),t0)

              end if

            end if

          end do
          end do

!$omp end do

        else

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

            if(land(i,j).lt.3) then

              tundp(i,j,1)=sstcst

            else

              tundp(i,j,1)=(ptbr(i,j,1)+ptp(i,j,1))                     &
     &          *exp(rddvcp*log(p0iv*(pbr(i,j,1)+pp(i,j,1))))

              tundp(i,j,1)=.5e0*(tundp(i,j,1)+(ptbr(i,j,2)+ptp(i,j,2))  &
     &          *exp(rddvcp*log(p0iv*(pbr(i,j,2)+pp(i,j,2)))))

              if(land(i,j).lt.10) then

                tundp(i,j,1)=min(tundp(i,j,1),t0)

              end if

            end if

          end do
          end do

!$omp end do

        end if

! -----

! Set the soil temperature.

!$omp do schedule(runtime)

        do k=2,nund
          ek(k)=exp(real(1-k)*dzgrd)
        end do

!$omp end do

        do k=2,nund

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

            if(land(i,j).lt.10) then

              tundp(i,j,k)=tundp(i,j,1)

            else

              tundp(i,j,k)=((tgdeep-tundp(i,j,1))*enkm1v)*ek(k)         &
     &          +(tundp(i,j,1)*enk-tgdeep)*enkm1v

            end if

          end do
          end do

!$omp end do

        end do

! -----

! Copy the past value to the present.

        if(advopt.le.3) then

          do k=1,nund

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1
              tund(i,j,k)=tundp(i,j,k)
            end do
            end do

!$omp end do

          end do

        end if

! -----

!! -----

! Reset the sea temperature.

      else if(sfcopt.gt.10) then

        if(advopt.le.3) then

          if(sfcdat(2:2).eq.'o') then

            do k=1,nund

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1

                if(land(i,j).lt.3) then

                  tund(i,j,k)=sst(i,j)
                  tundp(i,j,k)=sst(i,j)

                end if

              end do
              end do

!$omp end do

            end do

          else

            do k=1,nund

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1

                if(land(i,j).lt.3) then

                  tund(i,j,k)=sstcst
                  tundp(i,j,k)=sstcst

                end if

              end do
              end do

!$omp end do

            end do

          end if

        else

          if(sfcdat(2:2).eq.'o') then

            do k=1,nund

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1

                if(land(i,j).lt.3) then

                  tundp(i,j,k)=sst(i,j)

                end if

              end do
              end do

!$omp end do

            end do

          else

            do k=1,nund

!$omp do schedule(runtime) private(i,j)

              do j=1,nj-1
              do i=1,ni-1

                if(land(i,j).lt.3) then

                  tundp(i,j,k)=sstcst

                end if

              end do
              end do

!$omp end do

            end do

          end if

        end if

      end if

! -----

!$omp end parallel

!!! -----

      end subroutine s_initund

!-----7--------------------------------------------------------------7--

      end module m_initund
