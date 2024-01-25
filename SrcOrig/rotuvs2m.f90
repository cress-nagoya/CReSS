!***********************************************************************
      module m_rotuvs2m
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/03/17
!     Modification: 1999/03/25, 1999/05/20, 1999/06/28, 1999/07/05,
!                   1999/09/30, 1999/10/12, 2000/01/05, 2000/01/17,
!                   2000/07/05, 2001/01/09, 2001/03/13, 2001/04/15,
!                   2001/05/29, 2002/04/02, 2002/09/09, 2003/04/30,
!                   2003/05/19, 2003/12/12, 2004/03/05, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     rotate the x and the y components of velocity from the latitude
!     and longitude grid to the projected grid.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_getiname
      use m_getrname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rotuvs2m, s_rotuvs2m

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rotuvs2m

        module procedure s_rotuvs2m

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic cos
      intrinsic sin
      intrinsic tan
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rotuvs2m(fpmpopt,fpnspol,fptlon,                     &
     &                      cpj,nid,njd,nkd,londat,udat,vdat)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpnspol
                       ! Formal parameter of unique index of nspol

      integer, intent(in) :: fptlon
                       ! Formal parameter of unique index of tlon

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

      real, intent(in) :: cpj(1:7)
                       ! Map projection parameters

      real, intent(in) :: londat(1:nid,1:njd)
                       ! Longitude in data

! Input and output variables

      real, intent(inout) :: udat(1:nid,1:njd,1:nkd)
                       ! x components of velocity in data

      real, intent(inout) :: vdat(1:nid,1:njd,1:nkd)
                       ! y components of velocity in data

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer nspol    ! Option for projected region

      real tlon        ! True longitude

      real rpol        ! real(nspol)

      real d2rcpj      ! cpj(4) x d2r

! Internal private variables

      integer id       ! Array index in x direction
      integer jd       ! Array index in y direction
      integer kd       ! Array index in z direction

      real dlon        ! Distance of longitude

      real cosdl       ! cos(dlon)
      real sindl       ! rpol x sin(dlon)
      real tandl       ! tan(dlon)
      real cosiv       ! 1.0 / cos(dlon)

      real utmp        ! Temporary variable of x components of velocity
      real vtmp        ! Temporary variable of y components of velocity

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpmpopt,mpopt)
      call getiname(fpnspol,nspol)
      call getrname(fptlon,tlon)

! -----

! Set the common used variables.

      rpol=real(nspol)

      d2rcpj=cpj(4)*d2r

! -----

!! Rotate the x and the y components of velocity from the data sphere to
!! the model grid.

!$omp parallel default(shared) private(kd)

! Rotate the x and the y components of velocity with the Polar
! Stereographic projection method.

      if(mpopt.eq.1) then

        do kd=1,nkd

!$omp do schedule(runtime) private(id,jd,dlon,cosdl,sindl,utmp,vtmp)

          do jd=1,njd
          do id=1,nid

            if(udat(id,jd,kd).gt.lim34n                                 &
     &        .and.vdat(id,jd,kd).gt.lim34n) then

              dlon=(londat(id,jd)-tlon)*d2r

              cosdl=cos(dlon)
              sindl=rpol*sin(dlon)

              utmp=udat(id,jd,kd)
              vtmp=vdat(id,jd,kd)

              udat(id,jd,kd)=cosdl*utmp-sindl*vtmp
              vdat(id,jd,kd)=sindl*utmp+cosdl*vtmp

            end if

          end do
          end do

!$omp end do

        end do

! -----

! Rotate the x and the y components of velocity with the Lambert
! Conformal Conic projection method.

      else if(mpopt.eq.2) then

        do kd=1,nkd

!$omp do schedule(runtime) private(id,jd,dlon,cosdl,sindl,utmp,vtmp)

          do jd=1,njd
          do id=1,nid

            if(udat(id,jd,kd).gt.lim34n                                 &
     &        .and.vdat(id,jd,kd).gt.lim34n) then

              dlon=d2rcpj*(londat(id,jd)-tlon)

              cosdl=cos(dlon)
              sindl=rpol*sin(dlon)

              utmp=udat(id,jd,kd)
              vtmp=vdat(id,jd,kd)

              udat(id,jd,kd)=cosdl*utmp-sindl*vtmp
              vdat(id,jd,kd)=sindl*utmp+cosdl*vtmp

            end if

          end do
          end do

!$omp end do

        end do

! -----

! Rotate the x and the y components of velocity without any projection
! method.

      else if(mpopt.eq.4) then

        do kd=1,nkd

!$omp do schedule(runtime) private(id,jd,dlon,tandl,cosiv)

          do jd=1,njd
          do id=1,nid

            if(udat(id,jd,kd).gt.lim34n                                 &
     &        .and.vdat(id,jd,kd).gt.lim34n) then

              dlon=(londat(id,jd)-tlon)*d2r

              tandl=rpol*(dlon)
              cosiv=1.e0/cos(dlon)

              udat(id,jd,kd)=udat(id,jd,kd)-tandl*vdat(id,jd,kd)
              vdat(id,jd,kd)=cosiv*vdat(id,jd,kd)

            end if

          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_rotuvs2m

!-----7--------------------------------------------------------------7--

      end module m_rotuvs2m
