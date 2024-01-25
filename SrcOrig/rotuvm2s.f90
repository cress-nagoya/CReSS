!***********************************************************************
      module m_rotuvm2s
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/03/17
!     Modification: 1999/03/25, 1999/05/20, 1999/06/28, 1999/07/05,
!                   1999/09/30, 1999/10/12, 2000/01/05, 2000/01/17,
!                   2000/07/05, 2001/01/09, 2001/03/13, 2001/04/15,
!                   2001/05/29, 2001/06/29, 2002/04/02, 2002/09/09,
!                   2003/04/30, 2003/05/19, 2003/12/12, 2004/03/05,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     rotate the x and the y components of velocity from the projected
!     grid to the latitude and longitude grid.

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

      public :: rotuvm2s, s_rotuvm2s

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rotuvm2s

        module procedure s_rotuvm2s

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic cos
      intrinsic sin
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rotuvm2s(fpmpopt,fpnspol,fptlon,                     &
     &                      istr,iend,jstr,jend,kstr,kend,              &
     &                      cpj,imin,imax,jmin,jmax,kmin,kmax,lon,u,v)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpnspol
                       ! Formal parameter of unique index of nspol

      integer, intent(in) :: fptlon
                       ! Formal parameter of unique index of tlon

      integer, intent(in) :: istr
                       ! Minimum do loops index in x direction

      integer, intent(in) :: iend
                       ! Maximum do loops index in x direction

      integer, intent(in) :: jstr
                       ! Minimum do loops index in y direction

      integer, intent(in) :: jend
                       ! Maximum do loops index in y direction

      integer, intent(in) :: kstr
                       ! Minimum do loops index in z direction

      integer, intent(in) :: kend
                       ! Maximum do loops index in z direction

      integer, intent(in) :: imin
                       ! Minimum array index in x direction

      integer, intent(in) :: imax
                       ! Maximum array index in x direction

      integer, intent(in) :: jmin
                       ! Minimum array index in y direction

      integer, intent(in) :: jmax
                       ! Maximum array index in y direction

      integer, intent(in) :: kmin
                       ! Minimum array index in z direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      real, intent(in) :: cpj(1:7)
                       ! Map projection parameters

      real, intent(in) :: lon(imin:imax,jmin:jmax)
                       ! Longitude

! Input and output variables

      real, intent(inout) :: u(imin:imax,jmin:jmax,kmin:kmax)
                       ! x components of velocity

      real, intent(inout) :: v(imin:imax,jmin:jmax,kmin:kmax)
                       ! y components of velocity

! Internal shared variables

      integer mpopt    ! Option for map projection
      integer nspol    ! Option for projected region

      real tlon        ! True longitude

      real rpol        ! real(nspol)

      real d2rcpj      ! cpj(4) x d2r

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real dlon        ! Distance of longitude

      real cosdl       ! cos(dlon)
      real sindl       ! rpol x sin(dlon)

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

!! Rotate the x and the y components of velocity from the model grid to
!! the data sphere.

!$omp parallel default(shared) private(k)

! Rotate the x and the y components of velocity with the Polar
! Stereographic projection method.

      if(mpopt.eq.1) then

        do k=kstr,kend

!$omp do schedule(runtime) private(i,j,dlon,cosdl,sindl,utmp,vtmp)

          do j=jstr,jend
          do i=istr,iend

            if(u(i,j,k).gt.lim34n.and.v(i,j,k).gt.lim34n) then

              dlon=(lon(i,j)-tlon)*d2r

              cosdl=cos(dlon)
              sindl=rpol*sin(dlon)

              utmp=u(i,j,k)
              vtmp=v(i,j,k)

              u(i,j,k)=cosdl*utmp+sindl*vtmp
              v(i,j,k)=cosdl*vtmp-sindl*utmp

            end if

          end do
          end do

!$omp end do

        end do

! -----

! Rotate the x and the y components of velocity with the Lambert
! Conformal Conic projection method.

      else if(mpopt.eq.2) then

        do k=kstr,kend

!$omp do schedule(runtime) private(i,j,dlon,cosdl,sindl,utmp,vtmp)

          do j=jstr,jend
          do i=istr,iend

            if(u(i,j,k).gt.lim34n.and.v(i,j,k).gt.lim34n) then

              dlon=d2rcpj*(lon(i,j)-tlon)

              cosdl=cos(dlon)
              sindl=rpol*sin(dlon)

              utmp=u(i,j,k)
              vtmp=v(i,j,k)

              u(i,j,k)=cosdl*utmp+sindl*vtmp
              v(i,j,k)=cosdl*vtmp-sindl*utmp

            end if

          end do
          end do

!$omp end do

        end do

! -----

! Rotate the x and the y components of velocity without any projection
! method.

      else if(mpopt.eq.4) then

        do k=kstr,kend

!$omp do schedule(runtime) private(i,j,dlon,cosdl,sindl)

          do j=jstr,jend
          do i=istr,iend

            if(u(i,j,k).gt.lim34n.and.v(i,j,k).gt.lim34n) then

              dlon=(lon(i,j)-tlon)*d2r

              cosdl=cos(dlon)
              sindl=rpol*sin(dlon)

              u(i,j,k)=u(i,j,k)+sindl*v(i,j,k)
              v(i,j,k)=cosdl*v(i,j,k)

            end if

          end do
          end do

!$omp end do

        end do

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_rotuvm2s

!-----7--------------------------------------------------------------7--

      end module m_rotuvm2s
