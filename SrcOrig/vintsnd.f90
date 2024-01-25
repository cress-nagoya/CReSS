!***********************************************************************
      module m_vintsnd
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/03/25, 1999/04/06, 1999/06/14,
!                   1999/06/21, 1999/11/01, 1999/11/19, 2000/01/17,
!                   2000/04/18, 2001/01/15, 2001/05/29, 2001/10/17,
!                   2002/04/02, 2002/06/18, 2002/07/15, 2003/03/28,
!                   2003/04/30, 2003/05/19, 2003/06/27, 2003/12/12,
!                   2006/01/10, 2006/09/21, 2007/01/31, 2007/10/19,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     interpolate the sounding data to fine interval levels vartically.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_copy1d
      use m_getcname
      use m_inichar
      use m_vint11

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vintsnd, s_vintsnd

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vintsnd

        module procedure s_vintsnd

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_vintsnd(fpsndtyp,nsnd,zsnd,ltmp1,nlev,u1d,v1d,pt1d,  &
     &                     qv1d,z1d)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsndtyp
                       ! Formal parameter of unique index of sndtyp

      integer, intent(in) :: nsnd
                       ! Sounding data dimension

      integer, intent(in) :: nlev
                       ! Horizontally averaged vertical dimension

      real, intent(in) :: zsnd(1:nsnd)
                       ! z physical coordinates in sounding data

! Input and output variables

      real, intent(inout) :: u1d(1:nlev)
                       ! Horizontally averaged x components of velocity

      real, intent(inout) :: v1d(1:nlev)
                       ! Horizontally averaged y components of velocity

      real, intent(inout) :: pt1d(1:nlev)
                       ! Horizontally averaged potential temrerature

      real, intent(inout) :: qv1d(1:nlev)
                       ! Horizontally averaged water vapor mixing ratio

! Output variable

      real, intent(out) :: z1d(1:nlev)
                       ! Horizontally averaged z physical coordinates

! Internal shared variables

      character(len=108) sndtyp
                       ! Control flag of sounding data type

      real nlm1iv      ! 1.0 / real(nlev - 1)

      real, intent(inout) :: ltmp1(1:nsnd)
                       ! Temporary array

! Internal private variable

      integer kl       ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(sndtyp)

! -----

! Get the required namelist variable.

      call getcname(fpsndtyp,sndtyp)

! -----

!! Interpolate the sounding data to fine interval levels vartically.

      if(sndtyp(1:1).eq.'z') then

! Set the common used variables.

        nlm1iv=1.e0/real(nlev-1)

        z1d(1)=zsnd(1)
        z1d(nlev)=zsnd(nsnd)

! -----

! Calculate the interpolated zeta coordinates.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(kl)

        do kl=2,nlev-1
         z1d(kl)=z1d(1)+nlm1iv*real(kl-1)*(z1d(nlev)-z1d(1))
        end do

!$omp end do

!$omp end parallel

! -----

! Perform interpolation.

        pt1d(nlev)=pt1d(nsnd)
        qv1d(nlev)=qv1d(nsnd)

        u1d(nlev)=u1d(nsnd)
        v1d(nlev)=v1d(nsnd)

        call copy1d(1,nsnd,pt1d,ltmp1)

        call vint11(nlev-1,z1d,pt1d,nsnd,zsnd,ltmp1)

        call copy1d(1,nsnd,qv1d,ltmp1)

        call vint11(nlev-1,z1d,qv1d,nsnd,zsnd,ltmp1)

        call copy1d(1,nsnd,u1d,ltmp1)

        call vint11(nlev-1,z1d,u1d,nsnd,zsnd,ltmp1)

        call copy1d(1,nsnd,v1d,ltmp1)

        call vint11(nlev-1,z1d,v1d,nsnd,zsnd,ltmp1)

! -----

      end if

!! -----

! Copy the array zsnd to the z1d.

      if(sndtyp(1:1).eq.'p') then

        call copy1d(1,nsnd,zsnd,z1d)

      end if

! -----

      end subroutine s_vintsnd

!-----7--------------------------------------------------------------7--

      end module m_vintsnd
