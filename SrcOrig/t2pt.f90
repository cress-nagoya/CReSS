!***********************************************************************
      module m_t2pt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/20
!     Modification: 1999/07/05, 1999/09/30, 1999/10/12, 2000/01/17,
!                   2000/03/23, 2000/07/05, 2001/01/15, 2001/03/13,
!                   2002/04/02, 2003/04/30, 2003/05/19, 2004/04/15,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     convert the temperature to the potential temperature.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: t2pt, s_t2pt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface t2pt

        module procedure s_t2pt

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_t2pt(nid,njd,nkd,pdat,ptdat)
!***********************************************************************

! Input variables

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

      real, intent(in) :: pdat(1:nid,1:njd,1:nkd)
                       ! Pressure in data

! Input and output variable

      real, intent(inout) :: ptdat(1:nid,1:njd,1:nkd)
                       ! Potential temperature in data

! Internal shared variable

      real rddvcp      ! rd / cp

! Internal private variables

      integer id       ! Array index in x direction
      integer jd       ! Array index in y direction
      integer kd       ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Calculate rd / cp.

      rddvcp=rd/cp

! -----

! Convert the temperature to the potential temperature.

!$omp parallel default(shared) private(kd)

      do kd=1,nkd

!$omp do schedule(runtime) private(id,jd)

        do jd=1,njd
        do id=1,nid
          ptdat(id,jd,kd)                                               &
     &      =ptdat(id,jd,kd)*exp(rddvcp*log(p0/pdat(id,jd,kd)))
        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_t2pt

!-----7--------------------------------------------------------------7--

      end module m_t2pt
