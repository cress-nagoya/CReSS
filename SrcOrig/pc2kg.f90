!***********************************************************************
      module m_pc2kg
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/12/18
!     Modification: 2001/01/15, 2001/03/13, 2001/04/15, 2001/05/29,
!                   2001/06/29, 2001/08/07, 2001/10/18, 2001/12/11,
!                   2002/01/07, 2002/04/02, 2002/04/09, 2002/11/20,
!                   2003/03/28, 2003/04/30, 2003/05/19, 2004/03/22,
!                   2004/04/15, 2004/09/10, 2005/01/31, 2005/02/10,
!                   2005/08/05, 2007/10/19, 2008/05/02, 2008/07/01,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     convert the relative humidity to the water vapor mixing ratio.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: pc2kg, s_pc2kg

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface pc2kg

        module procedure s_pc2kg

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
      subroutine s_pc2kg(nid,njd,nkd,pdat,ptdat,qvdat)
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

      real, intent(in) :: ptdat(1:nid,1:njd,1:nkd)
                       ! Potential temperature in data

! Input and output variable

      real, intent(inout) :: qvdat(1:nid,1:njd,1:nkd)
                       ! Water vapor mixing ratio in data

! Internal shared variables

      real es001       ! 0.01 x es0

      real rddvcp      ! rd / cp

      real p0iv        ! 1.0 / p0

! Internal private variables

      integer id       ! Array index in x direction
      integer jd       ! Array index in y direction
      integer kd       ! Array index in z direction

      real t           ! Temperature

      real ea          ! Pertial vapor pressure

!-----7--------------------------------------------------------------7--

! Set the common used variables.

      es001=.01e0*es0

      rddvcp=rd/cp

      p0iv=1.e0/p0

! -----

! Convert the relative humidity to the water vapor mixing ratio.

!$omp parallel default(shared) private(kd)

      do kd=1,nkd

!$omp do schedule(runtime) private(id,jd,t,ea)

        do jd=1,njd
        do id=1,nid
          t=ptdat(id,jd,kd)*exp(rddvcp*log(p0iv*pdat(id,jd,kd)))

          if(t.gt.tlow) then

            ea=es001*qvdat(id,jd,kd)*exp(17.269e0*(t-t0)/(t-35.86e0))

            qvdat(id,jd,kd)=epsva*ea/(pdat(id,jd,kd)-ea)

          else

            ea=es001*qvdat(id,jd,kd)*exp(21.875e0*(t-t0)/(t-7.66e0))

            qvdat(id,jd,kd)=epsva*ea/(pdat(id,jd,kd)-ea)

          end if

        end do
        end do

!$omp end do

      end do

!$omp end parallel

! -----

      end subroutine s_pc2kg

!-----7--------------------------------------------------------------7--

      end module m_pc2kg
