!***********************************************************************
      module m_ininame
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2003/11/05
!     Modification: 2006/12/04, 2007/01/05, 2007/01/20, 2007/01/31,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2008/10/10,
!                   2009/01/30, 2009/02/27, 2013/01/28

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     initialize the table to archive namelist variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comerr
      use m_comname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: ininame, s_ininame

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface ininame

        module procedure s_ininame

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
      subroutine s_ininame
!***********************************************************************

! Internal shared variables

      integer in       ! Namelist table index

      integer ierr     ! Index of table to check namelist error

! Internal private variable

      integer in_sub   ! Substitute for in

!-----7--------------------------------------------------------------7--

!!! Initialize the table to archive namelist variables.

! For the module file, m_comerr.

      do ierr=1,nerr

        write(errlst(ierr)(1:14),'(a14)') '              '

      end do

! -----

!! For the module file, m_comname.

! For the extra defined variables.

      iname(-1)=0
      iname(0)=0

! -----

! For the character variables.

      do in=1,ncn

        call inichar(cname(in))
        call inichar(rcname(in))

      end do

! -----

! For the integer and real variables.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(in_sub)

      do in_sub=1,nin
        iname(in_sub)=0
        riname(in_sub)=0
      end do

!$omp end do

!$omp do schedule(runtime) private(in_sub)

      do in_sub=1,nrn
        rname(in_sub)=0.e0
        rrname(in_sub)=0.e0
      end do

!$omp end do

!$omp end parallel

! -----

!! -----

!!! -----

      end subroutine s_ininame

!-----7--------------------------------------------------------------7--

      end module m_ininame
