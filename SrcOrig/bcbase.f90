!***********************************************************************
      module m_bcbase
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/04/06, 1999/07/05, 1999/08/03, 1999/08/18,
!                   1999/08/23, 1999/09/30, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/03/08, 2000/07/05, 2001/04/15,
!                   2001/05/29, 2001/12/11, 2002/04/02, 2002/06/18,
!                   2002/08/15, 2003/04/30, 2003/05/19, 2003/11/05,
!                   2003/12/12, 2004/01/09, 2004/03/05, 2004/04/15,
!                   2004/08/20, 2005/02/10, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the bottom and the top boundary conditions for the base state
!     variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bcyclex
      use m_bcycley
      use m_combuf
      use m_comindx
      use m_comphy
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getiname
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

      public :: bcbase, s_bcbase

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface bcbase

        module procedure s_bcbase

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log
      intrinsic mod

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_bcbase(fpsmtopt,ni,nj,nk,zph8s,                      &
     &                    ubr,vbr,pbr,ptbr,qvbr,rbr,pibr,ptvbr)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsmtopt
                       ! Formal parameter of unique index of smtopt

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zph8s(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates at scalar points

! Input and output variables

      real, intent(inout) :: ubr(0:ni+1,0:nj+1,1:nk)
                       ! Base state x components of velocity

      real, intent(inout) :: vbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state y components of velocity

      real, intent(inout) :: pbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state pressure

      real, intent(inout) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(inout) :: qvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state water vapor mixing ratio

      real, intent(inout) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(inout) :: pibr(0:ni+1,0:nj+1,1:nk)
                       ! Base state Exner function

      real, intent(inout) :: ptvbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state virtual potential temperature

! Internal shared variables

      integer smtopt   ! Option for numerical smoothing

      integer nkm1     ! nk - 1
      integer nkm2     ! nk - 2

      real gdvcp2      ! 2.0 x g / cp

      real cpdvrd      ! cp / rd

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpsmtopt,smtopt)

! -----

! Set the common used variables

      nkm1=nk-1
      nkm2=nk-2

      gdvcp2=2.e0*g/cp

      cpdvrd=cp/rd

! -----

!! Set the bottom and the top boundary conditions for the base state
!! variables.

!$omp parallel default(shared)

! Set the bottom and the top boundary conditions for the base state
! velocity.

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=1,ni
        ubr(i,j,1)=ubr(i,j,2)
        ubr(i,j,nkm1)=ubr(i,j,nkm2)
      end do
      end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

      do j=1,nj
      do i=0,ni
        vbr(i,j,1)=vbr(i,j,2)
        vbr(i,j,nkm1)=vbr(i,j,nkm2)
      end do
      end do

!$omp end do

! -----

! Set the bottom and the top boundary conditions for the base state
! potential temperature and water vapor mixing ratio.

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=0,ni
        ptbr(i,j,1)=ptbr(i,j,2)
        ptbr(i,j,nkm1)=ptbr(i,j,nkm2)

        qvbr(i,j,1)=qvbr(i,j,2)
        qvbr(i,j,nkm1)=qvbr(i,j,nkm2)

      end do
      end do

!$omp end do

! -----

! Set the bottom and the top boundary conditions for the base state
! virtual potential temperature.

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=0,ni
        ptvbr(i,j,1)=ptvbr(i,j,2)
        ptvbr(i,j,nkm1)=ptvbr(i,j,nkm2)
      end do
      end do

!$omp end do

! -----

! Set the bottom and the top boundary conditions for the base state
! exnar function.

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=0,ni
        pibr(i,j,1)=pibr(i,j,2)                                         &
     &   +gdvcp2*(zph8s(i,j,2)-zph8s(i,j,1))/(ptvbr(i,j,1)+ptvbr(i,j,2))

        pibr(i,j,nkm1)=pibr(i,j,nkm2)-gdvcp2*(zph8s(i,j,nkm1)           &
     &   -zph8s(i,j,nkm2))/(ptvbr(i,j,nkm2)+ptvbr(i,j,nkm1))

      end do
      end do

!$omp end do

! -----

! Set the bottom and the top boundary conditions for the base state
! pressure and the base state density.

!$omp do schedule(runtime) private(i,j)

      do j=0,nj
      do i=0,ni
        pbr(i,j,1)=p0*exp(cpdvrd*log(pibr(i,j,1)))
        pbr(i,j,nkm1)=p0*exp(cpdvrd*log(pibr(i,j,nkm1)))

        rbr(i,j,1)=pbr(i,j,1)/(rd*ptvbr(i,j,1)*pibr(i,j,1))
        rbr(i,j,nkm1)=pbr(i,j,nkm1)/(rd*ptvbr(i,j,nkm1)*pibr(i,j,nkm1))

      end do
      end do

!$omp end do

! -----

!$omp end parallel

!! -----

!! Exchange the value in the case the 4th order calculation is
!! performed.

      if(mod(smtopt,10).eq.2.or.mod(smtopt,10).eq.3) then

! Exchange the value in x direction.

        call s_putbufsx(idwbc,idebc,'all',4,ni-3,ni,nj,nk,ubr,1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'all',0,ni+1,ni,nj,nk,ubr,1,1,rbuf)

        call s_putbufgx(idwbc,idebc,'all',4,ni-3,ni,nj,nk,ubr,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',0,ni+1,ni,nj,nk,ubr,1,1,rbuf)

        call bcyclex(idwbc,idebc,4,0,ni-3,ni+1,ni,nj,nk,ubr)

! -----

! Exchange the value in y direction.

        call s_putbufsy(idsbc,idnbc,'all',4,nj-3,ni,nj,nk,vbr,1,1,sbuf)

        call s_shiftsy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'all',0,nj+1,ni,nj,nk,vbr,1,1,rbuf)

        call s_putbufgy(idsbc,idnbc,'all',4,nj-3,ni,nj,nk,vbr,1,1,sbuf)

        call s_shiftgy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'all',0,nj+1,ni,nj,nk,vbr,1,1,rbuf)

        call bcycley(idsbc,idnbc,4,0,nj-3,nj+1,ni,nj,nk,vbr)

! -----

      end if

!! -----

      end subroutine s_bcbase

!-----7--------------------------------------------------------------7--

      end module m_bcbase
