!***********************************************************************
      module m_warmblk
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/11/01
!     Modification: 2000/01/17, 2000/03/08, 2000/04/18, 2000/06/01,
!                   2000/07/05, 2000/12/18, 2001/04/15, 2001/06/29,
!                   2001/10/18, 2002/01/15, 2003/02/13, 2003/05/19,
!                   2003/12/12, 2004/03/22, 2004/05/31, 2004/06/10,
!                   2004/08/01, 2004/09/01, 2004/09/25, 2004/12/17,
!                   2006/01/10, 2006/02/13, 2006/04/03, 2006/05/12,
!                   2006/09/30, 2008/05/02, 2008/07/01, 2008/08/25,
!                   2009/02/27

!     Author      : Satoki Tsujino
!     Modification: 2024/12/25

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate the source amounts of the bulk warm cloud micro physics.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_collc2r
      use m_comindx
      use m_convc2r
      use m_dmpcld
      use m_evapr2v
      use m_fallqr
      use m_termqr

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: warmblk, s_warmblk

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface warmblk

        module procedure s_warmblk

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
      subroutine s_warmblk(dtb,ni,nj,nk,jcb,ptbr,rbr,rst,rbv,pi,p,      &
     &                     ptpp,qvp,qcp,qrp,ptpf,qvf,qcf,qrf,prr,       &
     &                     urq,tmp1)
!***********************************************************************

! Input variables

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dtb
                       ! Large time steps interval

      real, intent(in) :: jcb(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian

      real, intent(in) :: ptbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state potential temperature

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: rst(0:ni+1,0:nj+1,1:nk)
                       ! Base state density x Jacobian

      real, intent(in) :: rbv(0:ni+1,0:nj+1,1:nk)
                       ! Inverse of base state density

      real, intent(in) :: pi(0:ni+1,0:nj+1,1:nk)
                       ! Exnar function

      real, intent(in) :: p(0:ni+1,0:nj+1,1:nk)
                       ! Pressure

      real, intent(in) :: ptpp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at past

      real, intent(in) :: qvp(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at past

      real, intent(in) :: qcp(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio at past

      real, intent(in) :: qrp(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio at past

! Input and output variables

      real, intent(inout) :: ptpf(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation at future

      real, intent(inout) :: qvf(0:ni+1,0:nj+1,1:nk)
                       ! Water vapor mixing ratio at future

      real, intent(inout) :: qcf(0:ni+1,0:nj+1,1:nk)
                       ! Cloud water mixing ratio at future

      real, intent(inout) :: qrf(0:ni+1,0:nj+1,1:nk)
                       ! Rain water mixing ratio at future

      real, intent(inout) :: prr(0:ni+1,0:nj+1,1:2)
                       ! Precipitation and accumulation for rain

! Internal shared variables

      real, intent(inout) :: urq(0:ni+1,0:nj+1,1:nk)
                       ! Terminal velocity of rain water

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer :: i, j, k
      real :: tmpdmp(0:ni+1,0:nj+1,1:nk)

!-----7--------------------------------------------------------------7--

! Calculate the terminal velocity of the rain water.

      call termqr(idthresq,ni,nj,nk,rbr,rbv,qrp,urq)

! -----

!! Calculate the source amounts of the bulk warm cloud micro physics.

! Calculate the auto correction rate from the cloud water to the rain
! water.

      call convc2r(dtb,ni,nj,nk,qcp,qcf,qrf)

! -----

! Calculate the collection rate between the cloud water and the rain
! water.

      call collc2r(idthresq,dtb,ni,nj,nk,qcp,qrp,qcf,qrf)

! -----

! Setting dqdt to ptpp

      do k=1,nk
         do j=0,nj+1
            do i=0,ni+1
               tmpdmp(i,j,k)=ptpp(i,j,k)
            end do
         end do
      end do

      call s_dmpcld( ni, nj, nk, tmpdmp, 'dqdt ', 'o', 1.0e0 )

! -----

! Calculate the evaporation rate from the rain water to the water vapor.

      call evapr2v(idthresq,dtb,ni,nj,nk,ptbr,rbr,pi,p,ptpp,qvp,qrp,    &
     &             ptpf,qvf,qrf)

! -----

! Setting ptpf-ptpf' to dqdt
! Setting dqdt+dqadt' to dqadt

      call s_dmpcld( ni, nj, nk, ptpf, 'dqdt ', 'm', 1.0e0 )
      call s_dmpcld( ni, nj, nk, tmpdmp, 'dqdt ', 'r', 1.0e0 )
      call s_dmpcld( ni, nj, nk, tmpdmp, 'dqadt', 'p', 1.0e0 )

! -----

!! -----

! Perform the fall out of the rain water.

      call fallqr(iddz,dtb,ni,nj,nk,jcb,rbr,rst,urq,qrf,prr,tmp1)

! -----

      end subroutine s_warmblk

!-----7--------------------------------------------------------------7--

      end module m_warmblk
