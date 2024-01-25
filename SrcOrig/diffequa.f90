!***********************************************************************
      module m_diffequa
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2000/06/01
!     Modification: 2000/07/05, 2000/12/18, 2001/02/13, 2001/04/15,
!                   2001/05/29, 2001/06/06, 2001/08/07, 2001/11/20,
!                   2002/04/02, 2002/06/06, 2003/02/13, 2003/04/30,
!                   2003/05/19, 2003/10/31, 2003/11/05, 2003/12/12,
!                   2004/03/05, 2004/08/01, 2004/08/31, 2004/09/01,
!                   2004/09/25, 2005/02/10, 2006/05/12, 2006/08/08,
!                   2006/12/04, 2007/01/05, 2007/01/20, 2007/01/20,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/11/13, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     solve the parabolic pertial differential equation.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bc4news
      use m_bcphi
      use m_bcycle
      use m_chkitr
      use m_combuf
      use m_comindx
      use m_commath
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getrname
      use m_putbufgx
      use m_putbufgy
      use m_putbufsx
      use m_putbufsy
      use m_shiftgx
      use m_shiftgy
      use m_shiftsy
      use m_shiftsx

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: diffequa, s_diffequa

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface diffequa

        module procedure s_diffequa

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_diffequa(fpdx,fpdy,fpdz,iteps,alpha1,alpha2,itcnt,   &
     &                      ni,nj,nk,known,phi,dphi)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdx
                       ! Formal parameter of unique index of dx

      integer, intent(in) :: fpdy
                       ! Formal parameter of unique index of dy

      integer, intent(in) :: fpdz
                       ! Formal parameter of unique index of dz

      integer, intent(in) :: ni
                       ! Model dimensionn in x direction

      integer, intent(in) :: nj
                       ! Model dimensionn in y direction

      integer, intent(in) :: nk
                       ! Model dimensionn in z direction

      real, intent(in) :: iteps
                       ! Value of convergence of iteration

      real, intent(in) :: alpha1
                       ! Weighting coeffient

      real, intent(in) :: alpha2
                       ! Weighting coeffient

      real, intent(in) :: known(0:ni+1,0:nj+1,1:nk)
                       ! Known quantity

! Input and output variables

      integer, intent(inout) :: itcnt
                       ! Iteration count

      real, intent(inout) :: phi(0:ni+1,0:nj+1,1:nk)
                       ! Optional solved variable

! Internal shared variables

      real dx          ! Grid distance in x direction
      real dy          ! Grid distance in y direction
      real dz          ! Grid distance in z direction

      real dxivsq      ! dxiv x dxiv
      real dyivsq      ! dyiv x dyiv

      real adzivs      ! (alpha1 / alpha2 x dziv)^2

      real dtitr       ! Time interval of iteration

      real itcon       ! Control flag of continuation of iteration

      real, intent(inout) :: dphi(0:ni+1,0:nj+1,1:nk)
                       ! Difference of phi in time interval

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in y direction

      real phi2        ! 2.0 x phi

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getrname(fpdx,dx)
      call getrname(fpdy,dy)
      call getrname(fpdz,dz)

! -----

! Set the common used variables.

      dxivsq=1.e0/(dx*dx)
      dyivsq=1.e0/(dy*dy)

      adzivs=alpha1/(alpha2*dz)
      adzivs=adzivs*adzivs

! -----

! Get the time interval of iteration.

      dtitr=min(.3e0*dx*dx,.3e0*dy*dy,.3e0/adzivs)

! -----

!! Solve the parabolic pertial differential equation.

      iterate: do

! Inclimemt the iteration counter.

        itcnt=itcnt+1

! -----

! Calculate the Laplacian terms and add the known quantity terms in the
! right hand and get new phi.

!$omp parallel default(shared) private(k)

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j,phi2)

          do j=2,nj-2
          do i=2,ni-2
            phi2=2.e0*phi(i,j,k)

            dphi(i,j,k)=known(i,j,k)                                    &
     &        +(phi(i+1,j,k)+phi(i-1,j,k)-phi2)*dxivsq                  &
     &        +(phi(i,j+1,k)+phi(i,j-1,k)-phi2)*dyivsq                  &
     &        +(phi(i,j,k+1)+phi(i,j,k-1)-phi2)*adzivs

            phi(i,j,k)=phi(i,j,k)+dtitr*dphi(i,j,k)

          end do
          end do

!$omp end do

        end do

!$omp end parallel

! -----

!! Exchange the value horizontally.

! Exchange the value horizontally between sub domain.

        call s_putbufsx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,phi,1,1,sbuf)

        call s_shiftsx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufsx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,phi,1,1,rbuf)

        call s_putbufsy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,phi,1,1,sbuf)

        call s_shiftsy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

        call s_getbufsy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,phi,1,1,rbuf)

! -----

! Exchange the value horizontally between group domain.

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,phi,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,phi,1,1,rbuf)

        call s_putbufgy(idsbc,idnbc,'all',2,nj-2,ni,nj,nk,phi,1,1,sbuf)

        call s_shiftgy(idsbc,idnbc,'all',ni,nk,1,sbuf,rbuf)

        call s_getbufgy(idsbc,idnbc,'all',1,nj-1,ni,nj,nk,phi,1,1,rbuf)

        call s_putbufgx(idwbc,idebc,'all',2,ni-2,ni,nj,nk,phi,1,1,sbuf)

        call s_shiftgx(idwbc,idebc,'all',nj,nk,1,sbuf,rbuf)

        call s_getbufgx(idwbc,idebc,'all',1,ni-1,ni,nj,nk,phi,1,1,rbuf)

! -----

!! -----

! Set the periodic boundary conditions.

        call bcycle(idwbc,idebc,idsbc,idnbc,                            &
     &              2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,phi)

! -----

! Set the boundary conditions at the four corners.

        call bc4news(idwbc,idebc,idsbc,idnbc,1,ni-1,1,nj-1,ni,nj,nk,phi)

! -----

! Set the bottom and top boundary conditions.

        call bcphi(ni,nj,nk,phi)

! -----

! Check the convergence of the iteration.

        call chkitr('common',2,ni-2,2,nj-2,2,nk-2,itcon,ni,nj,nk,dphi)

! -----

! Stop iteration.

        if(itcon.lt.iteps) then

          exit iterate

        end if

! -----

      end do iterate

!! -----

      end subroutine s_diffequa

!-----7--------------------------------------------------------------7--

      end module m_diffequa
