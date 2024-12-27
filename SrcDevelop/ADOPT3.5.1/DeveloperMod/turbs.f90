!***********************************************************************
      module m_turbs
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/07/05
!     Modification: 1999/07/21, 1999/07/28, 1999/08/03, 1999/09/30,
!                   1999/10/12, 1999/11/01, 2000/01/17, 2000/11/17,
!                   2001/06/06, 2002/01/21, 2002/04/02, 2003/01/04,
!                   2003/03/13, 2003/03/21, 2003/04/30, 2003/05/19,
!                   2003/11/28, 2004/02/01, 2004/06/10, 2006/11/06,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2011/08/09, 2013/01/28, 2013/02/13

!     Author      : Satoki Tsujino
!     Modification: 2024/12/25

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     calculate optional scalar turbulent mixing.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getcname
      use m_getiname
      use m_getrname
      use m_comtub

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: turbs, s_turbs

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface turbs

        module procedure s_turbs

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
      subroutine s_turbs(fptrnopt,fpmpopt,fpmfcopt,fpdxiv,fpdyiv,fpdziv,& 
     &                   fpdmpvar,ni,nj,nk,dflag,j31,j32,               &
     &                   jcb8u,jcb8v,mf,rmf,rmf8u,rmf8v,h1,h2,h3,sfrc,  &
     &                   tmp1,tmp2,tmp3)
!***********************************************************************

      use m_comtub     ! adding by satoki

! Input variables

      integer, intent(in) :: fptrnopt
                       ! Formal parameter of unique index of trnopt

      integer, intent(in) :: fpmpopt
                       ! Formal parameter of unique index of mpopt

      integer, intent(in) :: fpmfcopt
                       ! Formal parameter of unique index of mfcopt

      integer, intent(in) :: fpdxiv
                       ! Formal parameter of unique index of dxiv

      integer, intent(in) :: fpdyiv
                       ! Formal parameter of unique index of dyiv

      integer, intent(in) :: fpdziv
                       ! Formal parameter of unique index of dziv

      integer, intent(in) :: fpdmpvar
                       ! Formal parameter of unique index of dmpvar

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      character(len=2), intent(in) :: dflag
                       ! flag of turbulence dump variable (scalar)

      real, intent(in) :: j31(0:ni+1,0:nj+1,1:nk)
                       ! z-x components of Jacobian

      real, intent(in) :: j32(0:ni+1,0:nj+1,1:nk)
                       ! z-y components of Jacobian

      real, intent(in) :: jcb8u(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at u points

      real, intent(in) :: jcb8v(0:ni+1,0:nj+1,1:nk)
                       ! Jacobian at v points

      real, intent(in) :: mf(0:ni+1,0:nj+1)
                       ! Map scale factors

      real, intent(in) :: rmf(0:ni+1,0:nj+1,1:4)
                       ! Related parameters of map scale factors

      real, intent(in) :: rmf8u(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at u points

      real, intent(in) :: rmf8v(0:ni+1,0:nj+1,1:3)
                       ! Related parameters of map scale factors
                       ! at v points

      real, intent(in) :: h1(0:ni+1,0:nj+1,1:nk)
                       ! x components of turbulent fluxes

      real, intent(in) :: h2(0:ni+1,0:nj+1,1:nk)
                       ! y components of turbulent fluxes

      real, intent(in) :: h3(0:ni+1,0:nj+1,1:nk)
                       ! z components of turbulent fluxes

! Input and output variable

      real, intent(inout) :: sfrc(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar forcing term

! Internal shared variables

      integer trnopt   ! Option for terrain height setting
      integer mpopt    ! Option for map projection
      integer mfcopt   ! Option for map scale factor

      real dxiv        ! Inverse of dx
      real dyiv        ! Inverse of dy
      real dziv        ! Inverse of dz

      character(len=108) dmpvar
                       ! Control flag of dump variables

      real, intent(inout) :: tmp1(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp2(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

      real, intent(inout) :: tmp3(0:ni+1,0:nj+1,1:nk)
                       ! Temporary array

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fptrnopt,trnopt)
      call getiname(fpmpopt,mpopt)
      call getiname(fpmfcopt,mfcopt)
      call getrname(fpdxiv,dxiv)
      call getrname(fpdyiv,dyiv)
      call getrname(fpdziv,dziv)
      call getcname(fpdmpvar,dmpvar)

! -----

      if(dmpvar(16:16).eq.'o'.or.dmpvar(16:16).eq.'+')then

        if(dflag(1:2)=='pt')then

          turbpt=sfrc

        !else if(dflag(1:2)=='qv')then

        !  turbqv=sfrc

        end if

      end if

! -----

! Calculate the scalar turbulent mixing.

!$omp parallel default(shared) private(k)

      if(trnopt.eq.0) then

        if(mfcopt.eq.0) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              tmp1(i,j,k)=jcb8u(i,j,k)*h1(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              tmp2(i,j,k)=jcb8v(i,j,k)*h2(i,j,k)
            end do
            end do

!$omp end do

          end do

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              sfrc(i,j,k)=sfrc(i,j,k)+((h3(i,j,k+1)-h3(i,j,k))*dziv     &
     &          +((tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv                      &
     &          +(tmp2(i,j+1,k)-tmp2(i,j,k))*dyiv))
            end do
            end do

!$omp end do

          end do

        else

          if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

            if(mpopt.eq.0.or.mpopt.eq.10) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  tmp1(i,j,k)=jcb8u(i,j,k)*h1(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  tmp2(i,j,k)=rmf8v(i,j,2)*jcb8v(i,j,k)*h2(i,j,k)
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  tmp1(i,j,k)=rmf8u(i,j,2)*jcb8u(i,j,k)*h1(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  tmp2(i,j,k)=jcb8v(i,j,k)*h2(i,j,k)
                end do
                end do

!$omp end do

              end do

            end if

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                sfrc(i,j,k)=sfrc(i,j,k)+((h3(i,j,k+1)-h3(i,j,k))*dziv   &
     &            +mf(i,j)*((tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv            &
     &            +(tmp2(i,j+1,k)-tmp2(i,j,k))*dyiv))
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                tmp1(i,j,k)=rmf8u(i,j,2)*jcb8u(i,j,k)*h1(i,j,k)
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                tmp2(i,j,k)=rmf8v(i,j,2)*jcb8v(i,j,k)*h2(i,j,k)
              end do
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                sfrc(i,j,k)=sfrc(i,j,k)+((h3(i,j,k+1)-h3(i,j,k))*dziv   &
     &            +rmf(i,j,1)*((tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv         &
     &            +(tmp2(i,j+1,k)-tmp2(i,j,k))*dyiv))
              end do
              end do

!$omp end do

            end do

          end if

        end if

      else

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-1
            tmp1(i,j,k)=j31(i,j,k)*(h1(i,j,k-1)+h1(i,j,k))
          end do
          end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-1
          do i=2,ni-2
            tmp2(i,j,k)=j32(i,j,k)*(h2(i,j,k-1)+h2(i,j,k))
          end do
          end do

!$omp end do

        end do

        do k=2,nk-1

!$omp do schedule(runtime) private(i,j)

          do j=2,nj-2
          do i=2,ni-2
            tmp3(i,j,k)=h3(i,j,k)+.25e0                                 &
     &        *((tmp1(i,j,k)+tmp1(i+1,j,k))+(tmp2(i,j,k)+tmp2(i,j+1,k)))
          end do
          end do

!$omp end do

        end do

        if(mfcopt.eq.0) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-1
              tmp1(i,j,k)=jcb8u(i,j,k)*h1(i,j,k)
            end do
            end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-1
            do i=2,ni-2
              tmp2(i,j,k)=jcb8v(i,j,k)*h2(i,j,k)
            end do
            end do

!$omp end do

          end do

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=2,nj-2
            do i=2,ni-2
              sfrc(i,j,k)=sfrc(i,j,k)+((tmp3(i,j,k+1)-tmp3(i,j,k))*dziv &
     &          +((tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv                      &
     &          +(tmp2(i,j+1,k)-tmp2(i,j,k))*dyiv))
            end do
            end do

!$omp end do

          end do

        else

          if(mpopt.eq.0.or.mpopt.eq.5.or.mpopt.eq.10) then

            if(mpopt.eq.0.or.mpopt.eq.10) then

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  tmp1(i,j,k)=jcb8u(i,j,k)*h1(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  tmp2(i,j,k)=rmf8v(i,j,2)*jcb8v(i,j,k)*h2(i,j,k)
                end do
                end do

!$omp end do

              end do

            else

              do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-2
                do i=2,ni-1
                  tmp1(i,j,k)=rmf8u(i,j,2)*jcb8u(i,j,k)*h1(i,j,k)
                end do
                end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

                do j=2,nj-1
                do i=2,ni-2
                  tmp2(i,j,k)=jcb8v(i,j,k)*h2(i,j,k)
                end do
                end do

!$omp end do

              end do

            end if

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                sfrc(i,j,k)=sfrc(i,j,k)                                 &
     &            +((tmp3(i,j,k+1)-tmp3(i,j,k))*dziv                    &
     &            +mf(i,j)*((tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv            &
     &            +(tmp2(i,j+1,k)-tmp2(i,j,k))*dyiv))
              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-1
                tmp1(i,j,k)=rmf8u(i,j,2)*jcb8u(i,j,k)*h1(i,j,k)
              end do
              end do

!$omp end do

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-1
              do i=2,ni-2
                tmp2(i,j,k)=rmf8v(i,j,2)*jcb8v(i,j,k)*h2(i,j,k)
              end do
              end do

!$omp end do

            end do

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

              do j=2,nj-2
              do i=2,ni-2
                sfrc(i,j,k)=sfrc(i,j,k)                                 &
     &            +((tmp3(i,j,k+1)-tmp3(i,j,k))*dziv                    &
     &            +rmf(i,j,1)*((tmp1(i+1,j,k)-tmp1(i,j,k))*dxiv         &
     &            +(tmp2(i,j+1,k)-tmp2(i,j,k))*dyiv))
              end do
              end do

!$omp end do

            end do

          end if

        end if

      end if

      if(dmpvar(16:16).eq.'o'.or.dmpvar(16:16).eq.'+')then

        if(dflag(1:2).eq.'pt')then

!$omp do schedule(runtime) private(i,j,k)

          do k=1,nk-1
          do j=1,nj-1
          do i=1,ni-1
            turbpt(i,j,k)=sfrc(i,j,k)-turbpt(i,j,k)
          end do
          end do
          end do

!$omp end do

!        else if(dflag(1:2).eq.'qv')then

!!$omp do schedule(runtime) private(i,j,k)

!          do k=1,nk-1
!          do j=1,nj-1
!          do i=1,ni-1
!            turbqv(i,j,k)=sfrc(i,j,k)-turbqv(i,j,k)
!          end do
!          end do
!          end do

!!$omp end do

        end if

      end if

!$omp end parallel

! -----

      end subroutine s_turbs

!-----7--------------------------------------------------------------7--

      end module m_turbs
