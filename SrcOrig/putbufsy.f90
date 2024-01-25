!***********************************************************************
      module m_putbufsy
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/25
!     Modification: 1999/03/25, 1999/04/06, 1999/07/05, 1999/08/18,
!                   1999/08/23, 1999/10/07, 1999/10/12, 1999/11/01,
!                   1999/11/19, 2000/01/17, 2000/04/18, 2000/07/05,
!                   2001/05/29, 2002/04/02, 2002/06/06, 2003/04/30,
!                   2003/05/19, 2004/08/20, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/10/19, 2008/05/02, 2008/08/25,
!                   2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     put the sending buffer in y direction for sub domain.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_getiname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: putbufsy, s_putbufsy

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface putbufsy

        module procedure s_putbufsy

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
      subroutine s_putbufsy(fpsbc,fpnbc,fproc,jssnd,jnsnd,ni,nj,kmax,   &
     &                      var,ib,nb,sbufy)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: jssnd
                       ! Array index of south sended value

      integer, intent(in) :: jnsnd
                       ! Array index of north sended value

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: kmax
                       ! Maximum array index in z direction

      integer, intent(in) :: ib
                       ! Exchanging variables number

      integer, intent(in) :: nb
                       ! Number of exchanging variables

      real, intent(in) :: var(0:ni+1,0:nj+1,1:kmax)
                       ! Optional exchanged variable

! Output variable

      real, intent(out) :: sbufy(0:ni+1,1:kmax,1:2*nb)
                       ! Sending buffer in y direction

! Internal shared variables

      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

! Internal private variables

      integer i        ! Array index in x direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)

! -----

!! Put the sending buffer in y direction.

      if(njsub.ge.2) then

        if(fproc(1:3).eq.'all'.or.(sbc.eq.-1.and.nbc.eq.-1)) then

!$omp parallel default(shared) private(k)

! Fill in the sending buffer with the value in the south halo regions.

          if(sbc.eq.-1) then

            if(fproc(1:3).eq.'all') then

              do k=1,kmax

!$omp do schedule(runtime) private(i)

                do i=0,ni+1
                  sbufy(i,k,ib)=var(i,jssnd,k)
                end do

!$omp end do

              end do

            else if(fproc(1:3).eq.'bnd') then

              if(jsub.eq.0) then

                do k=1,kmax

!$omp do schedule(runtime) private(i)

                  do i=0,ni+1
                    sbufy(i,k,ib)=var(i,jssnd,k)
                  end do

!$omp end do

                end do

              end if

            end if

          else if(sbc.ne.-1) then

            if(jsub.ne.0) then

              do k=1,kmax

!$omp do schedule(runtime) private(i)

                do i=0,ni+1
                  sbufy(i,k,ib)=var(i,jssnd,k)
                end do

!$omp end do

              end do

            end if

          end if

! -----

! Fill in the sending buffer with the value in the north halo regions.

          if(nbc.eq.-1) then

            if(fproc(1:3).eq.'all') then

              do k=1,kmax

!$omp do schedule(runtime) private(i)

                do i=0,ni+1
                  sbufy(i,k,nb+ib)=var(i,jnsnd,k)
                end do

!$omp end do

              end do

            else if(fproc(1:3).eq.'bnd') then

              if(jsub.eq.njsub-1) then

                do k=1,kmax

!$omp do schedule(runtime) private(i)

                  do i=0,ni+1
                    sbufy(i,k,nb+ib)=var(i,jnsnd,k)
                  end do

!$omp end do

                end do

              end if

            end if

          else if(nbc.ne.-1) then

            if(jsub.ne.njsub-1) then

              do k=1,kmax

!$omp do schedule(runtime) private(i)

                do i=0,ni+1
                  sbufy(i,k,nb+ib)=var(i,jnsnd,k)
                end do

!$omp end do

              end do

            end if

          end if

! -----

!$omp end parallel

        end if

      end if

!! -----

      end subroutine s_putbufsy

!-----7--------------------------------------------------------------7--

      end module m_putbufsy
