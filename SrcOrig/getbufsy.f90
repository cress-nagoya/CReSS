!***********************************************************************
      module m_getbufsy
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/25
!     Modification: 1999/03/25, 1999/04/06, 1999/07/05, 1999/08/18,
!                   1999/08/23, 1999/10/07, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/04/18, 2000/07/05, 2001/05/29,
!                   2002/04/02, 2002/06/06, 2003/04/30, 2003/05/19,
!                   2004/08/20, 2006/12/04, 2007/01/05, 2007/01/20,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the receiving buffer in y direction for sub domain.

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

      public :: getbufsy, s_getbufsy

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getbufsy

        module procedure s_getbufsy

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
      subroutine s_getbufsy(fpsbc,fpnbc,fproc,jsrcv,jnrcv,ni,nj,kmax,   &
     &                      var,ib,nb,rbufy)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: jsrcv
                       ! Array index of south received value

      integer, intent(in) :: jnrcv
                       ! Array index of north received value

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

      real, intent(in) :: rbufy(0:ni+1,1:kmax,1:2*nb)
                       ! Receiving buffer in y direction

! Output variable

      real, intent(out) :: var(0:ni+1,0:nj+1,1:kmax)
                       ! Optional exchanged variable

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

!! Get the receiving buffer in y direction.

      if(njsub.ge.2) then

        if(fproc(1:3).eq.'all'.or.(sbc.eq.-1.and.nbc.eq.-1)) then

!$omp parallel default(shared) private(k)

! Fill in the south halo regions with the received value.

          if(sbc.eq.-1) then

            if(fproc(1:3).eq.'all') then

              do k=1,kmax

!$omp do schedule(runtime) private(i)

                do i=0,ni+1
                  var(i,jsrcv,k)=rbufy(i,k,ib)
                end do

!$omp end do

              end do

            else if(fproc(1:3).eq.'bnd') then

              if(jsub.eq.0) then

                do k=1,kmax

!$omp do schedule(runtime) private(i)

                  do i=0,ni+1
                    var(i,jsrcv,k)=rbufy(i,k,ib)
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
                  var(i,jsrcv,k)=rbufy(i,k,ib)
                end do

!$omp end do

              end do

            end if

          end if

! -----

! Fill in the north halo regions with the received value.

          if(nbc.eq.-1) then

            if(fproc(1:3).eq.'all') then

              do k=1,kmax

!$omp do schedule(runtime) private(i)

                do i=0,ni+1
                  var(i,jnrcv,k)=rbufy(i,k,nb+ib)
                end do

!$omp end do

              end do

            else if(fproc(1:3).eq.'bnd') then

              if(jsub.eq.njsub-1) then

                do k=1,kmax

!$omp do schedule(runtime) private(i)

                  do i=0,ni+1
                    var(i,jnrcv,k)=rbufy(i,k,nb+ib)
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
                  var(i,jnrcv,k)=rbufy(i,k,nb+ib)
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

      end subroutine s_getbufsy

!-----7--------------------------------------------------------------7--

      end module m_getbufsy
