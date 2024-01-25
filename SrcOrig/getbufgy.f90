!***********************************************************************
      module m_getbufgy
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/12/04
!     Modification: 2007/01/05, 2007/01/20, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the receiving buffer in y direction for group domain.

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

      public :: getbufgy, s_getbufgy

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getbufgy

        module procedure s_getbufgy

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
      subroutine s_getbufgy(fpsbc,fpnbc,fproc,jsrcv,jnrcv,ni,nj,kmax,   &
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

      if(njgrp.ge.2.and.(sbc.ne.-1.and.nbc.ne.-1)) then

        if(fproc(1:3).eq.'all'.or.(sbc.eq.1.and.nbc.eq.1)) then

!$omp parallel default(shared) private(k)

! Fill in the south halo regions with the received value.

          if(jsub.eq.0) then

            if(sbc.eq.1) then

              if(fproc(1:3).eq.'all') then

                do k=1,kmax

!$omp do schedule(runtime) private(i)

                  do i=0,ni+1
                    var(i,jsrcv,k)=rbufy(i,k,ib)
                  end do

!$omp end do

                end do

              else if(fproc(1:3).eq.'bnd') then

                if(jgrp.eq.0) then

                  do k=1,kmax

!$omp do schedule(runtime) private(i)

                    do i=0,ni+1
                      var(i,jsrcv,k)=rbufy(i,k,ib)
                    end do

!$omp end do

                  end do

                end if

              end if

            else if(sbc.ne.1) then

              if(ebs.eq.0) then

                do k=1,kmax

!$omp do schedule(runtime) private(i)

                  do i=0,ni+1
                    var(i,jsrcv,k)=rbufy(i,k,ib)
                  end do

!$omp end do

                end do

              end if

            end if

          end if

! -----

! Fill in the north halo regions with the received value.

          if(jsub.eq.njsub-1) then

            if(nbc.eq.1) then

              if(fproc(1:3).eq.'all') then

                do k=1,kmax

!$omp do schedule(runtime) private(i)

                  do i=0,ni+1
                    var(i,jnrcv,k)=rbufy(i,k,nb+ib)
                  end do

!$omp end do

                end do

              else if(fproc(1:3).eq.'bnd') then

                if(jgrp.eq.njgrp-1) then

                  do k=1,kmax

!$omp do schedule(runtime) private(i)

                    do i=0,ni+1
                      var(i,jnrcv,k)=rbufy(i,k,nb+ib)
                    end do

!$omp end do

                  end do

                end if

              end if

            else if(nbc.ne.1) then

              if(ebn.eq.0) then

                do k=1,kmax

!$omp do schedule(runtime) private(i)

                  do i=0,ni+1
                    var(i,jnrcv,k)=rbufy(i,k,nb+ib)
                  end do

!$omp end do

                end do

              end if

            end if

          end if

! -----

!$omp end parallel

        end if

      end if

!! -----

      end subroutine s_getbufgy

!-----7--------------------------------------------------------------7--

      end module m_getbufgy
