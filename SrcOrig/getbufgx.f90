!***********************************************************************
      module m_getbufgx
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/12/04
!     Modification: 2007/01/05, 2007/01/20, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the receiving buffer in x direction for group domain.

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

      public :: getbufgx, s_getbufgx

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getbufgx

        module procedure s_getbufgx

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
      subroutine s_getbufgx(fpwbc,fpebc,fproc,iwrcv,iercv,ni,nj,kmax,   &
     &                      var,ib,nb,rbufx)
!***********************************************************************

! Input variables

      character(len=3), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: iwrcv
                       ! Array index of west received value

      integer, intent(in) :: iercv
                       ! Array index of east received value

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

      real, intent(in) :: rbufx(0:nj+1,1:kmax,1:2*nb)
                       ! Receiving buffer in x direction

! Output variable

      real, intent(out) :: var(0:ni+1,0:nj+1,1:kmax)
                       ! Optional exchanged variable

! Internal shared variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundaty conditions

! Internal private variables

      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)

! -----

!! Get the receiving buffer in x direction.

      if(nigrp.ge.2.and.(wbc.ne.-1.and.ebc.ne.-1)) then

        if(fproc(1:3).eq.'all'.or.(wbc.eq.1.and.ebc.eq.1)) then

!$omp parallel default(shared) private(k)

! Fill in the west halo regions with the received value.

          if(isub.eq.0) then

            if(wbc.eq.1) then

              if(fproc(1:3).eq.'all') then

                do k=1,kmax

!$omp do schedule(runtime) private(j)

                  do j=0,nj+1
                    var(iwrcv,j,k)=rbufx(j,k,ib)
                  end do

!$omp end do

                end do

              else if(fproc(1:3).eq.'bnd') then

                if(igrp.eq.0) then

                  do k=1,kmax

!$omp do schedule(runtime) private(j)

                    do j=0,nj+1
                      var(iwrcv,j,k)=rbufx(j,k,ib)
                    end do

!$omp end do

                  end do

                end if

              end if

            else if(wbc.ne.1) then

              if(ebw.eq.0) then

                do k=1,kmax

!$omp do schedule(runtime) private(j)

                  do j=0,nj+1
                    var(iwrcv,j,k)=rbufx(j,k,ib)
                  end do

!$omp end do

                end do

              end if

            end if

          end if

! -----

! Fill in the east halo regions with the received value.

          if(isub.eq.nisub-1) then

            if(ebc.eq.1) then

              if(fproc(1:3).eq.'all') then

                do k=1,kmax

!$omp do schedule(runtime) private(j)

                  do j=0,nj+1
                    var(iercv,j,k)=rbufx(j,k,nb+ib)
                  end do

!$omp end do

                end do

              else if(fproc(1:3).eq.'bnd') then

                if(igrp.eq.nigrp-1) then

                  do k=1,kmax

!$omp do schedule(runtime) private(j)

                    do j=0,nj+1
                      var(iercv,j,k)=rbufx(j,k,nb+ib)
                    end do

!$omp end do

                  end do

                end if

              end if

            else if(ebc.ne.1) then

              if(ebe.eq.0) then

                do k=1,kmax

!$omp do schedule(runtime) private(j)

                  do j=0,nj+1
                    var(iercv,j,k)=rbufx(j,k,nb+ib)
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

      end subroutine s_getbufgx

!-----7--------------------------------------------------------------7--

      end module m_getbufgx
