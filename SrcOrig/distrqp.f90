!***********************************************************************
      module m_distrqp
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/09/09
!     Modification: 2002/12/02, 2003/04/30, 2003/05/19, 2003/12/12,
!                   2004/08/20, 2004/09/01, 2006/02/13, 2006/09/21,
!                   2007/05/07, 2007/06/27, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2011/09/22, 2013/01/28,
!                   2013/02/11, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     distribute the observed precipitation to the bulk categories.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar
      use m_temparam

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: distrqp, s_distrqp

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface distrqp

        module procedure s_distrqp

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_distrqp(fpdatype_rdr,fpcphopt,fphaiopt,fpthresq,     &
     &                     ni,nj,nk,nqw,nqi,rbr,qwtr,qice,              &
     &                     qprdr,qwrdr,qirdr)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpdatype_rdr
                       ! Formal parameter of unique index of datype_rdr

      integer, intent(in) :: fpcphopt
                       ! Formal parameter of unique index of cphopt

      integer, intent(in) :: fphaiopt
                       ! Formal parameter of unique index of haiopt

      integer, intent(in) :: fpthresq
                       ! Formal parameter of unique index of thresq

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nqi
                       ! Number of categories of ice hydrometeor

      real, intent(in) :: rbr(0:ni+1,0:nj+1,1:nk)
                       ! Base state density

      real, intent(in) :: qwtr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor

      real, intent(in) :: qice(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor

! Input and output variable

      real, intent(inout) :: qprdr(0:ni+1,0:nj+1,1:nk)
                       ! Precipitation mixing ratio of radar data
                       ! at marked time

! Output variables

      real, intent(out) :: qwrdr(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Water hydrometeor of radar data at marked time

      real, intent(out) :: qirdr(0:ni+1,0:nj+1,1:nk,1:nqi)
                       ! Ice hydrometeor of radar data at marked time

! Internal shared variables

      character(len=108) datype_rdr
                       ! Control flag of radar data type

      integer cphopt   ! Option for cloud micro physics
      integer haiopt   ! Option for additional hail processes

      real thresq      ! Minimum threshold value of mixing ratio

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real qpsum       ! Total precipitation

      real qpsumv      ! Inverse of qpsum

!-----7--------------------------------------------------------------7--

! Initialize character variable.

      call inichar(datype_rdr)

! -----

! Get the required namelist variables.

      call getcname(fpdatype_rdr,datype_rdr)
      call getiname(fpcphopt,cphopt)
      call getiname(fphaiopt,haiopt)
      call getrname(fpthresq,thresq)

! -----

!!! Distribute the observed precipitation to the bulk categories.

!$omp parallel default(shared) private(k)

! Change measurement.

      if(datype_rdr(1:1).eq.'r') then

        do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

          do j=1,nj-1
          do i=1,ni-1

            if(qprdr(i,j,k).lt.qpmin.and.qprdr(i,j,k).gt.lim34n) then
              qprdr(i,j,k)=0.e0
            end if

            qprdr(i,j,k)=qprdr(i,j,k)/rbr(i,j,k)

          end do
          end do

!$omp end do

        end do

      end if

! -----

!! Perform distribution.

      if(abs(cphopt).lt.10) then

! Distribute the observed precipitation to the rain mixing ratio.

        if(abs(cphopt).eq.1) then

          do k=2,nk-2

!$omp do schedule(runtime) private(i,j)

            do j=1,nj-1
            do i=1,ni-1

              if(qprdr(i,j,k).le.lim34n) then

                qwrdr(i,j,k,2)=lim35n

              else

!ORIG           if(qwtr(i,j,k,2).gt.thresq) then
!ORIG             qwrdr(i,j,k,2)=qprdr(i,j,k)
!ORIG           else
!ORIG             qwrdr(i,j,k,2)=0.e0
!ORIG           end if

                qwrdr(i,j,k,2)=qprdr(i,j,k)

              end if

            end do
            end do

!$omp end do

          end do

! -----

! Distribute the observed precipitation to the rain, snow, graupel and
! hail mixing ratio.

        else if(abs(cphopt).ge.2) then

          if(haiopt.eq.0) then

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j,qpsum,qpsumv)

              do j=1,nj-1
              do i=1,ni-1

                if(qprdr(i,j,k).le.lim34n) then

                  qwrdr(i,j,k,2)=lim35n
                  qirdr(i,j,k,2)=lim35n
                  qirdr(i,j,k,3)=lim35n

                else

                  qpsum=qwtr(i,j,k,2)+qice(i,j,k,2)+qice(i,j,k,3)

                  if(qpsum.gt.thresq) then

                    qpsumv=1.e0/qpsum

                    qwrdr(i,j,k,2)=qprdr(i,j,k)*qwtr(i,j,k,2)*qpsumv
                    qirdr(i,j,k,2)=qprdr(i,j,k)*qice(i,j,k,2)*qpsumv
                    qirdr(i,j,k,3)=qprdr(i,j,k)*qice(i,j,k,3)*qpsumv

                  else

!ORIG               qwrdr(i,j,k,2)=0.e0
!ORIG               qirdr(i,j,k,2)=0.e0
!ORIG               qirdr(i,j,k,3)=0.e0

                    qwrdr(i,j,k,2)=qprdr(i,j,k)
                    qirdr(i,j,k,2)=0.e0
                    qirdr(i,j,k,3)=0.e0

                  end if

                end if

              end do
              end do

!$omp end do

            end do

          else

            do k=2,nk-2

!$omp do schedule(runtime) private(i,j,qpsum,qpsumv)

              do j=1,nj-1
              do i=1,ni-1

                if(qprdr(i,j,k).le.lim34n) then

                  qwrdr(i,j,k,2)=lim35n
                  qirdr(i,j,k,2)=lim35n
                  qirdr(i,j,k,3)=lim35n
                  qirdr(i,j,k,4)=lim35n

                else

                  qpsum=qwtr(i,j,k,2)                                   &
     &              +qice(i,j,k,2)+qice(i,j,k,3)+qice(i,j,k,4)

                  if(qpsum.gt.thresq) then

                    qpsumv=1.e0/qpsum

                    qwrdr(i,j,k,2)=qprdr(i,j,k)*qwtr(i,j,k,2)*qpsumv
                    qirdr(i,j,k,2)=qprdr(i,j,k)*qice(i,j,k,2)*qpsumv
                    qirdr(i,j,k,3)=qprdr(i,j,k)*qice(i,j,k,3)*qpsumv
                    qirdr(i,j,k,4)=qprdr(i,j,k)*qice(i,j,k,4)*qpsumv

                  else

!ORIG               qwrdr(i,j,k,2)=0.e0
!ORIG               qirdr(i,j,k,2)=0.e0
!ORIG               qirdr(i,j,k,3)=0.e0
!ORIG               qirdr(i,j,k,4)=0.e0

                    qwrdr(i,j,k,2)=qprdr(i,j,k)
                    qirdr(i,j,k,2)=0.e0
                    qirdr(i,j,k,3)=0.e0
                    qirdr(i,j,k,4)=0.e0

                  end if

                end if

              end do
              end do

!$omp end do

            end do

          end if

        end if

! -----

      end if

!! -----

!$omp end parallel

!!! -----

      end subroutine s_distrqp

!-----7--------------------------------------------------------------7--

      end module m_distrqp
