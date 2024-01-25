!***********************************************************************
      module m_emitqt
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/05/31
!     Modification: 2004/09/10, 2007/01/31, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2011/05/16, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     emit the tracer mixing ratio from user specified location.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_comindx
      use m_commath
      use m_getiname
      use m_getrname
      use m_getxy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: emitqt, s_emitqt

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface emitqt

        module procedure s_emitqt

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic cos
      intrinsic real
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_emitqt(fpqt0opt,fpqt0num,fpqt0rx,fpqt0ry,fpqt0rz,    &
     &                    fpqt0cx,fpqt0cy,fpqt0cz,fpqt0ds,fpqtdt,       &
     &                    ni,nj,nk,zph,qtfrc,xs,ys)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpqt0opt
                       ! Formal parameter of unique index of qt0opt

      integer, intent(in) :: fpqt0num
                       ! Formal parameter of unique index of qt0num

      integer, intent(in) :: fpqt0rx
                       ! Formal parameter of unique index of qt0rx

      integer, intent(in) :: fpqt0ry
                       ! Formal parameter of unique index of qt0ry

      integer, intent(in) :: fpqt0rz
                       ! Formal parameter of unique index of qt0rz

      integer, intent(in) :: fpqt0cx
                       ! Formal parameter of unique index of qt0cx

      integer, intent(in) :: fpqt0cy
                       ! Formal parameter of unique index of qt0cy

      integer, intent(in) :: fpqt0cz
                       ! Formal parameter of unique index of qt0cz

      integer, intent(in) :: fpqt0ds
                       ! Formal parameter of unique index of qt0ds

      integer, intent(in) :: fpqtdt
                       ! Formal parameter of unique index of qtdt

      integer, intent(in) :: ni
                       ! Model demension in x direction

      integer, intent(in) :: nj
                       ! Model demension in y direction

      integer, intent(in) :: nk
                       ! Model demension in z direction

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

! Input and output variable

      real, intent(inout) :: qtfrc(0:ni+1,0:nj+1,1:nk)
                       ! Forcing term in tracer equation

! Internal shared variables

      integer qt0opt   ! Option for initial tracer location

      integer qt0num   ! Number of buble shaped initial tracer

      real qt0rx       ! Tracer radius or length in x direction
      real qt0ry       ! Tracer radius or length in y direction
      real qt0rz       ! Tracer radius or length in z direction

      real qt0cx       ! Center or origin in x coordinates of tracer
      real qt0cy       ! Center or origin in y coordinates of tracer
      real qt0cz       ! Center or origin in z coordinates of tracer

      real qt0ds       ! Distance between each tracer buble

      real qtdt        ! Emitted tracer intensity

      real cc05        ! 0.5 x cc

      real rxiv        ! 1.0 / qt0rx
      real ryiv        ! 1.0 / qt0ry
      real rziv        ! 1.0 / qt0rz

      real, intent(inout) :: xs(0:ni+1)
                       ! x coordinates at scalar points

      real, intent(inout) :: ys(0:nj+1)
                       ! y coordinates at scalar points

      real ctr(1:256)  ! Temporary variable

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer iqt      ! Index of number of buble

      real str         ! Temporary variable

      real a           ! Temporary variable
      real b           ! Temporary variable
      real c           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variables.

      call getiname(fpqt0opt,qt0opt)
      call getiname(fpqt0num,qt0num)
      call getrname(fpqt0rx,qt0rx)
      call getrname(fpqt0ry,qt0ry)
      call getrname(fpqt0rz,qt0rz)
      call getrname(fpqt0cx,qt0cx)
      call getrname(fpqt0cy,qt0cy)
      call getrname(fpqt0cz,qt0cz)
      call getrname(fpqt0ds,qt0ds)
      call getrname(fpqtdt,qtdt)

! -----

!! Emit the tracer mixing ratio from user specified location.

      if(qt0opt.eq.1.or.qt0opt.eq.2) then

! Get the x and the y coordinates at the scalar points.

        call getxy(iddx,iddy,'xx',0,ni+1,0,nj+1,xs,ys)

! -----

! Set the common used variables.

        cc05=.5e0*cc

        rxiv=1.e0/qt0rx
        ryiv=1.e0/qt0ry
        rziv=1.e0/qt0rz

! -----

! Emit the buble shaped tracer.

!$omp parallel default(shared) private(k,iqt)

        if(qt0opt.eq.1) then

!$omp do schedule(runtime)

          do iqt=1,qt0num
            ctr(iqt)=qt0cx+(real(iqt-1)+.5e0*real(1-qt0num))*qt0ds
          end do

!$omp end do

          do iqt=1,qt0num

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,str,a,b,c)

              do j=1,nj-1
              do i=1,ni-1
                a=rxiv*(xs(i)-ctr(iqt))
                b=ryiv*(ys(j)-qt0cy)
                c=rziv*(.5e0*(zph(i,j,k)+zph(i,j,k+1))-qt0cz)

                str=sqrt((a*a+b*b)+c*c)

                if(str.lt.1.e0) then

                  a=cos(cc05*str)

                  qtfrc(i,j,k)=qtfrc(i,j,k)+qtdt*a*a

                end if

              end do
              end do

!$omp end do

            end do

          end do

        else if(qt0opt.eq.2) then

!$omp do schedule(runtime)

          do iqt=1,qt0num
            ctr(iqt)=qt0cy+(real(iqt-1)+.5e0*real(1-qt0num))*qt0ds
          end do

!$omp end do

          do iqt=1,qt0num

            do k=1,nk-1

!$omp do schedule(runtime) private(i,j,str,a,b,c)

              do j=1,nj-1
              do i=1,ni-1
                a=rxiv*(xs(i)-qt0cx)
                b=ryiv*(ys(j)-ctr(iqt))
                c=rziv*(.5e0*(zph(i,j,k)+zph(i,j,k+1))-qt0cz)

                str=sqrt((a*a+b*b)+c*c)

                if(str.lt.1.e0) then

                  a=cos(cc05*str)

                  qtfrc(i,j,k)=qtfrc(i,j,k)+qtdt*a*a

                end if

              end do
              end do

!$omp end do

            end do

          end do

        end if

!$omp end parallel

! -----

      end if

!! -----

      end subroutine s_emitqt

!-----7--------------------------------------------------------------7--

      end module m_emitqt
