!***********************************************************************
      module m_getqt0
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2004/05/31
!     Modification: 2004/08/01, 2004/09/10, 2005/02/10, 2006/12/04,
!                   2007/01/05, 2007/01/20, 2007/01/31, 2007/10/19,
!                   2008/01/11, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2011/05/16, 2013/01/28, 2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the initial tracer.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_bcycle
      use m_combuf
      use m_comindx
      use m_commath
      use m_commpi
      use m_getbufgx
      use m_getbufgy
      use m_getbufsx
      use m_getbufsy
      use m_getiname
      use m_getrname
      use m_getxy
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

      public :: getqt0, s_getqt0

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getqt0

        module procedure s_getqt0

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic cos
      intrinsic sin
      intrinsic mod
      intrinsic real
      intrinsic sqrt

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_getqt0(fpqt0opt,fpqt0num,fpqt0,fpqt0rx,fpqt0ry,      &
     &                    fpqt0rz,fpqt0cx,fpqt0cy,fpqt0cz,fpqt0ds,      &
     &                    ni,nj,nk,zph,qt,xs,ys)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpqt0opt
                       ! Formal parameter of unique index of qt0opt

      integer, intent(in) :: fpqt0num
                       ! Formal parameter of unique index of qt0num

      integer, intent(in) :: fpqt0
                       ! Formal parameter of unique index of qt0

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

      integer, intent(in) :: ni
                       ! Model demension in x direction

      integer, intent(in) :: nj
                       ! Model demension in y direction

      integer, intent(in) :: nk
                       ! Model demension in z direction

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

! Output variable

      real, intent(out) :: qt(0:ni+1,0:nj+1,1:nk)
                       ! Tracer mixing ratio

! Internal shared variables

      integer qt0opt   ! Option for initial tracer location

      integer qt0num   ! Number of buble shaped initial tracer

      integer nx       ! Composite model dimension in x direction
      integer ny       ! Composite model dimension in y direction

      integer icmin    ! Minimum array index in x direction
                       ! in composite model dimension

      integer icmax    ! Maximum array index in x direction
                       ! in composite model dimension

      integer jcmin    ! Minimum array index in y direction
                       ! in composite model dimension

      integer jcmax    ! Maximum array index in y direction
                       ! in composite model dimension

      integer ic       ! Array index in x direction
                       ! in composite model dimension

      integer jc       ! Array index in y direction
                       ! in composite model dimension

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      integer ix       ! Parameter to make random noize

      real qt0         ! Magnitude or amplitude of tracer

      real qt0rx       ! Tracer radius or length in x direction
      real qt0ry       ! Tracer radius or length in y direction
      real qt0rz       ! Tracer radius or length in z direction

      real qt0cx       ! Center or origin in x coordinates of tracer
      real qt0cy       ! Center or origin in y coordinates of tracer
      real qt0cz       ! Center or origin in z coordinates of tracer

      real qt0ds       ! Distance between each tracer buble

      real cc05        ! 0.5 x cc

      real qt05        ! 0.5 x qt0

      real rxiv        ! 1.0 / qt0rx
      real ryiv        ! 1.0 / qt0ry
      real rziv        ! 1.0 / qt0rz

      real ccrxiv      ! cc / qt0rx
      real ccryiv      ! cc / qt0ry
      real ccrziv      ! cc / qt0rz

      real qt0zl       ! Lowest height of located tracer
      real qt0zh       ! Highest height of located tracer

      real zph8s       ! z physical coordinates at scalar points

      real, intent(inout) :: xs(0:ni+1)
                       ! x coordinates at scalar points

      real, intent(inout) :: ys(0:nj+1)
                       ! y coordinates at scalar points

      real ctr(1:256)  ! Temporary variable

! Internal private variables

      integer iqt      ! Index of number of buble

      integer i_sub    ! Substitute for i
      integer j_sub    ! Substitute for j
      integer k_sub    ! Substitute for k

      real zph8s_sub   ! Substitute for zph8s

      real str         ! Temporary variable

      real a           ! Temporary variable
      real b           ! Temporary variable
      real c           ! Temporary variable

!-----7--------------------------------------------------------------7--

! Get the required namelist variable.

      call getiname(fpqt0opt,qt0opt)

! -----

!! Set the buble shaped initial tracer to the array qt.

      if(qt0opt.eq.1.or.qt0opt.eq.2) then

! Get the required namelist variables.

        call getiname(fpqt0num,qt0num)
        call getrname(fpqt0,qt0)
        call getrname(fpqt0rx,qt0rx)
        call getrname(fpqt0ry,qt0ry)
        call getrname(fpqt0rz,qt0rz)
        call getrname(fpqt0cx,qt0cx)
        call getrname(fpqt0cy,qt0cy)
        call getrname(fpqt0cz,qt0cz)
        call getrname(fpqt0ds,qt0ds)

! -----

! Get the x and the y coordinates at the scalar points.

        call getxy(iddx,iddy,'xx',0,ni+1,0,nj+1,xs,ys)

! -----

! Set the common used variables.

        cc05=.5e0*cc

        rxiv=1.e0/qt0rx
        ryiv=1.e0/qt0ry
        rziv=1.e0/qt0rz

! -----

! Get the buble shaped initial tracer to the array qt.

!$omp parallel default(shared) private(k_sub,iqt)

        if(qt0opt.eq.1) then

!$omp do schedule(runtime)

          do iqt=1,qt0num
            ctr(iqt)=qt0cx+(real(iqt-1)+.5e0*real(1-qt0num))*qt0ds
          end do

!$omp end do

          do iqt=1,qt0num

            do k_sub=1,nk-1

!$omp do schedule(runtime) private(i_sub,j_sub,str,a,b,c)

              do j_sub=1,nj-1
              do i_sub=1,ni-1
                a=rxiv*(xs(i_sub)-ctr(iqt))
                b=ryiv*(ys(j_sub)-qt0cy)

                c=rziv*(.5e0*(zph(i_sub,j_sub,k_sub)                    &
     &            +zph(i_sub,j_sub,k_sub+1))-qt0cz)

                str=sqrt((a*a+b*b)+c*c)

                if(str.lt.1.e0) then

                  a=cos(cc05*str)

                  qt(i_sub,j_sub,k_sub)=qt0*a*a

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

            do k_sub=1,nk-1

!$omp do schedule(runtime) private(i_sub,j_sub,str,a,b,c)

              do j_sub=1,nj-1
              do i_sub=1,ni-1
                a=rxiv*(xs(i_sub)-qt0cx)
                b=ryiv*(ys(j_sub)-ctr(iqt))

                c=rziv*(.5e0*(zph(i_sub,j_sub,k_sub)                    &
     &            +zph(i_sub,j_sub,k_sub+1))-qt0cz)

                str=sqrt((a*a+b*b)+c*c)

                if(str.lt.1.e0) then

                  a=cos(cc05*str)

                  qt(i_sub,j_sub,k_sub)=qt0*a*a

                end if

              end do
              end do

!$omp end do

            end do

          end do

        end if

!$omp end parallel

! -----

!! -----

!! Set the sine curved initial tracer to the array qt.

      else if(qt0opt.eq.3.or.qt0opt.eq.4) then

! Get the required namelist variables.

        call getrname(fpqt0,qt0)
        call getrname(fpqt0rx,qt0rx)
        call getrname(fpqt0ry,qt0ry)
        call getrname(fpqt0rz,qt0rz)
        call getrname(fpqt0cx,qt0cx)
        call getrname(fpqt0cy,qt0cy)
        call getrname(fpqt0cz,qt0cz)

! -----

! Get the x and the y coordinates at the scalar points.

        call getxy(iddx,iddy,'xx',0,ni+1,0,nj+1,xs,ys)

! -----

! Set the common used variables.

        qt05=.5e0*qt0

        ccrxiv=cc/qt0rx
        ccryiv=cc/qt0ry
        ccrziv=cc/qt0rz

        qt0zl=qt0cz-qt0rz
        qt0zh=qt0zl+2.e0*qt0rz

! -----

! Get the sine curved initial tracer to the array qt.

!$omp parallel default(shared) private(k_sub)

        if(qt0opt.eq.3) then

          do k_sub=1,nk-1

!$omp do schedule(runtime) private(i_sub,j_sub,zph8s_sub)

            do j_sub=1,nj-1
            do i_sub=1,ni-1

              zph8s_sub=.5e0                                            &
     &          *(zph(i_sub,j_sub,k_sub)+zph(i_sub,j_sub,k_sub+1))

              if(zph8s_sub.ge.qt0zl.and.zph8s_sub.le.qt0zh) then

                qt(i_sub,j_sub,k_sub)=qt05*sin(ccrxiv*(xs(i_sub)-qt0cx))&
     &            *(1.e0+cos(ccrziv*(zph8s_sub-qt0cz)))

              end if

            end do
            end do

!$omp end do

          end do

        else if(qt0opt.eq.4) then

          do k_sub=1,nk-1

!$omp do schedule(runtime) private(i_sub,j_sub,zph8s_sub)

            do j_sub=1,nj-1
            do i_sub=1,ni-1

              zph8s_sub=.5e0                                            &
     &          *(zph(i_sub,j_sub,k_sub)+zph(i_sub,j_sub,k_sub+1))

              if(zph8s_sub.ge.qt0zl.and.zph8s_sub.le.qt0zh) then

                qt(i_sub,j_sub,k_sub)=qt05*sin(ccryiv*(ys(j_sub)-qt0cy))&
     &            *(1.e0+cos(ccrziv*(zph8s_sub-qt0cz)))

              end if

            end do
            end do

!$omp end do

          end do

        end if

!$omp end parallel

! -----

!! -----

!! Set the random initial tracer to the array qt.

      else if(qt0opt.eq.5) then

! Get the required namelist variables.

        call getrname(fpqt0,qt0)
        call getrname(fpqt0rz,qt0rz)
        call getrname(fpqt0cz,qt0cz)

! -----

! Set the common used variables.

        nx=(ni-3)*nigrp*nisub+3
        ny=(nj-3)*njgrp*njsub+3

        icmin=(ni-3)*nisub*igrp+(ni-3)*isub
        icmax=(ni-3)*nisub*igrp+(ni-3)*isub+ni
        jcmin=(nj-3)*njsub*jgrp+(nj-3)*jsub
        jcmax=(nj-3)*njsub*jgrp+(nj-3)*jsub+nj

        ix=0

        ccrziv=cc/qt0rz

        qt0zl=qt0cz-qt0rz
        qt0zh=qt0zl+2.e0*qt0rz

! -----

! Get the random initial tracer to the array qt.

        do k=1,nk-1
        do jc=1,ny-1
        do ic=1,nx-1

          ix=ix+1
          ix=mod(97*mod(ix,65536),65536)

          if((ic.gt.icmin.and.ic.lt.icmax)                              &
     &      .and.(jc.gt.jcmin.and.jc.lt.jcmax)) then

            i=ic-icmin
            j=jc-jcmin

            zph8s=.5e0*(zph(i,j,k)+zph(i,j,k+1))

            if(zph8s.ge.qt0zl.and.zph8s.le.qt0zh) then

              qt(i,j,k)=qt0*(i65536*real(ix)-.5e0)                      &
     &          *(1.e0+cos(ccrziv*(zph8s-qt0cz)))

            end if

          end if

        end do
        end do
        end do

! -----

      end if

!! -----

! Exchange the value horizontally.

      call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qt,1,1,sbuf)

      call s_shiftsx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

      call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qt,1,1,rbuf)

      call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,qt,1,1,sbuf)

      call s_shiftsy(idsbc,idnbc,'bnd',ni,nk,1,sbuf,rbuf)

      call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,qt,1,1,rbuf)

      call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,qt,1,1,sbuf)

      call s_shiftgx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

      call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,qt,1,1,rbuf)

      call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,qt,1,1,sbuf)

      call s_shiftgy(idsbc,idnbc,'bnd',ni,nk,1,sbuf,rbuf)

      call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,qt,1,1,rbuf)

! -----

! Set the periodic boundary conditions.

      call bcycle(idwbc,idebc,idsbc,idnbc,                              &
     &            2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,qt)

! -----

      end subroutine s_getqt0

!-----7--------------------------------------------------------------7--

      end module m_getqt0
