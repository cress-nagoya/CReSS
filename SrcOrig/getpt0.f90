!***********************************************************************
      module m_getpt0
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/01/20
!     Modification: 1999/01/25, 1999/03/25, 1999/04/06, 1999/05/10,
!                   1999/06/07, 1999/07/05, 1999/07/23, 1999/07/26,
!                   1999/08/03, 1999/08/18, 1999/08/23, 1999/09/06,
!                   1999/09/30, 1999/10/07, 1999/10/12, 1999/11/01,
!                   1999/11/19, 1999/12/06, 2000/01/17, 2000/04/18,
!                   2000/07/05, 2000/12/18, 2001/03/13, 2001/05/29,
!                   2001/06/06, 2001/09/13, 2001/11/14, 2002/04/02,
!                   2002/06/06, 2002/06/18, 2002/08/15, 2003/04/30,
!                   2003/05/19, 2003/09/01, 2003/10/31, 2003/11/05,
!                   2004/05/07, 2004/05/31, 2004/08/01, 2004/09/10,
!                   2005/02/10, 2006/01/10, 2006/12/04, 2007/01/05,
!                   2007/01/20, 2007/01/31, 2007/10/19, 2008/01/11,
!                   2008/05/02, 2008/08/25, 2009/02/27, 2013/01/28,
!                   2012/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the initial potential temperature perturbation.

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
      use m_setcst3d
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

      public :: getpt0, s_getpt0

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface getpt0

        module procedure s_getpt0

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
      subroutine s_getpt0(fppt0opt,fppt0num,fpptp0,fppt0rx,fppt0ry,     &
     &                    fppt0rz,fppt0cx,fppt0cy,fppt0cz,fppt0ds,      &
     &                    ni,nj,nk,zph,ptp,xs,ys)
!***********************************************************************

! Input variables

      integer, intent(in) :: fppt0opt
                       ! Formal parameter of unique index of pt0opt

      integer, intent(in) :: fppt0num
                       ! Formal parameter of unique index of pt0num

      integer, intent(in) :: fpptp0
                       ! Formal parameter of unique index of ptp0

      integer, intent(in) :: fppt0rx
                       ! Formal parameter of unique index of pt0rx

      integer, intent(in) :: fppt0ry
                       ! Formal parameter of unique index of pt0ry

      integer, intent(in) :: fppt0rz
                       ! Formal parameter of unique index of pt0rz

      integer, intent(in) :: fppt0cx
                       ! Formal parameter of unique index of pt0cx

      integer, intent(in) :: fppt0cy
                       ! Formal parameter of unique index of pt0cy

      integer, intent(in) :: fppt0cz
                       ! Formal parameter of unique index of pt0cz

      integer, intent(in) :: fppt0ds
                       ! Formal parameter of unique index of pt0ds

      integer, intent(in) :: ni
                       ! Model demension in x direction

      integer, intent(in) :: nj
                       ! Model demension in y direction

      integer, intent(in) :: nk
                       ! Model demension in z direction

      real, intent(in) :: zph(0:ni+1,0:nj+1,1:nk)
                       ! z physical coordinates

! Output variable

      real, intent(out) :: ptp(0:ni+1,0:nj+1,1:nk)
                       ! Potential temperature perturbation

! Internal shared variables

      integer pt0opt   ! Option for initial potential temperature

      integer pt0num   ! Number of buble shaped
                       ! initial potential temperature perturbation

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

      real ptp0        ! Magnitude or amplitude
                       ! of potential temperature perturbation

      real pt0rx       ! Potential temperature perturbation radius
                       ! or length in x direction

      real pt0ry       ! Potential temperature perturbation radius
                       ! or length in y direction

      real pt0rz       ! Potential temperature perturbation radius
                       ! or length in z direction

      real pt0cx       ! Center or origin in x coordinates
                       ! of potential temperature perturbation

      real pt0cy       ! Center or origin in y coordinates
                       ! of potential temperature perturbation

      real pt0cz       ! Center or origin in z coordinates
                       ! of potential temperature perturbation

      real pt0ds       ! Distance between each perturbation buble

      real cc05        ! 0.5 x cc

      real ptp05       ! 0.5 x ptp0

      real rxiv        ! 1.0 / pt0rx
      real ryiv        ! 1.0 / pt0ry
      real rziv        ! 1.0 / pt0rz

      real ccrxiv      ! cc / pt0rx
      real ccryiv      ! cc / pt0ry
      real ccrziv      ! cc / pt0rz

      real pt0zl       ! Lowest height
                       ! of potential temperature perturbation

      real pt0zh       ! Highest height
                       ! of potential temperature perturbation

      real zph8s       ! z physical coordinates at scalar points

      real, intent(inout) :: xs(0:ni+1)
                       ! x coordinates at scalar points

      real, intent(inout) :: ys(0:nj+1)
                       ! y coordinates at scalar points

      real ctr(1:256)  ! Temporary variable

! Internal private variables

      integer ipt      ! Index of number of buble

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

      call getiname(fppt0opt,pt0opt)

! -----

! Fill in the array ptp with 0 in the case there is no initial potential
! temperature perturbation.

      if(pt0opt.eq.0) then

        call setcst3d(0,ni+1,0,nj+1,1,nk,0.e0,ptp)

! -----

!! Set the buble shaped initial potential temperature perturbation to
!! the array ptp.

      else if(pt0opt.eq.1.or.pt0opt.eq.2) then

! Get the required namelist variables.

        call getiname(fppt0num,pt0num)
        call getrname(fpptp0,ptp0)
        call getrname(fppt0rx,pt0rx)
        call getrname(fppt0ry,pt0ry)
        call getrname(fppt0rz,pt0rz)
        call getrname(fppt0cx,pt0cx)
        call getrname(fppt0cy,pt0cy)
        call getrname(fppt0cz,pt0cz)
        call getrname(fppt0ds,pt0ds)

! -----

! Get the x and the y coordinates at the scalar points.

        call getxy(iddx,iddy,'xx',0,ni+1,0,nj+1,xs,ys)

! -----

! Set the common used variables.

        cc05=.5e0*cc

        rxiv=1.e0/pt0rx
        ryiv=1.e0/pt0ry
        rziv=1.e0/pt0rz

! -----

! Get the buble shaped initial potential temperature perturbation to the
! array ptp.

!$omp parallel default(shared) private(k_sub,ipt)

        if(pt0opt.eq.1) then

!$omp do schedule(runtime)

          do ipt=1,pt0num
            ctr(ipt)=pt0cx+(real(ipt-1)+.5e0*real(1-pt0num))*pt0ds
          end do

!$omp end do

          do ipt=1,pt0num

            do k_sub=1,nk-1

!$omp do schedule(runtime) private(i_sub,j_sub,str,a,b,c)

              do j_sub=1,nj-1
              do i_sub=1,ni-1
                a=rxiv*(xs(i_sub)-ctr(ipt))
                b=ryiv*(ys(j_sub)-pt0cy)

                c=rziv*(.5e0*(zph(i_sub,j_sub,k_sub)                    &
     &            +zph(i_sub,j_sub,k_sub+1))-pt0cz)

                str=sqrt((a*a+b*b)+c*c)

                if(str.lt.1.e0) then

                  a=cos(cc05*str)

                  ptp(i_sub,j_sub,k_sub)=ptp0*a*a

                end if

              end do
              end do

!$omp end do

            end do

          end do

        else if(pt0opt.eq.2) then

!$omp do schedule(runtime)

          do ipt=1,pt0num
            ctr(ipt)=pt0cy+(real(ipt-1)+.5e0*real(1-pt0num))*pt0ds
          end do

!$omp end do

          do ipt=1,pt0num

            do k_sub=1,nk-1

!$omp do schedule(runtime) private(i_sub,j_sub,str,a,b,c)

              do j_sub=1,nj-1
              do i_sub=1,ni-1
                a=rxiv*(xs(i_sub)-pt0cx)
                b=ryiv*(ys(j_sub)-ctr(ipt))

                c=rziv*(.5e0*(zph(i_sub,j_sub,k_sub)                    &
     &            +zph(i_sub,j_sub,k_sub+1))-pt0cz)

                str=sqrt((a*a+b*b)+c*c)

                if(str.lt.1.e0) then

                  a=cos(cc05*str)

                  ptp(i_sub,j_sub,k_sub)=ptp0*a*a

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

!! Set the sine curved initial potential temperature perturbation to the
!! array ptp.

      else if(pt0opt.eq.3.or.pt0opt.eq.4) then

! Get the required namelist variables.

        call getrname(fpptp0,ptp0)
        call getrname(fppt0rx,pt0rx)
        call getrname(fppt0ry,pt0ry)
        call getrname(fppt0rz,pt0rz)
        call getrname(fppt0cx,pt0cx)
        call getrname(fppt0cy,pt0cy)
        call getrname(fppt0cz,pt0cz)

! -----

! Get the x and the y coordinates at the scalar points.

        call getxy(iddx,iddy,'xx',0,ni+1,0,nj+1,xs,ys)

! -----

! Set the common used variables.

        ptp05=.5e0*ptp0

        ccrxiv=cc/pt0rx
        ccryiv=cc/pt0ry
        ccrziv=cc/pt0rz

        pt0zl=pt0cz-pt0rz
        pt0zh=pt0zl+2.e0*pt0rz

! -----

! Get the sine curved initial potential temperature perturbation to the
! array ptp.

!$omp parallel default(shared) private(k_sub)

        if(pt0opt.eq.3) then

          do k_sub=1,nk-1

!$omp do schedule(runtime) private(i_sub,j_sub,zph8s_sub)

            do j_sub=1,nj-1
            do i_sub=1,ni-1

              zph8s_sub=.5e0                                            &
     &          *(zph(i_sub,j_sub,k_sub)+zph(i_sub,j_sub,k_sub+1))

              if(zph8s_sub.ge.pt0zl.and.zph8s_sub.le.pt0zh) then

                ptp(i_sub,j_sub,k_sub)                                  &
     &            =ptp05*sin(ccrxiv*(xs(i_sub)-pt0cx))                  &
     &            *(1.e0+cos(ccrziv*(zph8s_sub-pt0cz)))

              end if

            end do
            end do

!$omp end do

          end do

        else if(pt0opt.eq.4) then

          do k_sub=1,nk-1

!$omp do schedule(runtime) private(i_sub,j_sub,zph8s_sub)

            do j_sub=1,nj-1
            do i_sub=1,ni-1

              zph8s_sub=.5e0                                            &
     &          *(zph(i_sub,j_sub,k_sub)+zph(i_sub,j_sub,k_sub+1))

              if(zph8s_sub.ge.pt0zl.and.zph8s_sub.le.pt0zh) then

                ptp(i_sub,j_sub,k_sub)                                  &
     &            =ptp05*sin(ccryiv*(ys(j_sub)-pt0cy))                  &
     &            *(1.e0+cos(ccrziv*(zph8s_sub-pt0cz)))

              end if

            end do
            end do

!$omp end do

          end do

        end if

!$omp end parallel

! -----

!! -----

!! Set the random initial potential temperature perturbation to the
!! array ptp.

      else if(pt0opt.eq.5) then

! Get the required namelist variables.

        call getrname(fpptp0,ptp0)
        call getrname(fppt0rz,pt0rz)
        call getrname(fppt0cz,pt0cz)

! -----

! Set the common used variables.

        nx=(ni-3)*nigrp*nisub+3
        ny=(nj-3)*njgrp*njsub+3

        icmin=(ni-3)*nisub*igrp+(ni-3)*isub
        icmax=(ni-3)*nisub*igrp+(ni-3)*isub+ni
        jcmin=(nj-3)*njsub*jgrp+(nj-3)*jsub
        jcmax=(nj-3)*njsub*jgrp+(nj-3)*jsub+nj

        ix=0

        ccrziv=cc/pt0rz

        pt0zl=pt0cz-pt0rz
        pt0zh=pt0zl+2.e0*pt0rz

! -----

! Get the random initial potential temperature perturbation to the array
! ptp.

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

            if(zph8s.ge.pt0zl.and.zph8s.le.pt0zh) then

              ptp(i,j,k)=ptp0*(i65536*real(ix)-.5e0)                    &
     &          *(1.e0+cos(ccrziv*(zph8s-pt0cz)))

            end if

          end if

        end do
        end do
        end do

! -----

      end if

!! -----

! Exchange the value horizontally.

      call s_putbufsx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptp,1,1,sbuf)

      call s_shiftsx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

      call s_getbufsx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptp,1,1,rbuf)

      call s_putbufsy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ptp,1,1,sbuf)

      call s_shiftsy(idsbc,idnbc,'bnd',ni,nk,1,sbuf,rbuf)

      call s_getbufsy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ptp,1,1,rbuf)

      call s_putbufgx(idwbc,idebc,'bnd',2,ni-2,ni,nj,nk,ptp,1,1,sbuf)

      call s_shiftgx(idwbc,idebc,'bnd',nj,nk,1,sbuf,rbuf)

      call s_getbufgx(idwbc,idebc,'bnd',1,ni-1,ni,nj,nk,ptp,1,1,rbuf)

      call s_putbufgy(idsbc,idnbc,'bnd',2,nj-2,ni,nj,nk,ptp,1,1,sbuf)

      call s_shiftgy(idsbc,idnbc,'bnd',ni,nk,1,sbuf,rbuf)

      call s_getbufgy(idsbc,idnbc,'bnd',1,nj-1,ni,nj,nk,ptp,1,1,rbuf)

! -----

! Set the periodic boundary conditions.

      call bcycle(idwbc,idebc,idsbc,idnbc,                              &
     &            2,1,ni-2,ni-1,2,1,nj-2,nj-1,ni,nj,nk,ptp)

! -----

      end subroutine s_getpt0

!-----7--------------------------------------------------------------7--

      end module m_getpt0
