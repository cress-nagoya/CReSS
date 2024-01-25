!***********************************************************************
      module m_rbcs0
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2002/08/15
!     Modification: 2002/10/31, 2003/03/21, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2003/06/27, 2003/11/05, 2003/11/28,
!                   2003/12/12, 2004/05/31, 2004/09/10, 2006/04/03,
!                   2006/09/21, 2006/12/04, 2007/01/05, 2007/01/31,
!                   2007/10/19, 2008/05/02, 2008/08/25, 2009/02/27,
!                   2009/11/13, 2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the radiative lateral boundary conditions for optional scalar
!     variable.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commpi
      use m_getcname
      use m_getiname
      use m_getrname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: rbcs0, s_rbcs0

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface rbcs0

        module procedure s_rbcs0

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic max

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_rbcs0(fplbcvar,fpwbc,fpebc,fpsbc,fpnbc,fpadvopt,     &
     &                   fplbnews,apl,dt,ni,nj,nk,s,sp,scpx,scpy,sf)
!***********************************************************************

! Input variables

      integer, intent(in) :: fplbcvar
                       ! Formal parameter of unique index of lbcvar

      integer, intent(in) :: fpwbc
                       ! Formal parameter of unique index of wbc

      integer, intent(in) :: fpebc
                       ! Formal parameter of unique index of ebc

      integer, intent(in) :: fpsbc
                       ! Formal parameter of unique index of sbc

      integer, intent(in) :: fpnbc
                       ! Formal parameter of unique index of nbc

      integer, intent(in) :: fpadvopt
                       ! Formal parameter of unique index of advopt

      integer, intent(in) :: fplbnews
                       ! Formal parameter of unique index of lbnews

      integer, intent(in) :: apl
                       ! Pointer of lbcvar

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: dt
                       ! Time steps interval

      real, intent(in) :: s(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at present

      real, intent(in) :: sp(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at past

      real, intent(in) :: scpx(1:nj,1:nk,1:2)
                       ! Phase speed of optional scalar variable
                       ! on west and east boundary

      real, intent(in) :: scpy(1:ni,1:nk,1:2)
                       ! Phase speed of optional scalar variable
                       ! on south and north boundary

! Input and output variable

      real, intent(inout) :: sf(0:ni+1,0:nj+1,1:nk)
                       ! Optional scalar variable at future

! Internal shared variables

      character(len=108) lbcvar
                       ! Control flag of
                       ! lateral boundary forced variables

      integer wbc      ! Option for west boundary conditions
      integer ebc      ! Option for east boundary conditions
      integer sbc      ! Option for south boundary conditions
      integer nbc      ! Option for north boundary conditions

      integer advopt   ! Option for advection scheme

      integer nim1     ! ni - 1
      integer nim2     ! ni - 2
      integer njm1     ! nj - 1
      integer njm2     ! nj - 2

      real lbnews      ! Boundary damping coefficient

      real dmpdt       ! lbnews x dt or 2.0 x lbnews x dt

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction
      integer k        ! Array index in z direction

      real gamma       ! Temporary variable

      real radwe       ! Temporary variable
      real radsn       ! Temporary variable

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(lbcvar)

! -----

! Get the required namelist variables.

      call getcname(fplbcvar,lbcvar)
      call getiname(fpwbc,wbc)
      call getiname(fpebc,ebc)
      call getiname(fpsbc,sbc)
      call getiname(fpnbc,nbc)
      call getiname(fpadvopt,advopt)
      call getrname(fplbnews,lbnews)

! -----

! Set the common used variables.

      nim1=ni-1
      nim2=ni-2
      njm1=nj-1
      njm2=nj-2

      if(advopt.le.3) then

        if(lbcvar(apl:apl).eq.'o') then
          dmpdt=2.e0*lbnews*dt
        else
          dmpdt=0.e0
        end if

      else

        if(lbcvar(apl:apl).eq.'o') then
          dmpdt=lbnews*dt
        else
          dmpdt=0.e0
        end if

      end if

! -----

!! Set the radiative lateral boundary conditions.

!$omp parallel default(shared)

! Set the boundary conditions at the four corners.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(advopt.le.3) then

         if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.sbc.ge.4) then

!$omp do schedule(runtime) private(k,radwe,radsn)

           do k=2,nk-2
             radwe=(s(2,1,k)-sp(1,1,k))*scpx(1,k,1)/(1.e0-scpx(1,k,1))
             radsn=(s(1,2,k)-sp(1,1,k))*scpy(1,k,1)/(1.e0-scpy(1,k,1))

             sf(1,1,k)                                                  &
     &         =max(sp(1,1,k)-2.e0*(radwe+radsn)-dmpdt*sp(1,1,k),0.e0)

           end do

!$omp end do

         end if

         if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.sbc.ge.4) then

!$omp do schedule(runtime) private(k,radwe,radsn)

           do k=2,nk-2
             radwe=(s(nim2,1,k)-sp(nim1,1,k))                           &
     &         *scpx(1,k,2)/(1.e0+scpx(1,k,2))

             radsn=(s(nim1,2,k)-sp(nim1,1,k))                           &
     &         *scpy(nim1,k,1)/(1.e0-scpy(nim1,k,1))

             sf(nim1,1,k)=max(sp(nim1,1,k)+2.e0*(radwe-radsn)           &
     &         -dmpdt*sp(nim1,1,k),0.e0)

           end do

!$omp end do

         end if

        else

         if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.sbc.ge.4) then

!$omp do schedule(runtime) private(k,radwe,radsn)

           do k=2,nk-2
             radwe=(sp(2,1,k)-sp(1,1,k))*scpx(1,k,1)
             radsn=(sp(1,2,k)-sp(1,1,k))*scpy(1,k,1)

             sf(1,1,k)=max(sp(1,1,k)-(radwe+radsn)-dmpdt*sp(1,1,k),0.e0)

           end do

!$omp end do

         end if

         if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.sbc.ge.4) then

!$omp do schedule(runtime) private(k,radwe,radsn)

           do k=2,nk-2
             radwe=(sp(nim2,1,k)-sp(nim1,1,k))*scpx(1,k,2)
             radsn=(sp(nim1,2,k)-sp(nim1,1,k))*scpy(nim1,k,1)

             sf(nim1,1,k)                                               &
     &         =max(sp(nim1,1,k)+(radwe-radsn)-dmpdt*sp(nim1,1,k),0.e0)

           end do

!$omp end do

         end if

        end if

      end if

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(advopt.le.3) then

         if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.nbc.ge.4) then

!$omp do schedule(runtime) private(k,radwe,radsn)

           do k=2,nk-2
             radwe=(s(2,njm1,k)-sp(1,njm1,k))                           &
     &         *scpx(njm1,k,1)/(1.e0-scpx(njm1,k,1))

             radsn=(s(1,njm2,k)-sp(1,njm1,k))                           &
     &         *scpy(1,k,2)/(1.e0+scpy(1,k,2))

             sf(1,njm1,k)=max(sp(1,njm1,k)-2.e0*(radwe-radsn)           &
     &         -dmpdt*sp(1,njm1,k),0.e0)

           end do

!$omp end do

         end if

         if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.nbc.ge.4) then

!$omp do schedule(runtime) private(k,radwe,radsn)

           do k=2,nk-2
             radwe=(s(nim2,njm1,k)-sp(nim1,njm1,k))                     &
     &         *scpx(njm1,k,2)/(1.e0+scpx(njm1,k,2))

             radsn=(s(nim1,njm2,k)-sp(nim1,njm1,k))                     &
     &         *scpy(nim1,k,2)/(1.e0+scpy(nim1,k,2))

             sf(nim1,njm1,k)=max(sp(nim1,njm1,k)+2.e0*(radwe+radsn)     &
     &         -dmpdt*sp(nim1,njm1,k),0.e0)

           end do

!$omp end do

         end if

        else

         if(ebw.eq.1.and.isub.eq.0.and.wbc.ge.4.and.nbc.ge.4) then

!$omp do schedule(runtime) private(k,radwe,radsn)

           do k=2,nk-2
             radwe=(sp(2,njm1,k)-sp(1,njm1,k))*scpx(njm1,k,1)
             radsn=(sp(1,njm2,k)-sp(1,njm1,k))*scpy(1,k,2)

             sf(1,njm1,k)                                               &
     &         =max(sp(1,njm1,k)-(radwe-radsn)-dmpdt*sp(1,njm1,k),0.e0)

           end do

!$omp end do

         end if

         if(ebe.eq.1.and.isub.eq.nisub-1.and.ebc.ge.4.and.nbc.ge.4) then

!$omp do schedule(runtime) private(k,radwe,radsn)

           do k=2,nk-2
             radwe=(sp(nim2,njm1,k)-sp(nim1,njm1,k))*scpx(njm1,k,2)
             radsn=(sp(nim1,njm2,k)-sp(nim1,njm1,k))*scpy(nim1,k,2)

             sf(nim1,njm1,k)=max(sp(nim1,njm1,k)+(radwe+radsn)          &
     &         -dmpdt*sp(nim1,njm1,k),0.e0)

           end do

!$omp end do

         end if

        end if

      end if

! -----

! Set the west boundary conditions.

      if(ebw.eq.1.and.isub.eq.0) then

        if(wbc.ge.4) then

          if(advopt.le.3) then

!$omp do schedule(runtime) private(j,k,gamma)

            do k=2,nk-2
            do j=2,nj-2
              gamma=2.e0*scpx(j,k,1)/(1.e0-scpx(j,k,1))

              sf(1,j,k)=max(sp(1,j,k)-gamma*(s(2,j,k)-sp(1,j,k))        &
     &          -dmpdt*sp(1,j,k),0.e0)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=2,nj-2
              sf(1,j,k)=max(sp(1,j,k)-scpx(j,k,1)*(sp(2,j,k)-sp(1,j,k)) &
     &          -dmpdt*sp(1,j,k),0.e0)
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Set the east boundary conditions.

      if(ebe.eq.1.and.isub.eq.nisub-1) then

        if(ebc.ge.4) then

          if(advopt.le.3) then

!$omp do schedule(runtime) private(j,k,gamma)

            do k=2,nk-2
            do j=2,nj-2
              gamma=2.e0*scpx(j,k,2)/(1.e0+scpx(j,k,2))

              sf(nim1,j,k)                                              &
     &          =max(sp(nim1,j,k)+gamma*(s(nim2,j,k)-sp(nim1,j,k))      &
     &          -dmpdt*sp(nim1,j,k),0.e0)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(j,k)

            do k=2,nk-2
            do j=2,nj-2
              sf(nim1,j,k)=max(sp(nim1,j,k)+scpx(j,k,2)                 &
     &          *(sp(nim2,j,k)-sp(nim1,j,k))-dmpdt*sp(nim1,j,k),0.e0)
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Set the south boundary conditions.

      if(ebs.eq.1.and.jsub.eq.0) then

        if(sbc.ge.4) then

          if(advopt.le.3) then

!$omp do schedule(runtime) private(i,k,gamma)

            do k=2,nk-2
            do i=2,ni-2
              gamma=2.e0*scpy(i,k,1)/(1.e0-scpy(i,k,1))

              sf(i,1,k)=max(sp(i,1,k)-gamma*(s(i,2,k)-sp(i,1,k))        &
     &          -dmpdt*sp(i,1,k),0.e0)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=2,ni-2
              sf(i,1,k)=max(sp(i,1,k)-scpy(i,k,1)*(sp(i,2,k)-sp(i,1,k)) &
     &          -dmpdt*sp(i,1,k),0.e0)
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

! Set the north boundary conditions.

      if(ebn.eq.1.and.jsub.eq.njsub-1) then

        if(nbc.ge.4) then

          if(advopt.le.3) then

!$omp do schedule(runtime) private(i,k,gamma)

            do k=2,nk-2
            do i=2,ni-2
              gamma=2.e0*scpy(i,k,2)/(1.e0+scpy(i,k,2))

              sf(i,njm1,k)                                              &
     &          =max(sp(i,njm1,k)+gamma*(s(i,njm2,k)-sp(i,njm1,k))      &
     &          -dmpdt*sp(i,njm1,k),0.e0)

            end do
            end do

!$omp end do

          else

!$omp do schedule(runtime) private(i,k)

            do k=2,nk-2
            do i=2,ni-2
              sf(i,njm1,k)=max(sp(i,njm1,k)+scpy(i,k,2)                 &
     &          *(sp(i,njm2,k)-sp(i,njm1,k))-dmpdt*sp(i,njm1,k),0.e0)
            end do
            end do

!$omp end do

          end if

        end if

      end if

! -----

!$omp end parallel

!! -----

      end subroutine s_rbcs0

!-----7--------------------------------------------------------------7--

      end module m_rbcs0
