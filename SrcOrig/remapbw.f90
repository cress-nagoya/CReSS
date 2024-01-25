!***********************************************************************
      module m_remapbw
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2006/08/08
!     Modification: 2006/09/30, 2007/10/19, 2008/01/11, 2008/05/02,
!                   2008/08/25, 2009/02/27, 2011/08/18, 2013/01/28,
!                   2013/02/13

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     remap the shifted water mass and concentrations for original
!     distributed water bins.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_commath
      use m_comphy

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: remapbw, s_remapbw

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface remapbw

        module procedure s_remapbw

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
      subroutine s_remapbw(k,nstr,ni,nj,nk,nqw,nnw,nqws,nnws,bmw,rbmw,  &
     &                     dbmw,bmws,mws,nws,mwbin,nwbin,bmwsl,bmwsr,   &
     &                     nwtd,nw0)
!***********************************************************************

! Input variables

      integer, intent(in) :: k
                       ! Array index in z direction

      integer, intent(in) :: nstr
                       ! Start index of distributed bin

      integer, intent(in) :: ni
                       ! Model dimension in x direction

      integer, intent(in) :: nj
                       ! Model dimension in y direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      integer, intent(in) :: nqw
                       ! Number of categories of water hydrometeor

      integer, intent(in) :: nnw
                       ! Number of categories of water concentration

      integer, intent(in) :: nqws
                       ! Number of categories
                       ! of shifted water hydrometeor

      integer, intent(in) :: nnws
                       ! Number of categories
                       ! of shifted water concentrations

      real, intent(in) :: bmw(1:nqw+1,1:3)
                       ! Mass at water bin boundary [g]

      real, intent(in) :: rbmw(1:nqw,1:2)
                       ! Related parameters of bmw

      real, intent(in) :: dbmw(1:nqw)
                       ! Differential between adjacent water bins [g]

      real, intent(in) :: bmws(0:ni+1,0:nj+1,1:nqws+1)
                       ! Shifted mass at water bin boundary [g]

      real, intent(in) :: mws(0:ni+1,0:nj+1,1:nqws)
                       ! Shifted water mass [g/cm^3]

      real, intent(in) :: nws(0:ni+1,0:nj+1,1:nnws)
                       ! Shifted water concentrations [1/cm^3]

! Input and output variables

      real, intent(inout) :: mwbin(0:ni+1,0:nj+1,1:nk,1:nqw)
                       ! Total water mass [g/cm^3]

      real, intent(inout) :: nwbin(0:ni+1,0:nj+1,1:nk,1:nnw)
                       ! Water concentrations [1/cm^3]

! Internal shared variables

      integer nqwp1    ! nqw + 1

      real, intent(inout) :: bmwsl(0:ni+1,0:nj+1)
                       ! Shifted mass at left water bin boundary

      real, intent(inout) :: bmwsr(0:ni+1,0:nj+1)
                       ! Shifted mass at right water bin boundary

      real, intent(inout) :: nwtd(0:ni+1,0:nj+1)
                       ! Tendency between shifted water bins

      real, intent(inout) :: nw0(0:ni+1,0:nj+1)
                       ! Interception of shifted water concentrations

! Internal private variables

      integer i        ! Array index in x direction
      integer j        ! Array index in y direction

      integer n        ! Array index in bin categories
      integer ns       ! Array index in shifted bin categories

      real dbmwiv      ! Inverse of differential of shifted water bin

      real cmw         ! Shifted mass at
                       ! center coordinates between adjacent water bins

      real cnw         ! Shifted concentrations at
                       ! center coordinates between adjacent water bins

      real nleft       ! Shifted concentrations
                       ! at left water bin boundary

      real nright      ! Shifted concentrations
                       ! at right water bin boundary

      real trdbmw      ! Differential between adjacent water bins
                       ! for triangular distributions

      real db          ! Integrated range
                       ! for current processed water bin

      real b05         ! Water mass at center coordinates
                       ! between current processed water bins

      real mwbr        ! Mean water mass

!-----7--------------------------------------------------------------7--

! Set the common used variable.

      nqwp1=nqw+1

! -----

!!!! Remap the shifted water mass and concentrations for original
!!!! distributed water bins and adjust the mean water mass to be
!!!! between their boundaries.

!$omp parallel default(shared) private(n,ns)

!!! Remap the shifted water mass and concentrations for original
!!! distributed water bins.

      do ns=1,nqws

! Set the shifted bin parameters.

!$omp do schedule(runtime)                                              &
!$omp&   private(i,j,dbmwiv,cmw,cnw,nleft,nright,trdbmw)

        do j=1,nj-1
        do i=1,ni-1

          if(mws(i,j,ns).gt.0.e0) then

            dbmwiv=1.e0/(bmws(i,j,ns+1)-bmws(i,j,ns))

            cmw=.5e0*(bmws(i,j,ns)+bmws(i,j,ns+1))

            cnw=dbmwiv*nws(i,j,ns)

            nwtd(i,j)                                                   &
     &        =12.e0*(mws(i,j,ns)-cmw*nws(i,j,ns))*dbmwiv*dbmwiv*dbmwiv

            nw0(i,j)=cnw-cmw*nwtd(i,j)

            nleft=cnw+nwtd(i,j)*(bmws(i,j,ns)-cmw)
            nright=cnw+nwtd(i,j)*(bmws(i,j,ns+1)-cmw)

            bmwsl(i,j)=bmws(i,j,ns)
            bmwsr(i,j)=bmws(i,j,ns+1)

            if(nleft.lt.0.e0) then

              bmwsl(i,j)                                                &
     &          =3.e0*mws(i,j,ns)/nws(i,j,ns)-2.e0*bmws(i,j,ns+1)

              trdbmw=bmws(i,j,ns+1)-bmwsl(i,j)

              if(trdbmw.lt.dbmw0) then

                bmwsr(i,j)=bmwsl(i,j)

              else

                dbmwiv=1.e0/trdbmw

                nwtd(i,j)=2.e0*dbmwiv*dbmwiv*nws(i,j,ns)

                nw0(i,j)=-nwtd(i,j)*bmwsl(i,j)

              end if

            end if

            if(nright.lt.0.e0) then

              bmwsr(i,j)                                                &
     &          =3.e0*mws(i,j,ns)/nws(i,j,ns)-2.e0*bmws(i,j,ns)

              trdbmw=bmwsr(i,j)-bmws(i,j,ns)

              if(trdbmw.lt.dbmw0) then

                bmwsl(i,j)=bmwsr(i,j)

              else

                dbmwiv=1.e0/trdbmw

                nwtd(i,j)=-2.e0*dbmwiv*dbmwiv*nws(i,j,ns)

                nw0(i,j)=-nwtd(i,j)*bmwsr(i,j)

              end if

            end if

          end if

        end do
        end do

!$omp end do

! -----

!! Perform remapping.

! Distribute to smallest bin.

        if(nqws.ne.1) then

!$omp do schedule(runtime) private(i,j,db,b05)

          do j=1,nj-1
          do i=1,ni-1

            if(mws(i,j,ns).gt.0.e0) then

              db=bmwsr(i,j)-bmw(1,1)

              if(db.lt.dbmw(1)) then

                mwbin(i,j,k,1)=mwbin(i,j,k,1)+mws(i,j,ns)
                nwbin(i,j,k,1)=nwbin(i,j,k,1)+nws(i,j,ns)

              else

                db=bmw(2,1)-bmwsl(i,j)

                if(db.gt.0.e0) then

                  b05=.5e0*(bmwsl(i,j)+bmw(2,1))

                  mwbin(i,j,k,1)=mwbin(i,j,k,1)+db*(b05*nw0(i,j)        &
     &              +oned3*nwtd(i,j)*(4.e0*b05*b05-bmwsl(i,j)*bmw(2,1)))

                  nwbin(i,j,k,1)                                        &
     &              =nwbin(i,j,k,1)+db*(nw0(i,j)+b05*nwtd(i,j))

                end if

              end if

            end if

          end do
          end do

!$omp end do

        end if

! -----

! Distribute to largest bin.

!$omp do schedule(runtime) private(i,j,db,b05)

        do j=1,nj-1
        do i=1,ni-1

          if(mws(i,j,ns).gt.0.e0) then

            db=bmw(nqwp1,1)-bmwsl(i,j)

            if(db.lt.dbmw(nqw)) then

              mwbin(i,j,k,nqw)=mwbin(i,j,k,nqw)+mws(i,j,ns)
              nwbin(i,j,k,nqw)=nwbin(i,j,k,nqw)+nws(i,j,ns)

            else

              db=bmwsr(i,j)-bmw(nqw,1)

              if(db.gt.0.e0) then

                b05=.5e0*(bmw(nqw,1)+bmwsr(i,j))

                mwbin(i,j,k,nqw)=mwbin(i,j,k,nqw)+db*(b05*nw0(i,j)      &
     &            +oned3*nwtd(i,j)*(4.e0*b05*b05-bmw(nqw,1)*bmwsr(i,j)))

                nwbin(i,j,k,nqw)                                        &
     &            =nwbin(i,j,k,nqw)+db*(nw0(i,j)+b05*nwtd(i,j))

              end if

            end if

          end if

        end do
        end do

!$omp end do

! -----

! Distribute to the other bins.

        if(nstr.le.nqw-1) then

          do n=nstr,nqw-1

!$omp do schedule(runtime) private(i,j,db,b05)

            do j=1,nj-1
            do i=1,ni-1

              if(mws(i,j,ns).gt.0.e0) then

                db=bmw(n+1,1)-bmwsl(i,j)

                if(db.gt.0.e0) then

                  if(db.lt.dbmw(n)) then

                    if(bmw(n+1,1)-bmwsr(i,j).gt.0.e0) then

                      mwbin(i,j,k,n)=mwbin(i,j,k,n)+mws(i,j,ns)
                      nwbin(i,j,k,n)=nwbin(i,j,k,n)+nws(i,j,ns)

                    else

                      b05=.5e0*(bmwsl(i,j)+bmw(n+1,1))

                      mwbin(i,j,k,n)=mwbin(i,j,k,n)                     &
     &                  +db*(b05*nw0(i,j)+oned3*nwtd(i,j)               &
     &                  *(4.e0*b05*b05-bmwsl(i,j)*bmw(n+1,1)))

                      nwbin(i,j,k,n)                                    &
     &                  =nwbin(i,j,k,n)+db*(nw0(i,j)+b05*nwtd(i,j))

                    end if

                  else

                    db=bmwsr(i,j)-bmw(n,1)

                    if(db.gt.0.e0) then

                      if(db.lt.dbmw(n)) then

                        b05=.5e0*(bmw(n,1)+bmwsr(i,j))

                        mwbin(i,j,k,n)=mwbin(i,j,k,n)                   &
     &                    +db*(b05*nw0(i,j)+oned3*nwtd(i,j)             &
     &                    *(4.e0*b05*b05-bmw(n,1)*bmwsr(i,j)))

                        nwbin(i,j,k,n)                                  &
     &                    =nwbin(i,j,k,n)+db*(nw0(i,j)+b05*nwtd(i,j))

                      else

                        mwbin(i,j,k,n)=mwbin(i,j,k,n)+dbmw(n)           &
     &                    *(rbmw(n,1)*nw0(i,j)+rbmw(n,2)*nwtd(i,j))

                        nwbin(i,j,k,n)=nwbin(i,j,k,n)                   &
     &                    +dbmw(n)*(nw0(i,j)+rbmw(n,1)*nwtd(i,j))

                      end if

                    end if

                  end if

                end if

              end if

            end do
            end do

!$omp end do

          end do

        end if

! -----

!! -----

      end do

!!! -----

! Adjust the mean water mass to be between their boundaries.

      if(nqws.ne.1) then

!$omp do schedule(runtime) private(i,j,mwbr)

        do j=1,nj-1
        do i=1,ni-1

          if(mwbin(i,j,k,1).gt.0.e0.and.nwbin(i,j,k,1).gt.0.e0) then

            mwbr=mwbin(i,j,k,1)/nwbin(i,j,k,1)

            if(mwbr.lt.bmw(1,3)) then

              nwbin(i,j,k,1)=mwbin(i,j,k,1)/bmw(1,3)

            end if

            if(mwbr.gt.bmw(2,2)) then

              nwbin(i,j,k,1)=mwbin(i,j,k,1)/bmw(2,2)

            end if

          else

            mwbin(i,j,k,1)=0.e0
            nwbin(i,j,k,1)=0.e0

          end if

        end do
        end do

!$omp end do

      end if

      do n=nstr,nqw

!$omp do schedule(runtime) private(i,j,mwbr)

        do j=1,nj-1
        do i=1,ni-1

          if(mwbin(i,j,k,n).gt.0.e0.and.nwbin(i,j,k,n).gt.0.e0) then

            mwbr=mwbin(i,j,k,n)/nwbin(i,j,k,n)

            if(mwbr.lt.bmw(n,3)) then

              nwbin(i,j,k,n)=mwbin(i,j,k,n)/bmw(n,3)

            end if

            if(mwbr.gt.bmw(n+1,2)) then

              nwbin(i,j,k,n)=mwbin(i,j,k,n)/bmw(n+1,2)

            end if

          else

            mwbin(i,j,k,n)=0.e0
            nwbin(i,j,k,n)=0.e0

          end if

        end do
        end do

!$omp end do

      end do

! -----

!$omp end parallel

!!!! -----

      end subroutine s_remapbw

!-----7--------------------------------------------------------------7--

      end module m_remapbw
