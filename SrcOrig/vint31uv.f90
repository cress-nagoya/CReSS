!***********************************************************************
      module m_vint31uv
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2011/12/17
!     Modification: 2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     interpolate the 3 dimensional input variable to the 1 dimensional
!     flat plane for x and y components of velocity.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_getcname
      use m_getiname
      use m_inichar

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: vint31uv, s_vint31uv

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface vint31uv

        module procedure s_vint31uv

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
      subroutine s_vint31uv(fpetrvar_gpv,fprefsfc_gpv,ape,              &
     &                      nid,njd,nkd,zlow,zdat,vardat,nk,z,varef)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpetrvar_gpv
                       ! Formal parameter of unique index of etrvar_gpv

      integer, intent(in) :: fprefsfc_gpv
                       ! Formal parameter of unique index of refsfc_gpv

      integer, intent(in) :: ape
                       ! Pointer of etrvar_gpv

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zlow(1:nid,1:njd)
                       ! Lowest z physical coordinates in data

      real, intent(in) :: zdat(1:nid,1:njd,1:nkd)
                       ! z physical coordinates in data

      real, intent(in) :: vardat(1:nid,1:njd,1:nkd)
                       ! Optional variable in data

      real, intent(in) :: z(1:nk)
                       ! zeta coordinates

! Output variable

      real, intent(out) :: varef(1:nid,1:njd,1:nk)
                       ! Optional vertically interpolated variable

! Internal shared variables

      character(len=108) etrvar_gpv
                       ! Control flag of extrapolating method

      integer refsfc_gpv
                       ! Option for
                       ! surface data reference in interpolating

      integer nkdm1    ! nkd - 1

! Internal private variables

      integer id       ! Data array index in x direction
      integer jd       ! Data array index in y direction
      integer kd       ! Data array index in z direction

      integer k        ! Array index in z direction

      real dk          ! Distance in z direction
                       ! between flat plane and data points

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(etrvar_gpv)

! -----

! Get the required namelist variables.

      call getcname(fpetrvar_gpv,etrvar_gpv)
      call getiname(fprefsfc_gpv,refsfc_gpv)

! -----

! Set the common used variable.

      nkdm1=nkd-1

! -----

!!! Interpolate the variable to the flat plane vertically.

!$omp parallel default(shared) private(k,kd)

! Extrapolate the variable.

      do k=1,nk

!$omp do schedule(runtime) private(id,jd,dk)

        do jd=1,njd
        do id=1,nid

          if(zdat(id,jd,nkd).le.z(k)) then

            dk=(z(k)-zdat(id,jd,nkdm1))                                 &
     &        /(zdat(id,jd,nkd)-zdat(id,jd,nkdm1))

            varef(id,jd,k)                                              &
     &        =(1.e0-dk)*vardat(id,jd,nkdm1)+dk*vardat(id,jd,nkd)

          end if

        end do
        end do

!$omp end do

      end do

! -----

!! Interpolate the variable in the case the surface data is not
!! referenced.

      if(refsfc_gpv.eq.0) then

! Extrapolate the variable.

        if(etrvar_gpv(ape:ape).eq.'o') then

          do k=1,nk

!$omp do schedule(runtime) private(id,jd,dk)

            do jd=1,njd
            do id=1,nid

              if(zdat(id,jd,1).gt.z(k)) then

                dk=(z(k)-zdat(id,jd,1))/(zdat(id,jd,2)-zdat(id,jd,1))

                varef(id,jd,k)                                          &
     &            =(1.e0-dk)*vardat(id,jd,1)+dk*vardat(id,jd,2)

              end if

            end do
            end do

!$omp end do

          end do

        else

          do k=1,nk

!$omp do schedule(runtime) private(id,jd)

            do jd=1,njd
            do id=1,nid

              if(zdat(id,jd,1).gt.z(k)) then

                varef(id,jd,k)=vardat(id,jd,1)

              end if

            end do
            end do

!$omp end do

          end do

        end if

! -----

! Interpolate the variable.

        do kd=1,nkd-1

          do k=1,nk

!$omp do schedule(runtime) private(id,jd,dk)

            do jd=1,njd
            do id=1,nid

              if(zdat(id,jd,kd).le.z(k)                                 &
     &          .and.zdat(id,jd,kd+1).gt.z(k)) then

                dk=(z(k)-zdat(id,jd,kd))                                &
     &            /(zdat(id,jd,kd+1)-zdat(id,jd,kd))

                varef(id,jd,k)                                          &
     &            =(1.e0-dk)*vardat(id,jd,kd)+dk*vardat(id,jd,kd+1)

              end if

            end do
            end do

!$omp end do

          end do

        end do

! -----

!! -----

!! Interpolate the variable in the case the surface data is referenced.

      else if(refsfc_gpv.eq.1) then

! Extrapolate the variable.

        if(etrvar_gpv(ape:ape).eq.'o') then

          do kd=1,nkd-1

            do k=1,nk

!$omp do schedule(runtime) private(id,jd,dk)

              do jd=1,njd
              do id=1,nid

                if(zlow(id,jd).gt.z(k)) then

                  if(zdat(id,jd,kd).le.zlow(id,jd)                      &
     &              .and.zdat(id,jd,kd+1).gt.zlow(id,jd)) then

                    dk=(z(k)-zlow(id,jd))/(zdat(id,jd,kd+1)-zlow(id,jd))

                    varef(id,jd,k)                                      &
     &                =(1.e0-dk)*vardat(id,jd,1)+dk*vardat(id,jd,kd+1)

                  end if

                end if

              end do
              end do

!$omp end do

            end do

          end do

        else

          do kd=1,nkd-1

            do k=1,nk

!$omp do schedule(runtime) private(id,jd)

              do jd=1,njd
              do id=1,nid

                if(zlow(id,jd).gt.z(k)) then

                  if(zdat(id,jd,kd).le.zlow(id,jd)                      &
     &              .and.zdat(id,jd,kd+1).gt.zlow(id,jd)) then

                    varef(id,jd,k)=vardat(id,jd,kd+1)

                  end if

                end if

              end do
              end do

!$omp end do

            end do

          end do

        end if

! -----

! Interpolate the variable.

        do kd=1,nkd-1

          do k=1,nk

!$omp do schedule(runtime) private(id,jd,dk)

            do jd=1,njd
            do id=1,nid

              if(zlow(id,jd).le.z(k)) then

                if(zdat(id,jd,kd).le.z(k)                               &
     &            .and.zdat(id,jd,kd+1).gt.z(k)) then

                  if(zdat(id,jd,kd).le.zlow(id,jd)) then

                    varef(id,jd,k)=vardat(id,jd,kd+1)

                  else

                    if(kd.eq.1) then

                      varef(id,jd,k)=vardat(id,jd,2)

                    else

                      dk=(z(k)-zdat(id,jd,kd))                          &
     &                  /(zdat(id,jd,kd+1)-zdat(id,jd,kd))

                      varef(id,jd,k)=(1.e0-dk)*vardat(id,jd,kd)         &
     &                  +dk*vardat(id,jd,kd+1)

                    end if

                  end if

                end if

              end if

            end do
            end do

!$omp end do

          end do

        end do

! -----

      end if

!! -----

!$omp end parallel

!!! -----

      end subroutine s_vint31uv

!-----7--------------------------------------------------------------7--

      end module m_vint31uv
