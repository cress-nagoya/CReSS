!***********************************************************************
      module m_set1d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1998/12/28
!     Modification: 1999/01/20, 1999/06/14, 1999/06/21, 1999/08/23,
!                   1999/09/30, 1999/11/01, 1999/11/19, 2000/01/17,
!                   2000/02/14, 2000/07/05, 2000/12/18, 2001/01/15,
!                   2001/04/15, 2001/05/29, 2001/06/29, 2001/08/07,
!                   2001/10/18, 2001/11/20, 2001/12/11, 2002/04/02,
!                   2002/06/18, 2002/11/20, 2003/03/28, 2003/04/30,
!                   2003/05/19, 2004/01/09, 2004/03/22, 2004/04/15,
!                   2004/07/01, 2004/08/01, 2004/08/20, 2004/09/01,
!                   2004/09/10, 2005/02/10, 2006/01/10, 2006/02/13,
!                   2006/09/21, 2007/01/20, 2007/01/31, 2007/05/14,
!                   2007/06/27, 2007/10/19, 2008/05/02, 2008/07/01,
!                   2008/07/25, 2008/08/25, 2009/01/05, 2009/02/27,
!                   2009/11/13, 2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the horizontally averaged variables from the read out sounding
!     data.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_commath
      use m_commpi
      use m_comphy
      use m_copy1d
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_getrname
      use m_inichar
      use m_setcst1d
      use m_vextlow

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: set1d, s_set1d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface set1d

        module procedure s_set1d

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic abs
      intrinsic exp
      intrinsic log
      intrinsic max
      intrinsic min

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_set1d(fpsndtyp,fpzsnd0,fppsnd0,fpthresq,fmois,       &
     &                   nlev,z1d,u1d,v1d,pt1d,qv1d,p1d,t1d,ltmp1,ltmp2)
!***********************************************************************

! Input variables

      integer, intent(in) :: fpsndtyp
                       ! Formal parameter of unique index of sndtyp

      integer, intent(in) :: fpzsnd0
                       ! Formal parameter of unique index of zsnd0

      integer, intent(in) :: fppsnd0
                       ! Formal parameter of unique index of psnd0

      integer, intent(in) :: fpthresq
                       ! Formal parameter of unique index of thresq

      integer, intent(in) :: nlev
                       ! Horizontally averaged vertical dimension

! Input and output variables

      real, intent(inout) :: z1d(0:nlev)
                       ! Horizontally averaged z physical coordinates

      real, intent(inout) :: u1d(0:nlev)
                       ! Horizontally averaged x components of velocity

      real, intent(inout) :: v1d(0:nlev)
                       ! Horizontally averaged y components of velocity

      real, intent(inout) :: pt1d(0:nlev)
                       ! Horizontally averaged potential temrerature

      real, intent(inout) :: qv1d(0:nlev)
                       ! Horizontally averaged water vapor mixing ratio

! Output variables

      character(len=5), intent(out) :: fmois
                       ! Control flag of air moisture

      real, intent(out) :: p1d(0:nlev)
                       ! Horizontally averaged pressure

! Internal shared variables

      character(len=108) sndtyp
                       ! Control flag of sounding data type

      integer nlev25   ! nlev / 4

      integer stat     ! Runtime status

      real thresq      ! Minimum threshold value of mixing ratio

      real zsnd0       ! Reference height in sounding data
      real psnd0       ! Reference pressure in sounding data

      real gdvcp2      ! 2.0 x g / cp
      real gdvrd2      ! 2.0 x g / rd

      real rddvg       ! rd / g
      real rddvcp      ! rd / cp
      real cpdvrd      ! cp / rd

      real p0iv        ! 1.0 / p0

      real es001       ! 0.01 x es0

      real pmin        ! Pressure at lowest level

      real pimin       ! Exner function at lowest level

      real lpmin       ! Natural logarithm of pressure at lowest level

      real qvmax       ! Maximum water vapor mixing ratio

      real itcon       ! Control flag of continuation of iteration

      real qvveps      ! Minimum value of convergence of iteration

      real, intent(inout) :: t1d(1:nlev)
                       ! Horizontally averaged temperature

      real, intent(inout) :: ltmp1(1:nlev)
                       ! Temporary array

      real, intent(inout) :: ltmp2(1:nlev)
                       ! Temporary array

! Internal private variables

      integer kl       ! Array index in z direction

      real ea          ! Pertial vapor pressure

      real qvf         ! Iterated water vapor mixing ratio

      real dqv         ! Variations of water vapor mixing ratio

! Remark

!     pt1d: This variable is also temporary.

!-----7--------------------------------------------------------------7--

! Initialize the character variable.

      call inichar(sndtyp)

! -----

! Get the required namelist variables.

      call getcname(fpsndtyp,sndtyp)
      call getrname(fpzsnd0,zsnd0)
      call getrname(fppsnd0,psnd0)
      call getrname(fpthresq,thresq)

! -----

! Set the common used variables.

      nlev25=nlev/4

      gdvcp2=2.e0*g/cp
      gdvrd2=2.e0*g/rd

      rddvg=rd/g
      rddvcp=rd/cp
      cpdvrd=cp/rd

      p0iv=1.e0/p0

      es001=.01e0*es0

      qvveps=1.e-5

! -----

!! Set the control flag fmois.

! Initialize the processed variable, qvmax.

      qvmax=lim36n

! -----

! Check the air moisture.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(kl) reduction(max: qvmax)

      do kl=1,nlev
        qvmax=max(qv1d(kl),qvmax)
      end do

!$omp end do

!$omp end parallel

! -----

! Finally set the common control flag fmois.

      if(qvmax.gt.0.e0) then

        write(fmois(1:5),'(a5)') 'moist'

      else

        write(fmois(1:5),'(a5)') 'dry  '

      end if

! -----

!! -----

!!!! Set the sounding variables in the case the 1st column in the
!!!! sounding data is the z physical coordinates.

      if(sndtyp(1:1).eq.'z') then

! Set the common used variables.

        if(sndtyp(2:2).eq.'t') then

          ltmp1(1)=log(psnd0)

          call s_copy1d(1,nlev,pt1d(1),t1d)

        else if(sndtyp(2:2).eq.'p') then

          ltmp1(1)=exp(rddvcp*log(p0iv*psnd0))

        end if

        if(sndtyp(3:3).eq.'r') then

          call s_copy1d(1,nlev,qv1d(1),ltmp2)

          call s_setcst1d(1,nlev,0.e0,qv1d(1))

        end if

! -----

!!! Perform iteration to get the pressure, temperature or potential
!!! temperature and water vapor mixing ratio.

        iterate: do

! Initialize the processed variables.

          lpmin=lim36
          pimin=lim36

          itcon=0.e0

! -----

!! Calculate the natural logarithm of pressure and virtual potential
!! temperature in the case the 2nd column in the sounding data is the
!! temperature.

          if(sndtyp(2:2).eq.'t') then

! Calculate the natural logarithm of pressure and virtual potential
! temperature.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(kl)

            do kl=1,nlev
              pt1d(kl)=t1d(kl)*(1.e0+epsav*qv1d(kl))/(1.e0+qv1d(kl))
            end do

!$omp end do

!$omp single private(kl)

            do kl=2,nlev
              ltmp1(kl)=ltmp1(kl-1)                                     &
     &          -gdvrd2*(z1d(kl)-z1d(kl-1))/(pt1d(kl-1)+pt1d(kl))

              lpmin=min(ltmp1(kl),lpmin)

            end do

!$omp end single

!$omp end parallel

! -----

! If error occured, call the procedure destroy.

            if(lpmin.le.0.e0) then
              stat=1
            else
              stat=0
            end if

            call chkerr(stat)

            if(stat.lt.0) then

              if(mype.eq.-stat-1) then

                call destroy('set1d   ',5,'cont',7,'              ',14, &
     &                       101,stat)

              end if

              call cpondpe

              call destroy('set1d   ',5,'stop',1001,'              ',14,&
     &                     101,stat)

            end if

! -----

!! -----

!! Calculate the Exner function and virtual temperature in the case the
!! 2nd column in the sounding data is the potential temperature.

          else if(sndtyp(2:2).eq.'p') then

! Calculate the Exner function and virtual temperature.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(kl)

            do kl=1,nlev
              t1d(kl)=pt1d(kl)*(1.e0+epsav*qv1d(kl))/(1.e0+qv1d(kl))
            end do

!$omp end do

!$omp single private(kl)

            do kl=2,nlev
              ltmp1(kl)=ltmp1(kl-1)                                     &
     &          -gdvcp2*(z1d(kl)-z1d(kl-1))/(t1d(kl-1)+t1d(kl))

              pimin=min(ltmp1(kl),pimin)

            end do

!$omp end single

!$omp end parallel

! -----

! If error occured, call the procedure destroy.

            if(pimin.le.0.e0) then
              stat=1
            else
              stat=0
            end if

            call chkerr(stat)

            if(stat.lt.0) then

              if(mype.eq.-stat-1) then

                call destroy('set1d   ',5,'cont',7,'              ',14, &
     &                       101,stat)

              end if

              call cpondpe

              call destroy('set1d   ',5,'stop',1001,'              ',14,&
     &                     101,stat)

            end if

! -----

          end if

!! -----

! Calculate the pressure and temperature or potential temperature.

!$omp parallel default(shared)

          if(sndtyp(2:2).eq.'t') then

!$omp do schedule(runtime) private(kl)

            do kl=1,nlev
              p1d(kl)=exp(ltmp1(kl))
              pt1d(kl)=t1d(kl)*exp(rddvcp*log(p0/p1d(kl)))
            end do

!$omp end do

          else if(sndtyp(2:2).eq.'p') then

!$omp do schedule(runtime) private(kl)

            do kl=1,nlev
              p1d(kl)=p0*exp(cpdvrd*log(ltmp1(kl)))
              t1d(kl)=pt1d(kl)*ltmp1(kl)
            end do

!$omp end do

          end if

!$omp end parallel

! -----

! Calculate the variations of the water vapor mixing ratio and the
! iterated water vapor mixing ratio.

          if(sndtyp(3:3).eq.'r') then

!$omp parallel default(shared)

!$omp do schedule(runtime) private(kl,ea,qvf,dqv) reduction(max: itcon)

            do kl=1,nlev

              if(t1d(kl).gt.tlow) then

                ea=es001*ltmp2(kl)                                      &
     &            *exp(17.269e0*(t1d(kl)-t0)/(t1d(kl)-35.86e0))

                qvf=epsva*ea/(p1d(kl)-ea)

                dqv=abs((qvf+eps)/(qv1d(kl)+eps)-1.e0)

                qv1d(kl)=qvf

              else

                ea=es001*ltmp2(kl)                                      &
     &            *exp(21.875e0*(t1d(kl)-t0)/(t1d(kl)-7.66e0))

                qvf=epsva*ea/(p1d(kl)-ea)

                dqv=abs((qvf+eps)/(qv1d(kl)+eps)-1.e0)

                qv1d(kl)=qvf

              end if

              itcon=max(dqv,itcon)

            end do

!$omp end do

!$omp end parallel

          end if

! -----

! If the variations of the water vapor mixing ratio is greater than
! qvveps, perform the iteration again.

          if(itcon.lt.qvveps) then

            exit iterate

          end if

! -----

        end do iterate

!!! -----

!!!! -----

!!! Set the sounding variables in the case the 1st column in the
!!! sounding data is the pressure.

      else if(sndtyp(1:1).eq.'p') then

! Initialize the processed variables.

        pmin=lim36

        p1d(1)=z1d(1)
        z1d(1)=zsnd0

! -----

!! Set the pressure.

! Swap the variable and get the minimum value.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(kl)

        do kl=2,nlev25
          p1d(kl)=z1d(kl)
        end do

!$omp end do

!$omp do schedule(runtime) private(kl) reduction(min: pmin)

        do kl=1,nlev25
          pmin=min(p1d(kl),pmin)
        end do

!$omp end do

!$omp end parallel

! -----

! If error occured, call the procedure destroy.

        if(pmin.le.0.e0) then
          stat=1
        else
          stat=0
        end if

        call chkerr(stat)

        if(stat.lt.0) then

          if(mype.eq.-stat-1) then

            call destroy('set1d   ',5,'cont',7,'              ',14,101, &
     &                   stat)

          end if

          call cpondpe

          call destroy('set1d   ',5,'stop',1001,'              ',14,101,&
     &                   stat)

        end if

! -----

!! -----

!! Calcuate the potential temperature or temperature, water vapor mixing
!! ratio and z physical coordiantes.

!$omp parallel default(shared)

! Calculate the potential temperature in the case the 2nd column in the
! sounding data is the temperature.

        if(sndtyp(2:2).eq.'t') then

!$omp do schedule(runtime) private(kl)

          do kl=1,nlev25
            t1d(kl)=pt1d(kl)
            pt1d(kl)=t1d(kl)*exp(rddvcp*log(p0/p1d(kl)))
          end do

!$omp end do

! -----

! Calculate the temperature in the case the 2nd column in the sounding
! data is the potential temperature.

        else if(sndtyp(2:2).eq.'p') then

!$omp do schedule(runtime) private(kl)

          do kl=1,nlev25
            t1d(kl)=pt1d(kl)*exp(rddvcp*log(p0iv*p1d(kl)))
          end do

!$omp end do

        end if

! -----

! Calculate the water vapor mixing ratio.

        if(sndtyp(3:3).eq.'r') then

!$omp do schedule(runtime) private(kl,ea)

          do kl=1,nlev25

            if(t1d(kl).gt.tlow) then

              ea=es001*qv1d(kl)                                         &
     &          *exp(17.269e0*(t1d(kl)-t0)/(t1d(kl)-35.86e0))

              qv1d(kl)=epsva*ea/(p1d(kl)-ea)

            else

              ea=es001*qv1d(kl)                                         &
     &          *exp(21.875e0*(t1d(kl)-t0)/(t1d(kl)-7.66e0))

              qv1d(kl)=epsva*ea/(p1d(kl)-ea)

            end if

          end do

!$omp end do

        end if

! -----

! Calculate the z physical coordinates.

!$omp do schedule(runtime) private(kl)

        do kl=1,nlev25
          ltmp1(kl)=t1d(kl)*(1.e0+epsav*qv1d(kl))/(1.e0+qv1d(kl))
        end do

!$omp end do

!$omp single private(kl)

        do kl=2,nlev25
          z1d(kl)=z1d(kl-1)-rddvg*(ltmp1(kl-1)+ltmp1(kl))               &
     &      /(p1d(kl-1)+p1d(kl))*(p1d(kl)-p1d(kl-1))
        end do

!$omp end single

! -----

!$omp end parallel

!! -----

      end if

!!! -----

! Force the water vapor mixing ratio more than user specified threshold
! value.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(kl)

      do kl=1,nlev

        if(qv1d(kl).le.thresq) then
          qv1d(kl)=0.e0
        end if

      end do

!$omp end do

!$omp end parallel

! -----

! Perform extrapolation in the lowest layer.

      call vextlow('oooooo',z1d,u1d,v1d,p1d,pt1d,qv1d)

! -----

      end subroutine s_set1d

!-----7--------------------------------------------------------------7--

      end module m_set1d
