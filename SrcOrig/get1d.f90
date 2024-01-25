!***********************************************************************
      module m_get1d
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 1999/05/20
!     Modification: 1999/06/21, 1999/06/28, 1999/07/05, 1999/08/03,
!                   1999/08/23, 1999/09/30, 1999/10/12, 1999/11/01,
!                   2000/01/17, 2000/07/05, 2001/01/15, 2001/04/15,
!                   2001/05/29, 2001/06/29, 2001/07/13, 2002/01/21,
!                   2002/04/02, 2002/04/09, 2002/06/18, 2002/07/15,
!                   2002/08/15, 2002/09/09, 2003/04/30, 2003/05/19,
!                   2003/12/12, 2004/01/09, 2004/04/15, 2004/05/07,
!                   2004/07/01, 2004/08/01, 2004/08/20, 2004/09/01,
!                   2004/09/10, 2005/02/10, 2006/01/10, 2006/09/21,
!                   2007/01/20, 2007/01/31, 2007/03/23, 2007/05/14,
!                   2007/06/27, 2007/09/04, 2007/10/19, 2008/05/02,
!                   2008/08/25, 2008/12/11, 2009/01/05, 2009/02/27,
!                   2009/11/13, 2013/01/28, 2013/02/13, 2013/03/27

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     get the horizontally averaged variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_chkerr
      use m_comindx
      use m_commath
      use m_commpi
      use m_comphy
      use m_copy1d
      use m_cpondpe
      use m_destroy
      use m_getcname
      use m_inichar
      use m_setcst1d
      use m_vextlow
      use m_vextund
      use m_vint11

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: get1d, s_get1d

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface get1d

        module procedure s_get1d

      end interface

!-----7--------------------------------------------------------------7--

! Intrinsic procedure

      intrinsic exp
      intrinsic log
      intrinsic min
      intrinsic real

! External procedure

!     none

!-----7--------------------------------------------------------------7--

! Internal module procedure

      contains

!***********************************************************************
      subroutine s_get1d(fpgpvvar,fpetrvar_gpv,fproc,                   &
     &                   idstr,idend,jdstr,jdend,kref,nid,njd,nkd,      &
     &                   zdat,nk,z,ubdat,vbdat,pbdat,ptbdat,qvbdat,     &
     &                   z1d,u1d,v1d,p1d,pt1d,qv1d,tmp1)
!***********************************************************************

! Input variables

      character(len=6), intent(in) :: fproc
                       ! Control flag of processing type

      integer, intent(in) :: fpgpvvar
                       ! Formal parameter of unique index of gpvvar

      integer, intent(in) :: fpetrvar_gpv
                       ! Formal parameter of unique index of etrvar_gpv

      integer, intent(in) :: idstr
                       ! Minimum do loops index in x direction

      integer, intent(in) :: idend
                       ! Maximum do loops index in x direction

      integer, intent(in) :: jdstr
                       ! Minimum do loops index in y direction

      integer, intent(in) :: jdend
                       ! Maximum do loops index in y direction

      integer, intent(in) :: kref
                       ! Reference index

      integer, intent(in) :: nid
                       ! Data dimension in x direction

      integer, intent(in) :: njd
                       ! Data dimension in y direction

      integer, intent(in) :: nkd
                       ! Data dimension in z direction

      integer, intent(in) :: nk
                       ! Model dimension in z direction

      real, intent(in) :: zdat(1:nid,1:njd,1:nkd)
                       ! z physical coordinates in data

      real, intent(in) :: z(1:nk)
                       ! zeta coordinates

      real, intent(in) :: ubdat(1:nid,1:njd,1:nk)
                       ! Base state x components of velocity in data

      real, intent(in) :: vbdat(1:nid,1:njd,1:nk)
                       ! Base state y components of velocity in data

      real, intent(in) :: pbdat(1:nid,1:njd,1:nk)
                       ! Base state pressure in data

      real, intent(in) :: ptbdat(1:nid,1:njd,1:nk)
                       ! Base state potential temperature in data

      real, intent(in) :: qvbdat(1:nid,1:njd,1:nk)
                       ! Base state water vapor mixing ratio in data

! Output variables

      real, intent(out) :: z1d(0:4*nk-3)
                       ! Horizontally averaged z coordinates

      real, intent(out) :: u1d(0:4*nk-3)
                       ! Horizontally averaged x components of velocity

      real, intent(out) :: v1d(0:4*nk-3)
                       ! Horizontally averaged y components of velocity

      real, intent(out) :: p1d(0:4*nk-3)
                       ! Horizontally averaged pressure

      real, intent(out) :: pt1d(0:4*nk-3)
                       ! Horizontally averaged potential temperature

      real, intent(out) :: qv1d(0:4*nk-3)
                       ! Horizontally averaged water vapor mixing ratio

! Internal shared variables

      character(len=108) gpvvar
                       ! Control flag of input GPV data variables

      character(len=108) etrvar_gpv
                       ! Control flag of extrapolating method

      integer k        ! Array index in z direction

      integer kr4m3    ! 4 x kref - 3

      integer stat     ! Runtime status

      real nginv       ! Inverse of number of grid points

      real gdvcp2      ! 2.0 x g / cp

      real cpdvrd      ! cp / rd
      real rddvcp      ! rd / cp

      real p0iv        ! 1.0 / p0

      real dz1dl       ! 0.25 x (z(kref) - z(1)) / real(kref - 1)
      real dz1dh       ! 0.25 x (z(nk) - z(kref)) / real(nk - kref)

      real pref        ! Reference pressure

      real pimin       ! Lowest base state Exner function

      real, intent(inout) :: tmp1(1:nk)
                       ! Temporary array

! Internal private variables

      integer id       ! Data array index in x direction
      integer jd       ! Data array index in y direction

      integer k_sub    ! Substitute for k

      integer ngcnt    ! Number of grid points

      real ngcntv      ! 1.0 / real(ngcnt)

! Remark

!     p1d,pt1d: These variables are also temporary.

!-----7--------------------------------------------------------------7--

! Initialize the character variables.

      call inichar(gpvvar)
      call inichar(etrvar_gpv)

! -----

! Get the required namelist variables.

      call getcname(fpgpvvar,gpvvar)
      call getcname(fpetrvar_gpv,etrvar_gpv)

! -----

! Set the common used variables.

      kr4m3=4*kref-3

      if(idstr.le.idend) then
        nginv=1.e0/real((idend-idstr+1)*(jdend-jdstr+1))
      else
        nginv=1.e0/real((nid+idend-idstr+1)*(jdend-jdstr+1))
      end if

      gdvcp2=2.e0*g/cp

      cpdvrd=cp/rd
      rddvcp=rd/cp

      p0iv=1.e0/p0

! -----

! Fill in the 1 dimensional array with 0.

      pref=0.e0

      call s_setcst1d(1,nk,0.e0,u1d(1))
      call s_setcst1d(1,nk,0.e0,v1d(1))
      call s_setcst1d(1,nk,0.e0,pt1d(1))
      call s_setcst1d(1,nk,0.e0,qv1d(1))

! -----

!!!! Create the 1 dimensional variables.

!!! Create the 1 dimensional variables in vertical model dimension.

! Set the common used variables.

      dz1dl=.25e0*(z(kref)-z(1))/real(kref-1)
      dz1dh=.25e0*(z(nk)-z(kref))/real(nk-kref)

! -----

!! Get the 1 dimensional variables.

!$omp parallel default(shared)

! Calculate the interpolated zeta coordinates.

!$omp do schedule(runtime) private(k_sub)

      do k_sub=1,4*kref-3
        z1d(k_sub)=z(1)+real(k_sub-1)*dz1dl
      end do

!$omp end do

!$omp do schedule(runtime) private(k_sub)

      do k_sub=4*kref-2,4*nk-3
        z1d(k_sub)=z(kref)+real(k_sub-kr4m3)*dz1dh
      end do

!$omp end do

! -----

! Get the reference sumed pressure.

      if(idstr.le.idend) then

!$omp do schedule(runtime) private(id,jd) reduction(+: pref)

        do jd=jdstr,jdend
        do id=idstr,idend
          pref=pref+pbdat(id,jd,kref)
        end do
        end do

!$omp end do

      else

!$omp do schedule(runtime) private(id,jd) reduction(+: pref)

        do jd=jdstr,jdend
        do id=idstr,nid
          pref=pref+pbdat(id,jd,kref)
        end do
        end do

!$omp end do

!$omp do schedule(runtime) private(id,jd) reduction(+: pref)

        do jd=jdstr,jdend
        do id=1,idend
          pref=pref+pbdat(id,jd,kref)
        end do
        end do

!$omp end do

      end if

! -----

! Get the averaged x and y components of velocity and potential
! temperature.

      if(idstr.le.idend) then

        if(fproc(1:3).eq.'all') then

!$omp do schedule(runtime) private(id,jd,k_sub)

          do k_sub=1,nk

            do jd=jdstr,jdend
            do id=idstr,idend
              u1d(k_sub)=u1d(k_sub)+ubdat(id,jd,k_sub)
              v1d(k_sub)=v1d(k_sub)+vbdat(id,jd,k_sub)

              pt1d(k_sub)=pt1d(k_sub)+ptbdat(id,jd,k_sub)

            end do
            end do

            u1d(k_sub)=u1d(k_sub)*nginv
            v1d(k_sub)=v1d(k_sub)*nginv

            pt1d(k_sub)=pt1d(k_sub)*nginv

          end do

!$omp end do

        else if(fproc(1:6).eq.'refsfc') then

!$omp do schedule(runtime) private(id,jd,k_sub,ngcnt,ngcntv)

          do k_sub=1,kref-1

            ngcnt=0

            do jd=jdstr,jdend
            do id=idstr,idend

              if(z(k_sub).ge.zdat(id,jd,1)) then

                ngcnt=ngcnt+1

                u1d(k_sub)=u1d(k_sub)+ubdat(id,jd,k_sub)
                v1d(k_sub)=v1d(k_sub)+vbdat(id,jd,k_sub)

                pt1d(k_sub)=pt1d(k_sub)+ptbdat(id,jd,k_sub)

              end if

            end do
            end do

            if(ngcnt.eq.0) then

              u1d(k_sub)=lim35n
              v1d(k_sub)=lim35n

              pt1d(k_sub)=lim35n

            else

              ngcntv=1.e0/real(ngcnt)

              u1d(k_sub)=u1d(k_sub)*ngcntv
              v1d(k_sub)=v1d(k_sub)*ngcntv

              pt1d(k_sub)=pt1d(k_sub)*ngcntv

            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(id,jd,k_sub)

          do k_sub=kref,nk

            do jd=jdstr,jdend
            do id=idstr,idend
              u1d(k_sub)=u1d(k_sub)+ubdat(id,jd,k_sub)
              v1d(k_sub)=v1d(k_sub)+vbdat(id,jd,k_sub)

              pt1d(k_sub)=pt1d(k_sub)+ptbdat(id,jd,k_sub)

            end do
            end do

            u1d(k_sub)=u1d(k_sub)*nginv
            v1d(k_sub)=v1d(k_sub)*nginv

            pt1d(k_sub)=pt1d(k_sub)*nginv

          end do

!$omp end do

        end if

      else

        if(fproc(1:3).eq.'all') then

!$omp do schedule(runtime) private(id,jd,k_sub)

          do k_sub=1,nk

            do jd=jdstr,jdend
            do id=idstr,nid
              u1d(k_sub)=u1d(k_sub)+ubdat(id,jd,k_sub)
              v1d(k_sub)=v1d(k_sub)+vbdat(id,jd,k_sub)

              pt1d(k_sub)=pt1d(k_sub)+ptbdat(id,jd,k_sub)

            end do
            end do

            do jd=jdstr,jdend
            do id=1,idend
              u1d(k_sub)=u1d(k_sub)+ubdat(id,jd,k_sub)
              v1d(k_sub)=v1d(k_sub)+vbdat(id,jd,k_sub)

              pt1d(k_sub)=pt1d(k_sub)+ptbdat(id,jd,k_sub)

            end do
            end do

            u1d(k_sub)=u1d(k_sub)*nginv
            v1d(k_sub)=v1d(k_sub)*nginv

            pt1d(k_sub)=pt1d(k_sub)*nginv

          end do

!$omp end do

        else if(fproc(1:6).eq.'refsfc') then

!$omp do schedule(runtime) private(id,jd,k_sub,ngcnt,ngcntv)

          do k_sub=1,kref-1

            ngcnt=0

            do jd=jdstr,jdend
            do id=idstr,nid

              if(z(k_sub).ge.zdat(id,jd,1)) then

                ngcnt=ngcnt+1

                u1d(k_sub)=u1d(k_sub)+ubdat(id,jd,k_sub)
                v1d(k_sub)=v1d(k_sub)+vbdat(id,jd,k_sub)

                pt1d(k_sub)=pt1d(k_sub)+ptbdat(id,jd,k_sub)

              end if

            end do
            end do

            do jd=jdstr,jdend
            do id=1,idend

              if(z(k_sub).ge.zdat(id,jd,1)) then

                ngcnt=ngcnt+1

                u1d(k_sub)=u1d(k_sub)+ubdat(id,jd,k_sub)
                v1d(k_sub)=v1d(k_sub)+vbdat(id,jd,k_sub)

                pt1d(k_sub)=pt1d(k_sub)+ptbdat(id,jd,k_sub)

              end if

            end do
            end do

            if(ngcnt.eq.0) then

              u1d(k_sub)=lim35n
              v1d(k_sub)=lim35n

              pt1d(k_sub)=lim35n

            else

              ngcntv=1.e0/real(ngcnt)

              u1d(k_sub)=u1d(k_sub)*ngcntv
              v1d(k_sub)=v1d(k_sub)*ngcntv

              pt1d(k_sub)=pt1d(k_sub)*ngcntv

            end if

          end do

!$omp end do

!$omp do schedule(runtime) private(id,jd,k_sub)

          do k_sub=kref,nk

            do jd=jdstr,jdend
            do id=idstr,nid
              u1d(k_sub)=u1d(k_sub)+ubdat(id,jd,k_sub)
              v1d(k_sub)=v1d(k_sub)+vbdat(id,jd,k_sub)

              pt1d(k_sub)=pt1d(k_sub)+ptbdat(id,jd,k_sub)

            end do
            end do

            do jd=jdstr,jdend
            do id=1,idend
              u1d(k_sub)=u1d(k_sub)+ubdat(id,jd,k_sub)
              v1d(k_sub)=v1d(k_sub)+vbdat(id,jd,k_sub)

              pt1d(k_sub)=pt1d(k_sub)+ptbdat(id,jd,k_sub)

            end do
            end do

            u1d(k_sub)=u1d(k_sub)*nginv
            v1d(k_sub)=v1d(k_sub)*nginv

            pt1d(k_sub)=pt1d(k_sub)*nginv

          end do

!$omp end do

        end if

      end if

! -----

! Get the averaged water vapor mixing ratio.

      if(gpvvar(2:2).eq.'o') then

        if(idstr.le.idend) then

          if(fproc(1:3).eq.'all') then

!$omp do schedule(runtime) private(id,jd,k_sub)

            do k_sub=1,nk

              do jd=jdstr,jdend
              do id=idstr,idend
                qv1d(k_sub)=qv1d(k_sub)+qvbdat(id,jd,k_sub)
              end do
              end do

              qv1d(k_sub)=qv1d(k_sub)*nginv

            end do

!$omp end do

          else if(fproc(1:6).eq.'refsfc') then

!$omp do schedule(runtime) private(id,jd,k_sub,ngcnt)

            do k_sub=1,kref-1

              ngcnt=0

              do jd=jdstr,jdend
              do id=idstr,idend

                if(z(k_sub).ge.zdat(id,jd,1)) then

                  ngcnt=ngcnt+1

                  qv1d(k_sub)=qv1d(k_sub)+qvbdat(id,jd,k_sub)

                end if

              end do
              end do

              if(ngcnt.eq.0) then
                qv1d(k_sub)=lim35n
              else
                qv1d(k_sub)=qv1d(k_sub)/real(ngcnt)
              end if

            end do

!$omp end do

!$omp do schedule(runtime) private(id,jd,k_sub)

            do k_sub=kref,nk

              do jd=jdstr,jdend
              do id=idstr,idend
                qv1d(k_sub)=qv1d(k_sub)+qvbdat(id,jd,k_sub)
              end do
              end do

              qv1d(k_sub)=qv1d(k_sub)*nginv

            end do

!$omp end do

          end if

        else

          if(fproc(1:3).eq.'all') then

!$omp do schedule(runtime) private(id,jd,k_sub)

            do k_sub=1,nk

              do jd=jdstr,jdend
              do id=idstr,nid
                qv1d(k_sub)=qv1d(k_sub)+qvbdat(id,jd,k_sub)
              end do
              end do

              do jd=jdstr,jdend
              do id=1,idend
                qv1d(k_sub)=qv1d(k_sub)+qvbdat(id,jd,k_sub)
              end do
              end do

              qv1d(k_sub)=qv1d(k_sub)*nginv

            end do

!$omp end do

          else if(fproc(1:6).eq.'refsfc') then

!$omp do schedule(runtime) private(id,jd,k_sub,ngcnt)

            do k_sub=1,kref-1

              ngcnt=0

              do jd=jdstr,jdend
              do id=idstr,nid

                if(z(k_sub).ge.zdat(id,jd,1)) then

                  ngcnt=ngcnt+1

                  qv1d(k_sub)=qv1d(k_sub)+qvbdat(id,jd,k_sub)

                end if

              end do
              end do

              do jd=jdstr,jdend
              do id=1,idend

                if(z(k_sub).ge.zdat(id,jd,1)) then

                  ngcnt=ngcnt+1

                  qv1d(k_sub)=qv1d(k_sub)+qvbdat(id,jd,k_sub)

                end if

              end do
              end do

              if(ngcnt.eq.0) then
                qv1d(k_sub)=lim35n
              else
                qv1d(k_sub)=qv1d(k_sub)/real(ngcnt)
              end if

            end do

!$omp end do

!$omp do schedule(runtime) private(id,jd,k_sub)

            do k_sub=kref,nk

              do jd=jdstr,jdend
              do id=idstr,nid
                qv1d(k_sub)=qv1d(k_sub)+qvbdat(id,jd,k_sub)
              end do
              end do

              do jd=jdstr,jdend
              do id=1,idend
                qv1d(k_sub)=qv1d(k_sub)+qvbdat(id,jd,k_sub)
              end do
              end do

              qv1d(k_sub)=qv1d(k_sub)*nginv

            end do

!$omp end do

          end if

        end if

      end if

! -----

!$omp end parallel

!! -----

! Get the reference averaged pressure.

      pref=pref*nginv

! -----

!!! -----

! Extrapolate the base state variables in undefined layer.

      call s_vextund(idgpvvar,idetrvar_gpv,                             &
     &               nk,z,u1d(1),v1d(1),pt1d(1),qv1d(1))

! -----

! Perform interpolation vertically to fine plane.

      call s_copy1d(1,nk,u1d(1),tmp1)

      call s_vint11(4*nk-3,z1d(1),u1d(1),nk,z,tmp1)

      call s_copy1d(1,nk,v1d(1),tmp1)

      call s_vint11(4*nk-3,z1d(1),v1d(1),nk,z,tmp1)

      call s_copy1d(1,nk,pt1d(1),tmp1)

      call s_vint11(4*nk-3,z1d(1),pt1d(1),nk,z,tmp1)

      if(gpvvar(2:2).eq.'o') then

        call s_copy1d(1,nk,qv1d(1),tmp1)

        call s_vint11(4*nk-3,z1d(1),qv1d(1),nk,z,tmp1)

      end if

! -----

!!!! -----

!! Get the base state pressure.

! Calculate the base state virtual potential temperature.

!$omp parallel default(shared)

      if(gpvvar(2:2).eq.'o') then

!$omp do schedule(runtime) private(k_sub)

        do k_sub=1,4*nk-3

          pt1d(k_sub)=pt1d(k_sub)                                       &
     &      *(1.e0+epsav*qv1d(k_sub))/(1.e0+qv1d(k_sub))

        end do

!$omp end do

      end if

!$omp end parallel

! -----

! Calculate the base state Exner function.

      pimin=lim36

      p1d(4*kref-3)=exp(rddvcp*log(p0iv*pref))

      do k=4*(kref-1),1,-1
        p1d(k)=p1d(k+1)+gdvcp2*(z1d(k+1)-z1d(k))/(pt1d(k)+pt1d(k+1))

        pimin=min(p1d(k),pimin)

      end do

      do k=4*kref-2,4*nk-3
        p1d(k)=p1d(k-1)-gdvcp2*(z1d(k)-z1d(k-1))/(pt1d(k-1)+pt1d(k))

        pimin=min(p1d(k),pimin)

      end do

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

          call destroy('get1d   ',5,'cont',7,'              ',14,101,   &
     &                 stat)

        end if

        call cpondpe

        call destroy('get1d   ',5,'stop',1001,'              ',14,101,  &
     &               stat)

      end if

! -----

! Recalculate the base state potential temperature and get the base
! state pressure.

!$omp parallel default(shared)

      if(gpvvar(2:2).eq.'o') then

!$omp do schedule(runtime) private(k_sub)

        do k_sub=1,4*nk-3

          pt1d(k_sub)=pt1d(k_sub)                                       &
     &      *(1.e0+qv1d(k_sub))/(1.e0+epsav*qv1d(k_sub))

        end do

!$omp end do

      end if

!$omp do schedule(runtime) private(k_sub)

      do k_sub=1,4*nk-3
        p1d(k_sub)=p0*exp(cpdvrd*log(p1d(k_sub)))
      end do

!$omp end do

!$omp end parallel

! -----

!! -----

! Perform extrapolation in the lowest layer.

      call vextlow(etrvar_gpv(1:6),z1d,u1d,v1d,p1d,pt1d,qv1d)

! -----

      end subroutine s_get1d

!-----7--------------------------------------------------------------7--

      end module m_get1d
