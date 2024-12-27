!***********************************************************************
      module m_setname
!***********************************************************************

!     Author      : Sakakibara Atsushi
!     Date        : 2007/01/20
!     Modification: 2007/01/31, 2007/03/23, 2007/04/24, 2007/05/21,
!                   2007/07/30, 2007/09/04, 2007/10/19, 2008/01/11,
!                   2008/04/17, 2008/05/02, 2008/06/09, 2008/07/01,
!                   2008/08/25, 2008/10/10, 2008/12/11, 2009/01/05,
!                   2009/01/30, 2009/02/27, 2010/05/17, 2011/05/16,
!                   2011/08/09, 2011/08/18, 2011/09/22, 2011/11/10,
!                   2013/01/28, 2013/03/27

!     Author      : Satoki Tsujino
!     Modification: 2016/04/08

!-----7--1----+----2----+----3----+----4----+----5----+----6----+----7--

! In this module,
!     set the table to archive namelist variables.

!-----7--------------------------------------------------------------7--

! Module reference

      use m_castname
      use m_comindx
      use m_commpi
      use m_comname
      use m_defname

!-----7--------------------------------------------------------------7--

! Implicit typing

      implicit none

! Default access control

      private

! Exceptional access control

      public :: setname, s_setname

!-----7--------------------------------------------------------------7--

! Module variable

!     none

! Module procedure

      interface setname

        module procedure s_setname

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
      subroutine s_setname
!***********************************************************************

! Internal private variable

      integer iid      ! Index of do loops

!-----7--------------------------------------------------------------7--

!! Set the table to archive namelist variables.

      if(mype.eq.root) then

! For the character variables.

        cname(idexprim)(1:108)=exprim(1:108)

        cname(idcrsdir)(1:108)=crsdir(1:108)
        cname(iddatdir)(1:108)=datdir(1:108)

        cname(ididate)(1:12)=idate(1:12)

        cname(idlbcvar)(1:10)=lbcvar(1:10)

        cname(idgpvvar)(1:9)=gpvvar(1:9)
        cname(idnggvar)(1:8)=nggvar(1:8)
        cname(idexbvar)(1:9)=exbvar(1:9)
        cname(idlspvar)(1:10)=lspvar(1:10)
        cname(idvspvar)(1:10)=vspvar(1:10)

        cname(idrdrvar)(1:4)=rdrvar(1:4)
        cname(idngrvar)(1:5)=ngrvar(1:5)

        cname(idprvres)(1:108)=prvres(1:108)
        cname(idsfcdat)(1:3)=sfcdat(1:3)

        cname(idsndtyp)(1:3)=sndtyp(1:3)

        cname(iddmpvar)(1:18)=dmpvar(1:18)
        cname(idmxnvar)(1:12)=mxnvar(1:12)

        cname(iddatype_gpv)(1:2)=datype_gpv(1:2)
        cname(idetrvar_gpv)(1:7)=etrvar_gpv(1:7)

        cname(iddatype_rdr)(1:1)=datype_rdr(1:1)

        cname(idfltyp_uni)(1:3)=fltyp_uni(1:3)

! -----

! For the integer variables.

        iname(idzeropt)=0
        iname(idoneopt)=1

        iname(idwlngth)=wlngth
        iname(idsavmem)=savmem

        iname(idncexp)=ncexp

        iname(idnccrs)=nccrs
        iname(idncdat)=ncdat

        iname(idnumpe)=numpe

        iname(idxgroup)=xgroup
        iname(idygroup)=ygroup

        iname(idxsub)=xsub
        iname(idysub)=ysub

        iname(idxdim)=xdim
        iname(idydim)=ydim
        iname(idzdim)=zdim

        iname(idmpopt)=mpopt
        iname(idnspol)=nspol

        iname(idsthopt)=sthopt

        iname(idtrnopt)=trnopt

        iname(idwbc)=wbc
        iname(idebc)=ebc
        iname(idsbc)=sbc
        iname(idnbc)=nbc
        iname(idbbc)=bbc
        iname(idtbc)=tbc

        iname(idgsmopt)=gsmopt
        iname(idgsmcnt)=gsmcnt
        iname(idnggopt)=nggopt
        iname(idexbopt)=exbopt
        iname(idexbwid)=exbwid
        iname(idlspopt)=lspopt
        iname(idwdnews)=wdnews
        iname(idwdnorm)=wdnorm
        iname(idvspopt)=vspopt

        iname(idngropt)=ngropt

        iname(idsfcopt)=sfcopt
        iname(idlevpbl)=levpbl
        iname(idlevund)=levund
        iname(idncprv)=ncprv
        iname(idlnduse)=lnduse
        iname(iddstopt)=dstopt

        iname(idiniopt)=iniopt
        iname(idsnddim)=snddim
        iname(idmasopt)=masopt

        iname(idmovopt)=movopt

        iname(idpt0opt)=pt0opt
        iname(idpt0num)=pt0num

        iname(idgwmopt)=gwmopt
        iname(idimpopt)=impopt
        iname(idadvopt)=advopt

        iname(idsmtopt)=smtopt

        iname(idmfcopt)=mfcopt

        iname(idcoropt)=coropt

        iname(idcrvopt)=crvopt

        iname(idbuyopt)=buyopt

        iname(iddiaopt)=diaopt

        iname(iddivopt)=divopt

        iname(idcphopt)=cphopt
        iname(idhaiopt)=haiopt
        iname(idncbinw)=ncbinw
        iname(idqcgopt)=qcgopt

        iname(idaslopt)=aslopt

        iname(idtrkopt)=trkopt
        iname(idqt0opt)=qt0opt
        iname(idqt0num)=qt0num

        iname(idtubopt)=tubopt
        iname(idisoopt)=isoopt

        iname(iddmpfmt)=dmpfmt
        iname(iddmplev)=dmplev
        iname(iddmpmon)=dmpmon
        iname(idresopt)=resopt
        iname(idmxnopt)=mxnopt

        iname(idiwest)=iwest
        iname(idieast)=ieast
        iname(idjsouth)=jsouth
        iname(idjnorth)=jnorth

        iname(idmpopt_gpv)=mpopt_gpv
        iname(idnspol_gpv)=nspol_gpv

        iname(idxdim_gpv)=xdim_gpv
        iname(idydim_gpv)=ydim_gpv
        iname(idzdim_gpv)=zdim_gpv

        iname(idintopt_gpv)=intopt_gpv
        iname(idrotopt_gpv)=rotopt_gpv
        iname(idrefsfc_gpv)=refsfc_gpv

        iname(idmpopt_asl)=mpopt_asl
        iname(idnspol_asl)=nspol_asl

        iname(idxdim_asl)=xdim_asl
        iname(idydim_asl)=ydim_asl
        iname(idzdim_asl)=zdim_asl

        iname(idintopt_asl)=intopt_asl

        iname(idmpopt_rdr)=mpopt_rdr
        iname(idnspol_rdr)=nspol_rdr

        iname(idxdim_rdr)=xdim_rdr
        iname(idydim_rdr)=ydim_rdr
        iname(idzdim_rdr)=zdim_rdr

        iname(idrotopt_rdr)=rotopt_rdr

        iname(idmpopt_trn)=mpopt_trn
        iname(idnspol_trn)=nspol_trn

        iname(idxdim_trn)=xdim_trn
        iname(idydim_trn)=ydim_trn

        iname(idintopt_trn)=intopt_trn

        iname(idmpopt_lnd)=mpopt_lnd
        iname(idnspol_lnd)=nspol_lnd

        iname(idxdim_lnd)=xdim_lnd
        iname(idydim_lnd)=ydim_lnd

        iname(idintopt_lnd)=intopt_lnd
        iname(idnumctg_lnd)=numctg_lnd

        iname(idmpopt_sst)=mpopt_sst
        iname(idnspol_sst)=nspol_sst

        iname(idxdim_sst)=xdim_sst
        iname(idydim_sst)=ydim_sst

        iname(idmpopt_ice)=mpopt_ice
        iname(idnspol_ice)=nspol_ice

        iname(idxdim_ice)=xdim_ice
        iname(idydim_ice)=ydim_ice

        iname(idrmopt_uni)=rmopt_uni
        iname(iduniopt_uni)=uniopt_uni
        iname(idugroup_uni)=ugroup_uni

        iname(idxsub_rst)=xsub_rst
        iname(idysub_rst)=ysub_rst
        iname(idrmopt_rst)=rmopt_rst

! -----

! For the real variables.

        rname(idtlat1)=tlat1
        rname(idtlat2)=tlat2
        rname(idtlon)=tlon
        rname(iddisr)=disr

        rname(iddx)=dx
        rname(iddy)=dy
        rname(iddz)=dz
        rname(iddxiv)=dxiv
        rname(iddyiv)=dyiv
        rname(iddziv)=dziv
        rname(idulat)=ulat
        rname(idulon)=ulon
        rname(idriu)=riu
        rname(idrju)=rju

        rname(idzsfc)=zsfc
        rname(idzflat)=zflat
        rname(iddzmin)=dzmin
        rname(idlayer1)=layer1
        rname(idlayer2)=layer2

        rname(idmnthgh)=mnthgh(1)
        rname(idmnthgh+1)=mnthgh(2)
        rname(idmntwx)=mntwx
        rname(idmntwy)=mntwy
        rname(idmntcx)=mntcx
        rname(idmntcy)=mntcy

        rname(idstime)=stime
        rname(idetime)=etime

        rname(idlbnews)=lbnews
        rname(idlbnorm)=lbnorm
        rname(idgwave)=gwave

        rname(idgpvitv)=gpvitv
        rname(idaslitv)=aslitv
        rname(idgsmcoe)=gsmcoe
        rname(idnggcoe)=nggcoe
        rname(idnggdlt)=nggdlt
        rname(idnggstr)=nggstr
        rname(idnggend)=nggend
        rname(idnggc20)=nggc20
        rname(idexnews)=exnews
        rname(idexnorm)=exnorm
        rname(idlspsmt)=lspsmt
        rname(idlsnews)=lsnews
        rname(idlsnorm)=lsnorm
        rname(idvspgpv)=vspgpv
        rname(idvspbar)=vspbar
        rname(idbotgpv)=botgpv
        rname(idbotbar)=botbar

        rname(idrdritv)=rdritv
        rname(idngrcoe)=ngrcoe
        rname(idngrdlt)=ngrdlt
        rname(idngrstr)=ngrstr
        rname(idngrend)=ngrend
        rname(idngrc20)=ngrc20
        rname(idngraff)=ngraff

        rname(iddtgrd)=dtgrd
        rname(iddzgrd)=dzgrd
        rname(iddzsea)=dzsea
        rname(idtgdeep)=tgdeep
        rname(idgralbe)=gralbe
        rname(idgrbeta)=grbeta
        rname(idgrz0m)=grz0m
        rname(idgrz0h)=grz0h
        rname(idgrcap)=grcap
        rname(idgrnuu)=grnuu
        rname(idsstcst)=sstcst
        rname(idsstitv)=sstitv

        rname(idzsnd0)=zsnd0
        rname(idpsnd0)=psnd0
        rname(idmaseps)=maseps
        rname(idalpha1)=alpha1
        rname(idalpha2)=alpha2

        rname(idumove)=umove
        rname(idvmove)=vmove

        rname(idptp0)=ptp0
        rname(idpt0rx)=pt0rx
        rname(idpt0ry)=pt0ry
        rname(idpt0rz)=pt0rz
        rname(idpt0cx)=pt0cx
        rname(idpt0cy)=pt0cy
        rname(idpt0cz)=pt0cz
        rname(idpt0ds)=pt0ds

        rname(iddtbig)=dtbig
        rname(iddtsml)=dtsml
        rname(idgsdeps)=gsdeps
        rname(idweicoe)=weicoe
        rname(idfilcoe)=filcoe
        rname(iddtvcul)=dtvcul

        rname(idsmhcoe)=smhcoe
        rname(idsmvcoe)=smvcoe
        rname(idnlhcoe)=nlhcoe
        rname(idnlvcoe)=nlvcoe

        rname(idthresq)=thresq
        rname(iddtcmph)=dtcmph(1)
        rname(iddtcmph+1)=dtcmph(2)
        rname(iddtcmph+2)=dtcmph(3)
        rname(iddtcmph+3)=dtcmph(4)
        rname(iddtcmph+4)=dtcmph(5)
        rname(iddtcmph+5)=dtcmph(6)
        rname(idbbinw)=bbinw
        rname(idsbinw)=sbinw
        rname(ideledlt)=eledlt

        rname(idqt0)=qt0
        rname(idqt0rx)=qt0rx
        rname(idqt0ry)=qt0ry
        rname(idqt0rz)=qt0rz
        rname(idqt0cx)=qt0cx
        rname(idqt0cy)=qt0cy
        rname(idqt0cz)=qt0cz
        rname(idqt0ds)=qt0ds
        rname(idqtdt)=qtdt
        rname(idqtstr)=qtstr
        rname(idqtend)=qtend

        rname(iddmpitv)=dmpitv
        rname(idmonitv)=monitv
        rname(idresitv)=resitv
        rname(idmxnitv)=mxnitv

        rname(idtlat1_gpv)=tlat1_gpv
        rname(idtlat2_gpv)=tlat2_gpv
        rname(idtlon_gpv)=tlon_gpv

        rname(iddx_gpv)=dx_gpv
        rname(iddy_gpv)=dy_gpv
        rname(iddxiv_gpv)=dxiv_gpv
        rname(iddyiv_gpv)=dyiv_gpv
        rname(idulat_gpv)=ulat_gpv
        rname(idulon_gpv)=ulon_gpv
        rname(idriu_gpv)=riu_gpv
        rname(idrju_gpv)=rju_gpv

        rname(idtlat1_asl)=tlat1_asl
        rname(idtlat2_asl)=tlat2_asl
        rname(idtlon_asl)=tlon_asl

        rname(iddx_asl)=dx_asl
        rname(iddy_asl)=dy_asl
        rname(iddxiv_asl)=dxiv_asl
        rname(iddyiv_asl)=dyiv_asl
        rname(idulat_asl)=ulat_asl
        rname(idulon_asl)=ulon_asl
        rname(idriu_asl)=riu_asl
        rname(idrju_asl)=rju_asl

        rname(idtlat1_rdr)=tlat1_rdr
        rname(idtlat2_rdr)=tlat2_rdr
        rname(idtlon_rdr)=tlon_rdr

        rname(iddx_rdr)=dx_rdr
        rname(iddy_rdr)=dy_rdr
        rname(iddxiv_rdr)=dxiv_rdr
        rname(iddyiv_rdr)=dyiv_rdr
        rname(idulat_rdr)=ulat_rdr
        rname(idulon_rdr)=ulon_rdr
        rname(idriu_rdr)=riu_rdr
        rname(idrju_rdr)=rju_rdr

        rname(idrdrcoe_rdr)=rdrcoe_rdr
        rname(idrdrexp_rdr)=rdrexp_rdr

        rname(idtlat1_trn)=tlat1_trn
        rname(idtlat2_trn)=tlat2_trn
        rname(idtlon_trn)=tlon_trn

        rname(iddx_trn)=dx_trn
        rname(iddy_trn)=dy_trn
        rname(iddxiv_trn)=dxiv_trn
        rname(iddyiv_trn)=dyiv_trn
        rname(idulat_trn)=ulat_trn
        rname(idulon_trn)=ulon_trn
        rname(idriu_trn)=riu_trn
        rname(idrju_trn)=rju_trn

        rname(idtlat1_lnd)=tlat1_lnd
        rname(idtlat2_lnd)=tlat2_lnd
        rname(idtlon_lnd)=tlon_lnd

        rname(iddx_lnd)=dx_lnd
        rname(iddy_lnd)=dy_lnd
        rname(iddxiv_lnd)=dxiv_lnd
        rname(iddyiv_lnd)=dyiv_lnd
        rname(idulat_lnd)=ulat_lnd
        rname(idulon_lnd)=ulon_lnd
        rname(idriu_lnd)=riu_lnd
        rname(idrju_lnd)=rju_lnd

        rname(idtlat1_sst)=tlat1_sst
        rname(idtlat2_sst)=tlat2_sst
        rname(idtlon_sst)=tlon_sst

        rname(iddx_sst)=dx_sst
        rname(iddy_sst)=dy_sst
        rname(iddxiv_sst)=dxiv_sst
        rname(iddyiv_sst)=dyiv_sst
        rname(idulat_sst)=ulat_sst
        rname(idulon_sst)=ulon_sst
        rname(idriu_sst)=riu_sst
        rname(idrju_sst)=rju_sst

        rname(idtlat1_ice)=tlat1_ice
        rname(idtlat2_ice)=tlat2_ice
        rname(idtlon_ice)=tlon_ice

        rname(iddx_ice)=dx_ice
        rname(iddy_ice)=dy_ice
        rname(iddxiv_ice)=dxiv_ice
        rname(iddyiv_ice)=dyiv_ice
        rname(idulat_ice)=ulat_ice
        rname(idulon_ice)=ulon_ice
        rname(idriu_ice)=riu_ice
        rname(idrju_ice)=rju_ice

        rname(idflitv_uni)=flitv_uni

        rname(idflitv_rst)=flitv_rst

! -----

! For the table.

!$omp parallel default(shared)

!$omp do schedule(runtime) private(iid)

        do iid=0,numctg_lnd-1
          iname(idlnduse_lnd+iid)=lnduse_lnd(iid+1)
        end do

!$omp end do

!$omp do schedule(runtime) private(iid)

        do iid=0,numctg_lnd-1
          rname(idalbe_lnd+iid)=albe_lnd(iid+1)
        end do

!$omp end do

!$omp do schedule(runtime) private(iid)

        do iid=0,numctg_lnd-1
          rname(idbeta_lnd+iid)=beta_lnd(iid+1)
        end do

!$omp end do

!$omp do schedule(runtime) private(iid)

        do iid=0,numctg_lnd-1
          rname(idz0m_lnd+iid)=z0m_lnd(iid+1)
        end do

!$omp end do

!$omp do schedule(runtime) private(iid)

        do iid=0,numctg_lnd-1
          rname(idz0h_lnd+iid)=z0h_lnd(iid+1)
        end do

!$omp end do

!$omp do schedule(runtime) private(iid)

        do iid=0,numctg_lnd-1
          rname(idcap_lnd+iid)=cap_lnd(iid+1)
        end do

!$omp end do

!$omp do schedule(runtime) private(iid)

        do iid=0,numctg_lnd-1
          rname(idnuu_lnd+iid)=nuu_lnd(iid+1)
        end do

!$omp end do

!$omp end parallel

! -----

      end if

!! -----

! Broadcast the namelist variables.

      call castname

! -----

      end subroutine s_setname

!-----7--------------------------------------------------------------7--

      end module m_setname
