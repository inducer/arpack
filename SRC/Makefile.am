lib_LTLIBRARIES     = libarpack.la
libarpack_la_LDFLAGS = -version-info 2

libarpack_la_SOURCES = \
	debug.h stat.h version.h \
	\
        sgetv0.f slaqrb.f sstqrb.f ssortc.f ssortr.f sstatn.f sstats.f \
	snaitr.f snapps.f snaup2.f snaupd.f snconv.f sneigh.f sngets.f \
	ssaitr.f ssapps.f ssaup2.f ssaupd.f ssconv.f sseigt.f ssgets.f \
	sneupd.f sseupd.f ssesrt.f \
	\
        dgetv0.f dlaqrb.f dstqrb.f dsortc.f dsortr.f dstatn.f dstats.f \
	dnaitr.f dnapps.f dnaup2.f dnaupd.f dnconv.f dneigh.f dngets.f \
	dsaitr.f dsapps.f dsaup2.f dsaupd.f dsconv.f dseigt.f dsgets.f \
	dneupd.f dseupd.f dsesrt.f \
	\
        cnaitr.f cnapps.f cnaup2.f cnaupd.f cneigh.f cneupd.f cngets.f \
        cgetv0.f csortc.f cstatn.f \
        \
        znaitr.f znapps.f znaup2.f znaupd.f zneigh.f zneupd.f zngets.f \
        zgetv0.f zsortc.f zstatn.f \
	../UTIL/arscnd.f  ../UTIL/dvout.f   ../UTIL/iswap.f  ../UTIL/zmout.f \
	../UTIL/cmout.f   ../UTIL/icnteq.f  ../UTIL/ivout.f  ../UTIL/zvout.f \
	../UTIL/cvout.f   ../UTIL/icopy.f   ../UTIL/smout.f \
	../UTIL/dmout.f   ../UTIL/iset.f    ../UTIL/svout.f \
	\
	../LAPACK/cgbtf2.f ../LAPACK/cgbtrf.f ../LAPACK/cgbtrs.f	\
	../LAPACK/cgeqr2.f ../LAPACK/cgttrf.f ../LAPACK/cgttrs.f	\
	../LAPACK/clacgv.f ../LAPACK/clacon.f ../LAPACK/clacpy.f	\
	../LAPACK/cladiv.f ../LAPACK/clahqr.f ../LAPACK/clange.f	\
	../LAPACK/clanhs.f ../LAPACK/clarf.f ../LAPACK/clarfg.f		\
	../LAPACK/clarnv.f ../LAPACK/clartg.f ../LAPACK/clascl.f	\
	../LAPACK/claset.f ../LAPACK/classq.f ../LAPACK/claswp.f	\
	../LAPACK/clatrs.f ../LAPACK/cmach.f ../LAPACK/crot.f		\
	../LAPACK/ctrevc.f ../LAPACK/ctrexc.f ../LAPACK/ctrsen.f	\
	../LAPACK/ctrsyl.f ../LAPACK/cunm2r.f ../LAPACK/dgbtf2.f	\
	../LAPACK/dgbtrf.f ../LAPACK/dgbtrs.f ../LAPACK/dgeqr2.f	\
	../LAPACK/dgttrf.f ../LAPACK/dgttrs.f ../LAPACK/dlabad.f	\
	../LAPACK/dlacon.f ../LAPACK/dlacpy.f ../LAPACK/dladiv.f	\
	../LAPACK/dlae2.f ../LAPACK/dlaev2.f ../LAPACK/dlaexc.f		\
	../LAPACK/dlagtm.f ../LAPACK/dlahqr.f ../LAPACK/dlaln2.f	\
	../LAPACK/dlamch.f ../LAPACK/dlange.f ../LAPACK/dlanhs.f	\
	../LAPACK/dlanst.f ../LAPACK/dlanv2.f ../LAPACK/dlaptm.f	\
	../LAPACK/dlapy2.f ../LAPACK/dlapy3.f ../LAPACK/dlaran.f	\
	../LAPACK/dlarf.f ../LAPACK/dlarfg.f ../LAPACK/dlarfx.f		\
	../LAPACK/dlarnd.f ../LAPACK/dlarnv.f ../LAPACK/dlartg.f	\
	../LAPACK/dlaruv.f ../LAPACK/dlascl.f ../LAPACK/dlaset.f	\
	../LAPACK/dlasr.f ../LAPACK/dlasrt.f ../LAPACK/dlassq.f		\
	../LAPACK/dlaswp.f ../LAPACK/dlasy2.f ../LAPACK/dorm2r.f	\
	../LAPACK/dpttrf.f ../LAPACK/dpttrs.f ../LAPACK/dsteqr.f	\
	../LAPACK/dtrevc.f ../LAPACK/dtrexc.f ../LAPACK/dtrsen.f	\
	../LAPACK/dtrsyl.f ../LAPACK/dzsum1.f ../LAPACK/icmax1.f	\
	../LAPACK/ilaenv.f ../LAPACK/izmax1.f ../LAPACK/lsame.f		\
	../LAPACK/lsamen.f ../LAPACK/scsum1.f ../LAPACK/sgbtf2.f	\
	../LAPACK/sgbtrf.f ../LAPACK/sgbtrs.f ../LAPACK/sgeqr2.f	\
	../LAPACK/sgttrf.f ../LAPACK/sgttrs.f ../LAPACK/slabad.f	\
	../LAPACK/slacon.f ../LAPACK/slacpy.f ../LAPACK/sladiv.f	\
	../LAPACK/slae2.f ../LAPACK/slaev2.f ../LAPACK/slaexc.f		\
	../LAPACK/slagtm.f ../LAPACK/slahqr.f ../LAPACK/slaln2.f	\
	../LAPACK/slamch.f ../LAPACK/slange.f ../LAPACK/slanhs.f	\
	../LAPACK/slanst.f ../LAPACK/slanv2.f ../LAPACK/slaptm.f	\
	../LAPACK/slapy2.f ../LAPACK/slapy3.f ../LAPACK/slaran.f	\
	../LAPACK/slarf.f ../LAPACK/slarfg.f ../LAPACK/slarfx.f		\
	../LAPACK/slarnd.f ../LAPACK/slarnv.f ../LAPACK/slartg.f	\
	../LAPACK/slaruv.f ../LAPACK/slascl.f ../LAPACK/slaset.f	\
	../LAPACK/slasr.f ../LAPACK/slasrt.f ../LAPACK/slassq.f		\
	../LAPACK/slaswp.f ../LAPACK/slasy2.f ../LAPACK/sorm2r.f	\
	../LAPACK/spttrf.f ../LAPACK/spttrs.f ../LAPACK/ssteqr.f	\
	../LAPACK/strevc.f ../LAPACK/strexc.f ../LAPACK/strsen.f	\
	../LAPACK/strsyl.f ../LAPACK/xerbla.f ../LAPACK/xlaenv.f	\
	../LAPACK/zdscal.f ../LAPACK/zgbtf2.f ../LAPACK/zgbtrf.f	\
	../LAPACK/zgbtrs.f ../LAPACK/zgeqr2.f ../LAPACK/zgttrf.f	\
	../LAPACK/zgttrs.f ../LAPACK/zlacgv.f ../LAPACK/zlacon.f	\
	../LAPACK/zlacpy.f ../LAPACK/zladiv.f ../LAPACK/zlahqr.f	\
	../LAPACK/zlange.f ../LAPACK/zlanhs.f ../LAPACK/zlarf.f		\
	../LAPACK/zlarfg.f ../LAPACK/zlarnv.f ../LAPACK/zlartg.f	\
	../LAPACK/zlascl.f ../LAPACK/zlaset.f ../LAPACK/zlassq.f	\
	../LAPACK/zlaswp.f ../LAPACK/zlatrs.f ../LAPACK/zmach.f		\
	../LAPACK/zrot.f ../LAPACK/ztrevc.f ../LAPACK/ztrexc.f		\
	../LAPACK/ztrsen.f ../LAPACK/ztrsyl.f ../LAPACK/zunm2r.f
