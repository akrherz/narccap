#	Top-level Makefile for INTERPB

#	Macros, these should be generic for all machines

.IGNORE:

AR	=	ar ru
CD	=	cd
LN	=	ln -s
MAKE	=	make -i -f Makefile
RM	=	/bin/rm -f
RM_LIST	=	*.o *.M core *.kmo *.mod interpb

#	Targets for supported architectures

default:
	uname -a > .tmpfile
	@grep CRAY .tmpfile ; \
	if [ $$? = 0 ] ; then echo "Compiling for Cray"							; \
		echo "AR		=	$(AR)"				>  macros_interpb	; \
		echo "RM		=	$(RM)"				>> macros_interpb	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_interpb	; \
		echo "CD		=	$(CD)"				>> macros_interpb	; \
		echo "LN		=	$(LN)"				>> macros_interpb	; \
		echo "MAKE		=	$(MAKE)"			>> macros_interpb	; \
		echo "SHELL		=	/bin/sh"			>> macros_interpb	; \
		echo "TOUCH		=	touch"				>> macros_interpb	; \
		echo "CPP		=	cpp"				>> macros_interpb	; \
		echo "CPPFLAGS		=	-I. -C -P -DCRAY"		>> macros_interpb	; \
		echo "FC		=	f90"				>> macros_interpb	; \
		echo "FCFLAGS		=	-f free -x omp"			>> macros_interpb	; \
		echo "LDFLAGS		=	-Ca"				>> macros_interpb	; \
		echo "CCFLAGS		=	-DCRAY -I."			>> macros_interpb	; \
		echo "LOCAL_LIBRARIES	=	"				>> macros_interpb	; \
		( $(CD) src ; $(MAKE) all )								; \
	else grep OSF .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for Compaq"						; \
		echo "AR		=	$(AR)"				>  macros_interpb	; \
		echo "RM		=	$(RM)"				>> macros_interpb	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_interpb	; \
		echo "CD		=	$(CD)"				>> macros_interpb	; \
		echo "LN		=	$(LN)"				>> macros_interpb	; \
		echo "MAKE		=	$(MAKE)"			>> macros_interpb	; \
		echo "SHELL		=	/bin/sh"			>> macros_interpb	; \
		echo "TOUCH		=	touch"				>> macros_interpb	; \
		echo "CPP		=	/usr/bin/cpp"			>> macros_interpb	; \
		echo "CPPFLAGS		=	-I. -C -P -DDEC -DF95"		>> macros_interpb	; \
		echo "FC		=	f95"				>> macros_interpb	; \
		echo "FCFLAGS		=	-automatic -convert big_endian -fast -fpe -free -pipeline -O4 -std -tune ev56 -u" >> macros_interpb	; \
		echo "LDFLAGS		=	-fast -O4"			>> macros_interpb	; \
		echo "CCFLAGS		=	-DDEC -I."			>> macros_interpb	; \
		echo "LOCAL_LIBRARIES	=	"				>> macros_interpb	; \
		( $(CD) src ; $(MAKE) all )							; \
	else grep HP .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for HP"							; \
		echo "AR		=	$(AR)"				>  macros_interpb	; \
		echo "RM		=	$(RM)"				>> macros_interpb	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_interpb	; \
		echo "CD		=	$(CD)"				>> macros_interpb	; \
		echo "LN		=	$(LN)"				>> macros_interpb	; \
		echo "MAKE		=	$(MAKE)"			>> macros_interpb	; \
		echo "SHELL		=	/bin/sh"			>> macros_interpb	; \
		echo "TOUCH		=	touch"				>> macros_interpb	; \
		echo "CPP		=	/opt/langtools/lbin/cpp"	>> macros_interpb	; \
		echo "CPPFLAGS		=	-I. -C -P -DHP"			>> macros_interpb	; \
		echo "FC		=	f90"				>> macros_interpb	; \
		echo "FCFLAGS		=	+langlvl=90 +source=free"	>> macros_interpb	; \
		echo "LDFLAGS		=	" 				>> macros_interpb	; \
		echo "CCFLAGS		=	-DHP -I."			>> macros_interpb	; \
		echo "LOCAL_LIBRARIES	=	"				>> macros_interpb	; \
		( $(CD) src ; $(MAKE) all )							; \
	else grep AIX .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for IBM"							; \
		echo "AR		=	$(AR)"				>  macros_interpb	; \
		echo "RM		=	$(RM)"				>> macros_interpb	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_interpb	; \
		echo "CD		=	$(CD)"				>> macros_interpb	; \
		echo "LN		=	$(LN)"				>> macros_interpb	; \
		echo "MAKE		=	$(MAKE)"			>> macros_interpb	; \
		echo "SHELL		=	/bin/sh"			>> macros_interpb	; \
		echo "TOUCH		=	touch"				>> macros_interpb	; \
		echo "CPP		=	/usr/lib/cpp"			>> macros_interpb	; \
		echo "CPPFLAGS		=	-I. -C -P -DIBM"		>> macros_interpb	; \
		echo "FC		=	xlf90"				>> macros_interpb	; \
		echo "FCFLAGS		=	-qfree=f90 -qlanglvl=90pure -O3 -qarch=auto -qnosave -qmaxmem=-1 -q64 -qspillsize=20000 -Q"	>> macros_interpb	; \
		echo "LDFLAGS		=	-O3 -q64" 			>> macros_interpb	; \
		echo "CCFLAGS		=	-DIBM -I."			>> macros_interpb	; \
		echo "LOCAL_LIBRARIES	=	"				>> macros_interpb	; \
		( $(CD) src ; $(MAKE) all )							; \
	else grep Darwin .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for MAC"							; \
		echo "AR		=	$(AR)"				>  macros_interpb	; \
		echo "RM		=	$(RM)"				>> macros_interpb	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_interpb	; \
		echo "CD		=	$(CD)"				>> macros_interpb	; \
		echo "LN		=	$(LN)"				>> macros_interpb	; \
		echo "MAKE		=	$(MAKE)"			>> macros_interpb	; \
		echo "SHELL		=	/bin/sh"			>> macros_interpb	; \
		echo "TOUCH		=	touch"				>> macros_interpb	; \
		echo "CPP		=	/usr/bin/cpp"			>> macros_interpb	; \
		echo "CPPFLAGS		=	-I. -C -P -DIBM -xassembler-with-cpp" >> macros_interpb	; \
		echo "FC		=	xlf90"				>> macros_interpb	; \
		echo "FCFLAGS		=	-qfree=f90 -qarch=auto "	>> macros_interpb	; \
		echo "LDFLAGS		=	-Wl,-stack_size,10000000,-stack_addr,0xc0000000" >> macros_interpb	; \
		echo "CCFLAGS		=	-DNOUNDERSCORE -DIBM -I."	>> macros_interpb	; \
		echo "LOCAL_LIBRARIES	=	"				>> macros_interpb	; \
		( $(CD) src ; $(MAKE) all )							; \
	else grep IRIX .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for SGI"							; \
		echo "AR		=	$(AR)"				>  macros_interpb	; \
		echo "RM		=	$(RM)"				>> macros_interpb	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_interpb	; \
		echo "CD		=	$(CD)"				>> macros_interpb	; \
		echo "LN		=	$(LN)"				>> macros_interpb	; \
		echo "MAKE		=	$(MAKE)"			>> macros_interpb	; \
		echo "SHELL		=	/bin/sh"			>> macros_interpb	; \
		echo "TOUCH		=	touch"				>> macros_interpb	; \
		echo "CPP		=	/lib/cpp"			>> macros_interpb	; \
		echo "CPPFLAGS		=	-I. -C -P -DSGI"		>> macros_interpb	; \
		echo "FC		=	f90"				>> macros_interpb	; \
		echo "FCFLAGS		=	-freeform -n32 -O2 -I."		>> macros_interpb	; \
		echo "LDFLAGS		=	-n32 -O2 -lfastm"		>> macros_interpb	; \
		echo "CCFLAGS		=	-DSGI -I. -n32"			>> macros_interpb	; \
		echo "LOCAL_LIBRARIES	=	"				>> macros_interpb	; \
		( $(CD) src ; $(MAKE) all )							; \
	else grep SUN .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for SUN"							; \
		echo "AR		=	$(AR)"				>  macros_interpb	; \
		echo "RM		=	$(RM)"				>> macros_interpb	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_interpb	; \
		echo "CD		=	$(CD)"				>> macros_interpb	; \
		echo "LN		=	$(LN)"				>> macros_interpb	; \
		echo "MAKE		=	$(MAKE)"			>> macros_interpb	; \
		echo "SHELL		=	/bin/sh"			>> macros_interpb	; \
		echo "TOUCH		=	touch"				>> macros_interpb	; \
		echo "CPP		=	/usr/ccs/lib/cpp"		>> macros_interpb	; \
		echo "CPPFLAGS		=	-I. -C -P -DSUN"		>> macros_interpb	; \
		echo "FC		=	f90"				>> macros_interpb	; \
		echo "FCFLAGS		=	-ansi -free"			>> macros_interpb	; \
		echo "LDFLAGS		=	" 				>> macros_interpb	; \
		echo "CCFLAGS		=	-DSUN -I."			>> macros_interpb	; \
		echo "LOCAL_LIBRARIES	=	"				>> macros_interpb	; \
		( $(CD) src ; $(MAKE) all	) 							; \
	else grep Linux .tmpfile									; \
	if [ $$? = 0 ] ; then echo "Compiling for Linux"						; \
		echo "AR		=	$(AR)"				>  macros_interpb	; \
		echo "RM		=	$(RM)"				>> macros_interpb	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_interpb	; \
		echo "CD		=	$(CD)"				>> macros_interpb	; \
		echo "LN		=	$(LN)"				>> macros_interpb	; \
		echo "MAKE		=	$(MAKE)"			>> macros_interpb	; \
		echo "SHELL		=	/bin/sh"			>> macros_interpb	; \
		echo "TOUCH		=	touch"				>> macros_interpb	; \
		echo "CPP		=	/lib/cpp"			>> macros_interpb	; \
		echo "CPPFLAGS		=	-I. -C -P -DDEC -traditional"	>> macros_interpb	; \
		echo "FC		=	pgf90"				>> macros_interpb	; \
		echo "FCFLAGS		=	-Mfreeform -pc 32 -byteswapio"	>> macros_interpb	; \
		echo "LDFLAGS		=	" 				>> macros_interpb	; \
		echo "CCFLAGS		=	-DDEC -I."			>> macros_interpb	; \
		echo "LOCAL_LIBRARIES	=	"				>> macros_interpb	; \
		( $(CD) src ; $(MAKE) all	) 							; \
	else grep UXP .tmpfile										; \
	if [ $$? = 0 ] ; then echo "Compiling for Fujitsu"						; \
		echo "AR		=	$(AR)"				>  macros_interpb	; \
		echo "RM		=	$(RM)"				>> macros_interpb	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_interpb	; \
		echo "CD		=	$(CD)"				>> macros_interpb	; \
		echo "LN		=	$(LN)"				>> macros_interpb	; \
		echo "MAKE		=	$(MAKE)"			>> macros_interpb	; \
		echo "SHELL		=	/bin/sh"			>> macros_interpb	; \
		echo "TOUCH		=	touch"				>> macros_interpb	; \
		echo "CPP		=	/lib/cpp"			>> macros_interpb	; \
		echo "CPPFLAGS		=	-I. -C -P -DVPP -DBIT32"	>> macros_interpb	; \
		echo "FC		=	frt"				>> macros_interpb	; \
		echo "FCFLAGS		=	-Free -X9 -Am -sc -Kfast -Kfreealloc -Karraystack3"	>> macros_interpb	; \
		echo "LDFLAGS		=	-J" 				>> macros_interpb	; \
		echo "CCFLAGS		=	-DVPP -DBIT32 -I."		>> macros_interpb	; \
		echo "LOCAL_LIBRARIES	=	"				>> macros_interpb	; \
		( $(CD) src ; $(MAKE) all	) 							; \
	else echo "Do not know how to compile for the `cat .tmpfile` machine." 				; \
	fi ; \
	fi ; \
	fi ; \
	fi ; \
	fi ; \
	fi ; \
	fi ; \
	fi ; \
	fi
	( $(RM) interpb ; $(LN) src/interpb . )

intel:
	echo "Compiling for Linux with INTEL compiler"						; \
		echo "AR		=	$(AR)"				>  macros_interpb	; \
		echo "RM		=	$(RM)"				>> macros_interpb	; \
		echo "RM_LIST		=	$(RM_LIST)"			>> macros_interpb	; \
		echo "CD		=	$(CD)"				>> macros_interpb	; \
		echo "LN		=	$(LN)"				>> macros_interpb	; \
		echo "MAKE		=	$(MAKE)"			>> macros_interpb	; \
		echo "SHELL		=	/bin/sh"			>> macros_interpb	; \
		echo "TOUCH		=	touch"				>> macros_interpb	; \
		echo "CPP		=	/lib/cpp"			>> macros_interpb	; \
		echo "CPPFLAGS		=	-I. -C -P -DDEC -traditional"	>> macros_interpb	; \
		echo "FC		=	ifort"				>> macros_interpb	; \
		echo "FCFLAGS		=	-FR -pc 32 -convert big_endian"	>> macros_interpb	; \
		echo "LDFLAGS		=	" 				>> macros_interpb	; \
		echo "CCFLAGS		=	-DDEC -I."			>> macros_interpb	; \
		echo "LOCAL_LIBRARIES	=	"				>> macros_interpb	; \
		( $(CD) src ; $(MAKE) all	) 							; \
	( $(RM) interpb ; $(LN) src/interpb . )

clean:
	( $(CD) src ; $(MAKE) clean "CD = $(CD)" "RM = $(RM)" "RM_LIST = $(RM_LIST)" )
	$(RM) $(RM_LIST) .tmpfile
	$(RM) FILE_MMOUTP* MMOUTP_DOMAIN* REGRID_DOMAIN*
