#	Lower level Makefile for interpb

include ../configure.rules
include ../macros_interpb

.IGNORE:

SRC		=	interpb.F

OBJS		=	$(SRC:.F=.o)

MODULES		=	module_all_io.o module_date_pack.o \
			module_diags.o \
			module_header_data.o module_interp.o \
			module_map_utils.o module_phys_consts.o \
			module_util.o module_xyll.o

LITTLE_F	=	$(SRC:.F=.f) $(MODULES:.o=.f)

all:			interpb \
			$(OBJS) $(MODULES) \
			../configure.rules ./Makefile ../Makefile

interpb:		../configure.rules ./Makefile ../Makefile \
			$(OBJS) $(MODULES)
			$(FC) -o $@ $(LDFLAGS) $(OBJS) $(MODULES)

interpb.o:		../configure.rules ./Makefile ../Makefile \
			module_all_io.o module_date_pack.o \
			module_diags.o \
			module_header_data.o module_interp.o module_util.o \
			module_xyll.o

module_all_io.o:	../configure.rules ./Makefile ../Makefile \
			module_date_pack.o module_header_data.o module_util.o

module_date_pack.o:	../configure.rules ./Makefile ../Makefile

module_diags.o:		../configure.rules ./Makefile ../Makefile \
			module_phys_consts.o module_util.o

module_header_data.o:	../configure.rules ./Makefile ../Makefile

module_interp.o:	../configure.rules ./Makefile ../Makefile \
			module_phys_consts.o

module_phys_consts.o:	../configure.rules ./Makefile ../Makefile 

module_util.o:		../configure.rules ./Makefile ../Makefile

module_xyll.o:		../configure.rules ./Makefile ../Makefile \
			module_map_utils.o

source:			$(LITTLE_F)
			$(CPP) $(CPPFLAGS) interpb.F > interpb.f
			$(CPP) $(CPPFLAGS) module_all_io.F > module_all_io.f
			$(CPP) $(CPPFLAGS) module_date_pack.F > module_date_pack.f
			$(CPP) $(CPPFLAGS) module_diags.F > module_diags.f
			$(CPP) $(CPPFLAGS) module_header_data.F > module_header_data.f
			$(CPP) $(CPPFLAGS) module_interp.F > module_interp.f
			$(CPP) $(CPPFLAGS) module_map_utils.F > module_map_utils.f
			$(CPP) $(CPPFLAGS) module_phys_consts.F > module_phys_consts.f
			$(CPP) $(CPPFLAGS) module_util.F > module_util.f
			$(CPP) $(CPPFLAGS) module_xyll.F > module_xyll.f

clean:
			$(RM) $(OBJS) $(MODULES) $(LITTLE_F) interpb core
			$(RM) *.mod
			$(RM) *.kmo
			$(RM) *.ctmi
