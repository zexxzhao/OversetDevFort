###NAME := nrel_rot
###NAME := 5mw_fsi
###NAME := 5mw_old
NAME := quenching
### NAME:=post

##################################
# System definitions
##################################
include ../defs/maverick

##################################
#  Files
##################################
include file.list
#F77FILES = $(MDFILES) $(FILES)
#OFILES   = $(F77FILES:.f=.o) 
#OBJ      = $(F77FILES:%.f=obj/%.o)
#DOBJ     = $(F77FILES:%.f=obj/%.dbg.o)

OFILES = $(MODFILES:.f90=.o)\
         $(F90FILES:.f90=.o)\
         $(F77FILES:.f=.o)

OBJ = $(MODFILES:%.f90=obj/%.o)\
      $(F90FILES:%.f90=obj/%.o)\
      $(F77FILES:%.f=obj/%.o)

DOBJ = $(MODFILES:%.f90=obj/%.dbg.o)\
       $(F90FILES:%.f90=obj/%.dbg.o)\
       $(F77FILES:%.f=obj/%.dbg.o)

##################################
# Rules
##################################
default: $(DEST)/$(NAME).exe 
debug: $(DEST)/$(NAME).dbg.exe 

info:
	@echo "Host  :" $(host) 
	@echo "Defs  :" $(defs) 
	@echo "NAME  :" $(NAME)
	@echo "DEST  :" $(DEST) 
	@echo "CC    :" $(CC)  
	@echo "CFLAGS:" $(CFLAGS) 
	@echo "FC    :" $(FC)  
	@echo "FFLAGS:" $(FFLAGS) 
	@echo "DFLAGS:" $(DFLAGS) 
	@echo "LIBS  :" $(LIBS) 	
	@echo "MDF   :" $(MDFILES) 

$(DEST)/$(NAME).exe:$(OBJ) | $(DEST)
	rm -rf $(DEST)/$(NAME).exe
	$(FC)   $(OBJ) $(LIBS) -o  $(DEST)/$(NAME).exe
	@ls -lh $(DEST)/$(NAME).exe

$(DEST)/$(NAME).dbg.exe:  $(DOBJ) | $(DEST)
	rm -rf $(DEST)/$(NAME).dbg.exe
	$(FC)   $(DOBJ) $(LIBS) $(DFLAGS) -o  $(DEST)/$(NAME).dbg.exe
	@ls -lh $(DEST)/$(NAME).dbg.exe

obj:
	mkdir   obj 

$(DEST):	
	mkdir   $(DEST)
	
clean:
	rm -rf obj *.mod $(DEST)/$(NAME).exe $(DEST)/$(NAME).dbg.exe
	
obj/%.o: %.f90 | obj
	$(FC) $(FFLAGS) -c $*.f90 -o obj/$*.o

obj/%.o: %.f | obj
	$(FC) $(FFLAGS) -c $*.f -o obj/$*.o

# debug mode
obj/%.dbg.o: %.f90 | obj
	$(FC) $(DFLAGS) -c $*.f90 -o obj/$*.dbg.o

obj/%.dbg.o: %.f | obj 
	$(FC) $(DFLAGS) -c $*.f  -o obj/$*.dbg.o
