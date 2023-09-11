##################################
# System definitions
##################################
# include ../make.config
include ../defs/maverick
OBJ=${DIRECTORY}/genmesh/obj

####FFLAGS =$(DFFLAGS)  

#=========================================
all: doc tet2mesh patch2mesh patch2net partmesh restart2vtk restart2tec genmesh gmsh_mesh geomData average2tec

clean:
	/bin/rm -rf *.mod *.o obj $(DEST)/genmesh $(DEST)/tet2mesh $(DEST)/patch2mesh $(DEST)/partmesh $(DEST)/restart2vtk $(DEST)/restart2tec $(DEST)/geomData $(DEST)/average2tec

obj:
	mkdir   obj 

$(DEST):	
	mkdir   $(DEST)
	
info:
	@echo "Host  :" $(host) 
	@echo "Defs  :" $(defs) 
	@echo "DEST  :" $(DEST) 
	@echo "FC :" $(FC)  
	@echo "FFLAGS :" $(FFLAGS) 
	@echo "LIBS  :" $(LIBS) 

#=========================================
obj/%.o: %.f90 |obj
	$(FC) $(FFLAGS) $(FREE)  -c $*.f90 -o obj/$*.o

obj/%.o: %.f |obj
	$(FC) $(FFLAGS) $(FREE)  -c $*.f   -o obj/$*.o
      
#=========================================
genmesh: $(DEST)/genmesh 
	@ls -lh $(DEST)/genmesh 

gmsh_mesh: $(DEST)/gmsh_mesh 
	@ls -lh $(DEST)/gmsh_mesh 
	        
tet2mesh: $(DEST)/tet2mesh 
	@ls -lh $(DEST)/tet2mesh 
        
patch2mesh: $(DEST)/patch2mesh
	@ls -lh $(DEST)/patch2mesh 
        
patch2net: obj $(DEST)/patch2net
	@ls -lh $(DEST)/patch2net 
        
partmesh:   obj $(DEST)/partmesh 
	@ls -lh $(DEST)/partmesh 
        
restart2vtk: obj  $(DEST)/restart2vtk
	@ls -lh $(DEST)/restart2vtk

restart2tec: obj  $(DEST)/restart2tec
	@ls -lh $(DEST)/restart2tec
                       
geomData: obj  $(DEST)/geomData
	@ls -lh $(DEST)/geomData
          
average2tec: obj  $(DEST)/average2tec
	@ls -lh $(DEST)/average2tec             
#=========================================                
$(DEST)/genmesh: genmesh.sh | $(DEST)
	echo '#!/bin/bash'  > $(DEST)/genmesh
	echo DEST=$(DEST)   >> $(DEST)/genmesh
	echo METIS=$(METIS) >> $(DEST)/genmesh
	cat  genmesh.sh     >> $(DEST)/genmesh
	chmod +x $(DEST)/genmesh  

$(DEST)/gmsh_mesh: gmsh_mesh.sh | $(DEST)
	cp gmsh_mesh.sh     $(DEST)/gmsh_mesh
	chmod +x $(DEST)/gmsh_mesh  
	                           
$(DEST)/tet2mesh: obj/datatypes.o obj/mesh.o obj/input_FEM.o obj/tet2mesh.o obj/writeIEN.o obj/writeTEC.o obj/writeVTK.o  | $(DEST)
	$(FC) $(FFLAGS) $^ -o $(DEST)/tet2mesh

$(DEST)/patch2mesh: obj/datatypes.o obj/mesh.o obj/input_NURBS.o obj/mergeMeshes.o obj/removeDuplicateNodes.o obj/writeIEN.o obj/iso2hex.o obj/eval_shape_NURBS.o obj/writeVTK.o obj/patch2mesh.o | $(DEST)
	$(FC) $(FFLAGS) $^ -o $(DEST)/patch2mesh 
        
$(DEST)/patch2net: obj/datatypes.o obj/mesh.o obj/input_NURBS.o obj/mergeMeshes.o obj/patch2net.o obj/writeVTK.o | $(DEST)
	$(FC) $(FFLAGS) $^ -o $(DEST)/patch2net
        
$(DEST)/partmesh: obj/datatypes.o obj/mesh.o obj/partmesh.o | $(DEST)
	$(FC) $(FFLAGS) $^ -o $(DEST)/partmesh 

$(DEST)/restart2vtk: obj/datatypes.o obj/mesh.o obj/interpolate.o obj/eval_shape_NURBS.o obj/writeVTK.o obj/restart2vtk.o | $(DEST)
	$(FC) $(FFLAGS) $^ -o $(DEST)/restart2vtk 

$(DEST)/restart2tec: obj/datatypes.o obj/mesh.o obj/interpolate.o obj/eval_shape_NURBS.o obj/writeTEC.o obj/restart2tec.o | $(DEST)
	$(FC) $(FFLAGS) $^ -o $(DEST)/restart2tec

$(DEST)/geomData: obj/geomData.o | $(DEST)
	$(FC) $(FFLAGS) $^ -o $(DEST)/geomData 

$(DEST)/average2tec: obj/datatypes.o obj/mesh.o obj/writeTEC.o obj/average2tec.o | $(DEST)
	$(FC) $(FFLAGS) $^ -o $(DEST)/average2tec

doc:
	@echo "========================================================================="
	@echo "= "
	@echo "= This partition code handles internal surface elements"
	@echo "= as the slaves of the volume elements with bigger "
	@echo "= Element IDs. For example, if a triangle is shared "
	@echo "= by two tetrahedra with Element ID 101 and 102, respectively, "
	@echo "= F2E will save the element with ID 102 only."
	@echo "= "
	@echo "========================================================================="
.PHONY: doc
