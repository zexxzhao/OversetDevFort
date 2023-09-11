
#===============================================
# Generate monolithic meshfile if necessary
#===============================================
if [ $# -eq 2 ] ; then
  if [  -f mesh.dat ]; then
    echo "Mesh generation skipped"
    echo " - mesh.dat already exists"    
  elif [[ "$1" == *msh ]]; then
    rm -rf lwork.*.dat l2g.*.dat part.*.dat ien.dat* graph.dat.* hex.*.dat mesh.dat >& /dev/null
    echo "Converting gmsh tetrahedral mesh"
    $DEST/tet2mesh $1
  elif [[ "$1" == *list ]]; then 
    rm -rf lwork.*.dat l2g.*.dat part.*.dat ien.dat* graph.dat.* hex.*.dat mesh.dat >& /dev/null
    echo "Merging NURBS patches mesh"
    $DEST/patch2mesh $1    
  else
    echo "Usage: genmesh <file> <nproc>"
    echo " - *.msh for tetrahedral mesh in gmsh "       
    echo " - *.list for NURBS patch list file "  
  fi
 
  shift
fi

#===============================================
# Serial run is an exception
#===============================================
if [ "$1" == 1 ] ; then
  rm -rf lwork.*.dat l2g.*.dat part.*.dat  >& /dev/null 
  cp mesh.dat part.1.dat    
  echo "Serial run, no partitioning performed"
  exit
fi

#===============================================
# Generate partition using metis
#===============================================
if [ -f graph.dat ]; then
  partfile="graph.dat.part."$1 
   if [ ! -f $partfile ]; then
      rm -rf lwork.*.dat l2g.*.dat part.*.dat $partfile >& /dev/null
      mpmetis graph.dat $1  
   else
     echo "Metis skipped" 
     echo " - $partfile already exists"
  fi
elif [ -f ien.dat ]; then
  partfile="ien.dat.epart."$1
  if [ ! -f $partfile ]; then      
    rm -rf lwork.*.dat l2g.*.dat part.*.dat $partfile >& /dev/null
    # partition into ($1 * 2) parts
    mpmetis ien.dat $((2*$1))
    # mod every entry in ien.dat.epart.$((2*$1)) by $i and output to ien.dat.epart.$i
    awk -v n=$1 '{print $1%n}' ien.dat.epart.$((2*$1)) > ien.dat.epart.$1
    # mod every entry in ien.dat.npart.$((2*$1)) by $i and output to ien.dat.npart.$i
    awk -v n=$1 '{print $1%n}' ien.dat.npart.$((2*$1)) > ien.dat.npart.$1
    rm -rf ien.dat.epart.$((2*$1)) ien.dat.npart.$((2*$1))
    # mpmetis ien.dat $1
  else
    echo "Metis skipped" 
    echo " - "$partfile" already exists"
  fi
else
  echo "No file for metis to partition"
  echo " - ien.dat and graph.dat absent"
  exit
fi

if [ ! -f $partfile ]; then  
    echo "No partitioning file"
    echo " - $partfile absent" 
    exit 
fi
    
#===============================================
# Generate mesh partition vtk file
#===============================================
if [ -f mesh.vtk ] && [ -f $partfile ]; then
  cp mesh.vtk  part.$1.vtk
  
  CELLS=`wc -l $partfile | awk '{print $1}'`
  echo CELL_DATA $CELLS  >> part.$1.vtk
  echo FIELD FieldData 1 >> part.$1.vtk
  echo part 1 $CELLS int >> part.$1.vtk
  cat  $partfile         >> part.$1.vtk       
elif [ -f mesh.0.vtk ] && [ -f $partfile ]; then
  cp mesh.0.vtk  part.$1.vtk
  
  CELLS=`wc -l $partfile | awk '{print $1}'`
  echo CELL_DATA $CELLS  >> part.$1.vtk
  echo FIELD FieldData 1 >> part.$1.vtk
  echo part 1 $CELLS int >> part.$1.vtk
  cat  $partfile         >> part.$1.vtk   
else
  echo "No files for VTK partition viz"
  echo " - mesh.vtk,mesh.0.vtk or $partfile absent"
  exit
fi 

#===============================================
# Partition mesh 
#===============================================  
if [ ! -f mesh.dat ]; then 
    echo "No meshfile to partition" 
    echo " - mesh.dat absent"    
elif [ ! -f $partfile ]; then  
    echo "No partitioning file"
    echo " - $partfile absent"    
elif ((  `ls  lwork.*.dat l2g.*.dat part.*.dat|wc -l` != $((3*$1)) )); then
    rm -rf l2g.*.dat lwork.*.dat part.*.dat     
    $DEST/partmesh mesh.dat $partfile   
else
    echo "Partitioning skipped"
    echo " - Files already exist" 
fi

