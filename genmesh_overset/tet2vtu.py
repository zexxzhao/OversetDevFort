#!/usr/bin/env pvpython

"""Demo script that converts a legacy VTK file to a VTU file.
"""

from paraview import servermanager
import os
import sys
import time
#import subprocess


#executable="/home/ido/bin/tet2vtk/tet2vtk"
#executable="/home/ido/tetmesh/bin/vizGlobSolLowMem"
executable="/share/home/01146/akkerman/tet2vtk/tet2vtk"


# Conversion routine
def convert(infile, outfile):
        print "Converting ", infile 
        
	# create a built-in connection
	if not servermanager.ActiveConnection:
	    connection = servermanager.Connect()

	# create reader for legacy VTK files
	reader = servermanager.sources.LegacyVTKFileReader(FileNames=infile)

	# create VTU writer and connect it to the reader
	writer = servermanager.writers.XMLUnstructuredGridWriter(Input=reader,
                                                         	 DataMode=1,
                                                         	 FileName=outfile)

	# Trigger execution of pipeline
	writer.UpdatePipeline()


# Read input

proc=len([name for name in os.listdir('.') if "mesh." in name and ".dat" in name])

print "Proc = ", proc

print len(sys.argv)
print     sys.argv

if len(sys.argv)==4:
	begin=int(sys.argv[1])
	end=int(sys.argv[2])+1
	stride=int(sys.argv[3])
        ostride=stride
        rank=1
elif len(sys.argv)==3:
	begin=int(sys.argv[1])
	end=int(sys.argv[2])+1
	stride=1
        ostride=stride
        rank=1
elif len(sys.argv)==2:
	begin=int(sys.argv[1])
	end=begin+1
	stride=1 
        ostride=stride
        rank=1
elif len(sys.argv)==5:
        print "In 5 arg if"
        #visprocs=16
        #rank=int(sys.argv[1])
        visprocs=int(sys.argv[1])
        begin=int(sys.argv[2])
        end=int(sys.argv[3])+1
        ostride=int(sys.argv[4])
        stride=visprocs*ostride
       
        print visprocs,begin,end,ostride,stride  
        #print rank
        rank=-1
        rank=int(os.getenv('MPIRUN_RANK',str(rank)))
        rank=int(os.getenv('OMPI_COMM_WORLD_RANK',str(rank)))
        rank=int(os.getenv('OMPI_MCA_ns_nds_vpid',str(rank)))
        rank=int(os.getenv('PMI_ID',str(rank)))

        print rank,"  /  ",  visprocs      

        begin= begin + rank*ostride
else:
	sys.stderr.write('Usage: begin end(default=begin) stride(default=1)')
        sys.exit(1)               

null = open("/dev/null","w") 

# Loop over files
for idx in range(begin,end,stride):
	infile="solution."+str(idx)+".vtk"
	outfile="solution."+str(idx)+".vtu"


# Check for restart files
        for p in range(1,proc+1):
       	   while not os.path.isfile("restart."+str(idx)+"."+str(p)):       
              print  "Waiting for restart."+str(idx)+"."+str(p)
              time.sleep(5)
              
# Convert in legacy vtu file (in ASCII with vtk extension)        
        #subprocess.call(executable + " "+str(proc) +" "+str(idx), shell=True)#,stdout=null)
        os.system(executable + " "+str(proc) +" "+str(idx))
# Convert to binary
        if os.path.isfile(infile):     
           convert(infile, outfile)
        else:
           print "ERROR: ",infile, "does not exist" 

# Remove ACII vtk file if conversion succeeded          
        if os.path.isfile(outfile):     
           os.remove(infile)


# Check if all files are generated 
for idx in range(begin,end,ostride):
        outfile="solution."+str(idx)+".vtu"
        while not os.path.isfile(outfile):
              time.sleep(5)






