!------------------------------------------------------------------------
!                                                                        
!        Main routine to call all the subroutines                        
!                                                                        
!------------------------------------------------------------------------
program patch2net

  use class 
  implicit none

  character*10 inp
  type (mesh_class), allocatable ::  mesh(:)
  type (mesh_class) :: finalmesh 

   
  integer :: i,mn,NMESH,lfi
  character*30 :: listname
  character*30 :: patchname

  !==============================
  ! read in the run arguments
  !==============================
  i = iargc()	       
  if (i.eq.1) then 
    call getarg(1,listname)
  else
    write(*,*) 'Usage: patch2net  <listfile> '
    write(*,*) '   output: net.dat'
    stop
  endif
  
  !==============================
  ! read in the patches
  !==============================
  lfi = 44
  open (lfi, file=trim(listname), status='old')
  read (lfi,*) NMESH
  
  allocate(mesh(NMESH))
  do mn = 1, NMESH
      read (lfi,*) patchname
      call input_NURBS(mesh(mn),trim(patchname), .true.)
  enddo
  
  !==============================
  ! merge patches into one mesh
  !==============================  
  if (NMESH.gt.1) then
    call mergeMeshes(mesh, NMESH, finalmesh) 
  else
    finalmesh = mesh(1)
  endif
   
  do mn = 1, NMESH 
    call deallocMesh(mesh(mn))
  enddo 
  
  !==============================
  ! write mesh 
  !==============================      
  call writeMesh   (finalmesh, "net.dat")  
  call writeMeshVTK(finalmesh, "net.vtk")  
  
end program patch2net

