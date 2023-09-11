!------------------------------------------------------------------------
!                                                                        
!        Main routine to call all the subroutines                        
!                                                                        
!------------------------------------------------------------------------
program patch2mesh

  use class 
  implicit none

  character*10 inp
  type (mesh_class), allocatable ::  mesh(:)  
  type (patch_con_class) connectivity
  type (map_class) , allocatable ::  nodemap(:)
  type (mesh_class) :: finalmesh

   
  integer :: i,j,mn,NMESH,lfi,chk
  character*30 :: listname
  character*30 :: patchname
  real(8) :: eps
  logical :: ex

  !==============================
  ! read in the run arguments
  !==============================
  i = iargc()	       
  if (i.eq.1) then 
    call getarg(1,listname)
  else
    write(*,*) 'Usage: patch2mesh  <listfile> '
    write(*,*) '   output: mesh.dat graph.dat'
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
    call input_NURBS(mesh(mn),trim(patchname),.false.)
    !!call removeDuplBndNodes(mesh(mn))
  enddo
  
  !==============================
  ! merge patches into one mesh
  !==============================  
  if (NMESH.gt.0) then  
  
    read (lfi,*) patchname, connectivity%eps
    call getPatchConn (patchname, NMESH, mesh, connectivity)
    eps = connectivity%eps*connectivity%dsize
    do mn = 1, NMESH
      call removeDuplBndNodes(mesh(mn),eps)
    enddo
  
      
    read (lfi,*) patchname, connectivity%treshold
    call getPatchBndConn(patchname, NMESH, mesh, connectivity)
   
    allocate(nodemap(NMESH))
    call findNodeMap(NMESH,mesh,connectivity,nodemap)      
    
    call mergeMeshes(NMESH,mesh,nodemap,finalmesh) 
   
    !==============================
    ! Re arrange boundaries 
    ! Assumes a box shape alligned with mesh
    ! Assumes x=(0,0,0) is in the interior
    !==============================     
    call rearrangeBoundaries(finalmesh,eps)
  else
    call removeDuplBndNodes(mesh(1))
    finalmesh = mesh(1)
  endif 
  close(lfi)
   
  do mn = 1, NMESH 
    call deallocMesh(mesh(mn))
  enddo 
  

  !==============================
  ! write mesh and graph
  !==============================      
  call writeMesh(finalmesh, "mesh.dat")
        
  call iso2hex(finalMesh,eps)

end program patch2mesh 

