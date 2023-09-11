!======================================================================
! Main routine to call all the subroutines              
!======================================================================
program tet2mesh

  use class 
  implicit none

  type (mesh_class) :: mesh 
  integer :: i,b
  real(8), allocatable :: dg(:,:)
  character(len=10) :: inp
  character(len=30) :: fname 
  
  i = iargc()  
  if (i.eq.1) then 
    call getarg(1,fname)
  else
    write(*,*) 'Usage: tet2mesh  <meshfile> '
    write(*,*) '   output: mesh.dat ien.dat'
    stop
  end if

  call input_fem(mesh, fname)

  call writeMesh(mesh, "mesh.dat")  

  ! output surface mesh for blade
  do b = 1, mesh%NBOUND 
    write(*,*) b, mesh%bound(b)%Face_ID!, imesh%bound(b)%Face_OR(1)

    write(inp,'(I8)') mesh%bound(b)%Face_ID
    fname = 'bmesh.' // trim(adjustl(inp)) //'.dat'  
    call  writeBndMesh(mesh, fname, b)   
  end do

  call writeIEN(mesh, "ien.dat")

  call writeMeshVTK(mesh,'mesh.vtk')

  do b = 1, mesh%NBOUND 
    write(*,*) b, mesh%bound(b)%Face_ID!, imesh%bound(b)%Face_OR(1)

    write(inp,'(I8)') mesh%bound(b)%Face_ID
    fname = 'bmesh.' // trim(adjustl(inp)) //'.tec'  
    call writeBndMeshTEC(mesh, fname, b)   
  end do


!  if (mesh%nbound.ge.7) then
!    allocate(dg(mesh%NNODE,mesh%NSD))
!    dg = 0d0
!    call writeBndVTK(mesh,'obstacle.vtk',7,dg)
!  endif
  
!   do b = 1, mesh%NBOUND 
!     write(*,*) b, mesh%bound(b)%Face_ID!, imesh%bound(b)%Face_OR(1)

!     fname = 'bnd.' // trim(adjustl(inp))
!     write(inp,'(I8)') mesh%bound(b)%Face_ID
!     fname = 'bnd.' // trim(adjustl(inp)) //'.vtk'  
!     call  writeBndMeshVTK(mesh,fname,b)   
!   end do
      
end program tet2mesh 
