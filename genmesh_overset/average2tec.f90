!------------------------------------------------------------------------
! Main routine to call all the subroutines
! rflag : plot on reference config and relative vel
!------------------------------------------------------------------------
program average2tec
  
  use class
  
  implicit none
      
  type (mesh_class) :: mesh
  type (mesh_class) :: imesh
  type (map_class), allocatable :: L2GNODE(:) 
  integer :: i, NProc, UMF, solf
  integer :: pfile, itmp, intep, ii, ln, gn, bn, rflag

  integer, allocatable :: PNNODE(:)
  real(8), allocatable ::  dg(:,:),  ug(:,:),  ugm(:,:),  pg(:),  phig(:)
  real(8) :: rtmp, theta
  logical :: fexist
  character(len=10) :: cname
  character(len=30) :: fname
  character(len=10) :: inp

  theta = 0.0d0
  rflag = 0

  ! get main input        	   
  i = iargc()
  if (i == 1) then 
    call getarg(1,inp); read(inp,*) NProc
  else    
    write(*,*) 'Usage 1: average2tec NProc'    
    write(*,*) 'output: solution.avg.tec'
    stop
  end if

  ! Read global 2 local nodes
  allocate(PNNODE(NProc), L2GNODE(NProc))
  
  do i = 1, NProc
    ! Write local 2 global node map to file for postprocessing
    write(cname,'(I8)') i
    UMF = 88
    fname = 'l2g.' // trim(adjustl(cname)) //'.dat'
    write(*,*) "    Reading ", fname
    open(UMF, file = fname, status = 'old')    

    read(UMF,*)PNNODE(i)   
    allocate(L2GNODE(i)%node(PNNODE(i),1))
    do ln = 1, PNNODE(i)
      read(UMF,*) L2GNODE(i)%node(ln,1) 
    end do   
    
    close(UMF)
  end do

  ! Read global mesh
  write(*,*) "    Reading ", "mesh.dat"
  call readMesh(mesh, "mesh.dat")  
 
  allocate(dg  (mesh%NNODE,mesh%NSD))  
  allocate(ug  (mesh%NNODE,mesh%NSD))
  allocate(ugm (mesh%NNODE,mesh%NSD))  
  allocate(pg  (mesh%NNODE))  
  allocate(phig(mesh%NNODE))       

  dg = 0.0d0; ug = 0.0d0; ugm = 0.0d0; 
  pg = 0.0d0; phig = 0.0d0

  ! Read solution for each proc and put it in global matrix
  do i = 1, NProc      
    solf = 12    
    write(cname,'(I8)') i
    fname = 'avgsol.' // trim(adjustl(cname))             
    open(solf, file=fname, status='old')
      
    read(solf,*) 

    do ln = 1, PNNODE(i)   
      gn = L2GNODE(i)%node(ln,1)    
      read(solf,*) ug(gn,:)
    end do 
                
    do ln = 1, PNNODE(i) 
      gn = L2GNODE(i)%node(ln,1)         
      read(solf,*) pg(gn)  
    end do

    close(solf)
                  
  end do
                           
  ! Write global solution                    
  fname = "solution.avg.tec"
  call writeTEC(mesh, fname, dg, ug, ugm, pg, phig, theta, rflag)

  do bn = 1, mesh%NBOUND
    write(cname,'(I8)') mesh%bound(bn)%FACE_ID
    fname = 'sol.face.avg.'//trim(adjustl(cname))//'.tec'   
         
    call writeBndTEC(mesh, fname, bn, dg, ug, ugm, pg, theta, rflag)
  end do               
end program average2tec
