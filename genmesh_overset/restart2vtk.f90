!------------------------------------------------------------------------
!        Main routine to call all the subroutines
!------------------------------------------------------------------------
program restart2vtk
  
  use class
  
  implicit none
      
  type (mesh_class)  mesh
  type (mesh_class) imesh
  type (map_class), allocatable :: L2GNODE(:) 
  integer :: i,NProc,nstart, nstep, nskip,interp,istep,UMF,sol_file
  integer :: pfile,itmp,intep,ii,ln,gn,bn 

  integer, allocatable :: PNNODE(:)
  real(8), allocatable :: dg(:,:),ug(:,:),ugm(:,:),pg(:),phig(:) 
  real(8), allocatable :: idg(:,:),iug(:,:),iugm(:,:),ipg(:),iphig(:)
  real(8) :: rtmp
  character(len=10) :: cname
  character(len=30) :: fname
  character(len=10) :: inp
  logical :: fexist

  ! get main input        	   
  i = iargc()
  if (i == 2) then 
    call getarg(1,inp); read(inp,*) NProc
    call getarg(2,inp); read(inp,*) nstart 
    nstep = nstart 
    nskip = 1
    interp = -1    
  else if (i == 3) then 
    call getarg(1,inp); read(inp,*) NProc
    call getarg(2,inp); read(inp,*) nstart
    call getarg(3,inp); read(inp,*) interp  
    nstep = nstart 
    nskip = 1       
  else if (i == 4) then    
    call getarg(1,inp); read(inp,*) NProc
    call getarg(2,inp); read(inp,*) nstart
    call getarg(3,inp); read(inp,*) nstep
    call getarg(4,inp); read(inp,*) nskip
    interp = -1
  else if (i == 5) then    
    call getarg(1,inp); read(inp,*) NProc
    call getarg(2,inp); read(inp,*) nstart
    call getarg(3,inp); read(inp,*) nstep
    call getarg(4,inp); read(inp,*) nskip
    call getarg(5,inp); read(inp,*) interp 
  else    
    write(*,*) 'Usage 1: restart2vtk NProc step '    
    write(*,*) 'Usage 2: restart2vtk NProc step interpolation_mode'
    write(*,*) 'Usage 3: restart2vtk NProc start end skip '    
    write(*,*) 'Usage 4: restart2vtk NProc step  end skip interpolation_mode'
    write(*,*) '   output: solution.*.vtk'
    stop
  end if
                 

  ! Read global mesh
  write(*,*) "    Reading ", "mesh.dat"
  call readMesh(mesh, "mesh.dat")  

  ! Read global 2 local nodes
  allocate(PNNODE(NProc), L2GNODE(NProc))

  if (NProc .eq.1) then  
    PNNODE(1) = mesh%NNode   
    allocate(L2GNODE(1)%node(PNNODE(1),1))
    do ln = 1, PNNODE(1)
      L2GNODE(1)%node(ln,1) = ln 
    enddo 
  else
  do i = 1, NProc
    ! Write local 2 global node map to file for postprocessing
    write(cname,'(I8)') i
    UMF = 88
    fname = 'l2g.' // trim(adjustl(cname)) //'.dat'
    write(*,*) "    Reading ", fname
    open(UMF, file = fname, status = 'unknown')    

    read(UMF,*)PNNODE(i)   
    allocate(L2GNODE(i)%node(PNNODE(i),1))
    do ln = 1, PNNODE(i)
      read(UMF,*)L2GNODE(i)%node(ln,1) 
    enddo   
    
    close(UMF)
  end do
  endif

  ! Allocate for global mesh 
  allocate(dg(mesh%NNODE,mesh%NSD))  
  allocate(ug(mesh%NNODE,mesh%NSD))  
  allocate(ugm(mesh%NNODE,mesh%NSD))  
  allocate(pg(mesh%NNODE))  
  allocate(phig(mesh%NNODE))       


  ! Read interpolation mesh
  if (interp /= -1) then
    write(cname,'(I8)') interp
    fname = 'hex.' // trim(adjustl(cname)) //'.dat'
    write(*,*) "    Reading ", fname
    
    call readMesh(imesh,fname)  
 
    allocate(idg(imesh%NNODE,imesh%NSD))  
    allocate(iug(imesh%NNODE,imesh%NSD))  
    allocate(iugm(imesh%NNODE,imesh%NSD)) 
    allocate(ipg(imesh%NNODE))  
    allocate(iphig(imesh%NNODE))       
  end if

  ! Loop over time steps  
  do istep = nstart, nstep, nskip
    write(*,*) 'Step = ', istep
    
    ! Read solution for each proc and puti in global matrix
    do i = 1, NProc      
      sol_file = 12    
      write(cname,'(I8)') istep
      fname = 'restart.' // trim(adjustl(cname)) 
      write(cname,'(I8)') i
      fname = trim(fname) // '.' // trim(adjustl(cname))
      
      inquire(file=fname,exist=fexist)
      do while (.not.fexist) 
        inquire(file=fname,exist=fexist)
        call sleep(3)
      enddo  
            
      open(sol_file, file=fname, status='old')
      
      read(sol_file,*)
      
      read(sol_file,*)
      read(sol_file,*)
      read(sol_file,*)
         
      do ln = 1, PNNODE(i)   
        gn = L2GNODE(i)%node(ln,1)    
        read(sol_file,*) dg(gn,:)
      end do 
        
      do ln = 1, PNNODE(i)   
        gn = L2GNODE(i)%node(ln,1)    
        read(sol_file,*) ug(gn,:)
      end do 
      
      do ln = 1, PNNODE(i)   
        read(sol_file,*) ugm(gn,:)
      end do  
      do ln = 1, PNNODE(i)   
        read(sol_file,*) rtmp, rtmp, rtmp ! acg
      end do 
      do ln = 1, PNNODE(i)   
        read(sol_file,*) rtmp, rtmp, rtmp ! acgm
      end do              
                                      
      do ln = 1, PNNODE(i) 
        gn = L2GNODE(i)%node(ln,1)         
        read(sol_file,*) phig(gn)  
      end do
      
      do ln = 1, PNNODE(i)       
        read(sol_file,*) rtmp ! rphig 
      end do  
          
      do ln = 1, PNNODE(i) 
        gn = L2GNODE(i)%node(ln,1)         
        read(sol_file,*) pg(gn)  
      end do

      close(sol_file)
                  
    end do

    ! Interpolate solution if necessary
    if (interp /= -1) then
!      call interpolate(mesh,dg,ug,ugm,pg,phig,interp,        &
!                       imesh,idg,iug,iugm,ipg,iphig)
    end if
                           
    ! Write global solution                    
    write(cname,'(I8)') istep
    fname = 'solution.' // trim(adjustl(cname)) //'.vtk' 
    if (interp == -1) then 
      call writeVTK(mesh, fname, dg, ug, pg, phig)
    else
      call writeVTK(imesh, fname, idg, iug, ipg, iphig)
    end if

    do bn = 1, mesh%NBOUND

      if (mesh%bound(bn)%FACE_ID == 7) then
        write(cname,'(I8)') istep
        fname = 'sol.face.' // trim(adjustl(cname))
        write(cname,'(I8)') bn
        fname = trim(fname) // '.' // trim(adjustl(cname))//'.vtk'
        if (interp == -1) then 
          call writeBndVTK(mesh, fname, bn, dg)
        else
          call writeBndVTK(imesh, fname, bn, idg)
        end if
      end if
    end do

  end do  
               
end program restart2vtk
