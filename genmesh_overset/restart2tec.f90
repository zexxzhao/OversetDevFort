!------------------------------------------------------------------------
! Main routine to call all the subroutines
! rflag : plot on reference config and relative vel
!------------------------------------------------------------------------
program restart2tec
  
  use class
  
  implicit none
      
  type (mesh_class) :: mesh
  type (mesh_class) :: imesh
  type (map_class), allocatable :: L2GNODE(:) 
  integer :: i, j, k
  integer :: NProc, nstart, nstep, nskip, interp, istep, UMF, solf
  integer :: pfile, itmp, intep, ii, ln, gn, bn, rflag

  integer, allocatable :: PNNODE(:)
  real(8), allocatable ::  dg(:,:),  ug(:,:),  ugm(:,:),  pg(:), Tg(:),  phig(:) 
  real(8), allocatable :: rTg(:), rphig(:), acg(:, :)
  real(8), allocatable :: idg(:,:), iug(:,:), iugm(:,:), ipg(:), iphig(:), irphig(:)
  real(8), allocatable :: acgm(:, :), rhog(:), rho0g(:), mug(:), heg(:), hepg(:)
  real(8), allocatable :: var3d(:, :), var1d(:)
  real(8) :: rtmp, theta
  real(8) :: time
  integer :: Rstep
  logical :: fexist
  character(len=10) :: cname
  character(len=30) :: fname
  character(len=10) :: inp

  ! get main input        	   
  i = iargc()
  if (i == 2) then 
    call getarg(1,inp); read(inp,*) NProc
    call getarg(2,inp); read(inp,*) nstart 
    nstep  = nstart 
    nskip  = 1
    interp = -1
    rflag  = 0
  else if (i == 3) then 
    call getarg(1,inp); read(inp,*) NProc
    call getarg(2,inp); read(inp,*) nstart
    call getarg(3,inp); read(inp,*) interp  
    nstep = nstart 
    nskip = 1   
    rflag = 0    
  else if (i == 4) then    
    call getarg(1,inp); read(inp,*) NProc
    call getarg(2,inp); read(inp,*) nstart
    call getarg(3,inp); read(inp,*) nstep
    call getarg(4,inp); read(inp,*) nskip
    interp = -1
    rflag  = 0
  else if (i == 5) then    
    call getarg(1,inp); read(inp,*) NProc
    call getarg(2,inp); read(inp,*) nstart
    call getarg(3,inp); read(inp,*) nstep
    call getarg(4,inp); read(inp,*) nskip
    call getarg(5,inp); read(inp,*) interp 
    rflag = 0
  else if (i == 6) then    
    call getarg(1,inp); read(inp,*) NProc
    call getarg(2,inp); read(inp,*) nstart
    call getarg(3,inp); read(inp,*) nstep
    call getarg(4,inp); read(inp,*) nskip
    call getarg(5,inp); read(inp,*) interp 
    call getarg(6,inp); read(inp,*) rflag
  else    
    write(*,*) 'Usage 1: restart2tec NProc step'    
    write(*,*) 'Usage 2: restart2tec NProc step interp'
    write(*,*) 'Usage 3: restart2tec NProc start end skip '    
    write(*,*) 'Usage 4: restart2tec NProc step  end skip interp'
    write(*,*) 'Usage 5: restart2tec NProc step  end skip interp rot2ref'
    write(*,*) '   output: solution.*.tec'
    stop
  end if

  if (interp /= 0 .and. interp /= -1) then
    write(*,*) "!------------------------------------------------------------"
    write(*,*) "! Warning: only works for interp = 0 or -1, not ", interp
    write(*,*) "!------------------------------------------------------------"
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
    open(UMF, file = fname, status = 'unknown')    

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
  allocate(acg (mesh%NNODE,mesh%NSD))  
  allocate(pg  (mesh%NNODE))  
  allocate(phig(mesh%NNODE)) 
  allocate(rphig(mesh%NNODE)) 
  allocate(Tg  (mesh%NNODE))       
  allocate(rTg  (mesh%NNODE))       
  allocate(acgm(mesh%NNODE, mesh%NSD))
  ! allocate(rhog(mesh%NNODE))
  ! allocate(rho0g(mesh%NNODE))
  ! allocate(mug(mesh%NNODE))
  ! allocate(heg(mesh%NNODE))
  ! allocate(hepg(mesh%NNODE))
  allocate(var3d(mesh%NNODE,mesh%NSD))
  allocate(var1d(mesh%NNODE))


  ! Read interpolation mesh
  if (interp /= -1) then
    write(cname,'(I8)') interp
    fname = 'hex.' // trim(adjustl(cname)) //'.dat'
    write(*,*) "    Reading ", fname
    
    call readMesh(imesh,fname)  
 
    allocate(idg  (imesh%NNODE,imesh%NSD))  
    allocate(iug  (imesh%NNODE,imesh%NSD))  
    allocate(iugm (imesh%NNODE,imesh%NSD))  
    allocate(ipg  (imesh%NNODE))  
    allocate(iphig(imesh%NNODE))       
    allocate(irphig(imesh%NNODE))       
  end if

  ! Loop over time steps  
  do istep = nstart, nstep, nskip
    write(*,*) 'Step = ', istep
    
    ! Read solution for each proc and put it in global matrix
    do i = 1, NProc      
      solf = 12    
      write(cname,'(I8)') istep
      fname = 'restart.' // trim(adjustl(cname)) 
      write(cname,'(I8)') i
      fname = trim(fname) // '.' // trim(adjustl(cname))
      
      inquire(file=fname,exist=fexist)
      do while (.not.fexist) 
        inquire(file=fname,exist=fexist)
        call sleep(3)
      end do  
            
      open(solf, file=fname, status='old')
      
      ! read(solf) !Rstep
      ! read(solf) !time
      read(solf, *) Rstep, time
      ! read dg
      do k = 1, PNNODE(i)
        read(solf, *) (var3d(k, j), j=1,3)
      enddo
      do ln = 1, PNNODE(i)
        gn = L2GNODE(i)%node(ln,1)
        dg(gn, :) = var3d(ln, :)
      end do
       
      ! read ug
      do k = 1, PNNODE(i)
        read(solf, *) (var3d(k, j), j=1,3)
      enddo
      do ln = 1, PNNODE(i)
        gn = L2GNODE(i)%node(ln,1)
        ug(gn, :) = var3d(ln, :)
      end do

      ! read ugm
      do k = 1, PNNODE(i)
        read(solf, *) (var3d(k, j), j=1,3)
      enddo
      do ln = 1, PNNODE(i)
        gn = L2GNODE(i)%node(ln,1)
        ugm(gn, :) = var3d(ln, :)
      end do

      ! read acg 
      do k = 1, PNNODE(i)
        read(solf, *) (var3d(k, j), j=1,3)
      enddo
      do ln = 1, PNNODE(i)
        gn = L2GNODE(i)%node(ln,1)
        acg(gn, :) = var3d(ln, :)
      end do  

      ! read acgm 
      do k = 1, PNNODE(i)
        read(solf, *) (var3d(k, j), j=1,3)
      enddo
      do ln = 1, PNNODE(i)
        gn = L2GNODE(i)%node(ln,1)
        acgm(gn, :) = var3d(ln, :)
      end do

      ! read phi
      do k = 1, PNNODE(i)
        read(solf, *) var1d(k)
      enddo
      do ln = 1, PNNODE(i)
        gn = L2GNODE(i)%node(ln,1)
        phig(gn) = var1d(ln)
      end do

      ! read rphi

      do k = 1, PNNODE(i)
        read(solf, *) var1d(k)
      enddo
      do ln = 1, PNNODE(i)
        gn = L2GNODE(i)%node(ln,1)
        rphig(gn) = var1d(ln)
      end do
      
      ! read pg
      do k = 1, PNNODE(i)
        read(solf, *) var1d(k) 
      enddo
      do ln = 1, PNNODE(i) 
        gn = L2GNODE(i)%node(ln,1)         
        pg(gn) = var1d(ln)
      end do

      ! read Tg
      do k = 1, PNNODE(i)
        read(solf, *) var1d(k) 
      enddo
      do ln = 1, PNNODE(i) 
        gn = L2GNODE(i)%node(ln,1)         
        Tg(gn) = var1d(ln)
      end do

      ! read rTg
      do k = 1, PNNODE(i)
        read(solf, *) var1d(k) 
      enddo
      do ln = 1, PNNODE(i) 
        gn = L2GNODE(i)%node(ln,1)         
        rTg(gn) = var1d(ln)
      end do


      close(solf)
                  
    end do

    ! Interpolate solution if necessary
    if (interp /= -1) then
      call interpolate(mesh, dg, ug, ugm, pg, phig, interp, &
                       imesh, idg, iug, iugm, ipg, iphig)
    end if
                           
    ! Write global solution                    
    write(cname,'(I8)') istep
    !fname = 'solution.' // trim(adjustl(cname))
    !write(cname,'(I8)') rflag
    !fname = trim(fname)//'.'//trim(adjustl(cname))//'.tec'
    fname = 'solution.' // trim(adjustl(cname))//'.tec'


    if (interp == -1) then 
      call writeTEC(mesh, fname, ug, pg, Tg, phig, rphig, rflag)
    else
      call writeTEC(imesh, fname, iug, ipg, Tg, iphig, irphig, rflag)
    end if

    do bn = 1, mesh%NBOUND
      write(cname,'(I8)') istep
      fname = 'sol.face.'// trim(adjustl(cname))
      !write(cname,'(I8)') rflag
      !fname = trim(fname)//'.'//trim(adjustl(cname))
      write(cname,'(I8)') mesh%bound(bn)%FACE_ID
      fname = trim(fname)//'.'//trim(adjustl(cname))//'.tec'

      if (interp == -1) then
        call writeBndTEC(mesh,fname,bn,ug,pg,Tg,rflag)
      else
          !!!call writeBndVTK(imesh,fname,ibn,idg)
      end if
    end do

  end do  
               
end program restart2tec
