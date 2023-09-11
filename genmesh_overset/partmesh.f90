!============================================================
! program to partition mesh and generate lwork
!============================================================
program partmesh

  use class
  implicit none
  
  character(len=10) :: inp

  type (mesh_class) :: mesh  
  type (mesh_class), allocatable :: part(:)
  integer :: i, j, k, ier, UMF, NProc
  integer :: e, le, ge, n, n2, ln, gn, p1, p2, bn, fn
  integer, allocatable :: IEP(:)
  character(len=10) :: cname
  character(len=30) :: fname
  character(len=30) :: meshname
  character(len=30) :: partname
    
  type (map_class), allocatable :: L2GELEM(:)
  type (map_class), allocatable :: L2GNODE(:)  

  logical, allocatable :: flag(:)
  integer, allocatable :: dupl(:)
  
  integer :: Noverlap
  integer, allocatable ::ON2N(:),N2ON(:)
  integer, allocatable ::overlap(:,:) 
  integer, allocatable :: NSnodes(:,:)
  integer, allocatable :: task(:,:)
  integer :: ntask
  logical :: found
    
  type (map_class), allocatable :: senddata(:,:)  

  integer, allocatable :: G2LNODES(:), G2LELEMS(:), G2SNODES(:,:), &
                          G2SELEMS(:,:)
  
  !==============================
  ! read in the run arguments
  !==============================
  i = iargc()	       
  if (i.eq.2) then 
    call getarg(1,meshname)
    call getarg(2,partname)
  else
    write(*,*) 'Usage: partmesh  <meshfile> <partitionfile>'
    write(*,*) '   output: part.*.dat l2g.*.dat lwork.*.dat'
    stop
  endif

  !==============================
  ! read in the input mesh file
  !==============================
  call readMesh(mesh, meshname)

  !========================================
  ! read in the domain decomposition array
  !========================================
  allocate(IEP(mesh%NElem))

  UMF = 12
  write(*,*) "Reading in: ", partname 
  open(UMF, file = trim(partname), status = 'old')  
  do ge = 1, mesh%NElem
    read(UMF,*) IEP(ge)
  end do
  IEP = IEP + 1  ! proc ID started from 0, change to 1

  close(UMF)  

  NProc = maxval(IEP)
  
  !==============================
  ! Init final meshes
  !==============================
  allocate(part(NProc))
  do i = 1, NProc
    part(i)%NSD     = mesh%NSD
    part(i)%NSHLmax = mesh%NSHLmax
    part(i)%NSHLB   = mesh%NSHLB
    part(i)%NNODE   = 0
    part(i)%NELEM   = 0
    part(i)%NBOUND  = mesh%NBOUND
    part(i)%NPATCH  = mesh%NPATCH
  end do
  

  !========================================
  write(*,*) "Assigning elements to proc"
  !========================================
  ! get the number of total local elements
  do ge = 1, mesh%NElem
    i = IEP(ge)
    part(i)%NELEM = part(i)%NELEM + 1      
  end do
  
  allocate(L2GELEM(NProc))
  do i = 1, NProc
    write(*,*) "  Proc ", i, " ==> ", part(i)%NELEM 
    allocate(L2GELEM(i)%node(part(i)%NELEM,1))
    part(i)%NELEM = 0
  end do
  
  do ge = 1, mesh%NElem
    i = IEP(ge)
    part(i)%NELEM = part(i)%NELEM + 1 
    L2GELEM(i)%node(part(i)%NELEM,1) = ge      
  end do
    
  !========================================
  write(*,*) "Assigning nodes to proc"
  ! Keep track of duplicate (shared) nodes 
  !========================================  
  allocate(flag(mesh%NNODE),dupl(mesh%NNODE))
  allocate(L2GNODE(NProc))
  dupl = 0
  do i = 1, NProc
    flag = .false.
    do le = 1, part(i)%NElem
      ge = L2GELEM(i)%node(le,1) 
      do k = 1, mesh%NSHL(ge)
        gn = mesh%IEN(ge,k)       
        flag(gn) = .true.
      end do
    end do
    
    part(i)%NNODE = 0
    do gn = 1, mesh%NNODE
      if (flag(gn)) then
        dupl(gn) = dupl(gn) + 1
        part(i)%NNODE = part(i)%NNODE + 1
      endif  
    end do    
    
    write(*,*) "  Proc ", i, " ==> ", part(i)%NNODE
     
    allocate(L2GNODE(i)%node(part(i)%NNODE,1))
    part(i)%NNODE = 0
    do gn = 1, mesh%NNODE
      if (flag(gn)) then
        part(i)%NNODE = part(i)%NNODE + 1
        L2GNODE(i)%node(part(i)%NNODE,1) = gn
      end if
    end do    
    
    ! Write local 2 global node map to file for postprocessing
    write(cname,'(I8)') i
    fname = 'l2g.' // trim(adjustl(cname)) //'.dat'
    write(*,*) "    Writting ", fname
    open(UMF, file = fname, status = 'replace')    

    write(UMF,*) part(i)%NNODE    
    do ln = 1, part(i)%NNODE 
      write(UMF,*) L2GNODE(i)%node(ln,1) 
    end do   
    
    close(UMF)
    
  end do 
  
  !===================================================
  ! Translate duplicate nodes to more usefull format 
  write(*,*) "Find communication nodes"
  !===================================================    
  allocate(N2ON(mesh%NNODE))
  Noverlap = 0
  N2ON = -1
  do gn = 1, mesh%NNODE
    if(dupl(gn).ge.2) then
      Noverlap = Noverlap + 1    
      N2ON(gn) = Noverlap
    end if  
  end do  
  
  write(*,*) "Duplication info : ", Noverlap, MAXVAL(dupl)

  allocate(overlap(Noverlap,MAXVAL(dupl)))
  overlap  = -1  
  dupl     = 0
  do i = 1, NProc
    do ln = 1, part(i)%NNODE
      gn = L2GNODE(i)%node(ln,1)
      n2 = N2ON(gn)
      
      if (n2.ge.1) then
        dupl(gn) = dupl(gn) + 1
        overlap(n2,dupl(gn)) = i
      endif  
    enddo   
  enddo  

  
  !==============================
  ! Generate communication data 
  ! in terms of global nodes
  !==============================
  allocate(NSnodes(NProc,NProc))
  NSnodes = 0
  do n = 1, mesh%NNODE
    j = N2ON(n)
    if (j.ge.1) then
      p1 = overlap(j,1)
      do k = 2, dupl(n)
        p2 = overlap(j,k)
        NSnodes(p1,p2) = NSnodes(p1,p2) + 1
        NSnodes(p2,p1) = NSnodes(p2,p1) + 1       
      enddo  
    endif
  enddo
 
  allocate(senddata(NProc,NProc))
  do i = 1, NProc
    do j = 1, NProc
      allocate(senddata(i,j)%node(NSnodes(i,j),1))
    enddo
  enddo
    
  NSnodes = 0
  do n = 1, mesh%NNODE
    j = N2ON(n)
    if (j.ge.1) then
      p1 = overlap(j,1)
      do k = 2, dupl(n)
        p2 = overlap(j,k)
        NSnodes(p1,p2) = NSnodes(p1,p2) + 1
        NSnodes(p2,p1) = NSnodes(p2,p1) + 1 
        
        senddata(p1,p2)%node(NSnodes(p1,p2),1) = n
        senddata(p2,p1)%node(NSnodes(p2,p1),1) = n  
      enddo  
    endif
  enddo  
  
  !===========================================
  write(*,*) "Create communication pattern" 
  !===========================================  
 
  ntask = 0
  do p1 = 1, NProc
    do p2 = p1+1, NProc
      if (NSnodes(p1,p2).gt.0) then      
        ntask = ntask +1    
      endif 
    enddo
  enddo  
  write(*,*) "Maximum number of tasks : ",ntask
  
  
  allocate(task(ntask,NProc))
  task = -1  
  ntask = 0
  
  do p1 = 1, NProc
    do p2 = p1+1, NProc
      if (NSnodes(p1,p2).gt.0) then 

        ! Look for empty slots in task list     
        found = .false.    
        do i = 1, ntask
          if ((task(i,p1).eq.-1).and.(task(i,p2).eq.-1)) then
            task(i,p1) = p2
            task(i,p2) = p1
            found = .true.
            exit    
          endif
        enddo

        ! Otherwise add to new task    
        if (.not.found) then      
          ntask = ntask +1
          task(ntask,p1) = p2
          task(ntask,p2) = p1       
        endif 
      endif 
      
    enddo
  enddo 
  
  !=================================    
  write(*,*) "Actual number of tasks : ",ntask
  if (NProc.le.1116) then
    do i = 1, ntask
      write(*,'(256I8)') task(i,:)
    enddo
  endif
  
  !=================================
  write(*,*) "Create final meshes" 
  ! write out immediatly to save mem
  !=================================
  allocate(G2LNODES(mesh%NNODE))
  allocate(G2LELEMS(mesh%NELEM))  
  do i = 1, NProc

    !... Get global 2 local node map 
    G2LNODES = -1
    do ln = 1, part(i)%NNODE
      gn = L2GNODE(i)%node(ln,1)
      G2LNODES(gn) = ln
    end do  
     
    !... Get global 2 local elem map
    G2LELEMS = -1
    do le = 1, part(i)%NELEM
      ge = L2GELEM(i)%node(le,1)
      G2LELEMS(ge) = le
    end do

       
    !... Get local nodes
    allocate(part(i)%xg(part(i)%NNODE,mesh%NSD), &
             part(i)%NID(part(i)%NNODE))
    do le = 1, part(i)%NNODE
      ge = L2GNODE(i)%node(le,1)
      part(i)%xg(le,:) = mesh%xg(ge,:)
      part(i)%NID(le)  = mesh%NID(ge)
    end do  
     
    !... Generate generate local IEN
    allocate(part(i)%IEN(part(i)%NELEM,mesh%NSHLmax), &
             part(i)%NSHL(part(i)%NELEM), &
             part(i)%EID(part(i)%NELEM))
    do le = 1, part(i)%NELEM
      ge = L2GELEM(i)%node(le,1)
      part(i)%NSHL(le) = mesh%NSHL(ge)
      part(i)%EID(le)  = mesh%EID(ge)
      do j = 1, mesh%NSHL(ge)  
        gn = mesh%IEN(ge,j)
        ln = G2LNODES(gn)         
        part(i)%IEN(le,j)  = ln
      end do        
    end do    
    
    !... Generate local wg/EPID/EIJK
    if (mesh%NPATCH.gt.0) then 
      allocate(part(i)%EPID(part(i)%NELEM),part(i)%EIJK(part(i)%NELEM,mesh%NSD))
      do le = 1, part(i)%NELEM
        ge = L2GELEM(i)%node(le,1)

        part(i)%EPID(le)   = mesh%EPID(ge)
        part(i)%EIJK(le,:) = mesh%EIJK(ge,:)       
      end do
      
      allocate(part(i)%wg(part(i)%NNODE))
      do le = 1, part(i)%NNODE
        ge = L2GNODE(i)%node(le,1)
        part(i)%wg(le)   = mesh%wg(ge)
      end do  
    endif
            
    !... Copy patches 
    allocate(part(i)%patch(mesh%NPATCH)) 
    do j = 1, mesh%NPATCH
      part(i)%patch(j) = mesh%patch(j)
    end do 
      
    !... Copy boundaries, translate elements and nodes to local numbering 
    allocate(part(i)%bound(mesh%NBOUND)) 
    do bn = 1, mesh%NBOUND   
      part(i)%bound(bn)%FACE_ID = mesh%bound(bn)%FACE_ID   
      allocate(part(i)%bound(bn)%FACE_IEN(mesh%bound(bn)%NFACE,mesh%NSHLB))  
      allocate(part(i)%bound(bn)%F2E(mesh%bound(bn)%NFACE))
      allocate(part(i)%bound(bn)%FACE_OR(mesh%bound(bn)%NFACE))  
      part(i)%bound(bn)%NFACE = 0 
      do fn = 1, mesh%bound(bn)%NFACE       
        ge = mesh%bound(bn)%F2E(fn) 
        le = G2LELEMS(ge) 
        if (le.gt.0) then
          part(i)%bound(bn)%NFACE = part(i)%bound(bn)%NFACE + 1
          part(i)%bound(bn)%F2E(part(i)%bound(bn)%NFACE) = le
          part(i)%bound(bn)%FACE_OR(part(i)%bound(bn)%NFACE) = mesh%bound(bn)%FACE_OR(fn)
          do j = 1, mesh%NSHLB 
            gn = mesh%bound(bn)%FACE_IEN(fn,j)
            part(i)%bound(bn)%FACE_IEN(part(i)%bound(bn)%NFACE,j) = G2LNODES(gn) 
          end do
        end if
      end do
      
      allocate(part(i)%bound(bn)%BNODES(mesh%bound(bn)%NNODE))   
      part(i)%bound(bn)%NNODE = 0
      do j = 1, mesh%bound(bn)%NNODE      
        gn = mesh%bound(bn)%BNODES(j)
        ln = G2LNODES(gn)
        if (ln.gt.0) then
          part(i)%bound(bn)%NNODE = part(i)%bound(bn)%NNODE + 1
          part(i)%bound(bn)%BNODES(part(i)%bound(bn)%NNODE) = ln
        end if  
      end do
    end do 
   
    !... Write mesh    
    write(cname,'(I8)') i
    fname = 'part.' // trim(adjustl(cname)) //'.dat'
    write(*,*) "Writting ", fname
    call writeMesh(part(i), fname)

    !... Write lwork
    !... Translate senddata to local indices
    p1 = i
    write(cname,'(I8)') i
    fname = 'lwork.' // trim(adjustl(cname)) //'.dat'
    write(*,*) "Writting ", fname
    open(UMF, file = fname, status = 'replace')

    write(UMF,*) Ntask, 1 + 4*ntask + 2*sum(NSnodes(p1,:))
    do j = 1, Ntask
      p2 = task(j,p1) 
      if ( p2.eq.-1) then
        write(UMF,*) 100
        write(UMF,*) 0
        write(UMF,*) -1
        write(UMF,*) 0
      else
        write(UMF,*) 200               ! tag
        
        if (p1.gt.p2) then
          write(UMF,*) 0               ! send 
        else
          write(UMF,*) 1               ! recv
        endif
        
        write(UMF,*) p2                ! proc
        write(UMF,*) NSnodes(p1,p2)    ! number of msg's
        do n = 1, NSnodes(p1,p2)          
          gn = senddata(p1,p2)%node(n,1)
          ln = G2LNODES(gn) 
          write(UMF,*) ln              ! node
          write(UMF,*) 1
        end do
      endif
    enddo       
    close(UMF)
                 
  enddo
         

  !-------------------------------------------------------------------
  write(*,*) "Create partitioned boundary to shell node map" 
  !-------------------------------------------------------------------
  ! Build the mapping between partitioned boundary mesh node number
  ! and shell mesh node number.
  ! Shell mesh node is unpartitioned, and numbering starts from 1 for 
  ! each boundary.
  ! Partitioned boundary mesh: part(i)%bound(bn)%
  ! Unpartitioned  shell mesh: mesh%bound(bn)%
  ! To build l2s: Local boundary node -> Local volume node number ->
  ! Global volume node number -> Global boundary (shell) node number
  !-------------------------------------------------------------------
  allocate(G2SNODES(mesh%NBOUND,mesh%NNODE))
  G2SNODES = -1
  ! loop through the unpartitioned boundaries (shell meshes)
  do bn = 1, mesh%NBOUND
    ! loop through the shell node number
    do j = 1, mesh%bound(bn)%NNODE
      ! get the global unpartitioned volume mesh node number
      gn = mesh%bound(bn)%BNODES(j)
      G2SNODES(bn,gn) = j
    end do
  end do

  allocate(G2SELEMS(mesh%NBOUND,mesh%NELEM))
  G2SELEMS = -1
  ! loop through the unpartitioned boundaries (shell meshes)
  do bn = 1, mesh%NBOUND
    ! loop through the shell element number
    do j = 1, mesh%bound(bn)%NFACE
      ! get the global unpartitioned volume mesh element number
      gn = mesh%bound(bn)%F2E(j)
      G2SELEMS(bn,gn) = j
    end do
  end do


  ! loop through the partitions
  do i = 1, NProc

    write(cname,'(I8)') i
    fname = 'l2b.'// trim(adjustl(cname))//'.dat'
    write(*,*) "    Writting ", fname
    open(UMF, file = fname, status = 'replace')    

    ! for every partition, loop through all boundaries
    do bn = 1, mesh%NBOUND
      
      ! write out the boundary (face) ID, and how many nodes belong
      ! to this boundary. For example, if this partition does not
      ! have boundary number 3, then the node number should be 0.
      write(UMF,*) part(i)%bound(bn)%FACE_ID, part(i)%bound(bn)%NNODE, &
                   part(i)%bound(bn)%NFACE

      ! loop through the partitioned boundary node numbers
      do j = 1, part(i)%bound(bn)%NNODE
        ! for each node, get the node number corresponds to
        ! partitioned (local) volume mesh
        k = part(i)%bound(bn)%BNODES(j)

        ! check
        if (G2SNODES(bn,L2GNODE(i)%node(k,1)) <= 0) then
          write(*,*) "ERROR@Partition", i, ", Boundary", bn, ", Node", j
          stop
        end if

        ! output partitioned boundary node number and the corresponding
        ! boundary mesh node number
        ! From local volume mesh node number, you get global volume 
        ! mesh node number (L2GNODE), and based on that, you get
        ! the unpartitioned boundary (global, shell) node number
        write(UMF,*) k, G2SNODES(bn,L2GNODE(i)%node(k,1))
      end do

      ! loop through the partitioned boundary element numbers
      do j = 1, part(i)%bound(bn)%NFACE
        ! for each element, get the element number corresponds to
        ! partitioned (local) volume mesh
        k = part(i)%bound(bn)%F2E(j)

        ! check
        if (G2SELEMS(bn,L2GELEM(i)%node(k,1)) <= 0) then
          write(*,*) "ERROR@Partition", i, ", Boundary", bn, ", Face", j
          stop
        end if

        ! output partitioned boundary element number and the corresponding
        ! boundary mesh element number
        write(UMF,*) k, G2SELEMS(bn,L2GELEM(i)%node(k,1))
      end do

    end do
    
    close(UMF)
  end do

  deallocate(G2SNODES, G2SELEMS)

end program partmesh
