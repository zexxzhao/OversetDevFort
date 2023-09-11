subroutine removeDuplBndNodes(mesh,eps)

  use class 
  implicit none
  
  type (mesh_class) :: mesh 
  real(8) eps, eps2
  
  integer :: b,f,i,j,n1,n2,nu,e,n
  real(8) :: x1(mesh%NSD),x2(mesh%NSD)
  integer :: N2N(mesh%NNODE)
  integer :: collapse,treshold 
  
  integer :: NFACE, NBNODE  
  
  logical found(mesh%NNODE)
  
  integer, allocatable :: FACE_IEN(:,:)
    
  integer, allocatable :: F2E(:),FACE_OR(:) 
  integer, allocatable :: BNODES(:)   

  eps2 = eps*eps
  
!
!... Find duplicate nodes in Boundaries
! 
  do n1 = 1, mesh%NNODE
    N2N(n1) = n1 
  enddo  
  if (mesh%NSHLB.eq.4) then
    treshold = 2
  else if (mesh%NSHLB.eq.9) then  
    treshold = 6
  endif
  
  do b = 1, mesh%NBOUND
    NFACE=0
    allocate(FACE_IEN(mesh%bound(b)%NFace,mesh%NSHLB))
    allocate(F2E(mesh%bound(b)%NFace))
    allocate(FACE_OR(mesh%bound(b)%NFace))

    do f = 1,mesh%bound(b)%NFace
      collapse = 0
      do i = 1,mesh%NSHLB     
         n1 = mesh%bound(b)%FACE_IEN(f,i)
         x1 = mesh%xg(n1,:)
         do j = i+1,mesh%NSHLB     
           n2 = mesh%bound(b)%FACE_IEN(f,j)    
           x2 = mesh%xg(n2,:) - x1
         
           if (sum(x2*x2).lt.eps2) then             
             N2N(max(n1,n2)) = N2N(min(n1,n2))
             collapse = collapse + 1
           endif
         enddo    
      enddo

      if (collapse.lt.treshold ) then
        NFACE = NFACE + 1

        FACE_IEN(NFACE,:) = mesh%bound(b)%FACE_IEN(f,:)
        F2E(NFACE)        = mesh%bound(b)%F2E(f)
        FACE_OR(NFACE)    = mesh%bound(b)%FACE_OR(f)
      endif                  
    enddo  
          
    deallocate(mesh%bound(b)%FACE_IEN,mesh%bound(b)%F2E,mesh%bound(b)%FACE_OR)
    
    allocate(mesh%bound(b)%FACE_IEN(NFace,mesh%NSHLB))
    allocate(mesh%bound(b)%F2E(NFace))
    allocate(mesh%bound(b)%FACE_OR(NFace))
    mesh%bound(b)%NFACE    = NFACE
    mesh%bound(b)%FACE_IEN = FACE_IEN(1:NFACE,:)
    mesh%bound(b)%F2E      = F2E(1:NFACE) 
    mesh%bound(b)%FACE_OR  = FACE_OR(1:NFACE)
    
    if (NFACE.eq.0) then
      mesh%bound(b)%FACE_ID = -1
      mesh%bound(b)%NNODE   =  0 
      write(*,*) "   Remove boundary:",b
    endif

    deallocate(FACE_IEN,F2E,FACE_OR)
  enddo


!... Find duplicate node map
!   
  nu = 0
  do n1 = 1, mesh%NNODE
  
    if (N2N(n1).eq.n1) then
      nu = nu + 1
      N2N(n1) = nu
      mesh%xg(nu,:) = mesh%xg(n1,:)
      mesh%wg(nu)   = mesh%wg(n1)
    else
      N2N(n1) = N2N(N2N(n1))
    endif     
  enddo  
 
  if (mesh%NNODE.gt.nu) write(*,*) "     Cleaning nodes:", mesh%NNODE,nu,mesh%NNODE-nu
  mesh%NNODE=nu
  
  
  
!
!... Apply duplicate node map
!
  do e = 1, mesh%NELEM
    do i = 1,mesh%NSHLmax 
      mesh%IEN(e,i) = N2N(mesh%IEN(e,i)) 
    enddo
  enddo
  
  do b = 1, mesh%NBOUND
    found = .false.
    do f = 1,mesh%bound(b)%NFace
      do i = 1,mesh%NSHLB 
        n = N2N(mesh%bound(b)%FACE_IEN(f,i))  
        mesh%bound(b)%FACE_IEN(f,i) = n 
        found (n) = .true.
      enddo
    enddo
    
    
    allocate(BNODES(mesh%bound(b)%NNODE))
    NBNODE = 0
    do n = 1, mesh%NNODE
      if (found (n)) then      
        NBNODE = NBNODE +1
        BNODES(NBNODE) =  n
      endif
    enddo
   
    deallocate(mesh%bound(b)%BNODES)
    allocate(mesh%bound(b)%BNODES(NBNODE))
    mesh%bound(b)%BNODES = BNODES(1:NBNODE)
    mesh%bound(b)%NNODE  = NBNODE  
    deallocate(BNODES)
  enddo  
  
end subroutine 
