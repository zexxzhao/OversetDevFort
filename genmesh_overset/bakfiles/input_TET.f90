subroutine input_tet(mesh,fname)

  use class
      
  implicit none
  
  type (mesh_class) mesh
  character*(*) fname
                  
  integer :: i, j, ier,n,ii, meshf, itmp,ifac,iel,found
  integer :: flist(3), elist1(3), elist2(3), elist3(3)
  integer :: elist4(3)
     
  integer :: NFace,NFaceID, Nedge
  integer, allocatable ::  Face_IEN(:,:), FaceID(:)
  integer, allocatable ::  FID2BND(:)
  logical, allocatable ::  nflag(:)
  integer, allocatable ::  G2L(:)
  
  integer, allocatable :: INE(:,:)
  integer, allocatable :: cnt(:)
  
  real(8), allocatable :: dg(:,:)
          
  character*(20) :: ctmp
  mesh%NPATCH = 0
  mesh%NSHL   = 4
  mesh%NSHLb  = 3
  
  !==============================
  ! Read gmsh .mesh file
  !==============================
  meshf = 90
  open(meshf, file = fname, status = 'old')
  write(*,*) "Reading ", fname   
  
  read(meshf,*) ctmp, itmp  
  read(meshf,*) ctmp
  read(meshf,*) mesh%NSD
  write(*,*) "Dimension", mesh%NSD
  
  ! Vertices
  read(meshf,*) ctmp  
  read(meshf,*) mesh%NNode
  write(*,*) ctmp, " ", mesh%NNode
  allocate(mesh%xg(mesh%NNode,mesh%NSD))
  do i = 1, mesh%NNode
    read(meshf,*) ( mesh%xg(i,j), j = 1, mesh%NSD), itmp    
  end do      
  read(meshf,*) ctmp

  ! Edges
  if (ctmp.eq.'Edges') then  
    read(meshf,*) NEdge
    write(*,*) ctmp, " ", Nedge
    do i = 1, NEdge
      read(meshf,*) itmp,itmp,itmp
    end do    
    read(meshf,*) ctmp
  endif

  ! Triangles  (in temporary variables)
  if (ctmp.eq.'Triangles') then
    read(meshf,*) NFace
    write(*,*) ctmp, " ", NFace
    allocate(Face_IEN(NFACE,mesh%NSHLb))
    allocate(FaceID(NFace))  
    do i = 1, NFace
      read(meshf,*) (Face_IEN(i,j), j = 1, mesh%NSHLb),FaceID(i)  
    end do 
    read(meshf,*) ctmp
  endif
    
  ! Tetrahedra  
  read(meshf,*) mesh%NElem
  write(*,*) ctmp, " ", mesh%NElem  
  allocate(mesh%IEN(mesh%NElem,mesh%NSHL)) 
  do i = 1,mesh%NElem
    read(meshf,*) mesh%IEN(i,1),mesh%IEN(i,2),mesh%IEN(i,4),mesh%IEN(i,3), itmp
  end do 

  !==============================
  ! Convert face data in boundary data 
  !==============================
 
  ! Generate faceid 2 boundary map  
  NFaceID = maxval(FaceID)
  allocate(FID2BND(NFaceID))
  FID2BND = 0
  do i = 1, NFace
    FID2BND(FaceID(i)) = 1 
  end do 
  mesh%NBOUND = 0 
  do i = 1, NFaceID
    if (FID2BND(i).ge.1) then
      mesh%NBOUND = mesh%NBOUND  + 1
      FID2BND(i) = mesh%NBOUND
    endif        
  end do   
  
  allocate(mesh%bound(mesh%NBOUND))
  
  itmp = 0 
  do i = 1, NFaceID
    if (FID2BND(i).ge.1) then
      itmp = itmp  + 1
      
      mesh%bound(itmp)%FACE_ID = i
      
      mesh%bound(itmp)%NFACE = 0
      mesh%bound(itmp)%NNODE = 0      
    endif
  enddo 

 ! Get boundary faces    
  do i = 1, NFace
    j = FID2BND(FaceID(i))
    mesh%bound(j)%NFACE = mesh%bound(j)%NFACE + 1
  end do  
   
  write(*,*) "     Boundary     ID     NFACE" 
  do i = 1, mesh%NBOUND 
    write(*,*) i ,mesh%bound(i)%FACE_ID,mesh%bound(i)%NFACE 
    allocate(mesh%bound(i)%FACE_IEN(mesh%bound(i)%NFACE,mesh%NSHLB))
    allocate(mesh%bound(i)%FACE_OR (mesh%bound(i)%NFACE))
    allocate(mesh%bound(i)%F2E     (mesh%bound(i)%NFACE))
    mesh%bound(i)%NFACE = 0
  end do
  
  do i = 1, NFace
    j = FID2BND(FaceID(i))
    mesh%bound(j)%NFACE = mesh%bound(j)%NFACE + 1
    mesh%bound(j)%FACE_IEN(mesh%bound(j)%NFACE,:) = FACE_IEN(i,:)
  end do
  

  ! Find elements connected with faces
  write(*,*) "Get face elements " 
  allocate(cnt(mesh%NNODE))
  cnt = 0
  do iel = 1, mesh%NELEM
    do i = 1, mesh%NSHL
      n = mesh%IEN(iel,i)
      cnt(n) = cnt(n)+ 1 
    enddo
  enddo

  write(*,*) MAXVAL(cnt)
  allocate(INE(mesh%NNODE,MAXVAL(cnt)))
  
  cnt = 0
  do iel = 1, mesh%NELEM
    do i = 1, mesh%NSHL
      n = mesh%IEN(iel,i)
      cnt(n) = cnt(n)+ 1 
      INE(n,cnt(n)) = iel
    enddo
  enddo
  

  do i = 1, mesh%NBOUND 

    mesh%bound(i)%F2E = 0
    NFACE = 0
    do ifac = 1, mesh%bound(i)%NFACE 
      n = mesh%bound(i)%FACE_IEN(ifac,1)
      do ii = 1, cnt(n)
        iel = INE(n,ii)  
          
        found = 0
        do j = 1,3 
          if ((mesh%bound(i)%FACE_IEN(ifac,j).eq.mesh%IEN(iel,1)).or. &
              (mesh%bound(i)%FACE_IEN(ifac,j).eq.mesh%IEN(iel,2)).or. & 
              (mesh%bound(i)%FACE_IEN(ifac,j).eq.mesh%IEN(iel,3)).or. & 
              (mesh%bound(i)%FACE_IEN(ifac,j).eq.mesh%IEN(iel,4)) ) then              
          found = found + 1
          endif
        enddo
      
        if (found.eq.3) then
          mesh%bound(i)%F2E(ifac) = iel
	  NFACE = NFACE + 1
          exit
        endif                        
      enddo        
    enddo 
    
    write(*,*)  mesh%bound(i)%FACE_ID, mesh%bound(i)%NFACE ,NFACE, mesh%bound(i)%NFACE - NFACE
    
  enddo
  
  deallocate(INE,cnt)
 
  ! Reorder face connectivity 
  write(*,*) "Reorder faces" 
  do j = 1, mesh%NBOUND 
    write(*,*) "Reorder faces of bound",j
    do ifac = 1, mesh%bound(j)%NFACE

      iel = mesh%bound(j)%F2E(ifac)
      flist(:) = mesh%bound(j)%FACE_IEN(ifac,:)
      
      elist1(1) = mesh%IEN(iel,2)
      elist1(2) = mesh%IEN(iel,3)
      elist1(3) = mesh%IEN(iel,1)         

      elist2(1) = mesh%IEN(iel,4)
      elist2(2) = mesh%IEN(iel,3)  !2
      elist2(3) = mesh%IEN(iel,2)  !4       

      elist3(1) = mesh%IEN(iel,1)
      elist3(2) = mesh%IEN(iel,3)
      elist3(3) = mesh%IEN(iel,4)

      elist4(1) = mesh%IEN(iel,2)
      elist4(2) = mesh%IEN(iel,1)
      elist4(3) = mesh%IEN(iel,4)         

      ! Orientation 1  
      found = 0  
      do i = 1, 3
        if (  (flist(1).eq.elist1(i)).or. &
              (flist(2).eq.elist1(i)).or. &
              (flist(3).eq.elist1(i)) ) then
             found = found + 1
        endif
      enddo

      if (found.eq.3) then
         mesh%bound(j)%FACE_IEN(ifac,:) = elist1(:)
         mesh%bound(j)%FACE_OR(ifac) = 1
      endif 
     

      ! Orientation 2   
      found = 0  
      do i = 1, 3
        if (  (flist(1).eq.elist2(i)).or. &
              (flist(2).eq.elist2(i)).or. &
              (flist(3).eq.elist2(i)) ) then
             found = found + 1
        endif
      enddo

      if (found.eq.3) then
         mesh%bound(j)%FACE_IEN(ifac,:) = elist2(:)
         mesh%bound(j)%FACE_OR(ifac) = 2
      endif 
     

      ! Orientation 3   
      found = 0  
      do i = 1, 3
        if (  (flist(1).eq.elist3(i)).or. &
              (flist(2).eq.elist3(i)).or. &
              (flist(3).eq.elist3(i)) ) then
             found = found + 1
        endif
      enddo

      if (found.eq.3) then
         mesh%bound(j)%FACE_IEN(ifac,:) = elist3(:)
         mesh%bound(j)%FACE_OR(ifac) = 3
      endif 
     
      ! Orientation 4   
      found = 0  
      do i = 1, 3
        if (  (flist(1).eq.elist4(i)).or. &
              (flist(2).eq.elist4(i)).or. &
              (flist(3).eq.elist4(i)) ) then
             found = found + 1
        endif
      enddo

      if (found.eq.3) then
         mesh%bound(j)%FACE_IEN(ifac,:) = elist4(:)
         mesh%bound(j)%FACE_OR(ifac) = 4
      endif 
     
    enddo
  enddo

  ! Get face nodes
  write(*,*) "Get face nodes" 
  allocate(nflag(mesh%NNODE ))
  do i = 1, mesh%NBOUND 
    nflag = .false.
    do ifac = 1, mesh%bound(i)%NFACE
      do j = 1, mesh%NSHLB
         nflag(mesh%bound(i)%FACE_IEN(ifac,j)) = .true.
      enddo
    enddo
    
    itmp = 0
    do j = 1, mesh%NNODE
      if (nflag(j))  itmp = itmp + 1
    enddo
    mesh%bound(i)%NNODE = itmp
    allocate(mesh%bound(i)%BNODES(itmp))
    
    itmp = 0
    do j = 1, mesh%NNODE
      if (nflag(j))  then
        itmp = itmp + 1
        mesh%bound(i)%BNODES(itmp) = j
      endif  
    enddo        
  enddo    

  !==============================
  ! deallocate
  !==============================   
  deallocate(Face_IEN, FaceID)
  deallocate(FID2BND)
  deallocate(nflag)

end subroutine input_tet
