!======================================================================
! subroutine computes the volumn of a tetrahedron
!======================================================================
subroutine voltet(x1, x2, x3, x4, vol)
  implicit none
  real(8), intent(in)  :: x1(3), x2(3), x3(3), x4(3)
  real(8), intent(out) :: vol

  vol = x2(1)*x3(2)*x4(3) - x2(1)*x3(3)*x4(2) - x3(1)*x2(2)*x4(3) &
      + x3(1)*x2(3)*x4(2) + x4(1)*x2(2)*x3(3) - x4(1)*x2(3)*x3(2) &
      - x1(1)*x3(2)*x4(3) + x1(1)*x3(3)*x4(2) + x3(1)*x1(2)*x4(3) &
      - x3(1)*x1(3)*x4(2) - x4(1)*x1(2)*x3(3) + x4(1)*x1(3)*x3(2) &
      + x1(1)*x2(2)*x4(3) - x1(1)*x2(3)*x4(2) - x2(1)*x1(2)*x4(3) &
      + x2(1)*x1(3)*x4(2) + x4(1)*x1(2)*x2(3) - x4(1)*x1(3)*x2(2) &
      - x1(1)*x2(2)*x3(3) + x1(1)*x2(3)*x3(2) + x2(1)*x1(2)*x3(3) &
      - x2(1)*x1(3)*x3(2) - x3(1)*x1(2)*x2(3) + x3(1)*x1(3)*x2(2)

  vol = -vol/6.0d0
end subroutine voltet


!======================================================================
! This subroutine reads in the mesh in Gmsh standard format (.msh)
! and converts it into our format.
! It currently supports triangles, tetrahedron, and prism.
!======================================================================
subroutine input_fem(mesh, fname)

  use class      
  implicit none
  
  type(mesh_class) :: mesh
  character*(*) :: fname
         
  integer, allocatable :: Face_IEN(:,:), FaceID(:), FID2BND(:), G2L(:)
  integer, allocatable :: INE(:,:), IEN(:,:), E_Type(:), E_Flag(:)
  integer, allocatable :: cnt(:)  
  real(8), allocatable :: dg(:,:)
  logical, allocatable :: nflag(:)   

  integer :: i, j, k, ier, n, ii, meshf, itmp, ifac, iel, found
  integer :: NElem, NTet, NPrism, maxNSHL
  integer :: flist(3), elist1(3), elist2(3), elist3(3), elist4(3)
  integer :: NFace, NFaceID, Nedge, NSHL(20)
  real(8) :: vol, vol1, vol2, vol3          
  character(20) :: ctmp

  maxNSHL     = 6
  mesh%NPATCH = 0
! mesh%NSHL   = 4
  mesh%NSHLb  = 3
  
  mesh%NSD = 3

  ! Gmsh format (2-triangle,4-tet,6-prism)
  NSHL    = -1  
  NSHL(2) = 3
  NSHL(4) = 4
  NSHL(6) = 6

  write(*,*) "Dimension", mesh%NSD

  !------------------------------
  ! Read gmsh .mesh file
  !------------------------------
  meshf = 90
  open(meshf, file = fname, status = 'old')
  write(*,*) "Reading ", fname   

  ! Nodes
  do while (ctmp .ne.'$Nodes')  
    read(meshf,*) ctmp
  end do
  read(meshf,*) mesh%NNode
  write(*,*) ctmp, mesh%NNode
  allocate(mesh%xg(mesh%NNode,mesh%NSD), mesh%NID(mesh%NNode))
  mesh%NID = -1
  do i = 1, mesh%NNode
    read(meshf,*) itmp, (mesh%xg(i,j), j = 1, mesh%NSD)
  end do      

  ! Elements (Triangles+Tet+Prism in temporary array)  
  do while (ctmp .ne.'$Elements')  
    read(meshf,*) ctmp
  end do
  read(meshf,*) NElem   
  write(*,*) ctmp, NElem
  allocate(IEN(NElem,maxNSHL), E_Type(NElem), E_Flag(NElem))
  IEN = -1; E_Type = -1; E_Flag = -1
  do i = 1, NElem
    read(meshf,*) itmp, E_Type(i), itmp, itmp, E_Flag(i), &
                  (IEN(i,j), j = 1, NSHL(E_Type(i)))
  end do 
  close(meshf) 
  
  !---------------------------------------
  ! Split elements in interior/boundary
  !---------------------------------------
  NFace = 0; NTet = 0; NPrism = 0
  mesh%NElem = 0
  do i = 1, NElem
    if (E_Type(i) == 2) then
      NFace = NFace + 1
    else if (E_Type(i) == 4) then
      NTet = NTet + 1
    else if (E_Type(i) == 6) then
      NPrism = NPrism + 1
    else
      write(*,*) "Undefined element type", E_Type(i)
      stop
    end if
  end do

  mesh%NElem = NTet + NPrism
  write(*,*) "Elements(vol;tet;prism;tri):", mesh%NElem, NTet, NPrism, NFace

  allocate(mesh%IEN(mesh%NElem,maxNSHL), mesh%EID(mesh%NElem), &
           Face_IEN(NFace,mesh%NSHLb), FaceID(NFace), mesh%NSHL(NElem))

  NFace      = 0
  mesh%NElem = 0
  mesh%EID   = -1
  do i = 1, NElem   
    if (E_Type(i) == 2) then
      NFace = NFace + 1                         ! Number of Face element
      Face_IEN(NFACE,:) = IEN(i,1:mesh%NSHLb)   ! IEN is embeded into face-related IEN
      FaceID(NFACE)     = E_Flag(i)             ! Face ID following E_flag

    else if (E_Type(i)==4 .or. E_Type(i)==6) then  ! Element within volummesh
      mesh%NElem = mesh%NElem + 1    
      mesh%NSHL(mesh%NElem) = E_Type(i)   

      if (mesh%NSHL(mesh%NElem) == 4) then  
        mesh%IEN(mesh%NElem,1) = IEN(i,1)      ! 1,4,3,2 ==> works for pure tet
        mesh%IEN(mesh%NElem,2) = IEN(i,2)
        mesh%IEN(mesh%NElem,3) = IEN(i,4)      ! Switch ordering of tet to get positive DetJ
        mesh%IEN(mesh%NElem,4) = IEN(i,3)
        mesh%EID(mesh%NElem)   = E_Flag(i)

        ! setup node ID using element ID: if the node is shared by several
        ! elements, the node ID is always the larger element ID
        do k = 1, mesh%NSHL(mesh%NElem)
          if (mesh%NID(IEN(i,k)) < E_Flag(i)) then
            mesh%NID(IEN(i,k)) = E_Flag(i)
          end if
        end do

        call voltet(mesh%xg(mesh%IEN(mesh%NElem,1),1:3), &
                    mesh%xg(mesh%IEN(mesh%NElem,2),1:3), &
                    mesh%xg(mesh%IEN(mesh%NElem,3),1:3), &
                    mesh%xg(mesh%IEN(mesh%NElem,4),1:3), vol)
        if (vol <= 0d0) then
          write(*,*) "%Negative vol = ", vol, "in tet", mesh%NElem
!!!          stop
        end if

      else if (mesh%NSHL(mesh%NElem) == 6) then  
        mesh%IEN(mesh%NElem,:) = IEN(i,:)  
        mesh%EID(mesh%NElem)   = E_Flag(i)

        do k = 1, mesh%NSHL(mesh%NElem)
          if (mesh%NID(IEN(i,k)) < E_Flag(i)) then
            mesh%NID(IEN(i,k)) = E_Flag(i)
          end if
        end do

        call voltet(mesh%xg(mesh%IEN(mesh%NElem,1),1:3), &
                    mesh%xg(mesh%IEN(mesh%NElem,2),1:3), &
                    mesh%xg(mesh%IEN(mesh%NElem,6),1:3), &
                    mesh%xg(mesh%IEN(mesh%NElem,3),1:3), vol1)
        call voltet(mesh%xg(mesh%IEN(mesh%NElem,2),1:3), &
                    mesh%xg(mesh%IEN(mesh%NElem,4),1:3), &
                    mesh%xg(mesh%IEN(mesh%NElem,1),1:3), &
                    mesh%xg(mesh%IEN(mesh%NElem,6),1:3), vol2)
        call voltet(mesh%xg(mesh%IEN(mesh%NElem,4),1:3), &
                    mesh%xg(mesh%IEN(mesh%NElem,5),1:3), &
                    mesh%xg(mesh%IEN(mesh%NElem,6),1:3), &
                    mesh%xg(mesh%IEN(mesh%NElem,2),1:3), vol3)
        vol = vol1 + vol2 + vol3
        if (vol <= 0d0) then
          write(*,*) "%Negative vol = ", vol, "in prism", mesh%NElem
!!!          stop
        end if
      else 
        write(*,*) "Undefined shape function:", mesh%NSHL(mesh%NElem)
        stop
      end if

    else 
      write(*,*) "Undefined element type", E_Type(i)
      stop
    end if
  end do 

  mesh%NSHLmax = MAXVAL(mesh%NSHL)

  write(*,*) "Elements(vol;srf):", mesh%NElem, NFace

  if (mesh%NSHLmax > maxNSHL) then
    write(*,*) "ERROR: mesh%NSHLmax and maxNSHL don't match"
    write(*,*) "mesh%NSHLmax =", mesh%NSHLmax
    write(*,*) "maxNSHL      =", maxNSHL
    stop
  end if

  !---------------------------------------
  ! Convert face data in boundary data 
  !---------------------------------------
 
  ! Generate faceid 2 boundary map  
  NFaceID = maxval(E_Flag)

  allocate(FID2BND(NFaceID))
  FID2BND = 0
  do i = 1, NFace
    FID2BND(FaceID(i)) = 1 
  end do 
  mesh%NBOUND = 0 
  do i = 1, NFaceID
    if (FID2BND(i).ge.1) then
      mesh%NBOUND = mesh%NBOUND  + 1
      FID2BND(i)  = mesh%NBOUND
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
    end if
  end do 

  ! Get boundary faces    
  do i = 1, NFace
    j = FID2BND(FaceID(i))
    mesh%bound(j)%NFACE = mesh%bound(j)%NFACE + 1
  end do  

  write(*,*) "     Boundary     ID                    NFACE" 
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
    do i = 1, mesh%NSHL(iel)
      n = mesh%IEN(iel,i)
      if (n==0) write(*,*) iel, i , n
      cnt(n) = cnt(n)+ 1 
    enddo
  enddo

  allocate(INE(mesh%NNODE,maxval(cnt)))
  
  cnt = 0
  do iel = 1, mesh%NELEM
    do i = 1, mesh%NSHL(iel)
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
        if (mesh%NSHL(iel) == 4) then
          do j = 1, 3 
            if ((mesh%bound(i)%FACE_IEN(ifac,j) == mesh%IEN(iel,1)).or. &
                (mesh%bound(i)%FACE_IEN(ifac,j) == mesh%IEN(iel,2)).or. & 
                (mesh%bound(i)%FACE_IEN(ifac,j) == mesh%IEN(iel,3)).or. & 
                (mesh%bound(i)%FACE_IEN(ifac,j) == mesh%IEN(iel,4)) ) then              
                found = found + 1
            end if
          end do
        !! Do similar for prism
        else if (mesh%NSHL(iel) == 6) then
          do j = 1, 3
            if ((mesh%bound(i)%FACE_IEN(ifac,j) == mesh%IEN(iel,1)).or. &
                (mesh%bound(i)%FACE_IEN(ifac,j) == mesh%IEN(iel,2)).or. & 
                (mesh%bound(i)%FACE_IEN(ifac,j) == mesh%IEN(iel,3)).or. & 
                (mesh%bound(i)%FACE_IEN(ifac,j) == mesh%IEN(iel,4)).or. &              
                (mesh%bound(i)%FACE_IEN(ifac,j) == mesh%IEN(iel,5)).or. & 
                (mesh%bound(i)%FACE_IEN(ifac,j) == mesh%IEN(iel,6)) ) then
              found = found + 1
            end if      
          end do
        end if

        if (found == 3) then
          itmp = mesh%bound(i)%F2E(ifac)
          ! The facet is NOT assigned to an element before
          if(itmp == 0) then
            mesh%bound(i)%F2E(ifac) = iel
	        NFACE = NFACE + 1
          ! The facet has been assigned to an element
          ! and the new element has a greater ID than the old one
          else if(E_FLAG(itmp) < E_FLAG(iel)) then
            mesh%bound(i)%F2E(ifac) = iel
            ! assert FACE_ID == 5
            if(mesh%bound(i)%FACE_ID .ne. 5) write(*,*) "WTF: This internal facet is not 5"
          endif
        end if  
        
      end do        
    end do 
   
    write(*,*) "Face ID", "# of Boundary Face","# of Face"  
    write(*,*) mesh%bound(i)%FACE_ID, mesh%bound(i)%NFACE, NFACE, mesh%bound(i)%NFACE - NFACE

  end do

 
  deallocate(INE,cnt)
 

  ! Reorder face connectivity 
  write(*,*) "Reorder faces" 
  do j = 1, mesh%NBOUND 
    write(*,*) "Reorder faces of bound",j
    do ifac = 1, mesh%bound(j)%NFACE

      iel = mesh%bound(j)%F2E(ifac)      
      flist(:) = mesh%bound(j)%FACE_IEN(ifac,:)
     
      if (mesh%NSHL(iel)  == 4) then
      
        elist1(1) = mesh%IEN(iel,2)
        elist1(2) = mesh%IEN(iel,3)
        elist1(3) = mesh%IEN(iel,1)         

        elist2(1) = mesh%IEN(iel,4)
        elist2(2) = mesh%IEN(iel,3)
        elist2(3) = mesh%IEN(iel,2)       

        elist3(1) = mesh%IEN(iel,1)
        elist3(2) = mesh%IEN(iel,3)
        elist3(3) = mesh%IEN(iel,4)

        elist4(1) = mesh%IEN(iel,2)
        elist4(2) = mesh%IEN(iel,1)
        elist4(3) = mesh%IEN(iel,4)         

        ! Orientation 1  
        found = 0  
        do i = 1, 3
          if ( (flist(1) == elist1(i)).or. &
               (flist(2) == elist1(i)).or. &
               (flist(3) == elist1(i)) ) then
            found = found + 1
          end if
        end do

        if (found == 3) then
          mesh%bound(j)%FACE_IEN(ifac,:) = elist1(:)
          mesh%bound(j)%FACE_OR(ifac) = 1
        end if     

        ! Orientation 2   
        found = 0  
        do i = 1, 3
          if ( (flist(1) == elist2(i)).or. &
               (flist(2) == elist2(i)).or. &
               (flist(3) == elist2(i)) ) then
            found = found + 1
          end if
        end do

        if (found == 3) then
          mesh%bound(j)%FACE_IEN(ifac,:) = elist2(:)
          mesh%bound(j)%FACE_OR(ifac) = 2
        end if      

        ! Orientation 3   
        found = 0  
        do i = 1, 3
          if ( (flist(1) == elist3(i)).or. &
               (flist(2) == elist3(i)).or. &
               (flist(3) == elist3(i)) ) then
            found = found + 1
          end if
        end do

        if (found == 3) then
          mesh%bound(j)%FACE_IEN(ifac,:) = elist3(:)
          mesh%bound(j)%FACE_OR(ifac) = 3
        end if 
     
        ! Orientation 4   
        found = 0  
        do i = 1, 3
          if ( (flist(1) == elist4(i)).or. &
               (flist(2) == elist4(i)).or. &
               (flist(3) == elist4(i)) ) then
            found = found + 1
          end if
        end do

        if (found == 3) then
           mesh%bound(j)%FACE_IEN(ifac,:) = elist4(:)
           mesh%bound(j)%FACE_OR(ifac) = 4
        end if 

      ! Prism
      ! Need to know orientation of triangle nodes
      else if (mesh%NSHL(iel) == 6) then
      
        elist1(1) = mesh%IEN(iel,3)
        elist1(2) = mesh%IEN(iel,2)
        elist1(3) = mesh%IEN(iel,1)         

        elist2(1) = mesh%IEN(iel,4)
        elist2(2) = mesh%IEN(iel,5)  
        elist2(3) = mesh%IEN(iel,6)                

        ! Orientation 1  
        found = 0  
        do i = 1, 3
          if ( (flist(1) == elist1(i)).or. &
               (flist(2) == elist1(i)).or. &
               (flist(3) == elist1(i)) ) then
            found = found + 1
          end if
        end do

        if (found == 3) then
          mesh%bound(j)%FACE_IEN(ifac,:) = elist1(:)
          mesh%bound(j)%FACE_OR(ifac) = 1
        end if 
     
        ! Orientation 2   
        found = 0  
        do i = 1, 3
          if ( (flist(1) == elist2(i)).or. &
               (flist(2) == elist2(i)).or. &
               (flist(3) == elist2(i)) ) then
            found = found + 1
          end if
        end do

        if (found == 3) then
          mesh%bound(j)%FACE_IEN(ifac,:) = elist2(:)
          mesh%bound(j)%FACE_OR(ifac) = 2
        end if
     
      end if
    end do
  end do

  ! Get face nodes
  write(*,*) "Get face nodes" 
  allocate(nflag(mesh%NNODE ))
  do i = 1, mesh%NBOUND 
    nflag = .false.
    do ifac = 1, mesh%bound(i)%NFACE
      do j = 1, mesh%NSHLB
         nflag(mesh%bound(i)%FACE_IEN(ifac,j)) = .true.
      end do
    end do
    
    itmp = 0
    do j = 1, mesh%NNODE
      if (nflag(j))  itmp = itmp + 1
    end do
    mesh%bound(i)%NNODE = itmp
    allocate(mesh%bound(i)%BNODES(itmp))
    
    itmp = 0
    do j = 1, mesh%NNODE
      if (nflag(j))  then
        itmp = itmp + 1
        mesh%bound(i)%BNODES(itmp) = j        
      end if  
    end do        

  end do    

  !---------------------------------------
  ! deallocate
  !---------------------------------------
  deallocate(Face_IEN, FaceID, FID2BND, nflag)

end subroutine input_fem
