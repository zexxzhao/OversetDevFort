!------------------------------------------------------------------------
! Get Patch connectivity             
!------------------------------------------------------------------------
subroutine getPatchConn(conname,NMESH, mesh, connectivity)
      use class

      implicit none
      
      character*(*) conname
      integer NMESH
      type (mesh_class) mesh(NMESH)
      type (patch_con_class) connectivity             
      logical chk
  
      inquire(file=trim(conname), exist=chk)
  
      if (chk) then
        call readPatchConn(conname,connectivity)
      else
        call genPatchConn   (NMESH, mesh, connectivity)
        call writePatchConn (conname,connectivity)
      endif

end subroutine getPatchConn
!------------------------------------------------------------------------
subroutine genPatchConn(NMESH, mesh, connectivity)
      use class

      implicit none
      
      integer NMESH
      type (mesh_class) mesh(NMESH)
      type (patch_con_class) connectivity    
      integer :: ConNP(NMESH*(NMESH-1)/2,2)     
      real(8) :: boxmax(mesh(1)%NSD,NMESH)      
      real(8) :: boxmin(mesh(1)%NSD,NMESH)
      real(8) :: eps, dist(3)
      integer :: mn1,mn2
      
      !----------------------------------------------
      ! Loop over patches to get bounding boxes 
      !----------------------------------------------                
      do mn1 = 1, NMESH 
        boxmin(:,mn1) = MINVAL(mesh(mn1)%xg(:,:), DIM=1)
        boxmax(:,mn1) = MAXVAL(mesh(mn1)%xg(:,:), DIM=1) 
      enddo       
      
      !----------------------------------------------
      ! Compute tolerance based on total domain size
      !----------------------------------------------               
      dist = MAXVAL(boxmax(:,:), DIM=2)&
           - MINVAL(boxmin(:,:), DIM=2) 
      connectivity%dsize = sqrt(sum(dist*dist))
      eps = connectivity%dsize*connectivity%eps
        
      ConNP  = 0 
      connectivity%npc = 0      
      do mn1 = 1, NMESH 
        do mn2 = mn1+1, NMESH

          if ((boxmin(1,mn1).lt.(boxmax(1,mn2)+eps)).and. &
              (boxmin(2,mn1).lt.(boxmax(2,mn2)+eps)).and. &
              (boxmin(3,mn1).lt.(boxmax(3,mn2)+eps)).and. &
              (boxmax(1,mn1).gt.(boxmin(1,mn2)-eps)).and. & 
              (boxmax(2,mn1).gt.(boxmin(2,mn2)-eps)).and. &
              (boxmax(3,mn1).gt.(boxmin(3,mn2)-eps)) )  then  
             connectivity%npc = connectivity%npc + 1  
             ConNP(connectivity%npc,1) = mn1
             ConNP(connectivity%npc,2) = mn2
             write(*,*) "Connected :",mn1, mn2
          endif
      
        enddo    
      enddo 
       
      !----------------------------------------------
      ! Copy connectivity
      !----------------------------------------------                  
      allocate(connectivity%patch(connectivity%npc ,2))
      do mn1 = 1, connectivity%npc                
        connectivity%patch(mn1,:) = ConNP(mn1,:)
      enddo
            
                       

end subroutine genPatchConn   
!------------------------------------------------------------------------
subroutine writePatchConn(fname,connectivity)
      use class

      implicit none
      
      character*(*) fname
      type (patch_con_class) connectivity  
            
      integer :: cfid,c
                       
      write(*,*) "Write  connectivity: ", fname    
      cfid = 15
      open (cfid, file=trim(fname), status='new')
            
      write(cfid,*) connectivity%npc,connectivity%dsize  
      do c = 1, connectivity%npc                 
         write(cfid,*) connectivity%patch(c,:)
      enddo  
      
      close (cfid)
      
end subroutine writePatchConn  
!------------------------------------------------------------------------
subroutine readPatchConn(fname,connectivity)
      use class

      implicit none
      
      character*(*) fname
      type (patch_con_class)  connectivity  
            
      integer :: cfid,c
                       
      write(*,*) "Read  connectivity: ", fname    
      cfid = 15
      open (cfid, file=trim(fname), status='old')
                                  
      read(cfid,*) connectivity%npc,connectivity%dsize
      allocate(connectivity%patch(connectivity%npc ,2))  
      do c = 1, connectivity%npc                 
         read(cfid,*) connectivity%patch(c,:)
      enddo  
      
      close (cfid)
      
end subroutine readPatchConn  


!------------------------------------------------------------------------
! Get Patch boundary connectivity             
!------------------------------------------------------------------------
subroutine getPatchBndConn(conname, NMESH, mesh, connectivity)  
      use class

      implicit none
      
      character*(*) conname
      integer NMESH
      type (mesh_class) mesh(NMESH)
      type (patch_con_class)  connectivity 
      logical chk
  
      inquire(file=trim(conname), exist=chk)
  
      if (chk) then
        call readPatchBndConn (conname,connectivity )
      else
        call genPatchBndConn  (NMESH, mesh, connectivity)
        call writePatchBndConn(conname,connectivity)
      endif

end subroutine getPatchBndConn
!------------------------------------------------------------------------
subroutine genPatchBndConn(NMESH, mesh, connectivity)
      use class

      implicit none
      
      integer NMESH
      type (mesh_class) mesh(NMESH)
      type (patch_con_class)  connectivity 
     
      integer c,p1,p2,b1,b2,ln1,ln2,gn1,gn2,found,mode
      real(8) dist(mesh(1)%NSD),eps2
      
      integer bndconn(connectivity%npc*6*6,5)
      
      eps2 = (connectivity%dsize*connectivity%eps)**2
        
      connectivity%nbc = 0
        
      !----------------------------------------------
      ! Loop through patch connectivity
      !----------------------------------------------
      connectivity%nbc = 0
      do c = 1, connectivity%npc
         p1 = connectivity%patch(c,1)       
         p2 = connectivity%patch(c,2)
         write(*,*) "----------------------------------"
         write(*,*) ' Patch:', p1,p2

         ! loop through master boundary  nodes       
         do b1 = 1, mesh(p1)%NBOUND
         do b2 = 1, mesh(p2)%NBOUND

            found  = 0 
            
            do ln1 = 1, mesh(p1)%bound(b1)%NNODE
            do ln2 = 1, mesh(p2)%bound(b2)%NNODE
               gn1 = mesh(p1)%bound(b1)%BNODES(ln1) 
               gn2 = mesh(p2)%bound(b2)%BNODES(ln2)  
 
               dist = mesh(p1)%xg(gn1,:)-mesh(p2)%xg(gn2,:)           
          
               if (sum(dist*dist) < eps2 ) then                  
                 found = found + 1
               endif
          
            enddo
            enddo

            if (found.ge.mesh(1)%NSHLB) then                  
              !----------------------------------------------
              ! Find type of match
              !----------------------------------------------
              mode = -1  
                          
              if ((found.gt.mesh(p1)%bound(b1)%NNODE).or. &
                  (found.gt.mesh(p2)%bound(b2)%NNODE) ) then
                write(*,*)"     Found to many connecting points for bounds:", b1,b2 
                write(*,*) found,mesh(p1)%bound(b1)%NNODE,mesh(p2)%bound(b2)%NNODE  
                write(*,*)"  ***    -----------------------   ***" 
                write(*,*)"  ***    Adjust connectivity eps   ***" 
                write(*,*)"  ***    -----------------------   ***"         
                mode  = 1                                        
              else if( ((found.eq.mesh(p1)%bound(b1)%NNODE).and. &
                        (found.ge.mesh(p2)%bound(b2)%NNODE) ) .or. &
                       ((found.eq.mesh(p1)%bound(b1)%NNODE).and. &
                        (found.ge.mesh(p2)%bound(b2)%NNODE) ) ) then
                write(*,*)"   Full internal boundary:", b1,b2 
                mode = 1                                                             
              else if ((real(found)/real(mesh(p1)%bound(b1)%NNODE).gt.connectivity%treshold).or. &
                       (real(found)/real(mesh(p2)%bound(b2)%NNODE).gt.connectivity%treshold)) then 

                write(*,*)"   Partial internal boundary:", b1,b2                
                write(*,*) found,mesh(p1)%bound(b1)%NNODE
                write(*,*) found,mesh(p2)%bound(b2)%NNODE 
                write(*,*) real(found)/real(mesh(p1)%bound(b1)%NNODE)*100d0,"%" 
                write(*,*) real(found)/real(mesh(p2)%bound(b2)%NNODE)*100d0,"%"     
                mode = 0
              endif  
              
              if (mode.ge.0) then    
                connectivity%nbc = connectivity%nbc + 1
                bndconn(connectivity%nbc, :) = (/p1,p2,b1,b2,mode/)
              endif
            endif
            
         enddo
         enddo
       enddo 
       
      !----------------------------------------------
      ! Copy connectivity
      !----------------------------------------------                  
      allocate(connectivity%bnd(connectivity%nbc ,5))
           
      do c =1, connectivity%nbc    
         connectivity%bnd(c, :) = bndconn(c, :)
      enddo
      write(*,*) "----------------------------------"
                    
end subroutine genPatchBndConn  
 
!------------------------------------------------------------------------
subroutine writePatchBndConn(fname,connectivity)
      use class

      implicit none
      
      character*(*) fname
      type (patch_con_class) connectivity  
            
      integer :: cfid,c
                       
      write(*,*) "Write  connectivity: ", fname    
      cfid = 15
      open (cfid, file=trim(fname), status='new')
            
      write(cfid,*) connectivity%nbc   
      do c = 1, connectivity%nbc                 
         write(cfid,*) connectivity%bnd(c,:)
      enddo  
      
      close (cfid)
      
end subroutine writePatchBndConn  
!------------------------------------------------------------------------
subroutine readPatchBndConn(fname,connectivity)
      use class

      implicit none
      
      character*(*) fname
      type (patch_con_class)  connectivity  
            
      integer :: cfid,c
                       
      write(*,*) "Read  connectivity: ", fname    
      cfid = 15
      open (cfid, file=trim(fname), status='old')
                                  
      read(cfid,*) connectivity%nbc  
      allocate(connectivity%bnd(connectivity%nbc ,5))  
      do c = 1, connectivity%nbc                 
         read(cfid,*) connectivity%bnd(c,:)
      enddo  
      
      close (cfid)
      
end subroutine readPatchBndConn  
   
!------------------------------------------------------------------------
!                 
!------------------------------------------------------------------------
subroutine findNodeMap(NMESH,mesh,connectivity,nodemap)      
      use class

      implicit none

      integer NMESH
      type (mesh_class) mesh(NMESH)      
      type (patch_con_class)  connectivity 
      type (map_class)  nodemap(NMESH)

      integer c,p1,p2,b1,b2,ln1,ln2,gn1,gn2,f1,f2,j,mode,NNODE, NFACe
      real(8) dist(mesh(1)%NSD),eps2
      logical addFace
      logical, allocatable :: found1(:), found2(:)
      integer, allocatable :: Face_IEN(:,:), F2E(:), FACE_OR(:), BNODES(:)
     
      eps2 = (connectivity%dsize*connectivity%eps)**2
                   
      !---------------------------------------------- 
      ! Get local to local connectivity
      !----------------------------------------------            
      do p1 = 1, NMESH
        allocate(nodemap(p1)%node(mesh(p1)%NNODE,2))
        nodemap(p1)%node(:,1) = p1
        do gn1 = 1,  mesh(p1)%NNODE
          nodemap(p1)%node(gn1,2) = gn1
        enddo  
      enddo
      
      do c = 1, connectivity%nbc                 
        p1   = connectivity%bnd(c,1)                
        p2   = connectivity%bnd(c,2)                
        b1   = connectivity%bnd(c,3)                
        b2   = connectivity%bnd(c,4)                
        mode = connectivity%bnd(c,5)
        
        allocate(found1(mesh(p1)%NNODE))
        allocate(found2(mesh(p2)%NNODE))
        found1 = .false.
        found2 = .false.
                              
        write(*,*) "----------------------------------"
        write(*,*) "  Patches", p1, p2
        write(*,*) "  Bounds ", b1, b2
                                
        do ln1 = 1, mesh(p1)%bound(b1)%NNODE
        do ln2 = 1, mesh(p2)%bound(b2)%NNODE
           gn1 = mesh(p1)%bound(b1)%BNODES(ln1) 
           gn2 = mesh(p2)%bound(b2)%BNODES(ln2)  
 
           dist = mesh(p1)%xg(gn1,:)-mesh(p2)%xg(gn2,:)           
        
           if (sum(dist*dist) < eps2 ) then      
             found1(gn1) = .true.
             found2(gn2) = .true.                             
             if (nodemap(p2)%node(gn2,1).le.nodemap(p1)%node(gn1,1)) then           
               nodemap(p1)%node(gn1,:) = nodemap(p2)%node(gn2,:) 
             else
               nodemap(p2)%node(gn2,:) = nodemap(p1)%node(gn1,:)                  
             endif                          
           endif
        
        enddo
        enddo
        
        if (mode.eq.1) then             
          mesh(p1)%bound(b1)%FACE_ID = -1             
          mesh(p2)%bound(b2)%FACE_ID = -1   
        else
          call reduceBoundary(mesh(p1)%NSHLB, mesh(p1)%NNODE, mesh(p1)%bound(b1), found1) 
          call reduceBoundary(mesh(p2)%NSHLB, mesh(p2)%NNODE, mesh(p2)%bound(b2), found2)   
        endif         
        
        deallocate(found1,found2)                      
      enddo  
      write(*,*) "----------------------------------"

return
end 

!------------------------------------------------------------------------
!                 
!------------------------------------------------------------------------
subroutine reduceBoundary(NSHLB, NNODE, bound,found)   
      use class

      implicit none
      
      integer  NNODE,NSHLB
      type (bnd_class) bound
      logical  found(NNODE)

      integer NBNODE,NFACE,ln,gn,f,j
      logical addFace

      integer FACE_IEN(bound%NFACE,NSHLB), F2E(bound%NFACE), &
              FACE_OR(bound%NFACE), BNODES(bound%NNODE)

      ! Initilialize
      NFACE  = 0
      NBNODE = 0              
      
      ! Extract exterior boundary faces 
      do f = 1, bound%NFACE          
        addFace = .false.
        do j = 1,NSHLB  
          if (.not.found(bound%FACE_IEN(f,j))) addFace = .true.
        enddo
        if (addFace)then
          NFACE = NFACE + 1
          FACE_IEN(NFACE,:) = bound%FACE_IEN(f,:)
          F2E     (NFACE)   = bound%F2E(f)
          FACE_OR (NFACE)   = bound%FACE_OR(f)
        endif  
      enddo
      
      ! Tag nodes that are part of exterior boundary faces 
      found = .false.
      do f = 1, NFACE   
        do j = 1, NSHLB  
           found(FACE_IEN(f,j)) = .true.
        enddo
      enddo 
       
      ! Extract exterior boundary nodes    
      do ln = 1, bound%NNODE
        gn = bound%BNODES(ln)
        if (found(gn))then
          NBNODE = NBNODE + 1
          BNODES(NBNODE) = gn
        endif  
      enddo
      
      ! Copy newly generate boundary data 
      deallocate(bound%FACE_IEN)
      deallocate(bound%F2E) 
      deallocate(bound%FACE_OR) 
      deallocate(bound%BNODES)  
                
      bound%NFACE = NFACE
      bound%NNODE = NBNODE
      allocate(bound%FACE_IEN(NFACE,NSHLB))
      allocate(bound%F2E(NFACE)) 
      allocate(bound%FACE_OR(NFACE)) 
      allocate(bound%BNODES(NBNODE))   
      
      do f = 1,NFACE 
        bound%FACE_IEN(f,:) = FACE_IEN(f,:) 
        bound%F2E     (f)   = F2E     (f)   
        bound%FACE_OR (f)   = FACE_OR (f)   
      enddo 
      
      do ln = 1, NBNODE
        bound%BNODES(ln) = BNODES(ln)
      enddo 

return
end 

!------------------------------------------------------------------------
!                 
!------------------------------------------------------------------------
subroutine mergeMeshes(NMESH, mesh, nodemap, finalmesh)      
      use class

      implicit none

      integer NMESH
      type (mesh_class) mesh(NMESH)
      type (map_class)  nodemap(NMESH)
      type (mesh_class) finalmesh
      
                   
      integer :: NSD,NSHL,NSHLB
      integer :: TNNODE,NNODE,NELEM,NBOUND,NPATCH
      
      integer  :: mn,bn,pn,i,j,ln1,ln2,p1,p2      
      
      !---------------------------------------------- 
      !  Check combatibility of meshes    
      !----------------------------------------------  
      write(*,*) "Check patch dimension and shapefunction count"
      NSD   = mesh(1)%NSD
      NSHL  = mesh(1)%NSHLmax      
      NSHLB = mesh(1)%NSHLB    
      do mn = 1, NMESH
        if (mesh(mn)%NSD.ne.NSD) then
           write(*,*) "Patch :", mn
           write(*,*) "NSD :", NSD,mesh(mn)%NSD
        endif
        if (mesh(mn)%NSHLmax.ne.NSHL) then
           write(*,*) "Patch :", mn
           write(*,*) "NSHL :", NSHL,mesh(mn)%NSHL
        endif       
        if (mesh(mn)%NSHLB.ne.NSHLB) then
           write(*,*) "Patch :", mn
           write(*,*) "NSHLB :", NSHLB,mesh(mn)%NSHLB
        endif
      enddo 
      finalmesh%NSD  = NSD
      finalmesh%NSHL = NSHL        
      finalmesh%NSHLB = NSHLB   
      
      !---------------------------------------------- 
      !  Get local to global node map 
      !----------------------------------------------

      NNODE  = 0
      TNNODE = 0
      do p1 = 1, NMESH
        TNNODE = TNNODE + mesh(p1)%NNODE 
        do ln1 = 1, mesh(p1)%NNODE
  
          if (nodemap(p1)%node(ln1,1).eq.p1) then
            NNODE = NNODE + 1
            nodemap(p1)%node(ln1,1) = NNODE 
            nodemap(p1)%node(ln1,2) = 1
          else if (nodemap(p1)%node(ln1,1).lt.p1) then              
            p2  = nodemap(p1)%node(ln1,1)
            ln2 = nodemap(p1)%node(ln1,2)
            
            nodemap(p1)%node(ln1,1) = nodemap(p2)%node(ln2,1) 
            nodemap(p1)%node(ln1,2) = 0    
          else
            write(*,*) "Error in dof assignment"          
          endif
        enddo       
      enddo  

      write(*,*) "  Connecting Nodes : ", TNNODE , NNODE, TNNODE-NNODE   
         
      !----------------------------------------------  
      ! Get count of other mesh entities
      !----------------------------------------------                       
      NELEM  = 0
      NBOUND = 0
      NPATCH = 0
      
      do mn = 1, NMESH
        NELEM  = NELEM  + mesh(mn)%NELEM 
        NPATCH = NPATCH + mesh(mn)%NPATCH
        do bn = 1, mesh(mn)%NBOUND
          if  (mesh(mn)%bound(bn)%FACE_ID.ne.-1) then      
            NBOUND = NBOUND + 1
          endif 
        enddo   
      enddo 
             
      finalmesh%NNODE  = NNODE 
      finalmesh%NELEM  = NELEM 
      finalmesh%NBOUND = NBOUND
      finalmesh%NPATCH = NPATCH
      
      !----------------------------------------------  
      ! Allocate new mesh      
      !----------------------------------------------    
      allocate(finalmesh%xg(NNODE,NSD))
      allocate(finalmesh%wg(NNODE))
      allocate(finalmesh%IEN(NELEM,NSHL))
      allocate(finalmesh%EPID(NELEM))
      allocate(finalmesh%EIJK(NELEM,NSD))  
      
      allocate(finalmesh%bound(NBOUND))
      allocate(finalmesh%patch(NPATCH))  
      
      !----------------------------------------------    
      !  Merge meshes into new mesh (remove duplicate nodes)      
      !----------------------------------------------          
      NNODE  = 0
      NELEM  = 0
      NBOUND = 0
      NPATCH = 0
      do mn = 1, NMESH 
         write(*,*) "Merging mesh", mn
         !---------------------------------------------- 
         ! Add nodes to mesh
         !----------------------------------------------        
         write(*,*) "    Add nodes"       
         do i = 1,mesh(mn)%NNODE
           if (nodemap(mn)%node(i,2).ge.1)  then
             NNODE = NNODE + 1         
             finalmesh%xg(NNODE,:) = mesh(mn)%xg(i,:)
             finalmesh%wg(NNODE)   = mesh(mn)%wg(i)
           endif
         enddo
         
         !---------------------------------------------- 
         ! Add elements to mesh
         !----------------------------------------------   
         write(*,*) "    Add elements"                   
         do i = 1,mesh(mn)%NELEM
           NELEM = NELEM + 1
           do j = 1, NSHL
             finalmesh%IEN(NELEM,j) = nodemap(mn)%node(mesh(mn)%IEN(i,j),1)                              
           enddo
         
           finalmesh%EPID(NELEM)   = mesh(mn)%EPID(i)  + NPATCH     
           finalmesh%EIJK(NELEM,:) = mesh(mn)%EIJK(i,:) 
         enddo
         
         !---------------------------------------------- 
         ! Add boundaries to mesh
         !----------------------------------------------  
         write(*,*) "    Add boundaries"           
         do bn = 1, mesh(mn)%NBOUND
           !write(*,*) "  Local boundary ",bn,"/",mesh(mn)%NBOUND
           
           if  (mesh(mn)%bound(bn)%FACE_ID.ne.-1) then                          
             NBOUND = NBOUND + 1
             write(*,*) "    Merge boundary ",NBOUND,bn,mesh(mn)%bound(bn)%FACE_ID 
             finalmesh%bound(NBOUND)%FACE_ID = mesh(mn)%bound(bn)%FACE_ID  
             finalmesh%bound(NBOUND)%NFACE = mesh(mn)%bound(bn)%NFACE 
             finalmesh%bound(NBOUND)%NNODE = mesh(mn)%bound(bn)%NNODE
             
             allocate(finalmesh%bound(NBOUND)%FACE_IEN(mesh(mn)%bound(bn)%NFACE,NSHLB))
             allocate(finalmesh%bound(NBOUND)%F2E     (mesh(mn)%bound(bn)%NFACE))
             allocate(finalmesh%bound(NBOUND)%FACE_OR (mesh(mn)%bound(bn)%NFACE))
             allocate(finalmesh%bound(NBOUND)%BNODES  (mesh(mn)%bound(bn)%NNODE))
             
             do i = 1,mesh(mn)%bound(bn)%NFACE
               do j = 1, NSHLB
                 finalmesh%bound(NBOUND)%FACE_IEN(i,j) = nodemap(mn)%node(mesh(mn)%bound(bn)%FACE_IEN(i,j),1)
               enddo
                     
               finalmesh%bound(NBOUND)%F2E (i)    = mesh(mn)%bound(bn)%F2E(i) + NELEM - mesh(mn)%NELEM
               finalmesh%bound(NBOUND)%FACE_OR(i) = mesh(mn)%bound(bn)%FACE_OR(i)
             enddo
           
             do i = 1,mesh(mn)%bound(bn)%NNODE
               finalmesh%bound(NBOUND)%BNODES(i) = nodemap(mn)%node(mesh(mn)%bound(bn)%BNODES(i),1)  
             enddo
           
           endif       
         enddo

         !----------------------------------------------
         ! Add patches to mesh
         !----------------------------------------------       
         write(*,*) "    Add patches"     
         do pn = 1, mesh(mn)%NPATCH
           !write(*,*) "  Merge patch ",pn,"/",mesh(mn)%NPATCH
           NPATCH = NPATCH + 1
           finalmesh%patch(NPATCH)%P = mesh(mn)%patch(pn)%P
           finalmesh%patch(NPATCH)%Q = mesh(mn)%patch(pn)%Q 
           finalmesh%patch(NPATCH)%R = mesh(mn)%patch(pn)%R
           
           finalmesh%patch(NPATCH)%MCP = mesh(mn)%patch(pn)%MCP
           finalmesh%patch(NPATCH)%NCP = mesh(mn)%patch(pn)%NCP
           finalmesh%patch(NPATCH)%OCP = mesh(mn)%patch(pn)%OCP
           
           allocate(finalmesh%patch(NPATCH)%U_KNOT(mesh(mn)%patch(pn)%MCP + mesh(mn)%patch(pn)%P + 1))
           allocate(finalmesh%patch(NPATCH)%V_KNOT(mesh(mn)%patch(pn)%NCP + mesh(mn)%patch(pn)%Q + 1))
           allocate(finalmesh%patch(NPATCH)%W_KNOT(mesh(mn)%patch(pn)%OCP + mesh(mn)%patch(pn)%R + 1))
           
           finalmesh%patch(NPATCH)%U_KNOT = mesh(mn)%patch(pn)%U_KNOT
           finalmesh%patch(NPATCH)%V_KNOT = mesh(mn)%patch(pn)%V_KNOT
           finalmesh%patch(NPATCH)%W_KNOT = mesh(mn)%patch(pn)%W_KNOT
         enddo
                                                 
      enddo      
      
      finalmesh%NNODE  = NNODE 
      finalmesh%NELEM  = NELEM 
      finalmesh%NBOUND = NBOUND
      finalmesh%NPATCH = NPATCH

return
end subroutine mergeMeshes
!------------------------------------------------------------------------
subroutine rearrangeBoundaries(mesh,eps)

  use class
  
  implicit none

  type (mesh_class) mesh
  type (bnd_class), allocatable :: bound(:)  
  integer :: NFACE(7), b,f,e,j,NBOUND ,b2,n
  real(8) :: xl(mesh%NSHLB,mesh%NSD),eps
  logical :: flag(mesh%NNODE)
  real(8) :: gmax(3), gmin(3),lmax(3), lmin(3)
    
 
  NFACE = 0     
  gmax = maxval(mesh%xg,1)
  gmin = minval(mesh%xg,1)
  do b = 1, mesh%NBOUND
    !bounds(b)
    do f = 1, mesh%bound(b)%NFACE      
      do j = 1, mesh%NSHLB
        xl(j,:) = mesh%xg(mesh%bound(b)%FACE_IEN(f,j),:)
      enddo
      
      lmax = maxval(xl,1)
      lmin = minval(xl,1) 
       
      if      (lmax(1).lt.(gmin(1)+eps)) then
          b2 = 5  
      else if (lmin(1).gt.(gmax(1)-eps)) then       
          b2 = 3 
      else if (lmax(2).lt.(gmin(2)+eps)) then
          b2 = 2        
      else if (lmin(2).gt.(gmax(2)-eps)) then        
          b2 = 4    
      else if (lmax(3).lt.(gmin(3)+eps)) then             
          b2 = 1        
      else if (lmin(3).gt.(gmax(3)-eps)) then         
          b2 = 6      
      else
        b2 = 7         
      endif
     
      NFACE(b2) = NFACE(b2) + 1
    enddo        
  enddo
  write(*,*) "NFACE "     
  write(*,*) NFACE   

!... Allocate boundaries  
  if (NFACE(7).gt.0) then
    NBOUND = 7
  else  
    NBOUND = 6
  endif
  
  allocate(bound(NBOUND))
  do b = 1, NBOUND
    bound(b)%FACE_ID=b
    bound(b)%NFACE=NFACE(b)
    allocate(bound(b)%FACE_IEN(NFACE(b),mesh%NSHLB))    
    allocate(bound(b)%F2E(NFACE(b)))
    allocate(bound(b)%FACE_OR(NFACE(b)))    
  enddo

!... Asign faces to boundaries    
  NFACE = 0
  do b = 1, mesh%NBOUND
    do f = 1, mesh%bound(b)%NFACE      
      do j = 1, mesh%NSHLB
        xl(j,:) = mesh%xg(mesh%bound(b)%FACE_IEN(f,j),:)
      enddo
      lmax = maxval(xl,1)
      lmin = minval(xl,1) 
       
      if      (lmax(1).lt.(gmin(1)+eps)) then
          b2 = 5  
      else if (lmin(1).gt.(gmax(1)-eps)) then       
          b2 = 3 
      else if (lmax(2).lt.(gmin(2)+eps)) then
          b2 = 2        
      else if (lmin(2).gt.(gmax(2)-eps)) then        
          b2 = 4    
      else if (lmax(3).lt.(gmin(3)+eps)) then             
          b2 = 1        
      else if (lmin(3).gt.(gmax(3)-eps)) then         
          b2 = 6      
      else
          b2 = 7         
      endif
     
      NFACE(b2) = NFACE(b2) + 1

      bound(b2)%FACE_IEN(NFACE(b2),:) = mesh%bound(b)%FACE_IEN(f,:)
      bound(b2)%F2E     (NFACE(b2))   = mesh%bound(b)%F2E    (f)  
      bound(b2)%FACE_OR (NFACE(b2))   = mesh%bound(b)%FACE_OR(f)

    enddo        
  enddo  
  write(*,*) "NFACE "     
  write(*,*) NFACE  

!... Asign nodes to boundaries      

  do b = 1, NBOUND 
    bound(b)%NNODE=0
    flag = .true. 
    do f = 1,bound(b)%NFACE         
      do j = 1, mesh%NSHLB
        n = bound(b)%FACE_IEN(f,j)
        if (flag (n)) then
          flag (n) = .false.
          bound(b)%NNODE = bound(b)%NNODE + 1
        endif
      enddo
    enddo
    
    allocate(bound(b)%BNODES(bound(b)%NNODE))
    
    bound(b)%NNODE=0    
    flag = .true. 
    do f = 1,bound(b)%NFACE         
      do j = 1, mesh%NSHLB
        n = bound(b)%FACE_IEN(f,j)
        if (flag (n)) then
          flag (n) = .false.
          bound(b)%NNODE = bound(b)%NNODE + 1
          bound(b)%BNODES(bound(b)%NNODE) = n
        endif
      enddo
    enddo                
  enddo
   
  ! Dealloc 
  do b = 1, mesh%NBOUND
    deallocate(mesh%bound(b)%FACE_IEN)    
    deallocate(mesh%bound(b)%F2E)
    deallocate(mesh%bound(b)%FACE_OR)   
    deallocate(mesh%bound(b)%BNODES) 
  enddo   
  deallocate(mesh%bound)
  
  ! Copy
  mesh%NBOUND=NBOUND
  allocate(mesh%bound(NBOUND))
  mesh%bound = bound
     
return
end 
