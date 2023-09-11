      subroutine input_NURBS(mesh,fname, net)

      use class

      implicit none
      
      type (mesh_class) mesh
      character*(*) fname
      logical net
      
      integer NSD,NNODE,NELEM,NSHL,NSHLB
      integer :: P,Q,R
      integer :: MCP,NCP,OCP
      real(8), allocatable :: B_NET(:,:,:,:)
         
      integer :: i, j, k, l, meshf
      
      integer :: loop1,loop2,loop3,g, e, gtemp, ln,fn
      integer :: iel,face(6),node(6)

      write(*,*) "Open NURBS file:", fname           
      meshf = 11      
 
      allocate(mesh%patch(1))
      allocate(mesh%bound(6))
      
!     
!...  Read in preliminary information
!     
      open (meshf, file=fname, status='old')
      read (meshf,*)  NSD
      read (meshf,*)  P, Q, R      
      read (meshf,*)  MCP, NCP, OCP
      
      mesh%NSD = NSD      
      mesh%patch(1)%P = P
      mesh%patch(1)%Q = Q
      mesh%patch(1)%R = R
      mesh%patch(1)%MCP = MCP
      mesh%patch(1)%NCP = NCP
      mesh%patch(1)%OCP = OCP
      
!     
!...  Allocate arrays for knot vectors and for control net
!     
      allocate(mesh%patch(1)%U_KNOT(MCP + P + 1))    
      allocate(mesh%patch(1)%V_KNOT(NCP + Q + 1))
      allocate(mesh%patch(1)%W_KNOT(OCP + R + 1))
      allocate(B_NET(MCP,NCP,OCP,NSD+1))
      
!     
!...  Read knot vectors and control net
!
      read(meshf,*) (mesh%patch(1)%U_KNOT(i), i=1,MCP + P + 1)
      read(meshf,*) (mesh%patch(1)%V_KNOT(i), i=1,NCP + Q + 1)     
      read(meshf,*) (mesh%patch(1)%W_KNOT(i), i=1,OCP + R + 1)
      
      do k = 1,OCP
         do j = 1,NCP
            do i = 1,MCP                   
               read (meshf,*) (B_NET(i,j,k,l), l = 1,NSD+1)
            enddo
         enddo
      enddo

!
!... Reduce P,Q,R such that the NET remains
!    
      if (net) then 
        call knot2net(mesh%patch(1)%U_KNOT,MCP,P) 
        call knot2net(mesh%patch(1)%V_KNOT,NCP,Q)
        call knot2net(mesh%patch(1)%W_KNOT,OCP,R)        
        
        mesh%patch(1)%P = P
        mesh%patch(1)%Q = Q
        mesh%patch(1)%R = R              
      endif                   
!     
!...  Derived quantities
!      
      NNODE = MCP*NCP*OCP      
      NELEM = (MCP-P)*(NCP-Q)*(OCP-R)      
      NSHL  = (P+1)*(Q+1)*(R+1)
      NSHLB = (P+1)*(Q+1)
      
      mesh%NNODE  = NNODE
      mesh%NELEM  = NELEM
      mesh%NSHL   = NSHL
      mesh%NSHLB  = NSHLB
      mesh%NBOUND = 6
      mesh%NPATCH = 1      
!     
!...  Generate IEN,EPID,EIJK
!             
      allocate(mesh%IEN (NELEM,NSHL))       
      allocate(mesh%EPID(NELEM))      
      allocate(mesh%EIJK(NELEM,3))
      allocate(mesh%xg(NNODE,NSD))  
      allocate(mesh%wg(NNODE))  
!     
!...  Initialize matrices and variables
!     
      g = 0
      e = 0
      mesh%EPID = 1
!     
!...  Loop through control points assigning global node
!     numbers and filling out IEN and INN as we go
      do k = 1,OCP            
        do j = 1,NCP          
          do i = 1,MCP        
            g = g+1
            
            mesh%xg(g,:) = B_NET(i,j,k,1:3)
            mesh%wg(g)   = B_NET(i,j,k,4)
            
            if ( (mesh%patch(1)%U_KNOT(i).ne.mesh%patch(1)%U_KNOT(i+1)).and. &
                 (mesh%patch(1)%V_KNOT(j).ne.mesh%patch(1)%V_KNOT(j+1)).and. &
                 (mesh%patch(1)%W_KNOT(k).ne.mesh%patch(1)%W_KNOT(k+1))) then
        
               e = e +1
               
               mesh%EIJK(e,1) = i
               mesh%EIJK(e,2) = j
               mesh%EIJK(e,3) = k  
               
               do loop1 = 0,R
                  do loop2 = 0,Q
                     do loop3 = 0,P
                        gtemp = g - loop1*MCP*NCP -MCP*loop2-loop3
                        ln = (P+1)*(Q+1)*loop1+(P+1)*loop2+loop3+1
                        mesh%IEN(e,ln) = gtemp
                     enddo
                  enddo
               enddo
                           
            endif
          enddo
        enddo
      enddo   
      
      NELEM = e
      mesh%NELEM  = NELEM

!
!...  Generate Faces
!
      mesh%bound(1)%NFACE = (MCP-P)*(NCP-Q)            
      mesh%bound(2)%NFACE = (MCP-P)*(OCP-R)        
      mesh%bound(3)%NFACE = (NCP-Q)*(OCP-R)    
      mesh%bound(4)%NFACE = (MCP-P)*(OCP-R)
      mesh%bound(5)%NFACE = (NCP-Q)*(OCP-R)         
      mesh%bound(6)%NFACE = (MCP-P)*(NCP-Q)
           
      do i = 1,6
        allocate(mesh%bound(i)%FACE_IEN(mesh%bound(i)%NFACE,NSHLB))
        allocate(mesh%bound(i)%F2E(mesh%bound(i)%NFACE))   
        allocate(mesh%bound(i)%FACE_OR(mesh%bound(i)%NFACE)) 
        mesh%bound(i)%FACE_ID = i
      enddo     
   
!       Face Orientation scheme: 
!           1 - surface (u,v,1)
!           2 - surface (u,1,w)
!           3 - surface (MCP,v,w)
!           4 - surface (u,NCP,w)
!           5 - surface (1,v,w)
!           6 - surface (u,v,OCP)

      face = 0
      
      do iel=1, NELEM
         i = mesh%EIJK(iel,1)
         j = mesh%EIJK(iel,2)
         k = mesh%EIJK(iel,3)
         
         if (k.eq.(R+1)) then
           face(1) = face(1) + 1              
           mesh%bound(1)%F2E(face(1)) = iel
           mesh%bound(1)%FACE_OR = 1
                     
           fn = 0
           loop1 = R
           do loop2 = 0,Q
             do loop3 = 0,P            
                 ln = (P+1)*(Q+1)*loop1+(P+1)*loop2+loop3+1
                 fn = fn + 1
                 mesh%bound(1)%FACE_IEN(face(1),fn) = mesh%IEN(iel,ln) 
             enddo
           enddo
                         
               
         endif

         if (j.eq.(Q+1)) then
           face(2) = face(2) + 1             
           mesh%bound(2)%F2E(face(2)) = iel
           mesh%bound(2)%FACE_OR = 2

           fn = 0
           loop2 = Q
           do loop1 = 0,R
             do loop3 = 0,P                    
                 ln = (P+1)*(Q+1)*loop1+(P+1)*loop2+loop3+1
                 fn = fn + 1
                 mesh%bound(2)%FACE_IEN(face(2),fn) = mesh%IEN(iel,ln) 
             enddo
           enddo
         endif

         if (i.eq.(MCP)) then 
           face(3) = face(3) + 1             
           mesh%bound(3)%F2E(face(3)) = iel
           mesh%bound(3)%FACE_OR = 3

           fn = 0
           loop3 = 0
           do loop1 = 0,R
             do loop2 = 0,Q                    
                 ln = (P+1)*(Q+1)*loop1+(P+1)*loop2+loop3+1
                 fn = fn + 1
                 mesh%bound(3)%FACE_IEN(face(3),fn) = mesh%IEN(iel,ln) 
             enddo
           enddo
         
         endif

         if (j.eq.(NCP)) then
           face(4) = face(4) + 1             
           mesh%bound(4)%F2E(face(4)) = iel
           mesh%bound(4)%FACE_OR = 4

           fn = 0
           loop2 = 0
           do loop1 = 0,R
             do loop3 = 0,P                  
                 ln = (P+1)*(Q+1)*loop1+(P+1)*loop2+loop3+1
                 fn = fn + 1
                 mesh%bound(4)%FACE_IEN(face(4),fn) = mesh%IEN(iel,ln) 
             enddo
           enddo
         
         endif

         if (i.eq.(P+1)) then
           face(5) = face(5) + 1             
           mesh%bound(5)%F2E(face(5)) = iel
           mesh%bound(5)%FACE_OR = 5

           fn = 0             
           loop3 = P 
           do loop1 = 0,R
             do loop2 = 0,Q                   
                 ln = (P+1)*(Q+1)*loop1+(P+1)*loop2+loop3+1
                 fn = fn + 1
                 mesh%bound(5)%FACE_IEN(face(5),fn) = mesh%IEN(iel,ln) 
             enddo
           enddo
         
         endif
         
         if (k.eq.(OCP)) then
           face(6) = face(6) + 1             
           mesh%bound(6)%F2E(face(6)) = iel
           mesh%bound(6)%FACE_OR = 6

           fn = 0
           loop1 = 0
           do loop2 = 0,Q
             do loop3 = 0,P                    
                 ln = (P+1)*(Q+1)*loop1+(P+1)*loop2+loop3+1
                 fn = fn + 1
                 mesh%bound(6)%FACE_IEN(face(6),fn) = mesh%IEN(iel,ln) 
             enddo
           enddo
         
         endif
      
      enddo
      do i = 1,6
        if (mesh%bound(i)%NFACE .ne.face(i)) then
          mesh%bound(i)%NFACE = face(i)
        endif
      enddo         
      
!
!...  Generate Face nodes
!      
      mesh%bound(1)%NNODE = MCP*NCP            
      mesh%bound(2)%NNODE = MCP*OCP        
      mesh%bound(3)%NNODE = NCP*OCP    
      mesh%bound(4)%NNODE = MCP*OCP
      mesh%bound(5)%NNODE = NCP*OCP         
      mesh%bound(6)%NNODE = MCP*NCP    
         
      do i = 1,6
        allocate(mesh%bound(i)%BNODES(mesh%bound(i)%NNODE))
        mesh%bound(i)%FACE_OR = i
      enddo      
      
      g    = 0
      node = 0
!     
!...  Loop through control points assigning global node
!     numbers and filling out IEN and INN as we go
      do k = 1,OCP            
        do j = 1,NCP          
          do i = 1,MCP        
            g = g+1
            
            if (k.eq.1) then
              node(1) = node(1) + 1              
              mesh%bound(1)%BNODES(node(1)) = g         
            endif

            if (j.eq.1) then
              node(2) = node(2) + 1             
              mesh%bound(2)%BNODES(node(2)) = g
            endif

            if (i.eq.MCP) then 
              node(3) = node(3) + 1          
              mesh%bound(3)%BNODES(node(3)) = g
            endif

            if (j.eq.NCP) then
              node(4) = node(4) + 1          
              mesh%bound(4)%BNODES(node(4)) = g
            endif
   
            if (i.eq.1) then
              node(5) = node(5) + 1          
              mesh%bound(5)%BNODES(node(5)) = g
            endif
         
            if (k.eq.OCP) then
              node(6) = node(6) + 1          
              mesh%bound(6)%BNODES(node(6)) = g
            endif            
            
          enddo
        enddo
      enddo 
!
!...  Face check
!      
      !write(*,*) "Checking boundary nodes found.... "          
      do i = 1,6
        if (mesh%bound(i)%NNODE .ne.node(i)) then
          write(*,*) "Boundary = ",i
          write(*,*) "Nodes    = ",mesh%bound(i)%NNODE, node(i)
        endif        
      enddo 
                                 
      return
      end      
      
      
!============================================================
!
!============================================================

      subroutine   knot2net(KNOT,MCP,P)       
      
      implicit none
      
      integer :: MCP,P,i,j
      real(8) :: KNOT(MCP+P+1)
      
      if ( P.ge.2) then     
        write(*,'(128ES16.8)') KNOT
        
        KNOT(1:2)  = 0d0
             
        do i = 3,MCP
          KNOT(i) = real(i-2)/real(MCP-3)         
        enddo 
        
        KNOT(MCP + 1:MCP + 2) = 1d0          
        P = 1
        write(*,'(128ES16.8)') KNOT
      endif
            
      return
      end       
      
      
