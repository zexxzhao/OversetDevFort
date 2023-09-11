!------------------------------------------------------------------------
! Convert a multipatch NURBS mesh in linear/quad hex mesh                                                                                                                                                           
!------------------------------------------------------------------------
subroutine iso2hex(omesh,eps)

  use class 
  implicit none
  
  type (mesh_class) omesh
  type (mesh_class) imesh
  
  integer interp,pn,ni,nj,nk,i,j,k,e1,e2,ln1,ln2,gn,n,p1,p2, NGAUSS,found,ip  , gfid
  real(8) :: xl(omesh%NSHLmax,omesh%NSD),wl(omesh%NSHLmax),xi(omesh%NSD),dx(omesh%NSD),shlu(omesh%NSHLmax)
  real(8), allocatable :: gp(:,:)
  real(8) :: eps,eps2
  integer :: maxi, maxj,maxk 
  integer, allocatable :: IJKE(:,:,:,:)
  integer :: INE(omesh%NNODE,18*omesh%NSHLmax),cnt1(omesh%NELEM*27) 
  integer :: E2E(omesh%NELEM,7**omesh%NSD-1),cnt(omesh%NELEM) ,cnt2(omesh%NELEM) 
  integer :: E2Ew(omesh%NELEM,7**omesh%NSD-1)
  logical :: flag(omesh%NELEM)
  logical,allocatable :: flag2(:)
  character*30 fname
  character*30 inp
  integer b,f,o,e, i2f(6,8)

!
!... Invert EIJK to get IJKE 
!
  maxi = 0
  maxj = 0
  maxk = 0
  do pn = 1, omesh%NPATCH
    maxi = max(maxi, omesh%patch(pn)%MCP+1)
    maxj = max(maxj, omesh%patch(pn)%NCP+1)
    maxk = max(maxk, omesh%patch(pn)%OCP+1) 
  enddo
  
  write(*,*) "allocate IJKE array"
  write(*,*)omesh%NPATCH,maxi,maxj,maxk
  
  allocate(IJKE(omesh%NPATCH,maxi,maxj,maxk))
  IJKE = -1
  do e1 = 1, omesh%NELEM
    pn = omesh%EPID(e1)
    ni = omesh%EIJK(e1,1)
    nj = omesh%EIJK(e1,2)
    nk = omesh%EIJK(e1,3)   
    
    IJKE(pn,ni,nj,nk) = e1
  enddo

!
!... Use IJKE to generate the element to element (E2E) connectivity graph
!    For now the assumption is  that there are no cp-2 lines
!  
  write(*,*) "Generate element 2 element connectivity"     
  cnt = 0
  
  do e1 = 1, omesh%NELEM 
    pn = omesh%EPID(e1)
    ni = omesh%EIJK(e1,1)
    nj = omesh%EIJK(e1,2)
    nk = omesh%EIJK(e1,3)   
  
    do i = -1,1
      do j = -1,1
        do k = -1,1
          e2 = IJKE(pn,ni+i,nj+j,nk+k)
          if ((e2.ge.0).and.(e1.ne.e2)) then
            cnt(e1) = cnt(e1) + 1
            E2E(e1,cnt(e1)) = e2  
          endif
        enddo
      enddo
    enddo        
  enddo   
  
  write(*,*) "  Max element 2 element connectivity", maxval(cnt) 
  E2Ew = E2E    ! Store the patchwise E2E for later
  cnt2 = cnt
!
!... Invert relevant part of IEN to INE to get node to element map
!      
  
  cnt1 = 0
  do e1 = 1, omesh%NELEM
    if (cnt(e1).lt.26) then
      do i = 1, omesh%NSHLmax
        n = omesh%IEN(e1,i)
        cnt1(n) = cnt1(n)+ 1 
        INE(n,cnt1(n)) = e1
      enddo
    endif
  enddo
  
!
!... Use INE to generate the remaining part of element to element (E2E) connectivity graph
!  

  do e1 = 1, omesh%NELEM 
    if (cnt(e1).lt.26) then
      !p1 = omesh%EPID(e1)
    
      flag = .true.
      do i = 1, cnt(e1)
        e2 = E2E(e1,i) 
        flag(e2) = .false.
      enddo
      
      do i = 1, omesh%NSHLmax
        n = omesh%IEN(e1,i) 
        do j = 1, cnt1(n) 
          e2 = INE(n,j)
          !p2 = omesh%EPID(e2)                    
 
          !if ((p1.ne.p2).and.(flag(e2))) then           
           if (flag(e2)) then
            cnt(e1) = cnt(e1) + 1
            E2E(e1,cnt(e1)) = e2   
            flag(e2) = .false.  
          endif 

        enddo
      enddo      
  
    endif                    
  enddo  
   write(*,*) "  Max element 2 element connectivity", maxval(cnt)    

!
!... Loop over interpolation meshes
!   
  do interp = 1, 0, -1
!
!... Set interpolation
!      
    if (interp .eq.0) then
      NGAUSS      = 8    
      imesh%NSHLmax  = 8
      imesh%NSHLB = 4    
           
      i2f(1,1) = 1 
      i2f(1,2) = 2
      i2f(1,3) = 4   
      i2f(1,4) = 3

      i2f(2,1) = 1 
      i2f(2,2) = 2
      i2f(2,3) = 6   
      i2f(2,4) = 5

      i2f(3,1) = 2  
      i2f(3,2) = 4 
      i2f(3,3) = 8   
      i2f(3,4) = 6 
      
      i2f(4,1) = 3 
      i2f(4,2) = 4
      i2f(4,3) = 8   
      i2f(4,4) = 7  
      
      i2f(5,1) = 1 
      i2f(5,2) = 3
      i2f(5,3) = 7       
      i2f(5,4) = 5 
      
      i2f(6,1) = 5
      i2f(6,2) = 6   
      i2f(6,3) = 8      
      i2f(6,4) = 7      
                  
    else if (interp .eq.1) then     
      NGAUSS      = 20      
      imesh%NSHLmax  = 20
      imesh%NSHLB = 4      
       
      i2f(1,1) = 1 
      i2f(1,2) = 2
      i2f(1,3) = 3   
      i2f(1,4) = 4

      i2f(2,1) = 1 
      i2f(2,2) = 2
      i2f(2,3) = 6   
      i2f(2,4) = 5

      i2f(3,1) = 2  
      i2f(3,2) = 3 
      i2f(3,3) = 7   
      i2f(3,4) = 6 
      
      i2f(4,1) = 4 
      i2f(4,2) = 3
      i2f(4,3) = 7   
      i2f(4,4) = 8  
      
      i2f(5,1) = 1 
      i2f(5,2) = 4
      i2f(5,3) = 8       
      i2f(5,4) = 5 
      
      i2f(6,1) = 5
      i2f(6,2) = 6   
      i2f(6,3) = 7      
      i2f(6,4) = 8        
    else
      write(*,*) 'Wrong interpolation_mode: ', interp     
      write(*,*) ' choose:  0 / 1'   
      stop
    endif 
    allocate(gp(NGAUSS,omesh%NSD))
    call getInterpolationPoints(gp,NGAUSS,omesh%NSD)     
 !
!... init interpolation mesh
!           
    imesh%NSD    = omesh%NSD   
    imesh%NELEM  = omesh%NELEM 
    imesh%NNODE  = imesh%NELEM*imesh%NSHLmax
    imesh%NBOUND = omesh%NBOUND
    imesh%NPATCH = 0      
 
    allocate(imesh%xg(imesh%NNODE,imesh%NSD))
    allocate(imesh%IEN(imesh%NELEM,imesh%NSHLmax)) 
    allocate(imesh%bound(imesh%NBOUND))
 
    do i = 1, imesh%NBOUND
      imesh%bound(i)%NFACE = omesh%bound(i)%NFACE
      imesh%bound(i)%NNODE = imesh%bound(i)%NFACE*imesh%NSHLB
      allocate(imesh%bound(i)%FACE_IEN(imesh%bound(i)%NFACE,imesh%NSHLB))
      allocate(imesh%bound(i)%F2E(imesh%bound(i)%NFACE))
      allocate(imesh%bound(i)%FACE_OR(imesh%bound(i)%NFACE))
    
      allocate(imesh%bound(i)%BNODES(imesh%bound(i)%NNODE))
    enddo
  
!
!... Loop over element
!    Interpolate values and create IEN on the fly
!        
    flag = .false.
    imesh%IEN  = -1
    imesh%NNODE = 0
    eps2 = eps**2
    do e1 = 1, omesh%NELEM
      if (mod(e1,50000).eq.0) write(*,*) e1,"/",omesh%NELEM

      pn = omesh%EPID(e1)
      ni = omesh%EIJK(e1,1)
      nj = omesh%EIJK(e1,2)
      nk = omesh%EIJK(e1,3)  
   
      do j = 1, omesh%NSHLmax 
        xl(j,:) = omesh%xg(omesh%IEN(e1,j),:)
        wl(j)   = omesh%wg(omesh%IEN(e1,j))
      enddo

! Interpolate values   
      do ln1 = 1, NGAUSS
        call eval_shape_nurbs(gp(ln1,:),ni,nj,nk,omesh%patch(pn),xl,wl,shlu, omesh%NSD, omesh%NSHL)

        do j = 1, omesh%NSD
          xi(j) = sum(xl(:,j)*shlu) 
        enddo
 
! Find existing nodes in neighbor elements
        do j = 1, cnt(e1)
          e2 = E2E(e1,j)
          if (flag(e2)) then
            do ln2 = 1, NGAUSS
              gn = imesh%IEN(e2,ln2)
    
              dx = imesh%xg(gn,:) - xi
              if (sum(dx*dx).le.eps2) then
                imesh%IEN(e1,ln1)= gn            
                exit
              endif 
            enddo    
          endif
          if (imesh%IEN(e1,ln1).eq.gn) exit
        enddo 
    
! Add node if not found  
        if (imesh%IEN(e1,ln1).eq.-1) then
          imesh%NNODE = imesh%NNODE + 1
          imesh%IEN(e1,ln1) = imesh%NNODE
          imesh%xg(imesh%NNODE,:) = xi
        endif
      enddo 
  
      flag(e1)= .true.
    enddo
     
!
!...  Translate faces/nodes
!   
    allocate(flag2(imesh%NNODE))
    do b = 1, imesh%NBOUND 
      imesh%bound(b)%FACE_ID = omesh%bound(b)%FACE_ID 
      imesh%bound(b)%F2E     = omesh%bound(b)%F2E
      imesh%bound(b)%FACE_OR = omesh%bound(b)%FACE_OR 
      imesh%bound(b)%NNODE=0
      flag2 = .true. 
      do f = 1,imesh%bound(b)%NFACE         
        e = imesh%bound(b)%F2E(f)
        o = imesh%bound(b)%FACE_OR(f)
        do j = 1, imesh%NSHLB
          n = imesh%IEN(e, i2f(o,j))
          imesh%bound(b)%FACE_IEN(f,j) = n
          !write(*,*) b,f,e,o
          !write(*,*) n, omesh%NNODE,imesh%NNODE
          if (flag2 (n)) then
            flag2 (n) = .false.
            imesh%bound(b)%NNODE = imesh%bound(b)%NNODE + 1
            imesh%bound(b)%BNODES(imesh%bound(b)%NNODE) = n
          endif
        enddo
      enddo
                  
    enddo
    deallocate(flag2)  
!
!...  extract boundary nodes
!    
  
!
!... Write interpolation mesh
!    
    write(inp,'(I8)') interp
    fname = 'hex.' // trim(adjustl(inp)) //'.dat'      
    call writeMesh(imesh , fname)   
  
!
!... Write VTK file of interpolation mesh 
!      
    write(inp,'(I8)') interp
    fname = 'mesh.' // trim(adjustl(inp)) //'.vtk'    
    call writeMeshVTK(imesh, fname)
    
    do b = 1, imesh%NBOUND 
      !write(*,*) b, imesh%bound(b)%Face_ID!, imesh%bound(b)%Face_OR(1)
      write(inp,'(I8)') interp
      fname = 'bnd.' // trim(adjustl(inp))
      write(inp,'(I8)') b 
      fname = trim(fname) //'.' // trim(adjustl(inp)) //'.vtk'  
      call  writeBndMeshVTK(imesh,fname,b)   
    enddo
    
!
!... Deallocate interpolation mesh and Gauss points
!
    if (interp.eq.1) call deallocMesh(imesh)
    deallocate(gp)
  enddo
  
!
!... Write IEN of linear hex for partitioning
! 
 call writeIEN  (imesh,'ien.dat')  

!
!... Invert relevant part of IEN to INE to get node to element map
!      
write(*,*)  imesh%NNODE, omesh%NNODE  
  cnt1 = 0
  do e1 = 1, imesh%NELEM
    if (cnt2(e1).lt.26) then
      do i = 1, imesh%NSHLmax
        n = imesh%IEN(e1,i)
        cnt1(n) = cnt1(n)+ 1 
        INE(n,cnt1(n)) = e1
      enddo
    endif
  enddo
    
!
!... Use INE to generate the element to element (E2E) connectivity graph
!  
  E2E = E2Ew    ! restore the patchwise E2E for later
  do e1 = 1, omesh%NELEM 
    if (cnt2(e1).lt.26) then
      p1 = omesh%EPID(e1)
    
      flag = .true.
      do i = 1, imesh%NSHLmax
        n = imesh%IEN(e1,i) 
        do j = 1, cnt1(n) 
          e2 = INE(n,j)
          p2 = omesh%EPID(e2)                    
 
          if ((p1.ne.p2).and.(flag(e2))) then           
            cnt2(e1) = cnt2(e1) + 1
            E2E(e1,cnt2(e1)) = e2   
            flag(e2) = .false.  
          endif 

        enddo
      enddo      
  
    endif                    
  enddo  
  write(*,*) "  Max element 2 element connectivity", maxval(cnt2)  
   
!
!... Get graph weight
!   
  allocate(flag2(omesh%NNODE))   
  do e1 = 1, omesh%NELEM
    flag2 = .false.
    do i = 1, omesh%NSHLmax
      n = omesh%IEN(e1,i) 
      flag2(n) = .true.
    enddo
        
    do j=1,cnt2(e1)
       e2 = E2E(e1,j)
       found = 0
       do i = 1, omesh%NSHLmax
         n = omesh%IEN(e2,i) 
         if (flag2(n)) found = found + 1
       enddo          
       E2Ew(e1,j) = found
    enddo     
  enddo 
  deallocate(flag2)       
!
!... write graph
!     
  write(*,*) "Write element graph" 
  gfid = 55
  open (gfid, file='graph.dat', status='unknown')  
  write (gfid, *) omesh%NELEM, sum(cnt2)/2,1
  do e1 = 1, omesh%NELEM 
    write(gfid,'(256I8)') (E2E(e1,j),E2Ew(e1,j),j=1,cnt2(e1))
  enddo    
  close(gfid) 







end subroutine  
