program geomData

  implicit none
      
  integer,parameter :: NSHL  = 4
  integer,parameter :: NSHLb = 3
  integer,parameter :: NGAUSS = 4 
  integer,parameter :: NSD =3  
      
  integer :: meshf, NNode, NElem,i,j,itmp,iel,igauss, NFace, NEdge
  integer :: nmap(4)
  real(8), allocatable :: xg(:,:)
  integer, allocatable :: IEN(:,:)
  real(8) :: shlu(NSHL),&
          shgradlu(NSHL,NSD),&
          shgradgu(NSHL,NSD)
  real(8) :: water_level,xl(NSHL,NSD),xi(NSD),DetJ
  real(8) :: gp(NGAUSS, NSD), gw(NGAUSS), rho, xcg(NSD)
  real(8) :: masse,vole,Qe(NSD),Ie(NSD,NSD)  
  real(8) :: mass,vol,Qt(NSD),It(NSD,NSD) 
  
  character(len=30) :: fname
  character(len=40) :: ctmp
  
  nmap(1)=1
  nmap(2)=2
  nmap(3)=4
  nmap(4)=3
  
!---------------------------------
! Read input
!---------------------------------
  call getarg(1,fname)  
  call getarg(2,ctmp)
  
  read(ctmp, *) water_level 

!---------------------------------
! Read mesh
!---------------------------------


  meshf = 11   
  open(meshf, file = trim(fname), status = 'old')
  write(*,*) "Reading ", fname   
  
  read(meshf,*) ctmp, itmp  
  read(meshf,*) ctmp
  read(meshf,*) itmp  !NSD
  
  ! Vertices
  read(meshf,*) ctmp  
  read(meshf,*) NNode
  write(*,*) ctmp, " ", NNode
  allocate(xg(NNode,NSD))
  do i = 1, NNode
    read(meshf,*) ( xg(i,j), j = 1, NSD), itmp    
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
    do i = 1, NFace
      read(meshf,*) itmp,itmp,itmp, itmp
    end do 
    read(meshf,*) ctmp
  endif
    
  ! Tetrahedra  
  read(meshf,*) NElem
  write(*,*) ctmp, " ", NElem  
  allocate(IEN(NElem,NSHL)) 
  do i = 1,NElem
    read(meshf,*) IEN(i,1),IEN(i,2),IEN(i,4),IEN(i,3), itmp
  end do 
      
!---------------------------------
! Loop over elements  
!---------------------------------  
  mass = 0d0
  vol  = 0d0  
  Qt   = 0d0   
  
  call genGPandGW(gp,gw, NGAUSS)
  
  do iel = 1, NElem
    do i = 1, NSHL
      xl(i,:) = xg(IEN(iel,i),:)
    end do
    masse = 0d0
    vole  = 0d0  
    Qe    = 0d0   
    do igauss = 1, NGAUSS
  
      call eval_SHAPE_3D(NSHL, NSD,gp(igauss,:),xl, &
         shlu,shgradgu,DetJ)
      
      do i=1,NSD   
        xi(i) = sum(xl(:,i)*shlu)   
      enddo

      if (xi(3).gt.water_level) then
        rho = 1d0
      else
        rho = 1d3
      endif  

      masse = masse + rho*DetJ*gw(igauss) 
      vole  = vole  + DetJ*gw(igauss)
      Qe    = Qe    + xi*DetJ*gw(igauss) 

         
    enddo ! gauss     
    
    mass = mass + masse
    vol  = vol  + vole 
    Qt   = Qt   + Qe      
  enddo ! elem  
  
  xcg = Qt/vol 
  It = 0d0
  do iel = 1, NElem
    do i = 1, NSHL
      xl(i,:) = xg(IEN(iel,i),:)
    end do
    Ie    = 0d0
    do igauss = 1, NGAUSS
  
      call eval_SHAPE_3D(NSHL, NSD,gp(igauss,:),xl, &
         shlu,shgradgu,DetJ)
      
      do i=1,NSD   
        xi(i) = sum(xl(:,i)*shlu) - xcg(i)  
      enddo
 

      do i=1,NSD
        Ie(i,i)   = Ie(i,i) + sum(xi*xi)*DetJ*gw(igauss)   
        do j=1,NSD    
          Ie(i,j)   = Ie(i,j) - xi(i)*xi(j)*DetJ*gw(igauss)      
        enddo 
      enddo   
         
    enddo ! gauss     
       
    It   = It   + Ie
         
  enddo ! elem  
    
  write(*,*)  "Vol  = ", vol 
  write(*,*)  "Mass = ", mass
  write(*,*) "--------------------------------------------"
  write(*,*)  "Xcg  = ", xcg 
  write(*,*) "--------------------------------------------"
  write(*,*) " I = "
  write(*,*)  It  
  write(*,*) "--------------------------------------------"
 end program  geomData  
 
 
 
 
       
      subroutine eval_SHAPE_3D(NSHLu,NSD,gp,xl,shlu,shgradgu,DetJ)


      implicit none
      
      ! Element number
      integer :: NSHLu, NSD

      ! u and v coordinates of integration point in parent element
      real(8) :: gp(NSD)
      
      ! Vector of Local basis function values at (u_hat, v_hat), local and 
      ! global gradients.
      real(8) :: xl(NSHLu,NSD)      
      real(8) :: shlu(NSHLu),shgradlu(NSHLu,NSD),shgradgu(NSHLu,NSD)
      real(8) :: dxdxi(NSD,NSD), dxidx(NSD,NSD),DetJ     
      
      ! NURBS coordinates, counters for loops
      integer :: niu, nju, nku, nip, njp, nkp, i, j, k, icount, aa
      
      ! temporary variables
      real(8) :: tmp

      ! Get basis functions and local ders
      call lintetshl(gp, NSHLu, NSD, shlu, shgradlu)

      
      ! Now calculate gradients.
      
      ! Calculate dx/dxi
      dxdxi = 0.0d0
      do icount = 1, NSHLu
            
         dxdxi(1,1) = dxdxi(1,1) + xl(icount,1)*shgradlu(icount,1)
         dxdxi(1,2) = dxdxi(1,2) + xl(icount,1)*shgradlu(icount,2)
         dxdxi(1,3) = dxdxi(1,3) + xl(icount,1)*shgradlu(icount,3)
         dxdxi(2,1) = dxdxi(2,1) + xl(icount,2)*shgradlu(icount,1)
         dxdxi(2,2) = dxdxi(2,2) + xl(icount,2)*shgradlu(icount,2)
         dxdxi(2,3) = dxdxi(2,3) + xl(icount,2)*shgradlu(icount,3)
         dxdxi(3,1) = dxdxi(3,1) + xl(icount,3)*shgradlu(icount,1)
         dxdxi(3,2) = dxdxi(3,2) + xl(icount,3)*shgradlu(icount,2)
         dxdxi(3,3) = dxdxi(3,3) + xl(icount,3)*shgradlu(icount,3)
        
      enddo      

      ! Compute inverse of deformation gradient
      dxidx = 0.0d0
      
      dxidx(1,1) = dxdxi(2,2) * dxdxi(3,3) - dxdxi(3,2) * dxdxi(2,3)
      dxidx(1,2) = dxdxi(3,2) * dxdxi(1,3) - dxdxi(1,2) * dxdxi(3,3)
      dxidx(1,3) = dxdxi(1,2) * dxdxi(2,3) - dxdxi(1,3) * dxdxi(2,2)

      tmp = 1d+0 / ( dxidx(1,1) * dxdxi(1,1)  &
                   + dxidx(1,2) * dxdxi(2,1)  &
                   + dxidx(1,3) * dxdxi(3,1) )

      dxidx(1,1) = dxidx(1,1) * tmp
      dxidx(1,2) = dxidx(1,2) * tmp
      dxidx(1,3) = dxidx(1,3) * tmp

      dxidx(2,1) = (dxdxi(2,3) * dxdxi(3,1) - dxdxi(2,1) * dxdxi(3,3)) * tmp
      dxidx(2,2) = (dxdxi(1,1) * dxdxi(3,3) - dxdxi(3,1) * dxdxi(1,3)) * tmp
      dxidx(2,3) = (dxdxi(2,1) * dxdxi(1,3) - dxdxi(1,1) * dxdxi(2,3)) * tmp
      dxidx(3,1) = (dxdxi(2,1) * dxdxi(3,2) - dxdxi(2,2) * dxdxi(3,1)) * tmp
      dxidx(3,2) = (dxdxi(3,1) * dxdxi(1,2) - dxdxi(1,1) * dxdxi(3,2)) * tmp
      dxidx(3,3) = (dxdxi(1,1) * dxdxi(2,2) - dxdxi(1,2) * dxdxi(2,1)) * tmp
      
      DetJ = 1d+0/tmp           ! Note that DetJ resides in common

      ! Global gradients
      do i = 1, NSHLu
        shgradgu(i,1) = shgradlu(i,1) * dxidx(1,1) + &
                        shgradlu(i,2) * dxidx(2,1) + &
                        shgradlu(i,3) * dxidx(3,1)
        shgradgu(i,2) = shgradlu(i,1) * dxidx(1,2) + &
                        shgradlu(i,2) * dxidx(2,2) + &
                        shgradlu(i,3) * dxidx(3,2) 
        shgradgu(i,3) = shgradlu(i,1) * dxidx(1,3) + &
                        shgradlu(i,2) * dxidx(2,3) + &
                        shgradlu(i,3) * dxidx(3,3)
      end do

      !DetJ = -1d+0*DetJ

      if (DetJ.lt.0d+0) then

	 write(*,*) "Warning - negative determinant Jacobian !"
	    write(*,*) "DetJ = ", DetJ
	    write(*,*) "DxDXI = ", dxdxi
	    write(*,*) "DxIDX = ", dxidx
	    !read(*,*)
	 
	 DetJ = -1d+0*DetJ
      endif
      
      end subroutine eval_SHAPE_3D

 
      subroutine lintetshl(gp,NSHL,NSD,shl,shgradl)

      implicit none
      
      integer :: NSHL, NSD
      real(8) :: gp(NSD), zi, eta, zeta
      real(8) :: shl(NSHL), shgradl(NSHL,NSD)
      
      zi   = gp(1)
      eta  = gp(2)
      zeta = gp(3)

      shl(1) = zi
      shl(2) = eta
      shl(3) = zeta
      shl(4) = 1.0d0-zi-eta-zeta

      shgradl = 0.0d0

      shgradl(1,1) =  1.0d0
      shgradl(2,2) =  1.0d0
      shgradl(3,3) =  1.0d0
      shgradl(4,1) = -1.0d0
      shgradl(4,2) = -1.0d0
      shgradl(4,3) = -1.0d0

      end subroutine lintetshl

      subroutine genGPandGW (gp,gw,NGAUSS)

      implicit none

      integer NGAUSS
      real(8) gp(NGAUSS,3), gw(NGAUSS)
    
      if (NGAUSS.eq.1) then

         gp(1,1) = 0.25d0  
         gp(1,2) = 0.25d0    
         gp(1,3) = 0.25d0

         gw(1) = 1d0/6d0

      elseif (NGAUSS.eq.4) then
      
         gp(1,1) = 0.5854101966249685d0  
         gp(1,2) = 0.1381966011250105d0  
         gp(1,3) = 0.1381966011250105d0

         gp(2,1) = 0.1381966011250105d0
         gp(2,2) = 0.5854101966249685d0
         gp(2,3) = 0.1381966011250105d0

         gp(3,1) = 0.1381966011250105d0
         gp(3,2) = 0.1381966011250105d0
         gp(3,3) = 0.5854101966249685d0
         
         gp(4,1) = 0.1381966011250105d0
         gp(4,2) = 0.1381966011250105d0
         gp(4,3) = 0.1381966011250105d0
  
         gw(1) = 0.2500000000000000d0/6d0
         gw(2) = 0.2500000000000000d0/6d0
         gw(3) = 0.2500000000000000d0/6d0
         gw(4) = 0.2500000000000000d0/6d0
      endif
      
      return
      end
      
