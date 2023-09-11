module class

!------------------------------------------------------------------------
!                 
!------------------------------------------------------------------------  
type NURBSpatch
  
  integer :: P, Q, R
  integer :: MCP, NCP, OCP
  real(8), allocatable :: U_KNOT(:), V_KNOT(:), W_KNOT(:)
  
end type NURBSpatch

!------------------------------------------------------------------------
!                 
!------------------------------------------------------------------------
type bnd_class

  integer :: FACE_ID

  integer :: NFACE  
  
  integer, allocatable :: FACE_IEN(:,:)
    
  integer, allocatable :: F2E(:)
  integer, allocatable :: FACE_OR(:)
      
  integer :: NNODE
  integer, allocatable :: BNODES(:)
      
end type bnd_class

!------------------------------------------------------------------------
!                 
!------------------------------------------------------------------------
type mesh_class

  integer :: NSD, NSHLmax, NSHLB  
  integer :: NNODE, NELEM, NBOUND, NPATCH

  real(8), allocatable :: xg(:,:)
  integer, allocatable :: IEN(:,:), EID(:), NSHL(:), NID(:)
  
  type(bnd_class),  allocatable :: bound(:) 
  type(NURBSpatch), allocatable :: patch(:) 
  
  real(8), allocatable :: wg(:)
  integer, allocatable :: EPID(:)
  integer, allocatable :: EIJK(:,:)    
  
end type mesh_class

!------------------------------------------------------------------------
!                 
!------------------------------------------------------------------------
type map_class
  
  integer, allocatable :: node(:,:)
  
end type map_class

type patch_con_class

  real(8)  dsize,eps,treshold
  
  integer  npc
  integer, allocatable :: patch(:,:)
  
  integer  nbc
  integer, allocatable :: bnd(:,:)
  
end type patch_con_class


end module
