!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
module class_def

  ! type NURBSpatch

  !   integer :: P, Q, R
  !   integer :: MCP, NCP, OCP
  !   real(8), allocatable :: U_KNOT(:), V_KNOT(:), W_KNOT(:)

  ! end type NURBSpatch
  implicit none
  type bnd_class

    integer :: FACE_ID

    integer :: NFACE

    integer, allocatable :: FACE_IEN(:, :)

    integer, allocatable :: F2E(:)
    integer, allocatable :: FACE_OR(:)
    integer, allocatable :: NSHLB(:)
    integer, allocatable :: NGAUSSB(:)

    integer :: NNODE
    integer, allocatable :: BNODES(:)

    ! mapping between partitioned (local) boundary node/element
    ! and unpartitioned boundary (shell) node/element
    integer, allocatable :: L2SNODE(:), L2SELEM(:)

  end type bnd_class

  type bmesh_class
    integer :: NNODE, NSD, FACE_ID
    real(8), allocatable :: xg(:, :)
    real(8), allocatable :: dg(:, :)
    real(8), allocatable :: ug(:, :)
    integer, allocatable :: owned(:)
  end type bmesh_class

end module class_def

!------------------------------------------------------------------------
! Module for defining shell types and variables
!------------------------------------------------------------------------
module defs_shell

  implicit none

  ! Declare type mesh
  type :: mesh
    ! degree in U and V for each patch
    integer, allocatable :: P(:), Q(:)

    ! patch type
    ! 1-blade; 0-bending strips; 2-shear web; ...
    integer, allocatable :: PTYPE(:)

    ! Knot vectors in u and v directions for each element
    real(8), allocatable :: U_KNOT(:, :), V_KNOT(:, :)

    ! Size of the knot vector for each elements (e.g. NUK=P+MCP+1)
    integer, allocatable :: NUK(:), NVK(:)

    ! Control Net
    ! B_NET is reference config, B_NET_D is current config.
    ! For the pre-bend case, reference config. could change, so
    ! we created B_NET_U for undeformed, original config.
    real(8), allocatable :: B_NET(:, :), B_NET_U(:, :), B_NET_D(:, :)

    ! Boundary condition indicator for global nodes and edges
    ! respectively
    integer, allocatable :: IBC(:, :)

    ! array to store force vectors on the wind turbine blade
    real(8), allocatable :: FORCE(:, :)

    ! Global connectivity array
    integer, allocatable :: IEN(:, :), INN(:, :)

    ! number of shape functions for every element
    integer, allocatable :: NSHL(:)

    ! Bezier extraction operator
    real(8), allocatable :: Ext(:, :, :)

    ! Array for closest points (and the corresponding element)
    integer, allocatable :: CLE(:, :)
    real(8), allocatable :: CLP(:, :, :)

    ! Array for element list
    integer, allocatable :: Elm_Close(:, :, :)
    integer :: NEL_Close, NEL_Close_t

    ! Array for element gauss point location and for buiding the
    ! element list based on radial location
    real(8) :: Elm_Size
    real(8), allocatable :: Elm_Loc(:, :, :)
    integer, allocatable :: RAD_ELM_LIST(:, :, :), RAD_ELM_NUM(:, :)

    integer :: NGAUSS, NNODE, NEL, maxNSHL

    ! Store node number of the tip
    integer :: TipLoc, TipLocTr

    ! array for Solution vectors
    real(8), allocatable :: dsh(:, :), dshold(:, :), &
                            ush(:, :), ushold(:, :), &
                            ash(:, :), ashold(:, :)

    ! Surface ID and Boundary number
    integer :: FaceID, iBound
  end type mesh

  ! Declare type mesh for multi-patch
  type :: mesh_mp
    ! degree in U and V for each patch
    integer, allocatable :: P(:), Q(:)

    ! number of control points in U and V for each patch
    ! (no need for T-spline)
    integer, allocatable :: MCP(:), NCP(:)

    ! Total number of control points and elements for each patch
    integer, allocatable :: NNODE(:), NEL(:)

    ! patch type
    ! 1-blade surface; 0-bending strips; 2-shear web; ...
    integer, allocatable :: PTYPE(:)

    ! Knot vectors in u and v directions for each patch
    real(8), allocatable :: U_KNOT(:, :), V_KNOT(:, :)

    ! Control Net
    real(8), allocatable :: B_NET(:, :, :)

    ! Boundary condition indicator for global nodes and edges
    ! respectively
    integer, allocatable :: IBC(:, :, :)

    ! array to store force vectors on the wind turbine blade
    real(8), allocatable :: FORCE(:, :, :)

    ! Mapping between patches and global reduced numbering
    ! e.g., MAP(1,2) = 3 means 1st patch, 2nd node points to
    !   global reduced node number 3
    integer, allocatable :: MAP(:, :)
  end type mesh_mp

  ! Declare type shell (for wind turbine blade)
  type :: shell_bld

    type(mesh_mp) :: mNRB
    type(mesh)    :: NRB

    type(mesh)    :: TSP, BEZ

    type(mesh)    :: FEM

    ! number of patches for Blade Surface (S). Matches the fluid mesh
    ! number of patches for blade structure (B). May include shear webs
    ! number of total patches including bending strips (T)
    integer :: NPS, NPB, NPT

    integer :: M_Flag, T_Flag
    real(8) :: RHSGtol, G_fact(3), RHSGNorm

    ! row, col, and total of nonzero entries for sparse structure
    integer, allocatable :: row(:), col(:)
    integer :: icnt

    ! The right hand side load vector G and left hand stiffness matrix K
    real(8), allocatable :: RHSG(:, :), LHSK(:, :), &
                            RHSG_EXT(:, :), RHSG_GRA(:, :)

!    ! Solution vectors
!    real(8), allocatable :: yg(:,:), dg(:,:), tg(:), &
!                            mg(:), dl(:,:)

    ! material matrix for composite
    integer :: NMat
!    real(8), allocatable :: matA(:,:,:), matB(:,:,:), matD(:,:,:)
    real(8), allocatable :: matA(:, :, :, :), matB(:, :, :, :), matD(:, :, :, :), Density(:, :), Thickness(:, :)

    ! number of newton iterations for shell
    integer, allocatable :: Nnewt(:)

    ! Torque computed on shell mesh
    real(8) :: Tq1, Tq2

    ! Blade rotation. 0 degree is the straight-up position
    real(8) :: BldRot

    integer :: bmap
  end type shell_bld

  ! Declare type shell (for non-matching boundaries)
  type :: shell_nmb

    type(mesh), allocatable :: FEM(:)

  end type shell_nmb
end module defs_shell
!----------------------------------------------------------------------
! Module for common parameters (e.g. pi)
!----------------------------------------------------------------------
module commonpars
  implicit none
  save

  real(8), parameter :: pi = 3.14159265358979323846264338328d0
  integer, parameter :: NRES = 2
  integer, parameter :: NSubdomain = 2
 ! Assemble field
  integer, parameter :: ASSEMBLE_FIELD_NONE = 0
  integer, parameter :: ASSEMBLE_FIELD_NS = 1
  integer, parameter :: ASSEMBLE_FIELD_LS = 2
  integer, parameter :: ASSEMBLE_FIELD_VOF = ASSEMBLE_FIELD_LS
  integer, parameter :: ASSEMBLE_FIELD_TEM = 4

  ! Assemble tensor
  integer, parameter :: ASSEMBLE_TENSOR_NONE = 0
  integer, parameter :: ASSEMBLE_TENSOR_SCALAR = 1
  integer, parameter :: ASSEMBLE_TENSOR_VEC = 2
  integer, parameter :: ASSEMBLE_TENSOR_MAT = 4
end module commonpars


!------------------------------------------------------------------------
!     Module for storing arrays and allocation routines
!------------------------------------------------------------------------
module aAdjKeep

  use class_def
  use commonpars
  ! use defs_shell

  implicit none
  save

  ! Mesh
  real(8), allocatable :: xg(:, :), wg(:)

  integer, allocatable :: IEN(:, :), NodeID(:)
  integer, allocatable :: EPID(:), EIJK(:, :)
  integer, allocatable :: ELM_ID(:)
  integer :: ELM_ID_set(NRES), bmesh_ID_set(NRES), SRC_ID_set(NRES)

  
  real(8) :: omega_x

  type(bnd_class), allocatable :: bound(:)

  type(bmesh_class) :: bmesh(2)

  !type(NURBSpatch), allocatable :: patch(:)

  ! Contraint flags
  integer, allocatable :: IPER(:)
  integer, allocatable :: IBC(:, :)
  ! logical :: IS_SOLID_NODE_ASSIGNED
  ! integer, allocatable :: IS_SOLID_NODE(:)

  ! Type flags
  integer, allocatable :: EL_TYP(:), D_FLAG(:), P_FLAG(:)

  ! Spars Struc
  integer, allocatable :: row(:), col(:)

  integer :: InNNODE1, InCELL
  integer, allocatable :: ien_inlet(:, :)
  real(8), allocatable :: xg_inlet(:, :)

  ! The right hand side load vector G and left hand stiffness matrix K
  real(8), allocatable :: RHSGU(:, :), RHSGM(:, :), RHSGP(:), &
                          RHSGLS(:), RHSGTEM(:)

  real(8), allocatable :: lhsK11(:, :), lhsK12(:, :), lhsK22(:, :), &
                          lhsG(:, :), lhsD1(:, :), lhsD2(:, :), &
                          lhsM(:), lhsLS(:), LHSPi(:, :), &
                          LHSPti(:, :), LHSKi(:, :), lhsMi(:), &
                          LHSLSi(:), LHSmass(:), &
                          LHSlsu(:, :), LHSPls(:), LHSlsP(:), LHSUls(:, :), &
                          LHSTem(:)

  real(8), allocatable :: lhsgq(:), rhsgq(:)
  ! Solution vectors
  real(8), allocatable :: dg(:, :), dgold(:, :), &
                          ug(:, :), ugold(:, :), &
                          ugm(:, :), ugmold(:, :), &
                          acg(:, :), acgold(:, :), &
                          acgm(:, :), acgmold(:, :), &
                          pg(:), pgold(:), &
                          phig(:), phigold(:), &
                          rphig(:), rphigold(:), &
                          rTg(:), rTgold(:), &
                          Tg(:), Tgold(:)
  real(8), allocatable :: uavg(:, :), pavg(:)

  ! Rigid body
  ! real(8) :: vbn0(3), vbn1(3)
  ! real(8) :: dbn0(3), dbn1(3)
  ! real(8) :: wbn0(3), wbn1(3)
  ! real(8) :: Rn0(3, 3), Rn1(3, 3)

  ! First P-K Stress
  ! real(8), allocatable :: FPKS(:, :, :)

  ! Array for Prism
  integer, allocatable :: ELMNSHL(:), ELMNGAUSS(:)

  ! global information for individual blades
  ! type(mesh) :: blade(3)

  ! real(8) :: Center_Rot(3)

end module aAdjKeep

!----------------------------------------------------------------------
! Module for common variables
!----------------------------------------------------------------------
module commonvars

  implicit none
  save

  ! FLAGS
  integer :: fem_flag
  logical :: grad_prop_in_diff

  integer::check_num, InNNODE
  integer::index_plane_p, index_plane_n, max_planes
  real(8)::Delt_x, x_inflow, remain_x, position_init
  integer::tmp_num1, tmp_num2
  real(8), dimension(1:3)::C_tmp1, C_tmp2
  real(8)::V_tmp1, V_tmp2
  real(8), allocatable :: inflow_velocity(:, :), inflow_pressure(:), inflow_rho(:)
  real(8), allocatable :: inflow_x(:), inflow_y(:), inflow_z(:)
  integer, allocatable :: inflow_bnodes(:)

  real(8), allocatable :: inflow_velocity_check(:, :)
  character(len=1000) :: filename_p, filename_n, filename_check1, filename_check2
  real(8), dimension(1:3)::Inflow_velocity_p, Inflow_velocity_n
  logical::turbulence_flag
  real(8)::interpolation_factor, inflow_pressure_p, inflow_pressure_n, inflow_rho_p, inflow_rho_n
  real(8)::top_flag
  ! Mesh
  integer :: NSD, NNODE, NELEM, NBOUND, NPATCH, NSHLBmax, icnt, &
             NBlade, maxNSHL
  logical :: iga
  integer :: use_vms
  real(8) :: DetJ, DetJb, DetJinv, hglob

  ! Time step
  real(8) :: Delt, Dtgl, rhoinf, beti, gami, alfi, almi, &
             mgam, ogam, time, &
             conv_time, mono_time, move_time, shel_time
  ! real(8) :: lambda
  integer :: Nstep, ifq, ifq_sh, ifq_tq, mono_iter

  ! Navier-Stokes solver
  real(8) :: mua, rhoa, muw, rhow
  real(8) :: cpa, cpw, kappaa, kappaw
  real(8) :: rhos, mus, cps, kappas

  real(8) :: htca, htcw
  real(8) :: Trm = 273d0+25d0
  real(8) :: Tref = 5d2 - 75d0
  real(8) :: NS_kdc_w, NS_kdc_a, fine_tau

  real(8) :: NS_GMRES_atol, NS_NL_Uatol, NS_NL_Patol
  real(8) :: NS_GMRES_tol, NS_NL_Utol, NS_NL_Ptol
  integer :: NS_GMRES_itermin, NS_GMRES_itermax, NS_NL_itermax, &
             NS_hess_flag

  ! LevelSet Convection solver
  real(8) :: LSC_kdc
  real(8) :: LSC_GMRES_atol, LSC_NL_atol
  real(8) :: LSC_GMRES_tol, LSC_NL_tol
  integer :: LSC_GMRES_itermin, LSC_GMRES_itermax, LSC_NL_itermax, &
             LSC_pred_step

  ! Temperature Advection-diffusion solver
  real(8) :: TEM_kdc
  real(8) :: TEM_GMRES_atol, TEM_NL_atol
  real(8) :: TEM_GMRES_tol, TEM_NL_tol
  integer :: TEM_GMRES_itermin, TEM_GMRES_itermax, TEM_NL_itermax

  ! Rigid body solver
  real(8) :: RB_NL_Ftol, RB_NL_Mtol


  ! Mesh solver
  real(8) :: Mesh_GMRES_tol, Mesh_NL_tol
  integer :: Mesh_GMRES_itermin, Mesh_GMRES_itermax
  ! Level Set Redistance solver
  real(8) :: LSRD_kdc, LSRD_penalty_fac
  real(8) :: LSRD_pseudoDTGL_fac
  real(8) :: LSRD_GMRES_tol, LSRD_NL_tol
  integer :: LSRD_GMRES_itermin, LSRD_GMRES_itermax, LSRD_NL_itermax

  ! Level Set massfix
  real(8) :: Mass_init
  real(8) :: Mass_NL_tol
  integer :: Mass_NL_itermax

  ! Rigid body
  real(8) :: massb, Ibhat(3, 3), xcg(3), Fb(3), Mb(3)
  integer :: VBC(3), MBC(3)

  ! Domain and Hull
  real(8) :: domain_top, domain_bottom, &
             domain_left, domain_right, &
             domain_front, domain_back, &
             water_level, water_depth, air_height, hull_length

  ! Wave generating wall
  real(8) :: wave_periode, wave_length, wave_amp, wave_angle

  ! Gravity
  real(8) :: gravity, gravvec(3)

  ! Interface
  real(8) :: mp_eps, beta_t, phi_inf, &
             phi_t, C_t, Re, Pe, Gr, Fr, &
             Ra, Sr, Pr, dphi_bg(3), cross_flag

  !real(8) :: mu, rho, dmudphi, drhodphi
  real(8) :: kappa
  !real(8) :: cp, hk
  real(8) :: Ts, c_cond, c_evap, lh
  real(8) :: C_DMDOT
  ! Setup
  real(8) :: Froude, Uin
  integer :: BCtype(99)

  integer, allocatable :: BCugType(:, :)
  real(8), allocatable :: BCugValu(:, :), phi_bg(:)
  integer, allocatable :: BCphigType(:), BCTgType(:)
  real(8), allocatable :: BCphigValu(:), BCTgValu(:)

  ! real(8) :: usettle
  ! Flags
  ! logical :: move, mono, conv, shel

  ! rotation
  real(8) :: theta, thetd, thedd, thetaOld, thetdOld, theddOld, &
             maxthetd, torque1, torque2, torque3, torque4, torque_SH, &
             moment1, moment2, force_trac(3), force_cons(3)

  logical :: solshell, nonmatch
end module commonvars

!----------------------------------------------------------------------
!     Module for mpi
!----------------------------------------------------------------------
module mpi
  implicit none
  save

  include "mpif.h"

  integer, parameter :: mpi_master = 0
  integer, parameter :: maxtask = 50
  integer, parameter :: maxseg = 15000

  integer :: numnodes, myid, mpi_err
  integer :: status(MPI_STATUS_SIZE)
  integer :: lstrideu, lstridep, numtask
  integer :: lfrontu, maxfrontu
  integer :: lfrontp, maxfrontp
  integer :: nlworku, nlworkp
  integer :: itag, iacc, iother, numseg, isgbeg, itkbeg, isgend
  integer :: sevsegtypeu(maxtask, 15)
  integer :: sevsegtypep(maxtask, 15)
  logical :: ismaster

  integer, allocatable :: ilworku(:), ilworkp(:)
end module mpi

