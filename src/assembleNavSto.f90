!======================================================================
!
!======================================================================
subroutine assembleNavSto(assemble_tensor_flag, subdomain_flag)
  use aAdjKeep
  use mpi
  use commonvars
  use commonpars

  implicit none

  integer, intent(in) :: assemble_tensor_flag ! assembe Jacobian mat or vec
  integer, intent(in) :: subdomain_flag ! assemble NS on specified subdomains

  real(8) :: dgAlpha(NNODE, NSD), &
             ugAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
             acgAlpha(NNODE, NSD), acgmAlpha(NNODE, NSD), &
             pgAlpha(NNODE)

  real(8) :: t1, t2

  !---------------------------------
  ! Alpha stage
  !---------------------------------
  acgAlpha = acgold + almi*(acg - acgold)
  acgmAlpha = acgmold + almi*(acgm - acgmold)
  ugAlpha = ugold + alfi*(ug - ugold)
  ugmAlpha = ugmold + alfi*(ugm - ugmold)
  dgAlpha = dgold + alfi*(dg - dgold)
  pgAlpha = pg

  !---------------------------------
  ! zero out RHS and LHS
  !---------------------------------

  if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
    RHSGu = 0.0d0
    RHSGp = 0.0d0
  end if

  if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
    LHSK11 = 0.0d0
    LHSG = 0.0d0
    LHSD1 = 0.0d0
    LHSM = 0.0d0
  endif

  if (myid .eq. 0) then
    call CPU_TIME(t1)
  endif
  call IntElmAss_NS_overset( &
    dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
    acgmAlpha, pgAlpha, assemble_tensor_flag, subdomain_flag)

  call FaceAssembly_NS_weak(&
    dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
    acgmAlpha, pgAlpha, assemble_tensor_flag, subdomain_flag)

  call FaceAssembly_NS_outflow(&
    dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
    acgmAlpha, pgAlpha, assemble_tensor_flag, subdomain_flag)
  if (myid .eq. 0) then
    call CPU_TIME(t2)
    write (*, *) "Total time assemble:", t2 - t1, "seconds"
  end if

  if (numnodes > 1 .and. iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      call commu(RHSGp, 1, 'in ')
      call commu(RHSGu, NSD, 'in ')
  end if

end subroutine assembleNavSto
!======================================================================
!
!======================================================================
subroutine IntElmAss_NS_overset(&
  dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
  acgmAlpha, pgAlpha, assemble_tensor_flag, subdomain_flag)
  use aAdjKeep
  use commonvars
  use commonpars
  use mpi
  implicit none

  real(8), intent(in) :: dgAlpha(NNODE, NSD), ugAlpha(NNODE, NSD), &
                         acgAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
                         acgmAlpha(NNODE, NSD), pgAlpha(NNODE)
  integer, intent(in) :: subdomain_flag ! assemble NS on specified subdomains
  integer, intent(in) :: assemble_tensor_flag ! assembe Jacobian mat or vec

  ! Local variables
  integer :: iel, igauss, aa, bb, i, j, hess_flag, idx, nshl, NGAUSS
  real(8) :: volm, vol_ex, uprime(NSD)
  real(8) :: tauM, tauC, tauBar, tauBar1, pri, k_dc, k_dc_phi, tauP, tauLS

  real(8), allocatable :: xKebe11(:, :, :), xGebe(:, :, :)
  real(8), allocatable :: xDebe1(:, :, :), xMebe(:, :)
  real(8), allocatable :: xLSebe(:, :), xLSUebe(:, :, :)
  real(8), allocatable :: xULSebe(:, :, :), xPLSebe(:, :)

  real(8), allocatable :: Rhsu(:, :), Rhsp(:), Rhsphi(:)

  real(8), allocatable :: shlu(:), shgradgu(:, :), shhessgu(:, :, :)

  real(8), allocatable :: xl(:, :)
  real(8), allocatable :: dl(:, :), wl(:)
  real(8), allocatable :: acml(:, :), uml(:, :)
  real(8), allocatable :: fl(:, :), ul(:, :), acl(:, :), pl(:)

  real(8), allocatable :: gp(:, :), gw(:)
  real(8), allocatable :: shconv(:), shconv_full(:)

  ! integer, allocatable :: ibc_loc(:, :)

  real(8) :: di(NSD), ui(NSD), aci(NSD), umi(NSD), acmi(NSD), &
             fi(NSD), xi(NSD), ddidxi(NSD, NSD), &
             dphidxi(NSD), duidxi(NSD, NSD), dphidxidxj(NSD, NSD), &
             duidxixj(NSD, NSD, NSD), duiolddxi(NSD, NSD), &
             dxidx(NSD, NSD), dpridxi(NSD), rLi(NSD)
  real(8) :: uadvi(NSD), uadvi_full(NSD)
  real(8) :: Ti, rTi, dTdxi(NSD), Tic
  real(8) :: ns_kdc

  real(8) :: Gij(NSD, NSD), Ginv(NSD, NSD)
  real(8) :: cfl, gcfl, cfl_loc
  real(8) :: Ftens(NSD, NSD), Stens(NSD, NSD), Ctens(NSD, NSD, NSD, NSD)
  real(8) :: res_phic_tmp1, res_phic_tmp2
  real(8) :: Qi, Si(3, 3), Omegai(3, 3)
  logical :: is_target

  real(8) :: rhoi, mui, cpi, hki
  real(8) :: divu, lambda, rm(NSD), tmp1(NSD)
  real(8) :: fact1, fact2

  real(8) :: DetJgw

  cfl = 0.0d0

  fact1 = almi
  fact2 = alfi * delt * gami


  NGAUSS = -1
  NSHL = -1

  ! loop over elements
  do iel = 1, NELEM
!  write(*,*) "iel:",iel, NELEM
    if (NSHL /= ELMNSHL(iel)) then

      if (NSHL >= 0) then
        deallocate (xKebe11, xGebe)                      ! 2
        deallocate (xDebe1, xMebe)                       ! 4
        deallocate (xLSebe, xLSUebe)                     ! 6
        deallocate (xULSebe, xPLSebe)                    ! 8
        deallocate (Rhsu, Rhsp, Rhsphi)                  ! 11
        deallocate (shlu, shgradgu, shhessgu)            ! 14
        deallocate (xl)                                  ! 15
        deallocate (dl, wl)                              ! 17
        deallocate (acml, uml)                           ! 19
        deallocate (fl, ul, acl, pl)                     ! 23
        deallocate (gp, gw)                              ! 29
        deallocate (shconv, shconv_full)                 ! 31
      end if

      NSHL = ELMNSHL(iel)
      NGAUSS = ELMNGAUSS(iel)
      allocate (xKebe11(NSD * NSD, NSHL, NSHL), xGebe(NSD, NSHL, NSHL))
      allocate (xDebe1(NSD, NSHL, NSHL), xMebe(NSHL, NSHL))
      allocate (xLSebe(NSHL, NSHL), xLSUebe(NSD, NSHL, NSHL))
      allocate (xULSebe(NSD, NSHL, NSHL), xPLSebe(NSHL, NSHL))
      allocate (Rhsu(NSD, NSHL), Rhsp(NSHL), Rhsphi(NSHL))
      allocate (shlu(NSHL), shgradgu(NSHL, NSD), shhessgu(NSHL, NSD, NSD))
      allocate (xl(NSHL, NSD))
      allocate (dl(NSHL, NSD), wl(NSHL))
      allocate (acml(NSHL, NSD), uml(NSHL, NSD))
      allocate (fl(NSHL, NSD), ul(NSHL, NSD), acl(NSHL, NSD), pl(NSHL))
      allocate (gp(NGAUSS, NSD), gw(NGAUSS))
      allocate (shconv(NSHL), shconv_full(NSHL))
      ! get Gaussian points and weights
      call genGPandGW(gp, gw, NGAUSS)
    end if


    ! find out if this is a target element
    is_target = .false.
    if(iand(subdomain_flag, 1) > 0) then
      is_target = is_target .or. ELM_ID(iel) == ELM_ID_SET(1)
    endif
    if(iand(subdomain_flag, 2) > 0) then
      is_target = is_target .or. ELM_ID(iel) == ELM_ID_SET(2)
    endif


    fl = 0.0d0

    ! Get local solution arrays
    do i = 1, NSHL
      idx = IEN(iel, i)
      xl(i, :) = xg(idx, :)
      dl(i, :) = dgAlpha(idx, :)
      wl(i) = wg(idx)
      ul(i, :) = ugAlpha(idx, :)
      acl(i, :) = acgAlpha(idx, :)
      uml(i, :) = ugmAlpha(idx, :)
      acml(i, :) = acgmAlpha(idx, :)
      pl(i) = pgAlpha(idx)
    end do


    ! initialize local stiffness matrix
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      xKebe11 = 0.0d0
      xGebe = 0.0d0
      xDebe1 = 0.0d0
      xMebe = 0.0d0

    end if
    ! initialize local load vector
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      Rhsu = 0.0d0
      Rhsp = 0.0d0
    end if
    ! Loop over integration points (NGAUSS in each direction)
    do igauss = 1, NGAUSS

      ! Get Element Shape functions and their gradients
      ! initialize
      shlu(:) = 0.0d0
      shgradgu(:, :) = 0.0d0
      shhessgu(:, :, :) = 0.0d0
      ! hess_flag = .false.
      call eval_shape(nshl, iel, gp(igauss, :), xl, dl, wl, shlu, &
                      shgradgu, shhessgu, dxidx, Gij, Ginv, 0)
      DetJgw = DetJ * gw(igauss)

      call e3int_qr(NSHL, NSD, NSD, shlu, xl, xi(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, dl, di(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, ul, ui(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, acl, aci(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, uml, umi(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, acml, acmi(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, fl, fi(1))
      call e3int_qr(NSHL, NSD, 1, shlu, pl, pri)

      call e3int_qr_grad(NSHL, NSD, NSD, shgradgu, dl, ddidxi)
      call e3int_qr_grad(NSHL, NSD, NSD, shgradgu, ul, duidxi)
      call e3int_qr_grad(NSHL, NSD, 1, shgradgu, pl, dpridxi)

      call e3int_qr_hess(NSHL, NSD, NSD, shhessgu, ul, duidxixj)

      rhoi = rhoa
      mui = mua
      ns_kdc = NS_kdc_a

      ! ALE Advective Velocity
      uadvi(:) = ui(:) - umi(:)
      divu = duidxi(1, 1) + duidxi(2, 2) + duidxi(3, 3)

      fi(:) = gravvec(:) * rhoi
      rm(:) = rhoi * (aci(:) + matmul(duidxi, uadvi(:)) ) - fi(:)

      ! calculate residual for NS
      rLi(:) = rm(:) + dpridxi(:) - mui * (duidxixj(:, 1, 1) + duidxixj(:, 2, 2) + duidxixj(:, 3, 3))
      ! calculate the STABLIZATION
      call e3STAB_3D(NSD, Gij, Delt, uadvi, rhoi, mui, tauM, tauP, tauC)

      tauBar = 0.0d0

      call e3DC_beta2(NSD, ns_kdc, Gij, rLi, k_dc)

      if(.true.) then
        call e3CFL(NSD, uadvi, Gij, Delt, cfl_loc)
        cfl = max(cfl, cfl_loc)
      endif
      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0 .and. is_target) then
      call e3LHS_3D_fluid_Old( &
        nshl, rhoi, mui, ui, umi, aci, pri, duidxi, &
        dpridxi, rLi, tauM, tauP, tauC, tauBar, ns_kdc, gw(igauss), &
        shlu, shgradgu, shhessgu, xKebe11, xGebe, xDebe1, xMebe)
      end if

      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0 .and. is_target) then
        call e3RHS_3D_fluid( &
          NSHL, rhoi, mui, ui, aci, umi, acmi, uadvi, pri, &
          rLi, fi, duidxi, ddidxi, tauM, tauP, &
          tauC, tauBar, ns_kdc, gw(igauss), shlu, &
          shgradgu, shhessgu, Rhsu, Rhsp)
      end if

    end do

    ! Apply Dirichlet BCs
    call BCLhs_3D(nshl, iel, xKebe11, xGebe, xDebe1, &
                  xMebe, Rhsu, Rhsp)

    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      call FillSparseMat_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe)
    end if
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      do bb = 1, NSHL
        RHSGu(IEN(iel, bb), :) = RHSGu(IEN(iel, bb), :) + Rhsu(:, bb)
        RHSGp(IEN(iel, bb)) = RHSGp(IEN(iel, bb)) + Rhsp(bb)
      enddo
    end if

  end do

  deallocate (xKebe11, xGebe)                      ! 2
  deallocate (xDebe1, xMebe)                       ! 4
  deallocate (xLSebe, xLSUebe)                     ! 6
  deallocate (xULSebe, xPLSebe)                    ! 8
  deallocate (Rhsu, Rhsp, Rhsphi)                  ! 11
  deallocate (shlu, shgradgu, shhessgu)            ! 14
  deallocate (xl)                                  ! 15
  deallocate (dl, wl)                              ! 17
  deallocate (acml, uml)                           ! 19
  deallocate (fl, ul, acl, pl)                     ! 23
  deallocate (gp, gw)                              ! 29
  deallocate (shconv, shconv_full)                 ! 31
  if (.true.) then
    ! Find largest CFL-number and output to screen
    if (numnodes > 1) then
      call MPI_ALLReduce(cfl, gcfl, 1, MPI_DOUBLE_PRECISION, &
                         MPI_MAX, MPI_COMM_WORLD, mpi_err)
    else
      gcfl = cfl
    end if
    if (ismaster) then
      write (*, '(40("-"))')
      write (*, *) "    CFL = ", gcfl
    end if
  end if

 
end subroutine IntElmAss_NS_overset
!======================================================================
!
!======================================================================

