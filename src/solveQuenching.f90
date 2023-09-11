!======================================================================
!
!======================================================================
subroutine assembleQuenching(assemble_tensor_flag, assemble_field_flag)
  use aAdjKeep
  use mpi
  use commonvars
  use commonpars

  implicit none

  integer, intent(in) :: assemble_tensor_flag ! assembe Jacobian mat or vec
  integer, intent(in) :: assemble_field_flag ! assemble NS + LS/VOF + Tem

  real(8) :: dgAlpha(NNODE, NSD), &
             ugAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
             acgAlpha(NNODE, NSD), acgmAlpha(NNODE, NSD), &
             pgAlpha(NNODE), phigAlpha(NNODE), rphigAlpha(NNODE), &
             rTgAlpha(NNODE), TgAlpha(NNODE)

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

  phigAlpha = phigold + alfi*(phig - phigold)
  rphigAlpha = rphigold + almi*(rphig - rphigold)

  TgAlpha = Tgold + alfi*(Tg - Tgold)
  rTgAlpha = rTgold + almi*(rTg - rTgold)

  !---------------------------------
  ! zero out RHS and LHS
  !---------------------------------

  if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
    RHSGu = 0.0d0
    RHSGp = 0.0d0
    RHSGls = 0.0d0
    RHSGtem = 0.0d0
  end if

  if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
    LHSK11 = 0.0d0
    LHSG = 0.0d0
    LHSD1 = 0.0d0
    LHSM = 0.0d0

    LHSLS = 0.0d0
    LHSULS = 0.0d0
    LHSLSU = 0.0d0
    LHSPLS = 0.0d0

    LHSTem = 0.0d0
  endif

  if (myid .eq. 0) then
    call CPU_TIME(t1)
  endif
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF) > 0) then

    call IntElmAss_NSVOF_Quenching_STAB( &
      dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
      acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
      TgAlpha, rTgAlpha, &
      assemble_tensor_flag)

    call FaceAssembly_NS_weak_CF(&
      dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
      acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
      assemble_tensor_flag)

    call FaceAssembly_NS_outflow(&
      dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
      acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
      assemble_tensor_flag)

  end if
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_TEM) > 0) then
    call IntElmAss_Tem_Quenching_STAB( &
      dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
      acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
      TgAlpha, rTgAlpha, &
      assemble_tensor_flag)
    call FaceAssembly_ConvectiveH(&
         dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
         acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
         TgAlpha, assemble_tensor_flag)

  end if

  if (myid .eq. 0) then
    call CPU_TIME(t2)
    write (*, *) "Total time assemble:", t2 - t1, "seconds"
  end if

  if (numnodes > 1 .and. iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
    if(iand(assemble_field_flag, ASSEMBLE_FIELD_NS) > 0) then
      call commu(RHSGp, 1, 'in ')
      call commu(RHSGu, NSD, 'in ')
    endif
    if(iand(assemble_field_flag, ASSEMBLE_FIELD_VOF) > 0) then
      call commu(RHSGls, 1, 'in ')
    endif
    if(iand(assemble_field_flag, ASSEMBLE_FIELD_TEM) > 0) then
      call commu(RHSGTem, 1, 'in ')
    endif
  end if

end subroutine assembleQuenching
!======================================================================
!
!======================================================================

subroutine IntElmAss_NSVOF_Quenching_STAB(&
  dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
  acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
  TgAlpha, rTgAlpha, &
  assemble_tensor_flag)
  use aAdjKeep
  use commonvars
  use commonpars
  use mpi
  implicit none


  integer, intent(in) :: assemble_tensor_flag
  real(8), intent(in) :: dgAlpha(NNODE, NSD), ugAlpha(NNODE, NSD), &
                         acgAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
                         acgmAlpha(NNODE, NSD), pgAlpha(NNODE), &
                         phigAlpha(NNODE), rphigAlpha(NNODE), &
                         TgAlpha(NNODE), rTgAlpha(NNODE)
                         ! Local variables
  integer :: iel, igauss, aa, bb, i, j, hess_flag, idx, nshl, NGAUSS

  real(8) :: volm, vol_ex, uprime(NSD)
  real(8) :: tauM, tauC, tauBar, tauBar1, pri, k_dc, k_dc_phi, tauP, tauLS

  !real(8), allocatable :: xKebe11(:, :, :), xGebe(:, :, :), xDebe1(:, :, :), &
  !                        xMebe(:, :), Rhsu(:, :), Rhsm(:, :), Rhsp(:), Rhsq(:), Rhsl(:), &
  !                        xLSebe(:, :), xLSUebe(:, :, :), xPLSebe(:, :), &
  !                        xULSebe(:, :, :), Rhsphi(:)

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
  real(8), allocatable :: phil(:), rphil(:)
  real(8), allocatable :: Tl(:), rTl(:)

  real(8), allocatable :: gp(:, :), gw(:)
  real(8), allocatable :: shconv(:), shconv_full(:)

  ! integer, allocatable :: ibc_loc(:, :)

  real(8) :: di(NSD), ui(NSD), aci(NSD), umi(NSD), acmi(NSD), rphii, &
             fi(NSD), xi(NSD), ddidxi(NSD, NSD), phii, &
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
  logical :: is_fluid

  real(8) :: rhoi, mui, cpi, hki
  real(8) :: mdot, vdot, phic, dmdphii, drhodphi, dmudphi
  real(8) :: divu, lambda, rm(NSD), tmp1(NSD)
  real(8) :: fact1, fact2

  real(8) :: DetJgw

  volm = 0.0d0
  vol_ex = 0.0d0
  cfl = 0.0d0

  fact1 = almi
  fact2 = alfi * delt * gami


  NGAUSS = -1
  NSHL = -1
!  rhsgq = 0d0
!  lhsgq = 0d0

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
        deallocate (phil, rphil)                         ! 25
        deallocate (Tl, rTl)                             ! 27
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
      allocate (phil(NSHL), rphil(NSHL))
      allocate (Tl(NSHL), rTl(NSHL))
      allocate (gp(NGAUSS, NSD), gw(NGAUSS))
      allocate (shconv(NSHL), shconv_full(NSHL))
      ! get Gaussian points and weights
      call genGPandGW(gp, gw, NGAUSS)
    end if

    is_fluid = ELM_ID(iel) == 102

    fl = 0.0d0
!    do i = 1, NSHL
!      fl(i,3) = -1.0d0/(Fr**2.0d0)
!    end do

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
      phil(i) = phigAlpha(idx)
      rphil(i) = rphigAlpha(idx)
      Tl(i) = TgAlpha(idx)
      rTl(i) = rTgAlpha(idx)
    end do


    ! initialize local stiffness matrix
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      xKebe11 = 0.0d0
      xGebe = 0.0d0
      xDebe1 = 0.0d0
      xMebe = 0.0d0

      xLSebe = 0.0d0
      xLSUebe = 0.0d0
      xULSebe = 0.0d0
      xPLSebe = 0.0d0

    end if
    ! initialize local load vector
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      Rhsu = 0.0d0
      Rhsp = 0.0d0
      Rhsphi = 0.0d0
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

      !call e3int_fluid(nshl, xl, dl, ul, acl, uml, acml, &
      !                 pl, fl, phil, Tl, rTl, shlu, shgradgu, &
      !                 shgradgu, shhessgu, dxidx, Ginv, &
      !                 di, ui, aci, umi, acmi, pri, fi, &
      !                 ddidxi, duidxi, duidxixj, dpridxi, &
      !                 phi, dphidxi, dphidxidxj, rphil, dphidti, xi, &
      !                 rTi, Ti, dTdxi, rhoi, mui, cpi, hki, &
      !                 rLi, res_phic_tmp1)
      call e3int_qr(NSHL, NSD, NSD, shlu, xl, xi(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, dl, di(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, ul, ui(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, acl, aci(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, uml, umi(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, acml, acmi(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, fl, fi(1))
      call e3int_qr(NSHL, NSD, 1, shlu, pl, pri)
      call e3int_qr(NSHL, NSD, 1, shlu, phil, phii)
      call e3int_qr(NSHL, NSD, 1, shlu, rphil, rphii)
      call e3int_qr(NSHL, NSD, 1, shlu, Tl, Ti)
      call e3int_qr(NSHL, NSD, 1, shlu, rTl, rTi)

      call e3int_qr_grad(NSHL, NSD, NSD, shgradgu, dl, ddidxi)
      call e3int_qr_grad(NSHL, NSD, NSD, shgradgu, ul, duidxi)
      call e3int_qr_grad(NSHL, NSD, 1, shgradgu, pl, dpridxi)
      call e3int_qr_grad(NSHL, NSD, 1, shgradgu, phil, dphidxi)
      call e3int_qr_grad(NSHL, NSD, 1, shgradgu, Tl, dTdxi)

      call e3int_qr_hess(NSHL, NSD, NSD, shhessgu, ul, duidxixj)
      call e3int_qr_hess(NSHL, NSD, 1, shhessgu, phil, dphidxidxj)

        ! call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(ismaster) write(*,*) "Evaluate on qr", igauss, NGAUSS
      phic = max(min(phii, 1.0d0), 0.0d0)
      ! phic = phii
      ! Tic = max(min(Ti, 773d0), 348d0)
      Tic = Ti
      if(Ti > Ts) then
        mdot = c_evap * (1-phic) * rhow * (Tic - Ts) / Ts
        dmdphii = -c_evap * rhow * (Tic - Ts) / Ts * C_DMDOT
      else
        mdot = c_cond * (phic) * rhoa * (Tic - Ts) / Ts
        dmdphii = c_cond * rhoa * (Tic - Ts) / Ts * C_DMDOT
      endif
      vdot = mdot / rhoa - mdot / rhow

      if(is_fluid) then
        call prop_interp(rhow, rhoa, phic, rhoi)
        call prop_interp(muw, mua, phic, mui)
      else ! solid
        rhoi = rhos
        mui = mus
      endif
      lambda = -2d0 / 3d0 * mui
      drhodphi = rhoa - rhow
      dmudphi = mua - muw
      if(phii < 0d0 .or. phii > 1d0) then
        drhodphi = 0d0
        dmudphi = 0d0
        dmdphii = 0d0
      endif

      ! ALE Advective Velocity
      uadvi(:) = ui(:) - umi(:)
      divu = duidxi(1, 1) + duidxi(2, 2) + duidxi(3, 3)

      fi(:) = gravvec(:) * rhoi
      rm(1) = rhoi * (aci(1) + sum(duidxi(1, :) * uadvi(:))) - fi(1)
      rm(2) = rhoi * (aci(2) + sum(duidxi(2, :) * uadvi(:))) - fi(2)
      rm(3) = rhoi * (aci(3) + sum(duidxi(3, :) * uadvi(:))) - fi(3)

      ! calculate residual for NS and VOF
        rLi(:) = rhoi * aci(:) + &
        rhoi * (duidxi(:, 1) * uadvi(1) + duidxi(:, 2) * uadvi(2) + duidxi(:, 3) * uadvi(3)) + &
        dpridxi(:) - fi(:) - &
        mui * (duidxixj(:, 1, 1) + duidxixj(:, 2, 2) + duidxixj(:, 3, 3)) - &
        (mui+lambda) * (duidxixj(1, 1, :) + duidxixj(2, 2, :) + duidxixj(3, 3, :))

      ! Dissusion with Jinhui on Sep 1, 2023:
      ! The gradient of properies are not zero in multiphase flow. So we should not ignore
      ! the viscous term in the residual using linear finite element. Even though it might
      ! have minimal effect on the solution, it is still better to include it. Meanwhile,
      ! introduce a flag to control this option.
  
      !if(grad_prop_in_diff) then
      !  rLi(:) = rLi(:) - matmul(duidxi, dphidxi) * dmudphi - &
      !    matmul(dphidxi, duidxi) * dmudphi + &
      !    2d0/3d0 * dmudphi * divu * dphidxi(:)
      !endif

      res_phic_tmp1 = rphii + sum(uadvi(:) * dphidxi(:)) &
        + phii * vdot - mdot / rhoa
      !  + mdot / (rhoa*rhow) * rhoi

      ! calculate the STABLIZATION
      call e3STAB_3D_NSVOF(NSD, Gij, Delt, uadvi, rhoi, mui, tauM, tauP, tauC, tauLS)
      uprime(:) = -tauM * rLi(:)

      do aa = 1, NSHL
        shconv(aa) = sum(uadvi(:) * shgradgu(aa, :))
      enddo

      tauBar = 0.0d0
      !if(USE_TAUBAR) then
      !  call e3STAB_3D_NSVOF_TAUBAR(NSD, Gij, uadvi, uprime, tauBar)
      !endif

      k_dc = 0.0d0
      if(abs(NS_kdc_w) + abs(NS_kdc_a)  > 0.0d0) then
        call prop_interp(NS_kdc_w, NS_kdc_a, phic, ns_kdc)
        call e3DC_beta2(NSD, ns_kdc, Gij, rLi, k_dc)
      endif

      k_dc_phi = 0.0d0
      if(abs(LSC_kdc) > 0.0d0) then
        ! call e3DC_scalar(NSD, LSC_kdc, Gij, res_phic_tmp1, dphidxi, k_dc_phi)
        k_dc_phi = LSC_kdc * abs(res_phic_tmp1) / norm2(Gij)
      endif
      ! write(20000+myid,*) xi,res_phic_tmp1, k_dc_phi
        ! call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(ismaster) write(*,*) "Evaluate KDC", igauss, NGAUSS

      if(.true.) then
        call e3CFL(NSD, uadvi, Gij, Delt, cfl_loc)
        cfl = max(cfl, cfl_loc)
      endif
        ! call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(ismaster) write(*,*) "Evaluate CFL", igauss, NGAUSS
      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0 .and. is_fluid) then
        do bb = 1, NSHL
          do aa = 1, NSHL
            do i = 1, NSD
              do j = 1, NSD
                xKebe11((i-1)*NSD+j, aa, bb) = xKebe11((i-1)*NSD+j, aa, bb) &
                  + fact2 * shgradgu(aa, j) * shgradgu(bb, i) * (mui + k_dc) * DetJgw &      ! From viscosity
                  + fact2 * shgradgu(aa, i) * shgradgu(bb, j) * (tauC+lambda) * DetJgw       ! From LSIC and viscosity 
              enddo
                xKebe11((i-1)*NSD+i, aa, bb) = xKebe11((i-1)*NSD+i, aa, bb) &
                  + fact1 * shlu(aa) * rhoi * shlu(bb) * DetJgw &                            ! From time derivative
                  + fact2 * shlu(aa) * rhoi * shconv(bb) * DetJgw &                          ! From convection
                  + fact2 * (mui + k_dc) * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJgw & ! From viscosity
                  + fact1 * rhoi * shconv(aa) * tauM * rhoi * shlu(bb) * DetJgw &            ! From SUPG time derivative
                  + fact2 * rhoi * shconv(aa) * tauM * rhoi * shconv(bb) * DetJgw            ! From SUPG convection
            enddo
            xGebe(:, aa, bb) = xGebe(:, aa, bb) &
              - fact2 * shgradgu(aa, :) * shlu(bb) * DetJgw &                                ! From <-p, div(w)>
              + fact2 * rhoi * shconv(aa) * tauM * shgradgu(bb, :) * DetJgw                  ! From SUPG <-p, div(w)> 

            xDebe1(:, aa, bb) = xDebe1(:, aa, bb) &
              + fact2 * shlu(aa) * shgradgu(bb, :) * DetJgw &                                ! From <q, div(u)> 
              + shgradgu(aa, :) * tauM * rhoi * fact1 * shlu(bb) * DetJgw &                  ! From PSPG time derivative 
              + shgradgu(aa, :) * tauM * rhoi * fact2 * shconv(bb) * DetJgw                  ! From PSPG convection

            xMebe(aa, bb) = xMebe(aa, bb)  &
              + fact2 * tauM * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJgw               ! From PSPG <tauM grad(q), grad(p)>

            xLSebe(aa, bb) = xLSebe(aa, bb)  &
              + (shlu(aa) + tauls * shconv(aa)) * fact1 * shlu(bb) * DetJgw &                ! From time derivative 
              + (shlu(aa) + tauls * shconv(aa)) * fact2 * shconv(bb) * DetJgw                ! From convection
            ! xLSebe(aa, bb) = xLSebe(aa, bb) &                                                ! From mdot/(rhoa*rhow)*rhoi 
            !   - (shlu(aa) + tauls * shconv(aa)) * fact2 * drhodphi * mdot * shlu(bb) / (rhoa*rhow) * DetJgw
            ! xLSebe(aa, bb) = xLSebe(aa, bb) &
            !   - (shlu(aa) + tauls * shconv(aa)) * fact2 * dmdphii * rhoi * shlu(bb) / (rhoa*rhow) * DetJgw
            xLSebe(aa, bb) = xLSebe(aa, bb)  + &
              (shlu(aa) + tauls * shconv(aa)) * fact2 * shlu(bb) * vdot * DetJgw + &
              (shlu(aa) + tauls * shconv(aa)) * fact2 * phii * (1/rhoa-1/rhow) * dmdphii * shlu(bb) * DetJgw - &
              (shlu(aa) + tauls * shconv(aa)) * fact2 * dmdphii * shlu(bb) / rhoa * DetJgw
            xLSebe(aa, bb) = xLSebe(aa, bb) &
              + fact2 * k_dc_phi * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJgw           ! From kdc 
            xLSebe(aa, bb) = xLSebe(aa, bb) &                                                ! From LSIC 
              + fact2 * dmdphii * (1/rhoa-1/rhow) * shlu(aa) * tauC * dmdphii * (1/rhoa-1/rhow) * shlu(bb) * DetJgw


            xULSebe(:, aa, bb) = xULSebe(:, aa, bb) &                                        ! From LSIC 
              + fact2 * shgradgu(aa, :) * tauC * dmdphii * shlu(bb) * (1/rhoa-1/rhow) * DetJgw

           


            xPLSebe(aa, bb) = xPLSebe(aa, bb) &                                              ! From <-q, vdot>
              - fact2 * shlu(aa) * dmdphii * shlu(bb) * (1/rhoa - 1/rhow) * DetJgw  
            
            xLSUebe(:, aa, bb) = xLSUebe(:, aa, bb) &
              + fact2 * (shlu(aa) + tauls * shconv(aa)) * shlu(bb) * dphidxi(:) * DetJgw     ! From convection 
            xLSUebe(:, aa, bb) = xLSUebe(:, aa, bb) &                                        ! From LSIC  
              - dmdphii * (1/rhoa - 1/rhow) * fact2 * shlu(aa) * tauC * shgradgu(bb, :) * DetJgw
          end do
        end do
      end if

      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0 .and. is_fluid) then
        tmp1(:) = rhoi *aci(:) + &
          rhoi * uadvi(1) * duidxi(:, 1) + &
          rhoi * uadvi(2) * duidxi(:, 2) + &
          rhoi * uadvi(3) * duidxi(:, 3) - fi(:)
        do aa = 1, NSHL
          RHSu(:, aa) = RHSu(:, aa) - &
            shlu(aa) * tmp1(:) * DetJgw - &
            shgradgu(aa, :) * (-pri + tauC * divu) * DetJgw + &
            !!! Comment here
            shgradgu(aa, :) * tauC * vdot * DetJgw - &
            shgradgu(aa, :) * lambda * divu * DetJgw - & 
            (mui + k_dc) * sum(shgradgu(aa, :) * (duidxi(:, 1) + duidxi(1, :))) * DetJgw - &
            (mui + k_dc) * sum(shgradgu(aa, :) * (duidxi(:, 2) + duidxi(2, :))) * DetJgw - &
            (mui + k_dc) * sum(shgradgu(aa, :) * (duidxi(:, 3) + duidxi(3, :))) * DetJgw + &
            rhoi * shconv(aa) * uprime(:) * DetJgw
          RHSp(aa) = RHSp(aa) - &
            shlu(aa) * (divu - vdot) * DetJgw + &
            sum(shgradgu(aa, :) * uprime(:)) * DetJgw
          RHSphi(aa) = RHSphi(aa) - &
            (shlu(aa) + tauls * shconv(aa)) * res_phic_tmp1 * DetJgw - &
            k_dc_phi * sum(shgradgu(aa, :) * dphidxi(:)) * DetJgw + &
            !!! Comment here
            shlu(aa) * dmdphii * (1/rhoa-1/rhow) * tauC * (divu - vdot) * DetJgw
        enddo
      end if

    end do

    ! Apply Dirichlet BCs
    call BCLhs_3D(nshl, iel, xKebe11, xGebe, xDebe1, &
                  xMebe, Rhsu, Rhsp, &
                  xLSebe, xLSUebe, xULSebe, xPLSebe, Rhsphi)
    !    call MPI_barrier(MPI_COMM_WORLD, mpi_err)
    !if(ismaster) write(*,*) "Evaluate BC", NSHL
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      call FillSparseMat_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe, &
                            xLSebe, xLSUebe, xULSebe, xPLSebe)
    end if
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      call LocalToGlobalNSVOF_3D(RHSGu, RHSGp, RHSGls, &
                                 NNODE, NSD, NSHL, maxNSHL, &
                                 IEN(iel, :), rhsu, rhsp, rhsphi)
        ! call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(ismaster) write(*,*) "Evaluate add_local"
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
  deallocate (phil, rphil)                         ! 25
  deallocate (Tl, rTl)                             ! 27
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

end subroutine IntElmAss_NSVOF_Quenching_STAB


!======================================================================
!
!======================================================================
subroutine IntElmAss_Tem_Quenching_STAB(&
  dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
  acgmAlpha, pgAlpha, phigAlpha, dphidtgAlpha, &
  TgAlpha, rTgAlpha, &
  assemble_tensor_flag)
  use aAdjKeep
  use commonvars
  use commonpars
  use mpi

  implicit none

  integer, intent(in) :: assemble_tensor_flag
  real(8), intent(in) :: dgAlpha(NNODE, NSD), ugAlpha(NNODE, NSD), &
                         acgAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
                         acgmAlpha(NNODE, NSD), pgAlpha(NNODE), &
                         phigAlpha(NNODE), dphidtgAlpha(NNODE), &
                         TgAlpha(NNODE), rTgAlpha(NNODE)
  ! Local variables
  real(8), parameter :: damp = 0.5d0

  integer :: iel, igauss, aa, bb, i, j, hess_flag, idx, nshl, NGAUSS

  real(8) :: volm, vol_ex, uprime(NSD)
  real(8) :: tauM, tauC, tauBar, tauBar1, pri, kappa_str, k_dc, k_dc_phi, tauP, tauLS

  real(8), allocatable :: RHStem(:), xTebe(:, :)

  real(8), allocatable :: shlu(:), shgradgu(:, :), shhessgu(:, :, :)

  real(8), allocatable :: fl(:, :), dl(:, :), ul(:, :), wl(:), acl(:, :), &
                          uml(:, :), acml(:, :), pl(:), dlold(:, :), &
                          xl(:, :), dumb(:, :), phil(:), ulold(:, :), dphidtl(:)

  real(8), allocatable :: rTl(:), Tl(:)

  real(8), allocatable :: gp(:, :), gw(:)
  real(8), allocatable :: shconv(:), shconv_full(:), tmp(:), drdTi(:)

  real(8) :: di(NSD), ui(NSD), aci(NSD), umi(NSD), acmi(NSD), rphii, &
             fi(NSD), uadvi(NSD), xi(NSD), ddidxi(NSD, NSD), phii, &
             dphidxi(NSD), duidxi(NSD, NSD), dphidxidxj(NSD, NSD), &
             duidxixj(NSD, NSD, NSD), duiolddxi(NSD, NSD), &
             dxidx(NSD, NSD), dpridxi(NSD), rLi(NSD)
  real(8) :: uadvi_full(NSD)
  real(8) :: Gij(NSD, NSD), Ginv(NSD, NSD)
  real(8) :: cfl(2), gcfl(2), cfl_loc(2)
  real(8) :: Ftens(NSD, NSD), Stens(NSD, NSD), Ctens(NSD, NSD, NSD, NSD)
  !real(8) :: res_phic_tmp1, res_phic_tmp2
  real(8) :: Qi, Si(3, 3), Omegai(3, 3)

  !real(8) :: phi

  real(8) :: rTi, Ti, T0c, dTdxi(NSD), dTdxixj(NSD, NSD)
  real(8) :: res_tem1, res_tem2, tau_tem
  logical :: is_fluid
  real(8) :: fact1, fact2
  real(8) :: Se, mdot, vdot
  real(8) :: rhoi, mui, cpi, hki, rhocpi
  real(8) :: k_dc_tem
  real(8) :: dmdTi, phic, tmp1, dSedTi
  real(8) :: DetJgw

  fact1 = almi
  fact2 = alfi * gami * Delt

  volm = 0.0d0
  vol_ex = 0.0d0
  cfl = 0.0d0

  NGAUSS = -1
  NSHL = -1
!  rhsgq = 0d0
!  lhsgq = 0d0

  ! loop over elements
  do iel = 1, NELEM
!  write(*,*) "iel:",iel, NELEM
    if (NSHL /= ELMNSHL(iel)) then

      if (NSHL >= 0) then
        deallocate (RHStem, xTebe)
        deallocate (shlu, shgradgu, shhessgu)
        deallocate (fl, dl, ul, wl, acl)
        deallocate (uml, acml, pl, dlold)
        deallocate (xl, dumb, phil, ulold, dphidtl)
        deallocate (rTl, Tl)
        deallocate (gp, gw)
        deallocate (shconv, shconv_full, tmp, drdTi)
      end if

      NSHL = ELMNSHL(iel)
      NGAUSS = ELMNGAUSS(iel)
      allocate (RHStem(NSHL), xTebe(NSHL, NSHL))
      allocate (shlu(NSHL), shgradgu(NSHL, NSD), shhessgu(NSHL, NSD, NSD))
      allocate (fl(NSHL, NSD), dl(NSHL, NSD), ul(NSHL, NSD), wl(NSHL), acl(NSHL, NSD))
      allocate (uml(NSHL, NSD), acml(NSHL, NSD), pl(NSHL), dlold(NSHL, NSD))
      allocate (xl(NSHL, NSD), dumb(NSHL, NSD), phil(NSHL), ulold(NSHL, NSD), dphidtl(NSHL))
      allocate (rTl(NSHL), Tl(NSHL))
      allocate (gp(NGAUSS, NSD), gw(NGAUSS))
      allocate (shconv(NSHL), shconv_full(NSHL), tmp(NSHL), drdTi(NSHL))

      ! get Gaussian points and weights
      call genGPandGW(gp, gw, NGAUSS)
    end if
    is_fluid = ELM_ID(iel) == 102
    fl = 0.0d0

    ! Get local solution arrays
    do i = 1, NSHL
      idx = IEN(iel, i)
      xl(i, :) = xg(idx, :)
      dl(i, :) = dgAlpha(idx, :)
      dlold(i, :) = dgold(idx, :)
      wl(i) = wg(idx)
      ul(i, :) = ugAlpha(idx, :)
      ulold(i, :) = ugold(idx, :)
      acl(i, :) = acgAlpha(idx, :)
      uml(i, :) = ugmAlpha(idx, :)
      acml(i, :) = acgmAlpha(idx, :)
      pl(i) = pgAlpha(idx)
      phil(i) = phigAlpha(idx)
      dphidtl(i) = dphidtgAlpha(idx)
      Tl(i) = TgAlpha(idx)
      rTl(i) = rTgAlpha(idx)
    end do

    ! initialize local stiffness matrix
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      xTebe(:, :) = 0.0d0
    end if
    ! initialize local load vector
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      RHStem(:) = 0.0d0
    end if
    ! Loop over integration points (NGAUSS in each direction)
    do igauss = 1, NGAUSS

      ! Get Element Shape functions and their gradients
      ! initialize
      shlu = 0.0d0
      shgradgu = 0.0d0
      shhessgu = 0.0d0
  
      call eval_shape(nshl, iel, gp(igauss, :), xl, dl, wl, shlu, &
                      shgradgu, shhessgu, dxidx, Gij, Ginv, 0)

      DetJgw = DetJ * gw(igauss)
      ! Fluid
      call e3int_qr(NSHL, NSD, 1, shlu, phil, phii)
      call e3int_qr(NSHL, NSD, 1, shlu, rTl, rTi)
      call e3int_qr(NSHL, NSD, 1, shlu, Tl, Ti)
      call e3int_qr(NSHL, NSD, NSD, shlu, xl, xi(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, acl, aci(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, ul, ui(1))
      call e3int_qr(NSHL, NSD, NSD, shlu, uml, umi(1))

      call e3int_qr_grad(NSHL, NSD, 1, shgradgu, pl, dpridxi)
      call e3int_qr_grad(NSHL, NSD, 1, shgradgu, Tl, dTdxi)
      call e3int_qr_grad(NSHL, NSD, NSD, shgradgu, ul, duidxi)

      call e3int_qr_hess(NSHL, NSD, 1, shhessgu, Tl, dTdxixj)
      call e3int_qr_hess(NSHL, NSD, NSD, shhessgu, ul, duidxixj)

      phic = max(min(phii, 1.0d0), 0.0d0)
      ! phic = phii
      ! T0c = max(min(ti, 773d0), 348d0)
      T0c = Ti
      if(Ti > Ts) then
        mdot = c_evap * (1-phic) * rhow * (T0c - Ts) / Ts
        dmdTi = c_evap * (1-phic) * rhow / Ts * c_dmdot
      else
        mdot = c_cond * (phic) * rhoa * (T0c - Ts) / Ts
        dmdTi = c_cond * (phic) * rhoa / Ts * c_dmdot
      endif
      vdot = mdot / rhoa - mdot / rhow

      !if(Ti < 348d0 .or. Ti > 773d0) dmdTi = 0d0

      if(is_fluid) then
        call prop_interp(rhow, rhoa, phic, rhoi)
        call prop_interp(muw, mua, phic, mui)
        call prop_interp(cpw, cpa, phic, cpi)
        call prop_interp(kappaw, kappaa, phic, hki)
        call prop_interp(rhow * cpw, rhoa * cpa, phic, rhocpi)
      else ! solid
        rhoi = rhos
        mui = mus
        cpi = cps
        hki = kappas
        rhocpi = rhos * cps
        ! T0c = max(min(ti, 773d0), 348d0)
        ! T0c = T0c - 273.0
        ! rhoi = 2.79980834d+03 - 1.99121134d-01 * T0c
        ! cpi = 8.39321693d+02 + 4.57341083d-01 * T0c + 7.63296031d-04 * T0c ** 2
        ! hki = 1.62493037d+02 + 1.83277791d-01 * T0c - 5.31720761d-04 * T0c ** 2 + 2.46439207d-07 * T0c ** 3
        ! rhocpi = rhoi * cpi

      endif

      ! ALE Advective Velocity
      uadvi(:) = ui(:) - umi(:)
      ! call a function to get residual for Tem
      ! call e3int_rLi(NSD, rhoi, mui, aci, duidxi, uadvi, dpridxi, fi, duidxixj, rLi)


      ! call e3STAB_3D_NSVOF(NSD, Gij, Delt, uadvi, rhoi, mui, tauM, tauP, tauC, tauLS)
      ! uprime(:) = -tauM * rLi(:)
  
      Se = ((cpw - cpa) * (Ti - Ts) - Lh) * mdot 
      ! call e3int_restem(NSD, rTi, Ti, dTdxi, dTdxixj, uadvi, Se, rhocpi, hki, res_tem1)
      res_tem1 = rhocpi * (rTi + dot_product(uadvi(:), dTdxi(:))) - Se -&
        hki * (dTdxixj(1,1) + dTdxixj(2,2) + dTdxixj(3,3))

      call e3STAB_3D_TEM(NSD, Gij, uadvi, Delt, rhoi, cpi, hki, tau_tem)

      k_dc_tem = 0.0d0
      if(abs(Tem_kdc) > 0d0 .and. is_fluid) then
        k_dc_tem = Tem_kdc * abs(res_tem1) / norm2(Gij) / Tref
        ! if(k_dc_tem > 1d-3) then
        !   write(20000+myid,*) xi, k_dc_tem 
        ! endif
     endif

      ! uadvi_full(:) = uadvi(:) + uprime(:)
      ! do aa = 1, NSHL
      !   shconv(aa) = sum(uadvi(:)*shgradgu(aa, :))
      !   ! shconv_full(aa) = sum(uadvi_full(:)*shgradgu(aa, :))
      ! enddo
      shconv(:) = matmul(shgradgu(:, :), uadvi(:))
      ! tmp1 = rhocpi * vdot + (rhow*cpw-rhoa*cpa) * sum(uadvi_full(:) * dphidxi(:))
      dSedTi = (cpw - cpa) * mdot + ((cpw - cpa) * (Ti - Ts) - lh)* dmdTi
      ! drdTi(:) = rhocpi * (fact1 * shlu(:) + fact2 * shconv(:)) - fact2 * dSedTi * shlu(:)
      ! shconv_full(:) = shconv(:)
      tmp(:)  = shlu(:) + tau_tem * rhocpi * shconv(:)
      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
        ! e3LHS_3D_tem()
        ! xTebe(:, :) = 0.0d0
        if(is_fluid) then
          do aa = 1, NSHL
            do bb = 1,NSHL
              xTebe(aa, bb) = xTebe(aa, bb) &
                + shlu(aa) * rhocpi * (fact1 * shlu(bb) + fact2 * shconv(bb)) * DetJgw &
                - shlu(aa) * fact2 * dSedTi * shlu(bb) * DetJgw &
                + tau_tem * rhocpi * shconv(aa) * rhocpi &
                    * (fact1 * shlu(bb) + fact2 * shconv(bb)) * DetJgw &
                - tau_tem * rhocpi * shconv(aa) * fact2 * dSedTi * shlu(bb) * DetJgw &
                + fact2 * (hki + k_dc_tem) * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJgw
                
            enddo
          enddo
        else
          do aa = 1, NSHL
            do bb = 1,NSHL
              xTebe(aa, bb) = xTebe(aa, bb) &
                + shlu(aa) * rhocpi * fact1 * shlu(bb) * DetJgw &
                + fact2 * hki * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJgw
            enddo
          enddo
        endif
      end if

      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
        if(is_fluid) then
          RHSTem(:) = RHSTem(:) - shlu(:) * rhocpi * (rTi + dot_product(uadvi(:), dTdxi(:))) * DetJgw
          RHSTem(:) = RHSTem(:) + shlu(:) * Se * DetJgw
          RHSTem(:) = RHSTem(:) - rhocpi * tau_tem * shconv(:) * res_tem1 * DetJgw
          RHSTem(:) = RHSTem(:) - (hki + k_dc_tem) * shgradgu(:, 1) * dTdxi(1) * DetJgw
          RHSTem(:) = RHSTem(:) - (hki + k_dc_tem) * shgradgu(:, 2) * dTdxi(2) * DetJgw
          RHSTem(:) = RHSTem(:) - (hki + k_dc_tem) * shgradgu(:, 3) * dTdxi(3) * DetJgw
        else
          RHSTem(:) = RHSTem(:) - shlu(:) * rhocpi * rTi * DetJgw
          RHSTem(:) = RHSTem(:) - hki * shgradgu(:, 1) * dTdxi(1) * DetJgw 
          RHSTem(:) = RHSTem(:) - hki * shgradgu(:, 2) * dTdxi(2) * DetJgw 
          RHSTem(:) = RHSTem(:) - hki * shgradgu(:, 3) * dTdxi(3) * DetJgw 
        endif
      end if
    end do

    ! Apply Dirichlet BCs
    call BCLHS_tem(nshl, iel, xTebe, RHSTem)
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      ! Assemble global LHS Matrix
      call FillSparseMat_tem(nshl, iel, xTebe)
    end if
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      ! Assemble load RHS vector
      do aa = 1, NSHL
        RHSGTEM(IEN(iel, aa)) = RHSGTEM(IEN(iel, aa)) + Rhstem(aa)
      end do
    end if
  end do

  deallocate (RHStem, xTebe)
  deallocate (shlu, shgradgu, shhessgu)
  deallocate (fl, dl, ul, wl, acl)
  deallocate (uml, acml, pl, dlold)
  deallocate (xl, dumb, phil, ulold, dphidtl)
  deallocate (rTl, Tl)
  deallocate (gp, gw)
  deallocate (shconv, shconv_full, tmp, drdTi)


end subroutine IntElmAss_Tem_Quenching_STAB

