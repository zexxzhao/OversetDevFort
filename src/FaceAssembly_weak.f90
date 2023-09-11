!======================================================================
!
!======================================================================
subroutine FaceAssembly_NS_weak(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                                acgmAlpha, pgAlpha, &
                                assemble_tensor_flag, subdomain_flag)
  use aAdjKeep
  use mpi
  use commonvars
  use commonpars
  implicit none

  integer, intent(in) :: assemble_tensor_flag, subdomain_flag
  real(8), intent(in) :: dgAlpha(NNODE, NSD), ugAlpha(NNODE, NSD), &
                         acgAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
                         acgmAlpha(NNODE, NSD), pgAlpha(NNODE)

  ! Local variables
  real(8), allocatable :: shlu(:), shgradgu(:, :), shhessgu(:, :, :)

  real(8), allocatable :: dl(:, :), ul(:, :), acl(:, :), &
                          uml(:, :), acml(:, :), pl(:), &
                          xl(:, :), phil(:), rphil(:), wl(:), dphidtl(:), phi_bgl(:)
  real(8), allocatable :: xKebe11(:, :, :), xMebe(:, :), xGebe(:, :, :), &
                          xDebe1(:, :, :), xDebe2(:, :, :), &
                          xLSebe(:, :), xLSUebe(:, :, :), &
                          xULSebe(:, :, :), xPLSebe(:, :)
  real(8), allocatable :: Rhsu(:, :), Rhsp(:), Rhsm(:, :), Rhsphi(:)

  integer :: b, ifac, iel, igauss, i, j, k, hess_flag, aa, bb, nshl, outflow, traction

  real(8), allocatable :: gp(:, :), gw(:), mgp(:, :)
  real(8) :: ui(NSD), pri, duidxi(NSD, NSD), phii,phic, rphii, dphidxi(NSD), &
             Gij(NSD, NSD), Ginv(NSD, NSD), ti(NSD), nor(NSD), &
             bui(NSD), umi(NSD), dxidx(NSD, NSD), dphidx(NSD), xi(3)
  real(8) :: tmp1(NSD), tmp2(NSD, NSD)
  real(8) :: fact1, fact2, tauB, tauNor, gi(NSD), gphi, tauBLS
  real(8) :: gneg, wave_u, wave_phi
  real(8) :: e3(3)
  real(8) :: rhoi, mui
  integer :: NGAUSSf
  logical :: is_target
  e3 = 0d0
  e3(3) = -1d0
  ! "shb" will be the shape function array while "shbg" will hold the
  ! gradients of the shape functions
  fact1 = almi
  fact2 = alfi*gami*Delt

  outflow = 0 !not outflow boundary
  traction = 0 ! no traction BC

  gi = 0.0d0

  NSHL = -1

  ! loop over faces
  do b = 1, NBOUND

    if (BCugType(b, 1) /= 2 .and. &
        BCugType(b, 2) /= 2 .and. &
        BCugType(b, 3) /= 2) then
      cycle
    endif


    do ifac = 1, bound(b)%NFACE

      iel = bound(b)%F2E(ifac)

      if (NSHL /= ELMNSHL(iel)) then
        if (NSHL >= 0) then
          deallocate (shlu, shgradgu, shhessgu, &
                      dl, ul, acl, uml, acml, pl, xl, phil, rphil, wl, phi_bgl, &
                      xKebe11, xMebe, xGebe, xDebe1, xDebe2, &
                      Rhsu, Rhsp, Rhsm, xLSebe, xLSUebe, xULSebe, xPLSebe, Rhsphi)
          deallocate (gp, gw, mgp)
        end if

        NSHL = ELMNSHL(iel)
        NGAUSSf = bound(b)%NGAUSSB(ifac)

        allocate (shlu(NSHL), shgradgu(NSHL, NSD), shhessgu(NSHL, NSD, NSD), &
                  dl(NSHL, NSD), ul(NSHL, NSD), acl(NSHL, NSD), uml(NSHL, NSD), &
                  acml(NSHL, NSD), pl(NSHL), xl(NSHL, NSD), phil(NSHL), rphil(NSHL), &
                  wl(NSHL), &
                  xKebe11(NSD*NSD, NSHL, NSHL), xMebe(NSHL, NSHL), &
                  xGebe(NSD, NSHL, NSHL), xDebe1(NSD, NSHL, NSHL), &
                  xDebe2(NSD, NSHL, NSHL), &
                  Rhsu(NSD, NSHL), Rhsp(NSHL), Rhsm(NSD, NSHL), &
                  xLSebe(NSHL, NSHL), xLSUebe(NSD, NSHL, NSHL), &
                  xULSebe(NSD, NSHL, NSHL), xPLSebe(NSHL, NSHL), Rhsphi(NSHL), phi_bgl(NSHL))
        allocate (gp(NGAUSSf, 2), gw(NGAUSSf), mgp(NGAUSSf, 3))
      end if

      call genGPandGWb(gp, gw, NGAUSSf)
      call genGPMap(nshl, NGAUSSf, bound(b)%FACE_OR(ifac), iga, mgp)

      is_target = .false.
      if(iand(subdomain_flag, 1) > 0) then
        is_target = is_target .or. ELM_ID(iel) == ELM_ID_SET(1)
      endif
      if(iand(subdomain_flag, 2) > 0) then
        is_target = is_target .or. ELM_ID(iel) == ELM_ID_SET(2)
      endif

      ! Get local solution arrays
      do i = 1, NSHL
        j = IEN(iel, i)
        xl(i, :) = xg(j, :)
        dl(i, :) = dgAlpha(j, :)
        ul(i, :) = ugAlpha(j, :)
        acl(i, :) = acgAlpha(j, :)
        uml(i, :) = ugmAlpha(j, :)
        acml(i, :) = acgmAlpha(j, :)
        pl(i) = pgAlpha(j)
        wl(i) = wg(j)
      end do

      ! initialize local resistance matrix
      xKebe11 = 0.0d0
      xGebe = 0.0d0
      xDebe1 = 0.0d0
      xMebe = 0.0d0

      ! initialize local load vector
      Rhsu = 0.0d0
      Rhsp = 0.0d0

      do igauss = 1, NGAUSSf

        call eval_faceshape(nshl, iel, gp(igauss, :), mgp(igauss, :), &
                            bound(b)%FACE_OR(ifac), &
                            xl, dl, wl, shlu, shgradgu, dxidx, &
                            Gij, Ginv, nor)

        ! Interpolate
        call e3int_qr(NSHL, NSD, 1, shlu, pl, pri)
        call e3int_qr(NSHL, NSD, NSD, shlu, ul, ui(1))
        call e3int_qr(NSHL, NSD, NSD, shlu, uml, umi(1))
        call e3int_qr(NSHL, NSD, NSD, shlu, xl, xi(1))

        call e3int_qr_grad(NSHL, NSD, NSD, shgradgu, ul, duidxi)

        rhoi = rhoa
        mui = mua
        gi = umi


        if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0 .and. is_target) then
          call e3bSTAB_weak(tauB, tauNor, ui - umi, nor, dxidx, rhoi, mui)

          call e3bRHS_weak( &
            nshl, rhoi, mui, nor, tauB, tauNor, gw(igauss), &
            shlu, shgradgu, ui, umi, pri, duidxi, &
            gi, Rhsu, Rhsp, ti, tmp1, tmp2)

        end if

        if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0 .and. is_target) then
          call e3bLHS_weak( &
            nshl, rhoi, mui, ui, umi, duidxi, &
            tauB, tauNor, &
            gw(igauss), shlu, shgradgu, xKebe11, &
            xGebe, xDebe1, nor)
        end if

      end do
      call BCLhs_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe, Rhsu, Rhsp)
      if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
        call FillSparseMat_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe)
      endif
      if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
        do bb = 1, NSHL
          RHSGu(IEN(iel, bb), :) = RHSGu(IEN(iel, bb), :) + Rhsu(:, bb)
          RHSGp(IEN(iel, bb)) = RHSGp(IEN(iel, bb)) + Rhsp(bb)
        enddo
      endif

    end do


  end do

  if (NSHL >= 0) then
    deallocate (shlu, shgradgu, shhessgu, &
                dl, ul, acl, uml, acml, xl, phil, rphil, wl, &
                xKebe11, xMebe, xGebe, xDebe1, xDebe2, &
                Rhsu, Rhsp, Rhsm, Rhsphi, phi_bgl, xLSebe, xULSebe, xLSUebe, xPLSebe, &
                gp, gw, mgp)
  end if

end subroutine FaceAssembly_NS_weak
!======================================================================
!
!======================================================================
subroutine FaceAssembly_NS_outflow(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                                acgmAlpha, pgAlpha, &
                                assemble_tensor_flag, subdomain_flag)
  use aAdjKeep
  use mpi
  use commonvars
  use commonpars
  implicit none

  integer, intent(in) :: assemble_tensor_flag, subdomain_flag
  real(8), intent(in) :: dgAlpha(NNODE, NSD), ugAlpha(NNODE, NSD), &
                         acgAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
                         acgmAlpha(NNODE, NSD), pgAlpha(NNODE)

  ! Local variables
  real(8), allocatable :: shlu(:), shgradgu(:, :), shhessgu(:, :, :)

  real(8), allocatable :: dl(:, :), ul(:, :), acl(:, :), &
                          uml(:, :), acml(:, :), pl(:), &
                          xl(:, :), phil(:), rphil(:), wl(:), dphidtl(:), phi_bgl(:)
  real(8), allocatable :: xKebe11(:, :, :), xMebe(:, :), xGebe(:, :, :), &
                          xDebe1(:, :, :), xDebe2(:, :, :), &
                          xLSebe(:, :), xLSUebe(:, :, :), &
                          xULSebe(:, :, :), xPLSebe(:, :)
  real(8), allocatable :: Rhsu(:, :), Rhsp(:), Rhsm(:, :), Rhsphi(:)

  integer :: b, ifac, iel, igauss, i, j, k, hess_flag, aa, bb, nshl, outflow, traction
  integer :: ii

  real(8), allocatable :: gp(:, :), gw(:), mgp(:, :)
  real(8) :: ui(NSD), pri, duidxi(NSD, NSD), phii,phic, rphii, dphidxi(NSD), &
             Gij(NSD, NSD), Ginv(NSD, NSD), ti(NSD), nor(NSD), &
             bui(NSD), umi(NSD), dxidx(NSD, NSD), dphidx(NSD), xi(3)
  real(8), allocatable :: tmp2(:, :)
  real(8) :: tmp1(NSD)
  real(8) :: fact1, fact2, tauB, tauNor, gi(NSD), gphi, tauBLS
  real(8) :: gneg, wave_u, wave_phi
  real(8) :: gnor, unor, upos, uneg
  real(8) :: e3(3)
  real(8) :: rhoi, mui
  real(8) :: DetJgw
  integer :: NGAUSSf
  logical :: is_target
  e3 = 0d0
  e3(3) = -1d0
  ! "shb" will be the shape function array while "shbg" will hold the
  ! gradients of the shape functions
  fact1 = almi
  fact2 = alfi*gami*Delt

  outflow = 0 !not outflow boundary
  traction = 0 ! no traction BC

  gi = 0.0d0

  NSHL = -1

  ! loop over faces
  do b = 1, NBOUND

    if (BCugType(b, 1) /= 2 .and. &
        BCugType(b, 2) /= 2 .and. &
        BCugType(b, 3) /= 2) then
      cycle
    endif


    do ifac = 1, bound(b)%NFACE

      iel = bound(b)%F2E(ifac)

      if (NSHL /= ELMNSHL(iel)) then
        if (NSHL >= 0) then
          deallocate (shlu, shgradgu, shhessgu, &
                      dl, ul, acl, uml, acml, pl, xl, phil, rphil, wl, phi_bgl, &
                      xKebe11, xMebe, xGebe, xDebe1, xDebe2, &
                      Rhsu, Rhsp, Rhsm, xLSebe, xLSUebe, xULSebe, xPLSebe, Rhsphi)
          deallocate (gp, gw, mgp)
          deallocate (tmp2)
        end if

        NSHL = ELMNSHL(iel)
        NGAUSSf = bound(b)%NGAUSSB(ifac)

        allocate (shlu(NSHL), shgradgu(NSHL, NSD), shhessgu(NSHL, NSD, NSD), &
                  dl(NSHL, NSD), ul(NSHL, NSD), acl(NSHL, NSD), uml(NSHL, NSD), &
                  acml(NSHL, NSD), pl(NSHL), xl(NSHL, NSD), phil(NSHL), rphil(NSHL), &
                  wl(NSHL), &
                  xKebe11(NSD*NSD, NSHL, NSHL), xMebe(NSHL, NSHL), &
                  xGebe(NSD, NSHL, NSHL), xDebe1(NSD, NSHL, NSHL), &
                  xDebe2(NSD, NSHL, NSHL), &
                  Rhsu(NSD, NSHL), Rhsp(NSHL), Rhsm(NSD, NSHL), &
                  xLSebe(NSHL, NSHL), xLSUebe(NSD, NSHL, NSHL), &
                  xULSebe(NSD, NSHL, NSHL), xPLSebe(NSHL, NSHL), Rhsphi(NSHL), phi_bgl(NSHL))
        allocate (gp(NGAUSSf, 2), gw(NGAUSSf), mgp(NGAUSSf, 3))
        allocate (tmp2(NSHL, NSHL))
      end if

      call genGPandGWb(gp, gw, NGAUSSf)
      call genGPMap(nshl, NGAUSSf, bound(b)%FACE_OR(ifac), iga, mgp)

      is_target = .false.
      if(iand(subdomain_flag, 1) > 0) then
        is_target = is_target .or. ELM_ID(iel) == ELM_ID_SET(1)
      endif
      if(iand(subdomain_flag, 2) > 0) then
        is_target = is_target .or. ELM_ID(iel) == ELM_ID_SET(2)
      endif

      ! Get local solution arrays
      do i = 1, NSHL
        j = IEN(iel, i)
        xl(i, :) = xg(j, :)
        dl(i, :) = dgAlpha(j, :)
        ul(i, :) = ugAlpha(j, :)
        acl(i, :) = acgAlpha(j, :)
        uml(i, :) = ugmAlpha(j, :)
        acml(i, :) = acgmAlpha(j, :)
        pl(i) = pgAlpha(j)
        wl(i) = wg(j)
      end do

      ! initialize local resistance matrix
      xKebe11 = 0.0d0
      xGebe = 0.0d0
      xDebe1 = 0.0d0
      xMebe = 0.0d0

      ! initialize local load vector
      Rhsu = 0.0d0
      Rhsp = 0.0d0

      do igauss = 1, NGAUSSf

        call eval_faceshape(nshl, iel, gp(igauss, :), mgp(igauss, :), &
                            bound(b)%FACE_OR(ifac), &
                            xl, dl, wl, shlu, shgradgu, dxidx, &
                            Gij, Ginv, nor)
        DetJgw = DetJb*gw(igauss)

        ! Interpolate
        call e3int_qr(NSHL, NSD, 1, shlu, pl, pri)
        call e3int_qr(NSHL, NSD, NSD, shlu, ul, ui(1))
        call e3int_qr(NSHL, NSD, NSD, shlu, uml, umi(1))
        call e3int_qr(NSHL, NSD, NSD, shlu, xl, xi(1))

        call e3int_qr_grad(NSHL, NSD, NSD, shgradgu, ul, duidxi)

        rhoi = rhoa
        mui = mua
        gi = umi

        unor = sum((ui - umi)*nor)  ! u \cdot n
        upos = 0.5d0*(unor + abs(unor))
        uneg = 0.5d0*(unor - abs(unor))

        ! Absolute normal for the rest
        gnor = sum(gi*nor)
        unor = sum(ui*nor)  ! u \cdot n


        if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
          do aa = 1, NSHL
            RHSu(:, aa) = RHSu(:, aa) + &
              shlu(aa) * rhoi * uneg * (ui(:) - gi(:)) * DetJgw
          enddo

        end if

        if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
          do bb = 1, NSHL
            do aa = 1, NSHL
              tmp2(aa, bb) = - shlu(aa)*uneg*rhoi*shlu(bb)
            enddo
          enddo

          do bb = 1, NSHL
            do aa = 1, NSHL
              do ii = 1, NSD
                xKebe11((ii-1)*NSD+ii, aa, bb) = xKebe11((ii-1)*NSD+ii, aa, bb) +&
                    fact2 * tmp2(aa, bb) * DetJgw
              enddo
            enddo
          enddo
        end if

      end do
      call BCLhs_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe, Rhsu, Rhsp)
      if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
        call FillSparseMat_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe)
      endif
      if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
        do bb = 1, NSHL
          RHSGu(IEN(iel, bb), :) = RHSGu(IEN(iel, bb), :) + Rhsu(:, bb)
          RHSGp(IEN(iel, bb)) = RHSGp(IEN(iel, bb)) + Rhsp(bb)
        enddo
      endif

    end do


  end do

  if (NSHL >= 0) then
    deallocate (shlu, shgradgu, shhessgu, &
                dl, ul, acl, uml, acml, xl, phil, rphil, wl, &
                xKebe11, xMebe, xGebe, xDebe1, xDebe2, &
                Rhsu, Rhsp, Rhsm, Rhsphi, phi_bgl, xLSebe, xULSebe, xLSUebe, xPLSebe, &
                gp, gw, mgp)
    deallocate (tmp2)
  end if

end subroutine FaceAssembly_NS_outflow

!======================================================================
!
!======================================================================

