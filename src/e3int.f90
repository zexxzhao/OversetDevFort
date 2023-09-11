!======================================================================
!
!======================================================================
subroutine prop_interp(prop_a, prop_b, phi, prop_out)

  implicit none
  real(8), intent(in) :: prop_a, prop_b, phi
  real(8), intent(out) :: prop_out

  real(8) :: phi_clipped
  phi_clipped = max(0.0d0, min(1.0d0, phi))
  prop_out = prop_a * (1d0 - phi_clipped) + prop_b * phi_clipped
end subroutine prop_interp

!======================================================================
!
!======================================================================
subroutine e3int_qr(NSHL, NSD, NVAR, shlu, vall, vali)
  implicit none
  integer, intent(in) :: NSHL, NSD, NVAR
  real(8), intent(in) :: shlu(NSHL)
  real(8), intent(in) :: vall(NSHL, NVAR)
  real(8), intent(out) :: vali(NVAR)

  integer :: i
  do i = 1, NVAR
    vali(i) = sum(vall(:, i)*shlu(:))
  enddo
end subroutine e3int_qr

!======================================================================
!
!======================================================================
subroutine e3int_qr_grad(NSHL, NSD, NVAR, shgradlu, vall, gradvali)
  implicit none
  integer, intent(in) :: NSHL, NSD, NVAR
  real(8), intent(in) :: shgradlu(NSHL, NSD)
  real(8), intent(in) :: vall(NSHL, NVAR)
  real(8), intent(out) :: gradvali(NVAR, NSD)

  integer :: i, j
  do i = 1, NVAR
    do j = 1, NSD
      gradvali(i, j) = sum(vall(:, i)*shgradlu(:, j))
    enddo
  enddo
end subroutine e3int_qr_grad

!======================================================================
!
!======================================================================
subroutine e3int_qr_hess(NSHL, NSD, NVAR, shhesslu, vall, hessvali)
  implicit none
  integer, intent(in) :: NSHL, NSD, NVAR
  real(8), intent(in) :: shhesslu(NSHL, NSD, NSD)
  real(8), intent(in) :: vall(NSHL, NVAR)
  real(8), intent(out) :: hessvali(NVAR, NSD, NSD)

  integer :: i, j, k
  do i = 1, NVAR
    do j = 1, NSD
      do k = 1, NSD
        hessvali(i, j, k) = sum(vall(:, i)*shhesslu(:, j, k))
      enddo
    enddo
  enddo
end subroutine e3int_qr_hess
!======================================================================
!
!======================================================================
subroutine e3int_fluid(nshl, xl, dl, ul, acl, uml, acml, pl, fl, phil, &
                       shlu, shgradlu, shgradgu, shhessgu, dxidx, &
                       Ginv, di, ui, aci, umi, acmi, pri, fi, ddidxi, &
                       duidxi, duidxixj, dpridxi, phi, &
                       dphidxi, dphidxixj, dphidtl, dphidti, xi, &
                       rhoi, mui, &
                       rLi, res_phi)
  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl
  real(8), intent(in) :: Ginv(NSD, NSD)

  integer :: i, j, k
  real(8) :: shlu(NSHL), shgradlu(NSHL, NSD), &
             shgradgu(NSHL, NSD), shhessgu(NSHL, NSD, NSD), dxidx(NSD, NSD)
  real(8) :: dl(NSHL, NSD), ul(NSHL, NSD), xl(NSHL, NSD), &
             acl(NSHL, NSD), uml(NSHL, NSD), acml(NSHL, NSD), pl(NSHL), &
             fl(NSHL, NSD), phil(NSHL)
  real(8) :: di(NSD), ui(NSD), aci(NSD), umi(NSD), acmi(NSD), pri, &
             fi(NSD), ddidxi(NSD, NSD), &
             duidxi(NSD, NSD), duidxixj(NSD, NSD, NSD), &
             dpridxi(NSD), phi, dphidxi(NSD), rLi(NSD), xi(NSD)
  real(8) :: uadv(NSD), He0, rho0, h, dphidxi0(NSD), dphidxixj(NSD, NSD), &
             dphidtl(NSHL), dphidti
  real(8) :: rhoi, mui, res_phi

  xi = 0.0d0
  di = 0.0d0
  ui = 0.0d0
  aci = 0.0d0
  umi = 0.0d0
  acmi = 0.0d0
  pri = 0.0d0
  fi = 0.0d0
  phi = 0.0d0
  dphidti = 0.0d0

  ddidxi = 0.0d0
  duidxi = 0.0d0
  duidxixj = 0.0d0
  dpridxi = 0.0d0
  dphidxi = 0.0d0
  dphidxixj = 0.0d0

  rLi = 0.0d0

  ! Get quantities and their grads and hessians
  do i = 1, NSHL
    xi(:) = xi(:) + xl(i, :)*shlu(i)
    di(:) = di(:) + dl(i, :)*shlu(i)
    ui(:) = ui(:) + ul(i, :)*shlu(i)
    aci(:) = aci(:) + acl(i, :)*shlu(i)
    umi(:) = umi(:) + uml(i, :)*shlu(i)
    acmi(:) = acmi(:) + acml(i, :)*shlu(i)
    pri = pri + pl(i)*shlu(i)
    fi(:) = fi(:) + fl(i, :)*shlu(i)
    phi = phi + phil(i)*shlu(i)
    dpridxi(:) = dpridxi(:) + pl(i)*shgradgu(i, :)
    dphidxi(:) = dphidxi(:) + phil(i)*shgradgu(i, :)
    dphidti = dphidti + dphidtl(i)*shlu(i)
    do j = 1, NSD
      ddidxi(:, j) = ddidxi(:, j) + dl(i, :)*shgradgu(i, j)
      duidxi(:, j) = duidxi(:, j) + ul(i, :)*shgradgu(i, j)
      do k = 1, NSD
        duidxixj(:, j, k) = duidxixj(:, j, k) + ul(i, :)*shhessgu(i, j, k)
        dphidxixj(j, k) = dphidxixj(j, k) + phil(i)*shhessgu(i, j, k)
      end do
    end do
  end do

  rhoi = rhow * phi + rhoa * (1.0d0 - phi)
  mui = muw * phi + rhoa * (1.0d0 - phi)

!  ! NS-ALE PDE Residual
!  dphidxi0(1) = 0.0d0
!  dphidxi0(2) = 0.0d0
!  dphidxi0(3) = 1.0d0
!  call getElemSize(h,dxidx,dphidxi0,Ginv)
!  call getHeps(He0,water_level-xi(3)-di(3), h)
!  rho0 = (1.0d0-He0)*rhoa + He0*rhow

!  fi = fi + (rho-rho0)*gravvec
  fi(:) = rhoi * gravvec(:)
  uadv(:) = ui(:) - umi(:)

  rLi(:) = rhoi*aci(:) + &
           rhoi*(duidxi(:, 1)*uadv(1) + &
                 duidxi(:, 2)*uadv(2) + &
                 duidxi(:, 3)*uadv(3)) + &
           dpridxi(:) - fi(:) - &
           mui*(duidxixj(:, 1, 1) + duidxixj(:, 2, 2) + duidxixj(:, 3, 3)) - &
           mui*2d0/3d0*(duidxixj(1, 1, :) + duidxixj(2, 2, :) + duidxixj(3, 3, :))
end subroutine e3int_fluid

!======================================================================
!
!======================================================================

subroutine e3int_rLi(NSD, rhoi, mui, &
                    aci, duidxi, uadvi, &
                    dpridxi, fi, duidxixj, &
                    rLi)
  implicit none
  integer, intent(in) :: NSD
  real(8), intent(in) :: rhoi, mui
  real(8), intent(in) :: aci(NSD), duidxi(NSD, NSD), uadvi(NSD)
  real(8), intent(in) :: dpridxi(NSD), fi(NSD)
  real(8), intent(in) :: duidxixj(NSD, NSD, NSD)

  real(8), intent(out) :: rLi(NSD)
  rLi(:) = rhoi*aci(:) + &
           rhoi*(duidxi(:, 1)*uadvi(1) + &
                 duidxi(:, 2)*uadvi(2) + &
                 duidxi(:, 3)*uadvi(3)) + &
           dpridxi(:) - fi(:) - &
           mui*(duidxixj(:, 1, 1) + duidxixj(:, 2, 2) + duidxixj(:, 3, 3)) - &
           mui*2d0/3d0*(duidxixj(1, 1, :) + duidxixj(2, 2, :) + duidxixj(3, 3, :))
end subroutine e3int_rLi
!======================================================================
!
!======================================================================

subroutine e3int_resphi(NSD, rphii, phii, dphidxi, uadvi, mdot, divu, rhoa, resphi)
  implicit none
  integer, intent(in) :: NSD
  real(8), intent(in) :: rphii, phii, dphidxi(NSD)
  real(8), intent(in) :: uadvi(NSD)
  real(8), intent(in) :: mdot, divu, rhoa
  
  real(8), intent(out) :: resphi

  !real(8) :: vdot

  resphi = rphii + sum(uadvi(:) * dphidxi(:)) + phii * divu - mdot / rhoa
  ! resphi = rphii + sum(uadvi(:) * dphidxi(:)) - mdot / rhoa
end subroutine e3int_resphi
!======================================================================
!
!======================================================================

subroutine e3int_restem(NSD, rTi, Ti, dTdxi, dTdxixj, uadvi, Se, rhocpi, hki, restem)
  implicit none
  integer, intent(in) :: NSD
  real(8), intent(in) :: rTi, Ti, dTdxi(NSD), dTdxixj(NSD, NSD)
  real(8), intent(in) :: uadvi(NSD)
  real(8), intent(in) :: Se, rhocpi, hki
  real(8), intent(out) :: restem

  real(8) :: lapT
  integer :: i
  ! calculate laplacian of temperature
  lapT = 0.0d0
  do i = 1, NSD
    lapT = lapT + dTdxixj(i, i)
  enddo

  restem = rhocpi*(rTi + sum(uadvi(:)*dTdxi(:))) - Se - hki * lapT
end subroutine e3int_restem
!======================================================================
!
!======================================================================
subroutine e3int_struct(nshl, dl, ul, acl, fl, shlu, shgradgu, &
                        di, ui, aci, fi, ddidxi)
  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl

  integer :: i, j, k
  real(8) :: shlu(NSHL), shgradgu(NSHL, NSD)
  real(8) :: dl(NSHL, NSD), ul(NSHL, NSD), &
             acl(NSHL, NSD), fl(NSHL, NSD)
  real(8) :: di(NSD), ui(NSD), aci(NSD), fi(NSD), ddidxi(NSD, NSD)

  di = 0.0d0
  ui = 0.0d0
  aci = 0.0d0
  fi = 0.0d0
  ddidxi = 0.0d0

  ! Get quantities and their grads and hessians
  do i = 1, NSHL
    di(:) = di(:) + dl(i, :)*shlu(i)
    ui(:) = ui(:) + ul(i, :)*shlu(i)
    aci(:) = aci(:) + acl(i, :)*shlu(i)
    fi(:) = fi(:) + fl(i, :)*shlu(i)
    do j = 1, NSD
      ddidxi(:, j) = ddidxi(:, j) + dl(i, :)*shgradgu(i, j)
    end do
  end do

end subroutine e3int_struct

!======================================================================
!
!======================================================================
subroutine e3tensors_nonlin(NSD, ddidxi, mu, lambda, Ftens, Stens, Ctens)

  implicit none

  integer :: i, j, k, l, NSD
  real(8) :: ddidxi(NSD, NSD)
  real(8) :: Ftens(NSD, NSD), Stens(NSD, NSD), Ctens(NSD, NSD, NSD, NSD)
  real(8) :: Itens(NSD, NSD), Jdet, Finv(NSD, NSD), tmp, &
             Ct(NSD, NSD), CtInv(NSD, NSD), Jmod, trC, pt33, pt5, &
             fact1, fact2, fact3, mu, lambda, Ctemp(NSD, NSD, NSD, NSD)

  pt33 = 1.0d0/3.0d0
  pt5 = 0.5d0

  Itens = 0.0d0
  Itens(1, 1) = 1.0d0
  Itens(2, 2) = 1.0d0
  Itens(3, 3) = 1.0d0

  Ftens = Itens + ddidxi    ! F = I + du/dX

  Jdet = 0.0d0
  Finv = 0.0d0
  tmp = 0.0d0

  ! compute the inverse of deformation gradient and Jacobian
  Finv(1, 1) = Ftens(2, 2)*Ftens(3, 3) &
               - Ftens(3, 2)*Ftens(2, 3)
  Finv(1, 2) = Ftens(3, 2)*Ftens(1, 3) &
               - Ftens(1, 2)*Ftens(3, 3)
  Finv(1, 3) = Ftens(1, 2)*Ftens(2, 3) &
               - Ftens(1, 3)*Ftens(2, 2)
  tmp = 1.0d0/(Finv(1, 1)*Ftens(1, 1) &
               + Finv(1, 2)*Ftens(2, 1) &
               + Finv(1, 3)*Ftens(3, 1))
  Finv(1, 1) = Finv(1, 1)*tmp
  Finv(1, 2) = Finv(1, 2)*tmp
  Finv(1, 3) = Finv(1, 3)*tmp
  Finv(2, 1) = (Ftens(2, 3)*Ftens(3, 1) &
                - Ftens(2, 1)*Ftens(3, 3))*tmp
  Finv(2, 2) = (Ftens(1, 1)*Ftens(3, 3) &
                - Ftens(3, 1)*Ftens(1, 3))*tmp
  Finv(2, 3) = (Ftens(2, 1)*Ftens(1, 3) &
                - Ftens(1, 1)*Ftens(2, 3))*tmp
  Finv(3, 1) = (Ftens(2, 1)*Ftens(3, 2) &
                - Ftens(2, 2)*Ftens(3, 1))*tmp
  Finv(3, 2) = (Ftens(3, 1)*Ftens(1, 2) &
                - Ftens(1, 1)*Ftens(3, 2))*tmp
  Finv(3, 3) = (Ftens(1, 1)*Ftens(2, 2) &
                - Ftens(1, 2)*Ftens(2, 1))*tmp

  Jdet = 1.0d0/tmp       ! J
  Jmod = Jdet**(-2.0d0/3.0d0) ! J^-2/3
!!!  Jmod = Jdet

  Ct = 0.0d0

  do i = 1, NSD
    do j = 1, NSD
      Ct(i, j) = sum(Ftens(:, i)*Ftens(:, j))
    end do
  end do

  ! compute the inverse of the F^T F tensor
  CtInv = 0.0d0
  tmp = 0.0d0
  trC = 0.0d0

  CtInv(1, 1) = Ct(2, 2)*Ct(3, 3) &
                - Ct(3, 2)*Ct(2, 3)
  CtInv(1, 2) = Ct(3, 2)*Ct(1, 3) &
                - Ct(1, 2)*Ct(3, 3)
  CtInv(1, 3) = Ct(1, 2)*Ct(2, 3) &
                - Ct(1, 3)*Ct(2, 2)
  tmp = 1.0d0/(CtInv(1, 1)*Ct(1, 1) &
               + CtInv(1, 2)*Ct(2, 1) &
               + CtInv(1, 3)*Ct(3, 1))
  CtInv(1, 1) = CtInv(1, 1)*tmp
  CtInv(1, 2) = CtInv(1, 2)*tmp
  CtInv(1, 3) = CtInv(1, 3)*tmp
  CtInv(2, 1) = (Ct(2, 3)*Ct(3, 1) &
                 - Ct(2, 1)*Ct(3, 3))*tmp
  CtInv(2, 2) = (Ct(1, 1)*Ct(3, 3) &
                 - Ct(3, 1)*Ct(1, 3))*tmp
  CtInv(2, 3) = (Ct(2, 1)*Ct(1, 3) &
                 - Ct(1, 1)*Ct(2, 3))*tmp
  CtInv(3, 1) = (Ct(2, 1)*Ct(3, 2) &
                 - Ct(2, 2)*Ct(3, 1))*tmp
  CtInv(3, 2) = (Ct(3, 1)*Ct(1, 2) &
                 - Ct(1, 1)*Ct(3, 2))*tmp
  CtInv(3, 3) = (Ct(1, 1)*Ct(2, 2) &
                 - Ct(1, 2)*Ct(2, 1))*tmp

  trC = Ct(1, 1) + Ct(2, 2) + Ct(3, 3) ! Trace of the F^T Ftensor

  ! Compute Second P-K tensor
  ! for Neo-Hookean Mat w/ slight compressibility

  Stens = 0.0d0

  Stens(:, :) = 0.5d0*lambda*(Jdet*Jdet - 1.0d0)*CtInv(:, :) + &
                mu*Jmod*(Itens(:, :) - pt33*trC*CtInv(:, :))

  ! Compute Nonlinear Tangent
  ! Stiffness Tensor C_{IJKL}

  Ctens = 0.0d0

  fact1 = lambda*Jdet*Jdet + 2d+0*pt33*pt33*mu*Jmod*trC
  fact2 = 2d+0*pt33*mu*Jmod
  fact3 = 2d+0*pt33*mu*Jmod*trC - lambda*(Jdet*Jdet - 1.0d0)
  fact3 = fact3*0.5d0

  do i = 1, NSD
    do j = 1, NSD
      do k = 1, NSD
        do l = 1, NSD
          Ctemp(i, j, k, l) = fact1*CtInv(i, j)*CtInv(k, l) - &
                              fact2*(Itens(i, j)*CtInv(k, l) + &
                                     CtInv(i, j)*Itens(k, l)) + &
                              fact3*(CtInv(i, k)*CtInv(j, l) + &
                                     CtInv(i, l)*CtInv(j, k))
        end do
      end do
    end do
  end do

  ! Build C_{ijJI} = F_{iK}C_{KIJL}F_{jL}
  do i = 1, NSD
    do j = 1, NSD
      do k = 1, NSD
        do l = 1, NSD
          Ctens(i, j, l, k) = Ftens(i, 1)*Ftens(j, 1)*Ctemp(1, k, l, 1) + &
                              Ftens(i, 1)*Ftens(j, 2)*Ctemp(1, k, l, 2) + &
                              Ftens(i, 1)*Ftens(j, 3)*Ctemp(1, k, l, 3) + &
                              Ftens(i, 2)*Ftens(j, 1)*Ctemp(2, k, l, 1) + &
                              Ftens(i, 2)*Ftens(j, 2)*Ctemp(2, k, l, 2) + &
                              Ftens(i, 2)*Ftens(j, 3)*Ctemp(2, k, l, 3) + &
                              Ftens(i, 3)*Ftens(j, 1)*Ctemp(3, k, l, 1) + &
                              Ftens(i, 3)*Ftens(j, 2)*Ctemp(3, k, l, 2) + &
                              Ftens(i, 3)*Ftens(j, 3)*Ctemp(3, k, l, 3)
        end do
      end do
    end do
  end do

end subroutine e3tensors_nonlin

