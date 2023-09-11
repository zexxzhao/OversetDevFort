!======================================================================
!
!======================================================================
subroutine eval_shape_nurbs(nshl, gp, ni, nj, nk, patch, xl, wl, &
                            shlu, shgradgu, shhessg, dxidx, &
                            Gij, Ginv, hess_flag)

  use commonvars
  use class_def
  implicit none

  type(NURBSpatch), intent(in) :: patch

  integer, intent(in) :: nshl, ni, nj, nk, hess_flag
  real(8), intent(in) :: gp(NSD), xl(NSHL, NSD), wl(NSHL)

  real(8), intent(out) :: shlu(NSHL), shgradgu(NSHL, NSD), &
                          shhessg(NSHL, NSD, NSD), dxidx(NSD, NSD), &
                          Gij(NSD, NSD), Ginv(NSD, NSD)

  integer :: Pu, Qu, Ru, i, j, k, aa, nders
  real(8) :: u_hat, v_hat, w_hat, du, dv, dw, da

  real(8) :: shgradlu(NSHL, NSD), shhessl(NSHL, 6), tempshl(NSHL), &
             tempshgradl(NSHL, NSD), tempshhessl(NSHL, 6)

  real(8) :: dxdxi(NSD, NSD), dxdxixj(NSD, 6), locLHS(6, 6)

  real(8) :: Nnu(2 + hess_flag, patch%P + 1), &
             Mmu(2 + hess_flag, patch%Q + 1), &
             Oou(2 + hess_flag, patch%R + 1)    ! Mu = Mmu, since mu is viscosity

  real(8) :: u, v, w, denom_sum, derv_sum_U, derv_sum_V, &
             derv_sum_W, derv_sum_UU, &
             derv_sum_UV, derv_sum_UW, derv_sum_VV, derv_sum_VW, &
             derv_sum_WW, tmp

  u_hat = gp(1)
  v_hat = gp(2)
  w_hat = gp(3)

  Pu = patch%P
  Qu = patch%Q
  Ru = patch%R

  ! Get u and v coordinates of integration point
  u = ((patch%U_KNOT(ni + 1) - patch%U_KNOT(ni))*u_hat + &
       patch%U_KNOT(ni + 1) + patch%U_KNOT(ni))/2.0d0
  v = ((patch%V_KNOT(nj + 1) - patch%V_KNOT(nj))*v_hat + &
       patch%V_KNOT(nj + 1) + patch%V_KNOT(nj))/2.0d0
  w = ((patch%W_KNOT(nk + 1) - patch%W_KNOT(nk))*w_hat + &
       patch%W_KNOT(nk + 1) + patch%W_KNOT(nk))/2.0d0

  ! Get knot span sizes
  du = patch%U_KNOT(ni + 1) - patch%U_KNOT(ni)
  dv = patch%V_KNOT(nj + 1) - patch%V_KNOT(nj)
  dw = patch%W_KNOT(nk + 1) - patch%W_KNOT(nk)
  da = du*dv*dw/8.0d0

  ! Evaluate 1D shape functions and derivatives each direction
  if (hess_flag == 0) then
    nders = 1       ! number of derivatives to take
  else
    nders = 2
  end if

  ! calculate in u,v,w dir.
  call dersbasisfuns(ni, Pu, patch%MCP, u, nders, patch%U_KNOT, Nnu)
  call dersbasisfuns(nj, Qu, patch%NCP, v, nders, patch%V_KNOT, Mmu)
  call dersbasisfuns(nk, Ru, patch%OCP, w, nders, patch%W_KNOT, Oou)

  ! Form basis functions and derivatives dR/du and dR/dv
  aa = 0
  denom_sum = 0.0d0
  derv_sum_U = 0.0d0
  derv_sum_V = 0.0d0
  derv_sum_W = 0.0d0
  shgradlu = 0.0d0
  shhessl = 0.0d0
  shhessg = 0.0d0
  tempshl = 0.0d0
  tempshgradl = 0.0d0
  tempshhessl = 0.0d0

  do k = 0, Ru
    do j = 0, Qu
      do i = 0, Pu
        aa = aa + 1

        ! basis functions
        shlu(aa) = Nnu(1, Pu + 1 - i)*Mmu(1, Qu + 1 - j)*Oou(1, Ru + 1 - k)*wl(aa)
        denom_sum = denom_sum + shlu(aa)

        ! derivatives (u, v and w)
        shgradlu(aa, 1) = Nnu(2, Pu + 1 - i)*Mmu(1, Qu + 1 - j)*Oou(1, Ru + 1 - k)* &
                          wl(aa)
        derv_sum_U = derv_sum_U + shgradlu(aa, 1)

        shgradlu(aa, 2) = Nnu(1, Pu + 1 - i)*Mmu(2, Qu + 1 - j)*Oou(1, Ru + 1 - k)* &
                          wl(aa)
        derv_sum_V = derv_sum_V + shgradlu(aa, 2)

        shgradlu(aa, 3) = Nnu(1, Pu + 1 - i)*Mmu(1, Qu + 1 - j)*Oou(2, Ru + 1 - k)* &
                          wl(aa)
        derv_sum_W = derv_sum_W + shgradlu(aa, 3)

      end do
    end do
  end do

  ! Divide through by denominator
  tempshl = shlu
  tempshgradl = shgradlu

  do i = 1, NSHL
    shgradlu(i, 1) = shgradlu(i, 1)/denom_sum - &
                     (shlu(i)*derv_sum_U)/(denom_sum**2)
    shgradlu(i, 2) = shgradlu(i, 2)/denom_sum - &
                     (shlu(i)*derv_sum_V)/(denom_sum**2)
    shgradlu(i, 3) = shgradlu(i, 3)/denom_sum - &
                     (shlu(i)*derv_sum_W)/(denom_sum**2)
    shlu(i) = shlu(i)/denom_sum
  end do

  ! Now calculate gradients.
  ! calculate dx/dxi
  dxdxi = 0.0d0
  do i = 1, NSHL
    do j = 1, NSD
      dxdxi(j, :) = dxdxi(j, :) + xl(i, j)*shgradlu(i, :)
    end do
  end do

  call get_inverse_3x3(dxdxi, dxidx, DetJ)

  do i = 1, NSHL
    do j = 1, NSD
      shgradgu(i, j) = sum(shgradlu(i, :)*dxidx(:, j))
    end do
  end do

  !----------------------------------------------
  ! 2nd derivatives
  !----------------------------------------------
  if (hess_flag == 1) then

    aa = 0
    derv_sum_UU = 0.0d0
    derv_sum_UV = 0.0d0
    derv_sum_UW = 0.0d0
    derv_sum_VV = 0.0d0
    derv_sum_VW = 0.0d0
    derv_sum_WW = 0.0d0
    do k = 0, Ru
      do j = 0, Qu
        do i = 0, Pu
          aa = aa + 1

          ! 2nd derivatives
          ! ,uu
          tempshhessl(aa, 1) = Nnu(3, Pu + 1 - i)*Mmu(1, Qu + 1 - j)* &
                               Oou(1, Ru + 1 - k)*wl(aa)
          derv_sum_UU = derv_sum_UU + tempshhessl(aa, 1)
          ! ,uv
          tempshhessl(aa, 2) = Nnu(2, Pu + 1 - i)*Mmu(2, Qu + 1 - j)* &
                               Oou(1, Ru + 1 - k)*wl(aa)
          derv_sum_UV = derv_sum_UV + tempshhessl(aa, 2)
          ! ,uw
          tempshhessl(aa, 3) = Nnu(2, Pu + 1 - i)*Mmu(1, Qu + 1 - j)* &
                               Oou(2, Ru + 1 - k)*wl(aa)
          derv_sum_UW = derv_sum_UW + tempshhessl(aa, 3)
          ! ,vv
          tempshhessl(aa, 4) = Nnu(1, Pu + 1 - i)*Mmu(3, Qu + 1 - j)* &
                               Oou(1, Ru + 1 - k)*wl(aa)
          derv_sum_VV = derv_sum_VV + tempshhessl(aa, 4)
          ! ,vw
          tempshhessl(aa, 5) = Nnu(1, Pu + 1 - i)*Mmu(2, Qu + 1 - j)* &
                               Oou(2, Ru + 1 - k)*wl(aa)
          derv_sum_VW = derv_sum_VW + tempshhessl(aa, 5)
          ! ,ww
          tempshhessl(aa, 6) = Nnu(1, Pu + 1 - i)*Mmu(1, Qu + 1 - j)* &
                               Oou(3, Ru + 1 - k)*wl(aa)
          derv_sum_WW = derv_sum_WW + tempshhessl(aa, 6)
        end do
      end do
    end do

    ! local Hessian
    do i = 1, NSHL
      shhessl(i, 1) = tempshhessl(i, 1)/denom_sum - &
                      tempshgradl(i, 1)*derv_sum_U/(denom_sum**2) - &
                      ((tempshgradl(i, 1)*derv_sum_U + tempshl(i)*derv_sum_UU)/ &
                       (denom_sum**2) - 2d+0*tempshl(i)*derv_sum_U*derv_sum_U/ &
                       (denom_sum**3))
      shhessl(i, 2) = tempshhessl(i, 2)/denom_sum - &
                      tempshgradl(i, 1)*derv_sum_V/(denom_sum**2) - &
                      ((tempshgradl(i, 2)*derv_sum_U + tempshl(i)*derv_sum_UV)/ &
                       (denom_sum**2) - 2d+0*tempshl(i)*derv_sum_U*derv_sum_V/ &
                       (denom_sum**3))
      shhessl(i, 3) = tempshhessl(i, 3)/denom_sum - &
                      tempshgradl(i, 1)*derv_sum_W/(denom_sum**2) - &
                      ((tempshgradl(i, 3)*derv_sum_U + tempshl(i)*derv_sum_UW)/ &
                       (denom_sum**2) - 2d+0*tempshl(i)*derv_sum_U*derv_sum_W/ &
                       (denom_sum**3))
      shhessl(i, 4) = tempshhessl(i, 4)/denom_sum - &
                      tempshgradl(i, 2)*derv_sum_V/(denom_sum**2) - &
                      ((tempshgradl(i, 2)*derv_sum_V + tempshl(i)*derv_sum_VV)/ &
                       (denom_sum**2) - 2d+0*tempshl(i)*derv_sum_V*derv_sum_V/ &
                       (denom_sum**3))
      shhessl(i, 5) = tempshhessl(i, 5)/denom_sum - &
                      tempshgradl(i, 2)*derv_sum_W/(denom_sum**2) - &
                      ((tempshgradl(i, 3)*derv_sum_V + tempshl(i)*derv_sum_VW)/ &
                       (denom_sum**2) - 2d+0*tempshl(i)*derv_sum_V*derv_sum_W/ &
                       (denom_sum**3))
      shhessl(i, 6) = tempshhessl(i, 6)/denom_sum - &
                      tempshgradl(i, 3)*derv_sum_W/(denom_sum**2) - &
                      ((tempshgradl(i, 3)*derv_sum_W + tempshl(i)*derv_sum_WW)/ &
                       (denom_sum**2) - 2d+0*tempshl(i)*derv_sum_W*derv_sum_W/ &
                       (denom_sum**3))
    end do

    ! global Hessian
    ! Second derivatives of the geometrical map

    dxdxixj = 0.0d0
    aa = 0
    do i = 1, NSHL
      do j = 1, NSD
        dxdxixj(j, :) = dxdxixj(j, :) + xl(i, j)*shhessl(i, :)
      end do
    end do

    ! RHS of the matrix equation for the second derivatives of bases.
    ! Reuse local hess. array

    shhessl(:, 1) = shhessl(:, 1) - shgradgu(:, 1)*dxdxixj(1, 1) - &
                    shgradgu(:, 2)*dxdxixj(2, 1) - &
                    shgradgu(:, 3)*dxdxixj(3, 1)

    shhessl(:, 2) = shhessl(:, 2) - shgradgu(:, 1)*dxdxixj(1, 2) - &
                    shgradgu(:, 2)*dxdxixj(2, 2) - &
                    shgradgu(:, 3)*dxdxixj(3, 2)

    shhessl(:, 3) = shhessl(:, 3) - shgradgu(:, 1)*dxdxixj(1, 3) - &
                    shgradgu(:, 2)*dxdxixj(2, 3) - &
                    shgradgu(:, 3)*dxdxixj(3, 3)

    shhessl(:, 4) = shhessl(:, 4) - shgradgu(:, 1)*dxdxixj(1, 4) - &
                    shgradgu(:, 2)*dxdxixj(2, 4) - &
                    shgradgu(:, 3)*dxdxixj(3, 4)

    shhessl(:, 5) = shhessl(:, 5) - shgradgu(:, 1)*dxdxixj(1, 5) - &
                    shgradgu(:, 2)*dxdxixj(2, 5) - &
                    shgradgu(:, 3)*dxdxixj(3, 5)

    shhessl(:, 6) = shhessl(:, 6) - shgradgu(:, 1)*dxdxixj(1, 6) - &
                    shgradgu(:, 2)*dxdxixj(2, 6) - &
                    shgradgu(:, 3)*dxdxixj(3, 6)

    ! LHS (6x6, same for every basis function)
    locLHS(1, 1) = dxdxi(1, 1)*dxdxi(1, 1)
    locLHS(1, 2) = 2d+0*dxdxi(1, 1)*dxdxi(2, 1)
    locLHS(1, 3) = 2d+0*dxdxi(1, 1)*dxdxi(3, 1)
    locLHS(1, 4) = dxdxi(2, 1)*dxdxi(2, 1)
    locLHS(1, 5) = 2d+0*dxdxi(2, 1)*dxdxi(3, 1)
    locLHS(1, 6) = dxdxi(3, 1)*dxdxi(3, 1)

    locLHS(2, 1) = dxdxi(1, 1)*dxdxi(1, 2)
    locLHS(2, 2) = dxdxi(1, 1)*dxdxi(2, 2) + dxdxi(1, 2)*dxdxi(2, 1)
    locLHS(2, 3) = dxdxi(1, 1)*dxdxi(3, 2) + dxdxi(1, 2)*dxdxi(3, 1)
    locLHS(2, 4) = dxdxi(2, 1)*dxdxi(2, 2)
    locLHS(2, 5) = dxdxi(2, 1)*dxdxi(3, 2) + dxdxi(2, 2)*dxdxi(3, 1)
    locLHS(2, 6) = dxdxi(3, 1)*dxdxi(3, 2)

    locLHS(3, 1) = dxdxi(1, 1)*dxdxi(1, 3)
    locLHS(3, 2) = dxdxi(1, 1)*dxdxi(2, 3) + dxdxi(1, 3)*dxdxi(2, 1)
    locLHS(3, 3) = dxdxi(1, 1)*dxdxi(3, 3) + dxdxi(1, 3)*dxdxi(3, 1)
    locLHS(3, 4) = dxdxi(2, 1)*dxdxi(2, 3)
    locLHS(3, 5) = dxdxi(2, 1)*dxdxi(3, 3) + dxdxi(2, 3)*dxdxi(3, 1)
    locLHS(3, 6) = dxdxi(3, 1)*dxdxi(3, 3)

    locLHS(4, 1) = dxdxi(1, 2)*dxdxi(1, 2)
    locLHS(4, 2) = 2d+0*dxdxi(1, 2)*dxdxi(2, 2)
    locLHS(4, 3) = 2d+0*dxdxi(1, 2)*dxdxi(3, 2)
    locLHS(4, 4) = dxdxi(2, 2)*dxdxi(2, 2)
    locLHS(4, 5) = 2d+0*dxdxi(2, 2)*dxdxi(3, 2)
    locLHS(4, 6) = dxdxi(3, 2)*dxdxi(3, 2)

    locLHS(5, 1) = dxdxi(1, 2)*dxdxi(1, 3)
    locLHS(5, 2) = dxdxi(1, 2)*dxdxi(2, 3) + dxdxi(1, 3)*dxdxi(2, 2)
    locLHS(5, 3) = dxdxi(1, 2)*dxdxi(3, 3) + dxdxi(1, 3)*dxdxi(3, 2)
    locLHS(5, 4) = dxdxi(2, 2)*dxdxi(2, 3)
    locLHS(5, 5) = dxdxi(2, 2)*dxdxi(3, 3) + dxdxi(2, 3)*dxdxi(3, 2)
    locLHS(5, 6) = dxdxi(3, 2)*dxdxi(3, 3)

    locLHS(6, 1) = dxdxi(1, 3)*dxdxi(1, 3)
    locLHS(6, 2) = 2d+0*dxdxi(1, 3)*dxdxi(2, 3)
    locLHS(6, 3) = 2d+0*dxdxi(1, 3)*dxdxi(3, 3)
    locLHS(6, 4) = dxdxi(2, 3)*dxdxi(2, 3)
    locLHS(6, 5) = 2d+0*dxdxi(2, 3)*dxdxi(3, 3)
    locLHS(6, 6) = dxdxi(3, 3)*dxdxi(3, 3)

    ! (6x6) - Gaussian elimination
    do k = 1, 6
      do i = k + 1, 6
        tmp = locLHS(i, k)/locLHS(k, k)
        do j = k + 1, 6
          locLHS(i, j) = locLHS(i, j) - tmp*locLHS(k, j)
        end do
        shhessl(:, i) = shhessl(:, i) - tmp*shhessl(:, k)
      end do
    end do

    do i = 6, 1, -1
      do j = i + 1, 6
        shhessl(:, i) = shhessl(:, i) - locLHS(i, j)*shhessl(:, j)
      end do
      shhessl(:, i) = shhessl(:, i)/locLHS(i, i)
    end do

    ! Assign to global hessian of basis functions
    shhessg(:, 1, 1) = shhessl(:, 1)
    shhessg(:, 1, 2) = shhessl(:, 2)
    shhessg(:, 1, 3) = shhessl(:, 3)

    shhessg(:, 2, 1) = shhessl(:, 2)
    shhessg(:, 2, 2) = shhessl(:, 4)
    shhessg(:, 2, 3) = shhessl(:, 5)

    shhessg(:, 3, 1) = shhessl(:, 3)
    shhessg(:, 3, 2) = shhessl(:, 5)
    shhessg(:, 3, 3) = shhessl(:, 6)
  end if

  !-- end 2nd derivatives -----------------------

  dxidx(1, :) = dxidx(1, :)*2.0d0/du
  dxidx(2, :) = dxidx(2, :)*2.0d0/dv
  dxidx(3, :) = dxidx(3, :)*2.0d0/dw

  dxdxi(:, 1) = dxdxi(:, 1)*du/2.0d0
  dxdxi(:, 2) = dxdxi(:, 2)*dv/2.0d0
  dxdxi(:, 3) = dxdxi(:, 3)*dw/2.0d0

  do j = 1, NSD
    do i = 1, NSD
      Gij(i, j) = sum(dxidx(:, i)*dxidx(:, j))
      Ginv(i, j) = sum(dxdxi(i, :)*dxdxi(j, :))
    end do
  end do

  DetJ = abs(DetJ)*da

end subroutine eval_shape_nurbs

!======================================================================
! "shb" will be the shape function array while "shbg" will hold the
! gradients of the shape functions
!======================================================================
subroutine eval_faceshape_nurbs(nshl, gp, faceor, ni, nj, nk, patch, &
                                xl, wl, shb, shbg, dxidx, &
                                Gij, Ginv, nor)
  use class_def
  use commonvars
  implicit none

  type(NURBSpatch), intent(in) :: patch

  real(8), intent(in) :: gp(NSD - 1), xl(NSHL, NSD), wl(NSHL)
  integer, intent(in) :: nshl, faceor, ni, nj, nk

  real(8), intent(out) :: shb(NSHL), shbg(NSHL, NSD), dxidx(NSD, NSD), &
                          nor(NSD), Gij(NSD, NSD), Ginv(NSD, NSD)

  ! Local Variables
  integer :: Pu, Qu, Ru, aa, i, j, k
  real(8) :: u_hat, v_hat, w_hat, tmp, nor_0(NSD), cofF(NSD, NSD), &
             du, dv, dw, da
  real(8) :: dxdxi(NSD, NSD), shbl(NSHL, NSD), ee, ff, gg
  real(8) :: N(2, patch%P + 1), M(2, patch%Q + 1), O(2, patch%R + 1)
  real(8) :: u, v, w, denom_sum, derv_sum_U, derv_sum_V, derv_sum_W

  ! Get NURBS coordinates for local node
  Pu = patch%P
  Qu = patch%Q
  Ru = patch%R

  du = patch%U_KNOT(ni + 1) - patch%U_KNOT(ni)
  dv = patch%V_KNOT(nj + 1) - patch%V_KNOT(nj)
  dw = patch%W_KNOT(nk + 1) - patch%W_KNOT(nk)

  ! The details are dependent of face orientation.
  nor_0 = 0.0d0

  if (faceor == 1) then
    u_hat = gp(1)
    v_hat = gp(2)
    w_hat = -1.0d0     ! -1.0 ensure w will be at beginning knot vector
    nor_0(3) = -1.0d0
    da = du*dv/4.0d0

  else if (faceor == 2) then
    u_hat = gp(1)
    v_hat = -1.0d0
    w_hat = gp(2)

    nor_0(2) = -1.0d0
    da = du*dw/4.0d0

  else if (faceor == 3) then
    u_hat = 1.0d0      ! 1.0 ensures u will be at end of knot vector
    v_hat = gp(1)
    w_hat = gp(2)

    nor_0(1) = 1.0d0
    da = dv*dw/4.0d0

  else if (faceor == 4) then
    u_hat = gp(1)
    v_hat = 1.0d0
    w_hat = gp(2)

    nor_0(2) = 1.0d0
    da = du*dw/4.0d0

  else if (faceor == 5) then
    u_hat = -1.0d0
    v_hat = gp(1)
    w_hat = gp(2)

    nor_0(1) = -1.0d0
    da = dv*dw/4.0d0

  else if (faceor == 6) then
    u_hat = gp(1)
    v_hat = gp(2)
    w_hat = 1.0d0

    nor_0(3) = 1.0d0
    da = du*dv/4.0d0

  else
    write (*, *) "ERROR: wrong face orientation"
    stop

  end if

  ! Find integration point in parameter space
  u = ((patch%U_KNOT(ni + 1) - patch%U_KNOT(ni))*u_hat + &
       patch%U_KNOT(ni + 1) + patch%U_KNOT(ni))/2d+0

  v = ((patch%V_KNOT(nj + 1) - patch%V_KNOT(nj))*v_hat + &
       patch%V_KNOT(nj + 1) + patch%V_KNOT(nj))/2d+0

  w = ((patch%W_KNOT(nk + 1) - patch%W_KNOT(nk))*w_hat + &
       patch%W_KNOT(nk + 1) + patch%W_KNOT(nk))/2d+0

  ! Evaluate 1D shape functions for velocity
  call dersbasisfuns(ni, Pu, patch%MCP, u, 1, patch%U_KNOT, N)
  call dersbasisfuns(nj, Qu, patch%NCP, v, 1, patch%V_KNOT, M)
  call dersbasisfuns(nk, Ru, patch%OCP, w, 1, patch%W_KNOT, O)

  ! Form basis functions for velocity
  aa = 0
  denom_sum = 0.0d0

  shb = 0.0d0
  derv_sum_U = 0.0d0
  derv_sum_V = 0.0d0
  derv_sum_W = 0.0d0
  shbl = 0.0d0
  shbg = 0.0d0

  do k = 0, Ru
    do j = 0, Qu
      do i = 0, Pu
        aa = aa + 1

        ! basis functions
        shb(aa) = N(1, Pu + 1 - i)*M(1, Qu + 1 - j)*O(1, Ru + 1 - k)*wl(aa)
        denom_sum = denom_sum + shb(aa)

        ! derivatives (u, v and w)
        shbl(aa, 1) = N(2, Pu + 1 - i)*M(1, Qu + 1 - j)*O(1, Ru + 1 - k)*wl(aa)
        derv_sum_U = derv_sum_U + shbl(aa, 1)

        shbl(aa, 2) = N(1, Pu + 1 - i)*M(2, Qu + 1 - j)*O(1, Ru + 1 - k)*wl(aa)
        derv_sum_V = derv_sum_V + shbl(aa, 2)

        shbl(aa, 3) = N(1, Pu + 1 - i)*M(1, Qu + 1 - j)*O(2, Ru + 1 - k)*wl(aa)
        derv_sum_W = derv_sum_W + shbl(aa, 3)

      end do
    end do
  end do

  ! Divide through by denominator
  do i = 1, NSHL
    shbl(i, 1) = shbl(i, 1)/denom_sum - &
                 (shb(i)*derv_sum_U)/(denom_sum**2)
    shbl(i, 2) = shbl(i, 2)/denom_sum - &
                 (shb(i)*derv_sum_V)/(denom_sum**2)
    shbl(i, 3) = shbl(i, 3)/denom_sum - &
                 (shb(i)*derv_sum_W)/(denom_sum**2)
    shb(i) = shb(i)/denom_sum
  end do

  ! Now we calculate the face Jacobian
  dxdxi(:, :) = 0.0d0
  do i = 1, NSHL
    do j = 1, NSD
      dxdxi(j, :) = dxdxi(j, :) + xl(i, j)*shbl(i, :)
    end do
  end do

  ! Calculate Cofactor

  cofF(1, 1) = dxdxi(2, 2)*dxdxi(3, 3) - dxdxi(2, 3)*dxdxi(3, 2)
  cofF(1, 2) = -(dxdxi(2, 1)*dxdxi(3, 3) - dxdxi(2, 3)*dxdxi(3, 1))
  cofF(1, 3) = dxdxi(2, 1)*dxdxi(3, 2) - dxdxi(2, 2)*dxdxi(3, 1)

  cofF(2, 1) = -(dxdxi(1, 2)*dxdxi(3, 3) - dxdxi(1, 3)*dxdxi(3, 2))
  cofF(2, 2) = dxdxi(1, 1)*dxdxi(3, 3) - dxdxi(1, 3)*dxdxi(3, 1)
  cofF(2, 3) = -(dxdxi(1, 1)*dxdxi(3, 2) - dxdxi(1, 2)*dxdxi(3, 1))

  cofF(3, 1) = dxdxi(1, 2)*dxdxi(2, 3) - dxdxi(1, 3)*dxdxi(2, 2)
  cofF(3, 2) = -(dxdxi(1, 1)*dxdxi(2, 3) - dxdxi(1, 3)*dxdxi(2, 1))
  cofF(3, 3) = dxdxi(1, 1)*dxdxi(2, 2) - dxdxi(1, 2)*dxdxi(2, 1)

  ! Jacobian of Face Mapping
  if ((faceor == 1) .or. (faceor == 6)) then
    ee = dxdxi(1, 1)**2 + dxdxi(2, 1)**2 + dxdxi(3, 1)**2
    ff = dxdxi(1, 1)*dxdxi(1, 2) + dxdxi(2, 1)*dxdxi(2, 2) + &
         dxdxi(3, 1)*dxdxi(3, 2)
    gg = dxdxi(1, 2)**2 + dxdxi(2, 2)**2 + dxdxi(3, 2)**2

  else if ((faceor == 2) .or. (faceor == 4)) then
    ee = dxdxi(1, 1)**2 + dxdxi(2, 1)**2 + dxdxi(3, 1)**2
    ff = dxdxi(1, 1)*dxdxi(1, 3) + dxdxi(2, 1)*dxdxi(2, 3) + &
         dxdxi(3, 1)*dxdxi(3, 3)
    gg = dxdxi(1, 3)**2 + dxdxi(2, 3)**2 + dxdxi(3, 3)**2

  else if ((faceor == 3) .or. (faceor == 5)) then
    ee = dxdxi(1, 2)**2 + dxdxi(2, 2)**2 + dxdxi(3, 2)**2
    ff = dxdxi(1, 2)*dxdxi(1, 3) + dxdxi(2, 2)*dxdxi(2, 3) + &
         dxdxi(3, 2)*dxdxi(3, 3)
    gg = dxdxi(1, 3)**2 + dxdxi(2, 3)**2 + dxdxi(3, 3)**2
  end if

  DetJb = sqrt(ee*gg - ff**2) ! Jacobian of face mapping

  !------------------------------------------------------------------

  ! Computation of normal
  if (faceor == 1) then
    nor(1) = dxdxi(2, 2)*dxdxi(3, 1) - dxdxi(3, 2)*dxdxi(2, 1)
    nor(2) = dxdxi(3, 2)*dxdxi(1, 1) - dxdxi(1, 2)*dxdxi(3, 1)
    nor(3) = dxdxi(1, 2)*dxdxi(2, 1) - dxdxi(2, 2)*dxdxi(1, 1)
  else if (faceor == 2) then
    nor(1) = dxdxi(2, 1)*dxdxi(3, 3) - dxdxi(3, 1)*dxdxi(2, 3)
    nor(2) = dxdxi(3, 1)*dxdxi(1, 3) - dxdxi(1, 1)*dxdxi(3, 3)
    nor(3) = dxdxi(1, 1)*dxdxi(2, 3) - dxdxi(2, 1)*dxdxi(1, 3)
  else if (faceor == 3) then
    nor(1) = dxdxi(2, 2)*dxdxi(3, 3) - dxdxi(3, 2)*dxdxi(2, 3)
    nor(2) = dxdxi(3, 2)*dxdxi(1, 3) - dxdxi(1, 2)*dxdxi(3, 3)
    nor(3) = dxdxi(1, 2)*dxdxi(2, 3) - dxdxi(2, 2)*dxdxi(1, 3)
  else if (faceor == 4) then
    nor(1) = -dxdxi(2, 1)*dxdxi(3, 3) + dxdxi(3, 1)*dxdxi(2, 3)
    nor(2) = -dxdxi(3, 1)*dxdxi(1, 3) + dxdxi(1, 1)*dxdxi(3, 3)
    nor(3) = -dxdxi(1, 1)*dxdxi(2, 3) + dxdxi(2, 1)*dxdxi(1, 3)
  else if (faceor == 5) then
    nor(1) = -dxdxi(2, 2)*dxdxi(3, 3) + dxdxi(3, 2)*dxdxi(2, 3)
    nor(2) = -dxdxi(3, 2)*dxdxi(1, 3) + dxdxi(1, 2)*dxdxi(3, 3)
    nor(3) = -dxdxi(1, 2)*dxdxi(2, 3) + dxdxi(2, 2)*dxdxi(1, 3)
  else if (faceor == 6) then
    nor(1) = -dxdxi(2, 2)*dxdxi(3, 1) + dxdxi(3, 2)*dxdxi(2, 1)
    nor(2) = -dxdxi(3, 2)*dxdxi(1, 1) + dxdxi(1, 2)*dxdxi(3, 1)
    nor(3) = -dxdxi(1, 2)*dxdxi(2, 1) + dxdxi(2, 2)*dxdxi(1, 1)
  end if

  tmp = sqrt(sum(nor*nor))
  nor = nor/tmp

  ! compute the inverse of deformation gradient
  ! 1/ Det

  tmp = 1.0d0/(dxdxi(1, 1)*dxdxi(2, 2)*dxdxi(3, 3) - &
               dxdxi(1, 1)*dxdxi(2, 3)*dxdxi(3, 2) - &
               dxdxi(2, 1)*dxdxi(1, 2)*dxdxi(3, 3) + &
               dxdxi(2, 1)*dxdxi(1, 3)*dxdxi(3, 2) + &
               dxdxi(3, 1)*dxdxi(1, 2)*dxdxi(2, 3) - &
               dxdxi(3, 1)*dxdxi(1, 3)*dxdxi(2, 2))

  ! F^(-1) = 1/Det * (cofF)^T
  do i = 1, NSD
    do j = 1, NSD
      dxidx(i, j) = cofF(j, i)
    end do
  end do

  dxidx(:, :) = dxidx(:, :)*tmp

  do i = 1, NSHL
    shbg(i, 1) = shbl(i, 1)*dxidx(1, 1) + &
                 shbl(i, 2)*dxidx(2, 1) + &
                 shbl(i, 3)*dxidx(3, 1)
    shbg(i, 2) = shbl(i, 1)*dxidx(1, 2) + &
                 shbl(i, 2)*dxidx(2, 2) + &
                 shbl(i, 3)*dxidx(3, 2)
    shbg(i, 3) = shbl(i, 1)*dxidx(1, 3) + &
                 shbl(i, 2)*dxidx(2, 3) + &
                 shbl(i, 3)*dxidx(3, 3)
  end do

  dxidx(1, :) = dxidx(1, :)*2.0d0/du
  dxidx(2, :) = dxidx(2, :)*2.0d0/dv
  dxidx(3, :) = dxidx(3, :)*2.0d0/dw

  dxdxi(:, 1) = dxdxi(:, 1)*du/2.0d0
  dxdxi(:, 2) = dxdxi(:, 2)*dv/2.0d0
  dxdxi(:, 3) = dxdxi(:, 3)*dw/2.0d0

  do j = 1, NSD
    do i = 1, NSD
      Gij(i, j) = sum(dxidx(:, i)*dxidx(:, j))
      Ginv(i, j) = sum(dxdxi(i, :)*dxdxi(j, :))
    end do
  end do

  DetJb = DetJb*da

end subroutine eval_faceshape_nurbs

!======================================================================
! This subroutine evaluates the basis functions and first derivatives
! functions at a given parameter value u.
!
! Algorithm from Piegl, Les. "The NURBS Book". Springer-Verlag:
! Berlin 1995; pp. 72-73.
!======================================================================
subroutine dersbasisfuns(i, pl, ml, u, nders, u_knotl, ders)

  implicit none

  ! knot span, degree of curve, number of control points, counters
  integer :: i, pl, ml, j, r, k, j1, j2, s1, s2, rk, pk, nders

  ! parameter value, vector of knots, derivative matrix
  real(8) :: u, u_knotl(pl + ml + 1), ders(nders + 1, pl + 1), ndu(pl + 1, pl + 1), d
  real(8) :: left(pl + 1), right(pl + 1), saved, temp, a(2, pl + 1)

  ndu(1, 1) = 1.0d0
  do j = 1, pl
    left(j + 1) = u - u_knotl(i + 1 - j)
    right(j + 1) = u_knotl(i + j) - u
    saved = 0.0d0
    do r = 0, j - 1
      ndu(j + 1, r + 1) = right(r + 2) + left(j - r + 1)
      temp = ndu(r + 1, j)/ndu(j + 1, r + 1)
      ndu(r + 1, j + 1) = saved + right(r + 2)*temp
      saved = left(j - r + 1)*temp
    end do
    ndu(j + 1, j + 1) = saved
  end do

  ! load basis functions
  do j = 0, pl
    ders(1, j + 1) = ndu(j + 1, pl + 1)
  end do

  ! compute derivatives
  do r = 0, pl ! loop over function index
    s1 = 0
    s2 = 1         ! alternate rows in array a
    a(1, 1) = 1.0d0

    ! loop to compute kth derivative
    do k = 1, nders
      d = 0.0d0
      rk = r - k
      pk = pl - k
      if (r >= k) then
        a(s2 + 1, 1) = a(s1 + 1, 1)/ndu(pk + 2, rk + 1)
        d = a(s2 + 1, 1)*ndu(rk + 1, pk + 1)
      end if
      if (rk >= -1) then
        j1 = 1
      else
        j1 = -rk
      end if
      if ((r - 1) <= pk) then
        j2 = k - 1
      else
        j2 = pl - r
      end if
      do j = j1, j2
        a(s2 + 1, j + 1) = (a(s1 + 1, j + 1) - a(s1 + 1, j))/ndu(pk + 2, rk + j + 1)
        d = d + a(s2 + 1, j + 1)*ndu(rk + j + 1, pk + 1)
      end do
      if (r <= pk) then
        a(s2 + 1, k + 1) = -a(s1 + 1, k)/ndu(pk + 2, r + 1)
        d = d + a(s2 + 1, k + 1)*ndu(r + 1, pk + 1)
      end if
      ders(k + 1, r + 1) = d
      j = s1
      s1 = s2
      s2 = j      ! switch rows
    end do
  end do

  ! Multiply through by the correct factors
  r = pl
  do k = 1, nders
    do j = 0, pl
      ders(k + 1, j + 1) = ders(k + 1, j + 1)*r
    end do
    r = r*(pl - k)
  end do

end subroutine dersbasisfuns
