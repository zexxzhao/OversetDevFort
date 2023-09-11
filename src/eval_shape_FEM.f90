!======================================================================
!
!======================================================================
subroutine eval_shape_fem(nshl, gp, xl, dl, shlu, shgradgu, dxidx, &
                          Gij, Ginv)
  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in)  :: nshl
  real(8), intent(in)  :: gp(NSD), dl(NSHL, NSD), xl(NSHL, NSD)
  real(8), intent(out) :: shlu(NSHL), shgradgu(NSHL, NSD), dxidx(NSD, NSD), &
                          Gij(NSD, NSD), Ginv(NSD, NSD)
  integer :: i, j
  real(8) :: dxdxi(NSD, NSD), shgradlu(NSHL, NSD), tmp

  ! Get basis functions and local ders
  if (nshl == 4) then
    call lintetshl(gp, NSHL, NSD, shlu, shgradlu)
  else if (nshl == 6) then
    call linprishl(gp, NSHL, NSD, shlu, shgradlu)
  else
    write (*, *) "Undefined nshl in eval_shape_fem"
    stop
  end if
!!!   call linhexshl(gp, NSHL, NSD, shlu, shgradlu)

  ! Now calculate gradients.

  ! Calculate dx/dxi
  dxdxi = 0.0d0
  do i = 1, NSD
    do j = 1, NSD
      dxdxi(i, j) = sum((xl(:, i) + dl(:, i))*shgradlu(:, j))
    end do
  end do

  call get_inverse_3x3(dxdxi, dxidx, DetJ)

  ! Global gradients
  do i = 1, NSHL
    do j = 1, NSD
      shgradgu(i, j) = sum(shgradlu(i, :)*dxidx(:, j))
    end do
  end do

  do j = 1, NSD
    do i = 1, NSD
      Gij(i, j) = sum(dxidx(:, i)*dxidx(:, j))
      Ginv(i, j) = sum(dxdxi(i, :)*dxdxi(j, :))
    end do
  end do

  if (DetJ < 0.0d0) then
    write (*, *) "Warning - negative determinant Jacobian !"
    write (*, *) "Element Type:", nshl
    write (*, *) "DetJ = ", DetJ
    write (*, *) "DxDXI = ", dxdxi
    write (*, *) "DxIDX = ", dxidx
  end if
end subroutine eval_shape_fem

!======================================================================
!
!======================================================================
subroutine eval_faceshape_tet(nshl, bgp, igp, for, xl, dl, shlu, &
                              shgradgu, dxidx, lnor)
  use aAdjKeep
  use commonvars
  use mpi

  implicit none

  integer, intent(in) :: nshl
  real(8), intent(in) :: bgp(2), igp(NSD)

  real(8) :: dl(NSHL, NSD), xl(NSHL, NSD)
  real(8) :: shlu(NSHL), shgradgu(NSHL, NSD), dxidx(NSD, NSD), lnor(NSD)
  integer :: for

  real(8) :: bshlu(NSHLbmax), bshgradlu(NSHLbmax, 2)
  real(8) :: bxdl(NSHLbmax, NSD)
  real(8) :: dlxdxi(NSD, 2)
  real(8) :: diff(NSD)
  real(8) :: Gij(NSD, NSD), Ginv(NSD, NSD)
  integer :: i, dd

  if (for == 1) then
    bxdl(1, :) = xl(2, :) + dl(2, :)
    bxdl(2, :) = xl(3, :) + dl(3, :)
    bxdl(3, :) = xl(1, :) + dl(1, :)
  else if (for == 2) then
    bxdl(1, :) = xl(4, :) + dl(4, :)
    bxdl(2, :) = xl(3, :) + dl(3, :)
    bxdl(3, :) = xl(2, :) + dl(2, :)
  else if (for == 3) then
    bxdl(1, :) = xl(1, :) + dl(1, :)
    bxdl(2, :) = xl(3, :) + dl(3, :)
    bxdl(3, :) = xl(4, :) + dl(4, :)
  else if (for == 4) then
    bxdl(1, :) = xl(2, :) + dl(2, :)
    bxdl(2, :) = xl(1, :) + dl(1, :)
    bxdl(3, :) = xl(4, :) + dl(4, :)
  end if

  ! Get basis functions and local ders
  call lintrishl(bgp, NSHLbmax, NSD, bshlu, bshgradlu)

  ! Now calculate gradients.
  ! Calculate dx/dxi
  dlxdxi = 0.0d0
  do i = 1, NSHLbmax
    do dd = 1, NSD
      dlxdxi(dd, :) = dlxdxi(dd, :) + bxdl(i, dd)*bshgradlu(i, :)
    end do
  end do

  lnor = 0.0d0
  lnor(1) = dlxdxi(2, 1)*dlxdxi(3, 2) - dlxdxi(3, 1)*dlxdxi(2, 2)
  lnor(2) = dlxdxi(3, 1)*dlxdxi(1, 2) - dlxdxi(1, 1)*dlxdxi(3, 2)
  lnor(3) = dlxdxi(1, 1)*dlxdxi(2, 2) - dlxdxi(2, 1)*dlxdxi(1, 2)

  DetJb = sqrt(sum(lnor*lnor))
  lnor(:) = lnor(:)/DetJb

  if (DetJb < 1.0d-15) write (*, *) DetJb

!  ! Compute "interior" part
!   igp = 0d0
!   call eval_SHAPE_tet(nshl,igp,xl,dl,shlu,shgradgu,dxidx,Gij,Ginv)

!   do i = 1, NSD
!      diff(i) = sum(bxdl(:,i)*bshlu)-sum((xl(:,i) + dl(:,i))*shlu)
!   enddo

!   do i = 1, NSD
!      igp(i) =  sum(dxidx(i,:)*diff)
!   end do

  call eval_SHAPE_fem(nshl, igp, xl, dl, shlu, shgradgu, dxidx, &
                      Gij, Ginv)

end subroutine eval_faceshape_tet

!======================================================================
!
!======================================================================
subroutine linhexshl(gp, NSHL, NSD, shl, shgradl)
  implicit none

  integer, intent(in) :: NSHL, NSD
  real(8) :: gp(NSD), zi, eta, zeta, eighth
  real(8) :: shl(NSHL), shgradl(NSHL, NSD)

  zi = gp(1)
  eta = gp(2)
  zeta = gp(3)

  eighth = 0.125d0

  shl(1) = (zi - 1d+0)*(eta + 1d+0)*(zeta - 1d+0)*eighth
  shl(2) = -(zi + 1d+0)*(eta + 1d+0)*(zeta - 1d+0)*eighth
  shl(3) = (zi + 1d+0)*(eta + 1d+0)*(zeta + 1d+0)*eighth
  shl(4) = -(zi - 1d+0)*(eta + 1d+0)*(zeta + 1d+0)*eighth
  shl(5) = -(zi - 1d+0)*(eta - 1d+0)*(zeta - 1d+0)*eighth
  shl(6) = (zi + 1d+0)*(eta - 1d+0)*(zeta - 1d+0)*eighth
  shl(7) = -(zi + 1d+0)*(eta - 1d+0)*(zeta + 1d+0)*eighth
  shl(8) = (zi - 1d+0)*(eta - 1d+0)*(zeta + 1d+0)*eighth

  shgradl(1, 1) = (eta + 1d+0)*(zeta - 1d+0)*eighth
  shgradl(2, 1) = -(eta + 1d+0)*(zeta - 1d+0)*eighth
  shgradl(3, 1) = (eta + 1d+0)*(zeta + 1d+0)*eighth
  shgradl(4, 1) = -(eta + 1d+0)*(zeta + 1d+0)*eighth
  shgradl(5, 1) = -(eta - 1d+0)*(zeta - 1d+0)*eighth
  shgradl(6, 1) = (eta - 1d+0)*(zeta - 1d+0)*eighth
  shgradl(7, 1) = -(eta - 1d+0)*(zeta + 1d+0)*eighth
  shgradl(8, 1) = (eta - 1d+0)*(zeta + 1d+0)*eighth

  shgradl(1, 2) = (zi - 1d+0)*(zeta - 1d+0)*eighth
  shgradl(2, 2) = -(zi + 1d+0)*(zeta - 1d+0)*eighth
  shgradl(3, 2) = (zi + 1d+0)*(zeta + 1d+0)*eighth
  shgradl(4, 2) = -(zi - 1d+0)*(zeta + 1d+0)*eighth
  shgradl(5, 2) = -(zi - 1d+0)*(zeta - 1d+0)*eighth
  shgradl(6, 2) = (zi + 1d+0)*(zeta - 1d+0)*eighth
  shgradl(7, 2) = -(zi + 1d+0)*(zeta + 1d+0)*eighth
  shgradl(8, 2) = (zi - 1d+0)*(zeta + 1d+0)*eighth

  shgradl(1, 3) = (zi - 1d+0)*(eta + 1d+0)*eighth
  shgradl(2, 3) = -(zi + 1d+0)*(eta + 1d+0)*eighth
  shgradl(3, 3) = (zi + 1d+0)*(eta + 1d+0)*eighth
  shgradl(4, 3) = -(zi - 1d+0)*(eta + 1d+0)*eighth
  shgradl(5, 3) = -(zi - 1d+0)*(eta - 1d+0)*eighth
  shgradl(6, 3) = (zi + 1d+0)*(eta - 1d+0)*eighth
  shgradl(7, 3) = -(zi + 1d+0)*(eta - 1d+0)*eighth
  shgradl(8, 3) = (zi - 1d+0)*(eta - 1d+0)*eighth
end subroutine linhexshl

!======================================================================
!
!======================================================================
subroutine lintetshl(gp, NSHL, NSD, shl, shgradl)
  implicit none

  integer :: NSHL, NSD
  real(8) :: gp(NSD), zi, eta, zeta
  real(8) :: shl(NSHL), shgradl(NSHL, NSD)

  zi = gp(1)
  eta = gp(2)
  zeta = gp(3)

  shl(1) = zi
  shl(2) = eta
  shl(3) = zeta
  shl(4) = 1.0d0 - zi - eta - zeta

  shgradl = 0.0d0

  shgradl(1, 1) = 1.0d0
  shgradl(2, 2) = 1.0d0
  shgradl(3, 3) = 1.0d0
  shgradl(4, 1) = -1.0d0
  shgradl(4, 2) = -1.0d0
  shgradl(4, 3) = -1.0d0
end subroutine lintetshl

!======================================================================
!
!======================================================================
subroutine linquadshl(gp, NSHL, NSD, shl, shgradl)
  implicit none

  integer :: NSHL, NSD
  real(8) :: gp(2), zi, eta, fourth
  real(8) :: shl(NSHL), shgradl(NSHL, 2)

  zi = gp(1)
  eta = gp(2)

  fourth = 0.25d0

  shl(1) = (zi - 1d+0)*(eta - 1d+0)*fourth
  shl(2) = -(zi + 1d+0)*(eta - 1d+0)*fourth
  shl(3) = (zi + 1d+0)*(eta + 1d+0)*fourth
  shl(4) = -(zi - 1d+0)*(eta + 1d+0)*fourth

  shgradl(1, 1) = (eta - 1d+0)*fourth
  shgradl(2, 1) = -(eta - 1d+0)*fourth
  shgradl(3, 1) = (eta + 1d+0)*fourth
  shgradl(4, 1) = -(eta + 1d+0)*fourth

  shgradl(1, 2) = (zi - 1d+0)*fourth
  shgradl(2, 2) = -(zi + 1d+0)*fourth
  shgradl(3, 2) = (zi + 1d+0)*fourth
  shgradl(4, 2) = -(zi - 1d+0)*fourth
end subroutine linquadshl

!======================================================================
! shape functions for linear triangle
!======================================================================
subroutine lintrishl(gp, NSHL, NSD, shl, shgradl)
  implicit none
  integer, intent(in)  :: NSHL, NSD
  real(8), intent(in)  :: gp(2)
  real(8), intent(out) :: shl(NSHL), shgradl(NSHL, 2)
  real(8) :: zi, eta

  zi = gp(1)
  eta = gp(2)

  shl(1) = zi
  shl(2) = eta
  shl(3) = 1.0d0 - zi - eta

  shgradl = 0.0d0

  shgradl(1, 1) = 1.0d0
  shgradl(2, 2) = 1.0d0

  shgradl(3, 1) = -1.0d0
  shgradl(3, 2) = -1.0d0
end subroutine lintrishl

!======================================================================
! subroutine computes the volumn of a tetrahedron
!======================================================================
subroutine voltet(x1, x2, x3, x4, vol)
  implicit none
  real(8), intent(in)  :: x1(3), x2(3), x3(3), x4(3)
  real(8), intent(out) :: vol

  vol = x2(1)*x3(2)*x4(3) - x2(1)*x3(3)*x4(2) - x3(1)*x2(2)*x4(3) &
        + x3(1)*x2(3)*x4(2) + x4(1)*x2(2)*x3(3) - x4(1)*x2(3)*x3(2) &
        - x1(1)*x3(2)*x4(3) + x1(1)*x3(3)*x4(2) + x3(1)*x1(2)*x4(3) &
        - x3(1)*x1(3)*x4(2) - x4(1)*x1(2)*x3(3) + x4(1)*x1(3)*x3(2) &
        + x1(1)*x2(2)*x4(3) - x1(1)*x2(3)*x4(2) - x2(1)*x1(2)*x4(3) &
        + x2(1)*x1(3)*x4(2) + x4(1)*x1(2)*x2(3) - x4(1)*x1(3)*x2(2) &
        - x1(1)*x2(2)*x3(3) + x1(1)*x2(3)*x3(2) + x2(1)*x1(2)*x3(3) &
        - x2(1)*x1(3)*x3(2) - x3(1)*x1(2)*x2(3) + x3(1)*x1(3)*x2(2)

  vol = vol/6.0d0
end subroutine voltet

!======================================================================
!
!======================================================================
subroutine artri(v1, v2, v3, w1, w2, w3, thck, vol)

  implicit none
  real(8) :: v1, v2, v3, w1, w2, w3, thck, vol

  vol = sqrt((v2*w3 - w2*v3)**2 + (v1*w3 - w1*v3)**2 + (v1*w2 - w1*v2)**2)
  vol = thck*vol/2d0

end subroutine artri

!======================================================================
!
!======================================================================
subroutine eval_faceshape_pri(nshl, bgp, igp, for, xl, dl, shlu, shgradgu, &
                              dxidx, lnor)
  use aAdjKeep
  use commonvars
  use mpi

  implicit none

  integer, intent(in) :: nshl
  real(8), intent(in) :: bgp(2), igp(NSD)

  real(8) :: dl(NSHL, NSD), xl(NSHL, NSD)
  real(8) :: shlu(NSHL), shgradgu(NSHL, NSD), dxidx(NSD, NSD), lnor(NSD)
  integer :: for

  real(8) :: bshlu(NSHLbmax), bshgradlu(NSHLbmax, 2)
  real(8) :: bxdl(NSHLbmax, NSD)
  real(8) :: dlxdxi(NSD, 2)
  real(8) :: diff(NSD)
  real(8) :: Gij(NSD, NSD), Ginv(NSD, NSD)
  integer :: i, j, dd

!!!   igp = 0d0
  if (for == 1) then
    bxdl(1, :) = xl(3, :) + dl(3, :)
    bxdl(2, :) = xl(2, :) + dl(2, :)
    bxdl(3, :) = xl(1, :) + dl(1, :)
!!!       igp(3) = -1d0
  else if (for == 2) then
    bxdl(1, :) = xl(4, :) + dl(4, :)
    bxdl(2, :) = xl(5, :) + dl(5, :)
    bxdl(3, :) = xl(6, :) + dl(6, :)
!!!       igp(3) = 1d0
!!!       write(*,*) "Top-down prism"
  end if

  ! Get basis functions and local ders
  call lintrishl(bgp, NSHLbmax, NSD, bshlu, bshgradlu)

  ! Now calculate gradients.
  ! Calculate dx/dxi
  dlxdxi = 0.0d0
  do i = 1, NSHLbmax
    do dd = 1, NSD
      dlxdxi(dd, :) = dlxdxi(dd, :) + bxdl(i, dd)*bshgradlu(i, :)
    end do
  end do

  lnor = 0.0d0
  lnor(1) = dlxdxi(2, 1)*dlxdxi(3, 2) - dlxdxi(3, 1)*dlxdxi(2, 2)
  lnor(2) = dlxdxi(3, 1)*dlxdxi(1, 2) - dlxdxi(1, 1)*dlxdxi(3, 2)
  lnor(3) = dlxdxi(1, 1)*dlxdxi(2, 2) - dlxdxi(2, 1)*dlxdxi(1, 2)

  DetJb = sqrt(sum(lnor*lnor))
  lnor(:) = lnor(:)/DetJb

  if (DetJb < 1.0d-15) write (*, *) DetJb

!   ! Compute "interior" part
!   !igp = 0d0
!   call eval_SHAPE_fem(nshl,igp,xl,dl,shlu,shgradgu,dxidx,Gij,Ginv)
!
!   do i = 1, NSD
!     diff(i) = sum(bxdl(:,i)*bshlu)-sum((xl(:,i) + dl(:,i))*shlu)
!   enddo
!
!   do while ((sum(diff**2)).gt.1d-10)
!
!     do i = 1, NSD
!        igp(i) = igp(i)  + sum(dxidx(i,:)*diff)
!     end do
!
!     igp(1) = max(0d0, min(igp(1),1d0))
!     igp(2) = max(0d0, min(igp(2),1d0))
!     !! 4th constraint!!!
!     igp(3) = max(-1d0,min(igp(3),1d0))
!
  call eval_SHAPE_fem(nshl, igp, xl, dl, shlu, shgradgu, dxidx, Gij, Ginv)
!
!     do i = 1, NSD
!       diff(i) = sum(bxdl(:,i)*bshlu)-sum((xl(:,i) + dl(:,i))*shlu)
!     enddo
!   enddo

  if (DetJ < 0.0d0) then
    write (*, *) "******** BOUNDARY ****** "
    write (*, *) "Element type: Prism   "
    write (*, *) "igp = ", igp
    do i = 1, NSD
      diff(i) = sum(bxdl(:, i)*bshlu) - sum((xl(:, i) + dl(:, i))*shlu)
    end do
    write (*, *) "diff = ", diff
    do i = 1, NSD
      diff(i) = sum(bxdl(:, i)*bshlu)
    end do
    write (*, *) "x2D = ", diff

    do i = 1, NSD
      diff(i) = sum((xl(:, i) + dl(:, i))*shlu)
    end do
    write (*, *) "x3D = ", diff

    write (*, *) "bdxl = "
    write (*, *) bxdl

    write (*, *) "bshlu = ", bshlu

    write (*, *) "xl = "
    write (*, *) xl

    write (*, *) "shlu = ", shlu
    DetJ = -1d+0*DetJ
  end if
end subroutine eval_faceshape_pri

!=======================================================================
! shape function for prism
!=======================================================================
subroutine linprishl(gp, NSHL, NSD, shl, shgradl)
  implicit none

  integer, intent(in) :: NSHL, NSD
  real(8) :: gp(NSD), zi, eta, zeta, sixth
  real(8) :: shl(NSHL), shgradl(NSHL, NSD)

  zi = gp(1)
  eta = gp(2)
  zeta = gp(3)

  shl(1) = 0.5d0*zi*(1d+0 - zeta)
  shl(2) = 0.5d0*eta*(1d+0 - zeta)
  shl(3) = 0.5d0*(1d+0 - zi - eta)*(1d+0 - zeta)
  shl(4) = 0.5d0*zi*(1d+0 + zeta)
  shl(5) = 0.5d0*eta*(1d+0 + zeta)
  shl(6) = 0.5d0*(1d+0 - zi - eta)*(1d+0 + zeta)

  shgradl = 0.0d0
  shgradl(1, 1) = 0.5d0*(1d+0 - zeta)
  shgradl(2, 1) = 0.0d0
  shgradl(3, 1) = -0.5d0*(1d+0 - zeta)
  shgradl(4, 1) = 0.5d0*(1d+0 + zeta)
  shgradl(5, 1) = 0.0d0
  shgradl(6, 1) = -0.5d0*(1d+0 + zeta)

  shgradl(1, 2) = 0.0d0
  shgradl(2, 2) = 0.5d0*(1d+0 - zeta)
  shgradl(3, 2) = -0.5d0*(1d+0 - zeta)
  shgradl(4, 2) = 0.0d0
  shgradl(5, 2) = 0.5d0*(1d+0 + zeta)
  shgradl(6, 2) = -0.5d0*(1d+0 + zeta)

  shgradl(1, 3) = -0.5d0*zi
  shgradl(2, 3) = -0.5d0*eta
  shgradl(3, 3) = -0.5d0*(1d+0 - zi - eta)
  shgradl(4, 3) = 0.5d0*zi
  shgradl(5, 3) = 0.5d0*eta
  shgradl(6, 3) = 0.5d0*(1d+0 - zi - eta)
end subroutine linprishl
