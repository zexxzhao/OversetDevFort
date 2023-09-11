!======================================================================
!
!======================================================================
subroutine e3STAB_3D(NSD, Gij, dt, uadvi, rhoi, mui, &
                     tauM, tauP, tauC)
  ! use aAdjKeep
  ! use commonvars
  implicit none

  integer, intent(in)  :: NSD
  real(8), intent(in)  :: Gij(NSD, NSD), uadvi(NSD),&
                          dt, rhoi, mui
  real(8), intent(out) :: tauM, tauC, tauP

  real(8) :: taua, taut, taud, m_k, nu
  integer :: i, j

  ! dtfact = 4.0d0

!!!  dtfact = 16.0d0
!!!  dtfact = 64.0d0
!!!  dtfact = 256.0d0

  m_k = 12.0d0

  ! get tauM, tauC, cfl(1)
  taua = 0.0d0
  do j = 1, NSD
    do i = 1, NSD
      taua = taua + uadvi(i)*Gij(i, j)*uadvi(j)
    end do
  end do

  nu = mui/rhoi
  taud = nu * nu * m_k * sum(Gij(:, :) * Gij(:, :))
  taut = 4.0d0 / (dt * dt)

  ! Get tauM
  tauM = 1.0d0/(rhoi*sqrt(taut + taua + taud))
  tauP = tauM

  ! Get tauC
  tauC = rhoi*sqrt(taua + taud)/(Gij(1, 1) + Gij(2, 2) + Gij(3, 3))
end subroutine e3STAB_3D


!======================================================================
! Stabilization parameters (tauM, tauC) for VMS
!======================================================================
subroutine e3STAB_3D_OLD(Gij, Ginv, AD_VEL_L, AD_VEL_L1, rLi, rhoi, mui, &
                     tauM, tauP, tauLS, &
                     tauC, tauBar, tauBar1, uprime, cfl)
  use aAdjKeep
  use commonvars
  implicit none

  real(8), intent(in)  :: Gij(NSD, NSD), Ginv(NSD, NSD), &
                          AD_VEL_L(NSD), AD_VEL_L1(NSD), rLi(NSD), &
                          rhoi, mui
  real(8), intent(out) :: tauM, tauC, tauP, tauBar, tauBar1, &
                          cfl(2), uprime(NSD), tauLS

  real(8) :: taua, taua1, taut, taud, m_k, dtfact, nu, taus, gij2
  real(8) :: a1, a2, a3
  integer :: i, j

  dtfact = 4.0d0

!!!  dtfact = 16.0d0
!!!  dtfact = 64.0d0
!!!  dtfact = 256.0d0

  m_k = 12d0
!!!  m_k    = 36.0d0!*real(max(Pu**2,Qu**2,Ru**2))

  cfl = 0.0d0
!  ! get cfl(2)
!  taua  = 0.0d0
!  do j = 1, NSD
!    do i = 1, NSD
!      taua = taua + AD_VEL_L(i)*Ginv(i,j)*AD_VEL_L(j)
!    end do
!  end do
!  cfl(2)  = Delt*sum(AD_VEL_L*AD_VEL_L)/sqrt(taua+1d-8)

  ! get tauM, tauC, cfl(1)
  taua = 0.0d0
  taut = 0.0d0
  taud = 0.0d0
  gij2 = 0.0d0
  do j = 1, NSD
    do i = 1, NSD
      taua = taua + AD_VEL_L(i)*Gij(i, j)*AD_VEL_L(j)
      taua1 = taua1 + AD_VEL_L1(i)*Gij(i, j)*AD_VEL_L1(j)
      gij2 = gij2 + Gij(i, j)*Gij(i, j)
    end do
  end do

  cfl(1) = Delt*sqrt(taua)

  nu = mui/rhoi
  taud = nu*nu*m_k*gij2
  taut = dtfact*Dtgl*Dtgl

  taus = taua + taud

  ! Get tauM - so far a scalar
  tauM = 1.0d0/(rhoi*sqrt(taus + taut))
  tauLS = 1.0d0/(sqrt(taua1 + taut + kappa*kappa*m_k*gij2))
  tauP = tauM

  ! Get tauC
  tauC = rhoi*sqrt(taus)/(Gij(1, 1) + Gij(2, 2) + Gij(3, 3))
!!!  tauC = rho*sqrt(taus)/sqrt(gij2)
!!!  tauC = 1.0d0/(tauM*(Gij(1,1)+Gij(2,2)+Gij(3,3)))

  ! Tayfun's tauC
!!!  tauC = rho**2*tauM*sum(AD_VEL_L*AD_VEL_L)

  uprime(:) = -tauM*rLi(:)

  ! Get TauBar
  tauBar = 0.0d0
!  do i = 1, NSD
!    do j = 1, NSD
!      tauBar = tauBar + uprime(i)*Gij(i,j)*uprime(j)
!    end do
!  end do
!  tauBar = 1.0d0/sqrt(tauBar+1.0d-15)
  a1 = taut + taua + taud
  a3 = taut + taua1 + kappa*kappa*m_k*gij2
  a2 = dtfact*Dtgl
  tauBar1 = 0.0d0
  !tauBar1 = -a2/(a1*sqrt(a3) + a3*sqrt(a1))
  !tauBar1 = tauBar1*cross_flag
!  tauBar1 = 0d0
end subroutine e3STAB_3D_OLD

!======================================================================
!
!======================================================================
subroutine e3DC_beta3(ui, duidxi, Gij, Ginv, res, dxidx, tauM, kdc)

  use aAdjKeep
  use commonvars
  implicit none
  real(8), intent(in) :: Gij(NSD, NSD), Ginv(NSD, NSD), tauM
  real(8) :: ui(NSD), duidxi(NSD, NSD), res(NSD), He, dxidx(NSD, NSD)
  integer :: i, j, k
  real(8) :: tmp1, tmp2, resnorm, kdc, h2, gunrm
  real(8) :: gabsu(NSD), absu

  He = 1d0!!!(rho - rhoa)/(rhow - rhoa)

  !  Meshsize
  tmp1 = 0d0
  tmp2 = 0d0
  do j = 1, NSD
    do i = 1, NSD
      tmp1 = tmp1 + Gij(i, j)*Gij(i, j)
      tmp2 = tmp2 + Ginv(i, j)*Ginv(i, j)
    end do
  end do

  h2 = 1d0/sqrt(Gij(1, 1) + Gij(2, 2) + Gij(3, 3))
  resnorm = sqrt(sum(res*res))

  kdc = h2*h2*resnorm

end subroutine e3DC_beta3

!======================================================================
!
!======================================================================
subroutine e3DC_beta2(NSD, ns_kdc, Gij, res, kdc)

  !use aAdjKeep
  !use commonvars

  implicit none

  integer, intent(in) :: NSD
  integer, intent(in) :: ns_kdc
  real(8), intent(in) :: Gij(NSD, NSD), res(NSD)
  real(8), intent(out) :: kdc

  kdc = ns_kdc * sqrt(sum(res(:)*res(:)) / sum(Gij(:, :)*Gij(:, :)))

end subroutine e3DC_beta2
!======================================================================
!
!======================================================================
subroutine e3DC_scalar(NSD, c_kdc, Gij, res, gradc, kdc)

  !use aAdjKeep
  !use commonvars

  implicit none

  integer, intent(in) :: NSD
  real(8), intent(in) :: c_kdc, Gij(NSD, NSD)
  real(8), intent(in) :: res, gradc(NSD)
  real(8) :: tmp, kdc
  integer :: i, j

  tmp = 0d0

  do i = 1, NSD
    do j = 1, NSD
      tmp = tmp + Gij(i, j)*gradc(i)*gradc(j)
    enddo
  enddo

  kdc = c_kdc * abs(res) / sqrt(tmp + 1d-12)
  ! kdc = c_kdc * abs(res) / sqrt(sum(Gij(:,:) * Gij(:, :)))

end subroutine e3DC_scalar

!======================================================================
!
!======================================================================
! subroutine e3DC_shakibNS(ui, duidxi, Ginv, res, dxidx, KAPPA_DC)
! 
!   use aAdjKeep
!   use commonvars
!   implicit none
! 
!   real(8), intent(in) :: Ginv(NSD, NSD)
!   real(8) :: ui(NSD), duidxi(NSD, NSD), res(NSD), He, dxidx(NSD, NSD)
!   real(8) :: KAPPA_DC(NSD, NSD)
!   integer :: i, j
!   real(8) :: resnorm, temp(NSD)
! 
!   He = (rho - rhoa)/(rhow - rhoa)
! 
!   ! Grad norm in reference coordinates
!   temp = 0.0d0
!   do j = 1, NSD
!     do i = 1, NSD
!       temp = temp + duidxi(:, i)*Ginv(i, j)*duidxi(:, j)
!     end do
!   end do
! 
!   ! Residual norm
!   resnorm = sqrt(sum(res(:)*res(:)))
! 
!   ! Kdc
!   KAPPA_DC = Ginv*(He*NS_kdc_w + (1d0 - He)*NS_kdc_a)* &
!              resnorm/(sqrt(sum(temp)) + 1d-3)
! 
! end subroutine e3DC_shakibNS

!======================================================================
! Stabilization parameter (tauB) for weak BC
!======================================================================
subroutine e3bSTAB_weak(tauB, tauNor, ui, nor, dxidx, rhoi, mui)

  use aAdjKeep
  use commonvars
  use mpi
  implicit none

  real(8), intent(in)  :: ui(NSD), nor(NSD), dxidx(NSD, NSD)
  real(8), intent(in)  :: rhoi, mui

  real(8), intent(out) :: tauB, tauNor

  real(8) :: temp(NSD), hn, unor, upos, uneg
  integer :: i

  do i = 1, NSD
    temp(i) = sum(dxidx(i, :)*nor)
  end do

!!!  hn = 1.0d0/(sqrt(sum(temp*temp)))  ! hn for tets
  hn = 2.0d0/(sqrt(sum(temp*temp)))  ! hn for other elements

  tauB = 4.0d0*mui/hn
  ! tauBLS = 4.0d0*kappa/hn

  unor = sum(ui*nor(:))
  upos = 0.5d0*(unor + abs(unor))
  uneg = 0.5d0*(unor - abs(unor))

  tauNor = tauB !+ rho*upos

end subroutine e3bSTAB_weak

!======================================================================
! Stabilization parameters (tauM, tauC) for VMS
!======================================================================
subroutine e3STAB_3D_NSVOF(NSD, Gij, dt, uadvi, rhoi, mui, &
                           tauM, tauP, tauC, tauLS)
  ! use aAdjKeep
  ! use commonvars
  implicit none

  integer, intent(in)  :: NSD
  real(8), intent(in)  :: Gij(NSD, NSD), uadvi(NSD),&
                          dt, rhoi, mui
  real(8), intent(out) :: tauM, tauC, tauP, tauLS

  real(8) :: taua, taut, taud, m_k, nu
  integer :: i, j

  ! dtfact = 4.0d0

!!!  dtfact = 16.0d0
!!!  dtfact = 64.0d0
!!!  dtfact = 256.0d0

  m_k = 12.0d0

  ! get tauM, tauC, cfl(1)
  taua = 0.0d0
  do j = 1, NSD
    do i = 1, NSD
      taua = taua + uadvi(i)*Gij(i, j)*uadvi(j)
    end do
  end do

  nu = mui/rhoi
  taud = nu * nu * m_k * sum(Gij(:, :) * Gij(:, :))
  taut = 4.0d0 / (dt * dt)

  ! Get tauM
  tauM = 1.0d0/(rhoi*sqrt(taut + taua + taud))
  tauLS = 1.0d0/(sqrt(taut + taua))
  tauP = tauM

  ! Get tauC
  tauC = rhoi*sqrt(taua + taud)/(Gij(1, 1) + Gij(2, 2) + Gij(3, 3))
end subroutine e3STAB_3D_NSVOF

!======================================================================
! Stabilization parameter tauBar for Tem
!======================================================================
subroutine e3STAB_3D_TEM(NSD, Gij, uadvi, dt, rhoi, cpi, hki, &
                         tau_tem)
! use aAdjKeep
! use commonvars
implicit none

integer, intent(in)  :: NSD
real(8), intent(in)  :: Gij(NSD, NSD), uadvi(NSD)
real(8), intent(in)  :: dt, rhoi, cpi, hki
real(8), intent(out) :: tau_tem

real(8) :: m_k, k
integer :: i, j

m_k = 12.0d0

tau_tem = 0.0d0
do j = 1, NSD
  do i = 1, NSD
    tau_tem = tau_tem + uadvi(i)*Gij(i, j)*uadvi(j)
  end do
end do
tau_tem = tau_tem + 4.0d0 / (dt * dt)

k = hki / (rhoi * cpi)
tau_tem = tau_tem + k * k * m_k * sum(Gij(:, :) * Gij(:, :))

! Get tauM
tau_tem = 1.0d0/(rhoi * cpi * sqrt(tau_tem))
end subroutine e3STAB_3D_TEM
!======================================================================
! Stabilization parameter tauBar for NSVOF
!======================================================================
subroutine e3STAB_3D_NSVOF_TAUBAR(NSD, Gij, uadvi, uprime, tauBar)
  implicit none

  integer, intent(in)  :: NSD
  real(8), intent(in)  :: Gij(NSD, NSD), uadvi(NSD), uprime(NSD)

  real(8), intent(out) :: tauBar

  integer :: i, j

  ! Get TauBar
  tauBar = 0.0d0
  do i = 1, NSD
    do j = 1, NSD
      tauBar = tauBar + uprime(i)*Gij(i,j)*uprime(j)
    end do
  end do
  tauBar = 1.0d0/sqrt(tauBar+1.0d-15)
end subroutine e3STAB_3D_NSVOF_TAUBAR

!======================================================================
! Calculate the CFL number
!======================================================================
subroutine e3CFL(NSD, uadvi, Gij, dt, cfl)
  implicit none
  integer, intent(in)  :: NSD
  real(8), intent(in)  :: uadvi(NSD), Gij(NSD, NSD)
  real(8), intent(in)  :: dt

  real(8), intent(out) :: cfl

  integer :: i, j

  cfl = 0.0d0
  do i = 1, NSD
    do j = 1, NSD
      cfl = cfl + uadvi(i)*Gij(i, j)*uadvi(j)
    end do
  end do
  cfl = dt*sqrt(cfl)

end subroutine e3CFL
