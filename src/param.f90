subroutine readinlet()
  use aAdjKeep
  use commonvars
  use commonpars
  use mpi
  implicit none

  integer :: inode, iel
  open (11, file='bmesh.1.dat', status='old')
  read (11, *)
  read (11, *) InNNODE1, InCELL
  if (ismaster) write (*, *) InNNODE1, InCELL
  allocate (xg_inlet(InNNODE1, 3))
  allocate (ien_inlet(InCELL, 3))
  do inode = 1, InNNODE1
    read (11, *) xg_inlet(inode, :)
  end do
  do iel = 1, InCELL
    read (11, *) ien_inlet(iel, :)
  end do
  close (11)

  allocate (inflow_velocity(InNNODE1, 3), inflow_x(InNNODE1), inflow_y(InNNODE1), inflow_z(InNNODE1), &
            inflow_bnodes(InNNODE1), inflow_velocity_check(InNNODE1, 3), inflow_pressure(InNNODE1), inflow_rho(InNNODE1))
end subroutine readinlet

subroutine calarea(xa, xb, xc, area)
  implicit none
  real(8), intent(in) :: xa(3), xb(3), xc(3)
  real(8), intent(out) :: area
  real(8) :: A(2, 2)
  A(1, 1) = xb(2) - xa(2)
  A(1, 2) = xb(3) - xa(3)
  A(2, 1) = xc(2) - xa(2)
  A(2, 2) = xc(3) - xa(3)
  area = 0.5d0*abs(A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))
end subroutine calarea
!======================================================================
! subroutine to get the inflow boundary conditions by Jinhui
!======================================================================
subroutine getinflow_domain(xt, utmp, ptmp, rhotmp)

  use aAdjKeep
  use commonvars
  use commonpars
  use mpi

  implicit none
  real(8), intent(in) :: xt(3)
  real(8), intent(out) :: utmp(3), ptmp, rhotmp
  integer :: num_plane
  real(8) :: areaa, areab, areac, areasum
  integer :: inode, iel, iel_found, i
  real(8) :: xa(3), xb(3), xc(3)
  real(8) :: r, s
  do iel = 1, InCELL
    xa(:) = xg_inlet(ien_inlet(iel, 1), :)
    xb(:) = xg_inlet(ien_inlet(iel, 2), :)
    xc(:) = xg_inlet(ien_inlet(iel, 3), :)
    areaa = 0d0
    areab = 0d0
    areac = 0d0
    areasum = 0d0
    call calarea(xa, xb, xc, areasum)
    call calarea(xt, xa, xb, areac)
    call calarea(xt, xa, xc, areab)
    call calarea(xt, xb, xc, areaa)
    iel_found = -1
    ! write(*,*)  areaa,  areab,  areac, areasum
    if (abs(areaa + areab + areac - areasum) < 1d-7) then
      iel_found = iel
      r = areaa/areasum
      s = areab/areasum
      exit
    end if
  end do
  if (iel > InCELL) then
    write (*, *) "did not find the elem:", myid, xt(:)
  end if
  if (ismaster) then
    write (*, *) iel_found, abs(areaa + areab + areac - areasum)
    write (*, *) r, s
  end if
  open (100000, file='inflow.1.dat', status='old')
  do i = 1, InNNODE1
    read (100000, *) inflow_x(i), inflow_y(i), inflow_z(i), &
      inflow_velocity(i, 1), inflow_velocity(i, 2), inflow_velocity(i, 3), &
      inflow_pressure(i), inflow_rho(i)
  end do

  utmp(:) = inflow_velocity(ien_inlet(iel_found, 1), :)*r &
            + inflow_velocity(ien_inlet(iel_found, 2), :)*s &
            + inflow_velocity(ien_inlet(iel_found, 3), :)*(1 - r - s)

  ptmp = inflow_pressure(ien_inlet(iel_found, 1))*r &
         + inflow_pressure(ien_inlet(iel_found, 2))*s &
         + inflow_pressure(ien_inlet(iel_found, 3))*(1 - r - s)

  rhotmp = inflow_rho(ien_inlet(iel_found, 1))*r &
           + inflow_rho(ien_inlet(iel_found, 2))*s &
           + inflow_rho(ien_inlet(iel_found, 3))*(1 - r - s)

  close (100000)

end subroutine getinflow_domain

!======================================================================
! subroutine to get the inflow boundary conditions by Jinhui
!======================================================================
subroutine getinflow(istep, Rstep)

  use aAdjKeep
  use commonvars
  use commonpars
  use mpi

  implicit none

  integer, intent(in)::istep, Rstep
  integer :: num_plane
  integer :: i
!   For the first inflow
!   ********************

  if (istep > 0) then

    max_planes = 2048

    num_plane = floor(x_inflow/Delt_x)
    index_plane_p = mod(num_plane, max_planes) + 1

    remain_x = x_inflow - num_plane*Delt_x
    interpolation_factor = remain_x/Delt_x

    if (index_plane_p < 10) then
      write (filename_p, "(A7,I1,A4)") "inflow.", index_plane_p, ".dat"
    end if
    if ((index_plane_p >= 10) .and. (index_plane_n < 100)) then
      write (filename_p, "(A7,I2,A4)") "inflow.", index_plane_p, ".dat"
    end if
    if ((index_plane_p >= 100) .and. (index_plane_p < 1000)) then
      write (filename_p, "(A7,I3,A4)") "inflow.", index_plane_p, ".dat"
    end if
    if (index_plane_p >= 1000) then
      write (filename_p, "(A7,I4,A4)") "inflow.", index_plane_p, ".dat"
    end if
    if (ismaster) then
      write (*, *) "The previous plane"
      print *, trim(filename_p)
    end if

    index_plane_n = index_plane_p + 1

    if (index_plane_n < 10) then
      write (filename_n, "(A7,I1,A4)") "inflow.", index_plane_n, ".dat"
    end if
    if ((index_plane_n >= 10) .and. (index_plane_n < 100)) then
      write (filename_n, "(A7,I2,A4)") "inflow.", index_plane_n, ".dat"
    end if
    if ((index_plane_n >= 100) .and. (index_plane_n < 1000)) then
      write (filename_n, "(A7,I3,A4)") "inflow.", index_plane_n, ".dat"
    end if
    if (index_plane_n > 1000) then
      write (filename_n, "(A7,I4,A4)") "inflow.", index_plane_n, ".dat"
    end if

    if (ismaster) then
      write (*, *) "The next plane"
      print *, trim(filename_n)
      write (*, *) "num_plane", num_plane
      write (*, *) "x_inflow,", x_inflow
      write (*, *) "remain_x", remain_x
      write (*, *) "Delt_x", Delt_x
      write (*, *) "interpolation_factor", interpolation_factor
      write (*, *) "InNNODE1:", InNNODE1
      write (*, *) "---------------------------------------"
    end if

    ! open (100000, file=filename_p, status='old', iostat=ios)
    ! open (100001, file=filename_n, status='old', iostat=ios)
    open (100000, file=filename_p, status='old')
    open (100001, file=filename_n, status='old')
    do i = 1, InNNODE1
      !  read(100000,*) Inflow_velocity_p(1),Inflow_velocity_p(2),Inflow_velocity_p(3)
      !  read(100001,*) Inflow_velocity_n(1),Inflow_velocity_n(2),Inflow_velocity_n(3)
      read (100000, *) inflow_x(i), inflow_y(i), inflow_z(i), inflow_velocity_p(1), &
        inflow_velocity_p(2), inflow_velocity_p(3), inflow_pressure_p, inflow_rho_p
      read (100001, *) inflow_x(i), inflow_y(i), inflow_z(i), inflow_velocity_n(1), &
        inflow_velocity_n(2), inflow_velocity_n(3), inflow_pressure_n, inflow_rho_n
      inflow_velocity(i, :) = (1 - interpolation_factor)*inflow_velocity_p + interpolation_factor*inflow_velocity_n
      inflow_rho_p = inflow_rho_p + 1
      inflow_rho_n = inflow_rho_n + 1
      inflow_pressure(i) = (1 - interpolation_factor)*inflow_pressure_p + interpolation_factor*inflow_pressure_n
      inflow_rho(i) = (1 - interpolation_factor)*inflow_rho_p + interpolation_factor*inflow_rho_n
    end do
    close (100000)
    close (100001)
  end if

end subroutine getinflow

!======================================================================
! subroutine to read in the parameters defined in 'param.dat'
!======================================================================
subroutine getparam()

  use aAdjKeep
  use commonvars
  use commonpars
  use mpi
  implicit none

  integer :: inurbs, i, j, bb
  character(len=30) :: fname(2), fnum(3)

  ! Time stepping
  call iread("Nstep", Nstep)
  call iread("ifq", ifq)
  call rread("Delt", Delt)
  call rread("rhoinf", rhoinf)
  ! Physics
  call rread("DENS_AIR", rhoa)
  call rread("DENS_WATER", rhow)
  call rread("VISC_AIR", mua)
  call rread("VISC_WATER", muw)
  call rread("CP_AIR", cpa)
  call rread("CP_WATER", cpw)
  call rread("HK_AIR", kappaa)
  call rread("HK_WATER", kappaw)
  call rread("DENS_SOLID", rhos)
  call rread("CP_SOLID", cps)
  call rread("HK_SOLID", kappas)
  call rread("T_SAT", Ts)
  call rread("LATENT_HEAT", lh)
  call rread("C_COND", c_cond)
  call rread("C_EVAP", c_evap)

  call rread("HTC_AIR", htca)
  call rread("HTC_WATER", htcw)

  call rread('TRM', Trm)
  call rread('TREF', Tref)

  call rread('OMEGA_X', omega_x)

  call iread("SUBDOMAIN_ID1", ELM_ID_SET(1))
  call iread("SUBDOMAIN_ID2", ELM_ID_SET(2))
  call iread("BMESH_ID1", BMESH_ID_SET(1))
  call iread("BMESH_ID2", BMESH_ID_SET(2))
  call iread("SRC_ID1", SRC_ID_SET(1))
  call iread("SRC_ID2", SRC_ID_SET(2))
  mp_eps = 3.0d0
  gravity = -9.81d0 ! -1.0d0/(Fr**2.0d0)

  gravvec(1) = 0.0d0
  gravvec(2) = 0.0d0
  gravvec(3) = gravity!Ra/(Pr*Re**2.0d0)

  kappa = 1.0d0/(Sr*Re)!2.15d-4 !For Advection-Diffusion Equation
  kappa = 0d0
!  beta_t = 9.0d-4!3.4d-3
!  phi_inf = 20.0d0

  ! call rread("hull_length", hull_length)
  ! call rread("water_level", water_level)
  ! call rread("position_init", position_init)
  ! water_depth = water_level - domain_bottom
  ! air_height = domain_top - water_level

  ! Problem setup
!!!  call rread("U_in",   Uin)
!!!  call rread("Froude", Froude)
!!!  Uin = Froude*sqrt(gravity*hull_length)
!!!  if (ismaster) write(*,"(a20,x,' = ',x,ES12.4)") "Uin", Uin

!  call iread("VBC(1)", VBC(1))
!  call iread("VBC(2)", VBC(2))
!  call iread("VBC(3)", VBC(3))

!  call iread("MBC(1)", MBC(1))
!  call iread("MBC(2)", MBC(2))
!  call iread("MBC(3)", MBC(3))

!  call iread("BCtype1", BCtype(1))
!  call iread("BCtype2", BCtype(2))
!  call iread("BCtype3", BCtype(3))
!  call iread("BCtype4", BCtype(4))
!  call iread("BCtype5", BCtype(5))
!  call iread("BCtype6", BCtype(6))

!  call iread("BCtype7", BCtype(7))
!  BCtype(8:20) = BCtype(7)
  call iread("NBOUND", NBOUND)
  call iread("NSD", NSD)
  ! BC for "setBCs_CFD"
  allocate (BCugType(NBOUND, NSD), &
            BCugValu(NBOUND, NSD))
  BCugType = 0
  do bb = 1, NBOUND
    i = bb
    write (fnum(1), '(I4)') i
    fname(1) = 'BCugType'//trim(adjustl(fnum(1)))
    call iread3(fname(1), BCugType(i, 1:3))
  end do


  BCugValu = 0.0d0
  do bb = 1, NBOUND
    ! i = bound(bb)%FACE_ID
    i = bb
    do j = 1, NSD
      if (BCugType(i, j) == 1) then
        write (fnum(1), '(I4)') i
        write (fnum(2), '(I4)') j
        fname(1) = 'BCugValu('//trim(adjustl(fnum(1)))//',' &
                   //trim(adjustl(fnum(2)))//')'
        call rread(fname(1), BCugValu(i, j))
      end if
    end do
  end do

  allocate (BCphigType(NBOUND), BCphigValu(NBOUND))
  allocate (BCTgType(NBOUND), BCTgValu(NBOUND))
  BCphigType(:) = 0
  BCTgType(:) = 0
  do bb = 1, NBOUND
    i = bb
    write (fnum(1), '(I4)') i
    fname(1) = 'BCphigType'//trim(adjustl(fnum(1)))
    call iread(fname(1), BCphigType(i))
  end do

  do bb = 1, NBOUND
    i = bb
    if (BCphigType(i) /= 1) cycle
    write (fnum(1), '(I4)') i
    fname(1) = 'BCphigValu'//trim(adjustl(fnum(1)))
    call rread(fname(1), BCphigValu(i))
  end do

  do bb = 1, NBOUND
    i = bb
    write (fnum(1), '(I4)') i
    fname(1) = 'BCTgType'//trim(adjustl(fnum(1)))
    call iread(fname(1), BCTgType(i))
  end do

  do bb = 1, NBOUND
    i = bb
    if (BCTgType(i) /= 1) cycle
    write (fnum(1), '(I4)') i
    fname(1) = 'BCTgValu'//trim(adjustl(fnum(1)))
    call rread(fname(1), BCTgValu(i))
  end do
  ! call rread("usettle", usettle)
  ! Wave generating wall
  ! call rread("wave_amp", wave_amp)
  ! if (wave_amp > 0.0d0) then
    ! call rread("wave_length", wave_length)
    ! call rread("wave_angle", wave_angle)

    ! wave_angle = pi*wave_angle/180.0d0
    ! wave_length = 2.0d0*pi/wave_length
    ! wave_periode = sqrt(gravity*wave_length* &
    !                     dtanh(wave_length/water_level))

    ! if (ismaster) write (*, "(a20,x,' = ',x,ES12.4)") "wave_periode", &
    !  wave_periode

    ! call rread("domain_top", domain_top)
    ! call rread("domain_bottom", domain_bottom)
    ! call rread("domain_left", domain_left)
    ! call rread("domain_right", domain_right)
    ! call rread("domain_front", domain_front)
    ! call rread("domain_back", domain_back)

    ! water_depth = water_level - domain_bottom
  !else
    wave_amp = 0.0d0
    wave_length = 1.0d0
    wave_angle = 0.0d0
    wave_periode = 1.0d0

    domain_top = 0.0d0
    domain_bottom = 0.0d0
    domain_left = 0.0d0
    domain_right = 0.0d0
    domain_front = 0.0d0
    domain_back = 0.0d0
    water_depth = 1.0d0
  !end if

  move_time = 99.0d9

  call iread("fem_flag", fem_flag)
  call iread("USE_VMS", USE_VMS)
  call lread("GRAD_PROP_IN_DIFF", GRAD_PROP_IN_DIFF)
  !call iread("IGA", iga)
  iga = .false.

  ! Navier-stokes solver
  call rread("NS_kdc_w", NS_kdc_w)
  call rread("NS_kdc_a", NS_kdc_a)
  call rread("LSC_kdc", LSC_kdc)
  call rread("TEM_kdc", TEM_kdc)
  call rread("C_DMDOT", C_DMDOT)

  call rread("NS_GMRES_rtol", NS_GMRES_tol)
  call rread("NS_GMRES_atol", NS_GMRES_atol)
  call iread("NS_GMRES_itermax", NS_GMRES_itermax)
  call iread("NS_GMRES_itermin", NS_GMRES_itermin)

  ! call iread("NS_hessian", NS_hess_flag)

  ! Rigid body solver
  ! RB_NL_Ftol = 1.0d-3
  ! RB_NL_Mtol = 1.0d-3

  ! Mesh solver
  ! call rread("Mesh_NL_tol", Mesh_NL_tol)
  ! call rread("Mesh_GMRES_tol", Mesh_GMRES_tol)
  ! call iread("Mesh_GMRES_itermin", Mesh_GMRES_itermin)
  ! call iread("Mesh_GMRES_itermax", Mesh_GMRES_itermax)
  ! Level Set Convection solver

  call rread("LSC_GMRES_atol", LSC_GMRES_atol)
  call rread("LSC_GMRES_rtol", LSC_GMRES_tol)
  call iread("LSC_GMRES_itermax", LSC_GMRES_itermax)
  call iread("LSC_GMRES_itermin", LSC_GMRES_itermin)

  call rread("TEM_GMRES_atol", TEM_GMRES_atol)
  call rread("TEM_GMRES_rtol", TEM_GMRES_tol)
  call iread("TEM_GMRES_itermax", TEM_GMRES_itermax)
  call iread("TEM_GMRES_itermin", TEM_GMRES_itermin)

  call rread("NS_NL_Urtol", NS_NL_Utol)
  call rread("NS_NL_Prtol", NS_NL_Ptol)
  call rread("NS_NL_Uatol", NS_NL_Uatol)
  call rread("NS_NL_Patol", NS_NL_Patol)
  call iread("NS_NL_itermax", NS_NL_itermax)

  call rread("LSC_NL_atol", LSC_NL_atol)
  call rread("LSC_NL_rtol", LSC_NL_tol)


  call rread("TEM_NL_atol", TEM_NL_atol)
  call rread("TEM_NL_rtol", TEM_NL_tol)

  ! Level Set redistance solver
  ! call rread("LSRD_kdc", LSRD_kdc)
  ! call rread("LSRD_penalty_fac", LSRD_penalty_fac)
  ! call rread("LSRD_pseudoDTGL_fac", LSRD_pseudoDTGL_fac)

  ! call rread("LSRD_GMRES_tol", LSRD_GMRES_tol)
  ! call iread("LSRD_GMRES_itermax", LSRD_GMRES_itermax)

  ! call rread("LSRD_NL_tol", LSRD_NL_tol)
  ! call iread("LSRD_NL_itermax", LSRD_NL_itermax)
  ! LSRD_GMRES_itermin = 20

  ! Level Set massfix solver
  ! call rread("Mass_NL_tol", Mass_NL_tol)
  ! call iread("Mass_NL_itermax", Mass_NL_itermax)

  ! Wind turbine rotation
  ! call rread("Rotational_Vel", maxthetd)

  ! Build time integration parameters
  if ((rhoinf < 0.0d0) .or. (rhoinf > 1.0d0)) then ! backward Euler
    almi = 1.0d0
    alfi = 1.0d0
    gami = 1.0d0
    beti = 1.0d0
  else
    almi = (3.0d0 - rhoinf)/(1.0d0 + rhoinf)/2.0d0
    alfi = 1.0d0/(1.0d0 + rhoinf)

    gami = 0.5d0 + almi - alfi
    beti = 0.25d0*(1.0d0 + almi - alfi)*(1.0d0 + almi - alfi)
  end if
  ogam = 1.0d0/gami
  mgam = gami - 1.0d0

  if (ismaster) write (*, "(a20,x,' = ',x,ES12.4)") "almi", almi
  if (ismaster) write (*, "(a20,x,' = ',x,ES12.4)") "alfi", alfi
  if (ismaster) write (*, "(a20,x,' = ',x,ES12.4)") "gami", gami
  if (ismaster) write (*, "(a20,x,' = ',x,ES12.4)") "beti", beti
  Dtgl = 1.0d0/Delt

end subroutine getparam

!======================================================================
!
!======================================================================
subroutine rread(vname, val)
  use mpi
  implicit none

  character(len=*), intent(in)  :: vname
  real(8), intent(out) :: val
  character(len=20) :: buf

  call cread(vname, buf)
  read (buf, *) val

  if (ismaster) write (*, "(a20,x,' = ',x,ES12.4)") vname, val
end subroutine rread

!======================================================================
!
!======================================================================
subroutine iread(vname, val)
  use mpi
  implicit none

  character(len=*), intent(in)  :: vname
  integer, intent(out) :: val
  character(len=20) :: buf

  call cread(vname, buf)
  read (buf, *) val

  if (ismaster) write (*, "(a20,x,' = ',x,I12)") vname, val
end subroutine iread

!======================================================================
!
!======================================================================
subroutine lread(vname, val)
  use mpi
  implicit none

  character(len=*), intent(in)  :: vname
  logical, intent(out) :: val

  integer :: tmp

  call iread(vname, tmp)
  val = tmp /= 0

  if (ismaster) write (*, "(a20,x,' = ',x,I12)") vname, tmp
end subroutine lread

!======================================================================
!
!======================================================================
subroutine iread3(vname, val)
  use mpi
  implicit none

  character(len=*), intent(in)  :: vname
  integer, intent(out) :: val(3)
  character(len=20) :: buf

  call cread(vname, buf)
  read (buf, *) val(1:3)

  if (ismaster) write (*, "(a20,x,' = ',x,3I4)") vname, val
end subroutine iread3

!======================================================================
! subroutine to read in a variable name, check with input file
! and read in the corresponding value.
!======================================================================
subroutine cread(vname, val)
  use mpi
  implicit none

  character(len=*), intent(in)  :: vname
  character(len=20), intent(out) :: val
  character(len=40) :: buf, buf1
  integer :: paramf, ios, found, pos

  found = 0
  paramf = 45
  open (paramf, file='param.dat', status='old', iostat=ios)

  do while ((ios == 0) .and. (found == 0))
    read (paramf, '(A)', iostat=ios) buf

    ! find the position of character "="
    pos = scan(buf, '=')

    ! string before "=" will be checked with input
    buf1 = buf(1:pos - 1)

    ! string after "=" is the output value
    val = buf(pos + 1:)

    if (buf1 == vname) found = 1
  end do

  close (paramf)

  if (found == 0) then
    write (*, *) "Could not find ", vname, " in param.dat."
    call MPI_ABORT(MPI_COMM_WORLD, 911, mpi_err)
!!!    val=''
  end if
end subroutine cread
