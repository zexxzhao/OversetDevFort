
subroutine setBCs_CFD_init1(istep)
  use aAdjKeep
  use commonvars
  use commonpars
  use mpi
  implicit none
  integer, intent(in)::istep
  integer :: b, i, j, k, n, d, dir, ptmp, iface, global_index
  real(8) :: u(NSD), phi

  IBC = 0
  ogam = 1.0d0/gami
  mgam = gami - 1.0d0

  do b = 1, NBOUND
    if (bound(b)%FACE_ID == 1) then
      do j = 1, bound(b)%NNODE
        k = bound(b)%BNODES(j)
        global_index = bound(b)%L2SNODE(j)
        do d = 1, NSD
          IBC(k, d) = 1
        end do
      end do
    end if
  end do

end subroutine setBCs_CFD_init1

subroutine setBCs_CFD_init(istep)
  use aAdjKeep
  use commonvars
  use commonpars
  use mpi
  implicit none
  integer, intent(in)::istep
  integer :: b, i, j, k, n, d, dir, ptmp, iface, global_index
  real(8) :: u(NSD), phi

  IBC = 0
  ogam = 1.0d0/gami
  mgam = gami - 1.0d0

  do b = 1, NBOUND
    if (bound(b)%FACE_ID == 1) then
      do j = 1, bound(b)%NNODE
        k = bound(b)%BNODES(j)
        global_index = bound(b)%L2SNODE(j)
!         IBC(k,7) = 1
!         pgold = inflow_pressure(global_index)
!         IBC(k,8) = 1
!         phigold = inflow_rho(global_index)
        do d = 1, NSD
          IBC(k, d) = 1
          ugold(k, d) = inflow_velocity(global_index, d)
        end do
      end do
    end if
  end do

end subroutine setBCs_CFD_init

!======================================================================
!
!======================================================================
subroutine setBCs_CFD()
  use commonvars
  use aAdjkeep
  use mpi

  implicit none

  integer :: b, i, j, k, n, d, dir, ptmp, iface, global_index

  real(8) :: u(NSD), phi

  IBC = 0
  ogam = 1.0d0/gami
  mgam = gami - 1.0d0

  do b = 1, NBOUND

    do d = 1, NSD
      if (BCugType(b, d) == 1) then
        do j = 1, bound(b)%NNODE
          k = bound(b)%BNODES(j)
          IBC(k, d) = 1
          ug(k, d) = BCugValu(b, d)
          acg(k, d) = ((ug(k, d) - ugold(k, d))*Dtgl &
                       + (gami - 1.0d0)*acgold(k, d))/gami

        end do
      endif
    end do

!      if (bound(b)%FACE_ID == 1) then
!        do j = 1, bound(b)%NNODE
!          k = bound(b)%BNODES(j)
!          IBC(k,8) = 1
!          phig(k) = -1.0d0/(gravity*Fr**2.0d0)*xg(k,3)+1.0d0+ 1.0d0/(2.0d0*gravity*Fr**2.0d0)  ! Linear gradient on the inflow
!          rphig(k) = ( (phig(k)-phigold(k))*Dtgl  &
!                        + (gami-1.0d0)*rphigold(k)  )/gami
!        end do
!      end if

!     if (bound(b)%FACE_ID == 21) then
!        do j = 1, bound(b)%NNODE
!          k = bound(b)%BNODES(j)
!          IBC(k,8) = 1
!          phig(k) = -1.0d0/(gravity*Fr**2.0d0)*xg(k,3)+1.0d0+ 1.0d0/(2.0d0*gravity*Fr**2.0d0)  ! Linear gradient on the inflow
!          rphig(k) = ( (phig(k)-phigold(k))*Dtgl  &
!                        + (gami-1.0d0)*rphigold(k)  )/gami
!        end do
!      end if

  end do  ! end loop: all faces

end subroutine setBCs_CFD

!======================================================================
!
!======================================================================

subroutine setBCs_overset(subdomain_flag)
  use commonvars
  use aAdjkeep
  use mpi

  implicit none

  integer, intent(in) :: subdomain_flag
  integer :: b, i, j, k, n, d, dir, ptmp, iface, global_index

  call setBCs_CFD()

  ! check if zero out subdomain one
  if(iand(subdomain_flag, 1) == 0) then
    do i = 1, NNODE
      if(NodeID(i) == ELM_ID_SET(1)) then
        IBC(i, 1:4) = 1
      endif
    enddo
  endif
  if(iand(subdomain_flag, 2) == 0) then
    do i = 1, NNODE
      if(NodeID(i) == ELM_ID_SET(2)) then
        IBC(i, 1:4) = 1
      endif
    enddo
  endif

end subroutine setBCs_overset

!======================================================================
!
!======================================================================
subroutine setBCs_NSVOF()
  use commonvars
  use aAdjkeep
  use mpi

  implicit none

  integer :: b, i, j, k, n, d, dir, ptmp, iface, global_index

  real(8) :: u(NSD), phi

  !-----------------
  ! IBC(:, 1:3) - velocity
  ! IBC(:, 4) - pressure
  ! IBC(:, 5) - VOF
  ! IBC(:, 6) - Tem
  IBC = 0
  ogam = 1.0d0/gami
  mgam = gami - 1.0d0

  ! if (.not. is_solid_node_assigned) then
  !   if(ismaster) write(*,*) "Generating BCs"
  !   call get_solid_nodes_private()
  !   is_solid_node_assigned = .true.
  ! end if
  ! IBC(:, 1) = IS_SOLID_NODE(:)
  ! IBC(:, 2) = IS_SOLID_NODE(:)
  ! IBC(:, 3) = IS_SOLID_NODE(:)
  ! IBC(:, 4) = IS_SOLID_NODE(:)
  ! IBC(:, 5) = IS_SOLID_NODE(:)
  do i = 1, NNODE
    if(NodeID(i) == 101) then
      IBC(i, 1:5) = 1
    else
      IBC(i, 1:5) = 0
    endif
  enddo

  do b = 1, NBOUND
    do d = 1, NSD
      if (BCugType(bound(b)%FACE_ID, d) == 1) then
        do j = 1, bound(b)%NNODE
          k = bound(b)%BNODES(j)
          IBC(k, d) = 1
          ug(k, d) = BCugValu(b, d)
          acg(k, d) = ((ug(k, d) - ugold(k, d))*Dtgl &
                       + (gami - 1.0d0)*acgold(k, d))/gami

        end do

      ! else if (BCugType(bound(b)%FACE_ID, d) == 3) then
      !   do j = 1, bound(b)%NNODE
      !     k = bound(b)%BNODES(j)
      !     IBC(k, d) = 1
      !     ug(k, d) = ugm(k, d)
      !     acg(k, d) = ((ug(k, d) - ugold(k, d))*Dtgl &
      !                  + (gami - 1.0d0)*acgold(k, d))/gami
      !   end do
      end if
    end do

    if (BCphigType(bound(b)%FACE_ID) == 1) then
      do j = 1, bound(b)%NNODE
        k = bound(b)%BNODES(j)
        IBC(k, 5) = 1
        phig(k) = BCphigValu(b) ! further dev
        rphig(k) = ((phig(k) - phigold(k))*Dtgl &
                    + (gami - 1.0d0)*rphigold(k))/gami
      end do
    end if
  end do  ! end loop: all faces

end subroutine setBCs_NSVOF
!======================================================================
!
!======================================================================
subroutine setBCs_Tem()
  use commonvars
  use aAdjkeep
  use mpi

  implicit none

  integer :: b, i, j, k, n, d, dir, ptmp, iface, global_index

  real(8) :: u(NSD), phi

  !-----------------
  ! IBC(:, 1:3) - velocity
  ! IBC(:, 4) - pressure
  ! IBC(:, 5) - VOF
  ! IBC(:, 6) - Tem
  ! IBC(:, :) = 1
  IBC(:, 6) = 0
  do b = 1, NBOUND
    if (BCTgType(bound(b)%FACE_ID) == 1) then
      do j = 1, bound(b)%NNODE
        k = bound(b)%BNODES(j)
        IBC(k, 6) = 1
        Tg(k) = BCTgValu(b)
        rTg(k) = ((Tg(k) - Tgold(k))*Dtgl &
                  + (gami - 1.0d0)*rTgold(k))/gami
      end do
    end if
  end do  ! end loop: all faces

end subroutine setBCs_Tem

!======================================================================
! Get the rotation matrices
!======================================================================
subroutine get_Rmat(theta, thetd, thedd, Rmat, Rdot, Rddt)
  implicit none

  real(8), intent(in) :: theta, thetd, thedd
  real(8), intent(out) :: Rmat(3, 3), Rdot(3, 3), Rddt(3, 3)

  Rmat = 0.0d0
  Rmat(2, 2) = cos(theta) - 1.0d0
  Rmat(2, 3) = -sin(theta)
  Rmat(3, 2) = sin(theta)
  Rmat(3, 3) = cos(theta) - 1.0d0

  Rdot = 0.0d0
  Rdot(2, 2) = -sin(theta)*thetd
  Rdot(2, 3) = -cos(theta)*thetd
  Rdot(3, 2) = cos(theta)*thetd
  Rdot(3, 3) = -sin(theta)*thetd

  Rddt = 0.0d0
  Rddt(2, 2) = -cos(theta)*thetd**2 - sin(theta)*thedd
  Rddt(2, 3) = sin(theta)*thetd**2 - cos(theta)*thedd
  Rddt(3, 2) = -sin(theta)*thetd**2 + cos(theta)*thedd
  Rddt(3, 3) = -cos(theta)*thetd**2 - sin(theta)*thedd
end subroutine get_Rmat

!======================================================================
! Get the rotation matrices
!======================================================================
subroutine get_Rmat_X(theta, thetd, thedd, Rmat, Rdot, Rddt)
  implicit none

  real(8), intent(in) :: theta, thetd, thedd
  real(8), intent(out) :: Rmat(3, 3), Rdot(3, 3), Rddt(3, 3)

  Rmat = 0.0d0
  Rmat(2, 2) = cos(theta) - 1.0d0
  Rmat(2, 3) = -sin(theta)
  Rmat(3, 2) = sin(theta)
  Rmat(3, 3) = cos(theta) - 1.0d0

  Rdot = 0.0d0
  Rdot(2, 2) = -sin(theta)*thetd
  Rdot(2, 3) = -cos(theta)*thetd
  Rdot(3, 2) = cos(theta)*thetd
  Rdot(3, 3) = -sin(theta)*thetd

  Rddt = 0.0d0
  Rddt(2, 2) = -cos(theta)*thetd**2 - sin(theta)*thedd
  Rddt(2, 3) = sin(theta)*thetd**2 - cos(theta)*thedd
  Rddt(3, 2) = -sin(theta)*thetd**2 + cos(theta)*thedd
  Rddt(3, 3) = -cos(theta)*thetd**2 - sin(theta)*thedd
end subroutine get_Rmat_X

!======================================================================
! Set the BC for the case of rotating HAWT
!======================================================================
subroutine setMeshBCs_hawt
  use commonvars
  use aAdjkeep
  implicit none

  integer :: b, i, j, k
  real(8) :: rot(NSD)

  do b = 1, NBOUND
    ! apply the same mesh BC for inlet/side/outlet/rotor
    do i = 1, bound(b)%NNODE
      j = bound(b)%BNODES(i)
      IBC(j, 4:6) = 1
    end do
  end do

  !put IBC on the whole mesh
!  IBC(:,4:6) = 1
end subroutine setMeshBCs_hawt

!======================================================================
! Set the BC for the case of rotating 5MW with tower
!======================================================================
subroutine setMeshBCs_tower
  use commonvars
  use aAdjkeep
  implicit none

  integer :: b, i, j, k

  ! for the fixed outer box domain (including the tower)
  do i = 1, NNODE
    if (NodeID(i) == 31) then
      IBC(i, 4:6) = 1
    end if
  end do

  ! for the disk and the blade
  do b = 1, NBOUND

    if (((bound(b)%FACE_ID >= 21) .and. (bound(b)%FACE_ID <= 29)) .or. &
        ((bound(b)%FACE_ID >= 41) .and. (bound(b)%FACE_ID <= 49)) .or. &
        (bound(b)%FACE_ID == 11)) then

      do i = 1, bound(b)%NNODE
        j = bound(b)%BNODES(i)
        IBC(j, 4:6) = 1
      end do

    end if
  end do
end subroutine setMeshBCs_tower

!======================================================================
!
!======================================================================
subroutine getWave(u, phi, x, t)
  use commonvars
  use aAdjkeep
  use mpi
  implicit none

  real(8) :: u(NSD), phi, x(NSD), t
  real(8) :: kxwt, kzh, amp, xy

  amp = wave_periode*(wave_amp/sinh(wave_length*water_depth))
  xy = dcos(wave_angle)*x(1) + dsin(wave_angle)*x(2)

  kxwt = wave_length*(xy - domain_left - Uin*t) - wave_periode*time
  kzh = wave_length*min(x(3) - domain_bottom, water_depth)

  u(1) = amp*dcosh(kzh)*dcos(kxwt)*dcos(wave_angle) + Uin
  u(2) = amp*dcosh(kzh)*dcos(kxwt)*dsin(wave_angle)
  u(3) = amp*dsinh(kzh)*dsin(kxwt)
  phi = wave_amp*dcos(kxwt) + water_level - x(3)
end subroutine getWave

!======================================================================
!
!======================================================================
subroutine setBCs()
  use commonvars
  use aAdjkeep
  use mpi

  implicit none

  integer :: b, i, j, k, n, dir, ptmp
  real(8) :: u(NSD), phi

  IBC = 0
  ogam = 1.0d0/gami
  mgam = gami - 1.0d0
  do b = 1, NBOUND
    if ((bound(b)%FACE_ID >= 7) .and. (BCType(7) == 1)) then ! Hull
      do j = 1, bound(b)%NNODE
        i = bound(b)%BNODES(j)
        ug(i, :) = ugm(i, :)
        ugold(i, :) = ugmold(i, :)
        acg(i, :) = acgm(i, :)
        acgold(i, :) = acgmold(i, :)

        IBC(i, 1:3) = 1
      end do
    end if

    ! Sides
    if (((bound(b)%FACE_ID == 2) .and. (BCType(2) == 1)) .or. &
        ((bound(b)%FACE_ID == 4) .and. (BCType(4) == 1))) then
      do j = 1, bound(b)%NNODE
        i = bound(b)%BNODES(j)
        call getWave(u, phi, xg(i, :), time)

        if (((bound(b)%FACE_ID == 2) .and. (u(2) > 0.0d0)) .or. &
            ((bound(b)%FACE_ID == 4) .and. (u(2) < 0.0d0))) then
          ug(i, :) = u
          phig(i) = phi

          acg(i, :) = ((ug(i, :) - ugold(i, :))*Dtgl &
                       + (gami - 1.0d0)*acgold(i, :))/gami

          rphig(i) = ((phig(i) - phigold(i))*Dtgl &
                      + (gami - 1.0d0)*rphigold(i))/gami

          IBC(i, 1:3) = 1
          IBC(i, 8) = 1
        end if
      end do
    end if

    ! inflow
    if ((bound(b)%FACE_ID == 5) .and. (BCType(5) == 1)) then
      do j = 1, bound(b)%NNODE
        i = bound(b)%BNODES(j)
        call getWave(ug(i, :), phig(i), xg(i, :), time)

        acg(i, :) = ((ug(i, :) - ugold(i, :))*Dtgl &
                     + (gami - 1.0d0)*acgold(i, :))/gami

        rphig(i) = ((phig(i) - phigold(i))*Dtgl &
                    + (gami - 1.0d0)*rphigold(i))/gami

        IBC(i, 1:3) = 1
        IBC(i, 8) = 1
      end do
    end if

    ! Bottom
    if ((bound(b)%FACE_ID == 1) .and. (BCType(1) == 1)) then
      do j = 1, bound(b)%NNODE
        i = bound(b)%BNODES(j)
        ug(i, 3) = 0.0d0
        ugold(i, 3) = 0.0d0
        acg(i, 3) = 0.0d0
        acgold(i, 3) = 0.0d0

        IBC(i, 3) = 1
      end do
    end if

    ! Top
    if ((bound(b)%FACE_ID == 6) .and. (BCType(6) == 1)) then
      do j = 1, bound(b)%NNODE
        i = bound(b)%BNODES(j)
        ug(i, 2) = 0.0d0
        ugold(i, 2) = 0.0d0
        acg(i, 2) = 0.0d0
        acgold(i, 2) = 0.0d0

        IBC(i, 2) = 1
      end do
    end if

  end do

end subroutine setBCs

!======================================================================
!
!======================================================================
! subroutine setMeshBCs()
!   use commonvars
!   use aAdjkeep
! 
!   implicit none
! 
!   integer :: b, i, j, k
!   real(8) :: rot(NSD)
! 
!   do b = 1, NBOUND
! 
!     if (bound(b)%FACE_ID >= 7) then ! Object
!       do j = 1, bound(b)%NNODE
!         i = bound(b)%BNODES(j)
! 
!         rot = matmul(Rn1, (xg(i, :) - xcg))
!         dg(i, :) = dbn1 + rot - (xg(i, :) - xcg)
! 
!         acgm(i, :) = (1.0d0/beti)*(-(ugmold(i, :)/Delt) &
!                                    + ((dg(i, :) - dgold(i, :))/(Delt*Delt)) &
!                                    + (beti - 0.5d0)*acgmold(i, :))
! 
!         ugm(i, :) = ugmold(i, :) &
!                     + Delt*((1.0d0 - gami)*acgmold(i, :) + gami*acgm(i, :))
! 
!         IBC(i, 4:6) = 1
!       end do
!     end if
! 
!     ! Outflow/Inflow
!     if ((bound(b)%FACE_ID == 3) .or. (bound(b)%FACE_ID == 5)) then
!       do j = 1, bound(b)%NNODE
!         i = bound(b)%BNODES(j)
!         dg(i, 1:3) = 0.0d0
!         ugm(i, 1:3) = 0.0d0
!         acgm(i, 1:3) = 0.0d0
!         dgold(i, 1:3) = 0.0d0
!         ugmold(i, 1:3) = 0.0d0
!         acgmold(i, 1:3) = 0.0d0
! 
!         IBC(i, 4:6) = 1
!       end do
!     end if
! 
!     ! Sides
!     if ((bound(b)%FACE_ID == 2) .or. (bound(b)%FACE_ID == 4)) then
!       do j = 1, bound(b)%NNODE
!         i = bound(b)%BNODES(j)
!         dg(i, 1:3) = 0.0d0
!         ugm(i, 1:3) = 0.0d0
!         acgm(i, 1:3) = 0.0d0
!         dgold(i, 1:3) = 0.0d0
!         ugmold(i, 1:3) = 0.0d0
!         acgmold(i, 1:3) = 0.0d0
! 
!         IBC(i, 4:6) = 1
!       end do
!     end if
! 
!     ! Top/bottom
!     if ((bound(b)%FACE_ID == 1) .or. (bound(b)%FACE_ID == 6)) then
!       do j = 1, bound(b)%NNODE
!         i = bound(b)%BNODES(j)
!         dg(i, 1:3) = 0.0d0
!         ugm(i, 1:3) = 0.0d0
!         acgm(i, 1:3) = 0.0d0
!         dgold(i, 1:3) = 0.0d0
!         ugmold(i, 1:3) = 0.0d0
!         acgmold(i, 1:3) = 0.0d0
! 
!         IBC(i, 4:6) = 1
!       end do
!     end if
! 
!   end do
! 
! end subroutine setMeshBCs
