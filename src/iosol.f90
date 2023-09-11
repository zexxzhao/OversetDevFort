!======================================================================
! Output aerodynamic torque for HAWT
!======================================================================
subroutine writeTQ(istep)
  use commonvars
  use mpi
  implicit none

  integer, intent(in) :: istep
  integer :: ifile

  character(len=30) :: fname
  character(len=10) :: cname

  ifile = 15
  fname = 'torque'//trim(cname(istep))
  open (ifile, file=fname, status='replace')
  write (ifile, "(4ES17.8)") time, theta, torque1, torque2
  write (ifile, "(4ES17.8)") torque3, torque4, moment1, moment2
  write (ifile, "(4ES17.8)") force_trac(1:3), force_cons(1:3)
  close (ifile)
end subroutine writeTQ

!======================================================================
!
!======================================================================
subroutine writeRB(istep)

  use aAdjKeep
  use commonvars
  use mpi

  implicit none

  integer, intent(in) :: istep

  integer :: mfile, dfile, vfile, rfile, ffile

  character(len=30) :: fname
  character(len=10) :: cname

  if (ismaster) then

    write (*, *) "Open files"
    mfile = 10
    fname = 'mass'//trim(cname(istep))//'.dat'
    open (mfile, file=fname, status='replace')

    dfile = 11
    fname = 'disp'//trim(cname(istep))//'.dat'
    open (dfile, file=fname, status='replace')

    vfile = 12
    fname = 'vel'//trim(cname(istep))//'.dat'
    open (vfile, file=fname, status='replace')

    rfile = 13
    fname = 'rot'//trim(cname(istep))//'.dat'
    open (rfile, file=fname, status='replace')

    ffile = 14
    fname = 'force'//trim(cname(istep))//'.dat'
    open (ffile, file=fname, status='replace')

    ! write (dfile, "(4ES14.6)") time, dbn1
    ! write (vfile, "(7ES14.6)") time, vbn1, wbn1
    ! write (rfile, "(10ES14.6)") time, Rn1
    ! write (ffile, "(7ES14.6)") time, Fb, Mb

    close (mfile)
    close (dfile)
    close (vfile)
    close (rfile)
    close (ffile)

  end if

end subroutine writeRB

!======================================================================
!
!======================================================================
subroutine readStep(Rstep)

  use aAdjKeep
  use commonvars
  use mpi

  implicit none

  integer :: solf, i, j, n, ierr, Rstep
  real(8) :: xi(3), xi1, xi3, rr

  character(len=30) :: fname
  character(len=10) :: cname

  logical :: exists

  ! Read step file
  Rstep = 0

!  if (ismaster) then
  solf = 15
  open (solf, file='step.dat', status='old', iostat=ierr)
  if (ierr == 0) then
    read (solf, *) Rstep, x_inflow
    close (solf)
  end if
!  end if

!  ! Communicate step
!  if (numnodes.gt.1) then
!    call MPI_BCAST(Rstep, 1, MPI_INTEGER, mpi_master, &
!                   MPI_COMM_WORLD, mpi_err)
!  end if

end subroutine readStep

!======================================================================
! Read the solution
!======================================================================
subroutine readSol(Rstep)
  use aAdjKeep
  use commonvars
  use mpi
  implicit none

  integer :: solf, i, j, n, ierr, Rstep
  real(8) :: xi(3), xi1, xi3, rr

  character(len=30) :: fname
  character(len=10) :: cname

  logical :: exists

  solf = 15
  if (ismaster) write (*, *) "Reading step:", Rstep
  fname = trim('restart'//cname(Rstep))//cname(myid + 1)
  open (solf, file=fname, status='old')

  read (solf, *) Rstep, time

  ! read(solf,*) (vbn0(j), j=1,3)
  ! read(solf,*) (dbn0(j), j=1,3)
  ! read(solf,*) (wbn0(j), j=1,3)

  do i = 1, NNODE
    read (solf, *) (dgold(i, j), j=1, 3)
  end do
  do i = 1, NNODE
    read (solf, *) (ugold(i, j), j=1, 3)
  end do
  do i = 1, NNODE
    read (solf, *) (ugmold(i, j), j=1, 3)
  end do
  do i = 1, NNODE
    read (solf, *) (acgold(i, j), j=1, 3)
  end do
  do i = 1, NNODE
    read (solf, *) (acgmold(i, j), j=1, 3)
  end do
  do i = 1, NNODE
    read (solf, *) phigold(i)
  end do
  do i = 1, NNODE
    read (solf, *) rphigold(i)
  end do
  do i = 1, NNODE
    read (solf, *) pgold(i)
  end do
  do i = 1, NNODE
    read (solf, *) Tgold(i)
  end do
  do i = 1, NNODE
    read (solf, *) rTgold(i)
  end do

!  read(solf,*) (vbn0(j),  j=1,3)
!  read(solf,*) (dbn0(j),  j=1,3)
!  read(solf,*) (wbn0(j),  j=1,3)
!  read(solf,*) (Rn0(1,j), j=1,3)
!  read(solf,*) (Rn0(2,j), j=1,3)
!  read(solf,*) (Rn0(3,j), j=1,3)
!  read(solf,*)  Mass_init

!  read(solf,*) thetaOld, thetdOld, theddOld
!  thetaOld = 0.0d0; thetdOld = 0.0d0; theddOld = 0.0d0

!  read(solf,*) i
!  close(solf)

!  if (i /= Rstep)  write(*,*) "ERROR IN READING SOLUTION!!"

end subroutine readSol

!======================================================================
! Generate initial condition
!======================================================================
subroutine generateIC()

  use aAdjKeep
  use commonvars
  use mpi
  use defs_shell

  implicit none

  integer :: i, j, k, n, iel, b
  real(8) :: xi(3), xi1, xi3, rr
  real(8) :: utmp(3), ptmp, rhotmp
  real(8) :: elem_h, rtmp
  real(8) :: ug_rand(NNODE, 3)
  if (ismaster) write (*, *) "Generating initial condition"

!  call init_random_seed()
  ! call random_number(ug_rand)
  time = 0.0d0
  Mass_init = -1.0d0

  ! vbn0(1) = 0.0d0
  ! vbn0(2) = 0.0d0
  ! vbn0(3) = 0.0d0
  ! dbn0 = 0.0d0
  ! wbn0(1) = 0.0d0
  ! wbn0(2) = 0.0d0
  ! wbn0(3) = 0.0d0
  ! Rn0 = 0.0d0
  ! Rn0(1, 1) = 1.0d0
  ! Rn0(2, 2) = 1.0d0
  ! Rn0(3, 3) = 1.0d0

!!$  do i = 1, NNODE
!!$    call getWave(ugold(i,:), phigold(i), xg(i,:), 0.0d0)
!!$  end do
!!$
!!$  do i = 1, NBOUND
!!$    if (bound(i)%FACE_ID >= 7) then
!!$      do j = 1, bound(i)%NNODE
!!$        n = bound(i)%BNODES(j)
!!$        ugold(n,1) = 0d0
!!$      end do
!!$    end if
!!$  end do

  ugold = 0.0d0
!  call getinflow(0,0)
!  ugold(:,1) = BCugValu(1,1)
!  call getinflow(0)
!  call setBCs_CFD_init(0)

  ! do i = 1, NBOUND
  !   if (bound(i)%FACE_ID >= 21 .and. bound(i)%FACE_ID <= 23) then
  !     do j = 1, bound(i)%NNODE
  !       n = bound(i)%BNODES(j)
  !       ugold(n, 1) = 0.0d0
  !     end do
  !   end if
  ! end do

  acgold = 0.0d0
  pgold = 0.0d0
  rphigold = 0.0d0
  ! elem_h = 0.0005d0
  phigold(:) = 0d0
  phig(:) = 0d0
  ! ugold(:, 1) = 1d0
  ! ugold(:, 2) = 0d0
  ! ugold(:, 3) = 0d0
  do i = 1, NNODE
    if(NODEID(i) == ELM_ID_SET(2)) then
      ugold(i, 1) = 1d0
    endif
  enddo
  do b = 1,NBOUND
    if(bound(b)%FACE_ID /= 5) cycle
    do i = 1,bound(b)%NNODE
      n = bound(b)%BNODES(i)
      ugold(n, 1) = 0d0
    end do
  enddo

  ! call commu(ugold, 3, "out")

  !call getinflow(0,0)
  !call setBCs_CFD_init1(0)
  rTgold(:) = 0d0

  dgold = 0.0d0
  ugmold = 0.0d0
  acgmold = 0.0d0

  ! thetaOld = 0.0d0
  ! thetdOld = 0.0d0
  ! theddOld = 0.0d0

  ! vbn1 = vbn0
  ! dbn1 = dbn0
  ! wbn1 = wbn0
  ! Rn1 = Rn0

  ug = ugold
  acg = acgold
  pg = pgold

  phig = phigold
  rphig = rphigold
  Tg = Tgold
  rTg = rTgold

  dg = dgold
  ugm = ugmold
  acgm = acgmold

  theta = thetaOld
  thetd = thetdOld
  thedd = theddOld

  ! if (solshell) then
! !    SH%TSP%dshOld = 0.0d0
  !   !   SH%TSP%ushOld = 0.0d0
  !   !  SH%TSP%ashOld = 0.0d0

  !   SH%NRB%dshOld = 0.0d0
  !   SH%NRB%ushOld = 0.0d0
  !   SH%NRB%ashOld = 0.0d0

  !   SH%FEM%dshOld = 0.0d0
  !   SH%FEM%ushOld = 0.0d0
  !   SH%FEM%ashOld = 0.0d0

! !    SH%TSP%dsh = SH%TSP%dshOld
  !   !   SH%TSP%ush = SH%TSP%ushOld
  !   !  SH%TSP%ash = SH%TSP%ashOld

  !   SH%NRB%dsh = SH%NRB%dshOld
  !   SH%NRB%ush = SH%NRB%ushOld
  !   SH%NRB%ash = SH%NRB%ashOld

  !   SH%FEM%dsh = SH%FEM%dshOld
  !   SH%FEM%ush = SH%FEM%ushOld
  !   SH%FEM%ash = SH%FEM%ashOld
  ! end if
end subroutine generateIC

subroutine writeVelocity(istep)

  use aAdjKeep
  use commonvars
  use mpi

  implicit none

  integer :: solf, i, j, istep
  character(len=30) :: fname
  character(len=10) :: cname

  ! Output results
  solf = 98
  fname = trim('velocity'//cname(istep))//cname(myid + 1)
  open (solf, file=fname, status='replace')

  write (solf, *) istep, time

  do i = 1, NNODE
    write (solf, *) (ug(i, j), j=1, 3), phig(i)
  end do

  close (solf)
end subroutine writeVelocity

!======================================================================
! Output solution
!======================================================================
subroutine writeSol(istep)

  use aAdjKeep
  use commonvars
  use mpi

  implicit none

  integer :: solf, i, j, istep
  character(len=30) :: fname
  character(len=10) :: cname

  ! Output results
  solf = 98
  fname = trim('restart'//cname(istep))//cname(myid + 1)
  open (solf, file=fname, status='replace')

  write (solf, *) istep, time

!  write(solf,*) (vbn1(j), j=1,3)
!  write(solf,*) (dbn1(j), j=1,3)
!  write(solf,*) (wbn1(j), j=1,3)

  do i = 1, NNODE
    write (solf, *) (dg(i, j), j=1, 3)
  end do
  do i = 1, NNODE
    write (solf, *) (ug(i, j), j=1, 3)
  end do
  do i = 1, NNODE
    write (solf, *) (ugm(i, j), j=1, 3)
  end do
  do i = 1, NNODE
    write (solf, *) (acg(i, j), j=1, 3)
  end do
  do i = 1, NNODE
    write (solf, *) (acgm(i, j), j=1, 3)
  end do
  do i = 1, NNODE
    write (solf, *) phig(i)
  end do
  do i = 1, NNODE
    write (solf, *) rphig(i)
  end do
  do i = 1, NNODE
    write (solf, *) pg(i)
  end do
  do i = 1, NNODE
    write (solf, *) Tgold(i)
  end do
  do i = 1, NNODE
    write (solf, *) rTgold(i)
  end do

  ! write (solf, *) (vbn1(j), j=1, 3)
  ! write (solf, *) (dbn1(j), j=1, 3)
  ! write (solf, *) (wbn1(j), j=1, 3)
  ! write (solf, *) (Rn1(1, j), j=1, 3)
  ! write (solf, *) (Rn1(2, j), j=1, 3)
  ! write (solf, *) (Rn1(3, j), j=1, 3)
  ! write (solf, *) Mass_init

  ! write (solf, *) theta, thetd, thedd

  ! write (solf, *) istep, Delt

  close (solf)

  if (numnodes > 1) call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

  if (ismaster) then
    solf = 88
    open (solf, file='step.dat', status='replace')
    write (solf, *) istep, x_inflow
    close (solf)
  end if

end subroutine writeSol

!======================================================================
! Read shell solution
!======================================================================
subroutine readShellSol(istep, SH)
  use mpi
  use defs_shell

  implicit none

  type(shell_bld), intent(inout) :: SH
  integer, intent(in)    :: istep

  integer :: solf, i, j
  character(len=30) :: fname, iname(2)

  solf = 97

  ! T-Spline
!  write(iname(1),'(I30)') istep
!  write(iname(2),'(I30)') myid + 21

!  fname = 'sh.rest.tsp.'//trim(adjustl(iname(1)))//'.'&
!                        //trim(adjustl(iname(2)))
!  open(solf, file=fname, status='old')

!  do i = 1, SH%TSP%NNODE
!    read(solf,*) (SH%TSP%dshOld(i,j), j = 1, 3)
!  end do
!  do i = 1, SH%TSP%NNODE
!    read(solf,*) (SH%TSP%ushOld(i,j), j = 1, 3)
!  end do
!  do i = 1, SH%TSP%NNODE
!    read(solf,*) (SH%TSP%ashOld(i,j), j = 1, 3)
!  end do

!  close(solf)

  ! NURBS
  write (iname(1), '(I30)') istep
  write (iname(2), '(I30)') myid + 21

  fname = 'sh.rest.nrb.'//trim(adjustl(iname(1)))//'.' &
          //trim(adjustl(iname(2)))
  open (solf, file=fname, status='old')

  do i = 1, SH%NRB%NNODE
    read (solf, *) (SH%NRB%dshOld(i, j), j=1, 3)
  end do
  do i = 1, SH%NRB%NNODE
    read (solf, *) (SH%NRB%ushOld(i, j), j=1, 3)
  end do
  do i = 1, SH%NRB%NNODE
    read (solf, *) (SH%NRB%ashOld(i, j), j=1, 3)
  end do

  close (solf)

  ! FEM
  fname = 'sh.rest.fem.'//trim(adjustl(iname(1)))//'.' &
          //trim(adjustl(iname(2)))
  open (solf, file=fname, status='old')

  do i = 1, SH%FEM%NNODE
    read (solf, *) (SH%FEM%dshOld(i, j), j=1, 3)
  end do
  do i = 1, SH%FEM%NNODE
    read (solf, *) (SH%FEM%ushOld(i, j), j=1, 3)
  end do
  do i = 1, SH%FEM%NNODE
    read (solf, *) (SH%FEM%ashOld(i, j), j=1, 3)
  end do

  close (solf)

end subroutine readShellSol

!======================================================================
! Output shell solution
!======================================================================
subroutine writeShellSol(istep, SH)
  use commonvars
  use mpi
  use defs_shell

  implicit none

  type(shell_bld), intent(in) :: SH
  integer, intent(in) :: istep

  integer :: solf, i, j
  character(len=30) :: fname, iname(2)

!!!  if (ismaster) then

  solf = 98

  ! T-Spline
!    write(iname(1),'(I30)') istep
!    write(iname(2),'(I30)') myid + 21

!    fname = 'sh.rest.tsp.'//trim(adjustl(iname(1)))//'.'&
!                          //trim(adjustl(iname(2)))
!    open(solf, file=fname, status='replace')

!    do i = 1, SH%TSP%NNODE
!      write(solf,*) (SH%TSP%dsh(i,j), j = 1, 3)
!    end do
!    do i = 1, SH%TSP%NNODE
!      write(solf,*) (SH%TSP%ush(i,j), j = 1, 3)
!    end do
!    do i = 1, SH%TSP%NNODE
!      write(solf,*) (SH%TSP%ash(i,j), j = 1, 3)
!    end do

!    write(solf,'(60("="))')
  !   write(solf,*) istep, time, Delt
  !   write(solf,*) theta, thetd, thedd
  !   write(solf,*) SH%TSP%dsh(SH%TSP%TipLoc,3), SH%Tq2
  !   close(solf)

  ! NURBS
  write (iname(1), '(I30)') istep
  write (iname(2), '(I30)') myid + 21

  fname = 'sh.rest.nrb.'//trim(adjustl(iname(1)))//'.' &
          //trim(adjustl(iname(2)))
  open (solf, file=fname, status='replace')

  do i = 1, SH%NRB%NNODE
    write (solf, *) (SH%NRB%dsh(i, j), j=1, 3)
  end do
  do i = 1, SH%NRB%NNODE
    write (solf, *) (SH%NRB%ush(i, j), j=1, 3)
  end do
  do i = 1, SH%NRB%NNODE
    write (solf, *) (SH%NRB%ash(i, j), j=1, 3)
  end do

  write (solf, '(60("="))')
  write (solf, *) istep, time, Delt
  write (solf, *) theta, thetd, thedd
  write (solf, *) SH%NRB%dsh(SH%NRB%TipLoc, 3), SH%Tq2
  close (solf)
  ! FEM
  fname = 'sh.rest.fem.'//trim(adjustl(iname(1)))//'.' &
          //trim(adjustl(iname(2)))
  open (solf, file=fname, status='replace')

  do i = 1, SH%FEM%NNODE
    write (solf, *) (SH%FEM%dsh(i, j), j=1, 3)
  end do
  do i = 1, SH%FEM%NNODE
    write (solf, *) (SH%FEM%ush(i, j), j=1, 3)
  end do
  do i = 1, SH%FEM%NNODE
    write (solf, *) (SH%FEM%ash(i, j), j=1, 3)
  end do

  write (solf, '(60("="))')
  write (solf, *) istep, time, Delt
  write (solf, *) theta, thetd, thedd
  write (solf, *) SH%FEM%dsh(SH%FEM%TipLoc, 3), SH%Tq1
  close (solf)

!!!  end if

end subroutine writeShellSol

!======================================================================
! Output tip displacement
!======================================================================
subroutine writeTip(istep, SH)
  use commonvars
  use mpi
  use defs_shell
  implicit none

  type(shell_bld), intent(in) :: SH
  integer, intent(in) :: istep
  integer :: ifile

  character(len=30) :: fname, iname(2)

  write (iname(1), '(I30)') istep
  write (iname(2), '(I30)') myid + 21

  ifile = 15
  fname = 'tipdispLeading.'//trim(adjustl(iname(1)))//'.' &
          //trim(adjustl(iname(2)))
  open (ifile, file=fname, status='replace')
  write (ifile, "(4ES17.8)") time, theta, SH%NRB%dsh(SH%NRB%TipLoc, 3), SH%Tq2
  write (ifile, "(4ES17.8)") time, theta, SH%FEM%dsh(SH%FEM%TipLoc, 3), SH%Tq1
  close (ifile)

  ifile = 15
  fname = 'tipdispTrailing.'//trim(adjustl(iname(1)))//'.' &
          //trim(adjustl(iname(2)))
  open (ifile, file=fname, status='replace')
  write (ifile, "(4ES17.8)") time, theta, SH%NRB%dsh(SH%NRB%TipLocTr, 3), SH%Tq2
  write (ifile, "(4ES17.8)") time, theta, SH%FEM%dsh(SH%FEM%TipLocTr, 3), SH%Tq1
  close (ifile)
end subroutine writeTip

!======================================================================
! Read shell solution
!======================================================================
subroutine readShellSol_NM(istep, NM)
  use defs_shell
  implicit none

  type(shell_nmb), intent(inout) :: NM
  integer, intent(in)    :: istep

  integer :: solf, i, j, bb
  character(len=30) :: fname, iname

  solf = 97

  write (iname, '(I30)') istep

  fname = 'nm.rest.fem.'//trim(adjustl(iname))
  open (solf, file=fname, status='old')

  do bb = 1, 2
    do i = 1, NM%FEM(bb)%NNODE
      read (solf, *) (NM%FEM(bb)%dshOld(i, j), j=1, 3)
    end do
    do i = 1, NM%FEM(bb)%NNODE
      read (solf, *) (NM%FEM(bb)%ushOld(i, j), j=1, 3)
    end do
    do i = 1, NM%FEM(bb)%NNODE
      read (solf, *) (NM%FEM(bb)%ashOld(i, j), j=1, 3)
    end do
  end do

  close (solf)
end subroutine readShellSol_NM

!======================================================================
! Output shell solution
!======================================================================
subroutine writeShellSol_NM(istep, NM)
  use commonvars
  use defs_shell
  implicit none

  type(shell_nmb), intent(in) :: NM
  integer, intent(in) :: istep

  integer :: solf, i, j, bb
  character(len=30) :: fname, iname

  solf = 98

  write (iname, '(I30)') istep

  fname = 'nm.rest.fem.'//trim(adjustl(iname))
  open (solf, file=fname, status='replace')

  do bb = 1, 2
    do i = 1, NM%FEM(bb)%NNODE
      write (solf, *) (NM%FEM(bb)%dsh(i, j), j=1, 3)
    end do
    do i = 1, NM%FEM(bb)%NNODE
      write (solf, *) (NM%FEM(bb)%ush(i, j), j=1, 3)
    end do
    do i = 1, NM%FEM(bb)%NNODE
      write (solf, *) (NM%FEM(bb)%ash(i, j), j=1, 3)
    end do
  end do

  write (solf, '(60("="))')
  write (solf, *) istep, time, Delt
  write (solf, *) theta, thetd, thedd

  close (solf)
end subroutine writeShellSol_NM

!======================================================================
! Output averaged solution
!======================================================================
subroutine writeAvgSol(istep)
  use aAdjKeep
  use commonvars
  use mpi
  implicit none

  integer, intent(in) :: istep
  integer :: solf, i, j
  character(len=30) :: fname
  character(len=10) :: cname

  ! Output results
  solf = 98
  fname = "avgsol"//cname(myid + 1)
  open (solf, file=fname, status='replace')

  write (solf, *) istep, Delt

  do i = 1, NNODE
    write (solf, *) (uavg(i, j)/dble(istep), j=1, 3)
  end do
  do i = 1, NNODE
    write (solf, *) pavg(i)/dble(istep)
  end do

  close (solf)
end subroutine writeAvgSol

!======================================================================
! Read in averaged solution
! istep: number of steps that have been averaged
!======================================================================
subroutine readAvgSol(istep)
  use aAdjKeep
  use commonvars
  use mpi
  implicit none

  integer, intent(out) :: istep
  integer :: solf, i, j, exts
  character(len=30) :: fname
  character(len=10) :: cname

  solf = 98
  fname = "avgsol"//cname(myid + 1)
  open (solf, file=fname, status='old', iostat=exts)

  if (exts == 0) then
    read (solf, *) istep
    do i = 1, NNODE
      read (solf, *) (uavg(i, j), j=1, 3)
    end do
    do i = 1, NNODE
      read (solf, *) pavg(i)
    end do
  else
    istep = 0
    uavg = 0.0d0
    pavg = 0.0d0
  end if

  close (solf)

  uavg = uavg*dble(istep)
  pavg = pavg*dble(istep)
end subroutine readAvgSol

!!$!======================================================================
!!$!
!!$!======================================================================
!!$subroutine writeVTK(istep)
!!$
!!$  use aAdjKeep
!!$  use commonvars
!!$  use mpi
!!$
!!$  implicit none
!!$
!!$  integer :: nshl
!!$
!!$  character(len=30) :: fname
!!$  character(len=10) :: cname
!!$  integer :: ifil, i, j, istep
!!$
!!$  ifil = 99
!!$  fname = trim('solution'//cname(istep))//cname(myid+1)
!!$
!!$  open(ifil, file=trim(fname)//'.vtk', status='replace', &
!!$             form='formatted')
!!$
!!$  write(ifil,'(a)') '# vtk DataFile Version 3.0'
!!$  write(ifil,'(a)') 'vtk output'
!!$  write(ifil,'(a)') 'ASCII'
!!$  write(ifil,'(a)') 'DATASET UNSTRUCTURED_GRID'
!!$
!!$  write(ifil,'(a,x,I8,x,a)') 'POINTS ',NNODE, 'float'
!!$
!!$  do i = 1, NNODE
!!$    write(ifil,'(3E17.8)') (real(xg(i,j) + dg(i,j)), j = 1, NSD)
!!$  end do
!!$
!!$  write(ifil,'(a,x,I8,x,I8)') 'CELLS ',NELEM, NELEM*5
!!$  do i = 1, NELEM
!!$    write(ifil,'(5I8)') NSHL, (IEN(i,j)-1, j = 1, NSHL)
!!$  end do
!!$
!!$  write(ifil,'(a,x,I8)') 'CELL_TYPES',NELEM
!!$  do i = 1, NELEM
!!$    write(ifil,'(I8)') 10
!!$  enddo
!!$
!!$  write(ifil,'(a,x,I8)') 'POINT_DATA', NNODE
!!$
!!$  write(ifil,'(a)') 'VECTORS u float'
!!$  do i = 1, NNODE
!!$    write(ifil,'(3E17.8)') real(ug(i,1)),real(ug(i,2)),real(ug(i,3))
!!$  enddo
!!$
!!$  write(ifil,'(a)') 'SCALARS p float'
!!$  write(ifil,'(a)') 'LOOKUP_TABLE default'
!!$  do i = 1, NNODE
!!$    write(ifil,'(E17.8)') real(pg(i))
!!$  enddo
!!$
!!$  write(ifil,'(a)') 'SCALARS phi float'
!!$  write(ifil,'(a)') 'LOOKUP_TABLE default'
!!$  do i = 1, NNODE
!!$    write(ifil,'(E17.8)') real(phig(i))
!!$  enddo
!!$
!!$
!!$  close(ifil)
!!$
!!$end subroutine writeVTK
!!$
!!$
!!$
!!$!======================================================================
!!$!
!!$!======================================================================
!!$subroutine writeVTKXML(istep)
!!$
!!$  use aAdjKeep
!!$  use commonvars
!!$  use mpi
!!$
!!$  implicit none
!!$
!!$  integer :: nshl
!!$
!!$  character(len=30) :: fname
!!$  character(len=10) :: cname
!!$  character(len=10) :: cname2
!!$  integer :: ifil, i,j,istep
!!$
!!$  ! Open file
!!$  ifil = 99
!!$  fname = trim('solution' // cname(myid+1))  // cname (istep)
!!$
!!$  open(ifil, file=trim(fname)//'.vtu', status='replace', &
!!$             form='formatted')
!!$
!!$  ! Header
!!$  write(ifil,'(a)') '<?xml version="1.0"?>'
!!$  write(ifil,'(a)') '<VTKFile type="UnstructuredGrid" '// &
!!$                    'version="0.1" byte_order="LittleEndian">'
!!$  write(ifil,'(a)') '   <UnstructuredGrid>  '
!!$  write(ifil,'(a)') '     <Piece ' &
!!$                    //'NumberOfPoints="'//trim(cname2(NNODE))//'" ' &
!!$                    //'NumberOfCells="'//trim(cname2(NELEM))//'">  '
!!$
!!$  ! Node coordinates
!!$  write(ifil,'(a)') '<Points>'
!!$  write(ifil,'(a)') '  <DataArray NumberOfComponents="3" '
!!$  write(ifil,'(a)') '       type="Float32" format="ascii">'
!!$  do i = 1, NNODE
!!$    write(ifil,'(3E17.8)') (real(xg(i,j) + dg(i,j)), j = 1, NSD)
!!$  end do
!!$  write(ifil,'(a)') '  </DataArray>'
!!$  write(ifil,'(a)') '</Points>'
!!$
!!$  ! Element connectivity and type
!!$  write(ifil,'(a)') '<Cells>'
!!$  write(ifil,'(a)') '  <DataArray type="Int32" '
!!$  write(ifil,'(a)') '       Name="connectivity" format="ascii">'
!!$  do i = 1, NELEM
!!$    write(ifil,'(4I7)') (IEN(i,j)-1, j = 1, NSHL)
!!$  end do
!!$  write(ifil,'(a)') '  </DataArray>'
!!$
!!$  write(ifil,'(a)') '  <DataArray type="Int32" '
!!$  write(ifil,'(a)') '       Name="offsets" format="ascii">'
!!$  do i = 1, NELEM
!!$    write(ifil,'(I12)') NSHL*i
!!$  end do
!!$  write(ifil,'(a)') '  </DataArray>'
!!$
!!$  write(ifil,'(a)') '  <DataArray type="UInt8" '
!!$  write(ifil,'(a)') '       Name="types" format="ascii"> '
!!$  do i = 1, NELEM
!!$    write(ifil,'(I2)') 10
!!$  enddo
!!$  write(ifil,'(a)') '  </DataArray>'
!!$
!!$  write(ifil,'(a)') '</Cells>'
!!$
!!$  ! Data
!!$  write(ifil,'(a)') '  <PointData Scalars="scalars">'
!!$
!!$  write(ifil,'(a)') '  <DataArray type="Float32" Name="velocity" '
!!$  write(ifil,'(a)') '       NumberOfComponents="3" format="ascii">'
!!$  do i = 1, NNODE
!!$    write(ifil,'(3E17.8)') real(ug(i,1)),real(ug(i,2)),real(ug(i,3))
!!$  enddo
!!$  write(ifil,'(a)') '  </DataArray>'
!!$
!!$  write(ifil,'(a)') '  <DataArray type="Float32" Name="pressure" '
!!$  write(ifil,'(a)') '       format="ascii">'
!!$  do i = 1, NNODE
!!$    write(ifil,'(E17.8)') real(pg(i))
!!$  enddo
!!$  write(ifil,'(a)') '  </DataArray>'
!!$
!!$  write(ifil,'(a)') '  <DataArray type="Float32" Name="level set" '
!!$  write(ifil,'(a)') '       format="ascii">'
!!$  do i = 1, NNODE
!!$    write(ifil,'(E17.8)') real(phig(i))
!!$  enddo
!!$  write(ifil,'(a)') '  </DataArray>'
!!$
!!$
!!$  write(ifil,'(a)') '  </PointData>'
!!$  write(ifil,'(a)') '  <CellData Scalars="scalars">'
!!$
!!$  write(ifil,'(a)') '  </CellData>'
!!$
!!$  ! Footer
!!$  write(ifil,'(a)') '    </Piece>'
!!$  write(ifil,'(a)') '  </UnstructuredGrid>'
!!$  write(ifil,'(a)') '</VTKFile>'
!!$
!!$  close(ifil)
!!$  if (numnodes.gt.1) call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!!$
!!$  ! Gather file
!!$  if (ismaster) then
!!$    fname = trim('solution' //cname (istep))//'.pvtu'
!!$
!!$    open(ifil, file=trim(fname),status='replace', form='formatted')
!!$
!!$    write(ifil,'(a)') '<?xml version="1.0"?>'
!!$    write(ifil,'(a)') '<VTKFile type="PUnstructuredGrid" ' &
!!$                      //'version="0.1" byte_order="LittleEndian">'
!!$    write(ifil,'(a)') '  <PUnstructuredGrid GhostLevel="0">'
!!$
!!$    write(ifil,'(a)') '    <PPoints>'
!!$    write(ifil,'(a)') '      <PDataArray type="Float32" ' &
!!$                      //'NumberOfComponents="3"/>'
!!$    write(ifil,'(a)') '    </PPoints>'
!!$
!!$    write(ifil,'(a)') '    <PPointData Scalars="scalars">'
!!$    write(ifil,'(a)') '<PDataArray type="Float32"  ' &
!!$                      //'Name="velocity" NumberOfComponents="3"/>'
!!$    write(ifil,'(a)') '<PDataArray type="Float32" ' &
!!$                      //'Name="pressure"/>'
!!$    write(ifil,'(a)') '<PDataArray type="Float32" ' &
!!$                      //'Name="level set"/>'
!!$    write(ifil,'(a)') '    </PPointData>'
!!$
!!$    write(ifil,'(a)') '    <PCellData Scalars="scalars">'
!!$    write(ifil,'(a)') '    </PCellData>'
!!$
!!$    do i=1, numnodes
!!$      fname=trim(trim('solution'//cname(i))//cname(istep))//'.vtu'
!!$      write(ifil,'(a)') '    <Piece Source="'//trim(fname)//'"/>'
!!$    enddo
!!$
!!$    write(ifil,'(a)') '  </PUnstructuredGrid>'
!!$    write(ifil,'(a)') '</VTKFile>'
!!$    close(ifil)
!!$
!!$  endif
!!$
!!$end subroutine writeVTKXML
