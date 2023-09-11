!======================================================================
! Main routine to call all the subroutines
!======================================================================
program NURBScode

  use aAdjKeep
  use mpi
  use commonvars
  use defs_shell

  implicit none

  integer :: i, j, k, ii, istep, Rstep, nn, dd, ibld, avgstepold, avgstep
  real(8) :: Rmat(3, 3), Rdot(3, 3), Rddt(3, 3), &
             RmatOld(3, 3), RdotOld(3, 3), RddtOld(3, 3), &
             utmp(3), umtmp(3), Rtmp(3, 3), Tacc, ForceTemp(3)
  real(8), allocatable :: dshalpha(:, :)
  real(8), allocatable :: NRmat(:, :, :), NRdot(:, :, :), NRddt(:, :, :), &
                          NRmatOld(:, :, :), NRdotOld(:, :, :), NRddtOld(:, :, :)

  integer :: aa, bb, iel
  real(8) :: xs(3, 11), xl(4, 3)
  integer :: telem(11)
  real(8) :: n_tmp(4), basis_tmp(4, 11), v_tmp
  real(8) :: Tlocal(11), Tglobal(11), Tl(4)

  character(len=30) :: fname, iname, cname

  ! Initialize MPI
  call MPI_INIT(mpi_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numnodes, mpi_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpi_err)
  ismaster = myid .eq. mpi_master

!A  Tacc = 2.0d0


!!!  solshell = myid.eq.0

  ! flag for non-matching computation
  nonmatch = .false.
  if (ismaster) write (*, *) "Get run parameters"
  call getparam()
  ! Read mesh and MPI-communication Data
  if (ismaster) write (*, *) "Read mesh and communication data"
  call input(myid + 1)
  if (numnodes > 1) call ctypes()

  ! Generate Sparse Structures
  if (ismaster) write (*, *) "Generating sparse structure"
  call genSparStruc()

  ! Allocate Matrices and Vectors
  if (ismaster) write (*, *) "Allocating matrices and vectors"
  call allocMatVec()

  ! ! Read in restart files
  ! call readStep(Rstep)
  ! ! Get initial condition
  ! if (Rstep == 0) then
  !   call generateIC()
  !   call writeSol(Rstep)
  ! else
  !   call readSol(Rstep)
  ! end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  xs(:, 1) = (/12.5, -56.6584, 44.64270/)
  xs(:, 2) = (/12.5, -58.3866, -39.8476/)
  xs(:, 3) = (/12.5, -54.7381, -36.1991/)
  xs(:, 4) = (/12.5, -42.0646, 28.70476/)
  xs(:, 5) = (/12.5, -38.0321, 32.73725/)
  xs(:, 6) = (/12.5, -23.8224, -51.5610/)
  xs(:, 7) = (/12.5, -21.7101, -55.5935/)
  xs(:, 8) = (/12.5, 46.45816, -51.1770/)
  xs(:, 9) = (/12.5, 48.18637, -55.5935/)
  xs(:, 10) = (/12.5, -31.3852, 51.93958/)
  xs(:, 11) = (/12.5, -24.9745, 51.93958/)
  xs(:, :) = xs(:, :) * 1d-3
  telem(:) = 0
  basis_tmp = 0.0d0
  do iel = 1,NELEM
    do i = 1, size(xs, 2)
      n_tmp(:) = 0d0
      do aa = 1, 4
        xl(aa, :) = xg(ien(iel, aa), :)
      enddo
      call voltet_abs(xs(:, i), xl(2, :), xl(3, :), xl(4, :), n_tmp(1)) 
      call voltet_abs(xl(1, :), xs(:, i), xl(3, :), xl(4, :), n_tmp(2)) 
      call voltet_abs(xl(1, :), xl(2, :), xs(:, i), xl(4, :), n_tmp(3)) 
      call voltet_abs(xl(1, :), xl(2, :), xl(3, :), xs(:, i), n_tmp(4)) 
      call voltet_abs(xl(1, :), xl(2, :), xl(3, :), xl(4, :), v_tmp) 
      n_tmp(:) = n_tmp(:) / v_tmp

      if(abs(sum(n_tmp(:)) - 1d0) < 1d-9) then
        telem(i) = iel
        basis_tmp(:, i) = n_tmp
        write(*,*) "Found", i, "@(", myid, iel, ")N=", n_tmp
      endif
    enddo
  enddo
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! find xs in elements
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !------------------------------------------
  ! Loop over time steps
  !------------------------------------------
  call readStep(Nstep)
  do istep = 0, Nstep, ifq

    Rstep = istep
    call readSol(Rstep)
    if (ismaster) then
      write (*, '(60("="))')
      write (*, "(a,x,I8,x,ES14.6)") "Time Step Number:", istep, time
      write (*, '(60("="))')
    end if
    ! call readStep(istep)
    Tlocal(:) = 0d0
    do i = 1, size(xs, 2)
        iel = telem(i)
        if(iel == 0) cycle
        do aa = 1, 4
          Tl(aa) = Tgold(ien(iel, aa))
        enddo
        Tlocal(i) = dot_product(Tl(:), basis_tmp(:, i))
    enddo
    call MPI_Allreduce(Tlocal, Tglobal, size(Tlocal), MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
    Tglobal(:) = Tglobal(:) - 273d0
    if(ismaster) then
      write(636,*) time, Tglobal
      write(*,*) "Tl=", time, Tglobal
    endif
  end do

  !--------------------------------------------
  ! Deallocate Matrices and Vectors
  !--------------------------------------------
  if (ismaster) write (*, *) "Deallocating matrices and vectors"
  call deallocMatVec()
  ! Finalize MPI
  call MPI_FINALIZE(mpi_err)

end program NURBScode

subroutine voltet_abs(x1, x2, x3, x4, vol)
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

  vol = abs(-vol/6.0d0)

end subroutine voltet_abs
