!======================================================================
!
!======================================================================
subroutine SparseGMRES_up(col, row, IBC, IPER, D_FLAG, P_FLAG, &
                          rhsGu, rhsGp, solu, solp, &
                          lhsK11, lhsG, lhsD1, lhsM, &
                          icntu, Utol, Kspaceu, Kspaceu_mn, &
                          NNODZu, NSHLu, NSD) 
  use mpi  
  implicit none

  integer, intent(in) :: icntu, NNODZu, NSHLu, NSD, &
                         col(NNODZu+1), row(icntu), &
                         IBC(NNODZu,2*NSD+1), IPER(NNODZu), &
                         D_FLAG(NNODZu), P_FLAG(NNODZu), &
                         Kspaceu, Kspaceu_mn

  real(8), intent(in) :: rhsGu(NNODZu,NSD), rhsGp(NNODZu), &
                         lhsK11(NSD*NSD,icntu), lhsG(NSD,icntu), &
                         lhsD1(NSD,icntu), lhsM(icntu)

  real(8), intent(inout) :: solu(NNODZu,NSD), solp(NNODZu)

  integer :: n, i, j, k, iK, iKs, jK, lK, itask, is, lenseg
  real(8) :: tmp2, abctmp(9)
    
  real(8) :: lhsKBdiagu(NNODZu,NSD*NSD), lhsKBdiagm(NNODZu,NSD*NSD), &
             lhsMdiag(NNODZu)
  
  real(8) :: rhstmpu(NNODZu,NSD), rhstmpm(NNODZu,NSD), &
             rhstmpp(NNODZu), temp1u(NNODZu,NSD), &
             temp1m(NNODZu,NSD), temp1p(NNODZu)
  
  real(8) :: uBrgu(NNODZu,NSD,Kspaceu+1), uBrgm(NNODZu,NSD,Kspaceu+1), &
             uBrgp(NNODZu,Kspaceu+1)
  
  real(8) :: HBrg(Kspaceu+1,Kspaceu), eBrg(Kspaceu+1), &
             yBrg(Kspaceu), &
             Rcos(Kspaceu), Rsin(Kspaceu), rr, unorm, &
             epsnrm, beta, ercheck, tmp, rr0, Utol, &
             Binv(NSD,NSD), flrate, unorm_ref, rrglob

  ! Algorithm from Shakib/Johan's theses
  rhstmpu(:,:) = RHSGu(:,:)
  rhstmpp(:)   = RHSGp(:)
  
  ! Precondition Residual
  lhsKBdiagu = 0.0d0
  lhsMdiag   = 0.0d0
  
  ! reading from CSR matrix, only when found in CSR
  ! rest is initialized 0 anyway
  do i = 1, NNODZu     
    do j = col(i), col(i+1)-1
      n = row(j)
      if (n == i) then   
        lhsKBdiagu(i,:) = LHSK11(:,j)
        lhsMdiag(i)     = LHSM(j)      
      end if
    end do     
  end do  

  ! communicate the block-diagonal ! LGGL    
  if (numnodes .gt. 1) then
    call commu(lhsKBdiagu, NSD*NSD, 'in ')
    call commu(lhsMdiag,   1,       'in ')
!!!    call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
    call commu(lhsKBdiagu, NSD*NSD, 'out')
    call commu(lhsMdiag,   1,       'out')
  end if
  
  ! invert block-diagonal: momentum  
  Binv = 0.0d0
  do i = 1, NNODZu
    Binv(1,1) = lhsKBdiagu(i,5)*lhsKBdiagu(i,9) &
              - lhsKBdiagu(i,8)*lhsKBdiagu(i,6)
    Binv(1,2) = lhsKBdiagu(i,8)*lhsKBdiagu(i,3) &
              - lhsKBdiagu(i,2)*lhsKBdiagu(i,9)
    Binv(1,3) = lhsKBdiagu(i,2)*lhsKBdiagu(i,6) &
              - lhsKBdiagu(i,3)*lhsKBdiagu(i,5)

    tmp2 = Binv(1,1)*lhsKBdiagu(i,1) &
         + Binv(1,2)*lhsKBdiagu(i,4) &
         + Binv(1,3)*lhsKBdiagu(i,7)

    tmp = 1.0d0/tmp2
     
    Binv(1,1) = Binv(1,1)*tmp
    Binv(1,2) = Binv(1,2)*tmp
    Binv(1,3) = Binv(1,3)*tmp

    Binv(2,1) = (lhsKBdiagu(i,6)*lhsKBdiagu(i,7) &
               - lhsKBdiagu(i,4)*lhsKBdiagu(i,9)) * tmp
    Binv(2,2) = (lhsKBdiagu(i,1)*lhsKBdiagu(i,9) &
               - lhsKBdiagu(i,7)*lhsKBdiagu(i,3)) * tmp
    Binv(2,3) = (lhsKBdiagu(i,4)*lhsKBdiagu(i,3) &
               - lhsKBdiagu(i,1)*lhsKBdiagu(i,6)) * tmp
    Binv(3,1) = (lhsKBdiagu(i,4)*lhsKBdiagu(i,8) &
               - lhsKBdiagu(i,5)*lhsKBdiagu(i,7)) * tmp
    Binv(3,2) = (lhsKBdiagu(i,7)*lhsKBdiagu(i,2) &
               - lhsKBdiagu(i,1)*lhsKBdiagu(i,8)) * tmp
    Binv(3,3) = (lhsKBdiagu(i,1)*lhsKBdiagu(i,5) &
               - lhsKBdiagu(i,2)*lhsKBdiagu(i,4)) * tmp
    lhsKBdiagu(i,1) = Binv(1,1)
    lhsKBdiagu(i,2) = Binv(1,2)
    lhsKBdiagu(i,3) = Binv(1,3)
    lhsKBdiagu(i,4) = Binv(2,1)
    lhsKBdiagu(i,5) = Binv(2,2)
    lhsKBdiagu(i,6) = Binv(2,3)
    lhsKBdiagu(i,7) = Binv(3,1)
    lhsKBdiagu(i,8) = Binv(3,2)
    lhsKBdiagu(i,9) = Binv(3,3)
  end do

  do i = 1, NNODZu
    if (lhsMdiag(i) == 0.0d0) then
      lhsMdiag(i) = 1.0d0
    else   
      lhsMdiag(i) = 1.0d0/lhsMdiag(i)
    end if
  end do

  ! Precondition residual
  do i = 1, NNODZu
    uBrgu(i,1,1) = lhsKBdiagu(i,1)*rhstmpu(i,1) + &
                   lhsKBdiagu(i,2)*rhstmpu(i,2) + &
                   lhsKBdiagu(i,3)*rhstmpu(i,3)
    uBrgu(i,2,1) = lhsKBdiagu(i,4)*rhstmpu(i,1) + &
                   lhsKBdiagu(i,5)*rhstmpu(i,2) + &
                   lhsKBdiagu(i,6)*rhstmpu(i,3)
    uBrgu(i,3,1) = lhsKBdiagu(i,7)*rhstmpu(i,1) + &
                   lhsKBdiagu(i,8)*rhstmpu(i,2) + &
                   lhsKBdiagu(i,9)*rhstmpu(i,3)

    uBrgp(i,1) = rhstmpp(i)*lhsMdiag(i)
  end do
    
  rhstmpu(:,:) = uBrgu(:,:,1)
  rhstmpp(:)   = uBrgp(:,1)
  
  !----------------------------------------------  
  ! calculate norm of rhs
  rr = sum(rhstmpp(:)*rhstmpp(:))
  do i = 1, NSD
    rr = rr + sum(rhstmpu(:,i)*rhstmpu(:,i))
  end do

  rrglob = rr

  if (numnodes > 1) then
    call MPI_ALLREDUCE(rrglob, rr, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
  end if

  unorm = sqrt(rr)

!!!  if (ismaster) 
!!!     &     write(*,*) "After-Prec. Residual L_2 Norm is...", unorm

  iKs = 0
  
  ! set up tolerance
  epsnrm = Utol*unorm
  unorm_ref = 1.0d2/unorm
  
  ! set up RHS of the Hessenberg's problem
  eBrg(:) = 0.0d0
  eBrg(1) = unorm
  
  ! normalize the first Krylov vector
  uBrgu(:,:,1) = uBrgu(:,:,1) / unorm
  uBrgp(:,1)   = uBrgp(:,1)   / unorm
  
  ! loop through GMRES iterations
  do iK = 1, Kspaceu
  
    iKs = iK
     
    ! matrix-vector product
    temp1u(:,:) = uBrgu(:,:,iKs)
    temp1p(:)   = uBrgp(:,iKs)
    !----------------------------------------------------------------
    ! Periodicity (Slave = Master) - GL
    if (numnodes > 1) then
      call commu(temp1u, NSD, 'out')
      call commu(temp1p,   1, 'out')
    end if
    ! sp uncommeneted this to account for IBC     
  !  do i = 1, NNODZu
  !    if ((IBC(i,1)==3).or.(IBC(i,2)==3).or.(IBC(i,3)==3)) then
  !      temp1u(i,:) = temp1u(IPER(i),:) ! Slave = Master
  !    end if 
 
  !    if (IBC(i,7)==3) then
  !      temp1p(i)   = temp1p(IPER(i)) ! Slave = Master     
  !    end if
  !  end do

    ! Product
    call SparseProdUP_3D(col, row, lhsK11, lhsG, lhsD1, lhsM, &
                         temp1u, temp1p, uBrgu(:,:,iKs+1), &
                         uBrgp(:,iKs+1), D_FLAG, P_FLAG, &
                         NNODZu, NSHLu, icntu, NSD)
    !sp uncommented this part as well to account for IBC
!!!    ! Communicate to Masters, Zero out Slaves - LG 
  !  do i = 1, NNODZu    
  !  if ((IBC(i,1)==3).or.(IBC(i,2)==3).or.(IBC(i,3)==3)) then
  !     uBrgu(IPER(i),:,iKs+1) = uBrgu(IPER(i),:,iKs+1) + &
  !                     uBrgu(i,:,iKs+1) ! Master = Master+Slave
  !     uBrgu(i,:,iKs+1) = 0.0d0 ! Slave = zero
  !    end if

  !   if (IBC(i,7)==3) then
  !     uBrgp(IPER(i),iKs+1) = uBrgp(IPER(i),iKs+1) +   & 
  !                   uBrgp(i,iKs+1) ! Master = Master+Slave
  !     uBrgp(i,iKs+1) = 0d+0 ! Slave = zero
  !    end if
  !  end do
    
    if (numnodes > 1) then     
      call commu(uBrgu(:,:,iKs+1), NSD, 'in ')
      call commu(uBrgp(:,iKs+1),     1, 'in ')
    end if
  
    !----------------------------------------------------------------------

    rhstmpu(:,:) = uBrgu(:,:,iKs+1)
    rhstmpp(:)   = uBrgp(:,iKs+1)
     
    ! Precondition product
    do i = 1, NNODZu    
      uBrgu(i,1,iKs+1) = lhsKBdiagu(i,1)*rhstmpu(i,1) + &
                         lhsKBdiagu(i,2)*rhstmpu(i,2) + &
                         lhsKBdiagu(i,3)*rhstmpu(i,3)
      uBrgu(i,2,iKs+1) = lhsKBdiagu(i,4)*rhstmpu(i,1) + &
                         lhsKBdiagu(i,5)*rhstmpu(i,2) + &
                         lhsKBdiagu(i,6)*rhstmpu(i,3)
      uBrgu(i,3,iKs+1) = lhsKBdiagu(i,7)*rhstmpu(i,1) + &
                         lhsKBdiagu(i,8)*rhstmpu(i,2) + &
                         lhsKBdiagu(i,9)*rhstmpu(i,3)
        
      uBrgp(i,iKs+1) = rhstmpp(i)*lhsMdiag(i)
    end do
       
    ! orthogonalize and get the norm
    do jK = 1, iKs+1  

      if (jK == 1) then
     
        rr = sum(uBrgp(:,iKs+1)*uBrgp(:,1)) !{u_{i+1}*u_1} vector
        do i = 1, NSD
          rr = rr + sum(uBrgu(:,i,iKs+1)*uBrgu(:,i,1))
        end do
        rrglob = rr 
        if (numnodes > 1) then
          call MPI_ALLREDUCE(rrglob, rr, 1, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, &
                             MPI_COMM_WORLD, mpi_err)  
        end if
       
        beta = rr
!!!        if (ismaster) write(*,*) iKs, jK, beta   
      else
  
        ! project off jK-1 vector
        uBrgu(:,:,iKs+1) = uBrgu(:,:,iKs+1) - beta*uBrgu(:,:,jK-1)
        uBrgp(:,iKs+1)   = uBrgp(:,iKs+1)   - beta*uBrgp(:,jK-1)
     
        rr = sum(uBrgp(:,iKs+1)*uBrgp(:,jK)) !{u_{i+1}*u_j} vector
       
        do i = 1, NSD
          rr = rr + sum(uBrgu(:,i,iKs+1)*uBrgu(:,i,jK))
        end do
        rrglob = rr
        if (numnodes > 1) then
          call MPI_ALLREDUCE(rrglob, rr, 1, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, &
                             MPI_COMM_WORLD, mpi_err)
        end if
       
        beta = rr
!!!        if (ismaster) write(*,*) iKs, jK, beta   
      end if
    
      HBrg(jK,iKs) = beta ! put this in the Hessenberg Matrix

    end do    
     
    unorm = sqrt(beta)
    HBrg(iKs+1,iKs) = unorm ! this fills the 1 sub diagonal band
    
    ! normalize the Krylov vector
    uBrgu(:,:,iKs+1) = uBrgu(:,:,iKs+1) / unorm
    uBrgp(:,iKs+1)   = uBrgp(:,iKs+1)   / unorm
    
    ! construct and reduce the Hessenberg Matrix
    ! since there is only one subdiagonal we can use a Givens rotation
    ! to rotate off each subdiagonal AS IT IS FORMED. We do this because it
    ! allows us to check progress of solution and quit when satisfied.  Note
    ! that all future K vects will put a subdiagonal in the next column so
    ! there is no penalty to work ahead as  the rotation for the next vector
    ! will be unaffected by this rotation.
     
    ! H Y = E ========>   R_i H Y = R_i E
    do jK = 1, iKs-1
      tmp = Rcos(jK)*HBrg(jK,iKs) + Rsin(jK)*HBrg(jK+1,iKs)

      HBrg(jK+1,iKs) = -Rsin(jK)*HBrg(jK,  iKs) + &
                        Rcos(jK)*HBrg(jK+1,iKs)
      HBrg(jK,  iKs) = tmp
    end do
     
    tmp     = sqrt(HBrg(iKs,iKs)**2 + HBrg(iKs+1,iKs)**2)
    Rcos(iKs)   = HBrg(iKs,  iKs) / tmp
    Rsin(iKs)   = HBrg(iKs+1,iKs) / tmp
    HBrg(iKs,  iKs) = tmp
    HBrg(iKs+1,iKs) = 0.0d0

    ! rotate eBrg    R_i E
    tmp = Rcos(iKs)*eBrg(iKs) + Rsin(iKs)*eBrg(iKs+1)
    eBrg(iKs+1) = -Rsin(iKs)*eBrg(iKs) + Rcos(iKs)*eBrg(iKs+1)
    eBrg(iKs)   =  tmp

    ! check for convergence
    ercheck = abs(eBrg(iKs+1))

    if ((ercheck <= epsnrm).and.(iKs >= Kspaceu_mn)) exit
     ! sp commented this for visual ease   
    if (mod(iKs,20) == 0) then
      if (ismaster) then
        write(*,8000) iKs, ercheck, epsnrm, ercheck*unorm_ref
      end if
    end if 

  end do    ! end of GMRES iteration loop


  ! -- Solution -------------------------------------------------
  ! if converged or end of Krylov space
  ! solve for yBrg
  do jK = iKs, 1, -1
    yBrg(jK) = eBrg(jK)/HBrg(jK,jK)
    do lK = 1, jK-1
      eBrg(lK) = eBrg(lK) - yBrg(jK)*HBrg(lK,jK)
    end do
  end do
  
  ! update du, dp
  do jK = 1, iKs
    solu(:,:) = solu(:,:) + yBrg(jK)*uBrgu(:,:,jK)
    solp(:)   = solp(:)   + yBrg(jK)*uBrgp(:,jK)
  end do

  ! communicate solution - GL
  if (numnodes > 1) then 
    call commu(solu, NSD, 'out')
    call commu(solp,   1, 'out')
  end if

  !---------------------------------------------------------------
 !sp commented this
 ! if (ismaster) write(*,8500) iKs, ercheck*unorm_ref
  
  8000 format(3x,i4,') Residual, Goal, Reduction= (%):', &
              1x,e10.4,1x,e10.4,1x,F12.6)
  8500 format(10x,' Iterations:', 2x,i4,'  Reduction:', &
              2x,f10.6,"(%)",/)
    
end subroutine SparseGMRES_up
