      subroutine SparseGMRES_m(col, row,
     &     IBC, IPER, D_FLAG, P_FLAG, rhsGm,
     &     solm, lhsK22, icntu, Utol, Kspaceu,Kspaceu_mn,
     &     NNODZu, NSHLu, NSD)
      
      use mpi
      implicit none
      
      integer :: n, i, j, k, iK, iKs, jK, lK, Kspaceu,Kspaceu_mn, 
     &           count,icntu, NNODZu, NSHLu, NSD, itask, is, lenseg
      
      integer :: col(NNODZu+1), row(icntu)      
      integer :: IBC(NNODZu,2*NSD+1), IPER(NNODZu)      
      integer :: D_FLAG(NNODZu), P_FLAG(NNODZu)
      
      real(8) :: rhsGu(NNODZu,NSD), rhsGm(NNODZu,NSD), rhsGp(NNODZu),
     &           solu(NNODZu,NSD), solm(NNODZu,NSD), solp(NNODZu)
      
      real(8) :: lhsK11(NSD*NSD,icntu), lhsK12(NSD*NSD,icntu),
     &           lhsK22(NSD*NSD,icntu), lhsG(NSD,icntu), 
     &           lhsD1(NSD,icntu), lhsD2(NSD,icntu), lhsM(icntu)
      
      real(8) :: lhsKBdiagu(NNODZu,NSD*NSD), lhsKBdiagm(NNODZu,NSD*NSD),
     &           lhsMdiag(NNODZu)
      
      real(8) :: rhstmpu(NNODZu,NSD), rhstmpm(NNODZu,NSD),
     &           rhstmpp(NNODZu), temp1u(NNODZu,NSD), 
     &           temp1m(NNODZu,NSD), temp1p(NNODZu)
      
      real(8) :: uBrgu(NNODZu,NSD,Kspaceu+1),
     &           uBrgm(NNODZu,NSD,Kspaceu+1),
     &           uBrgp(NNODZu,Kspaceu+1)
      
      real(8) :: HBrg(Kspaceu+1,Kspaceu), eBrg(Kspaceu+1),
     &           yBrg(Kspaceu), 
     &           Rcos(Kspaceu), Rsin(Kspaceu), rr, unorm,
     &           epsnrm, beta, ercheck, tmp, rr0, Utol,
     &           Binv(NSD,NSD), flrate, unorm_ref,rrglob
      
      !...  Algorithm from Shakib/Johan's theses

      rhstmpm(:,:) = RHSGm(:,:)
      
      ! Precondition Resudual
      
      lhsKBdiagm = 0d+0
      
      do i = 1, NNODZu
        do j = col(i), col(i+1)-1
          n = row(j)
          if (n == i) then
               
            lhsKBdiagm(i,:) = LHSK22(:,j)
                              
            !if (lhsKBdiagm(i,1) == 0d+0) lhsKBdiagm(i,1) = 1.0d0
            !if (lhsKBdiagm(i,5) == 0d+0) lhsKBdiagm(i,5) = 1.0d0
            !if (lhsKBdiagm(i,9) == 0d+0) lhsKBdiagm(i,9) = 1.0d0                              
          end if
        end do
      end do
             
      ! communicate the block-diagonal ! LGGL

!      do i = 1, NNODZu          ! If a node is a slave
!        if ((IBC(i,4).eq.3).or.(IBC(i,5).eq.3).or.(IBC(i,6).eq.3)) then 
!          lhsKBdiagm(IPER(i),:) = lhsKBdiagm(IPER(i),:)
!     &                          + lhsKBdiagm(i,:) ! Ma = Ma+Sl
!          lhsKBdiagm(i,:) = 0.0d0
!        end if
!      end do

      if (numnodes .gt. 1) then         
         call commu(lhsKBdiagm, NSD*NSD, 'in ')
         call commu(lhsKBdiagm, NSD*NSD, 'out')         
      endif      
      
!      do i = 1, NNODZu          ! If a node is a slave
!        if ((IBC(i,4).eq.3).or.(IBC(i,5).eq.3).or.(IBC(i,6).eq.3)) then 
!          lhsKBdiagm(i,:) = lhsKBdiagm(IPER(i),:)
!        endif     
!      enddo
      
      ! invert block-diagonal: mesh
      Binv = 0d+0
      do i = 1, NNODZu
         
        Binv(1,1) = lhsKBdiagm(i,5) * lhsKBdiagm(i,9) 
     &            - lhsKBdiagm(i,8) * lhsKBdiagm(i,6)
        Binv(1,2) = lhsKBdiagm(i,8) * lhsKBdiagm(i,3) 
     &            - lhsKBdiagm(i,2) * lhsKBdiagm(i,9)
        Binv(1,3) = lhsKBdiagm(i,2) * lhsKBdiagm(i,6) 
     &            - lhsKBdiagm(i,3) * lhsKBdiagm(i,5)

        tmp = 1d+0 / ( Binv(1,1) * lhsKBdiagm(i,1) 
     &               + Binv(1,2) * lhsKBdiagm(i,4)  
     &               + Binv(1,3) * lhsKBdiagm(i,7) )

        Binv(1,1) = Binv(1,1) * tmp
        Binv(1,2) = Binv(1,2) * tmp
        Binv(1,3) = Binv(1,3) * tmp

        Binv(2,1) = (lhsKBdiagm(i,6) * lhsKBdiagm(i,7) 
     &             - lhsKBdiagm(i,4) * lhsKBdiagm(i,9)) * tmp
        Binv(2,2) = (lhsKBdiagm(i,1) * lhsKBdiagm(i,9) 
     &             - lhsKBdiagm(i,7) * lhsKBdiagm(i,3)) * tmp
        Binv(2,3) = (lhsKBdiagm(i,4) * lhsKBdiagm(i,3) 
     &             - lhsKBdiagm(i,1) * lhsKBdiagm(i,6)) * tmp
        Binv(3,1) = (lhsKBdiagm(i,4) * lhsKBdiagm(i,8) 
     &             - lhsKBdiagm(i,5) * lhsKBdiagm(i,7)) * tmp
        Binv(3,2) = (lhsKBdiagm(i,7) * lhsKBdiagm(i,2) 
     &             - lhsKBdiagm(i,1) * lhsKBdiagm(i,8)) * tmp
        Binv(3,3) = (lhsKBdiagm(i,1) * lhsKBdiagm(i,5) 
     &             - lhsKBdiagm(i,2) * lhsKBdiagm(i,4)) * tmp
         
        lhsKBdiagm(i,1) = Binv(1,1)
        lhsKBdiagm(i,2) = Binv(1,2)
        lhsKBdiagm(i,3) = Binv(1,3)
        lhsKBdiagm(i,4) = Binv(2,1)
        lhsKBdiagm(i,5) = Binv(2,2)
        lhsKBdiagm(i,6) = Binv(2,3)
        lhsKBdiagm(i,7) = Binv(3,1)
        lhsKBdiagm(i,8) = Binv(3,2)
        lhsKBdiagm(i,9) = Binv(3,3)

      end do      
            
      ! Precondition residual
      do i = 1, NNODZu
        uBrgm(i,1,1) = lhsKBdiagm(i,1)*rhstmpm(i,1) +
     &                 lhsKBdiagm(i,2)*rhstmpm(i,2) +
     &                 lhsKBdiagm(i,3)*rhstmpm(i,3)
        uBrgm(i,2,1) = lhsKBdiagm(i,4)*rhstmpm(i,1) +
     &                 lhsKBdiagm(i,5)*rhstmpm(i,2) +
     &                 lhsKBdiagm(i,6)*rhstmpm(i,3)
        uBrgm(i,3,1) = lhsKBdiagm(i,7)*rhstmpm(i,1) +
     &                 lhsKBdiagm(i,8)*rhstmpm(i,2) +
     &                 lhsKBdiagm(i,9)*rhstmpm(i,3)
      end do
        
      rhstmpm(:,:) = uBrgm(:,:,1)
      
      ! calculate norm of rhs
      rr = 0d+0      
      do i = 1, NSD
        rr = rr + sum(rhstmpm(:,i)*rhstmpm(:,i))
      end do
      rrglob=rr
      if (numnodes .gt. 1) then         
         call MPI_ALLREDUCE (rrglob, rr, 1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, 
     &        MPI_COMM_WORLD,mpi_err)         
      endif

      unorm = sqrt(rr)
      
      iKs    = 0
      
      ! set up tolerance
      epsnrm = Utol * unorm
      unorm_ref = 1d2/unorm
            
      ! set up RHS of the Hessenberg's problem
      eBrg(:) = 0d+0
      eBrg(1) = unorm
      
      ! normalize the first Krylov vector
      uBrgm(:,:,1) = uBrgm(:,:,1) / unorm
      
      ! loop through GMRES iterations     
      do 1000 iK = 1, Kspaceu
        iKs = iK
         
        ! matrix-vector product
        temp1m(:,:) = uBrgm(:,:,iKs)

c--------------------------------------------------------------------c             
        ! Periodicity (Slave = Master) - GL
        if (numnodes .gt. 1) call commu(temp1m,NSD,'out')
                  
!        do i = 1, NNODZu
!          if ((IBC(i,4)==3).or.(IBC(i,5)==3).or.(IBC(i,6)==3)) then
!            temp1m(i,:) = temp1m(IPER(i),:) ! Slave = Master 
!          end if
!        end do
        
        ! Product        
        call SparseProdM_3D(col, row, lhsK22, temp1m, uBrgm(:,:,iKs+1),
     &                      D_FLAG, P_FLAG, NNODZu, NSHLu, icntu, NSD)
     
        !Communicate to Masters, Zero out Slaves - LG         
!        do i = 1, NNODZu
!          if ((IBC(i,4)==3).or.(IBC(i,5)==3).or.(IBC(i,6)==3)) then
!               
!            uBrgm(IPER(i),:,iKs+1) = uBrgm(IPER(i),:,iKs+1) + 
!     &                               uBrgm(i,:,iKs+1) ! Master = Master+Slave
!            uBrgm(i,:,iKs+1) = 0d+0 ! Slave = zero
!               
!          end if          
!        end do
                  
        if (numnodes .gt. 1) call commu(uBrgm(:,:,iKs+1),NSD,'in ')         
c-------------------------------------------------------------------------c
        rhstmpm(:,:) = uBrgm(:,:,iKs+1)
         
        ! Precondition product
        do i = 1, NNODZu
          uBrgm(i,1,iKs+1) = lhsKBdiagm(i,1)*rhstmpm(i,1) +
     &                       lhsKBdiagm(i,2)*rhstmpm(i,2) +
     &                       lhsKBdiagm(i,3)*rhstmpm(i,3)
          uBrgm(i,2,iKs+1) = lhsKBdiagm(i,4)*rhstmpm(i,1) +
     &                       lhsKBdiagm(i,5)*rhstmpm(i,2) +
     &                       lhsKBdiagm(i,6)*rhstmpm(i,3)
          uBrgm(i,3,iKs+1) = lhsKBdiagm(i,7)*rhstmpm(i,1) +
     &                       lhsKBdiagm(i,8)*rhstmpm(i,2) +
     &                       lhsKBdiagm(i,9)*rhstmpm(i,3)
        end do
               
        ! orthogonalize and get the norm
        do jK = 1, iKs+1
          if (jK == 1) then
	  
            rr = 0d+0               
            do i = 1, NSD
              rr = rr + sum(uBrgm(:,i,iKs+1)*uBrgm(:,i,1))
            end do
            rrglob=rr                 
            if (numnodes .gt. 1) then
            call MPI_ALLREDUCE(rrglob,rr,1,MPI_DOUBLE_PRECISION,
     &                         MPI_SUM,  MPI_COMM_WORLD,mpi_err)
            endif
               
            beta = rr
               
          else

            ! project off jK-1 vector
            uBrgm(:,:,iKs+1)=uBrgm(:,:,iKs+1)-beta * uBrgm(:,:,jK-1)
     
            rr = 0d+0               
            do i = 1, NSD
              rr = rr + sum(uBrgm(:,i,iKs+1)*uBrgm(:,i,jK))
            end do
            rrglob =rr
            if (numnodes .gt. 1) then
              call MPI_ALLREDUCE(rrglob, rr, 1,
     &	                         MPI_DOUBLE_PRECISION,MPI_SUM, 
     &                           MPI_COMM_WORLD, mpi_err)
            end if
               
            beta = rr
     
          end if
            
          HBrg(jK,iKs) = beta ! put this in the Hessenberg Matrix
     
        end do            
         
        unorm           = sqrt(beta)
        HBrg(iKs+1,iKs) = unorm ! this fills the 1 sub diagonal band
        
        ! normalize the Krylov vector
        uBrgm(:,:,iKs+1) = uBrgm(:,:,iKs+1) / unorm
        
c.... construct and reduce the Hessenberg Matrix
c     since there is only one subdiagonal we can use a Givens rotation
c     to rotate off each subdiagonal AS IT IS FORMED. We do this because it
c     allows us to check progress of solution and quit when satisfied.  Note
c     that all future K vects will put a subdiagonal in the next column so
c     there is no penalty to work ahead as  the rotation for the next vector
c     will be unaffected by this rotation.
         
c     
c     H Y = E ========>   R_i H Y = R_i E
c              
        do jK = 1, iKs-1
          tmp =  Rcos(jK) * HBrg(jK,  iKs) +
     &           Rsin(jK) * HBrg(jK+1,iKs)
          HBrg(jK+1,iKs) = -Rsin(jK) * HBrg(jK,  iKs) +
     &                      Rcos(jK) * HBrg(jK+1,iKs)
          HBrg(jK,  iKs) =  tmp
        end do

        tmp            = sqrt(HBrg(iKs,iKs)**2 + HBrg(iKs+1,iKs)**2)
        Rcos(iKs)      = HBrg(iKs,  iKs) / tmp
        Rsin(iKs)      = HBrg(iKs+1,iKs) / tmp
        HBrg(iKs,  iKs)= tmp
        HBrg(iKs+1,iKs)= 0d+0
     
        ! rotate eBrg    R_i E
        tmp        = Rcos(iKs) * eBrg(iKs) + Rsin(iKs) * eBrg(iKs+1)
        eBrg(iKs+1)=-Rsin(iKs) * eBrg(iKs) + Rcos(iKs) * eBrg(iKs+1)
        eBrg(iKs)  = tmp
    
        ! check for convergence
        ercheck=abs(eBrg(iKs+1))
         
       ! if (mod(iKs,100) == 0) then               
       !   if (ismaster) then
!	    write(*,8000) iKs,ercheck, epsnrm, ercheck*unorm_ref
!          endif
!	end if 

        if ((ercheck <= epsnrm).and.(iKs.ge.Kspaceu_mn)) exit          
      ! end of GMRES iteration loop

 1000 continue

c.... ------------------------->   Solution   <------------------------
      ! if converged or end of Krylov space
      ! solve for yBrg
      do jK = iKs, 1, -1
        yBrg(jK) = eBrg(jK) / HBrg(jK,jK)
        do lK = 1, jK-1
          eBrg(lK) = eBrg(lK) - yBrg(jK) * HBrg(lK,jK)
        end do
      end do
      
      ! update Dy
      do jK = 1, iKs
        solm(:,:) = solm(:,:) + yBrg(jK) * uBrgm(:,:,jK)
      end do
      
c---------------------------------------------------------------------c
      ! communicate solution - GL
      if (numnodes .gt. 1) then 
         call commu(solm, NSD, 'out')
      endif
      
      do i = 1, NNODZu         
        if ((IBC(i,4)==3).or.(IBC(i,5)==3).or.(IBC(i,6)==3)) then
          solm(i,:) = solm(IPER(i),:) ! Slave = Master
        end if
      end do

      if (ismaster) write(*,8500) iKs,ercheck*unorm_ref
      
      return
      
 8000 format( 3x,i4,') Residual, Goal, reduction (%):',
     &        1x,ES10.4,1x,ES10.4,1x,F12.6)
 8500 format(/,10x,' Iterations:',2x,i4,'  Reduction:',2x,f10.6,"(%)",/)
     
      end subroutine SparseGMRES_m
