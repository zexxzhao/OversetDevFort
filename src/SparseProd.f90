!======================================================================
!
!======================================================================
subroutine SparseProdUP_3D(col, row, lhsK11, lhsG, lhsD1, lhsM,  &
                           rhstmpu, rhstmpp, prodtmpu, prodtmpp, &
                           D_FLAG, P_FLAG, &
                           NNODZu, NSHLu, icntu, NSD)
  implicit none
 
  integer, intent(in) :: col(NNODZu+1), row(NNODZu*27*NSHLu), &
                         NNODZu, NSHLu, icntu, NSD, &
                         D_FLAG(NNODZu), P_FLAG(NNODZu)

  real(8), intent(in) :: lhsK11(NSD*NSD,icntu), lhsG(NSD,icntu),  &
                         lhsD1(NSD,icntu), lhsM(icntu), &
                         rhstmpu(NNODZu,NSD), rhstmpp(NNODZu)

  real(8), intent(out) :: prodtmpu(NNODZu,NSD), prodtmpp(NNODZu)

  integer :: aa, bb, cc
  real(8) :: tmpvect(NNODZu,NSD)  
  real(8) :: tmp(4), pisave
  
  ! clear the vector
  prodtmpu = 0.0d0
  prodtmpp = 0.0d0
  
  do aa = 1, NNODZu     ! K*u
    
    tmp = 0.0d0
    do bb = col(aa), col(aa+1)-1
      cc = row(bb)
      
      tmp(1) = tmp(1) + LHSK11(1,bb)*rhstmpu(cc,1) + &
                        LHSK11(2,bb)*rhstmpu(cc,2) + &
                        LHSK11(3,bb)*rhstmpu(cc,3) + &
                        LHSG  (1,bb)*rhstmpp(cc)
      
      tmp(2) = tmp(2) + LHSK11(4,bb)*rhstmpu(cc,1) + &
                        LHSK11(5,bb)*rhstmpu(cc,2) + &
                        LHSK11(6,bb)*rhstmpu(cc,3) + &
                        LHSG  (2,bb)*rhstmpp(cc)
        
      tmp(3) = tmp(3) + LHSK11(7,bb)*rhstmpu(cc,1) + &
                        LHSK11(8,bb)*rhstmpu(cc,2) + &
                        LHSK11(9,bb)*rhstmpu(cc,3) + &
                        LHSG  (3,bb)*rhstmpp(cc)

      tmp(4) = tmp(4) + LHSD1(1,bb)*rhstmpu(cc,1) + &
                        LHSD1(2,bb)*rhstmpu(cc,2) + &
                        LHSD1(3,bb)*rhstmpu(cc,3) + &
                        LHSM   (bb)*rhstmpp(cc)      
    end do
    
    prodtmpu(aa,:) = prodtmpu(aa,:) + tmp(1:3)
    prodtmpp(aa)   = prodtmpp(aa)   + tmp(4)
    
  end do

end subroutine SparseProdUP_3D

  
!======================================================================
!
!======================================================================
subroutine SparseProdM_3D(col, row, &
       lhsK22, rhstmpm, prodtmpm, &
       D_FLAG, P_FLAG, &
       NNODZu, NSHLu, icntu, NSD)

  implicit none
  
  integer :: aa, bb, cc, NNODZu, NSHLu, icntu, NSD  
  integer :: col(NNODZu+1), row(NNODZu*27*NSHLu)
  integer :: D_FLAG(NNODZu), P_FLAG(NNODZu)
  
  real(8) :: lhsK11(NSD*NSD,icntu), lhsK12(NSD*NSD,icntu), &
       lhsK22(NSD*NSD,icntu), lhsG(NSD,icntu), lhsD1(NSD,icntu), &
       lhsD2(NSD,icntu), lhsM(icntu)

  real(8) :: rhstmpu(NNODZu,NSD), rhstmpm(NNODZu,NSD), &
       rhstmpp(NNODZu)

  real(8) :: prodtmpu(NNODZu,NSD), prodtmpm(NNODZu,NSD), &
             prodtmpp(NNODZu)

  real(8) :: tmpvect(NNODZu,NSD)  
  real(8) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, pisave
  
  ! clear the vector
     
  prodtmpm(1:NNODZu,:) = 0.0d0
  do aa = 1, NNODZu     ! K*u
    
    tmp1 = 0d+0
    tmp2 = 0d+0
    tmp3 = 0d+0
    do bb = col(aa), col(aa+1)-1
      cc = row(bb)
      
      tmp1 = tmp1 + LHSK22(1,bb)*rhstmpm(cc,1) + &
                    LHSK22(2,bb)*rhstmpm(cc,2) +&
                    LHSK22(3,bb)*rhstmpm(cc,3)

      tmp2 = tmp2 + LHSK22(4,bb)*rhstmpm(cc,1) + &
                    LHSK22(5,bb)*rhstmpm(cc,2) +&
                    LHSK22(6,bb)*rhstmpm(cc,3)

      tmp3 = tmp3 + LHSK22(7,bb)*rhstmpm(cc,1) + &
                    LHSK22(8,bb)*rhstmpm(cc,2) +&
                    LHSK22(9,bb)*rhstmpm(cc,3)
      
    end do
    
    prodtmpm(aa,1) = prodtmpm(aa,1) + tmp1
    prodtmpm(aa,2) = prodtmpm(aa,2) + tmp2
    prodtmpm(aa,3) = prodtmpm(aa,3) + tmp3
    
  end do

end subroutine SparseProdM_3D



!======================================================================
!
!======================================================================
subroutine SparseProd_LS( &
          lhsM, &
          col, row, &
          rhstmp, prodtmp, &
          NNODZu, NSHLu, icntu, NSD)
  
  implicit none
  
  integer NNODZu, NSHLu, icntu, NSD  
  real(8) lhsM(icntu)
  integer col(NNODZu+1), row(NNODZu*27*NSHLu)
  real(8) rhstmp(NNODZu),prodtmp(NNODZu)
  integer aa, bb, cc
  real(8) tmp
  
  ! clear the vector  
  prodtmp = 0d+0
  
  do aa = 1, NNODZu     ! K*u
     
     tmp = 0d+0

     do bb = col(aa), col(aa+1)-1

    cc = row(bb)
    tmp= tmp + lhsM(bb)*rhstmp(cc)
    
     enddo
     
     prodtmp(aa) = prodtmp(aa) + tmp
     
  enddo

end subroutine SparseProd_LS


!======================================================================
!
!======================================================================
subroutine SparseProd_NS_conv(col, row, &
       lhsK11, lhsG, lhsD1, lhsM, lhsLS,lhsLSu,LHSPls,LHSUls, &
       rhstmpu, rhstmpp, rhstmpls, &
       prodtmpu, prodtmpp, prodtmpls, &
       D_FLAG, P_FLAG, &
       NNODZu, NSHLu, icntu, NSD)

  implicit none
  
  integer :: aa, bb, cc, NNODZu, NSHLu, icntu, NSD
  
  integer :: col(NNODZu+1), row(NNODZu*27*NSHLu)

  integer :: D_FLAG(NNODZu), P_FLAG(NNODZu)
  
  real(8) :: lhsK11(NSD*NSD,icntu), lhsG(NSD,icntu),  &
             lhsD1(NSD,icntu), lhsls(icntu), lhsM(icntu), &
             lhsLSu(NSD,icntu),LHSPls(icntu),LHSUls(NSD,icntu)

  real(8) :: rhstmpu(NNODZu,NSD), rhstmpls(NNODZu), &
             rhstmpp(NNODZu)

  real(8) :: prodtmpu(NNODZu,NSD), prodtmpls(NNODZu), &
             prodtmpp(NNODZu)
 
  real(8) :: tmpvect(NNODZu,NSD)
  
  real(8) :: tmp(5), pisave
  
  ! clear the vector
  prodtmpu(1:NNODZu,:) = 0.0d0
  prodtmpp(1:NNODZu)   = 0.0d0
  prodtmpls(1:NNODZu)  = 0.0d0
  
  do aa = 1, NNODZu     ! K*u
    
    tmp = 0.0d0
    do bb = col(aa), col(aa+1)-1
      cc = row(bb)
      
      tmp(1) = tmp(1) + LHSK11(1,bb)*rhstmpu (cc,1) + &
                        LHSK11(2,bb)*rhstmpu (cc,2) + &
                        LHSK11(3,bb)*rhstmpu (cc,3) + &
                        LHSUls(1,bb)*rhstmpls(cc)   + &
                        LHSG  (1,bb)*rhstmpp (cc)
      
      tmp(2) = tmp(2) + LHSK11(4,bb)*rhstmpu (cc,1) + &
                        LHSK11(5,bb)*rhstmpu (cc,2) + &
                        LHSK11(6,bb)*rhstmpu (cc,3) + &
                        LHSUls(2,bb)*rhstmpls(cc)   + &
                        LHSG  (2,bb)*rhstmpp (cc)
        
      tmp(3) = tmp(3) + LHSK11(7,bb)*rhstmpu (cc,1) + &
                        LHSK11(8,bb)*rhstmpu (cc,2) + &
                        LHSK11(9,bb)*rhstmpu (cc,3) + &
                        LHSUls(3,bb)*rhstmpls(cc)   + &
                        LHSG  (3,bb)*rhstmpp (cc)

      tmp(4) = tmp(4) + LHSD1(1,bb)*rhstmpu (cc,1) + &
                        LHSD1(2,bb)*rhstmpu (cc,2) + &
                        LHSD1(3,bb)*rhstmpu (cc,3) + & 
                        LHSPls (bb)*rhstmpls(cc)   + &
                        LHSM   (bb)*rhstmpp (cc)     
      
      tmp(5) = tmp(5) + lhsLSu(1,bb)*rhstmpu (cc,1) + &
                        lhsLSu(2,bb)*rhstmpu (cc,2) + &
                        lhsLSu(3,bb)*rhstmpu (cc,3) + &
                        lhsLs   (bb)*rhstmpls(cc)
    end do
    
    prodtmpu (aa,:) = prodtmpu (aa,:) + tmp(1:3)
    prodtmpp (aa)   = prodtmpp (aa)   + tmp(4)
    prodtmpls(aa)   = prodtmpls(aa)   + tmp(5)
  end do

end subroutine SparseProd_NS_conv
