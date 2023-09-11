!======================================================================
!
!======================================================================
subroutine e3LHS_3D_fluid_Old(nshl, rhoi, mui, ui, umi, aci, pri, duidxi,       &
                              dpridxi, rLi, tauM, tauP, tauC,        &
                              tauBar, kap_dc, gwt, shgu, shgradgu,   &
                              shhessgu, xKebe11,                     &
                              xGebe, xDebe1, xMebe)    
  use aAdjKeep
  use commonvars
  implicit none  
  
  integer, intent(in) :: nshl
  
  real(8), intent(in) :: rhoi, mui
  real(8), intent(in) :: ui(NSD), umi(NSD), aci(NSD), pri, kap_dc, &
                         duidxi(NSD,NSD), dpridxi(NSD), rLi(NSD), &
                         shgu(NSHL), shgradgu(NSHL,NSD), &
                         shhessgu(NSHL,NSD,NSD), &
                         gwt, tauM, tauP, tauC, tauBar

  real(8), intent(inout) :: xKebe11(NSD*NSD,NSHL,NSHL), &
                            xGebe(NSD,NSHL,NSHL),&
                            xDebe1(NSD,NSHL,NSHL),& 
                            xMebe(NSHL,NSHL)

  integer :: aa, bb, i, j
 
  real(8) :: fact1, fact2, fact3,&
             tmp1(NSHL), tmp2(NSHL), tmp4(NSHL,NSHL),&
             advu1(NSD), advu2(NSD), mupkdc
  
  ! loop over local shape functions in each direction
  fact1 = almi
  fact2 = alfi*gami*Delt
  fact3 = alfi*beti*Delt*Delt
    
  mupkdc = mui + kap_dc

  tmp1 = 0.0d0
  tmp2 = 0.0d0
  tmp4 = 0.0d0

  advu1(:) = ui(:)-umi(:)
  advu2(:) = -tauM*rLi(:)
    
  tmp1(:) = rhoi*(advu1(1)*shgradgu(:,1) + &! Na,_j (u_j-um_j)
                 advu1(2)*shgradgu(:,2) + &
                 advu1(3)*shgradgu(:,3))
  
  tmp2(:) = advu2(1)*shgradgu(:,1) + &! Na,_i (-tauM*Li)
            advu2(2)*shgradgu(:,2) + &
            advu2(3)*shgradgu(:,3)
  
  
  do bb = 1, NSHL      ! Diagonal blocks of K11
    do aa = 1, NSHL
      
      tmp4(aa,bb) = fact1*(shgu(aa)*rhoi*shgu(bb)+&
                           tmp1(aa)*tauM*rhoi*shgu(bb)) +&
                           fact2*(shgu(aa)*tmp1(bb) +&
                           tmp1(aa)*tauM*tmp1(bb) +&
                           mupkdc*(shgradgu(aa,1)*shgradgu(bb,1) +&
                                   shgradgu(aa,2)*shgradgu(bb,2) +&
                                   shgradgu(aa,3)*shgradgu(bb,3)) +&
                           tmp2(aa)*rhoi*tauBar*tmp2(bb))
    end do
  end do
 
  ! Physics-Physics Interaction
  do bb = 1, NSHL
    do aa = 1, NSHL      
      xKebe11(1,aa,bb) = xKebe11(1,aa,bb) +&
        (tmp4(aa,bb) +&
         fact2*(shgradgu(aa,1)*mupkdc*shgradgu(bb,1) +&
                shgradgu(aa,1)*tauC*shgradgu(bb,1)))*DetJ*gwt
      
      xKebe11(5,aa,bb) = xKebe11(5,aa,bb) +&
        (tmp4(aa,bb) +&
         fact2*(shgradgu(aa,2)*mupkdc*shgradgu(bb,2) +&
                shgradgu(aa,2)*tauC*shgradgu(bb,2)))*DetJ*gwt
      
      xKebe11(9,aa,bb) = xKebe11(9,aa,bb) +&
        (tmp4(aa,bb) +&
         fact2*(shgradgu(aa,3)*mupkdc*shgradgu(bb,3) +&
                shgradgu(aa,3)*tauC*shgradgu(bb,3)))*DetJ*gwt
      
      xKebe11(2,aa,bb) = xKebe11(2,aa,bb) + &
        fact2*(shgradgu(aa,2)*mupkdc*shgradgu(bb,1) +&
               shgradgu(aa,1)*tauC*shgradgu(bb,2))*DetJ*gwt
      
      xKebe11(4,aa,bb) = xKebe11(4,aa,bb) + &
        fact2*(shgradgu(aa,1)*mupkdc*shgradgu(bb,2) +&
               shgradgu(aa,2)*tauC*shgradgu(bb,1))*DetJ*gwt
      
      xKebe11(3,aa,bb) = xKebe11(3,aa,bb) + &
        fact2*(shgradgu(aa,3)*mupkdc*shgradgu(bb,1) +&
               shgradgu(aa,1)*tauC*shgradgu(bb,3))*DetJ*gwt
      
      xKebe11(7,aa,bb) = xKebe11(7,aa,bb) + &
        fact2*(shgradgu(aa,1)*mupkdc*shgradgu(bb,3) +&
               shgradgu(aa,3)*tauC*shgradgu(bb,1))*DetJ*gwt
      
      xKebe11(6,aa,bb) = xKebe11(6,aa,bb) + &
        fact2*(shgradgu(aa,3)*mupkdc*shgradgu(bb,2) +&
               shgradgu(aa,2)*tauC*shgradgu(bb,3))*DetJ*gwt
      
      xKebe11(8,aa,bb) = xKebe11(8,aa,bb) + &
        fact2*(shgradgu(aa,2)*mupkdc*shgradgu(bb,3) +&
               shgradgu(aa,3)*tauC*shgradgu(bb,2))*DetJ*gwt
    end do    
  end do

  ! Physics-Mesh
  ! Divergence Matrix  
  do bb = 1, NSHL   
    do aa = 1, NSHL      
      xGebe(1,aa,bb) = xGebe(1,aa,bb) + &
        fact2*(-shgradgu(aa,1)*shgu(bb) +&
               tmp1(aa)*tauM*shgradgu(bb,1))*DetJ*gwt
      xGebe(2,aa,bb) = xGebe(2,aa,bb) + &
        fact2*(-shgradgu(aa,2)*shgu(bb) +&
               tmp1(aa)*tauM*shgradgu(bb,2))*DetJ*gwt
      xGebe(3,aa,bb) = xGebe(3,aa,bb) + &
        fact2*(-shgradgu(aa,3)*shgu(bb) +&
               tmp1(aa)*tauM*shgradgu(bb,3))*DetJ*gwt
    end do    
  end do

  ! Physics-Physics
  ! Divergence Matrix  
  do bb = 1, NSHL
    do aa = 1, NSHL      
      xDebe1(1,aa,bb) = xDebe1(1,aa,bb) +&
        (fact2*(shgu(aa)*shgradgu(bb,1)+&
                shgradgu(aa,1)*tauP*tmp1(bb)) +&
         fact1*(shgradgu(aa,1)*tauP*rhoi*shgu(bb)))*DetJ*gwt

      xDebe1(2,aa,bb) = xDebe1(2,aa,bb) +&
        (fact2*(shgu(aa)*shgradgu(bb,2)+&
                shgradgu(aa,2)*tauP*tmp1(bb)) +&
         fact1*(shgradgu(aa,2)*tauP*rhoi*shgu(bb)))*DetJ*gwt
      
      xDebe1(3,aa,bb) = xDebe1(3,aa,bb) +&
        (fact2*(shgu(aa)*shgradgu(bb,3)+&
                shgradgu(aa,3)*tauP*tmp1(bb)) +&
         fact1*(shgradgu(aa,3)*tauP*rhoi*shgu(bb)))*DetJ*gwt   
    end do    
  end do   

  ! Mass Matrix
  do bb = 1, NSHL
    do aa = 1, NSHL
      xMebe(aa,bb) = xMebe(aa,bb) +&
        fact2*tauP*(shgradgu(aa,1)*shgradgu(bb,1) +&
                    shgradgu(aa,2)*shgradgu(bb,2) +&
                    shgradgu(aa,3)*shgradgu(bb,3))*DetJ*gwt
    end do
  end do
  
end subroutine e3LHS_3D_fluid_Old


!======================================================================
! LHS for weak BC
!======================================================================
subroutine e3bLHS_weak(nshl, rhoi, mui, ui, umi, duidxi, tauB, tauNor, gwt, &
                       shlu, shgradgu, xKebe, xGebe, xDebe, nor)
  
  use aAdjKeep  
  use commonvars
  implicit none
  
  integer, intent(inout) :: nshl
  
  real(8), intent(in) :: rhoi, mui
  real(8), intent(in) :: ui(NSD), umi(NSD), duidxi(NSD,NSD), &
                         tauB, tauNor, gwt, nor(NSD), &
                         shlu(NSHL), shgradgu(NSHL,NSD)
  real(8), intent(inout) :: xKebe(NSD*NSD,NSHL,NSHL), &
                            xGebe(NSD,NSHL,NSHL), &
                            xDebe(NSD,NSHL,NSHL)

  integer :: aa, bb
  real(8) :: fact1, fact2, tmp1(NSHL), tmp2(NSHL,NSHL), &
             unor, umul, munor, gmul, uneg 
  
  ! loop over local shape functions in each direction

  fact1 = almi
  fact2 = alfi*gami*Delt

  tmp1 = 0.0d0
  tmp2 = 0.0d0

  tmp1(:) = shgradgu(:,1)*nor(1) + shgradgu(:,2)*nor(2) &
          + shgradgu(:,3)*nor(3) 

  unor = sum((ui-umi)*nor(:))  ! u \cdot n
  uneg = 0.5d0*(unor-abs(unor))
  munor = tauNor-tauB
  
  ! gmul =  1d0 => sym
  ! gmul = -1d0 => skew  
  gmul = 1.0d0
  do bb = 1, NSHL      ! Diagonal blocks of K
    do aa = 1, NSHL
    
      tmp2(aa,bb) = -shlu(aa)*mui*tmp1(bb) &
                  - gmul*tmp1(aa)*mui*shlu(bb) &
                  + shlu(aa)*tauB*shlu(bb) &
                  - shlu(aa)*uneg*rhoi*shlu(bb)
    end do
  end do


  do bb = 1, NSHL    
    do aa = 1, NSHL
      
      xKebe(1,aa,bb) = xKebe(1,aa,bb) + &
        fact2*( tmp2(aa,bb)                            &
              - shlu(aa)*mui*shgradgu(bb,1)*nor(1)      &
              - gmul*shgradgu(aa,1)*nor(1)*mui*shlu(bb) &
              + shlu(aa)*nor(1)*munor*nor(1)*shlu(bb) )*DetJb*gwt
       
      xKebe(5,aa,bb) = xKebe(5,aa,bb) + &
        fact2*(  tmp2(aa,bb)                           &
              - shlu(aa)*mui*shgradgu(bb,2)*nor(2)      &
              - gmul*shgradgu(aa,2)*nor(2)*mui*shlu(bb) &
              + shlu(aa)*nor(2)*munor*nor(2)*shlu(bb) )*DetJb*gwt
       
      xKebe(9,aa,bb) = xKebe(9,aa,bb) + &
        fact2*( tmp2(aa,bb)                            &
              - shlu(aa)*mui*shgradgu(bb,3)*nor(3)      &
              - gmul*shgradgu(aa,3)*nor(3)*mui*shlu(bb) &
              + shlu(aa)*nor(3)*munor*nor(3)*shlu(bb) )*DetJb*gwt
       
      xKebe(2,aa,bb) = xKebe(2,aa,bb) + &
        fact2*(-shlu(aa)*mui*shgradgu(bb,1)*nor(2)      &
              - gmul*shgradgu(aa,2)*nor(1)*mui*shlu(bb) &
              + shlu(aa)*nor(1)*munor*nor(2)*shlu(bb))*DetJb*gwt

      xKebe(4,aa,bb) = xKebe(4,aa,bb) +  &
        fact2*(-shlu(aa)*mui*shgradgu(bb,2)*nor(1)      &
              - gmul*shgradgu(aa,1)*nor(2)*mui*shlu(bb) &
              + shlu(aa)*nor(2)*munor*nor(1)*shlu(bb))*DetJb*gwt
       
      xKebe(3,aa,bb) = xKebe(3,aa,bb) +  &
        fact2*(-shlu(aa)*mui*shgradgu(bb,1)*nor(3)      &
              - gmul*shgradgu(aa,3)*nor(1)*mui*shlu(bb) &
              + shlu(aa)*nor(1)*munor*nor(3)*shlu(bb))*DetJb*gwt

      xKebe(7,aa,bb) = xKebe(7,aa,bb) +  &
        fact2*(-shlu(aa)*mui*shgradgu(bb,3)*nor(1)      &
              - gmul*shgradgu(aa,1)*nor(3)*mui*shlu(bb) &
              + shlu(aa)*nor(3)*munor*nor(1)*shlu(bb))*DetJb*gwt

      xKebe(6,aa,bb) = xKebe(6,aa,bb) +  &
        fact2*(-shlu(aa)*mui*shgradgu(bb,2)*nor(3)      &
              - gmul*shgradgu(aa,3)*nor(2)*mui*shlu(bb) &
              + shlu(aa)*nor(2)*munor*nor(3)*shlu(bb))*DetJb*gwt

      xKebe(8,aa,bb) = xKebe(8,aa,bb) + &
        fact2*(-shlu(aa)*mui*shgradgu(bb,3)*nor(2)      &
              - gmul*shgradgu(aa,2)*nor(3)*mui*shlu(bb) &
              + shlu(aa)*nor(3)*munor*nor(2)*shlu(bb))*DetJb*gwt

    end do    
  end do
  
  do bb = 1, NSHL   
    do aa = 1, NSHL
      xDebe(1,aa,bb) = xDebe(1,aa,bb) - &
                       fact2*shlu(aa)*shlu(bb)*nor(1)*DetJb*gwt
      xDebe(2,aa,bb) = xDebe(2,aa,bb) - &
                       fact2*shlu(aa)*shlu(bb)*nor(2)*DetJb*gwt
      xDebe(3,aa,bb) = xDebe(3,aa,bb) - &
                       fact2*shlu(aa)*shlu(bb)*nor(3)*DetJb*gwt
    end do
  end do 


  do bb = 1, NSHL
    do aa = 1, NSHL    
      xGebe(1,aa,bb) = xGebe(1,aa,bb) +  & 
                       fact2*shlu(aa)*shlu(bb)*nor(1)*DetJb*gwt
      xGebe(2,aa,bb) = xGebe(2,aa,bb) +  &
                       fact2*shlu(aa)*shlu(bb)*nor(2)*DetJb*gwt
      xGebe(3,aa,bb) = xGebe(3,aa,bb) +  &
                       fact2*shlu(aa)*shlu(bb)*nor(3)*DetJb*gwt    
    end do
  end do

end subroutine e3bLHS_weak


