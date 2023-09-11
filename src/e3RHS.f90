!======================================================================
!
!======================================================================
subroutine e3Rhs_3D_fluid(nshl, rhoi, mui, ui, aci, umi, acmi, uadvi, pri, &
                          rLi, fi, duidxi, ddidxi, tauM, tauP,  &
                          tauC, tauBar, kap_dc, gwt, shgu,      &
                          shgradgu, uprime, Rhsu, Rhsp)
  use aAdjKeep
  use commonvars
  implicit none
  
  integer, intent(in) :: nshl      
  real(8), intent(in) :: rhoi, mui
  real(8), intent(in) :: ui(NSD), aci(NSD), umi(NSD), acmi(NSD),  &
                         uadvi(NSD), pri, rLi(NSD), fi(NSD),      &
                         duidxi(NSD,NSD), ddidxi(NSD,NSD), gwt,   &
                         shgu(NSHL), shgradgu(NSHL,NSD), kap_dc,  &
                         tauM, tauP, tauC, tauBar, uprime(NSD)

  real(8), intent(inout) :: Rhsu(NSD,NSHL), Rhsp(NSHL)

  integer :: aa, bb, i, j, k
  real(8) :: fact1, fact2, mupkdc, divu, ptot, advu(NSD), &
             tmp1(NSD), tmp2(NSD,NSD), tmp4(NSD)

  tmp1 = 0.0d0; tmp2 = 0.0d0; tmp4 = 0.0d0; divu = 0.0d0

!!!  uprime(:) = -tauM*rLi(:)
      
  mupkdc = mui + kap_dc

  divu = duidxi(1,1) + duidxi(2,2) + duidxi(3,3)

  advu(:) = uadvi(:) + uprime(:)
          
  ptot = pri - tauC*divu

  tmp1(:) = rhoi*(aci(:) + advu(1)*duidxi(:,1) &
                        + advu(2)*duidxi(:,2) &
                        + advu(3)*duidxi(:,3)) - fi(:)
            
  tmp4(:) = uprime(1)*duidxi(:,1) + &
            uprime(2)*duidxi(:,2) + &
            uprime(3)*duidxi(:,3)
          
  tmp2(1,1) = -ptot + 2.0d0*mupkdc*duidxi(1,1)   &
            - rhoi*(uadvi(1)+uprime(1))*uprime(1) &
            + rhoi*uprime(1)*tauBar*tmp4(1)
      
  tmp2(1,2) = mupkdc*(duidxi(1,2) + duidxi(2,1)) &
            - rhoi*(uadvi(2)+uprime(2))*uprime(1) &
            + rhoi*uprime(2)*tauBar*tmp4(1)
      
  tmp2(1,3) = mupkdc*(duidxi(1,3) + duidxi(3,1)) &
            - rhoi*(uadvi(3)+uprime(3))*uprime(1) &
            + rhoi*uprime(3)*tauBar*tmp4(1)
      
  tmp2(2,1) = mupkdc*(duidxi(2,1) + duidxi(1,2)) &
            - rhoi*(uadvi(1)+uprime(1))*uprime(2) &
            + rhoi*uprime(1)*tauBar*tmp4(2) 
      
  tmp2(2,2) = -ptot + 2.0d0*mupkdc*duidxi(2,2)   &
            - rhoi*(uadvi(2)+uprime(2))*uprime(2) &
            + rhoi*uprime(2)*tauBar*tmp4(2)
      
  tmp2(2,3) = mupkdc*(duidxi(2,3) + duidxi(3,2)) &
            - rhoi*(uadvi(3)+uprime(3))*uprime(2) &
            + rhoi*uprime(3)*tauBar*tmp4(2)
      
  tmp2(3,1) = mupkdc*(duidxi(3,1) + duidxi(1,3)) &
            - rhoi*(uadvi(1)+uprime(1))*uprime(3) &
            + rhoi*uprime(1)*tauBar*tmp4(3)
     
  tmp2(3,2) = mupkdc*(duidxi(3,2) + duidxi(2,3)) &
            - rhoi*(uadvi(2)+uprime(2))*uprime(3) &
            + rhoi*uprime(2)*tauBar*tmp4(3)
      
  tmp2(3,3) = -ptot + 2.0d0*mupkdc*duidxi(3,3)   &
            - rhoi*(uadvi(3)+uprime(3))*uprime(3) &
            + rhoi*uprime(3)*tauBar*tmp4(3)
      
  ! Physics Residual   
  do aa = 1, NSHL     
    do i = 1, NSD
      Rhsu(i,aa) = Rhsu(i,aa) - &
                   ( shgu(aa)*tmp1(i) + &
                     sum(shgradgu(aa,:)*tmp2(i,:)) )*DetJ*gwt
    end do
  end do
      
  do aa = 1, NSHL
    Rhsp(aa) = Rhsp(aa) - &
               ( shgu(aa)*divu - &
                 sum(shgradgu(aa,:)*uprime(:)) )*DetJ*gwt
  end do

end subroutine e3Rhs_3D_fluid




!======================================================================
! RHS for weak BC
!======================================================================
subroutine e3bRHS_weak(nshl, rhoi, mui, nor, tauB, tauNor, gwt, &
                       shlu, shgradgu, ui, umi, pri, duidxi, gi, &
                       Rhsu, Rhsp, ti, tmp1, tmp2)  
  use aAdjKeep  
  use commonvars
  implicit none
  
  integer, intent(in) :: nshl

  real(8), intent(in)    :: rhoi, mui
  real(8), intent(inout) :: Rhsu(NSD,NSHL), Rhsp(NSHL)
  real(8), intent(inout)   :: ti(NSD), tmp1(NSD), tmp2(NSD,NSD)
  real(8), intent(in)    :: shlu(NSHL), shgradgu(NSHL,NSD), &
                            nor(NSD), tauB, tauNor, gwt, &
                            ui(NSD), umi(NSD), pri, duidxi(NSD,NSD), &
                            gi(NSD)

  integer :: aa, bb, i  
  real(8) :: fact1, fact2, upos, uneg, unor, tr, pt33, gmul, gnor
  
  tmp1 = 0.0d0
  tmp2 = 0.0d0
  
  ti(:) = -pri*nor(:) + mui*(duidxi(:,1)*nor(1) + &
                             duidxi(:,2)*nor(2) + &
                             duidxi(:,3)*nor(3) + &
                             duidxi(1,:)*nor(1) + &
                             duidxi(2,:)*nor(2) + &
                             duidxi(3,:)*nor(3))
  
  ! Relative normal velocity for convective term
  unor = sum((ui-umi)*nor)  ! u \cdot n
  upos = 0.5d0*(unor+abs(unor))
  uneg = 0.5d0*(unor-abs(unor))  

  ! Absolute normal for the rest
  gnor = sum(gi*nor)
  unor = sum(ui*nor)  ! u \cdot n
 
  tmp1(:) = - ti(:) &
            + tauB*(ui(:)-gi(:)) &           
            - uneg*rhoi*(ui(:)-gi(:))      &
            + (tauNor-tauB)*unor*nor(:)

  ! gmul =  1.0d0 => sym
  ! gmul = -1.0d0 => skew       
  gmul = 1.0d0 
  do aa = 1, NSD
    do bb = 1, NSD
      tmp2(aa,bb) = -gmul*mui*((ui(aa)-gi(aa))*nor(bb) &
                             +(ui(bb)-gi(bb))*nor(aa))
    end do
  end do

  ! gnor = 0.0d0
  do aa = 1, NSHL
    do i = 1, NSD
      Rhsu(i,aa) = Rhsu(i,aa) - &
                 ( shlu(aa)*tmp1(i) + &
                   sum(shgradgu(aa,:)*tmp2(i,:)) )*DetJb*gwt
    end do

    Rhsp(aa) = Rhsp(aa) + shlu(aa)*(unor-gnor)*DetJb*gwt     
  end do

end subroutine e3bRHS_weak


