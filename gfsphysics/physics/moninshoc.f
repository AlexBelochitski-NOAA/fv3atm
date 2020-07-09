!!!!!  ==========================================================  !!!!!
! subroutine 'moninshoc' computes pbl height and applies vertical diffusion
! using the coefficient provided by the SHOC scheme (from previous step)
! 2015-05-04 - Shrinivas Moorthi - original version based on monin
! 2018-03-21 - Shrinivas Moorthi - fixed a bug related to tke vertical diffusion
!                                  and gneralized the tke location in tracer array
! 2018-03-23 - Shrinivas Moorthi - used twice the momentum diffusion coefficient
!                                  for tke as in Deardorff (1980) - added tridi1
! 2019       - Alex Belochitski  - updated for SHOCv2. Addded options for PBL height
!                                  and Prandtl number calculations. Background diffusion
!                                  is now assumed to be minimum value of diffusion
!                                  coefficient.    

      subroutine moninshoc(ix,im,km,ntrac,ntcw,ncnd,dv,du,tau,rtg,
     &                     u1,v1,t1,q1,tkh,prnum,ntke,
     &                     psk,rbsoil,zorl,u10m,v10m,fm,fh,
     &                     tsea,heat,evap,stress,spd1,kpbl,
     &                     prsi,del,prsl,prslk,phii,phil,delt,
     &                     dusfc,dvsfc,dtsfc,dqsfc,dkt,hpbl,
     &                     kinver,xkzm_m,xkzm_h,xkzm_s,xkzminv, 
     &                     xkzo,  xkzmo, 
     &                     lprnt,ipr,me)
!
      use machine  , only : kind_phys
      use funcphys , only : fpvs
      use physcons, grav => con_g,    rd => con_rd,  cp => con_cp
     &,             hvap => con_hvap, fv => con_fvirt
      implicit none
!
!     arguments
!
      logical lprnt
      integer ipr, me , ix, im, km, ntrac, ntcw, ncnd, ntke
      integer, dimension(im) ::  kinver, kpbl
!
      real(kind=kind_phys) delt, xkzm_m, xkzm_h, xkzm_s, xkzminv
      real(kind=kind_phys), dimension(im,km) :: du, dv, tau, prnum
!
      real(kind=kind_phys), dimension(im,km,ntrac) :: rtg

      real(kind=kind_phys), dimension(ix,km)   :: u1, v1, t1, tkh 
     &,                                           prsl, del, phil, prslk
      real(kind=kind_phys), dimension(ix,km+1) :: prsi, phii
      real(kind=kind_phys), dimension(ix,km,ntrac) :: q1
      real(kind=kind_phys), dimension(im)          :: psk, rbsoil, zorl
     &,                                               spd1, u10m, v10m
     &,                                               fm, fh, tsea, hpbl
     &,                                               dusfc, dvsfc
     &,                                               dtsfc, dqsfc
!
!    locals
!
      integer i,is,k,kk,km1,kmpbl,kp1, ntloc
!
      logical  pblflg(im), sfcflg(im), flg(im),  flgup(im)

      real(kind=kind_phys), dimension(im) ::  evap, heat, phih, phim
     &,                     rbdn, rbup, sflux, z0, crb, zol, thermal
     &,                     stress, beta, tx1, 
     &                      ustar, wscaleu
!
      real(kind=kind_phys), dimension(im,km)  :: theta, thvx, zl, a1, ad
     &,                                          dt2odel
      real(kind=kind_phys), dimension(im,km-1):: xkzo, xkzmo, al, au
     &,                                          dku, dkt, rdzt
!
      real(kind=kind_phys) zi(im,km+1), a2(im,km*(ntrac+1))
!
      real(kind=kind_phys) dsdz2,  dsdzq,  dsdzt, dsig, dt2, rdt
     &,                    dtodsd, dtodsu, rdz,   tem,  tem1
     &,                    ttend,  utend, vtend,  qtend, vpert
     &,                    spdk2,  rbint, ri,     zol1, robn, bvf2
!
      real(kind=kind_phys), parameter ::  gravi=1.0/grav, zolcr=0.2,
     &                      zolcru=-0.5,  rimin=-100.,    sfcfrac=0.1,
     &                      crbcon=0.25,  crbmin=0.15,    crbmax=0.35,
     &                      qmin=1.e-8,   zfmin=1.e-8,    qlmin=1.e-12,
     &                      aphi5=5.,     aphi16=16.,     f0=1.e-4
     &,                     cont=cp/grav, conq=hvap/grav, conw=1.0/grav
     &,                     dkmin=0.0,    dkmax=1000.
!    &,                     dkin=0.0,    dkmax=1000.,    xkzminv=0.3
     &,                     gocp=grav/cp, prmin=0.25,     prmax=4.0
     &,                     vk=0.4, cfac=6.5



! Option for PBL height calculation
! hpbl_original      - a bulk critical Richardson number algorithm from SHOCv1 as implemented by Moorthi 
! hpbl_troen_mahrt86 - Troen and Mahrt (1986), default GFS algorithm sans enchancement 
!                      of Pr due to countergradient transport in unstable PBL
! hpbl_vogelezang_holtslag96 -   Vogelezang and Holtslag (1996), an update of Troen and 
!                      Mahrt (1986), sans enchancement of Pr due to countergradient 
!                      transport in unstable PBL
! hpbl_tke           - do not use

      integer, parameter :: hpbl_original = 1, hpbl_troen_mahrt86 = 2, 
     &                      hpbl_vogelezang_holtslag96 = 3, hpbl_tke =4

!      integer, parameter :: hpbl_algorithm = hpbl_original
      integer, parameter :: hpbl_algorithm = hpbl_troen_mahrt86
!      integer, parameter :: hpbl_algorithm = hpbl_vogelezang_holtslag96

! Ortions for Prandtl number calculation in stable BL

! Original GFS assumption is that Pr=1 in stable BL. It's based on the fact that
! measurments of Pr in SBL show substantial scatter around unity when plotted against
! M-O stability parameter z/L. However, dependence emerges when Pr is plotted against
! gradient Richardson number  (e.g. Li, D. 2019, Atmospheric Research, 216, 86-105)
      logical, parameter :: Pr_is_unity_in_SBL = .false.
! Original GFS used Pr-Ri dependence above PBL based on Kim and Mahrt, 1992, Tellus, 44A
      logical, parameter :: PrRi_Kim_Marht_92 = .false.
! Empirical function of Venayagamoorthy and Stretch, 2010, J. of Fluid Mechanics, 644
! For consistency with M-O derived values, Pr in neutral limit is assumed to be unity. 
      logical, parameter :: PrRi_Venayagamoorthy_Stretch_2010 = .true.
! 
!
!-----------------------------------------------------------------------
!
!     compute preliminary variables
!
      if (ix < im) stop
!
!     if (lprnt) write(0,*)' in moninshoc tsea=',tsea(ipr)
      dt2   = delt
      rdt   = 1. / dt2
      km1   = km - 1
      kmpbl =  km / 2 ! Look for PBL height only in the lower half of the atmosphere
!
      do k=1,km
        do i=1,im
          zi(i,k)      = phii(i,k) * gravi ! Geometric height of layer interface
          zl(i,k)      = phil(i,k) * gravi ! Geometric height of layer center
! del(i,k) is thickness of atmosperic layer in pressure units
          dt2odel(i,k) = dt2 / del(i,k)
        enddo
      enddo
      do i=1,im
         zi(i,km+1) = phii(i,km+1) * gravi
      enddo
!
      do k = 1,km1
        do i=1,im
          rdzt(i,k)  = 1.0 / (zl(i,k+1) - zl(i,k)) ! Thickness of a layer
          prnum(i,k) = 1.0
        enddo
      enddo

!Setup backgrond diffision
      do i=1,im
        prnum(i,km) = 1.0
        tx1(i) = 1.0 / prsi(i,1)
      enddo
      do k = 1,km1
        do i=1,im
          xkzo(i,k)  = 0.0   ! background diffusivity for heat 
          xkzmo(i,k) = 0.0   ! background diffusivity for momentum
          if (k < kinver(i)) then
! Vertical background diffusivity for heat and momentum
            tem1       = 1.0 - prsi(i,k+1) * tx1(i)
            tem1       = min(1.0, exp(-tem1 * tem1 * 10.0))
            xkzo(i,k)  = xkzm_h * tem1
            xkzmo(i,k) = xkzm_m * tem1
          endif
        enddo
      enddo


!     if (lprnt) then
!       print *,' xkzo=',(xkzo(ipr,k),k=1,km1)
!       print *,' xkzmo=',(xkzmo(ipr,k),k=1,km1)
!     endif
!
!  diffusivity in the inversion layer is set to be xkzminv (m^2/s)
!
!     do k = 1,kmpbl
!       do i=1,im
!         if(zi(i,k+1) > 250.) then
!           tem1 = (t1(i,k+1)-t1(i,k)) * rdzt(i,k)
!           if(tem1 > 1.e-5) then
!              xkzo(i,k)  = min(xkzo(i,k),xkzminv)
!           endif
!         endif
!       enddo
!     enddo
!
!
      do i = 1,im
! Surface roughness length, m
         z0(i)     = 0.01 * zorl(i) ! Convert cm to m
         kpbl(i)   = 1
         hpbl(i)   = zi(i,1)
         pblflg(i) = .true.
         sfcflg(i) = .true.
! if bulk Richardson number calculated between the surface and the first 
! atmospheric level is positive then BL is stable
         if(rbsoil(i) > 0.) sfcflg(i) = .false.
         dusfc(i)  = 0.
         dvsfc(i)  = 0.
         dtsfc(i)  = 0.
         dqsfc(i)  = 0.
      enddo
!
      do k = 1,km
        do i=1,im
          tx1(i) = 0.0
        enddo
        do kk=1,ncnd
          do i=1,im
! Total condesate mixing ratio
            tx1(i) = tx1(i) + max(q1(i,k,ntcw+kk-1), qlmin)
          enddo
        enddo
        do i = 1,im
! Potential temperature
          theta(i,k) = t1(i,k) * psk(i) / prslk(i,k)
! Virtual potential temperature
          thvx(i,k)  = theta(i,k)*(1.+fv*max(q1(i,k,1),qmin)-tx1(i))
        enddo
      enddo
!
!     if (lprnt) write(0,*)' heat=',heat(ipr),' evap=',evap(ipr)
      do i = 1,im
! Buoyancy flux at the surface
         sflux(i)  = heat(i) + evap(i)*fv*theta(i,1)
!         if (heat(i) > 0.5) write(0,*)' heat=',heat(i),' evap=',evap(i) 
! If PBL is not unstable, set pblflag to .false.
         if(.not.sfcflg(i) .or. sflux(i) <= 0.) pblflg(i)=.false.
! dt/dz
         beta(i)  = dt2 / (zi(i,2)-zi(i,1))
! Surface friction velocity scale
         ustar(i) = sqrt(stress(i))
      enddo
!
!  Compute PBL height
!
!     write(0,*)' IN moninbl u10=',u10m(1:5),' v10=',v10m(1:5)

      if ( hpbl_algorithm .ne. hpbl_tke) then 
       do i=1,im
          flg(i)  = .false.
          rbup(i) = rbsoil(i)

          if(pblflg(i) .or. 
     &       hpbl_algorithm .eq. hpbl_vogelezang_holtslag96) then
! Unstable or neutral PBL treatment in Troen and Mart (1986) and
! all BL types for Vogelezang and Holtslag (1996)
! Set temperature of the thermal to virtual temperature of the first atmosperic LAYER
            thermal(i) = thvx(i,1)
! Critical bulk Ri is set to a constant (currently 0.25)
            crb(i) = crbcon
          else
! Stable PBL in Troen and Mart (1986) with modifications from Vickers and  Mahrt (2004)
! Set temperature of the thermal to surface skin virtual temperature
            thermal(i) = tsea(i)*(1.+fv*max(q1(i,1,1),qmin))
! Empirical formula for critical Ri dependence on surface Rossby number in neutral 
! BL from  Vickers, D. and L. Mahrt, 2004: Evaluating Formulations of Stable 
! Boundary Layer Height. J. Appl. Meteor., 43, 1736–1749 
! See Eq 17
! Applicable only to surface bulk Ri. 
            tem   = max(1.0, sqrt(u10m(i)*u10m(i) + v10m(i)*v10m(i)))
! Surface Rossby number Ro = |V|/(f*z0). Note that the Coriolis parameter f0 in this 
! formula is a constant and is included  only for dimensional
! reasons, and there is no evidence that the depth of the
! stable boundary layer varies with f0 (Vickers and  Mahrt, 2004)
! Therefore f0 is set to a typical in midlatitudes value of 1e-4
            robn   = tem / (f0 * z0(i)) 
            tem1   = 1.e-7 * robn
            crb(i) = max(min(0.16 * (tem1 ** (-0.18)), crbmax), crbmin)
          endif
       enddo
! Set  PBL top to level where bulk Ri is > critical value.
       do k = 1 , kmpbl
         do i = 1, im
           if(.not.flg(i)) then
             rbdn(i) = rbup(i)
             if (hpbl_algorithm .eq. hpbl_vogelezang_holtslag96) then
                spdk2   = max((u1(i,k)-u1(i,1))**2 + 
     &                    (v1(i,k)-v1(i,1))**2 + 100.*ustar(i), 1.)
! Bulk Richardson number on a current layer center
                rbup(i) = (thvx(i,k)-thermal(i))*(phil(i,k)-phil(i,1))
     &                    / (thvx(i,1)*spdk2)
             else ! Troen and Mart (1986) and Moorthi's original implementation
                spdk2   = max((u1(i,k)*u1(i,k)+v1(i,k)*v1(i,k)), 1.)
! Surface bulk Richardson number on a current layer center
                rbup(i) = (thvx(i,k)-thermal(i))*phil(i,k)
     &                    / (thvx(i,1)*spdk2)
             endif
! Layer containing PBL top
             kpbl(i) = k
             flg(i)  = rbup(i) > crb(i)
           endif
         enddo
       enddo
! Interpolate to find the exact height of the PBL top
       do i = 1,im
         if(kpbl(i) > 1) then
           k = kpbl(i)
           if(rbdn(i) >= crb(i)) then
! I think only rbdn(i) == crb(i) is possible
             rbint = 0.
           elseif(rbup(i) <= crb(i)) then
! Possible when k == kmpbl
             rbint = 1.
           else
! Assume that bulk Rb changes linearly with height between layer centers
             rbint = (crb(i)-rbdn(i)) / (rbup(i)-rbdn(i))
           endif
! Interpolate PBL height between layer centers
           hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
! If the calcualted height is below the bottom interface of the current layer,
! set index of the layer containing PBL top to kpbl-1 
           if(hpbl(i) < zi(i,kpbl(i))) kpbl(i) = kpbl(i) - 1
         else
           hpbl(i) = zl(i,1)
           kpbl(i) = 1
         endif
       enddo
      endif 


! Compute PBL height based on TKE. 
! Experimental, do not use. 
!      
!      if ( hpbl_algorithm .eq. hpbl_tke) then 
!
!       do i=1,im
!          flg(i)  = .false.
!          flgup(i)  = .false.
! !         rbup(i) = q1(i,1,ntke)
!          rbdn(i) = q1(i,1,ntke)
!          crb(i) = 0.005
!       enddo
!
!       do k = 1, kmpbl
!         do i = 1, im
!           if(.not.flg(i)) then
!!             rbdn(i) = rbup(i)
!             rbup(i) = q1(i,k,ntke)
!!             kpbl(i) = k
!             flgup(i)  = rbdn(i) >= crb(i)
!!             if (flgup(i))  kpbl(i) = k
!             flg(i)  = flgup(i) .and. (rbup(i) < crb(i)) ! .or. (zl(i,k) >= 4000)
!             if (flg(i))  kpbl(i) = k
!             rbdn(i) = rbup(i)
!           endif
!         enddo
!       enddo
!       do i = 1,im
!         if(kpbl(i) > 1) then
!           k = kpbl(i)
!           if(rbdn(i) <= crb(i) ) then
!             rbint = 0.
!           elseif(rbup(i) >= crb(i)) then
!             rbint = 1.
!           else
!             rbint = (crb(i) - rbdn(i)) / (rbup(i)-rbdn(i))
!           endif
!           hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
!           if(hpbl(i) < zi(i,kpbl(i))) kpbl(i) = kpbl(i) - 1
!         else
!           hpbl(i) = zl(i,1)
!           kpbl(i) = 1
!           pblflg(i) = .false.
!         endif
!       enddo
!      endif

!
!  Compute similarity parameters 
!
      do i=1,im
! Monin-Obukhov stability parameter z/L, where L is the M-O length, at the 
! first atmosperic level
         zol(i) = max(rbsoil(i)*fm(i)*fm(i)/fh(i),rimin)
         if(sfcflg(i)) then
           zol(i) = min(zol(i),-zfmin)
         else
           zol(i) = max(zol(i),zfmin)
         endif
! Monin-Obukhov stability parameter at the top of the surface layer (=0.1*Hpbl)
         zol1 = zol(i)*sfcfrac*hpbl(i)/zl(i,1)
! Calculate values of stability functions  at the top of the surface layer
! Unlike Troen and Mahrt (1986), here we use stability functions from Dyer and Hicks (1970)
         if(sfcflg(i)) then
! Unstable PBL
           tem     = 1.0 / max(1. - aphi16*zol1, 1.0e-8)
           phih(i) = sqrt(tem)
           phim(i) = sqrt(phih(i))
         else
! Stable PBL
! Note that the surface layer scheme in GFS uses a different form of M-O stability
! functions for the stable PBL. They are a function of gradient Ri, not the M-O
! stability parameter z/L: phim = phih= 1 + 5*Ri
! Therefore, to assume here phim dependence on z/L would be inconsistent with the 
! surface scheme. However, we only need the phih/phim ratio which is a unity under 
! either formulation 
!           phim(i) = 1. + aphi5*zol1
           phim(i) = 1.
           phih(i) = phim(i)
         endif

        
         if (pblflg(i) .and. hpbl_algorithm .ne. hpbl_original) then
           vpert=sflux(i)*(ustar(i)**3. + 
     &     aphi16*vk*sfcfrac*hpbl(i)*sflux(i)*grav/theta(i,1))**(-1./3.)
            if (hpbl_algorithm .eq. hpbl_troen_mahrt86) vpert=6.5*vpert
            if (hpbl_algorithm .eq. hpbl_vogelezang_holtslag96) 
     &                                                  vpert=8.5*vpert
            vpert=max(3.,vpert)
            thermal(i)=thermal(i)+vpert
         endif
         if ( hpbl_algorithm .eq. hpbl_tke) then
! Calculate Prandtl number at the top of the surface layer
! Assume that this value is applicable in the rest of the PBL
          if (pblflg(i)) then
! Commented out countergradient transport correction for Pr in unstable BL. 
             tem = phih(i)/phim(i) !+ cfac*vk*sfcfrac     
          else
             tem = 1. ! phih(i)/phim(i)
          endif
         prnum(i,1) = min(prmin,max(prmax,tem))
         endif 
      enddo
!
!  enhance the pbl height by considering the thermal excess
!     
      if ( hpbl_algorithm .ne. hpbl_tke) then 
      do i=1,im
         flg(i)  = .true.
         if (pblflg(i)) then
           flg(i) = .false.
           rbup(i) = rbsoil(i)
         endif
      enddo
      do k = 2, kmpbl
        do i = 1, im
          if(.not.flg(i)) then
            rbdn(i) = rbup(i)
            if (hpbl_algorithm .eq. hpbl_vogelezang_holtslag96) then
               spdk2   = max((u1(i,k)-u1(i,1))**2 + 
     &              (v1(i,k)-v1(i,1))**2 + 100.*ustar(i), 1.)
! Bulk Richardson number on a current layer center
               rbup(i) = (thvx(i,k)-thermal(i))*(phil(i,k)-phil(i,1))
     &              / (thvx(i,1)*spdk2)
            else  ! Troen and Mart (1986) and Moorthi's original implementation
               spdk2   = max((u1(i,k)*u1(i,k)+v1(i,k)*v1(i,k)), 1.)
               rbup(i) = (thvx(i,k)-thermal(i)) * phil(i,k)
     &              / (thvx(i,1)*spdk2)
            endif
            kpbl(i) = k
            flg(i)  = rbup(i) > crb(i)
          endif
        enddo
      enddo
      do i = 1,im
        if (pblflg(i)) then
          k = kpbl(i)
          if(rbdn(i) >= crb(i)) then
            rbint = 0.
          elseif(rbup(i) <= crb(i)) then
            rbint = 1.
          else
            rbint = (crb(i)-rbdn(i)) / (rbup(i)-rbdn(i))
          endif
          if (k > 1) then
            hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
            if(hpbl(i) < zi(i,kpbl(i))) kpbl(i) = kpbl(i) - 1
!            if(kpbl(i) <= 1) then
!              pblflg(i) = .false.
!            endif
          else
            pblflg(i) = .false.
          endif
        endif
! Calculate Prandtl number at the top of the surface layer
! Assume that this value is applicable in the rest of the PBL
          if (pblflg(i)) then
! Commented out countergradient transport correction for Pr in unstable BL. 
          tem = phih(i)/phim(i) + cfac*vk*sfcfrac
        else
          tem = 1. !  phih(i)/phim(i)
        endif
        prnum(i,1) = min(prmin,max(prmax,tem))
      enddo
      endif

      do i = 1, im
        if(zol(i) > zolcr) then
          kpbl(i) = 1
        endif
      enddo


!
!  Compute Prandtl number 
!
      do k = 1, km1
        kp1 = k + 1
        do i=1,im
!old
!          if(k >= kpbl(i)) then
            rdz  = rdzt(i,k)
            tem  = u1(i,k) - u1(i,kp1)
            tem1 = v1(i,k) - v1(i,kp1)
            tem  = (tem*tem + tem1*tem1) * rdz * rdz
            bvf2 = (0.5*grav)*(thvx(i,kp1)-thvx(i,k))*rdz
!     &           / (t1(i,k)+t1(i,kp1))
! Use virtual potential temperature for consistency with the buoyancy term in the 
! TKE equation 
     &           / (thvx(i,k)+thvx(i,kp1)) 
! Gradient Ri
            ri   = max(bvf2/tem,rimin)
! new
! Free atmosphere
          if(k >= kpbl(i)) then
            if(ri < 0.) then ! unstable regime
              prnum(i,kp1) = 1.0
            else

               if (PrRi_Kim_Marht_92) then  
                  prnum(i,kp1) = min(1.0 + 2.1*ri, prmax)
               elseif (PrRi_Venayagamoorthy_Stretch_2010) then 
                  prnum(i,kp1) = min(exp(-3*ri) + 4*ri, prmax)
               endif 
               
            endif
          else ! if (k > 1) then !
! new
! PBL
           if (Pr_is_unity_in_SBL .or. pblflg(i)) then  
!             if (pblflg(i)) then
                prnum(i,kp1) = prnum(i,1)
             else

                if (PrRi_Kim_Marht_92) then  
                   prnum(i,kp1) = min(1.0 + 2.1*ri, prmax)
                elseif (PrRi_Venayagamoorthy_Stretch_2010) then 
                   prnum(i,kp1) = min(exp(-3*ri) + 4*ri, prmax)
                endif 
 
             endif
! old
!            prnum(i,kp1) = prnum(i,1)
          endif
!
!         prnum(i,kp1) = 1.0
          prnum(i,kp1) = max(prmin, min(prmax, prnum(i,kp1)))
          tem      = tkh(i,kp1) * prnum(i,kp1)
! old
!          dku(i,k) = max(min(tem+xkzmo(i,k),       dkmax), xkzmo(i,k))
!          dkt(i,k) = max(min(tkh(i,kp1)+xkzo(i,k), dkmax), xkzo(i,k))
! new
          dku(i,k) = max(min(tem,        dkmax), xkzmo(i,k))
          dkt(i,k) = max(min(tkh(i,kp1), dkmax), xkzo(i,k))
          
        enddo
      enddo
!
!     compute tridiagonal matrix elements for heat and moisture
!
      do i=1,im
         ad(i,1) = 1.
         a1(i,1) = t1(i,1)   + beta(i) * heat(i)
         a2(i,1) = q1(i,1,1) + beta(i) * evap(i)
      enddo
!     if (lprnt) write(0,*)' a1=',a1(ipr,1),' beta=',beta(ipr)
!    &,' heat=',heat(ipr), ' t1=',t1(ipr,1)

      ntloc = 1
      if(ntrac > 1) then
        is    = 0
        do k = 2, ntrac
          if (k /= ntke) then
            ntloc = ntloc + 1
            is = is + km
            do i = 1, im
              a2(i,1+is) = q1(i,1,k)
            enddo
          endif
        enddo
      endif
!
      do k = 1,km1
        kp1 = k + 1
        do i = 1,im
          dtodsd    = dt2odel(i,k)
          dtodsu    = dt2odel(i,kp1)
          dsig      = prsl(i,k)-prsl(i,kp1)
          rdz       = rdzt(i,k)
          tem1      = dsig * dkt(i,k) * rdz
          dsdz2     = tem1 * rdz
          au(i,k)   = -dtodsd*dsdz2
          al(i,k)   = -dtodsu*dsdz2
!
          ad(i,k)   = ad(i,k)-au(i,k)
          ad(i,kp1) = 1.-al(i,k)
          dsdzt     = tem1 * gocp
          a1(i,k)   = a1(i,k)   + dtodsd*dsdzt
          a1(i,kp1) = t1(i,kp1) - dtodsu*dsdzt
          a2(i,kp1) = q1(i,kp1,1)
!
        enddo
      enddo
!
      if(ntrac > 1) then
        is = 0
        do kk = 2, ntrac
          if (kk /= ntke) then
            is = is + km
            do k = 1, km1
              kp1 = k + 1
              do i = 1, im
                a2(i,kp1+is) = q1(i,kp1,kk)
              enddo
            enddo
          endif
        enddo
      endif
!
!     solve tridiagonal problem for heat and moisture
!
      call tridin(im,km,ntloc,al,ad,au,a1,a2,au,a1,a2)

!
!     recover tendencies of heat and moisture
!
      do  k = 1,km
         do i = 1,im
            ttend      = (a1(i,k)-t1(i,k))   * rdt
            qtend      = (a2(i,k)-q1(i,k,1)) * rdt
!            if (ttend < 0) then 
!               write(0,*) ' ttend=',ttend,' qtend=',qtend, " k=", k
!               write(0,*)  ' xkzo=',xkzo(i,:),' xkzmo=',xkzmo(i,:)
!            endif
            tau(i,k)   = tau(i,k)   + ttend
!     if(lprnt .and. i==ipr .and. k<11) write(0,*)' tau=',tau(ipr,k)
!    &,' ttend=',ttend,' a1=',a1(ipr,k),' t1=',t1(ipr,k)
            rtg(i,k,1) = rtg(i,k,1) + qtend
            dtsfc(i)   = dtsfc(i)   + cont*del(i,k)*ttend
            dqsfc(i)   = dqsfc(i)   + conq*del(i,k)*qtend
         enddo
      enddo
      if(ntrac > 1) then
        is = 0
        do kk = 2, ntrac
          if (kk /= ntke) then
            is = is + km
            do k = 1, km 
              do i = 1, im
                qtend = (a2(i,k+is)-q1(i,k,kk))*rdt
                rtg(i,k,kk) = rtg(i,k,kk) + qtend
              enddo
            enddo
          endif
        enddo
      endif
!
!     compute tridiagonal matrix elements for momentum
!
      do i=1,im
         ad(i,1) = 1.0 + beta(i) * stress(i) / spd1(i)
         a1(i,1) = u1(i,1)
         a2(i,1) = v1(i,1)
      enddo
!
      do k = 1,km1
        kp1 = k + 1
        do i=1,im
          dtodsd    = dt2odel(i,k)
          dtodsu    = dt2odel(i,kp1)
          dsig      = prsl(i,k)-prsl(i,kp1)
          rdz       = rdzt(i,k)
          tem1      = dsig*dku(i,k)*rdz
          dsdz2     = tem1 * rdz
          au(i,k)   = -dtodsd*dsdz2
          al(i,k)   = -dtodsu*dsdz2
!
          ad(i,k)   = ad(i,k) - au(i,k)
          ad(i,kp1) = 1.0 - al(i,k)
          a1(i,kp1) = u1(i,kp1)
          a2(i,kp1) = v1(i,kp1)
!
        enddo
      enddo

      call tridi2(im,km,al,ad,au,a1,a2,au,a1,a2)
!
!     recover tendencies of momentum
!
      do k = 1,km
        do i = 1,im
          utend    = (a1(i,k)-u1(i,k))*rdt
          vtend    = (a2(i,k)-v1(i,k))*rdt
          du(i,k)  = du(i,k)  + utend
          dv(i,k)  = dv(i,k)  + vtend
          dusfc(i) = dusfc(i) + conw*del(i,k)*utend
          dvsfc(i) = dvsfc(i) + conw*del(i,k)*vtend
        enddo
      enddo
!
      if (ntke > 0) then    ! solve tridiagonal problem for momentum and tke
!
!     compute tridiagonal matrix elements for tke
!
        do i=1,im
           ad(i,1) = 1.0
           a1(i,1) = q1(i,1,ntke)
        enddo
!
        do k = 1,km1
          kp1 = k + 1
          do i=1,im
            dtodsd    = dt2odel(i,k)
            dtodsu    = dt2odel(i,kp1)
            dsig      = prsl(i,k)-prsl(i,kp1)
            rdz       = rdzt(i,k)
!            tem1      = dsig*dku(i,k)*(rdz+rdz)
            tem1      = dsig*dku(i,k)*rdz
            dsdz2     = tem1 * rdz
            au(i,k)   = -dtodsd*dsdz2
            al(i,k)   = -dtodsu*dsdz2
!
            ad(i,k)   = ad(i,k) - au(i,k)
            ad(i,kp1) = 1.0 - al(i,k)
            a1(i,kp1) = q1(i,kp1,ntke)
          enddo
        enddo

        call tridi1(im,km,al,ad,au,a1,au,a1)
!
        do k = 1, km !     recover tendencies of tke
          do i = 1, im
            qtend = (a1(i,k)-q1(i,k,ntke))*rdt
            rtg(i,k,ntke) = rtg(i,k,ntke) + qtend
          enddo
        enddo
      endif
!
      return
      end
      subroutine tridi1(l,n,cl,cm,cu,r1,au,a1)
!
      use machine     , only : kind_phys
      implicit none
      integer             k,n,l,i
      real(kind=kind_phys) fk
!
      real(kind=kind_phys) cl(l,2:n),cm(l,n),cu(l,n-1),r1(l,n),         &
     &                     au(l,n-1),a1(l,n)
!
      do i=1,l
        fk      = 1./cm(i,1)
        au(i,1) = fk*cu(i,1)
        a1(i,1) = fk*r1(i,1)
      enddo
      do k=2,n-1
        do i=1,l
          fk      = 1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k) = fk*cu(i,k)
          a1(i,k) = fk*(r1(i,k)-cl(i,k)*a1(i,k-1))
        enddo
      enddo
      do i=1,l
        fk      = 1./(cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk*(r1(i,n)-cl(i,n)*a1(i,n-1))
      enddo
      do k=n-1,1,-1
        do i=1,l
          a1(i,k) = a1(i,k)-au(i,k)*a1(i,k+1)
        enddo
      enddo
!
      return
      end
