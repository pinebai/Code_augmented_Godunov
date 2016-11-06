module average

use solution
! 'solution.f90'

use physics
!' physics.f90'

use RiemannSolver
! 'RiemannSolver.f90'

implicit none


contains
!***************************************************************************************


!***************************************************************************************
subroutine average_vandp

implicit none

integer :: i,IW
real*8 :: RS,US,PS
real*8 :: vl_aver,vr_aver,pl_aver,pr_aver,entropyl,entropyr,uul_aver,uur_aver
real*8 :: rarefaction_rho_l,rarefaction_rho_r,shock_rho_l,shock_rho_r,contact_rho_l,contact_rho_r
real*8,dimension(-ighost+1:cells_number+ighost) :: temp_v_aver,temp_p_aver
real*8,dimension(-ighost+1:cells_number+ighost) :: R1,U1,P1,R2,U2,P2

real(8) :: u_star, p_star

do i=-ighost+1, cells_number+ighost 

 R1(i)=u(i,1)-d_rho(i)
 U1(i)=v_aver(i)-d_v(i)
 P1(i)=p_aver(i)-d_p(i)
   
 R2(i)=u(i,1)+d_rho(i)
 U2(i)=v_aver(i)+d_v(i)
 P2(i)=p_aver(i)+d_p(i)

enddo

temp_v_aver=v_aver; temp_p_aver=p_aver
v_aver=0.0d0; p_aver=0.0d0; entropy_increase=0.0d0

do i=1, cells_number

!***************************************************************************************
 rarefaction_rho_l=0.0d0;rarefaction_rho_r=0.0d0;shock_rho_l=0.0d0;shock_rho_r=0.0d0;contact_rho_l=0.0d0;contact_rho_r=0.0d0
 vl_aver=0.0d0;vr_aver=0.0d0;pl_aver=0.0d0;pr_aver=0.0d0;uul_aver=0.0d0;uur_aver=0.0d0;entropyl=0.0d0;entropyr=0.0d0
   
 CALL RIEMANN0(R2(i-1),U2(i-1),P2(i-1),R1(i),U1(i),P1(i),u_star,p_star,IW)
 if(iw.lt.4) then
  call aver_velocity_pressure0(vl_aver,vr_aver,pl_aver,pr_aver,entropyl,entropyr)    ! No vacuum.
  v_aver(i)=v_aver(i)+vr_aver
  p_aver(i)=p_aver(i)+pr_aver
  entropy_increase(i)=entropy_increase(i)+entropyr
 else
  call aver_velocity_pressure1(vl_aver,vr_aver,pl_aver,pr_aver)                      ! Vacuum.
  v_aver(i)=v_aver(i)+vr_aver
  p_aver(i)=p_aver(i)+pr_aver
 end if 

! print*,i,"зѓ",iw
!***************************************************************************************
 rarefaction_rho_l=0.0d0;rarefaction_rho_r=0.0d0;shock_rho_l=0.0d0;shock_rho_r=0.0d0;contact_rho_l=0.0d0;contact_rho_r=0.0d0
 vl_aver=0.0d0;vr_aver=0.0d0;pl_aver=0.0d0;pr_aver=0.0d0;uul_aver=0.0d0;uur_aver=0.0d0;entropyl=0.0d0;entropyr=0.0d0
   
 CALL RIEMANN0(R1(i),U1(i),P1(i),R2(i),U2(i),P2(i),u_star,p_star,IW)
 if(iw.lt.4) then
  call aver_velocity_pressure0(vl_aver,vr_aver,pl_aver,pr_aver,entropyl,entropyr)    ! No vacuum.
  v_aver(i)=v_aver(i)+vl_aver+vr_aver
  p_aver(i)=p_aver(i)+pl_aver+pr_aver
  entropy_increase(i)=entropy_increase(i)+entropyl+entropyr 
 else
  call aver_velocity_pressure1(vl_aver,vr_aver,pl_aver,pr_aver)                      ! Vacuum.
  v_aver(i)=v_aver(i)+vl_aver+vr_aver
  p_aver(i)=p_aver(i)+pl_aver+pr_aver
 end if 


! print*,i,"жа",iw
!***************************************************************************************
 rarefaction_rho_l=0.0d0;rarefaction_rho_r=0.0d0;shock_rho_l=0.0d0;shock_rho_r=0.0d0;contact_rho_l=0.0d0;contact_rho_r=0.0d0
 vl_aver=0.0d0;vr_aver=0.0d0;pl_aver=0.0d0;pr_aver=0.0d0;uul_aver=0.0d0;uur_aver=0.0d0;entropyl=0.0d0;entropyr=0.0d0
    
 CALL RIEMANN0(R2(i),U2(i),P2(i),R1(i+1),U1(i+1),P1(i+1),u_star,p_star,IW)
 if(iw.lt.4) then
  call aver_velocity_pressure0(vl_aver,vr_aver,pl_aver,pr_aver,entropyl,entropyr)    ! No vacuum.
  v_aver(i)=v_aver(i)+vl_aver
  p_aver(i)=p_aver(i)+pl_aver
  entropy_increase(i)=entropy_increase(i)+entropyl
 else
  call aver_velocity_pressure1(vl_aver,vr_aver,pl_aver,pr_aver)                      ! Vacuum.
  v_aver(i)=v_aver(i)+vl_aver
  p_aver(i)=p_aver(i)+pl_aver
 end if 

enddo

end subroutine average_vandp
!***************************************************************************************


!***************************************************************************************                                                  
subroutine aver_velocity_pressure0(vl_aver,vr_aver,pl_aver,pr_aver,entropyl,entropyr)
! Compute velocity and pressure integrals and entropy increase across shocks for non-vacuum case.

implicit none

real*8, intent(out) :: vl_aver, vr_aver, pl_aver, pr_aver, entropyl, entropyr

integer :: i
real*8 :: G1, G2, G3, G4, G5, G6, G7, G8, G9
real*8 :: DL, UL, PL, CL, DR, UR, PR, CR, PM, UM
real*8 :: SML, SMR
real*8 :: C,SHL,CML,STL,PML,SL,PMR,SR,SHR,CMR,STR
real*8 :: intl_rho,intl_v,intl_p,intll_rho_sonic,intll_v_sonic,intll_p_sonic,intlr_rho_sonic,intlr_v_sonic,intlr_p_sonic
real*8 :: intr_rho,intr_v,intr_p,intrl_rho_sonic,intrl_v_sonic,intrl_p_sonic,intrr_rho_sonic,intrr_v_sonic,intrr_p_sonic
real*8 :: rhoml,rhomr
real*8, dimension(-1:1) :: intp_v,intp_p 

COMMON /GAMMAS/ G1, G2, G3, G4, G5, G6, G7, G8, G9
COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR, PM, UM

entropyl=0.0d0; entropyr=0.0d0

!shock

sl=ul-cl*dsqrt(g2*pm/pl+g1) ! Left shock speed
sr=ur+cr*dsqrt(g2*pm/pr+g1) ! Right shock speed

!rarefaction

cml=cl*(pm/pl)**g1 ! Left middle sound speed
cmr=cr*(pm/pr)**g1 ! Right middle sound speed

shl=ul-cl  ! Left characteristic speed
stl=um-cml ! left middle characteristic speed

shr=ur+cr  ! Right characteristic speed
str=um+cmr ! Right middle characteristic speed

!The density between the shock and contact discontinuity
rhoml=dl*((pm/pl+g6)/(pm/pl*g6+1.0d0))
rhomr=dr*((pm/pr+g6)/(pm/pr*g6+1.0d0))
    
call compute_rarefaction_left(shl,stl,intl_rho,intl_v,intl_p)
call compute_rarefaction_left(shl,0.0d0,intll_rho_sonic,intll_v_sonic,intll_p_sonic)
call compute_rarefaction_left(0.0d0,stl,intlr_rho_sonic,intlr_v_sonic,intlr_p_sonic)
call compute_rarefaction_right(str,shr,intr_rho,intr_v,intr_p)
call compute_rarefaction_right(str,0.0d0,intrl_rho_sonic,intrl_v_sonic,intrl_p_sonic)
call compute_rarefaction_right(0.0d0,shr,intrr_rho_sonic,intrr_v_sonic,intrr_p_sonic)
    
if(um>0) then
! The contact discontinuity goes to the right.
     
 if(pm>pl) then
! The left wave is a shock.   
    
  if(sl>0.0d0) then
! The left shock goes to the right; therefore, the left half of the Riemann solution is occupied by the initial left state.
    
   vl_aver=0.25d0*ul                                      
   pl_aver=0.25d0*pl
   entropyr=entropyr+shock_entropy_increase(dl,ul,pl,rhoml,um,pm,sl)
      
   if(pm>pr) then
! The right wave is a shock and goes to the right.
      
    entropyr=entropyr+shock_entropy_increase(rhomr,um,pm,dr,ur,pr,sr)
    vr_aver=sl*dt/dx*ul+(sr-sl)*dt/dx*um+(0.25d0-sr*dt/dx)*ur
    pr_aver=sl*dt/dx*pl+(sr-sl)*dt/dx*pm+(0.25d0-sr*dt/dx)*pr
     
   else
! The right wave is a rarefaction wave and goes to the right.
    vr_aver=sl*dt/dx*ul+(um+cmr-sl)*dt/dx*um+(0.25d0-(ur+cr)*dt/dx)*ur+intr_v
    pr_aver=sl*dt/dx*pl+(um+cmr-sl)*dt/dx*pm+(0.25d0-(ur+cr)*dt/dx)*pr+intr_p
    
   endif
    
  else
! The left shock goes to the left.
   
   entropyl=entropyl+shock_entropy_increase(dl,ul,pl,rhoml,um,pm,sl)
   vl_aver=(0.25d0-dabs(sl)*dt/dx)*ul+dabs(sl)*dt/dx*um
   pl_aver=(0.25d0-dabs(sl)*dt/dx)*pl+dabs(sl)*dt/dx*pm
     
   if(pm>pr) then
! The right wave is a shock and goes to the right.
     
    entropyr=entropyr+shock_entropy_increase(rhomr,um,pm,dr,ur,pr,sr)
    vr_aver=sr*dt/dx*um+(0.25d0-sr*dt/dx)*ur
    pr_aver=sr*dt/dx*pm+(0.25d0-sr*dt/dx)*pr  
    
   else
! The right wave is a rarefaction wave and goes to the right.
    
    vr_aver=(um+cmr)*dt/dx*um+(0.25d0-(ur+cr)*dt/dx)*ur+intr_v
    pr_aver=(um+cmr)*dt/dx*pm+(0.25d0-(ur+cr)*dt/dx)*pr+intr_p
       
   endif
  endif
    
 else
! The left wave is a rarefaction wave.
    
  if((ul-cl)>0.0d0) then
! The left rarefaction wave goes to the right; therefore, the left half of the Riemann solution is occupied by the initial left state.
   
   vl_aver=0.25d0*ul
   pl_aver=0.25d0*pl
     
   if(pm>pr) then
! The right wave is a shock and goes to the right.
      
    entropyr=entropyr+shock_entropy_increase(rhomr,um,pm,dr,ur,pr,sr)
    vr_aver=(ul-cl)*dt/dx*ul+(sr-(um-cml))*dt/dx*um+(0.25d0-sr*dt/dx)*ur+intl_v
    pr_aver=(ul-cl)*dt/dx*pl+(sr-(um-cml))*dt/dx*pm+(0.25d0-sr*dt/dx)*pr+intl_p         
      
   else
! The right wave is a rarefaction wave.
         
    vr_aver=(ul-cl)*dt/dx*ul+((um+cmr)-(um-cml))*dt/dx*um+(0.25d0-(ur+cr)*dt/dx)*ur+intl_v+intr_v
    pr_aver=(ul-cl)*dt/dx*pl+((um+cmr)-(um-cml))*dt/dx*pm+(0.25d0-(ur+cr)*dt/dx)*pr+intl_p+intr_p       		    
          
   endif
           
  else
! The left edge of the left rarefaction wave goes to the left.
        
   if((um-cml)<0.0d0) then
! The right edge of the left rarefaction wave goes to the left; therefore, the whole left rarefaction wave goes to the left.
! Because the contact discontinuity now goes to the right, the left half of the Riemann solution consists of the left state, 
! the left rarefaction wave and part of the left middle state. 
     
    vl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*ul+dabs(um-cml)*dt/dx*um+intl_v 
    pl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*pl+dabs(um-cml)*dt/dx*pm+intl_p
      		 
    if(pm>pr) then
! The right wave is a shock. Since the contact discontinuity now goes to the right, the shock also goes to the right.
     
     entropyr=entropyr+shock_entropy_increase(rhomr,um,pm,dr,ur,pr,sr)
     vr_aver=sr*dt/dx*um+(0.25d0-sr*dt/dx)*ur
     pr_aver=sr*dt/dx*pm+(0.25d0-sr*dt/dx)*pr
    
    else
! The right wave is a rarefaction wave. Since the contact discontinuity now goes to the right, the rarefaction wave also goes to the right.
       
     vr_aver=(um+cmr)*dt/dx*um+(0.25d0-(ur+cr)*dt/dx)*ur+intr_v
     pr_aver=(um+cmr)*dt/dx*pm+(0.25d0-(ur+cr)*dt/dx)*pr+intr_p
      
    endif 
       
   else
!  The right edge of the left rarefaction wave goes to the right; therefore, the left rarefaction wave is sonic.
    
    vl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*ul+intll_v_sonic
    pl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*pl+intll_p_sonic
     
    if(pm>pr) then
! The right wave is a shock and goes to the right.
    
     entropyr=entropyr+shock_entropy_increase(rhomr,um,pm,dr,ur,pr,sr)
     vr_aver=(sr-dabs(um-cml))*dt/dx*um+(0.25d0-sr*dt/dx)*ur+intlr_v_sonic
     pr_aver=(sr-dabs(um-cml))*dt/dx*pm+(0.25d0-sr*dt/dx)*pr+intlr_p_sonic
    
    else
! The right wave is a rarefaction wave and goes to the right.
       
     vr_aver=((um+cmr)-(um-cml))*dt/dx*um+(0.25d0-(ur+cr)*dt/dx)*ur+intlr_v_sonic+intr_v
     pr_aver=((um+cmr)-(um-cml))*dt/dx*pm+(0.25d0-(ur+cr)*dt/dx)*pr+intlr_p_sonic+intr_p	
     
    endif
   endif
  endif
 end if
    
else
! The contact discontinuity goes to the left.
         
 if(pm>pr) then
! The right wave is a shock.
    
  if(sr<0.0d0) then
! The right shock goes to the left; therefore, the right half of the Riemann solution is occupied by the initial right state.
     
   vr_aver=0.25d0*ur    
   pr_aver=0.25d0*pr 
   entropyl=entropyl+shock_entropy_increase(rhomr,um,pm,dr,ur,pr,sr)
    
   if(pm>pl) then
! The left wave is a shock and goes to the left.
              
    entropyl=entropyl+shock_entropy_increase(dl,ul,pl,rhoml,um,pm,sl)
    vl_aver=(0.25d0-dabs(sl)*dt/dx)*ul+(dabs(sl)-dabs(sr))*dt/dx*um+dabs(sr)*dt/dx*ur 
    pl_aver=(0.25d0-dabs(sl)*dt/dx)*pl+(dabs(sl)-dabs(sr))*dt/dx*pm+dabs(sr)*dt/dx*pr 
     
   else
! The left wave is rarefaction wave and goes to the left.       
    vl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*ul+(dabs(um-cml)-dabs(sr))*dt/dx*um+dabs(sr)*dt/dx*ur+intl_v
    pl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*pl+(dabs(um-cml)-dabs(sr))*dt/dx*pm+dabs(sr)*dt/dx*pr+intl_p 
    
   endif
    
  else
! The right shock goes to the right
       
   entropyr=entropyr+shock_entropy_increase(rhomr,um,pm,dr,ur,pr,sr)
   vr_aver=sr*dt/dx*um+(0.25d0-sr*dt/dx)*ur    
   pr_aver=sr*dt/dx*pm+(0.25d0-sr*dt/dx)*pr
           
   if(pm>pl) then
! The left wave is a shock and goes to the left.
    
    entropyl=entropyl+shock_entropy_increase(dl,ul,pl,rhoml,um,pm,sl)
    vl_aver=(0.25d0-dabs(sl)*dt/dx)*ul+dabs(sl)*dt/dx*um 
    pl_aver=(0.25d0-dabs(sl)*dt/dx)*pl+dabs(sl)*dt/dx*pm
    
   else
! The left wave is a rarefaction wave and goes to the left.    
    
    vl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*ul+dabs(um-cml)*dt/dx*um+intl_v
    pl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*pl+dabs(um-cml)*dt/dx*pm+intl_p
    
   endif
  endif
   
 else
! The right wave is a rarefaction wave.
     
  if((ur+cr)<0.0d0) then
! The right rarefaction wave goes to the left as a whole; there, the right half of the Riemann solution is occupied by the initial right state.   
      
   vr_aver=0.25d0*ur    
   pr_aver=0.25d0*pr 
    
   if(pm>pl) then
! The left wave is a shock and goes to the left.    
    
    vl_aver=(0.25d0-dabs(sl)*dt/dx)*ul+(dabs(sl)-dabs(um+cmr))*dt/dx*um+dabs(ur+cr)*dt/dx*ur+intr_v
    pl_aver=(0.25d0-dabs(sl)*dt/dx)*pl+(dabs(sl)-dabs(um+cmr))*dt/dx*pm+dabs(ur+cr)*dt/dx*pr+intr_p	    
     
   else
! The left wave is a rarefaction wave and goes to the left.
     
    vl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*ul+(dabs(um-cml)-dabs(um+cmr))*dt/dx*um+dabs(ur+cr)*dt/dx*ur+intl_v+intr_v
    pl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*pl+(dabs(um-cml)-dabs(um+cmr))*dt/dx*pm+dabs(ur+cr)*dt/dx*pr+intl_p+intr_p
    
   endif	     
    
  else
! The right edge of the right rarefaction wave goes to the right.
    
   if((um+cmr)>0.0d0) then
! The left edge of the right rarefaction wave goes to the right; therefore, the whole right rarefaction wave goes to the right.
! Because the contact discontinuity now goes to the left, the right half of the Riemann solution consists of the rightt state, 
! the right rarefaction wave and part of the right middle state.         
    
    vr_aver=(um+cmr)*dt/dx*um+(0.25d0-(ur+cr)*dt/dx)*ur+intr_v  
    pr_aver=(um+cmr)*dt/dx*pm+(0.25d0-(ur+cr)*dt/dx)*pr+intr_p
       
    if(pm>pl) then
! The left wave is a shock and goes to the left.
      
     entropyl=entropyl+shock_entropy_increase(dl,ul,pl,rhoml,um,pm,sl)
     vl_aver=(0.25d0-dabs(sl)*dt/dx)*ul+dabs(sl)*dt/dx*um
     pl_aver=(0.25d0-dabs(sl)*dt/dx)*pl+dabs(sl)*dt/dx*pm
      
    else
! The left wave is a rarefaction wave and goes to the left.      
      
     vl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*ul+dabs(um-cml)*dt/dx*um+intl_v
     pl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*pl+dabs(um-cml)*dt/dx*pm+intl_p
        
    endif
     		
   else
! The left edge of the right rarefaction wave goes to the left; therefore, the right rarefaction wave is sonic. 
         
    vr_aver=(0.25d0-(ur+cr)*dt/dx)*ur+intrr_v_sonic    
    pr_aver=(0.25d0-(ur+cr)*dt/dx)*pr+intrr_p_sonic
     
    if(pm>pl) then
! The left wave is a shock and goes to the left.    
     
     entropyl=entropyl+shock_entropy_increase(dl,ul,pl,rhoml,um,pm,sl)
     vl_aver=(0.25d0-dabs(sl)*dt/dx)*ul+(dabs(sl)-dabs(um+cmr))*dt/dx*um+intrl_v_sonic   
     pl_aver=(0.25d0-dabs(sl)*dt/dx)*pl+(dabs(sl)-dabs(um+cmr))*dt/dx*pm+intrl_p_sonic   		    
      
    else
! The left wave is a rarefaction wave and goes to the left.     
      
     vl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*ul+(dabs(um-cml)-dabs(um+cmr))*dt/dx*um+intrl_v_sonic+intl_v   
     pl_aver=(0.25d0-dabs(ul-cl)*dt/dx)*pl+(dabs(um-cml)-dabs(um+cmr))*dt/dx*pm+intrl_p_sonic+intl_p    	
    
    endif  	      
   endif
  endif
 endif
end if

end  subroutine aver_velocity_pressure0
!***************************************************************************************


subroutine aver_velocity_pressure1(vl_aver,vr_aver,pl_aver,pr_aver)
! Compute velocity and pressure integrals and entropy increase across shocks for vacuum case.

implicit none

real*8, intent(out) :: vl_aver, vr_aver, pl_aver, pr_aver

integer :: i
real*8 :: G1, G2, G3, G4, G5, G6, G7, G8, G9
real*8 :: DL, UL, PL, CL, DR, UR, PR, CR, PM, UM
real*8 :: SML, SMR, chl, chr, ddd, vvv, ppp

COMMON /GAMMAS/ G1, G2, G3, G4, G5, G6, G7, G8, G9
COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR, PM, UM

SML = UL + 2.0D0*CL/G9; SMR = UR - 2.0D0*CR/G9
chl=ul-cl; chr=ur+cr 

vl_aver=0.0d0; vr_aver=0.0d0; pl_aver=0.0d0; pr_aver=0.0d0

! Computation in the left half cell.

if(chl.le.0.0d0) then
! The left edge of the left rarefaction wave goes to the left.
 vl_aver=vl_aver+(0.25d0+chl*dt/dx)*ul; pl_aver=pl_aver+(0.25d0+chl*dt/dx)*pl
   
 if(sml.le.0.0d0)then
! The right edge of the left rarefaction wave, also the left bound of the vacuum, goes to the left.
  call compute_rarefaction_left(chl*dt/dx,sml*dt/dx,ddd,vvv,ppp)
  vl_aver=vl_aver+vvv; pl_aver=pl_aver+ppp 
   
  if(smr.le.0.0d0) then
! The left edge of the right rarefaction wave, also the right bound of the vacuum, goes to the left.
   vl_aver=vl_aver+(smr-sml)*dt/dx*0.5d0*(sml+smr) ! Pressure is zero in vacuum.
       
   if(chr.le.0.0d0) then
! The right edge of the right rarefaction wave goes to the left.
    call compute_rarefaction_right(smr*dt/dx,chr*dt/dx,ddd,vvv,ppp)
    vl_aver=vl_aver+vvv-chr*dt/dx*ur; pl_aver=pl_aver+ppp-chr*dt/dx*pr     
   else
! The right edge of the right rarefaction wave goes to the right.
    call compute_rarefaction_right(smr*dt/dx,0.0d0,ddd,vvv,ppp)
    vl_aver=vl_aver+vvv; pl_aver=pl_aver+ppp
   end if
  
  else
! The left edge of the right rarefaction wave, also the right bound of the vacuum, goes to the right.
   vl_aver=vl_aver-sml*dt/dx*0.5d0*(sml+smr) ! Pressure is zero in vacuum.      
  end if   	         
   
 else
! The right edge of the left rarefaction wave, also the left bound of the vacuum, goes to the right.
  call compute_rarefaction_left(chl*dt/dx,0.0d0,ddd,vvv,ppp)
  vl_aver=vl_aver+vvv; pl_aver=pl_aver+ppp  
 end if
   
else
! The left edge of the left rarefaction wave goes to the right.
 vl_aver=vl_aver+0.25d0*ul; pl_aver=pl_aver+0.25d0*pl       
end if
   
! Computation in the right half cell.

if(chr.ge.0.0d0) then
! The right edge of the right rarefaction wave goes to the right.
 vr_aver=vr_aver+(0.25d0-chr*dt/dx)*ur; pr_aver=pr_aver+(0.25d0-chr*dt/dx)*pr
   
 if(smr.ge.0.0d0)then
! The left edge of the right rarefaction wave, also the right bound of the vacuum, goes to the right.
  call compute_rarefaction_right(smr*dt/dx,chr*dt/dx,ddd,vvv,ppp)
  vr_aver=vr_aver+vvv; pr_aver=pr_aver+ppp 
   
  if(sml.ge.0.0d0) then
! The right edge of the left rarefaction wave, also the left bound of the vacuum, goes to the right.
   vr_aver=vr_aver+(smr-sml)*dt/dx*0.5d0*(sml+smr) ! Pressure is zero in vacuum.
       
   if(chl.ge.0.0d0) then
! The left edge of the left rarefaction wave goes to the right.
    call compute_rarefaction_left(chl*dt/dx,sml*dt/dx,ddd,vvv,ppp)
    vr_aver=vr_aver+vvv+chl*dt/dx*ul; pr_aver=pr_aver+ppp+chl*dt/dx*pl     
   else
! The right edge of the right rarefaction wave goes to the right.
    call compute_rarefaction_left(0.0d0,sml*dt/dx,ddd,vvv,ppp)
    vr_aver=vr_aver+vvv; pr_aver=pr_aver+ppp
   end if
  
  else
! The right edge of the left rarefaction wave, also the left bound of the vacuum, goes to the left.
   vr_aver=vr_aver+smr*dt/dx*0.5d0*(sml+smr) ! Pressure is zero in vacuum.      
  end if   	         
   
 else
! The left edge of the right rarefaction wave, also the right bound of the vacuum, goes to the left.
  call compute_rarefaction_right(0.0d0,chr*dt/dx,ddd,vvv,ppp)
  vr_aver=vr_aver+vvv; pr_aver=pr_aver+ppp  
 end if
   
else
! The right edge of the right rarefaction wave goes to the left.
 vr_aver=vr_aver+0.25d0*ur; pr_aver=pr_aver+0.25d0*pr       
end if   
   
end subroutine aver_velocity_pressure1


!***************************************************************************************      
subroutine compute_rarefaction_left(a,b,value_intl_rho,value_intl_v,value_intl_p) 

integer :: i
real*8 :: G1, G2, G3, G4, G5, G6, G7, G8, G9
real*8 :: DL, UL, PL, CL, DR, UR, PR, CR, PM, UM
real*8 :: a,b,value_intl_rho,value_intl_v,value_intl_p
!real*8,dimension(-1:1) :: gp,intp_rho,intp_v,intp_p

COMMON /GAMMAS/ G1, G2, G3, G4, G5, G6, G7, G8, G9
COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR, PM, UM

value_intl_rho=dt/dx*dl*cl*((g5+g6/cl*(ul-a))**(1.0d0/g6)-(g5+g6/cl*(ul-b))**(1.0d0/g6))
value_intl_v=dt/dx*g5*0.5d0*((cl+g7*ul+b)**2.0d0-(cl+g7*ul+a)**2.0d0)
value_intl_p=dt/dx*(gamma+1.0d0)/(3.0d0*gamma-1.0d0)*pl*cl*((g5+g6/cl*(ul-a))**(g3+1.0d0)-(g5+g6/cl*(ul-b))**(g3+1.0d0))

end subroutine compute_rarefaction_left
!***************************************************************************************


!***************************************************************************************
subroutine compute_rarefaction_right(a,b,value_intr_rho,value_intr_v,value_intr_p) 

integer :: i
real*8 :: G1, G2, G3, G4, G5, G6, G7, G8, G9
real*8 :: DL, UL, PL, CL, DR, UR, PR, CR, PM, UM
real*8 :: a,b,value_intr_rho,value_intr_v,value_intr_p
real*8,dimension(-1:1) :: gp,intp_rho,intp_v,intp_p

COMMON /GAMMAS/ G1, G2, G3, G4, G5, G6, G7, G8, G9
COMMON /STATES/ DL, UL, PL, CL, DR, UR, PR, CR, PM, UM

value_intr_rho=dt/dx*dr*cr*((g5-g6/cr*(ur-b))**(1.0d0/g6)-(g5-g6/cr*(ur-a))**(1.0d0/g6))
value_intr_v=dt/dx*g5*0.5d0*((-cr+g7*ur+b)**2.0d0-(-cr+g7*ur+a)**2.0d0)
value_intr_p=dt/dx*(gamma+1.0d0)/(3.0d0*gamma-1.0d0)*pr*cr*((g5-g6/cr*(ur-b))**(g3+1.0d0)-(g5-g6/cr*(ur-a))**(g3+1.0d0))

end subroutine compute_rarefaction_right
!***************************************************************************************


!***************************************************************************************  
end module average