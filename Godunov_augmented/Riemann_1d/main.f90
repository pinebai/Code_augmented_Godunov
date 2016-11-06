
Program Main

use Riemann_solver

implicit none

real*8::pl, pr, rhol, rhor, vl, vr
real*8::pl_star, rhol_star, vl_star
real*8::pr_star, rhor_star, vr_star
character*11, dimension(3) :: wv_pt
logical :: vacuum
real(8) :: vacuum_velocity_l, vacuum_velocity_r

type(phys)::ur, ul, ur_star, ul_star, ull, urr, xxx

type(phys) :: fr_star, fr, fl, fl_star
real(8) :: sksp, ent

real(8) :: rhox, ux, px, sx, sxl, sxr, sr
real(8) :: aaa

real(8), dimension(0:2000) :: x, density, velocity, pressure, entropy_per_volume
integer :: i
real(8) :: left_left_wave_edge, right_left_wave_edge, left_right_wave_edge, right_right_wave_edge
real(8) :: entropy_l, entropy_r, entropy_star_l, entropy_star_r, sound_speed, sound_speed_l, sound_speed_r, gamma_l, gamma_r
real(8) :: time, shift

left_left_wave_edge=error_data; left_right_wave_edge=error_data
right_left_wave_edge=error_data; right_right_wave_edge=error_data
entropy_l=error_data; entropy_r=error_data; sound_speed_l=error_data; sound_speed_r=error_data
gamma_l=error_data; gamma_r=error_data

!The kind of gas on the left.

ull%type_of_fluid='Polytropic          '
nullify(ull%vdw); nullify(ull%bar); allocate(ull%pol)

!ull%type_of_fluid='Barotropic          '
!nullify(ull%vdw); nullify(ull%pol); allocate(ull%bar)

!ull%type_of_fluid='Van Der Waals       '
!nullify(ull%pol); nullify(ull%bar); allocate(ull%vdw)

!The kind of gas on the right.

urr%type_of_fluid='Polytropic          '
nullify(urr%vdw); nullify(urr%bar); allocate(urr%pol)

!urr%type_of_fluid='Barotropic          '
!nullify(urr%vdw); nullify(urr%pol); allocate(urr%bar)

!urr%type_of_fluid='Van Der Waals       '
!nullify(urr%bar); nullify(urr%pol); allocate(urr%vdw)
  

!rhol=1.1d0; vl=0.0d0; pl=2000.0d0; ull%bar%thermo%gamma=7.0d0; 
!ull%bar%thermo%beta=3000.0d0; ull%bar%thermo%rho0=1.0d0;
!ull%bar%thermo%p0=1.0d0
!rhor=1.1d0; vr=0.0d0; pr=1000.0d0; urr%bar%thermo%gamma=7.0d0; 
!urr%bar%thermo%beta=3000.0d0; urr%bar%thermo%rho0=1.0d0;
!urr%bar%thermo%p0=1.0d0

!rhol=5.0d1; vl=1.0d3; pl=1.0d+5; ull%vdw%thermo%gamma=1.4d0; 
!ull%vdw%thermo%a=5.0d0; ull%vdw%thermo%b=1.0d-3
!rhor=5.0d1; vr=0.0d0; pr=1.0d+5; urr%vdw%thermo%gamma=1.4d0; 
!urr%vdw%thermo%a=5.0d0; urr%vdw%thermo%b=1.0d-3

!rhol=2.27550882502480d0; vl=0.781011006283816d0; pl=2.56167433147329d0; ull%pol%thermo%gamma=1.4d0
!rhor=0.167d0; vr=0.0d0; pr=1.013d0; urr%pol%thermo%gamma=1.63d0; 

!rhol=5.0d+1; vl=0.0d0; pl=1.0d+5; ull%vdw%thermo%gamma=1.4d0; 
!ull%vdw%thermo%a=5.0d0; ull%vdw%thermo%b=1.0d-3
!rhor=1.0d0; vr=0.0d0; pr=1.0d0; urr%pol%thermo%gamma=1.4d0

!rhol=1.0d0; vl=6.3386d0; pl=1000.0d0; ull%bar%thermo%gamma=7.0d0; 
!ull%bar%thermo%beta=3000.0d0; ull%bar%thermo%rho0=1.0d0;
!ull%bar%thermo%p0=1.0d0
!rhor=5.0d+1; vr=0.0d0; pr=1.0d+5; urr%vdw%thermo%gamma=1.4d0; 
!urr%vdw%thermo%a=5.0d0; urr%vdw%thermo%b=1.0d-3

!rhol=5.0d+1; vl=0.0d0; pl=1.0d+5; ull%vdw%thermo%gamma=1.4d0; 
!ull%vdw%thermo%a=5.0d0; ull%vdw%thermo%b=1.0d-3
!rhor=1.0d0; vr=6.3386d0; pr=1000.0d0; urr%bar%thermo%gamma=7.0d0; 
!urr%bar%thermo%beta=3000.0d0; urr%bar%thermo%rho0=1.0d0;
!urr%bar%thermo%p0=1.0d0

!rhol=1.0d0; vl=6.3386d0; pl=1000.0d0; ull%bar%thermo%gamma=7.0d0; 
!ull%bar%thermo%beta=3000.0d0; ull%bar%thermo%rho0=1.0d0;
!ull%bar%thermo%p0=1.0d0
!rhor=0.125d0; vr=0.0d0; pr=0.1d0; urr%pol%thermo%gamma=1.2d0

!rhol=0.344568d0; vl=1.52872d0; pl=2.46610d0; ull%pol%thermo%gamma=1.4d0
!rhor=0.5d0; vr=0.0d0; pr=0.571d0; urr%pol%thermo%gamma=1.4d0
!rr%bar%thermo%beta=3000.0d0; urr%bar%thermo%rho0=1.0d0;
!urr%bar%thermo%p0=1.0d0

! The Sod's shock tube problem for gamma-law gas.
rhor=0.125d0; vr=0.0d0; pr=1.0d-1; urr%pol%thermo%gamma=1.4d0
rhol=1.0d0; vl=0.0d0; pl=1.0d0; ull%pol%thermo%gamma=1.4d0

!! The Sod's shock tube problem for gamma-law gas, left side right.
!rhol=0.125d0; vl=0.0d0; pl=1.0d-1; ull%pol%thermo%gamma=1.4d0
!rhor=1.0d0; vr=0.0d0; pr=1.0d0; urr%pol%thermo%gamma=1.4d0

! The left blast-wave Riemann problem.
!rhor=1.0d-1; vr=0.0d0; pr=1.0d-1; urr%pol%thermo%gamma=1.4d0
!rhol=1.0d0; vl=0.0d0; pl=1.0d0; ull%pol%thermo%gamma=1.4d0

!! Pure rarefaction wave with sonic point.
!rhol=1.0d0; vl=0.0d0; pl=2.0d0; ull%pol%thermo%gamma=1.4d0
!rhor=0.227490453638516d0; vr=2.14443337755679d0; pr=0.251639285328357d0; urr%pol%thermo%gamma=1.4d0

!rhol=0.299695567037000d0; vl=1.78477833876520d0; pl=0.375441565532103d0; urr%pol%thermo%gamma=1.4d0
!rhor=0.227490453639000d0; vr=2.12078000542392d0; pr=0.253649058379402d0; ull%pol%thermo%gamma=1.4d0

!! The modified Sod's shock tube problem for gamma-law gas.
!rhol=1.0d0; vl=0.75d0; pl=1.0d0; ull%pol%thermo%gamma=1.4d0
!rhor=0.125d0; vr=0.0d0; pr=1.0d-1; urr%pol%thermo%gamma=1.4d0

!rhol=115.0d0; vl=1000.0d0; pl=23084402.17d0; urr%pol%thermo%gamma=1.4d0
!rhor=1.290d-3; vr=0.0d0; pr=1.0d0; ull%pol%thermo%gamma=1.4d0

!rhol=19.237d0; vl=2000.0d0; pl=1.0d0; urr%pol%thermo%gamma=1.4d0
!rhor=19.237d0; vr=0.0d0; pr=1.0d0; ull%pol%thermo%gamma=1.4d0

! The Sod's shock tube problem for gamma-law gas.
!rhol=1.0d0; vl=1.0d0; pl=1.0d0; bl=0.0d0; gamml=1.4d0; rho0l=1.0d0
!rhor=0.125d0; vr=0.0d0; pr=1.0d0; br=0.0d0; gammr=1.4d0; rho0r=1.0d0

! The Sod's shock tube problem for gamma-law gas.
!rhol=1.0d0; vl=0.0d0; pl=1.0d0; bl=0.0d0; gamml=1.4d0; rho0l=1.0d0
!rhor=1.0d0; vr=0.0d0; pr=1.0d0; br=0.0d0; gammr=1.4d0; rho0r=1.0d0

!! The right-left Sod's shock tube problem for gamma-law gas.
!rhol=1.0d0; vl=0.0d0; pl=1.0d0; ull%pol%thermo%gamma=1.4d0
!rhor=0.125d0; vr=0.0d0; pr=1.0d-1; urr%pol%thermo%gamma=1.4d0

! The Sod's shock tube problem for gamma-law gas.
!rhol=1.0d0; vl=6.3386d0; pl=1000.0d0; ull%bar%thermo%gamma=7.0d0; 
!ull%bar%thermo%beta=3000.0d0; ull%bar%thermo%rho0=1.0d0
!rhor=1.0d0; vr=0.0d0; pr=1.0d0; gammr=7.0d0; br=3000.0d0; rho0r=1.0d0

! The Sod's shock tube problem for gamma-law gas.
!rhol=0.001d0; vl=0.0d0; pl=1000.0d0; gamml=1.4d0; bl=0.0d0; rho0l=1.0d0
!rhor=1.0d0; vr=0.0d0; pr=1.0d0; gammr=7.0d0; br=3000.0d0; rho0r=1.0d0

!! The Lax's shock tube problem for gamma-law gas.
!rhol=0.445d0; vl=0.698d0; pl=3.528d0; ull%pol%thermo%gamma=1.4d0
!rhor=0.5d0; vr=0.0d0; pr=0.571d0; urr%pol%thermo%gamma=1.4d0

!! The Tang-Liu shock tube problem for gamma-law gas.
!rhol=1000.0d0; vl=0.0d0; pl=1000.0d0; ull%pol%thermo%gamma=1.4d0
!rhor=1.0d0; vr=0.0d0; pr=1.0d0; urr%pol%thermo%gamma=1.4d0

!! The Vacuum tube problem for gamma-law gas.
!rhol=7.0d0; vl=-1.0d0; pl=0.2d0; ull%pol%thermo%gamma=1.4d0
!rhor=7.0d0; vr=1.0d0; pr=0.2d0; urr%pol%thermo%gamma=1.4d0

call phys_clon(ul,ull)
call tran3(rho=rhol,v=vl,p=pl,u=ul)
call phys_clon(ur,urr)
call tran3(rho=rhor,v=vr,p=pr,u=ur)

call riemann(ul,ur,ul_star,ur_star,wv_pt,vacuum)

if(.not.vacuum) then
  
 call tran5(rhol_star,vl_star,pl_star,ul_star)
 call tran5(rhor_star,vr_star,pr_star,ur_star)
  
! aaa=-dlog(0.01d0)
  
  fr_star=f(ur_star) 
  fr=f(ur)
  
  sksp=(df(fr_star)-df(fr))/(df(ur_star)-df(ur))
  sksp=(mf(fr_star)-mf(fr))/(mf(ur_star)-mf(ur))
  sksp=(ef(fr_star)-ef(fr))/(ef(ur_star)-ef(ur))
   
  fl=f(ul)  
  fl_star=f(ul_star)
   
! xxx%type_of_fluid='Polytropic          '
! nullify(xxx%vdw); nullify(xxx%bar); allocate(xxx%pol)
! xxx%pol%thermo%gamma=1.4d0

! xxx%pol%csq%value(1)=ur%pol%csq%value(1)-0.2d0*(fr%pol%csq%value(1)-fl%pol%csq%value(1))
! xxx%pol%csq%value(2)=ur%pol%csq%value(1)-0.2d0*(fr%pol%csq%value(2)-fl%pol%csq%value(2))
! xxx%pol%csq%value(3)=ur%pol%csq%value(3)-0.2d0*(fr%pol%csq%value(3)-fl%pol%csq%value(3))

! rhox=df(xxx) 
! ux=vf(xxx) 
! px=pf(xxx)
! sx=sf(xxx)
 
! sxl=sf(ul_star); sxr=sf(ur_star); sr=sf(ur)

! sksp=(df(fl_star)-df(fr))/(df(ul)-df(ul_star))
! sksp=(mf(fl_star)-mf(fr))/(mf(ul)-mf(ul_star))
! sksp=(ef(fl_star)-ef(fr))/(ef(ul)-ef(ul_star))

! ent=sf(ul)
! ent=sf(ul_star)

! aaa=(px/dexp(2.272d0))**(1.0d0/1.4d0)
! aaa=(px/dexp(0.409d0))**(1.0d0/1.4d0)

!aaa=vl_star-cf(ul_star)

! aaa=dlog(831.601d0)
! aaa=dlog(1.615d0)
  
 write(*,*)
  
 write(*,*) 'The left middle state is'
 write(*,'(3f25.14)') rhol_star,vl_star,pl_star
 write(*,*) 'The right middle state is'
 write(*,'(3f25.14)') rhor_star,vr_star,pr_star
  
! ent=dlog(460.0d0)
! ent=(dlog(5.4d0)-dlog(4.4d0))*1.4d0
! ent=(dlog(4.4d0)-dlog(1.7d0))*1.4d0
 
 write(*,*) 'Input the final time'
 read(*,'(f18.12)') time
 write(*,*) 'Inpute the position shift'
 read(*,'(f18.12)') shift

! Compute the left and right wave edges.
 if(pl_star.lt.pl) then
  left_left_wave_edge=(vl-cf(ul))*time; left_right_wave_edge=(vl_star-cf(ul_star))*time 
  entropy_l=pl/rhol**ull%pol%thermo%gamma; gamma_l=ull%pol%thermo%gamma; sound_speed_l=cf(ul)
  entropy_star_l=entropy_l
 else
  fl=f(ul); fl_star=f(ul_star)
  left_left_wave_edge=(df(fl_star)-df(fl))/(df(ul_star)-df(ul))*time; left_right_wave_edge=left_left_wave_edge
  entropy_l=pl/rhol**ull%pol%thermo%gamma; gamma_l=ull%pol%thermo%gamma; sound_speed_l=cf(ul)
  entropy_star_l=pl_star/rhol_star**ull%pol%thermo%gamma 
 end if
 if(pl_star.lt.pr) then
  right_right_wave_edge=(vr+cf(ur))*time; right_left_wave_edge=(vr_star+cf(ur_star))*time
  entropy_r=pr/rhor**urr%pol%thermo%gamma; gamma_r=urr%pol%thermo%gamma; sound_speed_r=cf(ur)
  entropy_star_r=entropy_r
 else
  fr=f(ur); fr_star=f(ur_star)
  right_right_wave_edge=(df(fr_star)-df(fr))/(df(ur_star)-df(ur))*time; right_left_wave_edge=right_right_wave_edge
  entropy_r=pr/rhor**urr%pol%thermo%gamma; gamma_r=urr%pol%thermo%gamma; sound_speed_r=cf(ur)
  entropy_star_r=pr_star/rhor_star**urr%pol%thermo%gamma 
 end if
 if(left_left_wave_edge.gt.left_right_wave_edge) call my_error_message
 if(left_right_wave_edge.gt.vl_star*time) call my_error_message
 if(vr_star*time.gt.right_left_wave_edge) call my_error_message
 if(right_left_wave_edge.gt.right_right_wave_edge) call my_error_message

 do i=0, 2000
  x(i)=dfloat(i)*0.5d-3
  if(x(i)-shift.le.left_left_wave_edge) then
   density(i)=rhol; velocity(i)=vl; pressure(i)=pl; entropy_per_volume(i)=rhol*dlog(entropy_l)
  else
   if(x(i)-shift.lt.left_right_wave_edge) then
    velocity(i)=(gamma_l-1.0d0)/(gamma_l+1.0d0)*(vl+2.0d0*sound_speed_l/(gamma_l-1.0d0)+2.0d0/(gamma_l-1.0d0)*(x(i)-shift)/time)
    sound_speed=(gamma_l-1.0d0)/(gamma_l+1.0d0)*(vl-(x(i)-shift)/time+2.0d0/(gamma_l-1.0d0)*sound_speed_l)
    density(i)=(sound_speed*sound_speed/gamma_l/entropy_l)**(1.0d0/(gamma_l-1.0d0))
    pressure(i)=entropy_l*density(i)**gamma_l
    entropy_per_volume(i)=density(i)*dlog(entropy_l)
   else
    if(x(i)-shift.lt.vl_star*time) then
     density(i)=rhol_star; velocity(i)=vl_star; pressure(i)=pl_star; entropy_per_volume(i)=density(i)*dlog(entropy_star_l)
    else
     if(x(i)-shift.le.right_left_wave_edge) then
      density(i)=rhor_star; velocity(i)=vl_star; pressure(i)=pr_star; entropy_per_volume(i)=density(i)*dlog(entropy_star_r)
     else
      if(x(i)-shift.lt.right_right_wave_edge) then
       velocity(i)=(gamma_r-1.0d0)/(gamma_r+1.0d0)*(vr-2.0d0*sound_speed_r/(gamma_r-1.0d0)+2.0d0/(gamma_r-1.0d0)*(x(i)-shift)/time)
       sound_speed=(gamma_r-1.0d0)/(gamma_r+1.0d0)*(-vr+(x(i)-shift)/time+2.0d0/(gamma_r-1.0d0)*sound_speed_r)
       density(i)=(sound_speed*sound_speed/gamma_r/entropy_r)**(1.0d0/(gamma_r-1.0d0))
       pressure(i)=entropy_r*density(i)**gamma_r 
       entropy_per_volume(i)=density(i)*dlog(entropy_r)
      else
       density(i)=rhor; velocity(i)=vr; pressure(i)=pr; entropy_per_volume(i)=density(i)*dlog(entropy_r)
      endif
     end if
    end if
   end if
  end if
 end do    

else
! At present the vacuum part is only for the $\gamma$-law gas.

 write(*,*) 'Input the final time'
 read(*,'(f18.12)') time
 write(*,*) 'Inpute the position shift'
 read(*,'(f18.12)') shift

 vacuum_velocity_l=vl+2.0d0*cf(ul)/(ull%pol%thermo%gamma-1.0d0)
 left_left_wave_edge=(vl-cf(ul))*time; left_right_wave_edge=vacuum_velocity_l*time 
 entropy_l=pl/rhol**ull%pol%thermo%gamma; gamma_l=ull%pol%thermo%gamma; sound_speed_l=cf(ul)
  
 vacuum_velocity_r=vr_star-2.0d0*cf(ur)/(ull%pol%thermo%gamma-1.0d0)
 right_right_wave_edge=(vr+cf(ur))*time; right_left_wave_edge=vacuum_velocity_r*time
 entropy_r=pr/rhor**urr%pol%thermo%gamma; gamma_r=urr%pol%thermo%gamma; sound_speed_r=cf(ur)
  
 do i=0, 2000
  x(i)=dfloat(i)*0.5d-3
  if(x(i)-shift.lt.left_left_wave_edge) then
   density(i)=rhol; velocity(i)=vl; pressure(i)=pl; entropy_per_volume(i)=rhol*dlog(entropy_l)
  else
   if(x(i)-shift.lt.left_right_wave_edge) then
    velocity(i)=(gamma_l-1.0d0)/(gamma_l+1.0d0)*(vl+2.0d0*sound_speed_l/(gamma_l-1.0d0)+2.0d0/(gamma_l-1.0d0)*(x(i)-shift)/time)
    sound_speed=(gamma_l-1.0d0)/(gamma_l+1.0d0)*(vl-(x(i)-shift)/time+2.0d0/(gamma_l-1.0d0)*sound_speed_l)
    density(i)=(sound_speed*sound_speed/gamma_l/entropy_l)**(1.0d0/(gamma_l-1.0d0))
    pressure(i)=entropy_l*density(i)**gamma_l
    entropy_per_volume(i)=density(i)*dlog(entropy_l)
   else
    if(x(i)-shift.le.right_left_wave_edge) then
	 density(i)=0.0d0; velocity(i)=0.5d0*(vacuum_velocity_l+vacuum_velocity_r); pressure(i)=0.0d0; entropy_per_volume(i)=0.0d0
    else
     if(x(i)-shift.lt.right_right_wave_edge) then
      velocity(i)=(gamma_r-1.0d0)/(gamma_r+1.0d0)*(vr-2.0d0*sound_speed_r/(gamma_r-1.0d0)+2.0d0/(gamma_r-1.0d0)*(x(i)-shift)/time)
      sound_speed=(gamma_r-1.0d0)/(gamma_r+1.0d0)*(-vr+(x(i)-shift)/time+2.0d0/(gamma_r-1.0d0)*sound_speed_r)
      density(i)=(sound_speed*sound_speed/gamma_r/entropy_r)**(1.0d0/(gamma_r-1.0d0))
      pressure(i)=entropy_r*density(i)**gamma_r 
      entropy_per_volume(i)=density(i)*dlog(entropy_r)
     else
      density(i)=rhor; velocity(i)=vr; pressure(i)=pr; entropy_per_volume(i)=density(i)*dlog(entropy_r)
     end if
    end if
   end if
  end if
 end do    	 

end if

open(1,file='d:\Godunov_augmented\show\exact\exact_solution.dat')
do i=0,2000
 write(1,'(5f20.8)') x(i), density(i), velocity(i), pressure(i), entropy_per_volume(i)
end do   	      
close(1)   

End program Main
 
