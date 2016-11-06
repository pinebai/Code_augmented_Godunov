module fluids

use Riemann_polytropic_Phi
use Riemann_barotropic_Phi
use Riemann_van_der_waals_Phi

implicit none

type phys
 character*20 :: type_of_fluid
 type(vdwaals), pointer :: vdw
 type(barotro), pointer :: bar
 type(polytro), pointer :: pol
end type phys

interface df
 module procedure density
end interface

interface mf
 module procedure momentum
end interface

interface ef
 module procedure energy
end interface

interface pf
 module procedure pressure
end interface

interface cf
 module procedure sound_speed
end interface
  
interface vf
 module procedure velocity
end interface

interface hf
 module procedure enthalpy
end interface

interface sf
 module procedure entropy
end interface

interface ief
 module procedure internal_energy
end interface

interface tran3
 module procedure dynamic_2_conserved
end interface

interface tran5
 module procedure conserved_2_dynamic
end interface

public  df, mf, ef, pf, cf, vf, hf, sf, ief, tran3, tran5
private density, momentum, energy, pressure, velocity, sound_speed, &
        entropy, enthalpy, internal_energy, dynamic_2_conserved, &
		conserved_2_dynamic, dynamic_2_conserved_polytropic, &
		dynamic_2_conserved_van_der_waals, density_rare_polytropic, &
        dynamic_2_conserved_barotropic, density_rare_barotropic, &
        density_rare_van_der_waals


contains


! Interfacing density calculation.
function density(u) result(c)

implicit none

type(phys), intent(in) :: u
real(8) :: c

select case(u%type_of_fluid)
 case('Van Der Waals       '); c=df(u%vdw)
 case('Polytropic          '); c=df(u%pol)
 case('Barotropic          '); c=df(u%bar)
 case default; call my_error_message
end select

end function density


! Interfacing momentum calculation.
function momentum(u) result(c)

implicit none

type(phys), intent(in) :: u
real(8) :: c

select case(u%type_of_fluid)
 case('Van Der Waals       '); c=mf(u%vdw)
 case('Polytropic          '); c=mf(u%pol)
 case('Barotropic          '); c=mf(u%bar)
 case default; call my_error_message
end select

end function momentum


! Interfacing energy calculation.
function energy(u) result(c)

implicit none

type(phys), intent(in) :: u
real(8) :: c

select case(u%type_of_fluid)
 case('Van Der Waals       '); c=ef(u%vdw)
 case('Polytropic          '); c=ef(u%pol)
 case('Barotropic          '); c=ef(u%bar)
 case default; call my_error_message
end select

end function energy


! Interfacing pressure calculation. 
function pressure(u) result(c)

implicit none

type(phys), intent(in) :: u
real*8 :: c

select case(u%type_of_fluid)
 case('Van Der Waals       '); c=pf(u%vdw)
 case('Polytropic          '); c=pf(u%pol)
 case('Barotropic          '); c=pf(u%bar)
 case default; call my_error_message
end select

end function pressure


! Interfacing sound speed calculation.
function sound_speed(u) result(c)

implicit none

type(phys), intent(in):: u
real*8:: c

select case(u%type_of_fluid)
 case('Van Der Waals       '); c=cf(u%vdw)
 case('Polytropic          '); c=cf(u%pol)
 case('Barotropic          '); c=cf(u%bar)
 case default; call my_error_message
end select

end function sound_speed


! Interfacing enthalpy calculation.
function enthalpy(u) result(c)

implicit none

type(phys), intent(in):: u
real*8 :: c
   
select case(u%type_of_fluid)
 case('Van Der Waals       '); c=hf(u%vdw)
 case('Polytropic          '); c=hf(u%pol)
 case('Barotropic          '); c=hf(u%bar)
 case default; call my_error_message
end select

end function enthalpy


! Interfacing velocity calculation.
function velocity(u) result(c)

implicit none

type(phys), intent(in):: u
real*8:: c
 
select case(u%type_of_fluid)
 case('Van Der Waals       '); c=vf(u%vdw)
 case('Polytropic          '); c=vf(u%pol)
 case('Barotropic          '); c=vf(u%bar)
 case default; call my_error_message
end select

end function velocity
 

! Interfacing entropy calculation.
function entropy(u) result(c)

implicit none

type(phys), intent(in):: u
real*8 :: c

select case(u%type_of_fluid)
 case('Van Der Waals       '); c=sf(u%vdw)
 case('Polytropic          '); c=sf(u%pol)
 case('Barotropic          '); c=sf(u%bar)
 case default; call my_error_message
end select  

end function entropy


! Interfacing internal_energy calculation.
function internal_energy(u) result(c)

implicit none

type(phys), intent(in) :: u
real(8) :: c

select case(u%type_of_fluid)
 case('Van Der Waals       '); c=ief(u%vdw)
 case('Polytropic          '); c=ief(u%pol)
 case('Barotropic          '); c=ief(u%bar)
 case default; call my_error_message
end select

end function internal_energy


! Interfacing tran3 subroutine.
subroutine dynamic_2_conserved(rho,v,p,u) 

implicit none

type(phys), intent(inout) :: u
real(8) :: rho, v, p

select case(u%type_of_fluid)
 case('Van Der Waals       ')
  call dynamic_2_conserved_van_der_waals(rho,v,p,u%vdw)
 case('Polytropic          ') 
  call dynamic_2_conserved_polytropic(rho,v,p,u%pol)
 case('Barotropic          ') 
  call dynamic_2_conserved_barotropic(rho,v,p,u%bar)
 case default; call my_error_message
end select

end subroutine dynamic_2_conserved


! Interfacing Phi function.
function phi(rho_side,p_side,u_side,p_star) result(c)

implicit none

type(phys), intent(in) :: u_side
real(8), intent(in) :: rho_side, p_side, p_star
real(8) :: c

select case(u_side%type_of_fluid)

 case('Polytropic          ')
  c=phi_polytropic(rho_side,p_side,u_side%pol,p_star)
 case('Barotropic          ')
  c=phi_barotropic(rho_side,p_side,u_side%bar,p_star)
 case('Van Der Waals       ')
  c=phi_van_der_waals(rho_side,p_side,u_side%vdw,p_star)
 case default; call my_error_message
end select

end function phi


! Interfacing Phi function.
function integer_chara(rho_side,p_side,u_side) result(c)

implicit none

type(phys), intent(in) :: u_side
real(8), intent(in) :: rho_side, p_side
real(8) :: c

select case(u_side%type_of_fluid)

 case('Polytropic          ')
  c=integer_chara_polytropic(rho_side,p_side,u_side%pol)
 case('Barotropic          ')
  c=integer_chara_barotropic(rho_side,p_side,u_side%bar)
 case('Van Der Waals       ')
  c=integer_chara_van_der_waals(rho_side,p_side,u_side%vdw)
 case default; call my_error_message
end select

end function integer_chara

! Interfacing density computation through rarefaction.
function density_rare(rho_side,p_side,u_side,p_star) result(c)

implicit none

type(phys), intent(in) :: u_side
real(8), intent(in) :: p_star, rho_side, p_side
real(8) :: c

select case(u_side%type_of_fluid)

 case('Polytropic          ')
  c=density_rare_polytropic(rho_side,p_side,u_side%pol,p_star)
 case('Barotropic          ')
  c=density_rare_barotropic(rho_side,p_side,u_side%bar,p_star)
 case('Van Der Waals       ')
  c=density_rare_van_der_waals(rho_side,p_side,u_side%vdw,p_star)
 case default; call my_error_message
end select

end function density_rare
 
  
! Tran5 calculates density, velocity, pressure, gamma, b and rho0 from a 
! given conserved array.
subroutine conserved_2_dynamic(rho,v,p,u)  

implicit none

type(phys), intent(in) :: u
real*8, intent(out) :: rho, v, p

rho=df(u)
v=vf(u)
p=pf(u)				   

end subroutine conserved_2_dynamic


subroutine phys_clon(u,u_sample)

implicit none

type(phys), intent(inout) :: u
type(phys), intent(in) :: u_sample

u%type_of_fluid=u_sample%type_of_fluid

select case(u%type_of_fluid)
 case('Van Der Waals       ')
  nullify(u%pol); nullify(u%bar)
  allocate(u%vdw)
  u%vdw%thermo=u_sample%vdw%thermo    
 case('Polytropic          ')
  nullify(u%vdw); nullify(u%bar)
  allocate(u%pol)
  u%pol%thermo=u_sample%pol%thermo
 case('Barotropic          ')
  nullify(u%pol); nullify(u%vdw)
  allocate(u%bar)
  u%bar%thermo=u_sample%bar%thermo
 case default; call my_error_message
end select

end subroutine phys_clon


!Compute the exact eigenvalue of the derivative function 
!of flux ie. the matrix A(u).
subroutine char_val(u,lambd)

implicit none

type(phys), intent(in):: u    
real*8, dimension(3), intent(out):: lambd
   
lambd(1)=vf(u)-cf(u)
lambd(2)=vf(u)
lambd(3)=vf(u)+cf(u)
  
end subroutine char_val


! Compute the pth eigenvalue of the derivative function,
!  ie. the matrix A(u)(u is a vector).
function lambd_p(u,p) result(lam_p)

implicit none

type(phys), intent(in):: u
integer, intent(in):: p 
real(8), dimension(3) :: lambd
real*8:: lam_p

call char_val(u,lambd)
lam_p=lambd(p)

end function lambd_p


! Flux function f of Euler_equation.
! u is a conserved vector of the flow, u(1)=density,u(2)=density*velocity
! u(3)=density*energy.
function f(u)

implicit none

type(phys), intent(in) :: u
type(phys) :: f

call phys_clon(f,u)

select case(u%type_of_fluid)
 case('Van Der Waals       ')
  f%vdw%csq%value(1)=mf(u)
  f%vdw%csq%value(2)=pf(u)+mf(u)**2.d0/df(u)
  f%vdw%csq%value(3)=vf(u)*(ef(u)+pf(u))
 case('Polytropic          ')
  f%pol%csq%value(1)=mf(u)
  f%pol%csq%value(2)=pf(u)+mf(u)**2.d0/df(u)
  f%pol%csq%value(3)=vf(u)*(ef(u)+pf(u))
 case('Barotropic          ')
  f%bar%csq%value(1)=mf(u)
  f%bar%csq%value(2)=pf(u)+mf(u)**2.d0/df(u)
  f%bar%csq%value(3)=vf(u)*(ef(u)+pf(u))
 case default; call my_error_message
end select

end function f
	 	

end module fluids