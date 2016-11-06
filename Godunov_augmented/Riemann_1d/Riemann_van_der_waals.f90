module Riemann_van_der_waals_Phi
! Phi function and computation density through a rarefaction.

use van_der_waals

implicit none

public  phi_van_der_waals, density_rare_van_der_waals


contains


function phi_van_der_waals(rho,p,u,p_star) 

implicit none

real(8), intent(in) :: rho, p, p_star
type(vdwaals), intent(in) :: u
real(8) :: c, rho_star
real(8) :: tao_star,x
real(8) :: gamma, a, b, epslon
real(8) :: phi_van_der_waals

epslon=1.0d-6
gamma=u%thermo%gamma
a=u%thermo%a
b=u%thermo%b

if(p_star/p.ge.(1.0d0-epslon)) then
 call root_tao_star(rho,p,gamma,a,b,p_star,tao_star)
 phi_van_der_waals=dsqrt(-(p_star-p)/(tao_star-1.0d0/rho))
else
 rho_star=root_func_rho(rho,p,gamma,a,b,p_star)
 phi_van_der_waals=-(p_star-p)/integr_func(rho,p,gamma,a,b,rho_star)
end if

end function phi_van_der_waals


function integer_chara_van_der_waals(rho,p,u)

implicit none

real(8), intent(in) :: rho, p
type(vdwaals), intent(in) :: u
real(8) :: integer_chara_van_der_waals
real(8) :: gamma, a, b

gamma=u%thermo%gamma
a=u%thermo%a
b=u%thermo%b

integer_chara_van_der_waals=integr_func(rho,p,gamma,a,b,0.0d0)

end function integer_chara_van_der_waals


function density_rare_van_der_waals(rho_side,p_side,u_side,p_star) result(c)

implicit none

real(8), intent(in) :: p_star, rho_side, p_side
type(vdwaals) :: u_side
real(8) :: c, c1, a, b
real(8) :: aa, gamma, beta, rho0

gamma=u_side%thermo%gamma
a=u_side%thermo%a
b=u_side%thermo%b

c=root_func_rho(rho_side,p_side,gamma,a,b,p_star)

end function density_rare_van_der_waals


!calculate the root of tao_star.
subroutine root_tao_star(rho,p,gamma,a,b,p_star,tao_star) 

implicit none

real*8, intent(in) :: rho, p, gamma, a, b, p_star
real*8, intent(out) :: tao_star
real*8 :: tao, epslon,x2
real*8 :: tao_a, tao_b, tao_c, tao_d
real*8 :: f_a, f_b, f_c, f_p, f_q
real*8 :: taoe, delta, pslon
complex(8) tao_star_1, tao_star_2, tao_star_3
complex(8) omega,x1

omega=(-0.5d0,0.866025403784439)
epslon=1.0d-6

tao=1/rho
tao_a=p_star/(gamma-1)+0.5d0*(p_star+p)
tao_b=-p_star*b/(gamma-1)-(p+a/(tao**2.0d0))* &
   (tao-b)/(gamma-1)+a/tao-0.5d0*(p_star+p)*tao
tao_c=a/(gamma-1)-a
tao_d=-a*b/(gamma-1)

f_a=tao_b/tao_a
f_b=tao_c/tao_a
f_c=tao_d/tao_a

f_p=-(f_a**2.0d0)/3.0d0+f_b
f_q=2.0d0*(f_a**3.0d0)/27.0d0-f_a*f_b/3.0d0+f_c

delta=(0.5d0*f_q)**2.0d0+(f_p/3.0d0)**3.0d0
x1=dsqrt(delta)
tao_star_1=(-0.5d0*f_q+dsqrt(delta))**(1.0d0/3.0d0)+(-0.5d0 &
           *f_q-dsqrt(delta))**(1.0d0/3.0d0)-f_a/3.0d0
tao_star_2=omega*(-0.5d0*f_q+dsqrt(delta))**(1.0d0/3.0d0)+ &
           (omega**2.0d0)*(-0.5d0*f_q-dsqrt(delta)) &
		   **(1.0d0/3.0d0)-f_a/3.0d0
tao_star_3=(omega**2.0d0)*(-0.5d0*f_q+dsqrt(delta)) &
           **(1.0d0/3.0d0)+omega*(-0.5d0*f_q-dsqrt(delta))** &
		   (1.0d0/3.0d0)-f_a/3.0d0

if(abs(aimag(tao_star_1)).le.epslon) then
  tao_star=tao_star_1
else if(abs(aimag(tao_star_2)).le.epslon) then
  tao_star=tao_star_2
else if(abs(aimag(tao_star_3)).le.epslon) then
  tao_star=tao_star_3
else
  call my_error_message
end if
  
x2=tao_a*(tao_star**3.0d0)+tao_b*(tao_star &
    **2.0d0)+tao_c*tao_star+tao_d

end subroutine root_tao_star


!calculate the intergal of func_state.
function integr_func(rho,p,gamm,a,b,rho_star)

implicit none

real*8, intent(in) :: rho,p,gamm,b,a,rho_star
real*8 :: x1
real*8 :: errabs, errest, error, errrel
real*8 :: pa, pb
real*8 :: func_state, integr_func
integer  irule
external func_state

pa=rho_star
pb=rho
errabs=1.0d-12
errrel=1.0d-12
irule=1

call dqdag(f,pa,pb,errabs,errrel,irule,x1,errest)
integr_func=x1

contains

function f(x)
real*8 :: f, x
f=func_state(rho,p,gamm,a,b,x)
end function f

end function integr_func


!calculate the root of function func_rho.
function root_func_rho(rho,p,gamm,a,b,p_star) 

implicit none

real*8, intent(in) :: rho,p,gamm,a,b,p_star
integer :: itmax,nroot
real*8 func_rho, root_func_rho
real*8 :: eps,errabs,errrel,eta,x1,xxx
parameter (nroot=1)
integer :: info(nroot)
real*8 :: x(nroot),xguess(nroot)
external func_rho

data xguess/30/

eps=1.0d-6
errabs=1.0d-6
errrel=1.0d-6
eta=1.0d-6
itmax=200

call dzreal(f,errabs,errrel,eps,eta,nroot,itmax,xguess,x1,info)
root_func_rho=x1
xxx=f(root_func_rho)

contains

function f(x)
real*8 :: f,x
f=func_rho(rho,p,gamm,a,b,p_star,x)
end function f

end function root_func_rho


!"x" in function func_f1 is rho_star.
function func_rho(rho,p,gamm,a,b,p_star,x)

implicit none

real*8, intent(in) :: rho,p,gamm,a,b,p_star,x
real*8:: func_rho

func_rho=(x*(1.0d0-b*rho)/(rho*(1.0d0-b*x))) &
    **gamm*(p+a*rho**2.0d0)-a*x**2.0d0-p_star

end function func_rho


!'x' in function func_state is variable 'rho' for integre.
function func_state(rho,p,gamm,a,b,x)

implicit none

real*8, intent(in) :: rho,p,gamm,a,b,x
real*8:: func_state

func_state=dsqrt((x*(1.0d0-b*rho)/(rho*(1.0d0-b*x)))**gamm &
          *(p+a*rho**2.0d0)*gamm/(x*(1.0d0-b*x))-2.0d0*a*x)/x

end function func_state


end module Riemann_van_der_waals_Phi