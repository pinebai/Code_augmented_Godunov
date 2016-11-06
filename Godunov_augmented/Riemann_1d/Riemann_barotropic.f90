module Riemann_barotropic_Phi
! Phi function and computation density through a rarefaction.

use barotropic

implicit none

public  phi_barotropic, density_rare_barotropic


contains


function phi_barotropic(rho,p,u,p_star) 

implicit none

real(8), intent(in) :: rho, p, p_star
type(barotro), intent(in) :: u
real(8) :: phi_barotropic
real(8) :: gamma, beta, epslon

epslon=1.0d-6
gamma=u%thermo%gamma
beta=u%thermo%beta

if((p_star+beta)/(p+beta).ge.1.0d0-epslon)then
 phi_barotropic=dsqrt((p+beta)*rho)*dsqrt((gamma+1.d0)*(p_star &
                +beta)/(2.0d0*(p+beta))+(gamma-1.d0)/2.0d0)
else
 phi_barotropic=(p_star-p)*(gamma-1.0d0)/(2.0d0*(((p_star+beta) &
                /(p+beta))**((gamma-1.d0)/(2.0d0*gamma))-1.0d0) &
				*dsqrt(gamma*(p+beta)/rho))
end if

end function phi_barotropic


function integer_chara_barotropic(rho,p,u)

implicit none

real(8), intent(in) :: rho, p
type(barotro), intent(in) :: u
real(8) :: integer_chara_barotropic
real(8) :: gamma, beta

gamma=u%thermo%gamma
beta=u%thermo%beta

integer_chara_barotropic=2.0d0*dsqrt(gamma*(p+beta)/rho)/(gamma-1.0d0)

end function integer_chara_barotropic


function density_rare_barotropic(rho_side,p_side,u_side,p_star) result(c)

implicit none

real(8), intent(in) :: p_star, rho_side, p_side
type(barotro), intent(in) :: u_side
real(8) :: c
real(8) :: aa, gamma, beta, rho0, p0

gamma=u_side%thermo%gamma
beta=u_side%thermo%beta
rho0=u_side%thermo%rho0
p0=u_side%thermo%p0

aa=(p_side+beta)*(rho0/rho_side)**gamma/(p0+beta)
c=((p_star+beta)*rho0**gamma/(aa*(p0+beta)))**(1.d0/gamma)

end function density_rare_barotropic


end module Riemann_barotropic_Phi