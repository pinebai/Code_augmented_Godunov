module barotropic
! Barotropic gas.

use conserved_quantities
use my_messages

implicit none

type baro_thermo
 real(8) :: gamma, rho0, beta, p0
end type baro_thermo

type barotro
 type(conserved) :: csq
 type(baro_thermo) :: thermo
end type barotro

interface df
 module procedure density
end interface

interface mf
 module procedure momentum
end interface

interface ef
 module procedure energy
end interface

interface ief
 module procedure internal_energy
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


public  df, mf, ef, pf, cf, vf, hf, sf, ief
private density, momentum, energy, pressure, velocity, sound_speed, &
        entropy, enthalpy, internal_energy


contains


! The density taken as the first component of the conserved quantities.
function density(u) result(c)

implicit none

type(barotro), intent(in) :: u
real(8) :: c

c=u%csq%value(1)

end function density


! The momentum taken as the second component of the conserved quantities.
function momentum(u) result(c)

implicit none

type(barotro), intent(in) :: u
real(8) :: c

c=u%csq%value(2)

end function momentum


! The total energy taken as the third component of the conserved quantities.
function energy(u) result(c)

implicit none

type(barotro), intent(in) :: u
real(8) :: c

c=u%csq%value(3)

end function energy


! The pressure calculated from the conservative quantities. 
function pressure(u) result(c)

implicit none

type(barotro), intent(in) :: u
real(8) :: c

real(8) :: gamma,rho0,beta

gamma=u%thermo%gamma
rho0=u%thermo%rho0
beta=u%thermo%beta

c=(gamma-1.0d0)*(u%csq%value(3)-0.5d0*u%csq%value(2)**2.d0 &
   /u%csq%value(1)+u%csq%value(1)*beta/rho0)-gamma*beta

end function pressure


! The sound speed computed from the conservative quantities.
function sound_speed(u) result(c)

implicit none

type(barotro), intent(in):: u
real(8) :: c

real*8:: gamma,rho0,beta

gamma=u%thermo%gamma
beta=u%thermo%beta

c=dsqrt(gamma*(pf(u)+beta)/u%csq%value(1)) 

end function sound_speed


! The enthalpy calculated from the conservative quantities.
function enthalpy(u) result(c)

implicit none

type(barotro), intent(in):: u
real*8:: c
   
c=(u%csq%value(3)+pf(u))/u%csq%value(1)

end function enthalpy


! The velocity calculated from the conservative quantities.
function velocity(u) result(c)

implicit none

type(barotro), intent(in):: u
real*8:: c
 
c=u%csq%value(2)/u%csq%value(1)

end function velocity
 

! The entropy computed from the conservative quantities.
function entropy(u) result(c)
implicit none
type(barotro), intent(in):: u
real(8) :: c

real*8:: gamma,rho0,beta,p0

p0=u%thermo%p0
gamma=u%thermo%gamma 
rho0=u%thermo%rho0
beta=u%thermo%beta

c=dlog((pf(u)+beta)*(rho0/u%csq%value(1))**gamma/(p0+beta))  

end function entropy


! The internal energy computed from the conservative quantities.
function internal_energy(u) result(c)

implicit none

type(barotro), intent(in):: u
real*8:: c

c=u%csq%value(3)-0.5d0*u%csq%value(2)**2.0d0/u%csq%value(1)

end function internal_energy


! Tran3 calculates a conserved array from given density, velocity, 
! pressure, gamma, beta and rho0.
subroutine dynamic_2_conserved_barotropic(rho,v,p,u)  

implicit none

real*8, intent(in) :: rho, v, p
type(barotro), intent(out) :: u 

real(8) :: gamma,rho0,beta

gamma=u%thermo%gamma
rho0=u%thermo%rho0
beta=u%thermo%beta

u%csq%value(1)=rho
u%csq%value(2)=v*rho
u%csq%value(3)=(p+gamma*beta)/(gamma-1.0d0)- &
              rho*beta/rho0+0.5d0*rho*v**2.0d0

end subroutine dynamic_2_conserved_barotropic


end module barotropic