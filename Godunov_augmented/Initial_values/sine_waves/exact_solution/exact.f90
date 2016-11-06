program exact_sine_wave

implicit none

real(8), dimension(0:2000) :: x, density, velocity, pressure, entropy_per_volume
real(8) :: time, pi, gamma
integer :: i

print*, 'Input time'
read(*,'(f20.8)') time

pi=4.0d0*datan(1.0d0); gamma=1.4d0

do i=0, 2000
 x(i)=dfloat(i)*0.0005d0
 density(i)=1.0d0+0.2d0*dsin(2.0d0*pi*(x(i)+time)) 
 velocity(i)=1.0d0; pressure(i)=1.0d0
 entropy_per_volume(i)=density(i)*dlog(pressure(i)/density(i)**gamma)
end do

open(1,file='d:\Godunov_augmented\show\exact\exact_solution.dat')
do i=0,2000
 write(1,'(5f20.8)') x(i), density(i), velocity(i), pressure(i), entropy_per_volume(i)
end do   	      
close(1)   

end program exact_sine_wave