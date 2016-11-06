program convergence_check

implicit none

real(8), dimension(:,:), allocatable :: numerical, exact
integer :: cells_number, i
real(8) :: accuracy_1, accuracy_m, xxx

print*, 'Input the cells number.'
read(*,*) cells_number

allocate(numerical(1:cells_number,3)); allocate(exact(1:cells_number,3))

open(2,file='d:\Godunov_augmented\input\solutions.dat')
do i=1, cells_number
 read(2,'(3f25.16)') numerical(i,1), numerical(i,2), numerical(i,3)
end do  
close(2)

open(3,file='d:\Godunov_augmented\output\solutions.dat')
do i=1, cells_number
 read(3,'(3f25.16)') exact(i,1), exact(i,2), exact(i,3)
end do  
close(3)

accuracy_1=0.0d0; accuracy_m=0.0d0
do i=cells_number/5+1, cells_number*7/10
 xxx=dabs(numerical(i,1)-exact(i,1))
 accuracy_1=accuracy_1+xxx
 if(accuracy_m.lt.xxx) accuracy_m=xxx
end do
accuracy_1=accuracy_1/dfloat(cells_number)

print*, accuracy_1
print*
print*, accuracy_m

end program convergence_check