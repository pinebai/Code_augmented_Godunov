module conserved_quantities
! The basic data structure of the fluid, the three conservative quantities,
! density, momentum and energy.

use my_messages

implicit none

type conserved
 real*8, dimension(3):: value
end type conserved

interface operator(+)
 module procedure add
end interface

interface operator(-)
 module procedure ded
end interface

private add, ded
  

contains


function add(u,v)
! Addition operation.
 
implicit none

type(conserved), intent(in) :: u,v
type(conserved) :: add

add%value=u%value+v%value

end function add


function ded(u,v)
! Deduction operation.

implicit none 

type(conserved), intent(in) :: u,v
type(conserved) :: ded

ded%value=u%value-v%value

end function ded


end module conserved_quantities