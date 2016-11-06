module boundary_conditions

use infinitive_boundary
! 'infinitive_boundary.f90'

use periodic_boundary
! 'periodic_boundary.f90'

use reflective_boundary
! 'reflective_boundary.f90'

implicit none

public  boundary_cell_average, boundary_HS

private infinitive_boundary_cell, infinitive_boundary_HS, periodic_boundary_cell, periodic_boundary_HS


contains


subroutine boundary_cell_average

implicit none

select case(boundary_type)

 case('infinitive  '); call infinitive_boundary_cell

 case('periodic    '); call periodic_boundary_cell

 case('reflective  '); call reflective_boundary_cell

end select

end subroutine boundary_cell_average


subroutine boundary_HS

implicit none

select case(boundary_type)

 case('infinitive  '); call infinitive_boundary_HS

 case('periodic    '); call periodic_boundary_HS

 case('reflective  '); call reflective_boundary_HS

end select

end subroutine boundary_HS


end module boundary_conditions