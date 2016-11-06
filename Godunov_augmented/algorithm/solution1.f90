module solution
! The definition of numerical solution. It includes the original quantities, augmented quantities, half-steps. 
! The numerical fluxes are also defined.

use grid_and_parameters

real(8), dimension(:,:), allocatable :: u
! Numerical solution, conservative variables, mass, momentum and energy.

real(8), dimension(:), allocatable :: uu, v_aver, p_aver
! Numerical solution, augmented variables, entropy per unit volume, velocity cell-average and pressure cell-average.

real(8), dimension(:), allocatable :: d_rho, d_v, d_p
! Half steps for mass, velocity and pressure.

real(8), dimension(:), allocatable :: entropy_increase 
! Entropy increase in cells where there are shocks.

real(8), dimension(:), allocatable :: d_rho_tvd, d_v_tvd, d_p_tvd
! TVD limiters for half steps for mass, velocity and pressure.

real(8), dimension(:,:), allocatable :: flux_u 
! Numerical fluxs for conservative variables.

real(8), dimension(:), allocatable :: flux_uu
! Numerical flux for entropy.   



end module solution