MODULE lj_module
  USE KINDS, ONLY: wp => dp
  IMPLICIT NONE

  PUBLIC :: Compute_force_lj, Compute_potential_energy

  CONTAINS

  ! Subroutine to compute LJ forces
  SUBROUTINE compute_force_lj(x, sigma, epsilon, f, epot)
     REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) :: x   ! positions (n,3)
     REAL (KIND=wp), INTENT(IN) :: sigma, epsilon       ! lj parameters
     REAL (KIND=wp), DIMENSION(:,:) INTENT(OUT) :: f    ! force (n,3)
     REAL (KIND=wp), INTENT(OUT) :: epot                 ! potential energy
     INTEGER :: n, i, j
     REAL (KIND=wp) :: rx, ry, rz, r
     REAL (KIND=wp) :: fx, fy, fz

     n = SIZE(x, 1)   ! particle's number
     f = 0.0_wp
     epot = 0.0_wp
  
   ! Loop over all pairs (i < j)
   DO i = 1, n - 1
   DO j = i+1, n

   ! Calculate components of distance
    rx = x(i,1) - x(j,1)
    ry = x(i,2) - x(j,2)
    rz = x(i,3) - x(j,3)

   ! Distance
   r = sqrt(rx*rx +ry*ry + rz*rz)

   ! Dont divide by zero 
   IF (r > 1.0e-10_wp) THEN

    ! Potential energy
    epot = epot + 4.0_wp * epsilon * (((sigma/r)**12) -((sigma/r)**6))


    ! Components of Force
    fx = -4.0_wp * epsilon * (1.0_wp / (r**2) * (-12.0_wp *((sigma / r)**))





  END MODULE lj_module

