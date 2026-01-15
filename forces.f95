MODULE forces
 
        USE kinds, ONLY: wp

        IMPLICIT NONE

CONTAINS

        SUBROUTINE compute_forces(n, pos, f, sigma, epsilon)
  ! ARGUMENTS
  INTEGER, INTENT(IN) :: n                            ! Number of atoms
  REAL (KIND=wp), DIMENSION(:,:), INTENT(IN) :: pos   ! Position
  REAL (KIND=wp), DIMENSION(:,:), INTENT(OUT) :: f    ! Force (Result)
  REAL (KIND=wp), INTENT(IN) :: sigma, epsilon       ! Constants


  ! VARIABLES
  INTEGER :: i, j
  REAL (KIND=wp) :: dx, dy, dz, distance, distance_sq
  REAL (KIND=wp) :: repulsion, attraction, force_mag

  ! Before starting
  f = 0.0_wp

  ! Loop every every atom through (Atom i)
  DO i = 1, n
  ! Loop through every other atom (Atom j)
   DO j = i + 1, n

   ! Step 1: Compute distance
   dx = pos(i, 1) - pos(j, 1)
   dy = pos(i, 2) - pos(j, 2)
   dz = pos(i, 3) - pos(j, 3)
  
   distance_sq = dx**2 + dy**2 + dz**2
   distance = sqrt(distance_sq)

   ! Step 2: Calculate push/pull strength

   ! The '12' part is the Repulsion
   ! From the formula: -12 * (sigma^12 / r^13)
   repulsion = -12.0_wp * (sigma**12) / (distance**13)

   ! The '6' part is the Attraction
   ! From the formula: 6 * (sigma^6 / r^7)
   attraction = 6.0_wp * (sigma**6) / (distance**7)

   ! On combining them together
   force_mag = 4.0_wp * epsilon * (repulsion + attraction)

   ! Step 4: Apply the Forces in X, Y, Z direction
   ! Multiply the force_mag (V'lj) by the direction (dx/distance)

   ! Updating Atom i ('The Action')
   f(i, 1) = f(i, 1) - (dx / distance) * force_mag
   f(i, 2) = f(i, 2) - (dy / distance) * force_mag
   f(i, 3) = f(i, 3) - (dz / distance) * force_mag

   ! Updating Atom j (The 'Reaction' - Equal and Opposite)
    f(i, 1) = f(j, 1) - (dx / distance) * force_mag
    f(j, 2) = f(j, 2) - (dy / distance) * force_mag
    f(j, 3) = f(j, 3) - (dz / distance) * force_mag
    END DO
   END DO

   END SUBROUTINE compute_forces

END MODULE forces
