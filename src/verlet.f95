PROGRAM verlet
  IMPLICIT NONE
  
  ! Parameters: Mass (m), Time Step (tau), Total Steps (N)
  REAL, PARAMETER :: m = 1.0, tau = 0.2
  INTEGER, PARAMETER ::  n = 600

  ! Constant Force Components
  REAL, PARAMETER :: fx = 0.0, fy = 0.1, fz = 0.0

  ! State Variables (Position and Velocity)
  REAL :: x = 0.0, y = 0.0, z = 0.0
  REAL :: vx = 0.0, vy = 0.0, vz = 0.0

  OPEN (UNIT=11, FILE="trajectory.f90",STATUS="replace", ACTION="write")

  WRITE (UNIT=11, FMT=*) a, b, c
  WRITE (UNIT=11, FMT=*) x(1), (2), x(3)


  ! Loop counter and Time
  INTEGER :: k
  real :: time_k = 0.0

  PRINT *, 'Time (s), X (m), Y (m), Z (m)'
  PRINT *, time_k, x, y, z

  ! Main Simulation Loop
  DO k = 1, n
  time_k = time_k + tau

  ! New Position
  x = x + tau * vx + (tau**2 / (2.0 * m)) * fx
  y = y + tau * vy + (tau**2 / (2.0 * m)) * fy
  z = z + tau * vz + (tau**2 / (2.0 * m)) * fz

  ! New Velocities
  vx = vx + (tau / m) * fx
  vy = vy + (tau / m) * fy
  vz = vz + (tau / m) * fz

  PRINT *, time_k, x, y, z

  END DO
  
END PROGRAM verlet

