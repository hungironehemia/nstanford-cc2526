PROGRAM verlet
  IMPLICIT NONE
  
  ! Parameters: Mass (m), Time Step (tau), Total Steps (N)
  REAL :: m, tau
  INTEGER ::  n

  ! Constant Force Components
  REAL :: fx, fy, fz

  ! State Variables (Position and Velocity)
  REAL :: x, y, z
  REAL :: vx, vy, vz
  REAL :: time_k

  ! Loop counter and Time
  INTEGER :: k

  ! Read input file
  OPEN (UNIT=13, FILE="verlet_input.txt", STATUS="old", ACTION="read")

  ! Each READ correspond to a line
  READ (UNIT=13, FMT=*) m, tau, n
  READ (UNIT=13, FMT=*) fx, fy, fz
  READ (UNIT=13, FMT=*) x, y, z
  READ (UNIT=13, FMT=*) vx, vy, vz
  READ (UNIT=13, FMT=*) time_k

  CLOSE(13)

  ! Write output file
  OPEN (UNIT=14, FILE="verlet_output.text", STATUS="replace", ACTION="write")

  ! Each WRITE corraspond to a line
  WRITE(UNIT=14, FMT=*) "Time (s), x (m), y (m), z (m)"
  WRITE(UNIT=14, FMT=*) time_k, x, y, z

 
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


  WRITE (UNIT=14, FMT=*) time_k, x, y, z

  END DO

  CLOSE(14)

  PRINT *, "Completed well"
  
END PROGRAM verlet




