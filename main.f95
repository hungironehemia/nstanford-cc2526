PROGRAM main
 USE kinds, ONLY: wp
 USE forces, ONLY: compute_forces
 IMPLICIT NONE

  ! Varibles from Input file
  INTEGER :: n, nk, k, i
  REAL (KIND=wp) :: tau, time, sigma, epsilon
  REAL (KIND=wp) :: mass_value

  ! Array for Atoms (Allocatable)
  REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: pos, vel, f
  REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: mass

  ! Unit Conversion (Bohr to Angiostrom)
  REAL (KIND=wp), PARAMETER :: bohr2ang = 0.529177_wp
  REAL (KIND=wp), PARAMETER :: amu_to_au = 1822.888485

  ! First Step: READ INPUT
  OPEN (UNIT=10, FILE='input_neon.txt', STATUS='OLD')

  READ (UNIT=10, FMT=*) nk, tau
  READ (UNIT=10, FMT=*) sigma, epsilon
  READ (UNIT=10, FMT=*) n

  PRINT *, "Reading input parameters"

  ! Memory Allocation for Arrays
  ALLOCATE(pos(n, 3))
  ALLOCATE(vel(n, 3))
  ALLOCATE(f(n, 3))
  ALLOCATE(mass(n))

  ! Read atom data (Mass, X, Y, Z, Vx, Vy, Vz)
  DO i = 1, n
  READ (UNIT=10, FMT=*) mass(i), pos(i,1), pos(i,2), pos(i,3), vel(i,1), vel(i,2), vel(i,3)

   ! Convert mass
   mass(i) = mass(i) * amu_to_au
   END DO

  CLOSE(10)

  ! Open output file 
  OPEN (UNIT=15, FILE='trajectories.xyz', STATUS='REPLACE')
  
  ! Second Step: Initialization
  time = 0.0_wp
  
  ! First force before moving
  CALL compute_forces(n, pos, f, sigma,epsilon)

  PRINT *, "Starting simulation for", nk, "steps"

  ! Third step: Velocity verlet (time loop)
  DO k = 1, nk

  ! Write current frame to XYZ file
  ! Format
  ! N
  ! Time
  ! Atom1 X Y Z
  WRITE (UNIT=15, FMT=*) n
  WRITE (UNIT=15, FMT=*) "Time = ", time
   DO i = 1, n
      ! Change values to Angstrom
      WRITE (UNIT=15, FMT='(A, 3F12.6)') "Ne", pos(i,1)*bohr2ang, pos(i,2)*bohr2ang, pos(i,3)*bohr2ang
    END DO

  ! Update velocity (Firts Half-step)
  ! v(t + 0.5*dt) = v(t) + 0.5 * (F/m) * dt
  ! Assuming all Neon atoms has same mass
  mass_value = mass(1)
  vel = vel + 0.5_wp * (f / mass_value) * tau

  ! Update position (Full step)
  ! r(t + dt) = r(t) + v(t + 0.5*dt) * dt
  pos = pos + vel * tau

  ! Update forces (New position -> New forces)
  CALL compute_forces(n, pos, f, sigma, epsilon)

  ! Update velocity (Second Half-Step)
  ! v(t + dt) = v(t + 0.5*dt) + 0.5 * (F_new/m) * dt
  vel = vel + 0.5_wp * (f / mass_value) * tau

  ! Advance Time
  time = time + tau

  END DO 

  PRINT *, "Output written to trajectories.xyz"
  CLOSE(15)

  ! Clean up memory
  DEALLOCATE(pos, vel, f, mass)

END PROGRAM main
