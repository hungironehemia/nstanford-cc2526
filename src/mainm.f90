PROGRAM main

 USE kinds, ONLY: wp
 USE neon_module, ONLY: sigma, epsilon, v_lj, dv_lj_dr, distance

 IMPLICIT NONE

 REAL (KIND=wp), PARAMETER :: amu_to_au = 1822.8888485_wp ! 1 amu in atomic mass units
 INTEGER :: nk          ! number of time steps
 REAL (KIND=wp) :: tau  ! time step
 INTEGER :: n           ! number of atoms
 REAL (KIND=wp), ALLOCATABLE :: mass(:) ! atomic mass (amu)
 REAL (KIND=wp), ALLOCATABLE :: x(:), y(:), z(:)      ! position
 REAL (KIND=wp), ALLOCATABLE :: vx(:), vy(:), vz(:)   ! velocities
 REAL (KIND=wp), ALLOCATABLE ::  fx(:), fy(:), fz(:)  ! forces
 INTEGER :: i, j, k, step
 REAL (KIND=wp) :: r, dx, dy, dz, force_mag, total_energy
 REAL (KIND=wp) :: kin_energy, pot_energy
 CHARACTER(LEN=100) :: input_file = "neon_input.txt"
 CHARACTER(LEN=100) :: output_file = "neon_output.xyz"

 ! Open input file [every READ correspond to a line]
 OPEN (UNIT=11, FILE="neon_input.txt", STATUS="old", ACTION="read")

 READ (UNIT=11, FMT=*) nk, tau
 READ (UNIT=11, FMT=*) sigma, epsilon
 READ (UNIT=11, FMT=*) n

 ! Arrays allocation
 allocate(mass(n), x(n), y(n), z(n), vx(n), vy(n), vz(n), fx(n), fy(n), fz(n))

 ! Initialize array to zero
 mass = 0.0_wp
 x = 0.0_wp; y = 0.0_wp; z = 0.0_wp
 vx = 0.0_wp; vy = 0.0_wp; vz = 0.0_wp
 fx = 0.0_wp fy = 0.0_wp; fz = 0.0_wp

 DO i = 1, n
 READ (UNIT=11, FMT=*) mass(i), x(i), y(i), z(i), vx(i), vy(i), vz(i)

 END DO

 CLOSE(11)


 ! Change masses into atomic unit from atomic mass unit
 mass = mass * amu_to_mu


 ! Output file [every WRITE correspond to a line]
 OPEN (UNIT=14, FILE="neon_output", STATUS="replace", ACTION="write")

 ! Calculate initial force
 CALL compute_force()

 ! Verlet velocity integration loop

 DO step = 1, nk

 ! First step: Update velocities (half step)
 DO i = 1, n

 vx(i) = vx(i) + 0.5_wp * tau * fx(i) / mass(i)
 vy(i) = vy(i) + 0.5_wp * tau * fy(i) / mass(i)
 vz(i) = vz(i) + 0.5_wp * tau * fz(i) / mass(i)

 END DO

 ! Second step: Update position (full step)
 x = x + tau * vx
 y = y + tau * vy
 z = z + tau * vz

 ! Third step: New forces
 CALL compute_forces()

 ! Fourth step: Update velocities second half step
 DO i = 1, n
 vx(i) = vx(i) + 0.5_wp * tau * fx(i) / mass(i)
 vy(i) = vy(i) + 0.5_wp * tau * fy(i) / mass(i)
 vz(i) = vz(i) + 0.5_wp * tau * fz(i) / mass(i)

 ! Write output file
 IF (mod(step, 10) == 0) then

 ! XYZ format
 WRITE (UNIT=14, FMT=*) n
 WRITE (UNIT=14, FMT=*) "Step", step

 DO i = 1, n

 ! Change Bohr to Angstrom for XYZ format (1 Bohr = 0.529177 Å)
  
 WRITE (UNIT=14, "(A2,3F12.6)") "Ne", x(i)*0.529177_wp, z(i)*0.5299177_wp

 END DO
 END IF

 CLOSE(14)

 PRINT *, "Well done", neon_output

 ! Subroutine to compute force
 SUBROUTINE compute_forces()
 
 ! Reset forces to zero
 fx = 0.0_wp
 fy = 0.0_wp
 fz = 0.0_wp

 pot_energy = 0.0_wp

 ! Loop over all unique pairs
 DO i = 1, n-1

 DO j = i+1, n

 ! Distance
 dx = x(i) - x(j)
 dy = y(i) - y(j)
 dz = z(i) - z(j)
 r = sqrt(dx*dx + dy*dy + dz*dz)

 IF (r > 0.0_wp) then ! Avoid division by zero
 ! Add to potential energy
 pot_energy = pot_energy + v_lj(r)

 ! Force magnitude
 force_mag = -dv_lj_dr(r) ! Negative of derivative gives force

 ! Force component (unit vector * force_mag)
 fx(i) = fx(i) + force_mag * dx/r
 fy(i) = fx(i) + force_mag * dy/r
 fz(i) = fz(i) + force_mag * dz/r

 ! Opposite force on j (Newton's 3rd law)
 fx(j) = fx(j) + force_mag * dx/r
 fy(j) = fx(j) + force_mag * dy/r
 fz(j) = fx(j) + force_mag * dz/r

 END IF

 END DO

 END DO

 END SUBROUTINE compute_forces


END PROGRAM main
