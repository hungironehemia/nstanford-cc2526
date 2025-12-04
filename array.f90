! Read from file

  program read
          implicit none

  OPEN (UNIT=11, FILE="verlet2.f90", STATUS="old", ACTION="read"

  READ (UNIT=11, FMT=*) a, b, c
  READ (UNIT=11, FMT=*) i, j, k  ![every READ corresponds to a line]
  READ (UNIT=11, FMT=*) x(1), x(2), x(3)


  end program

