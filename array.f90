 ! read

   program read

    OPEN (UNIT=12, FILE="verlet2.f90", STATUS="replace", ACTION="write")

  WRITE (UNIT=12, FMT=*) a, b, c
  WRITE (UNIT=12, FMT=*) x(1), x(2), x(3)


  end program
