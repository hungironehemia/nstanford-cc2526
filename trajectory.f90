PRINT *, a, b, c

  OPEN (UNIT="trajectory.f90", STATUS="old", ACTION="read")

  READ (UNIT=11, FMT=*) a, b, c
  READ (UNIT=11, FMT=*) i, j, k ![every READ correspond to a line]
  READ (UNIT=11, FMT=*) x(1), x(2), x(3)

  3.0 5.0 1.0E-2
  4 200 350
  0.0 0.0 0.0

CLOSE (11)

