MODULE neon_module

  USE kinds, ONLY: wp => dp

  IMPLICIT NONE


  PUBLIC :: C
  ! LJ PARAMETERS
  REAL(wp) :: sigma, epsilon

  CONTAINS

  ! Lennard-Jones potential
  FUNCTION v_lj(r) RESULT(v)
  REAL (KIND=wp), INTENT(IN) :: r
  REAL (KIND=wp) :: v
  
  v = 4.0_wp * epsilon * ((sigma/r)**12 - (sigma/r)**6)

  END FUNCTION v_lj

  ! Derivative of LJ potential
  FUNCTION dv_lj_dr(r) RESULT(dv)
  REAL (KIND=wp), INTENT(IN) :: r
  REAL (KIND=wp) :: dv

  dv = 4.0_wp * (-12.0_wp * sigma**12/r**13 + 6.0_wp * sigma**6/r**7)

  END FUNCTION dv_lj_dr

  ! Distance between two points
  FUNCTION distance(x1, y1, z1, x2, y2, z2) RESULT(r)
  REAL (KIND=wp), INTENT(IN) :: x1, y1, z1, x2, y2, z2
  REAL (KIND=wp) :: r
  REAL (KIND=wp) :: dx, dy, dz

  dx = x1 - x2
  dy = y1 - y2
  dz = z1 - z2

  r = sqrt(dx*dx + dy*dy + dz*dz)

  END FUNCTION distance


END MODULE neon_module
