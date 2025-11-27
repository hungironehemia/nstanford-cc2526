PROGRAM verlet
  IMPLICIT NONE
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND (p=13, r=300)
  INTEGER :: i, j, k
  INTEGER :: l, m
  CHARACTER (LEN=72) :: str1, str2, str3
  REAL (KIND=wp), DIMENSION(:), ALLOCATABLE :: x
  REAL (KIND=wp), DIMENSION(:,:), ALLOCATABLE :: xmat
  
  ALLOCATE ( x (300) )
  ALLOCATE ( xmat (200, 300) )
  x(3)=4.0_wp
  PRINT *, x(3)
  
END PROGRAM verlet


