MODULE morton_order_3D
  IMPLICIT NONE

  ! Lookup tables
  INTEGER(KIND = 8), ALLOCATABLE :: ijk_to_m(:,:,:)
  INTEGER(KIND = 8), ALLOCATABLE :: m_to_ijk(:,:)

CONTAINS

  !****************************************************************************!

  FUNCTION to_morton(i,j,k) RESULT(m)
    ! Converts 3-D indicies (i,j,k) into corresponding morton index (m)

    INTEGER(KIND = 8), INTENT(IN) :: i,j,k
    INTEGER(KIND = 8) :: m
    INTEGER :: n

    m = 0

    DO n = 0,15
       CALL MVBITS(i,n,1,m,3*n)
       CALL MVBITS(j,n,1,m,3*n + 1)
       CALL MVBITS(k,n,1,m,3*n + 2)
    ENDDO

  END FUNCTION to_morton

  !****************************************************************************!

  SUBROUTINE from_morton(m,i,j,k)
    ! Converts morton index (m) into correspoding 3-D indices (i,j,k)

    INTEGER(KIND = 8), INTENT(IN) :: m
    INTEGER(KIND = 8), INTENT(OUT) :: i,j,k
    INTEGER :: n

    i = 0
    j = 0
    k = 0

    DO n = 0,15
       CALL MVBITS(m,3*n,1,i,n)
       CALL MVBITS(m,3*n + 1,1,j,n)
       CALL MVBITS(m,3*n + 2,1,k,n)
    ENDDO

  END SUBROUTINE from_morton

  !****************************************************************************!

  SUBROUTINE morton_table(depth,ijk_to_m,m_to_ijk)
    ! Initialize morton lookup tables to convert between i,j,k and morton indices

    INTEGER, INTENT(IN) :: depth
    INTEGER(KIND = 8), DIMENSION(0:(2**depth - 1),0:(2**depth - 1),0:(2**depth - 1)), &
         INTENT(OUT) :: ijk_to_m
    INTEGER(KIND = 8), DIMENSION(0:((2**depth)**3 - 1),0:2), INTENT(OUT) :: m_to_ijk
    INTEGER(KIND = 8) :: i, j, k, m

    DO k = 0,2**depth - 1
       DO j = 0,2**depth - 1
          DO i = 0,2**depth - 1

             m = to_morton(i,j,k)

             ijk_to_m(i,j,k) = m
             m_to_ijk(m,0) = i
             m_to_ijk(m,1) = j
             m_to_ijk(m,2) = k

          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE morton_table

  !****************************************************************************!

  FUNCTION m_above(m,d1,d2) RESULT(m_d2)
    ! Returns morton index that index in level d1 would have in level d2
    !
    ! m = morton index at d1
    ! d1 = current depth
    ! d2 = depth shallower than d1

    INTEGER(KIND = 8), INTENT(IN) :: m
    INTEGER, INTENT(IN) :: d1,d2
    INTEGER(KIND = 8) :: m_d2

    m_d2 = FLOOR(REAL(m) / 8**(d1 - d2))

  END FUNCTION m_above

  !****************************************************************************!

  SUBROUTINE m_below(m,d1,d2,m_min,m_max)
    ! Returns range of morton indicies in level d2 that fall under index in d1
    !
    ! m = morton index at d1
    ! d1 = current depth
    ! d2 = depth deeper than d2

    INTEGER(KIND = 8), INTENT(IN) :: m
    INTEGER, INTENT(IN):: d1,d2
    INTEGER(KIND = 8), INTENT(OUT) :: m_min, m_max

    m_min = m*8**(d2 - d1)
    m_max = (m + 1)*8**(d2 - d1) - 1

  END SUBROUTINE m_below

  !****************************************************************************!

END MODULE morton_order_3D
