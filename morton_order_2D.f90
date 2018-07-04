MODULE morton_order_2D
  IMPLICIT NONE
  
  ! Lookup tables
  INTEGER, ALLOCATABLE :: ij_to_m(:,:)
  INTEGER, ALLOCATABLE :: m_to_ij(:,:)

CONTAINS
  !****************************************************************************!

  FUNCTION to_morton(i,j) RESULT(m)
    ! Converts 2-D indicies (i,j) into corresponding morton index (m)
    ! Using 4 byte integers ...
    ! Max value of i, j, k coordinates is 2^15
    ! Max value of morton coordinate using is 2^31
    ! Can use 8 byte if higher values needed

    INTEGER(KIND = 4), INTENT(IN) :: i,j
    INTEGER(KIND = 4) :: m
    INTEGER :: k

    m = 0

    DO k = 0,15

      CALL MVBITS(i,k,1,m,2*k)
      CALL MVBITS(j,k,1,m,2*k + 1)

    ENDDO

  END FUNCTION to_morton

  !****************************************************************************!

  SUBROUTINE from_morton(m,i,j)
    ! Converts morton index (m) into corresponding 2-D indicies (i,j)

    INTEGER(KIND = 4), INTENT(IN) :: m
    INTEGER(KIND = 4), INTENT(OUT) :: i,j
    INTEGER :: k

    i = 0
    j = 0

    DO k = 0,15

      CALL MVBITS(m,2*k,1,i,k)
      CALL MVBITS(m,2*k + 1,1,j,k)

    ENDDO

  END SUBROUTINE from_morton

  !****************************************************************************!

  SUBROUTINE morton_table(depth,ij_to_m,m_to_ij)
    ! Generates morton lookup tables to convert between i,j and morton indices

    INTEGER, INTENT(IN) :: depth
    INTEGER, DIMENSION(0:(2**depth - 1),0:(2**depth - 1)), INTENT(OUT) :: ij_to_m
    INTEGER, DIMENSION(0:((2**depth)**2 - 1),0:1), INTENT(OUT) :: m_to_ij
    INTEGER(KIND = 4) :: i, j, m

    DO j = 0,2**depth - 1
      DO i = 0,2**depth - 1

        m = to_morton(i,j)

        ij_to_m(i,j) = m
        m_to_ij(m,0) = i
        m_to_ij(m,1) = j

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

    INTEGER, INTENT(IN) :: m, d1, d2
    INTEGER :: m_d2

    m_d2 = FLOOR(REAL(m) / 4**(d1 - d2))

  END FUNCTION m_above

  !****************************************************************************!

  SUBROUTINE m_below(m,d1,d2,m_min,m_max)
    ! Returns range of morton indicies in level d2 that fall under index in d1
    !
    ! m = morton index at d1
    ! d1 = current depth
    ! d2 = depth deeper than d2

    INTEGER, INTENT(IN) :: m, d1, d2
    INTEGER, INTENT(OUT) :: m_min, m_max

    m_min = m*4**(d2 - d1)
    m_max = (m + 1)*4**(d2 - d1) - 1

  END SUBROUTINE m_below

END MODULE morton_order_2D
