MODULE trees

  ! Classes for building tree data structures

  USE morton_order

  IMPLICIT NONE

!******************************************************************************!

  TYPE node_ptr
    ! Pointer for derived type

     TYPE(node), POINTER :: p

  END TYPE node_ptr

!******************************************************************************!

  TYPE node
     INTEGER :: morton ! Morton index
     INTEGER :: depth  ! Depth relative to root
     INTEGER :: height ! Height of subtree
     INTEGER :: isleaf ! Denotes a leaf node
     INTEGER :: index  ! Index in final leaf list

     REAL(KIND = 8) :: value ! Average value in area represented by leaf

     TYPE(node_ptr), DIMENSION(0:3) :: children ! Children of node

  END TYPE node

!******************************************************************************!

CONTAINS

  !****************************************************************************!

  FUNCTION make_node(my_morton,my_depth,my_children) RESULT(my_node)
    ! Create branch-node (Acts as junction in tree)

    INTEGER, INTENT(IN) :: my_morton, my_depth
    TYPE(node), DIMENSION(0:3), INTENT(IN), TARGET :: my_children
    INTEGER :: i, max_height
    TYPE(node) :: my_node

    my_node % morton = my_morton
    my_node % depth = my_depth
    my_node % isleaf = 0

    max_height = my_children(0) % height

    DO i = 0,3

       my_node % children(i) % p => my_children(i)

       IF (my_children(i) % height > max_height) THEN
         max_height = my_children(i) % height
       ENDIF

    ENDDO

    my_node % height = max_height + 1

  END FUNCTION make_node

  !****************************************************************************!

  FUNCTION make_leaf(my_morton,my_depth,my_value) RESULT(my_leaf)
    ! Create leaf node (contains spatial information)

    INTEGER, INTENT(IN) :: my_morton, my_depth
    REAL(KIND = 8), INTENT(IN) :: my_value
    TYPE(node) :: my_leaf

    my_leaf % morton = my_morton
    my_leaf % depth = my_depth
    my_leaf % height = 0
    my_leaf % isleaf = 1
    my_leaf % value = my_value

  END FUNCTION make_leaf

  !****************************************************************************!

  FUNCTION find_leaf(root,l_morton,l_depth) RESULT(found)
    ! See if leaf is contained within tree

    TYPE(node), INTENT(IN) :: root
    INTEGER, INTENT(IN) :: l_morton, l_depth
    TYPE(node) :: my_node
    INTEGER :: ind, found

    my_node = root
    found = 0

    DO WHILE (my_node % height + my_node % depth >= l_depth .AND. &
      my_node % isleaf == 0)

      ind = m_above(l_morton,l_depth,my_node % depth + 1)
      ind = MOD(ind,4)
      my_node = my_node % children(ind) % p

    ENDDO

    IF (my_node % isleaf == 1) THEN
      found = 1
    ENDIF

  END FUNCTION find_leaf

  !****************************************************************************!

  SUBROUTINE list_leaves(root,num_leaves,num_nodes,leaf_morton,leaf_depth)
    ! Extract list of leaf nodes from quadtree

    TYPE(node), INTENT(IN) :: root
    INTEGER, INTENT(IN) :: num_leaves, num_nodes
    INTEGER, DIMENSION(num_leaves), INTENT(OUT) :: leaf_morton, leaf_depth
    TYPE(node) :: my_node
    TYPE(node), DIMENSION(num_nodes) :: node_list
    INTEGER :: i, ind_l, ind_n, add_n

    node_list(1) = root
    ind_l = 0
    ind_n = 1
    add_n = 1

    DO WHILE(ind_l < num_leaves)
      my_node = node_list(ind_n)
      ind_n = ind_n + 1

      DO i = 0,3
        IF (my_node % children(i) % p % isleaf == 1) THEN
          ind_l = ind_l + 1
          leaf_morton(ind_l) = my_node % children(i) % p % morton
          leaf_depth(ind_l) = my_node % children(i) % p % depth
        ELSE
          add_n = add_n + 1
          node_list(add_n) = my_node % children(i) % p
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE list_leaves

END MODULE trees
