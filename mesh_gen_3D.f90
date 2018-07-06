PROGRAM mesh_gen_3D

  USE array_io
  USE morton_order_3D ! Contains lookup tables ijk_to_m and m_to_ijk
  USE trees_3D

  IMPLICIT NONE

  ! I/O file variables
  CHARACTER(LEN = 125) :: input_file, output_file
  REAL(KIND = 8), ALLOCATABLE :: input_array(:,:,:), domain_array(:,:,:), mesh(:,:)
  REAL(KIND = 8) :: system_length, dx
  INTEGER :: input_M, input_N, input_L

  ! Morton indexing variables
  INTEGER(KIND = 8) :: m,next_m,m1,m2
  INTEGER :: x_min,x_max,y_min,y_max,z_min,z_max

  ! Tree Variables
  TYPE(node), ALLOCATABLE :: node_list(:)
  TYPE(node) :: next_push, octree
  INTEGER, DIMENSION(0:7) :: q
  INTEGER(KIND = 8), DIMENSION(0:151) :: neighbors
  INTEGER :: num_nodes, depth, d, all_leaf, balanced, n_num, n_depth
  INTEGER :: isthere, start_search
  INTEGER :: total_leaves, total_nodes
  INTEGER(KIND = 8), ALLOCATABLE :: leaf_morton(:)
  REAL(KIND = 8) :: val_avg, volume
  INTEGER, ALLOCATABLE :: ind_nodes(:), leaf_depth(:)

  ! Vertex variables
  INTEGER, ALLOCATABLE :: vert_ijk(:,:), vert_depth(:), vert_isthere(:)
  INTEGER, ALLOCATABLE :: N_neigh(:,:),S_neigh(:,:),E_neigh(:,:),W_neigh(:,:), &
       F_neigh(:,:), B_neigh(:,:)
  INTEGER :: num_vert, max_vert, i_0, j_0, k_0, offset, r, v_i, v_j, v_k, v_d, ind, special
  INTEGER, DIMENSION(8,3) :: corners

  ! Miscellaneous Variables
  INTEGER :: i,j,k,l

  !************************** Read in input file ******************************!

  ! Input file should be a square 2-D array of size (2^depth + 1)^2
  input_file = '/Users/kieranfitzmaurice/Documents/REU_2018/test_domain_3D.dat'
  output_file = '/Users/kieranfitzmaurice/Documents/REU_2018/refined_mesh_3D.dat'

  CALL size3D_binary(input_file,input_M,input_N,input_L)
  ALLOCATE(input_array(0:(input_M - 1),0:(input_N - 1),0:(input_L - 1)))
  input_array = input3D_binary(input_file,input_M,input_N,input_L)

  ! Side length of system
  system_length = 10.0

  ! Represent system by values at centroids
  ALLOCATE(domain_array(0:(input_M - 2), 0:(input_N - 2),0:(input_L - 2)))

  domain_array = (input_array(0:(input_M - 2),0:(input_N - 2),0:(input_L - 2)) &
       + input_array(1:(input_M - 1),0:(input_N - 2),0:(input_L - 2)) &
       + input_array(1:(input_M - 1),1:(input_N - 1),0:(input_L - 2)) &
       + input_array(0:(input_M - 2),1:(input_N - 1),0:(input_L - 2)) &
       + input_array(0:(input_M - 2),0:(input_N - 2),1:(input_L - 1)) &
       + input_array(1:(input_M - 1),0:(input_N - 2),1:(input_L - 1)) &
       + input_array(1:(input_M - 1),1:(input_N - 1),1:(input_L - 1)) &
       + input_array(0:(input_M - 2),1:(input_N - 1),1:(input_L - 1)))/8

  DEALLOCATE(input_array)

  !*************** Allocate memory based on size of system ********************!

  ! Level of refinement of input file
  depth = LOG(input_M - 1.0)/LOG(2.0)

  ! Initialize morton index lookup tables
  ! lookup tables are contained in module morton_order,
  ! and accessible to any subroutine or function using that module

  ALLOCATE(ijk_to_m(0:(2**depth - 1),0:(2**depth - 1),0:(2**depth - 1)))
  ALLOCATE(m_to_ijk(0:((2**depth)**3 - 1),0:2))
  CALL morton_table(depth,ijk_to_m,m_to_ijk)

  ALLOCATE(ind_nodes(0:depth))
  ! Total number of nodes we will need to allocate for
  num_nodes = 0

  DO i = depth,0,-1
     ind_nodes(i) = num_nodes
     num_nodes = num_nodes + (2**i)**3
  ENDDO

  ! List of all the nodes that will make up tree
  ALLOCATE(node_list(0:(num_nodes - 1)))

  ! Read in input elements as leaves at deepest level
  DO m = ind_nodes(depth),ind_nodes(depth - 1) - 1
     node_list(m) = make_leaf(m,depth, &
          domain_array(m_to_ijk(m,0),m_to_ijk(m,1),m_to_ijk(m,2)))
  ENDDO

  DEALLOCATE(domain_array)

  !*************** Coarsen areas where less resolution needed *****************!

  total_leaves = (2**depth)**3
  total_nodes = 0
  d = depth - 1

  DO WHILE (d >= 0)

     ! Read in current layer
     DO i = 0,(2**d)**3 - 1

        all_leaf = 1

        ! Read in indices of next 8 nodes of current layer
        DO j = 0,7
           q(j) = ind_nodes(d + 1) + 8*i + j
           IF (node_list(q(j)) % isleaf == 0) THEN
              all_leaf = 0
           ENDIF
        ENDDO

        next_m = i ! Morton index of next node to push back

        IF (all_leaf == 1) THEN ! If all leaves, see if can combine

           ! Average value in combined region
           val_avg = (node_list(q(0)) % value + node_list(q(1)) % value &
                + node_list(q(2)) % value + node_list(q(3)) % value &
                + node_list(q(4)) % value + node_list(q(5)) % value &
                + node_list(q(6)) % value + node_list(q(7)) % value)/8.0

           ! Volume of combined region
           volume = (system_length/2**d)**3

           ! Check if combining criteria are met
           IF (combine_criteria(d,volume,val_avg) == 1) THEN

              balanced = 1
              n_num = 0
              n_depth = d + 2 ! Can't have neighbors more than 1 level apart

              ! Don't check for neighbors deeper than deepest level
              IF (n_depth <= depth) THEN
                 ! Range of morton indices of points 2 levels below
                 ! that are "buried" by combined region
                 CALL m_below(next_m,d,n_depth,m1,m2)

                 x_min = m_to_ijk(m1,0)
                 y_min = m_to_ijk(m1,1)
                 z_min = m_to_ijk(m1,2)
                 x_max = m_to_ijk(m2,0)
                 y_max = m_to_ijk(m2,1)
                 z_max = m_to_ijk(m2,2)

                 !*** Build list of potential illegal neighbors morton indices ***

                 ! Cube faces
                 IF (x_min > 0) THEN
                    DO k = y_min,y_max
                       DO l = z_min,z_max
                          neighbors(n_num) = ijk_to_m(x_min - 1, k, l)
                          n_num = n_num + 1
                       ENDDO
                    ENDDO
                 ENDIF

                 IF (x_max < 2**n_depth - 1) THEN
                    DO k = y_min,y_max
                       DO l = z_min,z_max
                          neighbors(n_num) = ijk_to_m(x_max + 1, k, l)
                          n_num = n_num + 1
                       ENDDO
                    ENDDO
                 ENDIF

                 IF (y_min > 0) THEN
                    DO k = x_min,x_max
                       DO l = z_min,z_max
                          neighbors(n_num) = ijk_to_m(k, y_min - 1, l)
                          n_num = n_num + 1
                       ENDDO
                    ENDDO
                 ENDIF

                 IF (y_max < 2**n_depth - 1) THEN
                    DO k = x_min,x_max
                       DO l = z_min,z_max
                          neighbors(n_num) = ijk_to_m(k, y_max + 1, l)
                          n_num = n_num + 1
                       ENDDO
                    ENDDO
                 ENDIF

                 IF (z_min > 0) THEN
                    DO k = x_min,x_max
                       DO l = y_min,y_max
                          neighbors(n_num) = ijk_to_m(k, l, z_min - 1)
                          n_num = n_num + 1
                       ENDDO
                    ENDDO
                 ENDIF

                 IF (z_max < 2**n_depth - 1) THEN
                    DO k = x_min,x_max
                       DO l = y_min,y_max
                          neighbors(n_num) = ijk_to_m(k, l, z_max + 1)
                          n_num = n_num + 1
                       ENDDO
                    ENDDO
                 ENDIF

                 ! Cube edges
                 IF (x_min > 0) THEN

                    IF (y_min > 0) THEN
                       DO k = z_min,z_max
                          neighbors(n_num) = ijk_to_m(x_min - 1,y_min - 1,k)
                          n_num = n_num + 1
                       ENDDO
                    ENDIF

                    IF (z_min > 0) THEN
                       DO k = y_min,y_max
                          neighbors(n_num) = ijk_to_m(x_min - 1,k,z_min - 1)
                          n_num = n_num + 1
                       ENDDO
                    ENDIF

                    IF (y_max < 2**n_depth - 1) THEN
                       DO k = z_min,z_max
                          neighbors(n_num) = ijk_to_m(x_min - 1,y_max + 1,k)
                          n_num = n_num + 1
                       ENDDO
                    ENDIF

                    IF (z_max <= 2**n_depth - 2) THEN
                       DO k = y_min,y_max
                          neighbors(n_num) = ijk_to_m(x_min - 1,k,z_max + 1)
                          n_num = n_num + 1
                       ENDDO
                    ENDIF

                 ENDIF

                 IF (x_max < 2**n_depth - 1) THEN

                    IF (y_min > 0) THEN
                       DO k = z_min,z_max
                          neighbors(n_num) = ijk_to_m(x_max + 1,y_min - 1,k)
                          n_num = n_num + 1
                       ENDDO
                    ENDIF

                    IF (z_min > 0) THEN
                       DO k = y_min,y_max
                          neighbors(n_num) = ijk_to_m(x_max + 1,k,z_min - 1)
                          n_num = n_num + 1
                       ENDDO
                    ENDIF

                    IF (y_max < 2**n_depth - 1) THEN
                       DO k = z_min,z_max
                          neighbors(n_num) = ijk_to_m(x_max + 1,y_max + 1,k)
                          n_num = n_num + 1
                       ENDDO
                    ENDIF

                    IF (z_max < 2**n_depth - 1) THEN
                       DO k = y_min,y_max
                          neighbors(n_num) = ijk_to_m(x_max + 1,k,z_max + 1)
                          n_num = n_num + 1
                       ENDDO
                    ENDIF

                 ENDIF

                 IF (y_min > 0) THEN

                    IF (z_min > 0) THEN
                       DO k = x_min,x_max
                          neighbors(n_num) = ijk_to_m(k,y_min - 1,z_min - 1)
                          n_num = n_num + 1
                       ENDDO
                    ENDIF

                    IF (z_max < 2**n_depth - 1) THEN
                       DO k = x_min,x_max
                          neighbors(n_num) = ijk_to_m(k,y_min - 1,z_max + 1)
                          n_num = n_num + 1
                       ENDDO
                    ENDIF

                 ENDIF

                 IF (y_max < 2**n_depth - 1) THEN

                    IF (z_min > 0) THEN
                       DO k = x_min,x_max
                          neighbors(n_num) = ijk_to_m(k,y_max + 1,z_min - 1)
                          n_num = n_num + 1
                       ENDDO
                    ENDIF

                    IF (z_max < 2**n_depth - 1) THEN
                       DO k = x_min,x_max
                          neighbors(n_num) = ijk_to_m(k,y_max + 1,z_max + 1)
                          n_num = n_num + 1
                       ENDDO
                    ENDIF

                 ENDIF

                 ! Cube corners
                 IF (x_min > 0) THEN

                    IF (y_min > 0) THEN

                       IF (z_min > 0) THEN
                          neighbors(n_num) = ijk_to_m(x_min - 1, y_min - 1, z_min - 1)
                          n_num = n_num + 1
                       ENDIF
                       IF (z_max < 2**n_depth - 1) THEN
                          neighbors(n_num) = ijk_to_m(x_min - 1, y_min - 1, z_max + 1)
                          n_num = n_num + 1
                       ENDIF

                    ENDIF

                    IF (y_max < 2**n_depth - 1) THEN
                       IF (z_min > 0) THEN
                          neighbors(n_num) = ijk_to_m(x_min - 1, y_max + 1, z_min - 1)
                          n_num = n_num + 1
                       ENDIF
                       IF (z_max < 2**n_depth - 1) THEN
                          neighbors(n_num) = ijk_to_m(x_min - 1, y_max + 1, z_max + 1)
                          n_num = n_num + 1
                       ENDIF
                    ENDIF

                 ENDIF

                 IF (x_max < 2**n_depth - 1) THEN

                    IF (y_min > 0) THEN

                       IF (z_min > 0) THEN
                          neighbors(n_num) = ijk_to_m(x_max + 1, y_min - 1, z_min - 1)
                          n_num = n_num + 1
                       ENDIF
                       IF (z_max < 2**n_depth - 1) THEN
                          neighbors(n_num) = ijk_to_m(x_max + 1, y_min - 1, z_max + 1)
                          n_num = n_num + 1
                       ENDIF

                    ENDIF

                    IF (y_max < 2**n_depth - 1) THEN
                       IF (z_min > 0) THEN
                          neighbors(n_num) = ijk_to_m(x_max + 1, y_max + 1, z_min - 1)
                          n_num = n_num + 1
                       ENDIF
                       IF (z_max < 2**n_depth - 1) THEN
                          neighbors(n_num) = ijk_to_m(x_max + 1, y_max + 1, z_max + 1)
                          n_num = n_num + 1
                       ENDIF
                    ENDIF

                 ENDIF

                 !*** Check if illegal neighbors are present ***
                 DO k = 0,n_num - 1
                    isthere = 0
                    start_search = FLOOR(neighbors(k)/8.0) ! Index to start search

                    IF (node_list(ind_nodes(d + 1) + start_search) % isleaf == 0) THEN
                       isthere = find_leaf(node_list(ind_nodes(d + 1) + start_search), &
                            neighbors(k),n_depth)
                    ENDIF

                    ! If creates illegal neighbors, reject attempt to combine
                    IF (isthere == 1) THEN
                       balanced = 0
                       EXIT
                    ENDIF

                 ENDDO

              ENDIF

              IF (balanced == 1) THEN ! If mesh stays balanced, combine leaves
                 next_push = make_leaf(next_m,d,val_avg)
                 total_leaves = total_leaves - 7
              ELSE ! If mesh won't be balanced, create node instead
                 next_push = make_node(next_m,d,node_list(q(0):q(7)))
                 total_nodes = total_nodes + 1
              ENDIF

           ELSE ! If combining criteria not met, make node
              next_push = make_node(next_m,d,node_list(q(0):q(7)))
              total_nodes = total_nodes + 1
           ENDIF

        ELSE  ! If not all 8 are leaves, then make node
           next_push = make_node(next_m,d,node_list(q(0):q(7)))
           total_nodes = total_nodes + 1
        ENDIF

        ! Push to next level
        node_list(ind_nodes(d) + i) = next_push

     ENDDO

     ! Move on to next level of refinement
     d = d - 1

  ENDDO

  octree = node_list(ind_nodes(0))

  !******************* Extract List of leaves from tree ***********************!

  ALLOCATE(leaf_depth(total_leaves))
  ALLOCATE(leaf_morton(total_leaves))

  CALL list_leaves(octree,total_leaves,total_nodes,leaf_morton,leaf_depth)

  ! Free up memory no longer needed
  DEALLOCATE(node_list)

  !*************** Move from cell centers to cell vertices ********************!

  max_vert = 8 + 19*(total_leaves - 1)/7 ! Maximum # of vertices in worst case
  ALLOCATE(vert_isthere(0:(input_M**3 - 1))) ! Records if unique vertex is present or not
  ALLOCATE(vert_ijk(max_vert,3))  ! List of unique vertices
  ALLOCATE(vert_depth(max_vert)) ! Depth of leaf that vertex is attached to

  num_vert = 0
  vert_isthere = 0

  DO i = 1,total_leaves
     CALL m_below(leaf_morton(i),leaf_depth(i),depth,m1,m2)

     ! Difference in i or j index between different corners
     offset = 2**(depth - leaf_depth(i))

     i_0 = m_to_ijk(m1,0)
     j_0 = m_to_ijk(m1,1)
     k_0 = m_to_ijk(m1,2)

     ! Vertices of each leaf
     corners(1,1) = i_0
     corners(1,2) = j_0
     corners(1,3) = k_0

     corners(2,1) = i_0
     corners(2,2) = j_0 + offset
     corners(2,3) = k_0

     corners(3,1) = i_0
     corners(3,2) = j_0
     corners(3,3) = k_0 + offset

     corners(4,1) = i_0
     corners(4,2) = j_0 + offset
     corners(4,3) = k_0 + offset

     corners(5,1) = i_0 + offset
     corners(5,2) = j_0
     corners(5,3) = k_0

     corners(6,1) = i_0 + offset
     corners(6,2) = j_0 + offset
     corners(6,3) = k_0

     corners(7,1) = i_0 + offset
     corners(7,2) = j_0
     corners(7,3) = k_0 + offset

     corners(8,1) = i_0 + offset
     corners(8,2) = j_0 + offset
     corners(8,3) = k_0 + offset

     DO j = 1,8
        r = ijk_to_r(corners(j,1),corners(j,2),corners(j,3),input_M)
        IF (vert_isthere(r) == 0) THEN
           num_vert = num_vert + 1
           vert_ijk(num_vert,1) = corners(j,1)
           vert_ijk(num_vert,2) = corners(j,2)
           vert_ijk(num_vert,3) = corners(j,3)
           vert_depth(num_vert) = leaf_depth(i)
           vert_isthere(r) = num_vert
        ENDIF
     ENDDO

  ENDDO

  !********************* Get connectivity of nodes ****************************!

  ! List of neighbor(s) in given direction
  ALLOCATE(N_neigh(num_vert,0:4))
  ALLOCATE(S_neigh(num_vert,0:4))
  ALLOCATE(E_neigh(num_vert,0:4))
  ALLOCATE(W_neigh(num_vert,0:4))
  ALLOCATE(F_neigh(num_vert,0:4))
  ALLOCATE(B_neigh(num_vert,0:4))

  ! Each vertex starts out without neighbors
  N_neigh = 0
  S_neigh = 0
  E_neigh = 0
  W_neigh = 0
  F_neigh = 0
  B_neigh = 0

  DO i = 1,num_vert

     v_d = vert_depth(i)
     v_i = vert_ijk(i,1)
     v_j = vert_ijk(i,2)
     v_k = vert_ijk(i,3)

     !********************** Check for north neighbor **************************!
     special = 1
     DO d = v_d + 1, v_d - 1, -1

        IF (d > depth) THEN
           CYCLE
        ENDIF

        offset = 2**(depth - d)

        IF (v_i - offset >= 0) THEN
           r = ijk_to_r(v_i - offset, v_j, v_k, input_M)
           ind = vert_isthere(r)

           IF (ind > 0) THEN
              N_neigh(i,1) = ind
              N_neigh(i,0) = 1
              special = 0
              EXIT
           ENDIF

        ELSE
           special = 0
           EXIT
        ENDIF
     ENDDO

     IF (special == 1) THEN ! Special case where neighbors must be interpolated for
        ! 1st neighbor
        r = ijk_to_r(v_i - offset, v_j - offset/2, v_k - offset/2, input_M)
        ind = vert_isthere(r)
        N_neigh(i,1) = ind

        ! 2nd neighbor
        r = ijk_to_r(v_i - offset, v_j - offset/2, v_k + offset/2, input_M)
        ind = vert_isthere(r)
        N_neigh(i,2) = ind

        ! 3rd neighbor
        r = ijk_to_r(v_i - offset, v_j + offset/2, v_k - offset/2, input_M)
        ind = vert_isthere(r)
        N_neigh(i,3) = ind

        ! 4th neighbor
        r = ijk_to_r(v_i - offset, v_j + offset/2, v_k + offset/2, input_M)
        ind = vert_isthere(r)
        N_neigh(i,4) = ind

        ! Will have four "north" neighbors
        N_neigh(i,0) = 4
     ENDIF


     !********************** Check for south neighbor **************************!
     special = 1
     DO d = v_d + 1, v_d - 1, -1

        IF (d > depth) THEN
           CYCLE
        ENDIF

        offset = 2**(depth - d)

        IF (v_i + offset < input_M) THEN
           r = ijk_to_r(v_i + offset, v_j, v_k, input_M)
           ind = vert_isthere(r)

           IF (ind > 0) THEN
              S_neigh(i,1) = ind
              S_neigh(i,0) = 1
              special = 0
              EXIT
           ENDIF

        ELSE
           special = 0
           EXIT
        ENDIF
     ENDDO

     IF (special == 1) THEN ! Special case where neighbors must be interpolated for
        ! 1st neighbor
        r = ijk_to_r(v_i + offset, v_j - offset/2, v_k - offset/2, input_M)
        ind = vert_isthere(r)
        S_neigh(i,1) = ind

        ! 2nd neighbor
        r = ijk_to_r(v_i + offset, v_j - offset/2, v_k + offset/2, input_M)
        ind = vert_isthere(r)
        S_neigh(i,2) = ind

        ! 3rd neighbor
        r = ijk_to_r(v_i + offset, v_j + offset/2, v_k - offset/2, input_M)
        ind = vert_isthere(r)
        S_neigh(i,3) = ind

        ! 4th neighbor
        r = ijk_to_r(v_i + offset, v_j + offset/2, v_k + offset/2, input_M)
        ind = vert_isthere(r)
        S_neigh(i,4) = ind

        ! Will have four "south" neighbors
        S_neigh(i,0) = 4
     ENDIF

     !*********************** Check for east neighbor **************************!
     special = 1
     DO d = v_d + 1, v_d - 1, -1

        IF (d > depth) THEN
           CYCLE
        ENDIF

        offset = 2**(depth - d)

        IF (v_j + offset < input_M) THEN
           r = ijk_to_r(v_i, v_j + offset, v_k, input_M)
           ind = vert_isthere(r)

           IF (ind > 0) THEN
              E_neigh(i,1) = ind
              E_neigh(i,0) = 1
              special = 0
              EXIT
           ENDIF

        ELSE
           special = 0
           EXIT
        ENDIF
     ENDDO

     IF (special == 1) THEN ! Special case where neighbors must be interpolated for
        ! 1st neighbor
        r = ijk_to_r(v_i - offset/2, v_j + offset, v_k - offset/2, input_M)
        ind = vert_isthere(r)
        E_neigh(i,1) = ind

        ! 2nd neighbor
        r = ijk_to_r(v_i - offset/2, v_j + offset, v_k + offset/2, input_M)
        ind = vert_isthere(r)
        E_neigh(i,2) = ind

        ! 3rd neighbor
        r = ijk_to_r(v_i + offset/2, v_j + offset, v_k - offset/2, input_M)
        ind = vert_isthere(r)
        E_neigh(i,3) = ind

        ! 4th neighbor
        r = ijk_to_r(v_i + offset/2, v_j + offset, v_k + offset/2, input_M)
        ind = vert_isthere(r)
        E_neigh(i,4) = ind

        ! Will have four "east" neighbors
        E_neigh(i,0) = 4
     ENDIF

     !*********************** Check for west neighbor **************************!
     special = 1
     DO d = v_d + 1, v_d - 1, -1

        IF (d > depth) THEN
           CYCLE
        ENDIF

        offset = 2**(depth - d)

        IF (v_j - offset >= 0) THEN
           r = ijk_to_r(v_i, v_j - offset, v_k, input_M)
           ind = vert_isthere(r)

           IF (ind > 0) THEN
              W_neigh(i,1) = ind
              W_neigh(i,0) = 1
              special = 0
              EXIT
           ENDIF

        ELSE
           special = 0
           EXIT
        ENDIF
     ENDDO

     IF (special == 1) THEN ! Special case where neighbors must be interpolated for
        ! 1st neighbor
        r = ijk_to_r(v_i - offset/2, v_j - offset, v_k - offset/2, input_M)
        ind = vert_isthere(r)
        W_neigh(i,1) = ind

        ! 2nd neighbor
        r = ijk_to_r(v_i - offset/2, v_j - offset, v_k + offset/2, input_M)
        ind = vert_isthere(r)
        W_neigh(i,2) = ind

        ! 3rd neighbor
        r = ijk_to_r(v_i + offset/2, v_j - offset, v_k - offset/2, input_M)
        ind = vert_isthere(r)
        W_neigh(i,3) = ind

        ! 4th neighbor
        r = ijk_to_r(v_i + offset/2, v_j - offset, v_k + offset/2, input_M)
        ind = vert_isthere(r)
        W_neigh(i,4) = ind

        ! Will have four "west" neighbors
        W_neigh(i,0) = 4
     ENDIF

     !********************** Check for forward neighbor ************************!
     special = 1
     DO d = v_d + 1, v_d - 1, -1

        IF (d > depth) THEN
           CYCLE
        ENDIF

        offset = 2**(depth - d)

        IF (v_k - offset >= 0) THEN
           r = ijk_to_r(v_i, v_j, v_k - offset, input_M)
           ind = vert_isthere(r)

           IF (ind > 0) THEN
              F_neigh(i,1) = ind
              F_neigh(i,0) = 1
              special = 0
              EXIT
           ENDIF

        ELSE
           special = 0
           EXIT
        ENDIF
     ENDDO

     IF (special == 1) THEN ! Special case where neighbors must be interpolated for
        ! 1st neighbor
        r = ijk_to_r(v_i - offset/2, v_j - offset/2, v_k - offset, input_M)
        ind = vert_isthere(r)
        F_neigh(i,1) = ind

        ! 2nd neighbor
        r = ijk_to_r(v_i - offset/2, v_j + offset/2, v_k - offset, input_M)
        ind = vert_isthere(r)
        F_neigh(i,2) = ind

        ! 3rd neighbor
        r = ijk_to_r(v_i + offset/2, v_j - offset/2, v_k - offset, input_M)
        ind = vert_isthere(r)
        F_neigh(i,3) = ind

        ! 4th neighbor
        r = ijk_to_r(v_i + offset/2, v_j + offset/2, v_k - offset, input_M)
        ind = vert_isthere(r)
        F_neigh(i,4) = ind

        ! Will have four "forward" neighbors
        F_neigh(i,0) = 4
     ENDIF

     !********************* Check for backward neighbor ***********************!
     special = 1
     DO d = v_d + 1, v_d - 1, -1

        IF (d > depth) THEN
           CYCLE
        ENDIF

        offset = 2**(depth - d)

        IF (v_k + offset < input_M) THEN
           r = ijk_to_r(v_i, v_j, v_k + offset, input_M)
           ind = vert_isthere(r)

           IF (ind > 0) THEN
              B_neigh(i,1) = ind
              B_neigh(i,0) = 1
              special = 0
              EXIT
           ENDIF

        ELSE
           special = 0
           EXIT
        ENDIF
     ENDDO

     IF (special == 1) THEN ! Special case where neighbors must be interpolated for
        ! 1st neighbor
        r = ijk_to_r(v_i - offset/2, v_j - offset/2, v_k + offset, input_M)
        ind = vert_isthere(r)
        B_neigh(i,1) = ind

        ! 2nd neighbor
        r = ijk_to_r(v_i - offset/2, v_j + offset/2, v_k + offset, input_M)
        ind = vert_isthere(r)
        B_neigh(i,2) = ind

        ! 3rd neighbor
        r = ijk_to_r(v_i + offset/2, v_j - offset/2, v_k + offset, input_M)
        ind = vert_isthere(r)
        B_neigh(i,3) = ind

        ! 4th neighbor
        r = ijk_to_r(v_i + offset/2, v_j + offset/2, v_k + offset, input_M)
        ind = vert_isthere(r)
        B_neigh(i,4) = ind

        ! Will have four "backward" neighbors
        B_neigh(i,0) = 4

     ENDIF

  ENDDO

  !************************* Write mesh to output *****************************!

  DEALLOCATE(vert_isthere)
  ALLOCATE(mesh(num_vert,34))

  ! Read back in input array
  ALLOCATE(input_array(0:(input_M - 1),0:(input_N - 1),0:(input_L - 1)))
  input_array = input3D_binary(input_file,input_M,input_N,input_L)

  dx = system_length / 2**depth

  DO i = 1,num_vert
    mesh(i,1) = input_array(vert_ijk(i,1),vert_ijk(i,2),vert_ijk(i,3))
    mesh(i,2) = vert_ijk(i,1)*dx
    mesh(i,3) = vert_ijk(i,2)*dx
    mesh(i,4) = vert_ijk(i,3)*dx
    mesh(i,5:9) = N_neigh(i,0:4)
    mesh(i,10:14) = S_neigh(i,0:4)
    mesh(i,15:19) = E_neigh(i,0:4)
    mesh(i,20:24) = W_neigh(i,0:4)
    mesh(i,25:29) = F_neigh(i,0:4)
    mesh(i,30:34) = B_neigh(i,0:4)
  ENDDO

  CALL output2D_binary(output_file,mesh,num_vert,34)

  !************************* End of main program ******************************!


CONTAINS

  FUNCTION combine_criteria(d,volume,val) RESULT(c)

    INTEGER, INTENT(IN) :: d
    REAL(KIND = 8), INTENT(IN) :: volume, val
    INTEGER :: min_depth, c

    min_depth = 0

    IF (d >= min_depth) THEN
       IF (val*volume > 0.001) THEN
          c = 0 ! Do not combine
       ELSE
          c = 1 ! Combine
       ENDIF

    ELSE ! Do not combine
       c = 0
    ENDIF

  END FUNCTION combine_criteria

  FUNCTION ijk_to_r(i,j,k,L) RESULT(r)
    ! Converts i,j indices into index r in raster scan of space
    INTEGER, INTENT(IN) :: i,j,k,L
    INTEGER :: r

    r = i*L**2 + j*L + k

  END FUNCTION ijk_to_r

  SUBROUTINE r_to_ijk(r,L,i,j,k)
    INTEGER, INTENT(IN) :: r, L
    INTEGER, INTENT(OUT) :: i,j,k

    i = FLOOR(REAL(r)/L**2)
    j = FLOOR(REAL(r - i*L**2)/L)
    k = MOD(r,L)

  END SUBROUTINE r_to_ijk

END PROGRAM mesh_gen_3D
