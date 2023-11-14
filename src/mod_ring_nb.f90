!Brief: Bookkeeping for ring-ring neighborhood
!Principle:
!! ring_nb_pairs is a flat(1D) upper triangular matrix to store neighborness of all possible ring-ring pairs.
!! pair_to_index() maps any given ring-ring pair to a unique index of ring_nb_pairs.
!! index_to_pair() maps any given ring_nb_pairs index to the corresponding ring-ring pair.
!! var is a global variable that is incremented during every reset of ring_nb_pairs.
!! Elements in ring_nb_pairs that correspond to pairs that are neighbors equal var. Rest are not equal to var.
!! ring_nb_io stores the actual pairwise neighborlist in a compact form appropriate for I/O.
!! It basically stores the indices of elements of ring_nb_pairs that equal var (i.e. correspond to pairs).
!! pack_ring_nb generates ring_nb_io from ring_nb_pairs. unpack_ring_nb does the reverse.
!! Use of non-var and var in ring_nb_pairs instead of 0 and non-zero is to make searches O(m) [faster than O(m*m)].
module ring_nb
    implicit none

    integer, dimension(:), allocatable, private :: ring_nb_pairs
    integer, private :: var
    integer, dimension(:), allocatable :: ring_nb_io
    integer, dimension(:), allocatable, private :: index_to_pair_ref

contains

    subroutine allocate_ring_nb(m)
        integer, intent(in) :: m
        integer :: alloc_stat, npairs, l, q

        npairs = (m*m - m)/2 ! Number of all possible ring-ring pairs

        allocate (ring_nb_pairs(npairs), index_to_pair_ref(npairs), ring_nb_io(5*m), stat=alloc_stat)
        ! Assuming a max of 10 neighbors per ring. Avoiding double counting : 10*m/2 = 5*m elements for ring_nb_io
        if (alloc_stat /= 0) error stop 'Fatal: Problem while allocating ring_nb'

        ! One-time intialiation
        var = 0
        ring_nb_pairs = var
        ! Loop over all possible pairs to fix index_to_pair_ref array once and for all
        do l = 1, m - 1
            do q = l + 1, m
                index_to_pair_ref(pair_to_index(l, q)) = q
            end do
        end do
    end subroutine allocate_ring_nb

    !NOTE: The following subroutine needs to be executed by a single thread only
    subroutine reset_ring_nb()
        var = var + 1
    end subroutine reset_ring_nb

    integer function pair_to_index(l, q) result(i)
        integer, intent(in) :: l, q
        integer :: row, col

        if (l == q) error stop &
            'Fatal: Row number cannot be equal to column number in a strictly upper triangular matrix'
        row = min(l, q)
        col = max(l, q)
        i = ((col - 1)*(col - 2))/2 + row
    end function pair_to_index

    ! Gives ring indices l and q (l < q) from input ring_nb_pairs index i.
    pure subroutine index_to_pair(i, l, q)
        integer, intent(in) :: i
        integer, intent(out) :: l, q

        q = index_to_pair_ref(i)
        l = i - ((q - 1)*(q - 2))/2
    end subroutine index_to_pair

    subroutine assert_are_nb_rings(l, q)
        integer, intent(in) :: l, q
        integer :: i

        i = pair_to_index(l, q)
        ring_nb_pairs(i) = var
    end subroutine assert_are_nb_rings

    logical function are_nb_rings(l, q) result(check)
        integer, intent(in) :: l, q

        check = (ring_nb_pairs(pair_to_index(l, q)) == var)
    end function are_nb_rings

    subroutine pack_ring_nb()
        integer :: i

        if (size(ring_nb_io) < count(ring_nb_pairs == var)) error stop &
            'Fatal: Too small output buffer for holding the compressed ring-ring neighborhood information'

        ring_nb_io = pack(array=[(i, i=1, size(ring_nb_pairs))], mask=ring_nb_pairs == var, &
                          vector=[(0, i=1, size(ring_nb_io))])
    end subroutine pack_ring_nb

    subroutine unpack_ring_nb()
        integer :: i, valu

        call reset_ring_nb()
        i = 1
        do
            valu = ring_nb_io(i)
            if (valu == 0) exit
            ring_nb_pairs(valu) = var
            i = i + 1
        end do
    end subroutine unpack_ring_nb

end module ring_nb
