module grid_linked_list
       use shared

       implicit none

	   integer, protected:: w, ncell
       double precision, protected :: celli
	   integer, protected:: gridmapsiz
	   integer, dimension(:), allocatable:: gridmap
       integer, dimension(:), allocatable :: bead_nl_head, bead_nl_body ! bead-only neighborlist: head and body
       private :: grid_index
       
    contains
    
    !! Function To Give Cell Index !!
    pure integer function grid_index(ix, iy)
        integer, intent(in):: ix,iy

      	grid_index = 1 + mod(ix-1+w , w) + mod(iy-1+w , w) * w
    end function grid_index

	!!*** Subroutine to set up the cells and to make the list of neighbouring cells through PBC ***!!
	subroutine gridmaps()
	integer:: ix,iy,igridmap,alloc_stat
    double precision:: rcut
    double precision, parameter :: pi=dacos(-1.d0)

    rcut = l0/dsin(pi/n) + rc_adh ! ring diameter + interstitial space
    rcut = rcut/4 !TODO: What decides the divisor? Shouldn't this be a user-parameter or environment variable?
    ! ring diameter estimated from circumcirle of a regular n-gon with side l0
    w=int(box/rcut)
	celli = dble(w/box) ! celli is the inverse of cell length
	if((1.d0/celli).lt.rcut) error stop 'Grid size too small compared to ring diameter + interstitial space'
    ncell=w*w
    gridmapsiz=4*ncell
    
    allocate(gridmap(gridmapsiz), bead_nl_head(ncell), bead_nl_body(m*n), stat=alloc_stat)
        
     if(alloc_stat /= 0) error stop 'Problem while allocating grid_linked_list'
    
	!! Find Half The Nearest Neighbours Of Each Cell !!

	do iy=1,w
        	do ix=1,w

			
			igridmap = (grid_index(ix,iy) - 1) * 4
		
			gridmap(igridmap+1) = grid_index(ix+1,iy)
			gridmap(igridmap+2) = grid_index(ix+1,iy+1)
			gridmap(igridmap+3) = grid_index(ix,iy+1)
			gridmap(igridmap+4) = grid_index(ix-1,iy+1)
		end do
	end do
	end subroutine gridmaps


	!!*** Subroutine to make linked lists & head of chain arrays ***!!
	subroutine links()
	integer:: icell,l,i, bead_index
    double precision :: x3,y3
	
   	!! Zero Head Of Chain Array & List Array !!
    bead_nl_head = 0
    bead_nl_body = 0

	!! Sort All Beads !!

    do l=1,m	 !! Loop over rings
        do i=1,n !! Loop over beads

		    x3 = x(l,i) - box*floor(x(l,i)/box)
            y3 = y(l,i) - box*floor(y(l,i)/box)


			icell = 1 + int(x3 * celli) + int(y3 * celli) * w  ! determining the grid index for a particular bead

			
			bead_index = (l-1)*n+i ! Global serial of the bead : [1,mn]

			bead_nl_body(bead_index) = bead_nl_head(icell)
			bead_nl_head(icell) = bead_index
		end do	
	end do
	end subroutine links

end module grid_linked_list
