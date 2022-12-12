module grid_linked_list
       use shared

       implicit none

	   integer, protected:: w, ncell
       double precision, protected :: celli
	   integer, protected:: gridmapsiz
	   integer, dimension(:), allocatable:: gridmap
	   integer, dimension(:), allocatable:: headcell
       integer, dimension(:,:), allocatable:: headbead,listcell,listbead
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

    w=int(box/rcut)
    !TODO: rcut should be derived from number of beads per cell and l0
	celli = dble(w/box) !! celli is the inverse of cell length
	if((1.d0/celli).lt.rcut) error stop 'Grid size too small compared to interaction cut-off'
    ncell=w*w
    gridmapsiz=4*ncell
    
    allocate(gridmap(gridmapsiz), headcell(ncell), headbead(ncell,m),listcell(ncell,m),listbead(m,n), stat=alloc_stat)
        
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
	integer:: icell,l,i
    integer, dimension(size(headcell)):: headcell_dummy
    double precision, dimension(size(x,1),size(x,2)):: x3,y3
	
   	!! Zero Head Of Chain Array & List Array !!


    do icell=1,ncell
    headcell(icell)=0 !! head of chain array for the rings in a cell
	headcell_dummy(icell)=0
	listcell(icell,:)=0 !! List array for the rings in a cell
	headbead(icell,:)=0 !! head of chain array for the beads in a ring in the cell
    end do

	!! Sort All Beads !!

    do l=1,m	 !! Loop over rings
		listbead(l,:)=0
        do i=1,n !! Loop over beads

		        x3(l,i) = x(l,i) - box*floor(x(l,i)/box)
                        y3(l,i) = y(l,i) - box*floor(y(l,i)/box)


			icell = 1 + int(x3(l,i) * celli) + int(y3(l,i) * celli) * w  !! determining the cell index for a particular bead

			
			
			listcell(icell,l)     = headcell(icell)     !! The next lower ring index in that cell(icell)
			listbead(l,i)         = headbead(icell,l)   !! The next lower bead index of a part of a ring(l) in a cell(icell) 
			headbead(icell,l)     = i	            !! Highest bead index in the part of the ring(l) in a cell(icell)
			headcell_dummy(icell) = l

		end do	

		headcell = headcell_dummy                   !! Highest ring index in a cell(icell) 

	end do
	end subroutine links

end module grid_linked_list
