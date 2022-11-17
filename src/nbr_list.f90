	!!*** Subroutine to set up the cells (grids) and to make the list of neighbouring cells through PBC ***!!
	subroutine maps

        use map_arrays

	implicit none
	integer:: ix,iy,imap,icell

	!! Statement Function To Give Cell Index (in which grid, a bead belongs to) !!

	icell(ix,iy) = 1 + mod(ix-1+w , w) + mod(iy-1+w , w) * w

	!! Find Half The Nearest Neighbours Of Each Cell !!

	do iy=1,w
        	do ix=1,w

			
			imap = (icell(ix,iy) - 1) * 4
		
			map(imap+1) = icell(ix+1,iy)
			map(imap+2) = icell(ix+1,iy+1)
			map(imap+3) = icell(ix,iy+1)
			map(imap+4) = icell(ix-1,iy+1)

			!write(11,*)icell(ix,iy),map(imap+1),map(imap+2),map(imap+3),map(imap+4)
		end do
	end do	

	return
	end




	!!*** Subroutine to make linked lists & head of chain arrays ***!!
	subroutine links(j1,jf) 

	use position_arrays
	use list_arrays
	use map_arrays	

	implicit none
	integer:: j1,jf,icell,l,i,headcell_dummy(ncell)
	double precision:: cell,celli,x3(m,n),y3(m,n),delta
	
   	!! Zero Head Of Chain Array & List Array !!


	headcell      =0 !! head of chain array for the rings in a cell
	headcell_dummy=0
	listcell      =0 !! List array for the rings in a cell
	headbead      =0 !! head of chain array for the beads in a ring in the cell
	listbead      =0

	delta = 0.0d0

	celli = dble(w/box) !! celli is the inverse of cell length
	cell  = 1.0d0/celli !! cell is the cell length

	!if(j1.eq.1) write(*,*) 'cell length=',cell,'No.of cells=',ncell,'inverse length=',celli
	!write(*,*) 'cell length=',cell,'No.of cells=',ncell,'inverse length=',celli

	if(cell.lt.rcut) then
		stop 'cell size to small for cut off'
	end if

	!! Sort All Beads !!

	do l=1,m	 !! Loop over rings
		do i=1,n !! Loop over beads

		        x3(l,i) = x(l,i) - box*floor(x(l,i)/box)
                        y3(l,i) = y(l,i) - box*floor(y(l,i)/box)


8			icell = 1 + int(x3(l,i) * celli) + int(y3(l,i) * celli) * w  !! determining the cell index for a particular bead

			icell_memory(l,i) = icell
			
			
10			listcell(icell,l)     = headcell(icell)     !! The next lower ring index in that cell(icell)
			listbead(l,i)         = headbead(icell,l)   !! The next lower bead index of a part of a ring(l) in a cell(icell) 
			headbead(icell,l)     = i	            !! Highest bead index in the part of the ring(l) in a cell(icell)
			headcell_dummy(icell) = l

		end do	

		headcell = headcell_dummy                   !! Highest ring index in a cell(icell) 

	end do


	return
	end
