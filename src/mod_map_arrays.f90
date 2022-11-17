       module map_arrays
         
       use parameters

       implicit none
	integer,parameter:: w=int(box/rcut), ncell=w*w
	integer,parameter:: mapsiz=4*ncell
	integer:: map(mapsiz)
 
       end module map_arrays
