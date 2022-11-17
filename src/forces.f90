!! Subroutine for force calculation in a single cell 
	subroutine force(j1)

	use parameters	
	use position_arrays
	use forces

	implicit none
	integer:: i,j,l,j1
        double precision::l1,l2,dx1,dx2,dy1,dy2
   
        ! boundary condition
	
	do l=1,m
          x(l,0) = x(l,n)
          y(l,0) = y(l,n)
          x(l,n+1) = x(l,1)
          y(l,n+1) = y(l,1)
        end do
          

         
     
  	  do l=1,m             ! loop for cell no.
              do i=1,n         ! loop for beads in each cell

                   dx1 = x(l,i-1)-x(l,i)
                   dy1 = y(l,i-1)-y(l,i)
                   dx2 = x(l,i)-x(l,i+1)
                   dy2 = y(l,i)-y(l,i+1)

                   dx1 = dx1 - box*nint(dx1/box)
                   dx2 = dx2 - box*nint(dx2/box)
                   dy1 = dy1 - box*nint(dy1/box)             
                   dy2 = dy2 - box*nint(dy2/box)

                   l1 = dsqrt(dx1*dx1 + dy1*dy1)

                   l2 = dsqrt(dx2*dx2 + dy2*dy2)

                   fx(l,i)=k*(l1-l0)*dx1/l1 - k*(l2-l0)*dx2/l2 - & 
                           0.5d0*p*l0*(dy1/l1 + dy2/l2) 
                     

                   fy(l,i)=k*(l1-l0)*dy1/l1 - k*(l2-l0)*dy2/l2 + &  
                            0.5d0*p*l0*(dx1/l1 + dx2/l2)

               end do
	  end do
		 

	return
	end
    
    
    
     	!!*** Subroutine for the forces of interaction ***!!
	subroutine interaction(j1,frepx,frepy,dx,dy,r,icell,jcell,l)

	use parameters
        use position_arrays
	use forces
	use map_arrays
	use list_arrays
       
       
        implicit none
        integer:: i,j,j1,l,q
        double precision:: r,frepx,frepy,dx,dy,fadhx,fadhy,rlist,ti,tf
	integer:: icell,jcell,jcell0,nabor 
	
        	        f_intx=0.0d0
        	        f_inty=0.0d0
			frpx=0.0d0
			frpy=0.0d0
			fadx=0.0d0
			fady=0.0d0
              
        
        !! Loop Over All Cells !!

	do icell=1,ncell
			!write(*,*)'icell=',icell,'ncell=',ncell

			l=headcell(icell)  !! Highest ring index in a cell(icell)
			
50			if(l.ne.0) then
				!write(*,*)'icell=',icell,'l(headcell(icell)=',l

				i=headbead(icell,l) !! Highest bead index in the part of the ring(l) in a cell(icell)

100				if(i.ne.0) then
					!write(14,*) 'icell=',icell,'l(headcell(icell)=',l,'i(bead(icell,l)=',i

					q=listcell(icell,l)  !! The next lower ring index after l-th ring in a cell(icell) 

200					if((q.ne.0)) then !.and.(q.ne.l)) then
					       
						j=headbead(icell,q) !! The highest bead index in the q-th ring in a cell(icell) 

300						if(j.ne.0) then
											
							dx = x(q,j)-x(l,i)
							dy = y(q,j)-y(l,i)
							dx = dx - box*nint(dx/box)
							dy = dy - box*nint(dy/box)					        
							r = dsqrt(dx*dx + dy*dy)
							
							frepx=0.0d0
							frepy=0.0d0
							fadhx=0.0d0
							fadhy=0.0d0

                      					if(r.lt.rc_rep) then

				          			frepx = -k_rep*(rc_rep-r)*(dx)/r
				          			!fadhx = k_adh*(rc_adh-r)*(dx)/r
					  			frepy = -k_rep*(rc_rep-r)*(dy)/r
				          			!fadhy = k_adh*(rc_adh-r)*(dy)/r

				          			f_intx(l,i) = f_intx(l,i) + frepx !+ fadhx
				          			f_intx(q,j) = f_intx(q,j) - frepx !- fadhx

				          			f_inty(l,i) = f_inty(l,i) + frepy !+ fadhy 
				          			f_inty(q,j) = f_inty(q,j) - frepy !- fadhy

							! 	frpx(l,i) = frpx(l,i) + frepx
							!	frpx(q,j) = frpx(q,j) - frepx
							!	frpy(l,i) = frpy(l,i) + frepy
							!	frpy(q,j) = frpy(q,j) - frepy

							!	fadx(l,i) = fadx(l,i) + fadhx
							!	fadx(q,j) = fadx(q,j) - fadhx
							!	fady(l,i) = fady(l,i) + fadhy
							!	fady(q,j) = fady(q,j) - fadhy


                        				else if((r.le.rc_adh).and.(r.ge.rc_rep)) then

								fadhx = k_adh*(rc_adh-r)*(dx)/r
								fadhy = k_adh*(rc_adh-r)*(dy)/r

								f_intx(l,i) = f_intx(l,i) + fadhx
								f_intx(q,j) = f_intx(q,j) - fadhx

						                f_inty(l,i) = f_inty(l,i) + fadhy 
								f_inty(q,j) = f_inty(q,j) - fadhy

								
								!fadx(l,i) = fadx(l,i) + fadhx
								!fadx(q,j) = fadx(q,j) - fadhx
								!fady(l,i) = fady(l,i) + fadhy
								!fady(q,j) = fady(q,j) - fadhy

       							end if

							j=listbead(q,j) !! the next lower bead index in the q-th ring for the same cell(icell)
							goto 300
						end if

						q=listcell(icell,q) !! the next lower ring index in the same cell(icell)
						goto 200
					end if
	                      !! Loop Over Neighbouring Cells !!

				!	if(ncell.lt.6) goto 700

					jcell0 = 4*(icell-1)     
					
					do nabor=1,4
						jcell = map(jcell0 + nabor)  !!

	!! Loop Over All Rings & Beads Of The Neighbouring Cells!! 

						q=headcell(jcell)
																	
400						if(q.ne.0)then
						    
                                                     if(q.eq.l) goto 600
	
							j=headbead(jcell,q)

500							if(j.ne.0) then

							dx = x(q,j)-x(l,i)
							dy = y(q,j)-y(l,i)
							dx = dx - box*nint(dx/box)
							dy = dy - box*nint(dy/box)					        
							r = dsqrt(dx*dx + dy*dy)
							
							frepx=0.0d0
							frepy=0.0d0
							fadhx=0.0d0
							fadhy=0.0d0

                      					if(r.lt.rc_rep) then

				          			frepx = -k_rep*(rc_rep-r)*(dx)/r
				          			!fadhx = k_adh*(rc_adh-r)*(dx)/r
					  			frepy = -k_rep*(rc_rep-r)*(dy)/r
				          			!fadhy = k_adh*(rc_adh-r)*(dy)/r

				          			f_intx(l,i) = f_intx(l,i) + frepx !+ fadhx
				          			f_intx(q,j) = f_intx(q,j) - frepx !- fadhx

				          			f_inty(l,i) = f_inty(l,i) + frepy !+ fadhy 
				          			f_inty(q,j) = f_inty(q,j) - frepy !- fadhy

							! 	frpx(l,i) = frpx(l,i) + frepx
							!	frpx(q,j) = frpx(q,j) - frepx
							!	frpy(l,i) = frpy(l,i) + frepy
							!	frpy(q,j) = frpy(q,j) - frepy

							!	fadx(l,i) = fadx(l,i) + fadhx
							!	fadx(q,j) = fadx(q,j) - fadhx
							!	fady(l,i) = fady(l,i) + fadhy
							!	fady(q,j) = fady(q,j) - fadhy


                        				else if((r.le.rc_adh).and.(r.ge.rc_rep)) then

								fadhx = k_adh*(rc_adh-r)*(dx)/r
								fadhy = k_adh*(rc_adh-r)*(dy)/r

								f_intx(l,i) = f_intx(l,i) + fadhx
								f_intx(q,j) = f_intx(q,j) - fadhx

						                f_inty(l,i) = f_inty(l,i) + fadhy 
								f_inty(q,j) = f_inty(q,j) - fadhy

								
								!fadx(l,i) = fadx(l,i) + fadhx
								!fadx(q,j) = fadx(q,j) - fadhx
								!fady(l,i) = fady(l,i) + fadhy
								!fady(q,j) = fady(q,j) - fadhy

       							end if


								j=listbead(q,j)!! Considering the next bead of the current ring in the neighbour cell
								goto 500
							end if

600							q=listcell(jcell,q) !!Considering the next ring in the neighbour cell


							goto 400
						end if
					end do


700					i=listbead(l,i)  !! Considering the next bead of the current ring in the current cell
					goto 100
				end if

				l=listcell(icell,l)  !!Considering the next ring in the current cell  	
				goto 50
			end if
	end do		
					
	return
        end                

