module utilities
implicit none
private :: div_rem 
contains

!!!!!!!!/////// Gaussian random no. generator \\\\\\\\\\\\\\\\\\\\\\\\\\\
! Polar rejection method (Knop[1969]) related to Box-Muller transform
 Subroutine gasdev(g,mean,variance) 
      DOUBLE PRECISION,INTENT(IN)::variance,mean
      DOUBLE PRECISION,INTENT(OUT)::g(:)
      DOUBLE PRECISION:: fac,rsq,v1,v2
      DOUBLE PRECISION:: rands(2)
      INTEGER:: i, size_g
      
      size_g=size(g)
     
     harvest_g_array: DO i=1,size_g,2
        DO
        CALL RANDOM_NUMBER(rands)
        v1=2.0d0*rands(1)-1.0d0
        v2=2.0d0*rands(2)-1.0d0
        rsq=v1**2+v2**2
        if((rsq<1.0d0).AND.(rsq/=0.0d0))EXIT
        ENDDO
        fac=variance*DSQRT(-2.0d0*dlog(rsq)/(rsq))
        g(i)=v1*fac+mean
        if(i/=size_g) g(i+1)=v2*fac+mean 
     END DO harvest_g_array
      END subroutine gasdev



! Brief: This subroutine outputs the time spent in secs (cpu & wall clock) since its previous invocation
! Note: All arguments are optional and of type real
! Note: CPU usage = cpu x 100 % / wclock
! Note: #Threads = nint(cpu/wclock)

subroutine timestamp(cpu, wclock)
    real, intent(out), optional :: cpu, wclock
    integer, save :: sys_clock_count_prev
    integer :: sys_clock_count, sys_clock_max, sys_clock_rate, sys_clock_diff
    real, save :: cpu_sec_prev
    real :: cpu_sec
    integer :: call_count = 1 ! implicit save attribute

    call cpu_time(cpu_sec)
    
    if (present(cpu)) then
        if (call_count == 1) then
            cpu = 0.0
        else
            cpu = cpu_sec - cpu_sec_prev
         end if
    end if
    cpu_sec_prev = cpu_sec

    call system_clock(sys_clock_count, sys_clock_rate, sys_clock_max)

    if (present(wclock)) then
        if (call_count == 1) then
            sys_clock_diff = 0
        else
            sys_clock_diff = sys_clock_count - sys_clock_count_prev
        end if
        
        if (sys_clock_diff < 0) then
            wclock = real(sys_clock_diff + sys_clock_max) / sys_clock_rate
        else
            wclock = real(sys_clock_diff) / sys_clock_rate
        end if
    end if
    sys_clock_count_prev = sys_clock_count
     
    call_count = call_count + 1
end subroutine timestamp

    ! Transforms real seconds into human readable format
    pure character(len=15) function dhms(sec)
        real, intent(in) :: sec
        integer :: days,hrs,mins,secs
     
        secs=int(sec)
        call div_rem(secs, 3600*24, days)
        call div_rem(secs, 3600, hrs)
        call div_rem(secs, 60, mins)   
     
        write(dhms,"(i0,'d',1x,i0,'h',1x,i0,'m',1x,i0,'s')") days, hrs, mins, secs
    end function dhms
     
    pure subroutine div_rem(divided_becomes_remainder, divisor, quotient)
        integer, intent(inout) :: divided_becomes_remainder
        integer, intent(in) :: divisor
        integer, intent(out) :: quotient
    
        quotient = divided_becomes_remainder / divisor
        divided_becomes_remainder = mod(divided_becomes_remainder, divisor)
    end subroutine div_rem

end module utilities

module parameters
    implicit none
    public
    double precision, protected:: k=120.0d0      !  Single cell spring constant
    double precision, protected:: p=25.0d0       !  Single cell internal hydrostatic pressure coefficient
    double precision, protected:: l0=0.1d0       !  Single cell natural spring-length
    double precision, protected:: rc_adh=0.28d0  ! Adhesion interaction cut-off
    double precision, protected:: rc_rep=0.18d0  ! Repulsion interaction cut-off
    double precision, protected:: k_adh=0.001d0   !  Adhesion interaction strength
    double precision, protected:: k_rep=1000.0d0 !  Adhesion interaction strength
    double precision, parameter:: mean=0.0d0     !  Mean of the gaussian white noise
    double precision, protected:: var=0.05d0      !  Variance of the gaussian white noise
    double precision, protected:: Vo=0.05d0       !  Self propulsion of the beads
    double precision, protected:: c = 0.5d0         ! c is coeff. of viscous damping      
    double precision, protected:: dt=0.001d0   ! Integration timestep
    integer, protected:: tau_align=10 ! Tau for Vicsek alignment in multiples of dt
    integer, protected:: jf=10050000  !! No. of Iterations
    integer,protected:: n = 50    ! No. of beads
    integer,protected:: m = 256   ! No. of cell

    namelist /params/ dt, jf, Vo, k_adh, tau_align, var,m,n

    ! Ideally should not have parameter attribute. But that requires allocatable arrays with size : m,n, ncell etc.
    !TODO : Box should be derived and be present inside checkpoint and module state_vars
    double precision,parameter:: box = 46.0d0   !  Box length
    double precision,parameter:: rcut= 1.50d0 
    !TODO: radius should be derived from number of beads per cell and l0
    double precision,parameter:: radius= 1.0d0   ! Initial cell radius
    integer, parameter:: traj_dump_int=100 ! Trajectory file dump interval
    integer, parameter:: status_dump_int=100 ! Status file dump interval
    integer, parameter:: cpt_dump_int=50*traj_dump_int ! Checkpoint file dump interval

    contains
    
    subroutine assign_params(fname)
        character(len=*), intent(in) :: fname
        integer:: fd

        open(newunit=fd,file=fname, access='sequential', form='formatted', status='old', action='read', err=100)
            read(fd, nml=params, err=100, end=100)
        close(fd)
        
        return
100  error stop 'Problem with parameter input file : '//fname
    end subroutine assign_params
end module parameters

module state_vars
    implicit none
    public
    double precision, dimension(:,:), allocatable :: x,y ! Position coordinates of every bead
    double precision, dimension(:,:), allocatable :: mx,my ! Motility unit vector for every bead
    double precision, dimension(:,:), allocatable :: fx,fy ! Intracellular force
    double precision, dimension(:,:), allocatable :: f_rpx,f_rpy ! Intercellular steric repulsion / volume exclusion
    double precision, dimension(:,:), allocatable :: f_adx,f_ady ! Intercellular attraction / adhesion
    integer, dimension(:), allocatable :: prng_seeds ! Stores the state of the P(seudo) R(andom) N(um) G(enerator)

    contains

    integer function allocate_array_stat(cellnum, beadnum)
        integer, intent(in) :: cellnum, beadnum
        integer :: seeds_size

        call random_seed(size = seeds_size)
        
        allocate( &
            x(1:cellnum, 0:beadnum+1), y(1:cellnum, 0:beadnum+1), &
            mx(1:cellnum, 1:beadnum), my(1:cellnum, 1:beadnum), &
            fx(1:cellnum, 1:beadnum), fy(1:cellnum, 1:beadnum), &
            f_rpx(1:cellnum, 1:beadnum), f_rpy(1:cellnum, 1:beadnum), &
            f_adx(1:cellnum, 1:beadnum), f_ady(1:cellnum, 1:beadnum), &
            prng_seeds(seeds_size), &
            stat=allocate_array_stat )
        ! 0:beadnum+1 is to expand a circular array into linear one later in force
    end function allocate_array_stat
end module state_vars

module files
    use iso_fortran_env, only: err_fd => error_unit
    implicit none
    public
    ! File Names
    character(len=*), parameter :: traj_fname='traj.bin'
    character(len=*), parameter :: final_fname='final.xy'
    character(len=*), parameter :: params_fname='params.in'
    character(len=*), parameter :: status_fname='status.live'
    character(len=*), parameter :: cpt_fname='state.cpt'
    ! File Descriptors
    integer :: traj_fd, final_fd, status_fd, cpt_fd
    
    contains
    
    subroutine close_files()
        close(traj_fd, status='keep')
        close(final_fd, status='keep')
        close(status_fd, status='delete')
    end subroutine close_files
    
    subroutine log_this(msg)
        character(len=*), intent(in) :: msg
        character(len=8) :: curr_date
        character(len=10) :: curr_time
        call date_and_time(date=curr_date, time=curr_time)
        write(err_fd,'(/,a,1x,a,1x,a,1x,a)') curr_date(7:8) // '-' // curr_date(5:6) // '-' // curr_date(1:4), &
            curr_time(1:2) // ':' // curr_time(3:4) // ':' // curr_time(5:6), ' => ', msg
    end subroutine log_this
    
end module files

module shared
    use utilities
    use parameters
    use state_vars
    use files
    implicit none
    public
    double precision, parameter::  pi = dacos(-1.0d0)
end module shared

module grid_linked_list
       use shared

       implicit none

	   integer, protected:: w, ncell
	   integer, protected:: mapsiz
	   integer, dimension(:), allocatable:: map
	   integer, dimension(:), allocatable:: headcell
       integer, dimension(:,:), allocatable:: headbead,listcell,listbead,icell_memory
       private :: grid_index
       
    contains
    
    !! Function To Give Cell Index !!
    pure integer function grid_index(ix, iy)
        integer, intent(in):: ix,iy

      	grid_index = 1 + mod(ix-1+w , w) + mod(iy-1+w , w) * w
    end function grid_index

	!!*** Subroutine to set up the cells and to make the list of neighbouring cells through PBC ***!!
	subroutine maps()
	integer:: ix,iy,imap,alloc_stat

    w=int(box/rcut)
    ncell=w*w
    mapsiz=4*ncell
    
    allocate(map(mapsiz), headcell(ncell), headbead(ncell,m),listcell(ncell,m),listbead(m,n),icell_memory(m,n), &
        stat=alloc_stat)
        
     if(alloc_stat /= 0) error stop 'Problem while allocating grid_linked_list'
    
	!! Find Half The Nearest Neighbours Of Each Cell !!

	do iy=1,w
        	do ix=1,w

			
			imap = (grid_index(ix,iy) - 1) * 4
		
			map(imap+1) = grid_index(ix+1,iy)
			map(imap+2) = grid_index(ix+1,iy+1)
			map(imap+3) = grid_index(ix,iy+1)
			map(imap+4) = grid_index(ix-1,iy+1)
		end do
	end do
	end subroutine maps


	!!*** Subroutine to make linked lists & head of chain arrays ***!!
	subroutine links()
	integer:: icell,l,i
    integer, dimension(size(headcell)):: headcell_dummy
	double precision:: cell,celli,delta
    double precision, dimension(size(x,1),size(x,2)):: x3,y3
	
   	!! Zero Head Of Chain Array & List Array !!


	headcell      =0 !! head of chain array for the rings in a cell
	headcell_dummy=0
	listcell      =0 !! List array for the rings in a cell
	headbead      =0 !! head of chain array for the beads in a ring in the cell
	listbead      =0

	delta = 0.0d0

	celli = dble(w/box) !! celli is the inverse of cell length
	cell  = 1.0d0/celli !! cell is the cell length

	if(cell.lt.rcut) error stop 'Grid size too small compared to interaction cut-off'

	!! Sort All Beads !!

	do l=1,m	 !! Loop over rings
		do i=1,n !! Loop over beads

		        x3(l,i) = x(l,i) - box*floor(x(l,i)/box)
                        y3(l,i) = y(l,i) - box*floor(y(l,i)/box)


			icell = 1 + int(x3(l,i) * celli) + int(y3(l,i) * celli) * w  !! determining the cell index for a particular bead

			icell_memory(l,i) = icell
			
			
			listcell(icell,l)     = headcell(icell)     !! The next lower ring index in that cell(icell)
			listbead(l,i)         = headbead(icell,l)   !! The next lower bead index of a part of a ring(l) in a cell(icell) 
			headbead(icell,l)     = i	            !! Highest bead index in the part of the ring(l) in a cell(icell)
			headcell_dummy(icell) = l

		end do	

		headcell = headcell_dummy                   !! Highest ring index in a cell(icell) 

	end do
	end subroutine links

end module grid_linked_list

program ccd_run       ! Main Program Starts
	use shared
    use grid_linked_list

    implicit none

	integer:: l,i,j1,reclen,recnum
    real:: cpusec,wcsec
    logical:: another_run_is_live
    character(len=10):: buffer  ! Internal file

    call assign_params(params_fname)
    call log_this('Run parameters read in')

    write(*,'(/,a,/)') 'Parameters used:'
    write(*,nml=params) ! Dump all params on STDOUT

    call log_this('Creating all necessary arrays (such as positions) @ heap')
    if(allocate_array_stat(m,n) /= 0) error stop 'Array allocation in heap encountered some problem.'

    ! Check if another run is live and Open status file to manifest live status.
    ! This has to be done before the run changes any state / writes anything to memory.
    inquire(file=status_fname, exist=another_run_is_live)
    if(another_run_is_live) error stop 'Another run is going on'
    open(newunit=status_fd,file=status_fname, access='sequential', form='formatted',status='new', &
        asynchronous='yes', action='write')

    write(status_fd,*,asynchronous='yes') jf/status_dump_int
    flush(status_fd) 
    ! This doesn't guarantee data would be written to disk, but warrants it would be available to other processes

    ! Open trajectory file
    call log_this('Initiating trajectory file: '//traj_fname)
	inquire(iolength=reclen) j1, x, y, mx, my, fx, fy, f_rpx, f_rpy, f_adx, f_ady
    open(newunit=traj_fd,file=traj_fname, access='direct', recl=reclen, form='unformatted', status='replace', &
        asynchronous='yes', action='write')

    ! Initialize/Pre-run setup
    call log_this('Pre-run setups...')
    call initial      
	call maps
    call initial_angle
   
    call timestamp()
    
    call log_this('Starting the main run')
	timeseries: do j1=0,jf

    call links
	call force
	call interaction()

        traj_dump: if(mod(j1,traj_dump_int).eq.0) then
             recnum = j1/traj_dump_int + 1
             write(traj_fd, asynchronous='yes', rec=recnum) &
                j1, x, y, mx, my, fx, fy, f_rpx, f_rpy, f_adx, f_ady

            cpt_dump: if(mod(j1,cpt_dump_int).eq.0) then
                ! Complete trajectory-file dumps so far and flush output buffer
                wait(traj_fd)
                flush(traj_fd)

                ! Dumping to .cpt.tmp instead of *.cpt for now
                open(newunit=cpt_fd,file='.cpt.tmp', access='sequential', form='unformatted', &
                    status='replace', action='write')
                    call random_seed(get = prng_seeds)
                    write(cpt_fd) prng_seeds ! Saves current state of the PRNG. To be `put=` in `random_seeds` call...
                    write(cpt_fd) x,y,mx,my ! Saves current state of the physical system
                    write(cpt_fd) recnum ! Last record number of trajectory file, as of now
                close(cpt_fd)
                ! Atomically moving .cpt.tmp to *.cpt now
                call execute_command_line('mv '//'.cpt.tmp '//cpt_fname, wait=.false.)

                write(buffer,'(i0)') j1
                call log_this('Created checkpoint @ timestep = '//trim(buffer))
             end if cpt_dump
        end if traj_dump

        status_dump: if(mod(j1,status_dump_int).eq.0) then
            write(status_fd,*,asynchronous='yes')
            flush(status_fd)
        end if status_dump	

	call move_noise

	end do timeseries

    call timestamp(cpusec,wcsec)

    call log_this('Run complete...writing final config: '//final_fname)
        ! Open final config file
        open(newunit=final_fd, file=final_fname, access='sequential', form='formatted',status='replace', &
            action='write')
          final_dump: do l=1,M
            write(final_fd,'(a,1x,i0)') '#Cell:', l
            do i=1,N        
				   x(l,i) = x(l,i) - box*floor(x(l,i)/box)
				   y(l,i) = y(l,i) - box*floor(y(l,i)/box)
                write(final_fd,*) x(l,i),y(l,i)
            end do
            write(final_fd,'(a,1x,i0,/)') '#End_Cell:', l
          end do final_dump

    call close_files()

    call log_this('Done')

    write(err_fd,'(/,a)') 'Performance:'
    write(err_fd,'(a)') 'CPU = '//dhms(cpusec)
    write(err_fd,'(a)') 'Wallclock = '//dhms(wcsec)
    write(err_fd,'(a,1x,i0,/)') '# Threads = ', nint(cpusec/wcsec) 

end program ccd_run       ! Main Program Ends
	
 	!!*** Subroutine for intercellular forces of interaction ***!!
    ! Below `ring_a` and `ring_b` denote any ring pair within the same cell or grid
    ! `ring_a` and `ring_c` denote any ring pair within neighbouring cells/grids
	subroutine interaction()

	use shared
    use grid_linked_list
          
    implicit none
    integer:: i,j,l,q
    double precision:: r,frepx,frepy,dx,dy,fadhx,fadhy,rlist,factor
	integer:: icell,jcell,jcell0,nabor 
	
			f_rpx=0.0d0
			f_rpy=0.0d0
			f_adx=0.0d0
			f_ady=0.0d0
              
        
        !! Loop Over All Cells !!

	grids: do concurrent (icell=1:ncell)

			l=headcell(icell)  !! Highest ring index in a cell(icell)
			
        ring_a: do	 ! For explanation of ring_a see above
            if(l.eq.0) exit ring_a

				i=headbead(icell,l) !! Highest bead index in the part of the ring(l) in a cell(icell)

			ring_a_beads: do	
                if(i.eq.0) exit ring_a_beads

					q=listcell(icell,l)  !! The next lower ring index after l-th ring in a cell(icell) 

            ring_b: do
					if((q.eq.0)) exit ring_b
					       
						j=headbead(icell,q) !! The highest bead index in the q-th ring in a cell(icell) 

                ring_b_beads: do
						if(j.eq.0) exit ring_b_beads
											
							dx = x(q,j)-x(l,i)
							dy = y(q,j)-y(l,i)
							dx = dx - box*nint(dx/box)
							dy = dy - box*nint(dy/box)					        
							r = dsqrt(dx*dx + dy*dy)
							

                      					if(r.lt.rc_rep) then

				          		factor = (r-rc_rep)/r
				          		frepx = factor*dx
					  			frepy = factor*dy

				          			f_rpx(l,i) = f_rpx(l,i) + frepx 
				          			f_rpx(q,j) = f_rpx(q,j) - frepx 

				          			f_rpy(l,i) = f_rpy(l,i) + frepy 
				          			f_rpy(q,j) = f_rpy(q,j) - frepy 

                        				else if((r.le.rc_adh).and.(r.ge.rc_rep)) then

                                factor = (rc_adh-r)/r
								fadhx = factor*dx
								fadhy = factor*dy

								f_adx(l,i) = f_adx(l,i) + fadhx
								f_adx(q,j) = f_adx(q,j) - fadhx

						        f_ady(l,i) = f_ady(l,i) + fadhy 
								f_ady(q,j) = f_ady(q,j) - fadhy
								
       							end if

							j=listbead(q,j) !! the next lower bead index in the q-th ring for the same cell(icell)
						end do ring_b_beads

						q=listcell(icell,q) !! the next lower ring index in the same cell(icell)
					end do ring_b



	                      !! Loop Over Neighbouring Cells !!

					jcell0 = 4*(icell-1)     
					
					do nabor=1,4
						jcell = map(jcell0 + nabor)  !!

	!! Loop Over All Rings & Beads Of The Neighbouring Cells!! 

						q=headcell(jcell)
																	
                ring_c: do
						if(q.eq.0) exit ring_c
						    
                                                     if(q.eq.l) then
                                                        !!Considering the next ring in the neighbour cell
                                                        q=listcell(jcell,q)
                                                        cycle ring_c
                                                     end if
	
							j=headbead(jcell,q)

                    ring_c_beads: do
							if(j.eq.0) exit ring_c_beads

							dx = x(q,j)-x(l,i)
							dy = y(q,j)-y(l,i)
							dx = dx - box*nint(dx/box)
							dy = dy - box*nint(dy/box)					        
							r = dsqrt(dx*dx + dy*dy)
							

                      					if(r.lt.rc_rep) then

				          		factor = (r-rc_rep)/r
				          		frepx = factor*dx
					  			frepy = factor*dy


				          			f_rpx(l,i) = f_rpx(l,i) + frepx 
				          			f_rpx(q,j) = f_rpx(q,j) - frepx 

				          			f_rpy(l,i) = f_rpy(l,i) + frepy 
				          			f_rpy(q,j) = f_rpy(q,j) - frepy 

                        				else if((r.le.rc_adh).and.(r.ge.rc_rep)) then

                                factor = (rc_adh-r)/r
								fadhx = factor*dx
								fadhy = factor*dy

								f_adx(l,i) = f_adx(l,i) + fadhx
								f_adx(q,j) = f_adx(q,j) - fadhx

						        f_ady(l,i) = f_ady(l,i) + fadhy 
								f_ady(q,j) = f_ady(q,j) - fadhy
								
       							end if


								j=listbead(q,j)!! Considering the next bead of the current ring in the neighbour cell
							end do ring_c_beads

							q=listcell(jcell,q) !!Considering the next ring in the neighbour cell


						end do ring_c
					end do


					i=listbead(l,i)  !! Considering the next bead of the current ring in the current cell
				end do ring_a_beads

				l=listcell(icell,l)  !!Considering the next ring in the current cell  	
			end do ring_a
	end do grids		
	
    f_adx = k_adh*f_adx
    f_ady = k_adh*f_ady
    f_rpx = k_rep*f_rpx
    f_rpy = k_rep*f_rpy    

	return
        end subroutine interaction            

	

!!! Subroutine for random initial configurations
subroutine initial

use shared
 
implicit none
double precision:: a,b,a1,b1,g,dr,q,theta1
integer:: l,k1,i,t
      DOUBLE PRECISION:: rands(2)
      double precision, dimension(size(x,dim=1)) :: x1, y1


q = 0.5d0*dsqrt(m*7.3d0)
dr = 0.2d0
t = 0


    do
        CALL RANDOM_NUMBER(rands)
    a = q*2.0d0*rands(1)
    b = q*2.0d0*rands(2)
    if((a.ge.(radius+0.05d0)).and.(b.ge.(radius+0.05d0))) exit
    end do


   x1(1) = a
   y1(1) = b


do l=2,m
        lcell: do
    do
        CALL RANDOM_NUMBER(rands)
    a1 = q*2.0d0*rands(1)
    b1 = q*2.0d0*rands(2)
    if((a1.ge.(radius+0.05d0)).and.(b1.ge.(radius+0.05d0))) exit
    end do

   x1(l) = a1
   y1(l) = b1

    do k1=1,l-1

       g = dsqrt((x1(l)-x1(k1))**2 +(y1(l)-y1(k1))**2)

       if (g.lt.(2*radius+dr)) cycle lcell
    end do
        exit lcell
        end do lcell
end do

       theta1 = 0.0d0
       
  
     do l=1,m  
       do i=1,n
        CALL RANDOM_NUMBER(rands)
        
           x(l,i) = radius*cos(theta1) + 0.01d0*(2.0d0*rands(1)-1.0d0) + x1(l)
           y(l,i) = radius*sin(theta1) + 0.01d0*(2.0d0*rands(2)-1.0d0) + y1(l)
           theta1 = theta1 + 2.0d0*pi/n
                                                         
       end do
     end do

	return
	end subroutine initial


!! Subroutine for intracellular forces 
	subroutine force

	use shared	
	

	implicit none
	integer:: i,j,l
        double precision::l1,l2,dx1,dx2,dy1,dy2
   
        ! boundary condition
	
          x(:,0) = x(:,n)
          y(:,0) = y(:,n)
          x(:,n+1) = x(:,1)
          y(:,n+1) = y(:,1)
          

         
     
  	  do concurrent (l=1:m, i=1:n)

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
		 

	return
	end subroutine force


!! Subroutine for the movement step in time using Euler-Maruyama algo ( including noise)
       subroutine move_noise

       use shared
       

       implicit none
       integer:: i,j,l
       double precision:: vx,vy,wz
       double precision :: theta_x, theta_y, theta_sq_by_4 ! as in Dey arxiv: 1811.06450
       double precision :: noise(m,n), g(m*n)
          
 
             CALL gasdev(g,mean,var)
             noise=reshape(g, [m,n]) ! reshapes 1D array g into 2D

       do concurrent (l=1:m, i=1:n)
           vx = (fx(l,i) + f_adx(l,i) + f_rpx(l,i))*dt/c + Vo*mx(l,i)*dt
           vy = (fy(l,i) + f_ady(l,i) + f_rpy(l,i))*dt/c + Vo*my(l,i)*dt
           x(l,i) = x(l,i) + vx
           y(l,i) = y(l,i) + vy
            
      
            wz = (mx(l,i)*vy - my(l,i)*vx)/(tau_align*dt) + noise(l,i)
            theta_x = -my(l,i)*wz*dt
            theta_y = mx(l,i)*wz*dt
            theta_sq_by_4 = (theta_x*theta_x + theta_y*theta_y)/4.0d0
            
            ! Norm preserving rotation of m with ang vel w -> ang dispacement wz*dt
            mx(l,i) = ((1.0d0 - theta_sq_by_4)*mx(l,i) + theta_x)/(1.0d0 + theta_sq_by_4)
            my(l,i) = ((1.0d0 - theta_sq_by_4)*my(l,i) + theta_y)/(1.0d0 + theta_sq_by_4)
      end do

      return
      end subroutine move_noise


!! Subroutine for initializaion of noise-angle
       subroutine initial_angle           

        use shared

        implicit none

        double precision :: m_tmp,mx_tmp,my_tmp
        integer :: i,l
        DOUBLE PRECISION:: rands(2)
        
        do l=1,m
            do i=1,n
               CALL RANDOM_NUMBER(rands)

                mx_tmp = 2.0d0*rands(1) - 1.0d0
                my_tmp = 2.0d0*rands(2) - 1.0d0
                m_tmp = dsqrt(mx_tmp*mx_tmp + my_tmp*my_tmp)
                mx(l,i)=mx_tmp/m_tmp
                my(l,i)=my_tmp/m_tmp     
            end do
        end do

        return
        end subroutine initial_angle