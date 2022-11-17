        program many_cell       ! Main Program Starts

	use parameters
	use position_arrays
	use forces
	use map_arrays
	use list_arrays

	integer:: icell,jcell,l,i,j1,jf,q,j
	double precision:: r,cell,celli,dt,x2(m,n),y2(m,n),ti,tf,frepx,frepy,dx,dy,rint,h,sys_xcm,sys_ycm
	Character(len=len('With_Noise/Coords_shift_wrt_SysCom/Adhesions/Diff_noise_var/Var0.2/Kadh20/Ensmbls_for_alpha2/')) &
	& ::filepath

	filepath = 'With_Noise/Coords_shift_wrt_SysCom/Adhesions/Diff_noise_var/vo0.05_var0.05/'
	
	! open(14,file='Grid1.5_P2.5/Fin_5_10_4_highadh9_2.0_dt0.005_krep100_cnf_5.dat',status='unknown')
	 open(16,file=trim(adjustl(filepath))//'adh0.001_1.5_Noise(var0.05_v_0.05)P2.5_Kspr120_ite_1E7_cnf_2.dat',&
                   & status='unknown',position='append')
         !open(18,file=trim(adjustl(filepath))//'adh0.01_1.5_Noise(var0.5_v_0.2)P4.0_Kspr120_ite_1E7_cnf_1_Fin.dat',&
        ! & status='unknown')


        call cpu_time(ti)
	r=1.0d0
	dt=0.001d0
    tau_align=dt*10
	jf=10050000  !! No. of Iterations
	!jf=1
	call initial(r)   
          do l=1,M  
            do i=1,N
           
               write(16,*)l,i,x(l,i),y(l,i)
                
            end do
          end do
      
	call maps
        call initial_angle
	!write(*,*)'No. of boxes',ncell
	do j1=1,jf

	if(mod(j1,1).eq.0) then  !! Step to calculate all the coordinates w.r.t. System COM.
	sys_xcm = 0.0d0
	sys_ycm = 0.0d0
	do l=1,M
		do i=1,N
			sys_xcm = sys_xcm + x(l,i)
			sys_ycm = sys_ycm + y(l,i)
		end do
	end do
	sys_xcm = sys_xcm/(M*N)
	sys_ycm = sys_ycm/(M*N)
	x = x - sys_xcm
 	y = y - sys_ycm
	end if

        call links(j1,jf)
	call force(j1)
	call interaction(j1,frepx,frepy,dx,dy,rint,icell,jcell,l)

	if(j1.gt.50000) then
							
	call move_noise(dt)

        if(mod(j1,5000).eq.0) then

          do l=1,M  
            do i=1,N
           
               write(16,*)j1-50000,l,i,x(l,i),y(l,i)
                
            end do
          end do
        end if


	else

	call move_deterministic(dt)

        if(j1.eq.50000) then

          do l=1,M  
            do i=1,N
           
               write(16,*)1,l,i,x(l,i),y(l,i)
             
            end do
          end do
        end if


	end if


       if(j1.eq.jf) then
          do l=1,M  
            do i=1,N        
				   !x(l,i) = x(l,i) - box*floor(x(l,i)/box)
				   !y(l,i) = y(l,i) - box*floor(y(l,i)/box)
                !write(18,*)j1,l,i,x(l,i),y(l,i)	              
            end do
          end do
       end if
	
	end do

	call cpu_time(tf)

	write(*,*)'time=',tf-ti

        end program many_cell       ! Main Program Ends
