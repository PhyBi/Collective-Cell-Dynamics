!! ############################### Program For Distribution Of Various Geometrical Observables ###############################

	implicit none
	integer::i,j,jm,l,t,d,k
	integer,parameter::m=256,n=50
	double precision::Kadh
	double precision::x(m,n+1),y(m,n+1),A(m),mean_area(2000),sum_area,area_frac(2000),area_frac_compl,b,t1,t2,sum_mean_area,&
			  & sum_area_frac,timeavg_area(100),timeavg_area_frac(100),sum_timeavg_area,sum_sq_timeavg_area,&
			  & ensavg_area,ensavg_sq_area,std_dev_area,sum_timeavg_area_frac,ensavg_area_frac,&
			  & sum_sq_timeavg_area_frac,ensavg_sq_area_frac,std_dev_area_frac
	double precision::P(m),sum_perimeter,mean_perimeter(2000),sum_mean_perimeter,timeavg_perimeter(100),sum_timeavg_perimeter,&
			  & ensavg_perimeter,sum_sq_timeavg_perimeter,ensavg_sq_perimeter,std_dev_perimeter
	double precision::SI(m),sum_SI,mean_SI(2000),sum_mean_SI,timeavg_SI(100),sum_timeavg_SI,ensavg_SI,&
			  & sum_sq_timeavg_SI,ensavg_sq_SI,std_dev_SI
	double precision,parameter::box = 46.0d0
	character(len=100):: infile,outfile
	character(len=len('With_noise_coord_wrtSysCOM/Vo0.2_var0.6/Adh0.001_P2.5_vo0.2_sgma0.6_dt0.001/Geometry/'))::&
	& outpath,inpath,outpath2
														
	call cpu_time(t1)
	
	inpath = 'Var0.2/Kadh50/'
	outpath= 'Var0.2/Kadh50/Geometry/'	
	!outpath2= 'With_noise_coord_wrtSysCOM/Vo0.2_var0.4/Avg_Geometry/'

	!Kadh = 5d0	

	open(11,file=trim(adjustl(outpath))//'All_shape_indices_5E6_10ens_adh50_P2.5_var0.2_v_0.2_ite_1E7.dat',&
           & status='unknown',position='append',action='write')
	!open(12,file=trim(adjustl(outpath2))//'Avg_SI_adh20_P2.5_var1.5_v_0.2_ite_1E7.dat',&
         !  & status='unknown',position='append',action='write')


	jm=10
	do j=1,jm

	write (infile,"('adh50_1.5_Noise(var0.2_v_0.2)P2.5_Kspr120_ite_1E7_cnf_',I0,'.dat')") j
	
	open(10,file=trim(adjustl(inpath))//infile,status='unknown',action='read')

	t = 0

	do k=1,2001 
		if(k.ge.1001) then  !! data is considered from 5E6 iteration step.
			t=t+1
			A=0.0d0
			P=0.0d0
			SI=0.0d0		

			do l=1,m
		    		do i=1,n
					read(10,*)b,b,b,x(l,i),y(l,i)				
				end do
				x(l,n+1) = x(l,1)  !! Periodic boundary in a single cell
		 	        y(l,n+1) = y(l,1)
				do i=1,n
					A(l) = A(l) + 0.5d0*(x(l,i)*y(l,i+1) - x(l,i+1)*y(l,i)) 	 !! Area calculation of l-th cell
					P(l) = P(l) + dsqrt((x(l,i)-x(l,i+1))**2 + (y(l,i)-y(l,i+1))**2) !!Perimeter calculation of l-th cell
				end do
				A(l) = dabs(A(l))
				SI(l) = P(l)/dsqrt(A(l))             !! Shape Index
				write(11,*) SI(l)
			end do
		else

			do l=1,m
				do i=1,n
				read(10,*)b,b,b,b,b
				end do
			end do
		end if

	end do

	end do

	close(11)
	
	call cpu_time(t2)
        write(*,*)'time',t2-t1
	end
