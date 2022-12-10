!! ############################### Program For Various Geometrical Calculation #################################################

	implicit none
	integer::i,j,jm,l,t,d,k
	integer,parameter::m=256,n=50
	double precision::Kadh
	double precision::x(m,n+1),y(m,n+1),A(m),b,t1,t2,sum_area2,&
			  & sum_area,sum_sq_area,mean_area,mean_sq_area,var_area,&
			  & area_frac,sum_area_frac,sum_sq_area_frac,mean_area_frac,mean_sq_area_frac,var_area_frac
	double precision::P(m),sum_perimeter,sum_sq_perimeter,mean_perimeter,mean_sq_perimeter,var_perimeter
	double precision::SI(m),sum_SI,sum_sq_SI,mean_SI,mean_sq_SI,var_SI
	double precision,parameter::box = 46.0d0
	character(len=100):: infile,outfile
	character(len=len('With_noise_coord_wrtSysCOM/Adh0.001_P5_vo0.2_sgma0.5_dt0.001/Geometry_related/')):: outpath,inpath,&
														& outpath2
	call cpu_time(t1)
	
	inpath  = 'Var0.5/Kadh35/'
	outpath = 'Geometry_fluctuations/Var0.5/'
	outpath2= ''

	Kadh = 35d0	

	open(25,file=trim(adjustl(outpath))//'var_mean_area_perim_si_areafrac_5E6_1ens_adh35.dat',&
           & status='unknown',position='append',action='write') 
	!open(12,file=trim(adjustl(outpath2))//'Avg_area_frac_with_adh.dat',&
         !  & status='unknown',position='append',action='write')


	j=1
	sum_area=0.0d0
	sum_perimeter=0.0d0
	sum_SI=0.0d0
	sum_area_frac=0.0d0
	sum_sq_area=0.0d0	
	sum_sq_perimeter=0.0d0	
	sum_sq_SI=0.0d0
	sum_sq_area_frac=0.0d0

	!do j=1,jm

	write (infile,"('adh35_Noise(var0.5_v_0.2)p25_lo0.1_Kspr120_ite_1E7_cnf_',I0,'.dat')") j
	!write (outfile,"('Geometry_adh0.001_1.5_Noise(var0.5_v_0.2)P5_Kspr120_ite_1E7_cnf_',I0,'.dat')") j

	open(10,file=trim(adjustl(inpath))//infile,status='unknown',action='read')
             
	!open(11,file=trim(adjustl(outpath))//outfile,status='unknown',position='append',action='write')
	
	!d=5000
	t = 0
	do k=1,2001 
		if(k.ge.1001) then  !! data is considered from 5E6 iteration step.
			t=t+1
			sum_area2  = 0.0d0
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
		  	 	sum_area = sum_area + A(l)           !! Total area of all the cells at time t
				sum_area2 = sum_area2 + A(l)	     !! Total area of all the cells at time t to calculate area_frac
				sum_perimeter = sum_perimeter + P(l) !! Total perimeter of all the cells at time t
				sum_SI = sum_SI + SI(l)              !! Total shape-index of all the cells at time t     
				sum_sq_area = sum_sq_area + A(l)*A(l)
				sum_sq_perimeter = sum_sq_perimeter + P(l)*P(l)
				sum_sq_SI = sum_sq_SI + SI(l)*SI(l)
			end do
			!if(k.eq.20) write(*,*)A
			area_frac = sum_area2/(box*box)
			sum_area_frac = sum_area_frac + area_frac
			sum_sq_area_frac = sum_sq_area_frac + area_frac*area_frac
		else

			do l=1,m
				do i=1,n
				read(10,*)b,b,b,b,b
				end do
			end do
		end if

	end do

	mean_area = sum_area/(t*m)
	mean_sq_area = sum_sq_area/(t*m)
	var_area = mean_sq_area - mean_area*mean_area
	mean_perimeter = sum_perimeter/(t*m)
	mean_sq_perimeter = sum_sq_perimeter/(t*m)
	var_perimeter = mean_sq_perimeter - mean_perimeter*mean_perimeter
	mean_SI = sum_SI/(t*m)
	mean_sq_SI = sum_sq_SI/(t*m)
	var_SI = mean_sq_SI - mean_SI*mean_SI
	mean_area_frac = sum_area_frac/t
	mean_sq_area_frac = sum_sq_area_frac/t
	var_area_frac = mean_sq_area_frac - mean_area_frac*mean_area_frac

	write(25,*)Kadh,var_perimeter,mean_perimeter,var_area,mean_area,&
		  &var_SI,mean_SI,var_area_frac,mean_area_frac
	!close(11)
	!close(12)
	call cpu_time(t2)
        write(*,*)'time',t2-t1
	end
