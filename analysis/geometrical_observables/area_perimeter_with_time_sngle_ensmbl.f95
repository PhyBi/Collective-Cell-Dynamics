!! #################### Program For Evaluation Of Time-Series Of Various Geometrical Observables(single ensemble) ################

	implicit none
	integer::i,j,l,t,tini,tfin,tintrvl,d,k
	integer,parameter::m=256,n=50
	double precision::x(m,n+1),y(m,n+1),A(m),mean_area,sum_area,area_frac,area_frac_compl,b,t1,t2,xtot,ytot
	double precision::P(m),SI(m),sum_perimeter,mean_perimeter,sum_SI,mean_SI
	double precision,parameter::box = 46.0d0
	character(len=100):: infile,outfile
	character(len=len('Ks120_P2.5_Kadh0.001_Vo0.2_var0.2/CM_related/')):: outpath,inpath
	call cpu_time(t1)
	
	inpath  = 'Var0.6/Kadh30/'
	outpath = 'Geometry/Geometry_withTime/'

	j=1
!	do j=1,5

	write (infile,"('adh30_Noise(var0.6_v_0.2)p25_lo0.1_Kspr120_ite_1E7_cnf_',I0,'.dat')") j
	write (outfile,"('Area_perim_withtime_adh30_Noise(var0.6_v_0.2)p25_lo0.1_Kspr120_ite_1E7_cnf_',I0,'.dat')") j

	open(10,file=trim(adjustl(inpath))//infile,status='unknown',action='read')
             
	open(11,file=trim(adjustl(outpath))//outfile,status='unknown',position='append',action='write')
	
	d=5000
	do k=1,2001 
		if(k.gt.19) goto 20
		do l=1,m
			do i=1,n
			read(10,*)b,b,b,b,b
			end do
		end do
		goto 30
20		sum_area  = 0.0d0
		sum_perimeter = 0.0d0
		sum_SI = 0.0d0
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
			sum_perimeter = sum_perimeter + P(l) !! Total perimeter of all the cells at time t
			sum_SI = sum_SI + SI(l)              !! Total shape-index of all the cells at time t     
		end do
		!if(k.eq.20) write(*,*)A
		area_frac = sum_area/(box*box)
		area_frac_compl = 1-area_frac !! Complementary of the area-fraction(vacant space fraction)
		mean_area = sum_area/m        !! mean area of the 256 cells at a time instant t.
		mean_perimeter = sum_perimeter/m  !! mean perimeter of the 256 cells at a time instant t.
		mean_SI = sum_SI/m            !! mean shape-index of the 256 cells at a time instant t.  
		write(11,*)d,area_frac,mean_perimeter,mean_area,mean_SI
		!write(11,*)d,sum_perimeter,sum_area		
		d=d+5000                 
30	end do

!	end do

	call cpu_time(t2)
        write(*,*)'time',t2-t1
	end
