!! ############################### Program For Avg-Tension(over cells & beads) Calculation #########################

	implicit none
	integer::i,j,l,t,tini,tfin,tintrvl,d,k
	integer,parameter::m=256,n=50
	double precision::x(m,n+1),y(m,n+1),A(m),mean_area,sum_area,area_frac,area_frac_compl,b,t1,t2,xtot,ytot
	double precision::l_neighbors,dl,lo,sum_percell_dl,sum_allcells_dl,avg_allcells_dl
	double precision,parameter::box = 46.0d0
	character(len=100):: infile,outfile
	character(len=len('Ks120_P2.5_Kadh0.001_Vo0.2_var0.2/CM_related/')):: outpath,inpath
	call cpu_time(t1)
	
	inpath  = 'Var0.6/Kadh30/'
	outpath = 'Tension/Tension_with_time_single_ens/'

	lo=0.1d0   !! Equilibrium spring length
	j=1
!	do j=1,5

	write (infile,"('adh30_Noise(var0.6_v_0.2)p25_lo0.1_Kspr120_ite_1E7_cnf_',I0,'.dat')") j
	write (outfile,"('Avg_tension_withtime_1ens_adh30_Noise(var0.6_v_0.2)P2.5_Kspr120_ite_1E7_cnf_',I0,'.dat')") j

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
20		sum_allcells_dl  = 0.0d0
		do l=1,m
			sum_percell_dl=0.0d0
	    		do i=1,n
				read(10,*)b,b,b,x(l,i),y(l,i)				
			end do
			x(l,n+1) = x(l,1)  !! Periodic boundary in a single cell
         	        y(l,n+1) = y(l,1)
			do i=1,n
				l_neighbors = dsqrt((x(l,i+1)-x(l,i))**2 + (y(l,i+1)-y(l,i))**2) !! distance betwn two consequent beads
				dl = l_neighbors - lo 
				sum_percell_dl = sum_percell_dl + dl
			end do
			sum_percell_dl = sum_percell_dl/n
	  	 	sum_allcells_dl = sum_allcells_dl + sum_percell_dl           !! Total area of all the cells at time t
		end do
		avg_allcells_dl = sum_allcells_dl/m
		write(11,*)d,avg_allcells_dl
		d=d+5000                 
30	end do

!	end do

	call cpu_time(t2)
        write(*,*)'time',t2-t1
	end
