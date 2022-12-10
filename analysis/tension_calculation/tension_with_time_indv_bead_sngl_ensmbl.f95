!! ##################### Program For Time-Series Of Tension for individual beads (single ensemble) ##################

	implicit none
	integer::i,j,l,t,tini,tfin,tintrvl,d,k
	integer,parameter::m=256,n=50
	double precision::x(m,n+1),y(m,n+1),A(m),mean_area,sum_area,area_frac,area_frac_compl,b,t1,t2,xtot,ytot
	double precision::P(m),SI(m),sum_perimeter,mean_perimeter,sum_SI,mean_SI,l_neighbors,dl,lo
	double precision,parameter::box = 46.0d0
	character(len=100):: infile,outfile
	character(len=len('Ks120_P2.5_Kadh0.001_Vo0.2_var0.2/CM_related/')):: outpath,inpath
	call cpu_time(t1)
	
	inpath  = 'Var0.6/Kadh30/'
	outpath = 'Tension/'

	lo=0.1d0   !! Equilibrium spring length
	j=1
!	do j=1,5

	write (infile,"('adh30_Noise(var0.6_v_0.2)p25_lo0.1_Kspr120_ite_1E7_cnf_',I0,'.dat')") j
	write (outfile,"('Tension_withtime_frm_95E3_adh30_(var0.6_v_0.2)p25_lo0.1_Kspr120_cnf_',I0,'.dat')") j

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
		
20		do l=1,m
	    		do i=1,n
				read(10,*)b,b,b,x(l,i),y(l,i)				
			end do
			x(l,n+1) = x(l,1)  !! Periodic boundary in a single cell
         	        y(l,n+1) = y(l,1)
			do i=1,n
				l_neighbors = dsqrt((x(l,i+1)-x(l,i))**2 + (y(l,i+1)-y(l,i))**2) !! distance betwn two consequent beads
				dl = l_neighbors - lo 
				write(11,*)d,l,x(l,i),y(l,i),dl
			end do
		end do
		d=d+5000                 
30	end do

!	end do

	call cpu_time(t2)
        write(*,*)'time',t2-t1
	end
