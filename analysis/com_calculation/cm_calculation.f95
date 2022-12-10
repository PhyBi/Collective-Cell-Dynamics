!! ############################### Program For Center Of Mass Coordinate Calculation #################################################

	implicit none
	integer::i,j,l,t,tini,tfin,tintrvl,jm
	integer,parameter::m=256,n=50
	double precision::x(m,n),y(m,n),x_cm(m),y_cm(m),a,t1,t2,xtot,ytot
	character(len=100):: infile,outfile
	character(len=len('Ks120_P2.5_Kadh0.001_Vo0.2_var0.2/CM_related/')):: outpath,inpath
	call cpu_time(t1)
	
	jm = 10 !! # of ensembles
	inpath  = 'vo0.05_var0.05/'
	outpath = 'vo0.05_var0.05/CM_related/'

	!j=1
	do j=1,jm

	write (infile,"('adh0.001_1.5_Noise(var0.05_v_0.05)P2.5_Kspr120_ite_1E7_cnf_',I0,'.dat')") j
	write (outfile,"('CM_adh0.001_1.5_Noise(var0.05_v_0.05)P2.5_Kspr120_ite_1E7_cnf_',I0,'.dat')") j

	open(10,file=trim(adjustl(inpath))//infile,status='unknown',action='read')
             
	open(11,file=trim(adjustl(outpath))//outfile,status='unknown',position='append',action='write')
            
        
	tini=5000
	tfin=10000000
	tintrvl=5000

			      do l=1,m
                                     xtot=0.0d0
                                     ytot=0.0d0
				      do i=1,n
			      		      read(10,*)a,a,a,x(l,i),y(l,i)
					      xtot = xtot + x(l,i)
					      ytot = ytot + y(l,i)	
	
				      end do
				      xtot=xtot/n
				      ytot=ytot/n
                                      x_cm(l)=xtot
                                      y_cm(l)=ytot	
				      write(11,*)1,l,x_cm(l),y_cm(l)
			      end do 


	do t=tini,tfin,tintrvl

			    if(t.eq.tini) then
                               do l=1,m
                                 do i=1,n
                               
                                 end do
                               end do
			    end if

			      do l=1,m
                                     xtot=0.0d0
                                     ytot=0.0d0
				      do i=1,n
			      		      read(10,*)a,a,a,x(l,i),y(l,i)
					      xtot = xtot + x(l,i)
					      ytot = ytot + y(l,i)	
	
				      end do
				      xtot=xtot/n
				      ytot=ytot/n
                                      x_cm(l)=xtot
                                      y_cm(l)=ytot	
				      write(11,*)t,l,x_cm(l),y_cm(l)
			      end do 
	end do 

	end do

	call cpu_time(t2)
        write(*,*)'time',t2-t1
	end
