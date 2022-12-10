!######################################### Program For MSD Calculation ###############################################

	implicit none

	integer:: i,j,t,k,l,a,d,jm
	integer,parameter:: m=256
	double precision :: x1(0:20000,m),y1(0:20000,m),r_sq,r_sq_avg(20000,50),r_quad_avg(20000),r_quad,alpha,msd_avg
	character(len=90):: infile,outfile
	character(len=len('Test_dt0.001_highadh_kspr120/sgma3.0/Adh0.05_P2.5_Vo0.2/CM_related/MSD2/')):: outpath,inpath

	jm = 10  ! No. of ensembles

	inpath  = 'Var0.2/Kadh40/CM_related/'
	outpath = 'Var0.2/Kadh40/CM_related/MSD/'


	open(13,file=trim(adjustl(outpath))//'Stdy9E4_Avg10_MSD_CM_adh40_(var0.2_v_0.2)P2.5_ite_1E7.dat',&
           & status='unknown',position='append')

	do j=1,jm

	!if((j.eq.1).or.(j.eq.2).or.(j.eq.4).or.(j.eq.8).or.(j.eq.10)) goto 10
 	!if(j.eq.5) goto 10
	write (infile,"('CM_adh40_1.5_Noise(var0.2_v_0.2)P2.5_Kspr120_ite_1E7_cnf_',I0,'.dat')") j
	!write (outfile,"('Steady_CorrectMSD_CM_adh1.0_1.5_Noise(sigma3.0_v_0.2)P2.5_Kspr120_ite_10_7_raw_cnf_',I0,'.dat')") j
	!write (outfile,"('Stdy_msd_CM_adh4.0_1.5_Noise(sigma3.0_v_0.2)raw_cnf_',I0,'.dat')") j


        open(11,file=trim(adjustl(inpath))//infile,status='old',action='read')
   
       ! open(12,file=trim(adjustl(outpath))//outfile,status='unknown',position='append',action='write')
      


			t=0
			do k=1,2001 !200

			  if(k.gt.19) t=t+1  !k=11 step is the 50000 iterations with noise. we set it as t0 state.
                                             !k=21 step is the 1E5 iterations with noise. we set it as t0 state.
			 do i=1,m
			  read(11,*)l,a,x1(t,i),y1(t,i)
			 end do 

		       end do


				d = 5000
				do t=1,1982 ! tfinal = 2000 - 10 (10 data points for equilibration. t=10 is our t0)
                                            ! tfinal = 2000 - 20 (20 data points for equilibration. t=10 is our t0)
				    r_sq_avg(t,j) = 0.0d0

				  do i=1,m
				     !r_sq = (dsqrt(x1(t,i)**2 + y1(t,i)**2) - dsqrt(x1(0,i)**2 + y1(0,i)**2))**2

				     r_sq = (x1(t,i) - x1(0,i))**2 + (y1(t,i) - y1(0,i))**2

				     r_sq_avg(t,j) = r_sq_avg(t,j) + r_sq 

				  end do

				   r_sq_avg(t,j) = r_sq_avg(t,j)/m

				  ! write(12,*)d,r_sq_avg(t,j)
				  
				   d=d+5000
				end do

				
10	end do



        d=5000
	do t=1,1982
		msd_avg = 0.0d0
		do j=1,jm

			!if((j.eq.1).or.(j.eq.2).or.(j.eq.4).or.(j.eq.8).or.(j.eq.10)) goto 11
			!if(j.eq.5) goto 11
			msd_avg = msd_avg + r_sq_avg(t,j)

11		end do

		msd_avg = msd_avg/jm

		write(13,*)d,msd_avg
                d=d+5000
	end do
end
