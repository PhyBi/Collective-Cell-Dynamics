!! Subroutine for the movement step in time using Euler algo (No noise)
      subroutine move_deterministic(dt)

       use position_arrays
       use forces

       implicit none
       integer:: i,j,l,idum,j1
       double precision:: c,dt  
       
       c = 0.5d0       ! c is coeff. of viscous damping      
      
       do l=1,m


         do i=1,n

            x(l,i) = x(l,i) + (fx(l,i) + f_intx(l,i))*dt/c
            y(l,i) = y(l,i) + (fy(l,i) + f_inty(l,i))*dt/c
      

         end do

      end do

      return
      end



!! Subroutine for the movement step in time using Euler-Maruyama algo ( including noise)
       subroutine move_noise(dt)

       use position_arrays
       use forces

       implicit none
       integer:: i,j,l,idum,j1
       double precision:: c,dt,noise,tau
       double precision:: vx,vy,wz
       double precision :: theta_x, theta_y, theta_sq_by_4 ! as in Dey arxiv: 1811.06450
          
        c = 0.5d0       ! c is coeff. of viscous damping      
 
       do l=1,m


         do i=1,n

	     
           vx = (fx(l,i) + f_intx(l,i))*dt/c + Vo*mx(l,i)*dt
           vy = (fy(l,i) + f_inty(l,i))*dt/c + Vo*my(l,i)*dt       
           x(l,i) = x(l,i) + vx
           y(l,i) = y(l,i) + vy
            
            
            CALL SYSTEM_CLOCK(COUNT=idum)

            CALL gasdev(noise,idum,mean,var)

            wz = (mx(l,i)*vy - my(l,i)*vx)/tau_align + noise
            theta_x = -my(l,i)*wz*dt
            theta_y = mx(l,i)*wz*dt
            theta_sq_by_4 = (theta_x*theta_x + theta_y*theta_y)/4.0d0
            
            ! Norm preserving rotation of m with ang vel w -> ang dispacement wz*dt
            mx(l,i) = ((1.0d0 - theta_sq_by_4)*mx(l,i) + theta_x)/(1.0d0 + theta_sq_by_4)
            my(l,i) = ((1.0d0 - theta_sq_by_4)*my(l,i) + theta_y)/(1.0d0 + theta_sq_by_4)
            
            

         end do

      end do

      return
      end
