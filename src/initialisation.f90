!!! Subroutine for random initial configurations
subroutine initial(r)

use position_arrays
 
implicit none
double precision:: a,b,a1,b1,g,r,dr,q,pi,theta1,ran2
integer:: l,k1,i,t,idum

q = 0.5d0*dsqrt(m*7.3d0)
dr = 0.2d0
t = 0

CALL SYSTEM_CLOCK(COUNT=idum)

!call srand(653214)

4   a = q*2.0d0*ran2(idum)
    b = q*2.0d0*ran2(idum)
!4   a = q*2.0d0*rand(0)
!    b = q*2.0d0*rand(0)

   if((a.lt.(r+0.05d0)).or.(b.lt.(r+0.05d0))) goto 4

   x1(1) = a
   y1(1) = b


do l=2,m

5   a1 = q*2.0d0*ran2(idum)
    b1 = q*2.0d0*ran2(idum)
!5   a1 = q*2.0d0*rand(0)
!    b1 = q*2.0d0*rand(0)

    if((a1.lt.(r+0.05d0)).or.(b1.lt.(r+0.05d0))) goto 5

   x1(l) = a1
   y1(l) = b1

    do k1=1,l-1

       g = dsqrt((x1(l)-x1(k1))**2 +(y1(l)-y1(k1))**2)

       if (g.lt.(2*r+dr)) then
      !  t=t+1
      !  write(*,*)t
        goto 5
       end if
    end do

end do

       pi = 4.0d0*datan(1.0d0)
       theta1 = 0.0d0
       
  
     do l=1,m  
       do i=1,n
        
           x(l,i) = r*cos(theta1) + 0.01d0*(2.0d0*ran2(idum)-1.0d0) + x1(l)
           y(l,i) = r*sin(theta1) + 0.01d0*(2.0d0*ran2(idum)-1.0d0) + y1(l)
           !x(l,i) = r*cos(theta1) + x1(l) !+ 0.01d0*(2.0d0*rand(0)-1.0d0) 
           !y(l,i) = r*sin(theta1) + y1(l) !+ 0.01d0*(2.0d0*rand(0)-1.0d0)  
           
           theta1 = theta1 + 2.0d0*pi/n
                                                         
       end do
     end do

	return
	end     
    
    
    
    !! Subroutine for initializaion of the direction of motility vector
       subroutine initial_angle           

        use position_arrays

        implicit none

        double precision :: pi,ran2,m_tmp,mx_tmp,my_tmp
        integer :: idum,i,l

         pi = 4.0d0*datan(1.0d0)
       !  idum = 564326
        CALL SYSTEM_CLOCK(COUNT=idum)
         
        do l=1,m
            do i=1,n
                mx_tmp = 2.0d0*ran2(idum) - 1.0d0
                my_tmp = 2.0d0*ran2(idum) - 1.0d0
                m_tmp = dsqrt(mx_tmp*mx_tmp + my_tmp*my_tmp)
                mx(l,i)=mx_tmp/m_tmp
                my(l,i)=my_tmp/m_tmp     
            end do
        end do

        return
        end

