module analysis
    use state_vars
    use utilities, only: circular_next, eigen_2x2_mat
    implicit none
    
    contains
    
    pure subroutine cell_cm(xbeads, ybeads, xcm, ycm)
        double precision, dimension(:), intent(in) :: xbeads, ybeads
        double precision, intent(out) :: xcm, ycm
        integer :: nbeads
        
        nbeads = size(xbeads)
        xcm = sum(xbeads)/nbeads
        ycm = sum(ybeads)/nbeads
    end subroutine cell_cm
    
    pure subroutine cell_perimetry(xbeads, ybeads, area, perimeter)
        double precision, dimension(:), intent(in) :: xbeads, ybeads
        double precision, intent(out) :: area, perimeter
        integer :: nbeads, this_bead, next_bead
        double precision :: dx, dy

        nbeads=size(xbeads)
        area=0.0d0
        perimeter=0.0d0
        do this_bead=1,nbeads
            next_bead=circular_next(this_bead,+1,nbeads)
            
            area=area+xbeads(this_bead)*ybeads(next_bead)-xbeads(next_bead)*ybeads(this_bead)
            
            perimeter=perimeter+hypot(xbeads(this_bead)-xbeads(next_bead),ybeads(this_bead)-ybeads(next_bead))
        end do
        area=dabs(0.5d0*area)    
    end subroutine cell_perimetry

    pure subroutine cell_shape(xbeads, ybeads, major_axis, minor_axis)
        double precision, dimension(:), intent(in) :: xbeads, ybeads
        double precision, dimension(:), intent(out) :: major_axis, minor_axis
        double precision :: xcm, ycm, shape_tensor(2,2), eval1, eval2, evec1(2), evec2(2)
        integer :: nbeads
        
        call cell_cm(xbeads, ybeads, xcm, ycm)
        nbeads=size(xbeads)
        shape_tensor(1,1) = sum((xbeads-xcm)*(xbeads-xcm))/nbeads
        shape_tensor(2,2) = sum((ybeads-ycm)*(ybeads-ycm))/nbeads
        shape_tensor(2,1) = sum((ybeads-ycm)*(xbeads-xcm))/nbeads
        shape_tensor(1,2) = shape_tensor(2,1)
        
        call eigen_2x2_mat(shape_tensor,eval1,eval2,evec1,evec2)
        
        if(eval1 > eval2) then
            major_axis = evec1
            minor_axis = evec2
        else
            major_axis = evec2
            minor_axis = evec1            
        endif
    end subroutine cell_shape
    
    ! Hexatic/Bond-orientational order parameter: h.o.p
    ! Gives mod(sum(psi_6)/npairs)
    !TODO: parallelize (OMP) later
    subroutine psi_6(hop)
        use ring_nb, only: are_nb_rings
        double precision, intent(out) :: hop
        integer :: nrings, ring1, ring2 ! ring/cell index
        double precision :: re, im, abs_val, xcm_ring1, ycm_ring1, xcm_ring2, ycm_ring2
        complex :: hop_z_sum
        
        nrings = size(x,2)
        hop_z_sum = (0.0, 0.0)
        
        do ring1=1,nrings-1
            call cell_cm(x(:,ring1), y(:,ring1), xcm_ring1, ycm_ring1)
            do ring2=ring1+1,nrings
                if (are_nb_rings(ring1,ring2)) then
                    call cell_cm(x(:,ring2), y(:,ring2), xcm_ring2, ycm_ring2)
                    re = xcm_ring2 - xcm_ring1
                    im = ycm_ring2 - ycm_ring1
                    abs_val = abs(cmplx(re,im))
                    hop_z_sum = hop_z_sum + (cmplx(re,im)/abs_val)**6
                end if
            end do
        end do

        hop = hop_z_sum/count(ring_nb_io/=0)
    end subroutine psi_6


end module analysis
