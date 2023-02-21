!TODO: parallelize (OMP) later

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
    
    ! Hexatic/Bond-orientational order parameter: h.o.p. It is defined in two ways.
    ! hop1 is from Revalee et al., J. Chem. Phys. 128, 035102 (2008); https://doi.org/10.1063/1.2825300
    ! hop2 is from Loewe et al., Phy. Rev. Lett. 125(3):038003, 2020
    subroutine psi_6(nrings, hop1,hop2)
        use ring_nb, only: are_nb_rings
        integer, intent(in) :: nrings ! Number of rings/cells
        double precision, intent(out) :: hop1, hop2
        integer :: ring1, ring2 ! ring/cell index
        double precision :: re, im, xcm_ring1, ycm_ring1, xcm_ring2, ycm_ring2
        complex, dimension(nrings) :: hop_z_sum ! Stores the complex sum for every cell/ring
        integer, dimension(nrings) :: num_nb ! Number of nearest neighbors
        complex :: z ! Just to store any complex value
        
        hop_z_sum = (0.0, 0.0)
        num_nb = 0
        
        do ring1=1,nrings-1
            call cell_cm(x(:,ring1), y(:,ring1), xcm_ring1, ycm_ring1)
            do ring2=ring1+1,nrings
                if (are_nb_rings(ring1,ring2)) then
                    call cell_cm(x(:,ring2), y(:,ring2), xcm_ring2, ycm_ring2)
                    re = xcm_ring2 - xcm_ring1
                    im = ycm_ring2 - ycm_ring1
                    z = cmplx(re,im)**6 ; z = z/abs(z) ! gives e(i6theta)
                    hop_z_sum(ring1) = hop_z_sum(ring1) + z
                    num_nb(ring1) = num_nb(ring1) + 1
                    hop_z_sum(ring2) = hop_z_sum(ring2) + conjg(z)
                    num_nb(ring2) = num_nb(ring2) + 1
                end if
            end do
        end do

        hop1 = sum(abs(hop_z_sum/num_nb))/nrings
        hop2 = abs(sum(hop_z_sum/num_nb)/nrings)
    end subroutine psi_6


end module analysis
