! Help:Begin
! NOTE: This program requires the metadata, last checkpoint and trajectory. Outputs analysis dump.
! Usage: ccd_analysis <metadata file path>
! Help:End

program ccd_analysis
    use analysis
    use utilities, only: help_handler 
    implicit none
    
    character(len=:), allocatable :: metadata_fname
    integer :: metadata_fname_length, exitcode
    integer :: ring, rec_index
    double precision :: msd, shapeind, hexop1, hexop2, vicsekop, areafrac, tension
    double precision :: cell_area, cell_perim, sum_area, cell_tension
    double precision :: cell_vicsekop_x, cell_vicsekop_y, vicsekop_x, vicsekop_y
    
    call help_handler()

    ! Get the metadata file path
    call get_command_argument(1, length=metadata_fname_length, status=exitcode)
    if(exitcode /= 0) error stop 'Fatal: Pass metadata path as argument'
    allocate(character(len=metadata_fname_length) :: metadata_fname)
    call get_command_argument(1, metadata_fname)
    
    call init(metadata_fname)
    
    traj_records: do rec_index=1,recnum
        call traj_read(rec_index, timepoint)
        
        msd=0.d0
        shapeind=0.d0
        vicsekop_x=0.d0
        vicsekop_y=0.d0
        sum_area=0.d0
        tension=0.d0
        
        cells: do ring=1,nrings
            msd = msd + cell_sd(ring)
            
            call cell_perimetry(ring, cell_area, cell_perim, cell_tension)
            sum_area = sum_area + cell_area
            tension = tension + cell_tension
            shapeind = shapeind + cell_perim/dsqrt(cell_area)
            
            call cell_vicsekop(ring, cell_vicsekop_x,cell_vicsekop_y)
            vicsekop_x = vicsekop_x + cell_vicsekop_x
            vicsekop_y = vicsekop_y + cell_vicsekop_y
        end do cells
        
        msd = msd/nrings
        shapeind = shapeind/nrings
        areafrac = sum_area/(box*box)
        tension = tension/nrings
        vicsekop = hypot(vicsekop_x/nrings, vicsekop_y/nrings)
        
        call psi_6(nrings, hexop1,hexop2)
        
        call dump(timepoint, msd, shapeind, hexop1, hexop2, vicsekop, areafrac, tension)
    end do traj_records
    
end program ccd_analysis
