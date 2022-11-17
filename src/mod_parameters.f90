       module parameters

       implicit none
	integer,parameter:: n = 50    ! No. of beads
	integer,parameter:: m = 256   ! No. of cell
    integer,parameter:: s = n+1
    double precision,parameter:: box = 46.0d0   !  Box length
	double precision,parameter:: rcut= 1.50d0 
	double precision,parameter:: k=120.0d0      !  Single cell spring constant
	double precision,parameter:: p=25.0d0       !  Single cell internal hydrostatic pressure coefficient
	double precision,parameter:: l0=0.1d0       !  Single cell natural spring-length
	double precision,parameter:: rc_adh=0.28d0  ! Adhesion interaction cut-off
	double precision,parameter:: rc_rep=0.18d0  ! Repulsion interaction cut-off
	double precision,parameter:: k_adh=0.001d0   !  Adhesion interaction strength
	double precision,parameter:: k_rep=1000.0d0 !  Adhesion interaction strength
	double precision,parameter:: mean=0.0d0     !  Mean of the gaussian white noise
	double precision,parameter:: var=0.05d0      !  Variance of the gaussian white noise
	double precision,parameter:: Vo=0.05d0       !  Self propulsion of the beads
    double precision:: tau_align ! Tau for Vicsek alignment
       end module parameters

