! This module holds code dedicated to visualization of the 2D configurations in gnuplot.
! This is required to create special XY datafile (containing section based discontinuities and virtual beads
!! at box corners) so that gnuplot can handle cells that are broken into pieces (i.e. sections) due to folding
!! under PBC. Gnuplot should be able to produce solid fills and line plots of such cells from the datafile.

! Bragging right: This module could be easily added very late into the project thanks to the extensibility
!! offered by the code-base design and automated build mechanism.


!!!!!! Glossary to understand concepts mentioned below. All w.r.t a single cell !!!!!!

!**Arc : A continuous segment of the cellular periphery, comprised of consecutive beads

!**Lead : First bead in an arc when traversing the arc in ascending order of bead serial (i.e. anticlockwise).

!**Trail : Last bead in an arc when traversing the arc in ascending order of bead serial (i.e. anticlockwise).

!**Section : An isolated part of the cell. Cells unbroken by folding under PBC have only 1 section which is
! confined by a single arc consisting of all the beads in the cell. Cells broken by folding have multiple 
! sections. Each such section is an area enclosed by one or more arcs and one or more edges of the sim(ulation)
! box. If the section is confined by arc(s) and one box edge only, gnuplot can create the section from the arc(s)
! alone by creating a closed polygon from the beads. However, if the section contains two box edges, gnuplot also
! needs the vertex of those two edges along with the arc(s) to form the section polygon.

!**Virtual bead at a box corner : The vertex mentioned above. No bead may actually be positioned there. Yet we need
! to dump the vertex/corner position in the datafile for gnuplot. Hence, virtual.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module gnuplot
    implicit none
    private ! Most of the things here are only for the use of the public components below

    public :: gp_xy_dump

    ! Type containing an integer pair or 2-tuple
    type int_pair
        integer :: x = 0, y = 0

    contains

        private

        procedure :: int_pair_equival
        generic, public :: operator(==) => int_pair_equival
    end type int_pair

    ! Type to hold a continuous segment of a single cell periphery, i.e. an arc
    ! Also has a (singly) linked list nature to store an entire section of the cell
    type cell_arc
        integer :: lead = 0, trail = 0 ! Holds the indices of the two extreme beads on the arc
        type(cell_arc), pointer :: next => null() ! Points to the next arc in the same section
    end type cell_arc

    ! There can be at most 4 sections corresponding to the 4 sim boxes that meet at any corner
    ! A section is identified by a certain folding amount. This is because, originally, i.e.
    !! without folding, all beads of any given cell constitute a single unfragmented cell. The
    !! folding breaks the cell into sections dispersed in space, if originally it was spread
    !! across the edge(s) of the main/replical sim box.
    type cell_section
        type(int_pair) :: sig ! Signature: A section is identified by a certain folding amount
        type(cell_arc), pointer :: head => null() ! Head arc of the (linked) list of arcs within the section
    end type cell_section

contains

    ! This subroutine is similar to xy_dump from mod_files.f90.
    ! Includes additional logic to aid gnuplot create solid fills and line plots
    !! even for cells broken due to folding under PBC.

    ! Dumps xy file for any frame/timestep to be consumed by gnuplot
    ! This routine is threadsafe provided different threads use different `fname`s
    ! x and y are passed as arguments to aid threadsafety
    subroutine gp_xy_dump(fname, boxlen, x, y, title)
        character(len=*), intent(in) :: fname
        double precision, intent(in) :: boxlen
        double precision, dimension(:, :), intent(in) :: x, y
        character(len=*), intent(in), optional :: title
        integer :: fd, l

        open (newunit=fd, file=fname, access='sequential', form='formatted', status='replace', &
              action='write')
        if (present(title)) write (fd, '(a,1x,a)') '#Title:', title
        write (fd, '(a,1x,es23.16)') '#Box:', boxlen
        write (fd, '(a)') '#Column headers:'
        write (fd, '(a,4x,a,4x,a)') 'x', 'y', 'bead'

        cells: do l = 1, size(x, 2)
            ! Two consecutive blank records for separating datasets containing single cell info
            ! Basically creates the demarcator for gnuplot `index`s
            write (fd, '(/)')
            write (fd, '(a,1x,i0)') '#Cell:', l
            call dump_cell_xy(fd, boxlen, x(:, l), y(:, l))
            write (fd, '(a,1x,i0)') '#End_Cell:', l
        end do cells
        close (fd, status='keep')
    end subroutine gp_xy_dump

    ! Dumps folded structured XY for the given cellular configuration.
    ! Structured implies addition of dataset discontinuity (single blank record)
    !! and virtual beads at box corners as necessary for cells broken due to folding.
    subroutine dump_cell_xy(fd, boxlen, x, y)
        use utilities, only: circular_next
        integer, intent(in) :: fd
        double precision, intent(in) :: boxlen
        double precision, dimension(:), intent(in) :: x, y

        ! There can be at most 4 sections corresponding to the 4 sim boxes that meet at a corner.
        type(cell_section), dimension(4) :: secs
        type(int_pair) :: sec_sig ! current signature to determine which section
        integer :: secs_initialized, last_sec ! #sections initialized & index of last section encountered

        integer :: i ! bead serial
        type(cell_arc), pointer :: current_arc

        double precision :: vx, vy ! Holds coordinates of the virtual bead at box corners
        integer :: sec_extrm, last_extrm

        !!! Section construction begins

        secs_initialized = 0
        last_sec = 0

        ! Looping over beads in ascending order
        do i = 1, size(x)
            sec_sig = int_pair(floor(x(i)/boxlen), floor(y(i)/boxlen))
            if (last_sec > 0) then
                if (sec_sig == secs(last_sec)%sig) then
                    current_arc%trail = i
                    cycle
                end if
            end if
            last_sec = findloc(secs(1:secs_initialized)%sig == sec_sig, .true., DIM=1)
            if (last_sec == 0) then
                secs_initialized = secs_initialized + 1
                if (secs_initialized > 4) error stop 'Fatal: More than 4 sections'
                last_sec = secs_initialized
                secs(last_sec)%sig = sec_sig
            end if
            allocate (current_arc)
            current_arc = cell_arc(lead=i, next=secs(last_sec)%head)
            secs(last_sec)%head => current_arc
        end do

        !!! Section construction ends

        ! Structured dumping of folded coordinates:
        !! Section dumps are separated by a single blank record to signify discontinuity. Gnuplot must
        !! not join separate sections otherwise line plots and solid fills will be messed up.
        !! Mainly to aid solid fills, virtual beads at box corners are inserted as and when needed.

        sections: do last_sec = 1, secs_initialized
            sec_sig = secs(last_sec)%sig
            current_arc => secs(last_sec)%head

            ! Traversing the linked list below outputs arcs in reverse order than they were put in the list.
            !! Beads in a section may not maintain any order (ascending or descending) in this approach, hence
            !! messing up any possible continuity. To maintain continuity, the arc direction must also be reversed.
            !! This means swapping the lead and trail beads. This would make all the beads sorted in descending order.

            sec_extrm = current_arc%trail ! Holds first extreme bead (trail/lead) encountered in a section
            last_extrm = 0 ! Holds the last extreme point (trail/lead) encountered in the last arc
            linked_list: do
                if (.not. associated(current_arc)) exit linked_list

                if (last_extrm /= 0) then
                    ! Dump virtual bead at a sim box corner if needed
                    if (need_virtual(x(last_extrm), y(last_extrm), x(current_arc%trail), y(current_arc%trail), &
                                     boxlen, vx, vy)) then
                        write (fd, '(a)') '#Virtual bead follows:'
                        write (fd, '(es23.16,1x,es23.16)') vx, vy
                    end if
                end if

                ! Dump folded XY coordinates for the beads on the arc
                arc_beads: do i = current_arc%trail, current_arc%lead, -1 ! Note trail and lead have been swapped
                    write (fd, '(es23.16,1x,es23.16,1x,i0)') x(i) - sec_sig%x*boxlen, y(i) - sec_sig%y*boxlen, i
                end do arc_beads

                last_extrm = current_arc%lead
                current_arc => current_arc%next
            end do linked_list

            if (secs_initialized > 1) then
                if (last_extrm /= circular_next(sec_extrm, +1, size(x))) then
                    ! No discontinuity (so no need for virtual bead at corner) between beads circular_next to each other
                    if (need_virtual(x(last_extrm), y(last_extrm), x(sec_extrm), y(sec_extrm), boxlen, vx, vy)) then
                        write (fd, '(a)') '#Virtual bead follows:'
                        write (fd, '(es23.16,1x,es23.16)') vx, vy
                    end if
                end if
            end if

            write (fd, *) ! A blank record to signify discontinuity in gnuplot data set

        end do sections

    end subroutine dump_cell_xy

    ! Determine if any of the sim box corner is required to be added as a virtual bead.
    ! Also output the desired corner coordinates.
    ! Takes lead and trail coordinates as well as sim box length as inputs.
    ! Corner is the one that is on the left as we move from lead to trail along a
    !! lead-trail joining arc anticlockwise.
    logical function need_virtual(x_lead, y_lead, x_trail, y_trail, boxlen, corner_x, corner_y)
        double precision, intent(in) :: x_lead, y_lead, x_trail, y_trail, boxlen
        double precision, intent(out) :: corner_x, corner_y

        double precision :: x_lead_fold, y_lead_fold, x_trail_fold, y_trail_fold
        integer :: lead_axis, trail_axis ! Holds reference to the lead/trail's nearest X or Y axis
        ! Reference to the axes is defined as [1, 2, 3, 4] = [Y1, X1, Y2, X2]; e.g. 1 means Y1

        ! Fold the coordinates
        x_lead_fold = modulo(x_lead, boxlen)
        y_lead_fold = modulo(y_lead, boxlen)
        x_trail_fold = modulo(x_trail, boxlen)
        y_trail_fold = modulo(y_trail, boxlen)
        
        lead_axis = minloc(dabs([x_lead_fold, y_lead_fold, boxlen - x_lead_fold, boxlen - y_lead_fold]), DIM=1)
        trail_axis = minloc(dabs([x_trail_fold, y_trail_fold, boxlen - x_trail_fold, boxlen - y_trail_fold]), DIM=1)

        ! Ref. to X axes are even and that to Y are odd. Only odd+odd gives odd.
        need_virtual = mod(lead_axis + trail_axis, 2) /= 0
        ! True only if one of lead and trail is nearest to X and the other to Y axis.
        ! Otherwise, they wouldn't need virtual bead at all.

        ! Determine which corner of the box is appropriate to be inserted as a virtual bead, if needed
        if (need_virtual) then
            if (lead_axis == 2 .and. trail_axis == 1) then ! Lead on X1, Trail on Y1
                corner_x = 0.d0
                corner_y = 0.d0
            else if (lead_axis == 3 .and. trail_axis == 2) then ! Lead on Y2, Trail on X1
                corner_x = boxlen
                corner_y = 0.d0
            else if (lead_axis == 4 .and. trail_axis == 3) then ! Lead on X2, Trail on Y2
                corner_x = boxlen
                corner_y = boxlen
            else ! Lead on Y1, Trail on X2
                corner_x = 0.d0
                corner_y = boxlen
            end if
        end if
    end function need_virtual

    ! Defining equivalence between a pair of int_pair's
    logical elemental function int_pair_equival(val1, val2)
        class(int_pair), intent(in) :: val1, val2

        int_pair_equival = (val1%x == val2%x) .and. (val1%y == val2%y)
    end function int_pair_equival
end module gnuplot
