module utilities
    implicit none
    private :: div_rem
contains

!!!!!!!!/////// Gaussian (Std. Normal) random no. generator \\\\\\\\\\\\\\\\\\\\\\\\\\\
! Polar rejection method (Knop[1969]) related to Box-Muller transform
    subroutine stdnormal(g)
        double precision, intent(OUT) :: g(:)
        double precision :: fac, rsq, v1, v2
        double precision :: rands(2)
        integer :: i, size_g
        double precision, parameter :: small = epsilon(0.0d0)

        size_g = size(g)

        harvest_g_array: do i = 1, size_g, 2
            do
                call random_number(rands)
                v1 = 2.0d0*rands(1) - 1.0d0
                v2 = 2.0d0*rands(2) - 1.0d0
                rsq = v1**2 + v2**2
                if ((rsq < 1.0d0) .and. (rsq > small)) exit
            end do
            fac = DSQRT(-2.0d0*dlog(rsq)/(rsq))
            g(i) = v1*fac
            if (i /= size_g) g(i + 1) = v2*fac
        end do harvest_g_array
    end subroutine stdnormal

! Brief: This subroutine outputs the time spent in secs (cpu & wall clock) since its previous invocation
! Note: All arguments are optional and of type real
! Note: CPU usage = cpu x 100 % / wclock
! Note: #Threads = nint(cpu/wclock)

!TODO: Rename timestamp -> stopwatch. Add split capture mechanism: split_cpu, split_wclock.
! cpu and wclock would give time spent since first/resetting invocation
! Add logical reset flag as well
    subroutine timestamp(cpu, wclock)
        real, intent(out), optional :: cpu, wclock
        integer, save :: sys_clock_count_prev
        integer :: sys_clock_count, sys_clock_max, sys_clock_rate, sys_clock_diff
        real, save :: cpu_sec_prev
        real :: cpu_sec
        integer :: call_count = 1 ! implicit save attribute

        call cpu_time(cpu_sec)

        if (present(cpu)) then
            if (call_count == 1) then
                cpu = 0.0
            else
                cpu = cpu_sec - cpu_sec_prev
            end if
        end if
        cpu_sec_prev = cpu_sec

        call system_clock(sys_clock_count, sys_clock_rate, sys_clock_max)

        if (present(wclock)) then
            if (call_count == 1) then
                sys_clock_diff = 0
            else
                sys_clock_diff = sys_clock_count - sys_clock_count_prev
            end if

            if (sys_clock_diff < 0) then
                wclock = real(sys_clock_diff + sys_clock_max)/sys_clock_rate
            else
                wclock = real(sys_clock_diff)/sys_clock_rate
            end if
        end if
        sys_clock_count_prev = sys_clock_count

        call_count = call_count + 1
    end subroutine timestamp

    ! Transforms real seconds into human readable format
    pure character(len=15) function dhms(sec)
        real, intent(in) :: sec
        integer :: days, hrs, mins, secs

        secs = int(sec)
        call div_rem(secs, 3600*24, days)
        call div_rem(secs, 3600, hrs)
        call div_rem(secs, 60, mins)

        write (dhms, "(i0,'d',1x,i0,'h',1x,i0,'m',1x,i0,'s')") days, hrs, mins, secs
    end function dhms

    pure subroutine div_rem(divided_becomes_remainder, divisor, quotient)
        integer, intent(inout) :: divided_becomes_remainder
        integer, intent(in) :: divisor
        integer, intent(out) :: quotient

        quotient = divided_becomes_remainder/divisor
        divided_becomes_remainder = mod(divided_becomes_remainder, divisor)
    end subroutine div_rem

    ! Outputs the sha1 hash of any given file
    character(len=40) function sha1(fname)
        character(len=*), intent(in) :: fname
        character(len=*), parameter :: tmpfile = '.sha1.tmp'
        integer :: tmpunit, ios

        call execute_command_line("sha1sum "//fname//" 2>/dev/null > "//tmpfile)
        open (newunit=tmpunit, file=tmpfile, status='old', action='read')
        read (tmpunit, '(a)', iostat=ios) sha1
        if (ios /= 0) sha1 = ''
        close (tmpunit, status='delete')
    end function sha1

    function int_to_char(intarg)
        integer, intent(in) :: intarg
        character(len=floor(log10(real(intarg))) + 1) :: int_to_char
        write (int_to_char, '(i0)') intarg
    end function int_to_char

    ! Returns true if `flag` is present as a command line flag/option/argument
    ! Example: cmd_line_flag('--help')
    logical function cmd_line_flag(flag)
        character(len=*), intent(in) :: flag
        character(len=:), allocatable :: cmd_line
        integer :: cmd_line_length
        character(len=*), parameter :: delimiter = ' ' ! Flags in a command line are always delimited

        call get_command(length=cmd_line_length)
        allocate (character(len=cmd_line_length) :: cmd_line)
        call get_command(command=cmd_line)

        cmd_line_flag = index(cmd_line//delimiter, delimiter//flag//delimiter) /= 0
    end function cmd_line_flag

    ! Returns argument of a command line option `opt`. Format of cmd line : <cmd> --<opt>=<arg> ...
    ! Holds the value in `arg` and optionally its length in `length`. `length=0` means no `arg` present.
    ! If option is given multiple times, returns the argument for its last occurence only
    ! Example: cmd_line_opt('--box',boxlen_buffer)
    subroutine cmd_line_opt(opt, arg, length)
        character(len=*), intent(in) :: opt
        character(len=*), intent(out), optional :: arg ! optional in case only length is queried
        integer, intent(out), optional :: length ! Size of arg, optional
        character(len=:), allocatable :: cmd_line
        integer :: cmd_line_length, opt_start_index, opt_end_index, opt_length

        call get_command(length=cmd_line_length)
        allocate (character(len=cmd_line_length) :: cmd_line)
        call get_command(command=cmd_line)

        opt_start_index = index(cmd_line, ' '//opt//'=', back=.true.) + 1
        ! +1 above compensates for the leading ' ' in the substring
        if (opt_start_index == 1) then
            if (present(arg)) arg = ''
            if (present(length)) length = 0
            return
        end if
        opt_end_index = (opt_start_index - 1) + scan(cmd_line(opt_start_index:)//' ', ' ') - 1
        ! -1 above compensates for the trailing ' ' in the substring
        opt_length = len(opt) + 1 ! +1 is to take into account the delimiting '='
        if (present(arg)) arg = cmd_line(opt_start_index + opt_length:opt_end_index)
        if (present(length)) length = opt_end_index - (opt_start_index + opt_length) + 1
    end subroutine cmd_line_opt

    ! Returns `num`th command line argument (i.e. any argument that is not an option starting with - or --)
    ! Holds the value in `arg` and optionally its length in `length`. `length=0` means no `arg` present.
    ! Optionally, `total` holds the total number of such arguments in the command line
    subroutine cmd_line_arg(num, arg, length, total)
        integer, intent(in) :: num
        character(len=*), intent(out), optional :: arg ! optional in case only length is queried
        integer, intent(out), optional :: length ! Size of arg, optional
        integer, intent(out), optional :: total ! Holds the total number of cmd line non-option arguments
        integer :: counter, max_counter, arg_count, arg_length
        character :: prefix

        counter = 0 ! counts the number of all types of arguments including options
        arg_count = 0 ! counts the number of non-option arguments only
        max_counter = command_argument_count()
        ! Loop over all command line arguments
        do
            counter = counter + 1
            call get_command_argument(counter, prefix, arg_length)
            if (prefix /= '-') then
                arg_count = arg_count + 1
                if (arg_count == num) then
                    if (present(arg)) call get_command_argument(counter, arg)
                    if (present(length)) length = arg_length
                    if (present(total)) total = arg_count + (max_counter - counter)
                    return
                end if
            end if
        end do
        if (present(arg)) arg = ''
        if (present(length)) length = 0
        if (present(total)) total = arg_count
    end subroutine cmd_line_arg

    subroutine print_help()
        character(len=32) :: prog_name
        call get_command_argument(0, prog_name)
        call execute_command_line("helpdoc "//trim(prog_name))
    end subroutine print_help

    subroutine help_handler()
        if (cmd_line_flag('-h') .or. cmd_line_flag('--help')) then
            call print_help()
            stop
        end if
    end subroutine help_handler

    ! Returns next/previous element index in a circular array. Default start_index=1
    pure integer function circular_next(curr_index, stride, period, start_index)
        integer, intent(in) :: curr_index, stride, period
        integer, intent(in), optional :: start_index !Default: 1
        integer :: index_offset ! Offset w.r.t circular array starting with index 0

        if (present(start_index)) then
            index_offset = start_index
        else
            index_offset = 1
        end if

        circular_next = modulo(curr_index - index_offset + stride, period) + index_offset
    end function circular_next

    ! Returns eigenvalues and unit eigenvectors of a 2x2 matrix following the trace method
    ! Ref: https://en.wikipedia.org/wiki/Eigenvalue_algorithm#2.C3.972_matrices
    pure subroutine eigen_2x2_mat(matrix, eigenval1, eigenval2, eigenvec1, eigenvec2)
        double precision, dimension(2, 2), intent(in) :: matrix
        double precision, intent(out) :: eigenval1, eigenval2
        double precision, dimension(2), intent(out) :: eigenvec1, eigenvec2
        double precision :: trace, det, gap

        trace = matrix(1, 1) + matrix(2, 2)
        det = matrix(1, 1)*matrix(2, 2) - matrix(1, 2)*matrix(2, 1)
        gap = dsqrt(trace*trace - 4*det)
        eigenval1 = (trace + gap)/2
        eigenval2 = (trace - gap)/2

        eigenvec1 = [matrix(1, 1) - eigenval2, matrix(2, 1)]
        eigenvec2 = [matrix(1, 1) - eigenval1, matrix(2, 1)]

        eigenvec1 = eigenvec1/norm2(eigenvec1)
        eigenvec2 = eigenvec2/norm2(eigenvec2)
    end subroutine eigen_2x2_mat
end module utilities
