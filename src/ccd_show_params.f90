! -n | --no-check : Do not check parameters for consistency

program ccd_show_params
    use parameters
    use files
    use utilities, only: cmd_line_flag
    logical :: check_params_bool
    
    check_params_bool = .not. (cmd_line_flag('-n') .or. cmd_line_flag ('--no-check'))
    call assign_params(params_fname, check=check_params_bool)
    write(*,nml=params)
end program ccd_show_params
