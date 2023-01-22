! -n | --no-check : Do not check parameters for consistency

program ccd_show_params
    use parameters
    use files
    use utilities, only: cmd_line_flag
    
    call assign_params(params_fname, nocheck=cmd_line_flag('-n') .or. cmd_line_flag ('--no-check'))
    write(*,nml=params)
end program ccd_show_params
