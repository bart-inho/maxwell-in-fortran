program maxwell

    use iso_fortran_env
    use all

    implicit none

    type(everything) :: all
    integer :: nt = 2000, t = 0 ! Number of time steps

    call all%EverythingInit() ! Initialize everything

    do t = 1, nt ! Time loop
        
        call all%UpdateMagneticField() ! Update magnetic field
        call all%UpdateElectricField() ! Update electric field

        call all%GaussianSource(t) ! Add a Gaussian source

        if (mod(t, 20) == 0) then ! Write to file every 20 time steps
            call all%WriteToText()
        end if

    end do 

    call all%EverythingFini() ! Free memory

    print*, 'Done!' ! Print done

end program