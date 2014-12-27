module utils
    implicit none

    contains
        subroutine stop_error(msg) !{{{
        ! Aborts the program with nonzero exit code
        !
        ! The statement "stop msg" will return 0 exit code when compiled using
        ! gfortran. stop_error() uses the statement "stop 1" which returns an exit code
        ! 1 and a print statement to print the message.
        !
        ! Example
        ! -------
        !
        ! call stop_error("Invalid argument")

            character(len=*) :: msg ! Message to print on stdout
            print *, msg
            stop 1
        end subroutine stop_error !}}}

end module utils
