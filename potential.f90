module potential
    use matrix_routines

    implicit none

    private

    contains

        pure function potential_matrix(n,b,lambda,mu) result(v)
            real, intent(in) :: b, lambda, mu
            real, dimension(:,:), allocatable :: c, m, v

        end function potential_matrix

        pure function interaction_matrix(lambda, spinfield) result(p) !{{{
            integer, dimension(:), intent(in) :: spinfield
            integer :: i,n
            real, dimension(:,:), allocatable :: p

            n = size(spinfield)

            allocate(p(n,n))

            do i=1,n
                p(i,i) = spinfield(i)
            end do

            p = lambda * p
        end function interaction_matrix !}}}

        pure function construct_spin_field(n) result(s) !{{{
            ! This subroutine calls the Fortran intrinsic `random_number`
            ! intrinsic, assuming that the seed is externally initialized.
            integer, intent(in) :: n

            integer, dimension(n) :: s

            integer :: i
            real :: r

            do i=1,n
                call random_number(r)
                if (r < 0.5) then
                    s(i) = 0
                else
                    s(i) = 1
                end if
            end do
        end subroutine construct_spin_field !}}}
end module potential
