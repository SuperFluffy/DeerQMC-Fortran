module kinetic
    use, intrinsic :: iso_fortran_env
    use matrix_routines
    implicit none
    private

    public :: hopping_matrix, next_nearest

    interface next_nearest
        module procedure next_nearest_1d
        module procedure next_nearest_2d
    end interface next_nearest

    contains
        pure function next_nearest_1d(nx, phase) result (k) !{{{
            integer, intent(in) :: nx
            real, intent(in) :: phase

            real(real64), dimension(nx,nx) :: k

            k = eye(nx,+1) + eye(nx,-1)

            if (nx > 1) then
                k(1,nx) = phase
                k(nx,1) = phase
            end if
        end function next_nearest_1d !}}}

        pure function next_nearest_2d(nx,ny, xphase, yphase) result (k) !{{{
            integer, intent(in) :: nx, ny
            real, intent(in) :: xphase, yphase

            real(real64), dimension(:,:), allocatable :: k, k_x, k_y, um_x, um_y

            um_x = eye(nx)
            um_y = eye(ny)

            k_x = next_nearest_1d(nx,xphase)
            k_y = next_nearest_1d(ny,yphase)

            k = kronecker_product(um_y,k_x) + kronecker_product(k_y,um_x)
        end function next_nearest_2d !}}}

        pure function hopping_matrix(t,nx,ny,xphase,yphase) result(k) !{{{
            real, intent(in) :: t
            integer, intent(in) :: nx, ny
            real, intent(in) :: xphase, yphase

            real(real64), dimension(:,:), allocatable :: k

            k = -t * next_nearest_2d(nx,ny,xphase,yphase)

        end function hopping_matrix !}}}

end module kinetic
