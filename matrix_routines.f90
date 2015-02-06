module matrix_routines
    use types
    use utils
    implicit none

    private

    public :: diagonal, eye, kronecker
    public :: triangular_lower, triangular_upper
    public :: gemm_n, gemm_order

    interface diagonal !{{{
        module procedure i_diagonal_create
        module procedure s_diagonal_create
        module procedure d_diagonal_create
        module procedure q_diagonal_create
        module procedure i_diagonal_extract
        module procedure s_diagonal_extract
        module procedure d_diagonal_extract
        module procedure q_diagonal_extract
    end interface diagonal !}}}

    interface kronecker !{{{
        module procedure i_kronecker
        module procedure s_kronecker
        module procedure d_kronecker
        module procedure q_kronecker
    end interface kronecker !}}}

    interface triangular_upper !{{{
        module procedure i_triangular_upper
        module procedure s_triangular_upper
        module procedure d_triangular_upper
        module procedure q_triangular_upper
    end interface triangular_upper !}}}

    interface triangular_lower !{{{
        module procedure i_triangular_lower
        module procedure s_triangular_lower
        module procedure d_triangular_lower
        module procedure q_triangular_lower
    end interface triangular_lower !}}}

    interface gemm_n !{{{
        module procedure sgemm_n
        module procedure dgemm_n
    end interface gemm_n !}}}

    interface gemm_order !{{{
        module procedure sgemm_order
        module procedure dgemm_order
    end interface gemm_order !}}}

    contains

        ! diagonal {{{
        ! create diagonal {{{
        pure function i_diagonal_create(d) result(m) !{{{
            integer, dimension(:), intent(in) :: d

            integer, dimension(:,:), allocatable :: m
            integer :: i,n

            n = size(d)
            allocate(m(n,n))

            m = 0

            do i=1,n
                m(i,i) = d(i)
            end do
        end function i_diagonal_create !}}}

        pure function s_diagonal_create(d) result(m) !{{{
            real(sp), dimension(:), intent(in) :: d

            real(sp), dimension(:,:), allocatable :: m
            integer :: i,n

            n = size(d)
            allocate(m(n,n))

            m = 0.0

            do i=1,n
                m(i,i) = d(i)
            end do
        end function s_diagonal_create !}}}

        pure function d_diagonal_create(d) result(m) !{{{
            real(dp), dimension(:), intent(in) :: d

            real(dp), dimension(:,:), allocatable :: m
            integer :: i,n

            n = size(d)
            allocate(m(n,n))

            m = 0.0

            do i=1,n
                m(i,i) = d(i)
            end do
        end function d_diagonal_create !}}}

        pure function q_diagonal_create(d) result(m) !{{{
            real(qp), dimension(:), intent(in) :: d

            real(qp), dimension(:,:), allocatable :: m
            integer :: i,n

            n = size(d)
            allocate(m(n,n))

            m = 0.0

            do i=1,n
                m(i,i) = d(i)
            end do
        end function q_diagonal_create !}}}
        !}}}

        ! extract diagonal {{{
        pure function i_diagonal_extract(m,offset) result(d) !{{{
            integer, dimension(:,:), intent(in) :: m
            integer, intent(in), optional :: offset

            integer, dimension(:), allocatable :: d, sh
            integer :: i,l,n,os,u

            l = 1
            os = 0
            if (present(offset)) then
                os = offset
            end if


            sh = shape(m)
            n = min(sh(1),sh(2))
            !allocate(d(n))

            u = n

            if (abs(os) < n) then
                if (os > 0) then
                    u = n - os
                else if (offset < 0) then
                    l = 1 - os
                end if
            end if

            n = u-(l-1)
            allocate(d(n))
            d = 0

            do i=1, u
                d(i) = m(i+os,i+os)
            end do


        end function i_diagonal_extract !}}}

        pure function s_diagonal_extract(m) result(d) !{{{
            real(sp), dimension(:,:), intent(in) :: m

            real(sp), dimension(:), allocatable :: d
            integer, dimension(:), allocatable :: sh
            integer :: i,n

            sh = shape(m)
            n = min(sh(1),sh(2))
            allocate(d(n))

            d = 0.0

            do i=1,n
                d(i) = m(i,i)
            end do
        end function s_diagonal_extract !}}}

        pure function d_diagonal_extract(m) result(d) !{{{
            real(dp), dimension(:,:), intent(in) :: m

            real(dp), dimension(:), allocatable :: d
            integer, dimension(:), allocatable :: sh
            integer :: i,n

            sh = shape(m)
            n = min(sh(1),sh(2))
            allocate(d(n))

            d = 0.0

            do i=1,n
                d(i) = m(i,i)
            end do
        end function d_diagonal_extract !}}}

        pure function q_diagonal_extract(m) result(d) !{{{
            real(qp), dimension(:,:), intent(in) :: m

            real(qp), dimension(:), allocatable :: d
            integer, dimension(:), allocatable :: sh
            integer :: i,n

            sh = shape(m)
            n = min(sh(1),sh(2))
            allocate(d(n))

            d = 0.0

            do i=1,n
                d(i) = m(i,i)
            end do
        end function q_diagonal_extract !}}}
        !}}}
        !}}}

        pure function eye(n,offset) result(um) !{{{
            integer, intent(in) :: n
            integer, intent(in), optional :: offset

            integer, dimension(n,n) :: um

            integer :: i, l, u, os

            um = 0

            l = 1
            u = n
            os = 0

            if (present(offset)) then
                os = offset
            end if

            if (abs(os) < n) then
                if (os > 0) then
                    u = n - os
                else if (offset < 0) then
                    l = 1 - os
                end if

                do i=l, u
                    um(i, i+os) = 1
                end do
            end if

        end function eye !}}}

        ! kronecker {{{
        pure function i_kronecker(a,b) result(k) !{{{
            integer, dimension(:,:), intent(in) :: a, b

            integer, dimension(:), allocatable :: shape_a, shape_b
            integer :: nx, ny
            integer :: i,j

            integer :: x_low, x_upp, y_low, y_upp

            integer, dimension(:,:), allocatable :: k

            shape_a = shape(a)
            shape_b = shape(b)

            nx = shape_a(1) * shape_b(1)
            ny = shape_a(2) * shape_b(2)

            allocate(k(nx,ny))

            do j=1,shape_a(2)
                do i=1,shape_a(1)
                    x_low = (i-1)*shape_b(1) + 1
                    x_upp = i*shape_b(1)

                    y_low = (j-1)*shape_b(2) + 1
                    y_upp = j*shape_b(2)

                    k(x_low:x_upp,y_low:y_upp) = a(i,j) * b
                end do
            end do
        end function i_kronecker !}}}

        pure function s_kronecker(a,b) result(k) !{{{
            real(sp), dimension(:,:), intent(in) :: a, b

            integer, dimension(:), allocatable :: shape_a, shape_b
            integer :: nx, ny
            integer :: i,j

            integer :: x_low, x_upp, y_low, y_upp

            real(sp), dimension(:,:), allocatable :: k

            shape_a = shape(a)
            shape_b = shape(b)

            nx = shape_a(1) * shape_b(1)
            ny = shape_a(2) * shape_b(2)

            allocate(k(nx,ny))

            do j=1,shape_a(2)
                do i=1,shape_a(1)
                    x_low = (i-1)*shape_b(1) + 1
                    x_upp = i*shape_b(1)

                    y_low = (j-1)*shape_b(2) + 1
                    y_upp = j*shape_b(2)

                    k(x_low:x_upp,y_low:y_upp) = a(i,j) * b
                end do
            end do
        end function s_kronecker !}}}

        pure function d_kronecker(a,b) result(k) !{{{
            real(dp), dimension(:,:), intent(in) :: a, b

            integer, dimension(:), allocatable :: shape_a, shape_b
            integer :: nx, ny
            integer :: i,j

            integer :: x_low, x_upp, y_low, y_upp

            real(dp), dimension(:,:), allocatable :: k

            shape_a = shape(a)
            shape_b = shape(b)

            nx = shape_a(1) * shape_b(1)
            ny = shape_a(2) * shape_b(2)

            allocate(k(nx,ny))

            do j=1,shape_a(2)
                do i=1,shape_a(1)
                    x_low = (i-1)*shape_b(1) + 1
                    x_upp = i*shape_b(1)

                    y_low = (j-1)*shape_b(2) + 1
                    y_upp = j*shape_b(2)

                    k(x_low:x_upp,y_low:y_upp) = a(i,j) * b
                end do
            end do
        end function d_kronecker !}}}

        pure function q_kronecker(a,b) result(k) !{{{
            real(qp), dimension(:,:), intent(in) :: a, b

            integer, dimension(:), allocatable :: shape_a, shape_b
            integer :: nx, ny
            integer :: i,j

            integer :: x_low, x_upp, y_low, y_upp

            real(qp), dimension(:,:), allocatable :: k

            shape_a = shape(a)
            shape_b = shape(b)

            nx = shape_a(1) * shape_b(1)
            ny = shape_a(2) * shape_b(2)

            allocate(k(nx,ny))

            do j=1,shape_a(2)
                do i=1,shape_a(1)
                    x_low = (i-1)*shape_b(1) + 1
                    x_upp = i*shape_b(1)

                    y_low = (j-1)*shape_b(2) + 1
                    y_upp = j*shape_b(2)

                    k(x_low:x_upp,y_low:y_upp) = a(i,j) * b
                end do
            end do
        end function q_kronecker !}}}
        !}}}

        ! triangular_upper {{{
        pure function i_triangular_upper(m) result(r) !{{{
            integer, dimension(:,:), intent(in) :: m
            integer, dimension(:,:), allocatable :: r

            integer, dimension(:), allocatable :: sh
            integer :: i, j, lim

            sh = shape(m)
            ! This is currently broken in gfortran 4.9
            !allocate(r,source=m)
            ! Workaround:
            ! The statement r=m requires -standard-semantics in ifort as of 15.*
            !allocate(r(sh(1), sh(2)))
            r = m

            lim = min(sh(1),sh(2))

            do j = 1, lim
                do i = j+1,lim
                    r(i,j) = 0.0
                end do
            end do
        end function i_triangular_upper !}}}

        pure function s_triangular_upper(m) result(r) !{{{
            real(sp), dimension(:,:), intent(in) :: m
            real(sp), dimension(:,:), allocatable :: r

            integer, dimension(:), allocatable :: sh
            integer :: i, j, lim

            sh = shape(m)
            ! This is currently broken in gfortran 4.9
            !allocate(r,source=m)
            ! Workaround:
            ! The statement r=m requires -standard-semantics in ifort as of 15.*
            !allocate(r(sh(1), sh(2)))
            r = m

            lim = min(sh(1),sh(2))

            do j = 1, lim
                do i = j+1, lim
                    r(i,j) = 0.0
                end do
            end do
        end function s_triangular_upper !}}}

        pure function d_triangular_upper(m) result(r) !{{{
            real(dp), dimension(:,:), intent(in) :: m
            real(dp), dimension(:,:), allocatable :: r

            integer, dimension(:), allocatable :: sh
            integer :: i, j, lim

            sh = shape(m)
            ! This is currently broken in gfortran 4.9
            !allocate(r,source=m)
            ! Workaround:
            ! The statement r=m requires -standard-semantics in ifort as of 15.*
            !allocate(r(sh(1), sh(2)))
            r = m

            lim = min(sh(1),sh(2))

            do j = 1, lim
                do i = j+1, lim
                    r(i,j) = 0.0
                end do
            end do
        end function d_triangular_upper !}}}

        pure function q_triangular_upper(m) result(r) !{{{
            real(qp), dimension(:,:), intent(in) :: m
            real(qp), dimension(:,:), allocatable :: r

            integer, dimension(:), allocatable :: sh
            integer :: i, j, lim

            sh = shape(m)
            ! This is currently broken in gfortran 4.9
            !allocate(r,source=m)
            ! Workaround:
            ! The statement r=m requires -standard-semantics in ifort as of 15.*
            !allocate(r(sh(1), sh(2)))
            r = m

            lim = min(sh(1),sh(2))

            do j = 1, lim
                do i = j+1, lim
                    r(i,j) = 0.0
                end do
            end do
        end function q_triangular_upper !}}}
        !}}}

        ! triangular_lower {{{
        pure function i_triangular_lower(m) result(l) !{{{
            integer, dimension(:,:), intent(in) :: m
            integer, dimension(:,:), allocatable :: l

            integer, dimension(:), allocatable :: sh
            integer :: i, j, lim

            sh = shape(m)
            ! This is currently broken in gfortran 4.9
            !allocate(l,source=m)
            ! Workaround:
            ! The statement l=m requires -standard-semantics in ifort as of 15.*
            !allocate(l(sh(1), sh(2)))
            l = m

            lim = min(sh(1),sh(2))

            do j = 1, lim
                do i = 1, j-1
                    l(i,j) = 0.0
                end do
            end do
        end function i_triangular_lower !}}}

        pure function s_triangular_lower(m) result(l) !{{{
            real(sp), dimension(:,:), intent(in) :: m
            real(sp), dimension(:,:), allocatable :: l

            integer, dimension(:), allocatable :: sh
            integer :: i, j, lim

            sh = shape(m)
            ! This is currently broken in gfortran 4.9
            !allocate(l,source=m)
            ! Workaround:
            ! The statement l=m requires -standard-semantics in ifort as of 15.*
            !allocate(l(sh(1), sh(2)))
            l = m

            lim = min(sh(1),sh(2))

            do j = 1, lim
                do i = 1, j-1
                    l(i,j) = 0.0
                end do
            end do
        end function s_triangular_lower !}}}

        pure function d_triangular_lower(m) result(l) !{{{
            real(dp), dimension(:,:), intent(in) :: m
            real(dp), dimension(:,:), allocatable :: l

            integer, dimension(:), allocatable :: sh
            integer :: i, j, lim

            sh = shape(m)
            ! This is currently broken in gfortran 4.9
            !allocate(l,source=m)
            ! Workaround:
            ! The statement l=m requires -standard-semantics in ifort as of 15.*
            !allocate(l(sh(1), sh(2)))
            l = m

            lim = min(sh(1),sh(2))

            do j = 1, lim
                do i = 1, j-1
                    l(i,j) = 0.0
                end do
            end do
        end function d_triangular_lower !}}}

        pure function q_triangular_lower(m) result(l) !{{{
            real(qp), dimension(:,:), intent(in) :: m
            real(qp), dimension(:,:), allocatable :: l

            integer, dimension(:), allocatable :: sh
            integer :: i, j, lim

            sh = shape(m)
            ! This is currently broken in gfortran 4.9
            !allocate(l,source=m)
            ! Workaround:
            ! The statement l=m requires -standard-semantics in ifort as of 15.*
            !allocate(l(sh(1), sh(2)))
            l = m

            lim = min(sh(1),sh(2))

            do j = 1, lim
                do i = 1, j-1
                    l(i,j) = 0.0
                end do
            end do
        end function q_triangular_lower !}}}
        !}}}

        ! gemm_n {{{
    subroutine sgemm_n(m,ms) !{{{
        real(sp), dimension(:,:), intent(inout) :: m
        real(sp), dimension(:,:,:), intent(in) :: ms

        real(sp), dimension(:,:), allocatable :: work

        integer, dimension(2) :: sh_m
        integer, dimension(3) :: sh_ms
        integer :: i

        integer :: rows, columns, n_arrays

        sh_m = shape(m)
        rows = sh_m(1)
        columns = sh_m(2)

        sh_ms = shape(ms)

        if (rows /= columns) then
            call stop_error("sgemm_n: non-quadratic matrices passed.")
        else if (rows /= sh_ms(1)) then
            call stop_error("sgemm_n: row dimension of matrix does not agree with matrices")
        else if (columns /= sh_ms(2)) then
            call stop_error("sgemm_n: columns dimension of matrix does not agree with matrices")
        else
            n_arrays = sh_ms(3)

            allocate(work(rows,columns))

            do i=1, n_arrays ! call to gemm, such that work = work * ms(1) * ms(2) ... ms(n)
                call sgemm ('n', 'n', rows, columns, rows, 1.0, m, rows, ms(:,:,i), rows, 0.0, work, rows)
                m = work
            end do
        end if
    end subroutine sgemm_n !}}}

    subroutine dgemm_n(m,ms) !{{{
        real(dp), dimension(:,:), intent(inout) :: m
        real(dp), dimension(:,:,:), intent(in) :: ms

        real(dp), dimension(:,:), allocatable :: work

        integer, dimension(2) :: sh_m
        integer, dimension(3) :: sh_ms
        integer :: i

        integer :: rows, columns, n_arrays

        sh_m = shape(m)
        rows = sh_m(1)
        columns = sh_m(2)

        sh_ms = shape(ms)

        if (rows /= columns) then
            call stop_error("dgemm_n: non-quadratic matrices passed.")
        else if (rows /= sh_ms(1)) then
            call stop_error("dgemm_n: row dimension of matrix does not agree with matrices")
        else if (columns /= sh_ms(2)) then
            call stop_error("dgemm_n: columns dimension of matrix does not agree with matrices")
        else
            n_arrays = sh_ms(3)

            allocate(work(rows,columns))

            do i=1, n_arrays
                ! call to gemm, such that work = work * ms(1) * ms(2) ... ms(n)
                call dgemm ('n', 'n', rows, columns, rows, 1.0, m, rows, ms(:,:,i), rows, 0.0, work, rows)
                m = work
            end do
        end if
    end subroutine dgemm_n !}}}
    !}}} 

    ! gemm_order {{{
    subroutine sgemm_order(m,ms,order) !{{{
        real(sp), dimension(:,:), intent(inout) :: m
        real(sp), dimension(:,:,:), intent(in) :: ms

        integer, dimension(:), intent(in) :: order

        real(sp), dimension(:,:), allocatable :: work

        integer, dimension(2) :: sh_m
        integer, dimension(3) :: sh_ms
        integer :: i,o

        integer :: rows, columns

        o = size(order)

        sh_m = shape(m)
        rows = sh_m(1)
        columns = sh_m(2)

        sh_ms = shape(ms)

        if (rows /= columns) then
            call stop_error("dgemm_n: non-quadratic matrices passed.")
        else if (rows /= sh_ms(1)) then
            call stop_error("dgemm_n: row dimension of matrix does not agree with matrices")
        else if (columns /= sh_ms(2)) then
            call stop_error("dgemm_n: columns dimension of matrix does not agree with matrices")
        else
            allocate(work(rows,columns))

            do i=1, o
                ! call to gemm, such that work = work * ms(1) * ms(2) ... ms(n)
                call dgemm ('n', 'n', rows, columns, rows, 1.0, m, rows, ms(:,:,i), rows, 0.0, work, rows)
                m = work
            end do
        end if
    end subroutine sgemm_order !}}}

    subroutine dgemm_order(m,ms,order) !{{{
        real(dp), dimension(:,:), intent(inout) :: m
        real(dp), dimension(:,:,:), intent(in) :: ms

        integer, dimension(:), intent(in) :: order

        real(dp), dimension(:,:), allocatable :: work

        integer, dimension(2) :: sh_m
        integer, dimension(3) :: sh_ms
        integer :: i,o

        integer :: rows, columns

        o = size(order)

        sh_m = shape(m)
        rows = sh_m(1)
        columns = sh_m(2)

        sh_ms = shape(ms)

        if (rows /= columns) then
            call stop_error("dgemm_n: non-quadratic matrices passed.")
        else if (rows /= sh_ms(1)) then
            call stop_error("dgemm_n: row dimension of matrix does not agree with matrices")
        else if (columns /= sh_ms(2)) then
            call stop_error("dgemm_n: columns dimension of matrix does not agree with matrices")
        else
            allocate(work(rows,columns))

            do i=1, o
                ! call to gemm, such that work = work * ms(1) * ms(2) ... ms(n)
                call dgemm ('n', 'n', rows, columns, rows, 1.0, m, rows, ms(:,:,i), rows, 0.0, work, rows)
                m = work
            end do
        end if
    end subroutine dgemm_order !}}}
    ! }}}
end module matrix_routines
