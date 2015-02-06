module linalg
    use types
    use matrix_routines
    use utils
    implicit none

    private
    public :: qdr, rdq
    public :: stable_mm_qdr, stable_mm_rdq
    
    contains
        subroutine qdr(a, q, d, r, info) !{{{
            real(dp), dimension(:,:), intent(in) :: a
            real(dp), dimension(:,:), allocatable, intent(out) :: q,d,r
            integer, intent(out) :: info

            real(dp), dimension(:,:), allocatable :: tr
            real(dp), dimension(:), allocatable :: vd,vdi,tau, work, work_q

            integer, dimension(:), allocatable :: sh

            integer :: lwork, m, ws
            integer :: i,s

            q = a

            sh = shape(a)
            m = sh(1)

            if (m /= sh(2)) then
                call stop_error("qdr: non-square matrix detected.")
            end if

            allocate(tau(m))
            allocate(work_q(1))

            ! Compute the optimal work array: query
            call dgeqrf(m, m, q, m, tau, work_q, -1, info)
            ws = int(abs(work_q(1)))
            allocate(work(ws))

            ! Compute the QR factorization
            call dgeqrf(m, m, q, m, tau, work, ws, info)

            ! R = D R'
            tr = triangular_upper(q)

            vd = diagonal(tr)
            s = size(vd)
            allocate(vdi(s))
            do i=1,s
                vdi(i) = 1/vd(i)
            end do

            r = tr
            call dgemm('n','n',m,m,m,1.0_dp,diagonal(vdi),m,tr,m,0.0_dp,r,m)

            d = diagonal(vd)

            ! Compute the orthogonal matrix from the reflectors resulting from dgeqrf
            call dorgqr(m, m, m, q, m, tau, work_q, -1, info)

            ! Reallocate the work array if the optimal size is different from the QR factorization above.
            if (work_q(1) .ne. ws) then
                ws = work_q(1)
                deallocate(work)
                allocate(work(ws))
            end if

            call dorgqr(m, m, m, q, m, tau, work, ws, info)
        end subroutine qdr !}}}

        subroutine rdq(a, r, d, q, info) !{{{
            ! RDQ decomposition of the matrix A. Assumes that A is quadratic!
            real(dp), dimension(:,:), intent(in) :: a
            real(dp), dimension(:,:), allocatable, intent(out) :: r, d, q
            integer, intent(out) :: info

            real(dp), dimension(:,:), allocatable :: tr
            real(dp), dimension(:), allocatable :: vd,vdi,tau, work, work_q

            integer, dimension(:), allocatable :: sh

            integer :: lwork, m, n, mn, ws
            integer :: i,s

            q = a

            sh = shape(a)
            m = sh(1)
            n = sh(2)
            if (m /= n) then
                call stop_error("rdq: non-square matrix detected.")
            end if

            allocate(tau(m))
            allocate(work_q(1))
            ! Compute the optimal work array: query
            call dgerqf(m, m, q, m, tau, work_q, -1, info)
            ws = int(abs(work_q(1)))
            allocate(work(ws))

            ! Compute the RQ factorization
            call dgerqf(m, m, q, m, tau, work, ws, info)

            tr = triangular_upper(q)

            vd = diagonal(tr)
            s = size(vd)
            allocate(vdi(s))
            do i=1,s
                vdi(i) = 1/vd(i)
            end do

            r = tr
            call dgemm('n','n',m,m,m,1.0_dp,tr,m,diagonal(vdi),m,0.0_dp,r,m)

            d = diagonal(vd)

            ! Compute the orthogonal matrix from the reflectors resulting from dgerqf
            call dorgrq(m, m, m, q, m, tau, work_q, -1, info)

            ! Reallocate the work array if the optimal size is different from the RQ factorization above.
            if (work_q(1) .ne. ws) then
                ws = work_q(1)
                deallocate(work)
                allocate(work(ws))
            end if

            call dorgrq(m, m, m, q, m, tau, work, ws, info)
        end subroutine rdq !}}}

end module linalg
