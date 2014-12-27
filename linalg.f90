module linalg
    use types
    use matrix_routines
    implicit none

    public :: qdr
    private

    contains
        subroutine qdr(a, q, d, r, info) !{{{
            real(dp), dimension(:,:), intent(in) :: a
            real(dp), dimension(:,:), allocatable, intent(out) :: q,d,r
            integer, intent(out) :: info

            real(dp), dimension(:,:), allocatable :: tr
            real(dp), dimension(:), allocatable :: vd,vdi,tau, work, work_q

            integer, dimension(:), allocatable :: sh

            integer :: lda, lwork, k, m, n, mn, ws
            integer :: i,s

            q = a

            sh = shape(a)
            m = sh(1)
            n = sh(2)

            lda = m

            mn = min(m,n)
            allocate(tau(mn))
            allocate(work_q(1))
            ! Compute the optimal work array: query
            call dgeqrf(m, n, q, lda, tau, work_q, -1, info)
            ws = int(abs(work_q(1)))
            allocate(work(ws))

            ! Compute the QR factorization
            call dgeqrf(m, n, q, lda, tau, work, ws, info)

            ! R = D R'
            tr = triangular_upper(q)

            vd = diagonal(tr)
            s = size(vd)
            allocate(vdi(s))
            do i=1,s
                vdi(i) = 1/vd(i)
            end do

            r = tr
            call dgemm('n','n',mn,n,mn,1.0_dp,diagonal(vdi),mn,tr,mn,0.0_dp,r,mn)

            d = diagonal(vd)

            ! Compute the orthogonal matrix from the reflectors resulting from dgeqrf
            k = size(tau)
            call dorgqr(m, n, k, q, lda, tau, work_q, -1, info)
            if (work_q(1) .ne. ws) then
                ws = work_q(1)
                deallocate(work)
                allocate(work(ws))
            end if

            call dorgqr(m, n, k, q, lda, tau, work, ws, info)
        end subroutine qdr !}}}

        subroutine rdq(a, r, d, q, info) !{{{
            real(dp), dimension(:,:), intent(in) :: a
            real(dp), dimension(:,:), allocatable, intent(out) :: r, d, q
            integer, intent(out) :: info

            real(dp), dimension(:,:), allocatable :: tr
            real(dp), dimension(:), allocatable :: vd,vdi,tau, work, work_q

            integer, dimension(:), allocatable :: sh

            integer :: lda, lwork, k, m, n, mn, ws
            integer :: i,s

            q = a

            sh = shape(a)
            m = sh(1)
            n = sh(2)

            lda = m

            mn = min(m,n)
            allocate(tau(mn))
            allocate(work_q(1))
            ! Compute the optimal work array: query
            call dgerqf(m, n, q, lda, tau, work_q, -1, info)
            ws = int(abs(work_q(1)))
            allocate(work(ws))

            ! Compute the QR factorization
            call dgeqrf(m, n, q, lda, tau, work, ws, info)

            ! R = D R'
            tr = triangular_upper(q)

            vd = diagonal(tr)
            s = size(vd)
            allocate(vdi(s))
            do i=1,s
                vdi(i) = 1/vd(i)
            end do

            r = tr
            call dgemm('n','n',mn,n,mn,1.0_dp,diagonal(vdi),mn,tr,mn,0.0_dp,r,mn)

            d = diagonal(vd)

            ! Compute the orthogonal matrix from the reflectors resulting from dgeqrf
            k = size(tau)
            call dorgqr(m, n, k, q, lda, tau, work_q, -1, info)
            if (work_q(1) .ne. ws) then
                ws = work_q(1)
                deallocate(work)
                allocate(work(ws))
            end if

            call dorgqr(m, n, k, q, lda, tau, work, ws, info)
        end subroutine rdq !}}}

end module linalg
