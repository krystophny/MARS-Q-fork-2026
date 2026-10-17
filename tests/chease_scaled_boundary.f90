! Diagnostic fixture only; no real equilibrium, map, or scientific admission.
! Native routines: MARS-Q a4245da4996e5582024278ad0b68eeaa9a13c7cc.
! folded_band: public CHEASE b942066c64a6dfe1baaba546945473fb6db97828.
!
! The SPD matrix is independently prescribed diagonal storage. Every physical
! equation has diagonal 1e-12; inactive LIMITA equations have diagonal one and
! zero right-hand side. Native up/down numbering leaves row N physical and
! row N-1 inactive. Interior LDLT pivots equal their own original diagonal,
! while the final native comparison uses the preceding inactive unit pivot.
! Folded removal/permutation makes both the bulk and center matrices 1e-12 I.
!
! Build order: prec_const_fixture.f90, native_factor_a4245da.f,
! folded_band_b942066c.f90, this file. Retain the native default-real-8 and
! default-double-8 flags; use bounds checks. No BLAS library is required.
program scaled_boundary_oracle
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use folded_band, only: folded_factor, folded_solve, folded_reset, folded_active
    implicit none

    integer, parameter :: ns = 3, nt = 8
    integer, parameter :: n = 4 * (ns + 1) * nt, m = 4 * nt + 12
    integer, parameter :: mp = 4 * 640 + 12 ! Owning a4245 physical stride.
    integer, parameter :: centers(3) = [4 * nt - 2, 4 * nt - 1, 4 * nt]
    real(dp), parameter :: physical_diagonal = 1.0e-12_dp
    real(dp), parameter :: native_eps = 1.0e-10_dp
    real(dp), parameter :: tolerance = 128.0_dp * epsilon(1.0_dp)

    integer :: node_map((ns + 1) * nt), seen(nt)
    integer :: radial, angle, old, rank, i, info
    logical :: active(n)
    real(dp), allocatable :: original(:, :), factors(:, :)
    real(dp) :: diagonal(n), truth(n), rhs(n), solution(n)
    real(dp) :: solution_error, scaled_residual, relative_residual

    interface
        subroutine aldlt(a, eps, rows, band_width, stride, status)
            import dp
            integer, intent(in) :: rows, band_width, stride
            real(dp), intent(in) :: eps
            real(dp), intent(inout) :: a(*)
            integer, intent(out) :: status
        end subroutine aldlt
        subroutine lyv(a, u, rows, vector_size, band_width, stride)
            import dp
            integer, intent(in) :: rows, vector_size, band_width, stride
            real(dp), intent(in) :: a(*)
            real(dp), intent(inout) :: u(vector_size)
        end subroutine lyv
        subroutine dwy(a, u, rows, vector_size, band_width, stride)
            import dp
            integer, intent(in) :: rows, vector_size, band_width, stride
            real(dp), intent(in) :: a(*)
            real(dp), intent(inout) :: u(vector_size)
        end subroutine dwy
        subroutine ltxw(a, u, rows, vector_size, band_width, stride)
            import dp
            integer, intent(in) :: rows, vector_size, band_width, stride
            real(dp), intent(in) :: a(*)
            real(dp), intent(inout) :: u(vector_size)
        end subroutine ltxw
    end interface

    allocate(original(mp, n), factors(mp, n))

    ! Own NUPDWN equations, independently applied to every native radial ring.
    ! See a4245 CheaseMerge/chease.f:3988-3999; no folded-order indices used.
    do radial = 0, ns
        seen = 0
        do angle = 1, nt
            old = radial * nt + angle
            if (angle >= 2 .and. angle <= nt / 2 + 1) then
                node_map(old) = old + angle - 2
            else if (angle >= nt / 2 + 2) then
                node_map(old) = (radial + 2) * nt - 2 * angle + 3
            else
                node_map(old) = old
            end if
            rank = node_map(old) - radial * nt
            if (rank < 1 .or. rank > nt) error stop 'Node rank outside its ring'
            seen(rank) = seen(rank) + 1
        end do
        if (any(seen /= 1)) error stop 'Native node map is not a permutation'
    end do

    ! LIMITA layout: three center equations, all four interior components,
    ! and only components two/four at the outer boundary remain physical.
    active = .false.
    active(centers) = .true.
    do radial = 1, ns
        do angle = 1, nt
            old = 4 * (node_map(radial * nt + angle) - 1)
            if (radial < ns) then
                active(old + 1:old + 4) = .true.
            else
                active(old + 2) = .true.
                active(old + 4) = .true.
            end if
        end do
    end do
    if (count(active) /= (4 * ns - 2) * nt + 3) then
        error stop 'Unexpected number of physical equations'
    end if
    if (.not. active(n)) error stop 'Final native row must be physical'
    if (active(n - 1)) error stop 'Penultimate native row must be inactive'

    diagonal = merge(physical_diagonal, 1.0_dp, active)
    if (any(diagonal <= 0.0_dp)) error stop 'Prescribed matrix is not SPD'
    original = 0.0_dp
    original(1, :) = diagonal
    truth = 0.0_dp
    do i = 1, n
        if (active(i)) truth(i) = real(mod(i, 11) - 5, dp) / 8.0_dp
    end do
    rhs = diagonal * truth ! Independent diagonal matvec; inactive RHS is zero.

    ! No pivoting, updates or cancellations can change this diagonal matrix.
    ! All interior checks see d_i against native_eps * d_i. Only the final
    ! check sees 1e-12 against native_eps * 1, and must reject at native EPS.
    factors = original
    call aldlt(factors(1, 1), native_eps, n, m, mp, info)
    if (info /= -1) error stop 'Native scale-sensitive rejection not observed'
    if (any(factors /= original)) then
        error stop 'Diagonal native fixture changed during factorization'
    end if
    write(*, '(a,i0,a,es13.5)') &
        'native expected rejection: info=', info, ', EPS=', native_eps
    write(*, '(a,2es13.5)') &
        'final and preceding diagonals: ', diagonal(n), diagonal(n - 1)

    ! The exact public module invokes the same extracted native routines.
    ! Removing inactive units leaves homogeneous 1e-12 I bulk and Schur blocks.
    call folded_factor(original, ns, nt, n, m, mp, node_map, native_eps, info)
    if (info /= 0) error stop 'Folded factorization failed at the same EPS'
    if (.not. folded_active) error stop 'Folded factor state was not activated'
    if (any(original /= factors)) error stop 'Folded factor changed native input'
    solution = rhs
    call folded_solve(solution, info)
    if (info /= 0) error stop 'Folded solve failed'
    if (.not. all(ieee_is_finite(solution))) error stop 'Nonfinite folded solution'
    if (any(solution /= 0.0_dp .and. .not. active)) then
        error stop 'Inactive solution equations are not zero'
    end if
    solution_error = maxval(abs(solution - truth))
    scaled_residual = maxval(abs(diagonal * solution - rhs) / diagonal)
    relative_residual = maxval(abs(diagonal * solution - rhs)) / maxval(abs(rhs))
    if (solution_error > tolerance) error stop 'Folded known solution mismatch'
    if (scaled_residual > tolerance) error stop 'Folded scaled residual mismatch'
    if (relative_residual > tolerance) error stop 'Folded relative residual mismatch'
    write(*, '(a,3es13.5)') &
        'folded same-EPS solution/scaled/relative residual: ', &
        solution_error, scaled_residual, relative_residual
    call folded_reset()

    ! Positive control for all four extracted native routines at unchanged EPS.
    ! Keep the same constrained layout and zero inactive RHS, with physical
    ! diagonals also one; the resulting independently known matrix is I.
    factors = 0.0_dp
    factors(1, :) = 1.0_dp
    solution = truth
    call aldlt(factors(1, 1), native_eps, n, m, mp, info)
    if (info /= 0) error stop 'Native unit-scale positive control failed'
    call lyv(factors(1, 1), solution, n, n, m, mp)
    call dwy(factors(1, 1), solution, n, n, m, mp)
    call ltxw(factors(1, 1), solution, n, n, m, mp)
    if (.not. all(ieee_is_finite(solution))) error stop 'Nonfinite native control'
    if (maxval(abs(solution - truth)) > tolerance) then
        error stop 'Native unit-scale known solution mismatch'
    end if
    write(*, '(a)') 'PASS native final-scale rejection and folded same-EPS solve'
end program scaled_boundary_oracle
