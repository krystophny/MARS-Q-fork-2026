! Copyright 2026 ITP plasma group
! SPDX-License-Identifier: Apache-2.0
!
! Ported from public CHEASE b942066c64a6dfe1baaba546945473fb6db97828.
! The owning precision is KIND(1D0); factor/solve routines stay native.
! Active-variable permutation only: native assembly and reconstruction remain
! unchanged. Eliminate the folded angular bulk before the three center DOFs.
MODULE folded_band
  USE, INTRINSIC :: iso_fortran_env, ONLY: int64
  IMPLICIT NONE
  INTEGER, PARAMETER :: RKIND = KIND(1D0)
  PRIVATE
  PUBLIC :: folded_factor, folded_solve, folded_reset, folded_active
  LOGICAL, SAVE :: folded_active = .FALSE.
  INTEGER, ALLOCATABLE, SAVE :: index_map(:)
  REAL(RKIND), ALLOCATABLE, SAVE :: bulk(:,:), coupling(:,:), response(:,:)
  REAL(RKIND), SAVE :: schur(3,3)
  INTEGER, SAVE :: nbulk=0, width=0
CONTAINS
  SUBROUTINE folded_reset()
    folded_active = .FALSE.
    nbulk = 0
    width = 0
    IF (ALLOCATED(index_map)) DEALLOCATE(index_map)
    IF (ALLOCATED(bulk)) DEALLOCATE(bulk)
    IF (ALLOCATED(coupling)) DEALLOCATE(coupling)
    IF (ALLOCATED(response)) DEALLOCATE(response)
  END SUBROUTINE

  SUBROUTINE folded_factor(a, ns, nt, n, m, mp, node_map, eps, info)
    INTEGER, INTENT(IN) :: ns, nt, n, m, mp, node_map(:)
    REAL(RKIND), INTENT(IN) :: a(*), eps
    INTEGER, INTENT(OUT) :: info
    INTEGER :: block_size, angle, radial, component, old, new, i, j
    INTEGER :: ii, jj, lo, hi, k, rank, status
    INTEGER(int64) :: offset
    REAL(RKIND) :: h(3,3), value, svalue
    CALL folded_reset()
    info = -2
    IF (ns < 3 .OR. nt < 4 .OR. MOD(nt,2) /= 0) RETURN
    IF (INT(n,int64) /= 4_int64*(INT(ns,int64)+1)*INT(nt,int64)) RETURN
    IF (mp < m) RETURN
    IF (SIZE(node_map) < (ns+1)*nt) RETURN
    block_size = 4*ns-2
    nbulk = block_size*nt
    width = 8*ns+4
    ALLOCATE(index_map(n), bulk(width,nbulk), coupling(nbulk,3), &
             response(nbulk,3), STAT=status)
    IF (status /= 0) THEN
      CALL folded_reset()
      RETURN
    ENDIF
    index_map = 0
    bulk = 0._RKIND
    coupling = 0._RKIND
    h = 0._RKIND
    DO angle=1,nt
      ! The native center-ring permutation already folds the periodic cycle.
      rank = node_map(angle)-1
      IF (rank < 0 .OR. rank >= nt) RETURN
      DO radial=1,ns
        IF (node_map(radial*nt+angle) < 1) RETURN
        IF (node_map(radial*nt+angle) > n/4) RETURN
        old = 4*(node_map(radial*nt+angle)-1)
        IF (old < 0 .OR. old+4 > n) RETURN
        IF (radial < ns) THEN
          DO component=1,4
            new = rank*block_size+4*(radial-1)+component
            index_map(old+component) = new
          END DO
        ELSE
          new = rank*block_size+4*(ns-1)
          index_map(old+2) = new+1
          index_map(old+4) = new+2
        ENDIF
      END DO
    END DO
    DO k=1,3
      index_map(4*nt-3+k) = -k
    END DO
    ! Require a bijection; a malformed node map must fail before factorization.
    BLOCK
      INTEGER :: seen(nbulk)
      seen = 0
      DO i=1,n
        IF (index_map(i) > 0) seen(index_map(i)) = seen(index_map(i))+1
      END DO
      IF (ANY(seen /= 1)) RETURN
    END BLOCK
    DO i=1,n
      offset = INT(i-1,int64)*INT(mp,int64)
      ii = index_map(i)
      IF (ii == 0 .AND. a(offset+1) /= 1._RKIND) RETURN
      DO j=i,MIN(n,i+m-1)
        value = a(offset+j-i+1)
        IF (value == 0._RKIND) CYCLE
        jj = index_map(j)
        IF (ii == 0 .OR. jj == 0) THEN
          ! LIMITA installs unit inactive equations with no active coupling.
          IF (i /= j .OR. value /= 1._RKIND) RETURN
        ELSE IF (ii > 0 .AND. jj > 0) THEN
          lo = MIN(ii,jj)
          hi = MAX(ii,jj)
          IF (hi-lo+1 > width) RETURN
          bulk(hi-lo+1,lo) = value
        ELSE IF (ii < 0 .AND. jj < 0) THEN
          h(-ii,-jj) = value
          h(-jj,-ii) = value
        ELSE IF (ii > 0) THEN
          coupling(ii,-jj) = value
        ELSE
          coupling(jj,-ii) = value
        ENDIF
      END DO
    END DO
    CALL ALDLT(bulk,eps,nbulk,width,width,info)
    IF (info /= 0) RETURN
    response = coupling
    DO k=1,3
      CALL solve_bulk(response(:,k))
    END DO
    schur = 0._RKIND
    DO i=1,3
      DO j=i,3
        svalue = h(i,j)-DOT_PRODUCT(coupling(:,i),response(:,j))
        schur(j-i+1,i) = svalue
      END DO
    END DO
    CALL ALDLT(schur,eps,3,3,3,info)
    IF (info /= 0) RETURN
    folded_active = .TRUE.
  END SUBROUTINE

  SUBROUTINE solve_bulk(x)
    REAL(RKIND), INTENT(INOUT) :: x(:)
    CALL LYV(bulk,x,nbulk,nbulk,width,width)
    CALL DWY(bulk,x,nbulk,nbulk,width,width)
    CALL LTXW(bulk,x,nbulk,nbulk,width,width)
  END SUBROUTINE

  SUBROUTINE folded_solve(b, info)
    REAL(RKIND), INTENT(INOUT) :: b(:)
    INTEGER, INTENT(OUT) :: info
    REAL(RKIND), ALLOCATABLE :: y(:)
    REAL(RKIND) :: center(3)
    INTEGER :: i, k, status
    info = -2
    IF (.NOT. folded_active) RETURN
    IF (SIZE(b) < SIZE(index_map)) RETURN
    ALLOCATE(y(nbulk), STAT=status)
    IF (status /= 0) RETURN
    y = 0._RKIND
    center = 0._RKIND
    DO i=1,SIZE(index_map)
      k = index_map(i)
      IF (k > 0) THEN
        y(k) = b(i)
      ELSE IF (k < 0) THEN
        center(-k) = b(i)
      ELSE IF (b(i) /= 0._RKIND) THEN
        ! LIMITB must have zeroed every inactive native equation.
        RETURN
      ENDIF
    END DO
    CALL solve_bulk(y)
    DO k=1,3
      center(k) = center(k)-DOT_PRODUCT(coupling(:,k),y)
    END DO
    CALL LYV(schur,center,3,3,3,3)
    CALL DWY(schur,center,3,3,3,3)
    CALL LTXW(schur,center,3,3,3,3)
    y = y-MATMUL(response,center)
    DO i=1,SIZE(index_map)
      k = index_map(i)
      IF (k > 0) THEN
        b(i) = y(k)
      ELSE IF (k < 0) THEN
        b(i) = center(-k)
      ELSE
        b(i) = 0._RKIND
      ENDIF
    END DO
    info = 0
  END SUBROUTINE
END MODULE folded_band
