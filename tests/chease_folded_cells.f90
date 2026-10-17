! Public CHEASE b942066c physical-clique oracle, using owning EPS1e-10.
program folded_oracle
  use prec_const
  use folded_band
  implicit none
  integer :: ns, nt, n, m, mp, r, t, neighbor, q, d, i, j, k, count, info, rank
  integer :: ids(16), centers(3)
  integer, allocatable :: nodes(:), active(:)
  real(rkind), allocatable :: dense(:,:), band(:,:), original(:,:), truth(:), rhs(:), old_rhs(:)
  real(rkind) :: value, inactive_rhs(1)
  call folded_reset()
  inactive_rhs=123._rkind
  call folded_solve(inactive_rhs,info)
  if (info /= -2 .or. inactive_rhs(1) /= 123._rkind) &
    error stop "Inactive solve changed RHS or did not reject"
  do ns=3,5,2
    do nt=8,20,4
      n=4*(ns+1)*nt
      m=4*nt+12
      mp=m+2
      allocate(nodes((ns+1)*nt),active(n),dense(n,n),band(mp,n),original(mp,n), &
               truth(n),rhs(n),old_rhs(n))
      do r=0,ns
        do t=0,nt-1
          if(t == 0) then
            rank=0
          else if(t <= nt/2) then
            rank=2*t-1
          else
            rank=2*(nt-t)
          endif
          nodes(r*nt+t+1)=r*nt+rank+1
        enddo
      enddo
      active=0
      centers=[4*nt-2,4*nt-1,4*nt]
      active(centers)=1
      do r=1,ns
        do t=1,nt
          i=4*(nodes(r*nt+t)-1)
          if(r < ns) then
            active(i+1:i+4)=1
          else
            active(i+2)=1
            active(i+4)=1
          endif
        enddo
      enddo
      dense=0
      ! Construct independent physical cell cliques, including the periodic
      ! seam and both surviving boundary derivatives. No new-order indices.
      do r=1,ns-1
        do t=1,nt
          neighbor=mod(t,nt)+1
          count=0
          do q=r,r+1
            do k=1,2
              j=t
              if(k == 2) j=neighbor
              i=4*(nodes(q*nt+j)-1)
              do d=1,4
                if(active(i+d) == 0) cycle
                count=count+1
                ids(count)=i+d
              enddo
            enddo
          enddo
          do i=1,count
            do j=i+1,count
              call couple(ids(i),ids(j))
            enddo
          enddo
        enddo
      enddo
      do t=1,nt
        do d=1,4
          i=4*(nodes(nt+t)-1)+d
          do k=1,3
            call couple(i,centers(k))
          enddo
        enddo
      enddo
      do i=1,3
        do j=i+1,3
          call couple(centers(i),centers(j))
        enddo
      enddo
      band=0
      band(m+1:mp,:)=77
      do i=1,n
        dense(i,i)=1+sum(abs(dense(i,:)))
        do j=i,n
          if(j-i+1 > m) then
            if(dense(i,j) /= 0) error stop 'Fixture violates native band'
          else
            band(j-i+1,i)=dense(i,j)
          endif
        enddo
        truth(i)=sin(real(i,rkind))
        if(active(i) == 0) truth(i)=0
      enddo
      original=band
      call folded_factor(band,ns,nt,n,m,mp,nodes,1.0d-10,info)
      if(info /= 0 .or. .not. folded_active) error stop 'Folded factorization failed'
      if(any(band /= original)) error stop 'Native matrix modified'
      do k=1,3
        rhs=matmul(dense,truth)*k
        old_rhs=rhs
        call folded_solve(rhs,info)
        if(info /= 0) error stop 'Folded solve rejected valid native RHS'
        if(maxval(abs(rhs-k*truth)) > 1.0d-12) error stop 'Known physical solution differs'
        if(maxval(abs(matmul(dense,rhs)-old_rhs)) > 2.0d-12) error stop 'Dense residual fails'
      enddo
      rhs=0
      rhs(1)=1
      old_rhs=rhs
      call folded_solve(rhs,info)
      if(info /= -2 .or. any(rhs /= old_rhs)) error stop 'Invalid inactive RHS accepted'
      band(2,1)=1
      call folded_factor(band,ns,nt,n,m,mp,nodes,1.0d-10,info)
      if(info /= -2 .or. folded_active) error stop 'Inactive coupling accepted'
      call folded_reset()
      deallocate(nodes,active,dense,band,original,truth,rhs,old_rhs)
    enddo
  enddo
  print *, 'PASS physical cell cliques, periodic seam, center border, outer jets, multiple RHS and fail-closed guards'
contains
  subroutine couple(i,j)
    integer,intent(in) :: i,j
    real(rkind) :: value
    value=.01_rkind*cos(real(3*i+j,rkind))
    dense(i,j)=dense(i,j)+value
    dense(j,i)=dense(i,j)
  end subroutine
end program
