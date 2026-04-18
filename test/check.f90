program check
  use precision_mod
  use iterative_linear_system
implicit none

integer ::  i
integer :: n = 10000
real(dp), allocatable :: A(:, :), B(:), X(:)
real(dp) :: row_sum_offdiag, rand_val, diag_boost

allocate(A(n,n), B(n))
call random_number(A)
A = 10.0_dp * (A - 0.5_dp)

do i = 1, n
  row_sum_offdiag = sum(abs(A(i,:))) - abs(A(i,i))

  call random_number(rand_val)
  diag_boost = 1.0_dp + 5.0_dp * rand_val
  A(i,i) = sign(row_sum_offdiag + diag_boost, A(i,i))
  if (abs(A(i,i)) < 1.0e-12_dp) A(i,i) = row_sum_offdiag + diag_boost
end do

call random_number(B)

X = iterative_solve(A, B, "jacobi")
print *, "Jacobi residual: ", residual(A, B, X)
X = iterative_solve(A, B, "seidel")
print *, "Seidel residual: ", residual(A, B, X)
X = iterative_solve(A, B, "relaxation")
print *, "Relaxation residual: ", residual(A, B, X)

end program check
