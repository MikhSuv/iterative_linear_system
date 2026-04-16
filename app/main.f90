program main
  use precision_mod
  use iterative_linear_system
  implicit none
  real(dp), allocatable :: A(:, :), B(:), X(:)
  integer :: n, i 
  ! call read_linear_system("data.dat", A, B, n) 
  n = 100
  allocate(A(n, n))
  allocate(B(n))
  call RANDOM_NUMBER(A)
  call RANDOM_NUMBER(B)
  do i = 1, n
  A(i,i) = a(i,i) * 1000
  end do
  print *, A
  print *, is_diagonally_dominant(A)
  X = jacoby_solve(A, B)
  print *, residual(A, B, x)

end program main
