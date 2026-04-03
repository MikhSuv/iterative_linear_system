program main
  use precision_mod
  use iterative_linear_system
  implicit none
  real(dp), allocatable :: A(:, :), B(:)
  integer :: n
  call read_linear_system("data.dat", A, B, n) 
  print *, A
  print *, is_diagonally_dominant(A)

end program main
