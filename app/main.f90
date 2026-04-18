program main

  use precision_mod
  use iterative_linear_system

  implicit none

  character(len=256):: arg
  real(dp), allocatable :: A(:, :), B(:), X(:)
  integer :: n

  if (command_argument_count() < 1) then
    error stop "Usage: ./iterative_slover <jacobi|seidel|relaxation>"
  end if
  call get_command_argument(1, arg)
  call read_linear_system("data.dat", A, B, n) 
  X = iterative_solve(A, B, trim(arg))
  call write_result("result.dat", X)
  print *,"Модуль вектора невязки:", residual(A, B, X)
end program main
