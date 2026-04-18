module iterative_linear_system
  use precision_mod
  implicit none
  ! private
  integer, PARAMETER :: max_iter = 100000000
  ! public :: say_hello
contains
   subroutine read_linear_system(filename, A, B, n)
      character(len=*), intent(in) :: filename
      real(dp), allocatable, intent(out) :: A(:, :), B(:) ! СЛУ (A|B)
      integer, intent(out), optional :: n

      integer :: iunit, iostatus, i
      character(len=256) :: line ! для чтения первой строки

      open (newunit=iunit, file=filename, status='old', &
            action='read', iostat=iostatus)
      if (iostatus /= 0) then
         error stop 'Error occured while opening file'
      end if

      read (iunit, '(a)', iostat=iostatus) line

      if (iostatus /= 0) then
         error stop 'Error occured while reading line'
      end if
      read (line(2:), *) n

      ! Выделение памяти под матрицу A
      allocate (A(n, n))
      ! Выделение памяти под столбец B
      allocate (B(n))
      do i = 1, n
         read (iunit, *) A(i, :)
      end do

      do i = 1, n
         read (iunit, *) B(i)
      end do

      close (iunit)

   end subroutine read_linear_system

   subroutine write_result(filename, X)
      character(len=*), intent(in) :: filename
      real(dp), intent(in) ::  X(:) ! Столбец результата

      integer :: n
      integer :: ounit, i, iostatus

      n = size(X)
      open (newunit=ounit, file=filename, action='write', iostat=iostatus)
      if (iostatus /= 0) then
         error stop 'Error occured while opening file'
      end if

      write (ounit, '("# ", i0)') n

      do i = 1, n
         write (ounit, '(*(e23.15, 1x))') X(i)
      end do

      close (ounit)

   end subroutine write_result
  function is_diagonally_dominant(A) result(check)
    ! Проверяет, обладает ли матрица A диагоняльным преобладанием
    real(dp), intent(in) :: A(:, :) ! Матрица
    real(dp), allocatable :: diag_vals(:), row_summs(:)
    logical :: check
    integer :: n, i
    check = .false.
    n = size(A(:,1))

    allocate(diag_vals(n), row_summs(n))
    diag_vals = [(abs(A(i,i)), i = 1, n)]
    row_summs = sum(abs(A), dim=2) - diag_vals
    check = all(diag_vals > row_summs)
  end function is_diagonally_dominant

  function jacoby_solve(A, B) result(X)
    real(dp), intent(in) :: A(:, :), B(:)
    real(dp), allocatable :: X(:), X_new(:)
    real(dp), allocatable :: invD(:), Z(:,:), G(:)
    integer :: n, i, j
    
    if (.not. is_diagonally_dominant(A)) then
    print *, "WARDING: Not diagonally dominant matrix!"
    end if
    n = size(B)
    allocate(invD(n), G(n), Z(n,n))
    invD = 1.0_dp / [(A(i,i), i = 1, n)]
    G = B * invD
    X = G
    ! !$omp parallel do collapse(2) private(i, j) shared(A, invD, Z, G)
    do i = 1, n
      do j = 1, n
        if (i == j) then
        Z(i, j) = 0.0_dp
        else
          Z(i, j) = -invD(i) * A(i,j)
        end if
      end do
    end do
    ! !$omp end parallel do

    do i = 1, max_iter
      X_new = matmul(Z, X) + G
      if (norm2(X_new - X) < eps) then
        return
      end if
      X = X_new
    end do
    print *, "Ineration limit has been reached"
  end function jacoby_solve

  function residual(A, B, X) result(r)
    real(dp), intent(in) :: A(:, :), B(:), X(:)
    real(dp) :: r ! модуль вектора невязки

    r = norm2(matmul(A, X) - B)
  end function residual
end module iterative_linear_system
