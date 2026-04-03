module precision_mod

  implicit none
  private

  public :: dp, eps

  integer, parameter :: dp = selected_real_kind(15,307)
  real(dp), parameter :: eps = epsilon(1.0_dp)*100

end module precision_mod
