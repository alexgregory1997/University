module plot_mod
  implicit none
contains
  subroutine selection(ose,ose2,input)
    implicit none
    integer, parameter :: dp=selected_real_kind(15,300)
    integer, intent(out) :: input
    real(kind=dp), intent(out) :: ose, ose2
    print*, 'Select system to plot analytical solution'
    print*, '1. One atom per unit cell'
    print*, '2. Two identical atoms per unit cell'
    print*, '3. Two differing atoms per unit cell'
    read*, input
    select case (input)
    case (1)
      print*, 'please enter on-site energy'
      read*, ose
    case (2)
      print*, 'please enter on-site energy'
      read*, ose
    case (3)
      print*, 'please enter on-site energies'
      read*, ose
      read*, ose2
    case default
      print*, 'please select valid system'
    end select
  end subroutine selection

  function energy_oapuc(ose,t,k,a) !energy for one atom cell
    implicit none
    integer, parameter :: dp=selected_real_kind(15,300)
    real(kind=dp) :: ose, t, k, a
    real(kind=dp) :: energy_oapuc
    energy_oapuc=ose-2.0_dp*t*cos(k*a)
  end function energy_oapuc

  function energy_tapuc_pos(ose,ose2,t,k,a) !one solution of energy for two atom cell
    implicit none
    integer, parameter :: dp=selected_real_kind(15,300)
    real(kind=dp) :: ose, ose2, t, k, a
    real(kind=dp) :: energy_tapuc_pos
    energy_tapuc_pos=(ose+ose2+sqrt(((ose-ose2)**2.0_dp)+8.0_dp*(t**2.0_dp)*(1.0_dp+cos(k*a))))/2.0_dp
  end function energy_tapuc_pos

  function energy_tapuc_neg(ose,ose2,t,k,a) !other energy solution for two atom cell
    implicit none
    integer, parameter :: dp=selected_real_kind(15,300)
    real(kind=dp) :: ose, ose2, t, k, a
    real(kind=dp) :: energy_tapuc_neg
    energy_tapuc_neg=(ose+ose2-sqrt(((ose-ose2)**2.0_dp)+8.0_dp*(t**2.0_dp)*(1.0_dp+cos(k*a))))/2.0_dp
  end function energy_tapuc_neg
end module plot_mod
