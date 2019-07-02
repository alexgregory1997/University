program plotting
  use plot_mod
  implicit none
  integer, parameter :: dp=selected_real_kind(15,300)
  real(kind=dp), parameter :: pi=3.14159265359_dp
  integer :: ios, input
  real(kind=dp) :: ose, ose2, t, k, a, s

  t=1.0_dp
  k=0.0_dp
  a=1.0_dp
  s=0.05_dp

  !select analytical solution to plot
  call selection(ose,ose2,input)

  open(unit=50,file='plotting.dat',iostat=ios)
  if (ios.ne.0) stop 'error opening plotting.dat'

  do while (k.le.pi/a)
    select case (input)
    case (1)
      write(50,*) k, energy_oapuc(ose,t,k,a)
      k=k+s
    case (2)
      write(50,*) k, energy_tapuc_pos(ose,ose2,t,k,a), energy_tapuc_neg(ose,ose2,t,k,a)
      k=k+s
    case (3)
      write(50,*) k, energy_tapuc_pos(ose,ose2,t,k,a), energy_tapuc_neg(ose,ose2,t,k,a)
      k=k+s
    case default
      print*, 'Please enter valid system'
    end select
  end do

  close(unit=50,iostat=ios)
  if (ios.ne.0) stop 'error closing plotting.dat'

end program plotting
