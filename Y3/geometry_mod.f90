module geometry_mod
  implicit none
contains
  subroutine cartesian(i,n,ose,t,k,a,hamiltonian,option,col,row)
    implicit none
    integer, parameter :: dp=selected_real_kind(15,300)
    integer, intent(in) :: i, n, option, col, row
    real(kind=dp), intent(in) :: ose, t, k, a
    complex*16, dimension(:,:), allocatable, intent(inout) :: hamiltonian

    select case (option)
    case (1)
      hamiltonian(i,i)=hamiltonian(i,i)+ose+inter_left(t,k,a)+inter_right(t,k,a)
    case (2)
      hamiltonian(i,i)=hamiltonian(i,i)+ose
      if (i.eq.1) then
        hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
        hamiltonian(i,i+col-1)=hamiltonian(i,i+col-1)+inter_left(t,k,a)
      else if (i.eq.col) then
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
        hamiltonian(i,i-col+1)=hamiltonian(i,i-col+1)+inter_right(t,k,a)
      else
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
        hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
      end if
    case (3)
      hamiltonian(i,i)=hamiltonian(i,i)+ose
      if (i.eq.1) then
        hamiltonian(i,i)=hamiltonian(i,i)+inter_right(t,k,a)+inter_left(t,k,a)
        hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
      else if (i.eq.row) then
        hamiltonian(i,i)=hamiltonian(i,i)+inter_right(t,k,a)+inter_left(t,k,a)
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
      else
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
        hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
        hamiltonian(i,i)=hamiltonian(i,i)+inter_right(t,k,a)+inter_left(t,k,a)
      end if
    case (4)
      hamiltonian(i,i)=hamiltonian(i,i)+ose
      if (i.gt.1 .and. i.lt.col) then !top excluding corners
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
        hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
        hamiltonian(i,i+col)=hamiltonian(i,i+col)+intra(t)
      else if (i.gt.n-col+1 .and. i.lt.n) then !bottom excluding corners
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
        hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
        hamiltonian(i,i-col)=hamiltonian(i,i-col)+intra(t)
      else if (mod(i,col).eq.0 .and. i.gt.col .and. i.lt.n) then !right excluding corner
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
        hamiltonian(i,i+col)=hamiltonian(i,i+col)+intra(t)
        hamiltonian(i,i-col)=hamiltonian(i,i-col)+intra(t)
        hamiltonian(i,i-col+1)=hamiltonian(i,i-col+1)+inter_right(t,k,a)
      else if (mod(i,col).eq.1 .and. i.gt.1 .and. i.lt.n-col+1) then !left excluding corner
        hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
        hamiltonian(i,i+col)=hamiltonian(i,i+col)+intra(t)
        hamiltonian(i,i-col)=hamiltonian(i,i-col)+intra(t)
        hamiltonian(i,i+col-1)=hamiltonian(i,i+col-1)+inter_left(t,k,a)
      else if (i.eq.1) then !top left corner
        hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
        hamiltonian(i,i+col)=hamiltonian(i,i+col)+intra(t)
        hamiltonian(i,i+col-1)=hamiltonian(i,i+col-1)+inter_left(t,k,a)
      else if (i.eq.col) then !top right corner
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
        hamiltonian(i,i+col)=hamiltonian(i,i+col)+intra(t)
        hamiltonian(i,i-col+1)=hamiltonian(i,i-col+1)+inter_right(t,k,a)
      else if (i.eq.n-col+1) then !bottom left corner
        hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
        hamiltonian(i,i-col)=hamiltonian(i,i-col)+intra(t)
        hamiltonian(i,i+col-1)=hamiltonian(i,i+col-1)+inter_left(t,k,a)
      else if (i.eq.n) then !bottom right corner
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
        hamiltonian(i,i-col)=hamiltonian(i,i-col)+intra(t)
        hamiltonian(i,i-col+1)=hamiltonian(i,i-col+1)+inter_right(t,k,a)
      else !non-special case ie not on an edge
        hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
        hamiltonian(i,i+col)=hamiltonian(i,i+col)+intra(t)
        hamiltonian(i,i-col)=hamiltonian(i,i-col)+intra(t)
      end if
    case default
      stop 'Please select valid option'
    end select
  end subroutine cartesian

  subroutine zigzag(i,n,ose,t,k,a,hamiltonian)
    implicit none
    integer, parameter :: dp=selected_real_kind(15,300)
    integer, intent(in) :: i, n
    real(kind=dp), intent(in) :: ose, t, k, a
    complex*16, dimension(:,:), allocatable, intent(inout) :: hamiltonian

    hamiltonian(i,i)=hamiltonian(i,i)+ose
    if(i.eq.1) then
      hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)+inter_right(t,k,a)
    else if (i.eq.n) then
      if (mod(n,4).eq.2) then
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)+inter_left(t,k,a)
      else
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)+inter_right(t,k,a)
      end if
    else if (mod(i,4).eq.0) then
      hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)+inter_right(t,k,a)
      hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
    else if (mod(i,4).eq.1) then
      hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
      hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)+inter_right(t,k,a)
    else if (mod(i,4).eq.2) then
      hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)+inter_left(t,k,a)
      hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
    else if (mod(i,4).eq.3) then
      hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
      hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)+inter_left(t,k,a)
    end if
  end subroutine zigzag

  subroutine armchair(i,n,ose,t,k,a,hamiltonian)
    implicit none
    integer, parameter :: dp=selected_real_kind(15,300)
    integer, intent(in) :: i, n
    real(kind=dp), intent(in) :: ose, t, k, a
    complex*16, dimension(:,:), allocatable, intent(inout) :: hamiltonian

    hamiltonian(i,i)=hamiltonian(i,i)+ose
    if (i.eq.1) then
      hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
      hamiltonian(i,i+2)=hamiltonian(i,i+2)+intra(t)
    else if (i.eq.2) then
      hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
      hamiltonian(i,i+2)=hamiltonian(i,i+2)+intra(t)
    else if (i.eq.n) then
      if (mod(n,4).eq.0) then !fish tail
        hamiltonian(i,i-2)=hamiltonian(i,i-2)+intra(t)
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+inter_right(t,k,a)
      else if (mod(n,4).eq.2) then
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
        hamiltonian(i,i-2)=hamiltonian(i,i-2)+intra(t)
      end if
    else if (i.eq.n-1) then
      if (mod(n,4).eq.0) then !fish tail
        hamiltonian(i,i-2)=hamiltonian(i,i-2)+intra(t)
        hamiltonian(i,i+1)=hamiltonian(i,i+1)+inter_left(t,k,a)
      else if (mod(n,4).eq.2) then
        hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
        hamiltonian(i,i-2)=hamiltonian(i,i-2)+intra(t)
      end if
    else if(mod(i,2).eq.0) then !even
      if (mod(i,4).eq.0) then
        hamiltonian(i,i+2)=hamiltonian(i,i+2)+intra(t)
        hamiltonian(i,i-2)=hamiltonian(i,i-2)+intra(t)
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+inter_right(t,k,a)
      else
        hamiltonian(i,i+2)=hamiltonian(i,i+2)+intra(t)
        hamiltonian(i,i-2)=hamiltonian(i,i-2)+intra(t)
        hamiltonian(i,i-1)=hamiltonian(i,i-1)+intra(t)
      end if
    else !odd
      if (mod(i,4).eq.3) then
        hamiltonian(i,i+2)=hamiltonian(i,i+2)+intra(t)
        hamiltonian(i,i-2)=hamiltonian(i,i-2)+intra(t)
        hamiltonian(i,i+1)=hamiltonian(i,i+1)+inter_left(t,k,a)
      else
        hamiltonian(i,i+2)=hamiltonian(i,i+2)+intra(t)
        hamiltonian(i,i-2)=hamiltonian(i,i-2)+intra(t)
        hamiltonian(i,i+1)=hamiltonian(i,i+1)+intra(t)
      end if
    end if
  end subroutine armchair

  function intra(t)
    implicit none
    integer, parameter :: dp=selected_real_kind(15,300)
    real(kind=dp), intent(in) :: t
    real(kind=dp) :: intra
    intra=-t
  end function intra

  function inter_left(t,k,a)
    implicit none
    integer, parameter :: dp=selected_real_kind(15,300)
    real(kind=dp), intent(in) :: t, k, a
    complex(kind=dp) :: inter_left
    inter_left = -t*exp(-cmplx(0,1)*k*a)
  end function inter_left

  function inter_right(t,k,a)
    implicit none
    integer, parameter :: dp=selected_real_kind(15,300)
    real(kind=dp), intent(in) :: t, k, a
    complex(kind=dp) :: inter_right
    inter_right = -t*exp(cmplx(0,1)*k*a)
  end function inter_right

end module geometry_mod
