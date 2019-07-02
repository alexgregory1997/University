module switch_mod
  implicit none
contains
  subroutine selection(geometry,option,n,col,row)
    implicit none
    integer, intent(out) :: geometry, option, n, col, row

    print*, 'Please enter number corresponding to geometry:'
    print*, '1. Cartesian'
    print*, '2. Zigzag'
    print*, '3. Armchair'
    read*, geometry
    if (geometry.eq.1) then
      print*, 'Select unit cell type'
      print*, '1. Single atom per unit cell'
      print*, '2. Horizontal unit cell - max 1 atom in width'
      print*, '3. Verticle unit cell - max 1 atom length'
      print*, '4. All others with length and width greater than 1'
      read*, option
      select case(option)
      case (1)
        n=1
      case (2)
        row=1
        print*, 'Enter number of columns of unit cell'
        read*, col
        n=col*row
      case (3)
        col=1
        print*, 'Enter number of rows of unit cell'
        read*, row
        n=col*row
      case (4)
        print*, 'Enter number of columns of unit cell'
        read*, col
        print*, 'Enter number of rows of unit cell'
        read*, row
        n=col*row
      case default
        stop 'Please enter valid unit cell type'
      end select
    else
      print*, 'Please enter number of atoms per unit cell - must be multiple of 2'
      read*, n
      if (mod(n,2).ne.0) then
        stop 'Enter valid number of atoms per unit cell'
      end if
    end if
  end subroutine selection
end module switch_mod
