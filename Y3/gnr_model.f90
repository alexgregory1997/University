program gnr_model
use geometry_mod
use switch_mod
  implicit none
  integer, parameter :: dp=selected_real_kind(15,300)
  real(kind=dp), parameter :: pi=3.14159265359_dp
  real(kind=dp) :: ose, t, k, a
  complex*16, dimension(:,:), allocatable :: hamiltonian
  integer :: i, ierr, n, j, ios, col, row
  real(kind=dp) :: s
  !zheev variable decleration
  character*1 :: jobz, uplo
  integer :: info, lwork
  complex*16, dimension(:), allocatable :: work
  real(kind=dp), dimension(:), allocatable :: rwork
  real(kind=dp), dimension(:), allocatable :: eigenvalues
  !user inputs
  integer :: geometry, option

  call selection(geometry,option,n,col,row)

  !initialisation
  !condition related
  ose=0.0_dp
  t=1.0_dp
  a=1.0_dp
  !zheev related
  lwork=max(1,2*n-1)
  !loop related
  s=0.005_dp
  k=0.0_dp

  !allocate arrays
  allocate(hamiltonian(n,n),stat=ierr)
  if (ierr.ne.0) stop 'error allocating hamiltonian'
  !allocate zheev arrays
  allocate(eigenvalues(n),stat=ierr)
  if (ierr.ne.0) stop 'Error Allocating Eigenvalues'
  allocate(work(max(1,lwork)),stat=ierr)
  if (ierr.ne.0) stop 'Error Allocating work'
  allocate(rwork(max(1,3*n-2)),stat=ierr)
  if (ierr.ne.0) stop 'Error Allocating rwork'

  !initialse arrays
  hamiltonian=cmplx(0.0_dp,0.0_dp)

  open(unit=51,file='band.dat',iostat=ios)
  if (ios.ne.0) stop 'error opening band.dat'

  do while (k.le.pi/a)
    do i=1,n
      select case (geometry)
      case(1)
        call cartesian(i,n,ose,t,k,a,hamiltonian,option,col,row)
      case(2)
        call zigzag(i,n,ose,t,k,a,hamiltonian)
      case(3)
        call armchair(i,n,ose,t,k,a,hamiltonian)
      case default
        stop 'Please select valid geometry'
      end select
    end do
    call zheev('N','L',n,hamiltonian,n,eigenvalues,work,lwork,RWORK,INFO)
    write(51,*) k, eigenvalues
    k=k+s
    hamiltonian=cmplx(0.0_dp,0.0_dp)
  end do

  close(51,iostat=ios)
  if (ios.ne.0) stop 'error closing band.dat'

  !deallocate arrays
  deallocate(hamiltonian,stat=ierr)
  if (ierr.ne.0) stop 'error deallocating hamiltonian'
  !deallocate zheev arrays
  deallocate(eigenvalues,stat=ierr)
  if (ierr.ne.0) stop 'error deallocating eigenvalues'
  deallocate(work,stat=ierr)
  if (ierr.ne.0) stop 'error deallocating work'
  deallocate(rwork,stat=ierr)
  if (ierr.ne.0) stop 'error deallocating rwork'

end program gnr_model
