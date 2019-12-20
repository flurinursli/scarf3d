! this sample program is to test halo exchange

! program main
!
! use mpi
!
! implicit none
!
! integer                              :: i, j
! integer                              :: ierr, world_size, world_rank, topo
! integer                              :: down, up
! integer                              :: from_down, to_down, from_up, to_up
! integer                              :: tag_n, tag_s
! integer,              dimension(2)   :: npts, dims, coords, fs, fe
! integer,              dimension(2)   :: sizes, subsizes, starts
! integer,              dimension(4)   :: request
! integer,              dimension(mpi_status_size, 4) :: status
! integer, allocatable, dimension(:,:) :: gs, ge, n
! real,    allocatable, dimension(:,:) :: m1, m2
!
! !----------------------------------------------------
!
! call mpi_init(ierr)
!
! npts = [5, 7]
! dims = [2, 2]
!
! call mpi_comm_size(mpi_comm_world, world_size, ierr)
! call mpi_comm_rank(mpi_comm_world, world_rank, ierr)
!
! call mpi_cart_create(mpi_comm_world, 2, dims, [.false., .false.], .true., topo, ierr)
! call mpi_cart_coords(topo, world_rank, 2, coords, ierr)
!
! call mpi_rank2index(npts, dims, coords, fs, fe)
!
! do i = 0, world_size - 1
!    if (i .eq. world_rank) print*, world_rank, ' - ', fs(1), fe(1), ' - ', fs(2), fe(2)
!    call mpi_barrier(mpi_comm_world, ierr)
! enddo
!
! allocate(gs(2, 0:world_size-1), ge(2, 0:world_size-1))
!
! gs(:, world_rank) = fs
! ge(:, world_rank) = fe
!
! call mpi_allgather(mpi_in_place, 0, mpi_datatype_null, gs, 2, mpi_integer, mpi_comm_world, ierr)
! call mpi_allgather(mpi_in_place, 0, mpi_datatype_null, ge, 2, mpi_integer, mpi_comm_world, ierr)
!
! allocate(n(2, 0:world_size-1))
!
! ! points along each dimension for each process
! n = ge - gs + 1
!
! allocate(m1(n(1, world_rank), n(2, world_rank)))
!
! do j = 1, n(2, world_rank)
!    do i = 1, n(1, world_rank)
!       m1(i, j) = gs(1, world_rank) + i - 1 + (gs(2, world_rank) + j - 2)*npts(1)
!    enddo
! enddo
!
! do i = 0, world_size - 1
!    if (i .eq. world_rank) then
!       print*, world_rank
!       do j = 1, n(1, world_rank)
!          print*, m1(j,:)
!       enddo
!    endif
!    call mpi_barrier(mpi_comm_world, ierr)
! enddo
!
! ! allocate array with halo
! allocate(m2(0:n(1, world_rank) + 1, 0:n(2, world_rank) + 1))
!
! m2 = 0.
! m2(1:n(1, world_rank), 1:n(2, world_rank)) = m1
!
! ! loop over directions
! do j = 0, 1
!
! call mpi_cart_shift(topo, j, 1, down, up, ierr)
!
!    sizes = n(:, world_rank) + 2
!
!    ! exchange along N-S
!    if (j .eq. 0) then
!
!       ! "up" is process above, "down" is process below
!
!       subsizes = [1, n(2, world_rank)]
!
!       starts = [1, 1]
!       call mpi_type_create_subarray(2, sizes, subsizes, starts, mpi_order_fortran, mpi_real, to_down, ierr)
!       call mpi_type_commit(to_down, ierr)
!
!       starts = [0, 1]
!       call mpi_type_create_subarray(2, sizes, subsizes, starts, mpi_order_fortran, mpi_real, from_down, ierr)
!       call mpi_type_commit(from_down, ierr)
!
!       starts = [n(1, world_rank), 1]
!       call mpi_type_create_subarray(2, sizes, subsizes, starts, mpi_order_fortran, mpi_real, to_up, ierr)
!       call mpi_type_commit(to_up, ierr)
!
!       starts = [n(1, world_rank)+1, 1]
!       call mpi_type_create_subarray(2, sizes, subsizes, starts, mpi_order_fortran, mpi_real, from_up, ierr)
!       call mpi_type_commit(from_up, ierr)
!
!    ! exchange along E-W
!    elseif (j .eq. 1) then
!
!       ! "up" is process to the right, "down" is process to the left
!
!       subsizes = [n(1, world_rank) + 2, 1]
!
!       starts = [0, 1]
!       call mpi_type_create_subarray(2, sizes, subsizes, starts, mpi_order_fortran, mpi_real, to_down, ierr)
!       call mpi_type_commit(to_down, ierr)
!
!       starts = [0, 0]
!       call mpi_type_create_subarray(2, sizes, subsizes, starts, mpi_order_fortran, mpi_real, from_down, ierr)
!       call mpi_type_commit(from_down, ierr)
!
!       starts = [0, n(2, world_rank)]
!       call mpi_type_create_subarray(2, sizes, subsizes, starts, mpi_order_fortran, mpi_real, to_up, ierr)
!       call mpi_type_commit(to_up, ierr)
!
!       starts = [0, n(2, world_rank)+1]
!       call mpi_type_create_subarray(2, sizes, subsizes, starts, mpi_order_fortran, mpi_real, from_up, ierr)
!       call mpi_type_commit(from_up, ierr)
!
!    endif
!
!    ! as it is set now, all the "up" sides are not exchanged (this is useful to avoid repeating interpolation)
!    ! --> can remove one row/column in "m2"
!
!    ! exchange with "down" processor
!    !call mpi_sendrecv(m2, 1, to_down, down, 0, m2, 1, from_down, down, 0, mpi_comm_world, status(:,1), ierr)
!    call mpi_sendrecv(m2, 1, to_down, mpi_proc_null, 0, m2, 1, from_down, down, 0, mpi_comm_world, status(:,1), ierr)
!
!    ! exchange with "up" processor
!    !call mpi_sendrecv(m2, 1, to_up, up, 0, m2, 1, from_up, up, 0, mpi_comm_world, status(:,1), ierr)
!    call mpi_sendrecv(m2, 1, to_up, up, 0, m2, 1, from_up, mpi_proc_null, 0, mpi_comm_world, status(:,1), ierr)
!
!    call mpi_type_free(to_down, ierr)
!    call mpi_type_free(from_down, ierr)
!    call mpi_type_free(to_up, ierr)
!    call mpi_type_free(from_up, ierr)
!
! enddo
!
! do i = 0, world_size - 1
!    if (i .eq. world_rank) then
!       print*, world_rank
!       do j = 0, n(1, world_rank) + 1
!          print*, m2(j,:)
!       enddo
!    endif
!    call mpi_barrier(mpi_comm_world, ierr)
! enddo
!
!
! call mpi_barrier(mpi_comm_world, ierr)
! call mpi_finalize(ierr)
!
! end program main
!
! !------------------------------------------------------
!
! subroutine mpi_rank2index(npts, ntasks, rank, fs, fe)
!
! implicit none
!
! integer, dimension(2), intent(in)  :: npts
! integer, dimension(2), intent(in)  :: ntasks
! integer, dimension(2), intent(in)  :: rank
! integer, dimension(2), intent(out) :: fs, fe
! integer                            :: i
!
! !---------------------------------------------------------
!
! do i = 1, size(rank)
!   fs(i) = 1 + int( real(npts(i)) / real(ntasks(i)) * real(rank(i)) )
!   fe(i) = int( real(npts(i)) / real(ntasks(i)) * real(rank(i) + 1) )
! enddo
!
! end subroutine mpi_rank2index
