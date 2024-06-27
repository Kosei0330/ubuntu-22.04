program answer
!!! この解答では例えば,1(1)を解く際に

! 3

! 3 -1 2
! 5 1 -2
! 2 1 -1

! 1
! 2
! 3

!!! のようなdat形式のファイルを別途用意する、ということを想定する。
!!! 1(2)は、
! 4

! 3.141592653 1 1 1
! 1 -1.41421356 -1 -1.41421356
! -2 -1 1 1
! 1 -1 -4 2.82842712
    
! 11.141592653
! -7
! 3.17157288
! -36.41421356
!!! のようになる。
    
    implicit none
    real(8),allocatable :: a_raw(:,:), a(:,:)
    real(8),allocatable :: b(:), x(:)
    integer :: n,i,j,k

    read(*,*) n
    
    allocate(a_raw(n,n))
    read(*,*) a_raw
    allocate(b(n))
    read(*,*) b
    a =transpose(a_raw)

    allocate(x(n))

    do i=1,n-1
        if (a(i,i) == 0) then
            exit
        endif
        do j=i+1,n
            do k=i+1,n
                a(j,k) = a(j,k) - a(i,k) * (a(j,i) / a(i,i))
            enddo
            b(j) = b(j) - b(i) * (a(j,i) / a(i,i))
        enddo
    enddo

    x(n) = b(n) / a(n,n)
    do i= n-1,1,-1
        x(i) = b(i)
        do j = i+1, n
            x(i) = x(i) - a(i,j) * x(j)
        enddo
        x(i) = x(i) / a(i,i)
    enddo

    write(*,*) x

    deallocate(a_raw)
    deallocate(a)
    deallocate(b)
    deallocate(x)
endprogram answer