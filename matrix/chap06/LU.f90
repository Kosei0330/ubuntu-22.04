program example
    implicit none
    real(8), allocatable :: a_raw(:,:), a(:,:)
    real(8),allocatable :: x(:),y(:),z(:)
    integer :: n

    read(*,*) n
    allocate(a_raw(n,n))
    read(*,*) a_raw
    a = transpose(a_raw)

    allocate(x(n),y(n),z(n))

    !!! 初期ベクトルの設定
    x = 0
    x(1) = 1

    !!! aを、LとUが同時に存在するような行列に変換するsubroutineを使う。
    !!! Lだけを出力したりはしないので、注意。
    call LU(a,n,x,y,z)
    write(*,*) a

    stop
contains
    subroutine LU(a,n,x,y,z)
        implicit none
        integer, intent(in) :: n
        real(8), intent(inout) :: a(n,n), x(n)
        real(8), intent(out) :: y(n), z(n)
        real(8) :: c
        integer :: i,j,k
        
        
        !!! LU分解の実行
        do i=1,n-1
            if (a(i,i) == 0) then
                exit
            endif
            do j=i+1,n
                a(j,i) = a(j,i) / a(i,i)
                do k=i+1,n
                    a(j,k)=a(j,k)-a(i,k)*a(j,i)
                enddo
            enddo
        enddo
        !!! この時点ではaはLとUが同居した行列になっている。
        
    endsubroutine

endprogram example