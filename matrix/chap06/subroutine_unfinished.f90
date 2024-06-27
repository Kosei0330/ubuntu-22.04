program answer
    implicit none

    
    stop
contains
    subroutine LU(a)
        implicit none
        real(8) :: a(n,n)
        integer :: i,j,k,n
        
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

    !!! y(n) への値の代入
    y(1) = b(1)
    do i=2,n
        c=0   !!! 和を整理するための変数
        do j=1,i-1
            c = c + a(i,j)*y(j)
        enddo
        y(i) = b(i)-c
    enddo
    
    !!! x(n)への値の代入
    x(n) = y(n) / a(n,n)
    do i= n-1,1,-1
        c=0
        do j=i+1,n
            c = c + a(i,j)*x(j)
        enddo
        x(i) = (y(i)-c) / a(i,i)
    enddo

    !!! x(解)の出力
    write(*,*) x

    endsubroutine
endprogram answer