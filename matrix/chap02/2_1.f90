program answer
    !!! この解答では、別のファイル(dat形式など)に行列やベクトルの情報を与えて
    !!! それを読み込むことで答えを返す、という想定をしている。
    
    implicit none
    real(8),allocatable :: a_raw(:,:), a(:,:)
    !!! b(:),x(:) は省いている。2-1では使わないためである。
    integer :: n,i,j,k

    read(*,*) n
    
    allocate(a_raw(n,n))
    read(*,*) a_raw
    a =transpose(a_raw)

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

    !!! 最初に与えたa(i,j)が
    ! 3 -1 2
    ! 5 1 -2
    ! 2 1 -1
    !!! だったとして(2_1の状況)、この時点でa(i,j)は、
    ! 3     -1     2
    ! 1.667 2.667 -5.333
    ! 0.667 0.625 1
    !!! のように、LとUが同居した形になっている。

!!! 下三角行列の出力(1行目は、1,0,0,0…確定なので、別の操作にする必要があると思ったが、不要)
!!! do j=1,0 はループの中に入らないという仕様のおかげ、なはず。
    
    write(*,*) '下三角行列Lは、'
    do i=1,n
        do j=1,i-1   !!! i=1のみこのループには入らない。
            write(*, '(f12.6)', advance='no') a(i,j)
        enddo
        write(*, '(i12)', advance='no') 1
        do j=i+1,n
            write(*, '(i12)', advance='no') 0
        enddo
        write(*,*)
    enddo

    write(*,*)

!!! 上三角行列の出力
    
    write(*,*) '上三角行列Lは、'
    do i=1,n
        do j=1,i-1
            write(*, '(i12)', advance='no') 0
        enddo
        do j=i,n
            write(*, '(f12.6)', advance='no') a(i,j)
        enddo
        write(*,*)
    enddo

    deallocate(a_raw)
    deallocate(a)

endprogram answer