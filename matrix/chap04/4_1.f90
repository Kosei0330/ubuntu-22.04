program answer
    !!! この解答では、別のファイル(txt形式など)に行列やベクトルの情報を与えて
    !!! それを読み込むことで答えを返す、という想定をしている。
        
    implicit none
    real(8),allocatable :: a_raw(:,:), a(:,:)
    real(8),allocatable :: x(:), y(:)
    integer :: n,i,j,k
    real(8),parameter :: epsilon=0.001  !!! 打ち切り誤差
    real(8) :: gamma1, gamma2   !!! 固有値を求めるための変数
    real(8) :: c !!! 絶対値を求めるための変数
    
    read(*,*) n
        
    allocate(a_raw(n,n))
    read(*,*) a_raw
    a = transpose(a_raw)
    allocate(y(n))
    allocate(x(n))
    x = 0
    x(1) = 1   !!! 初期値を与える

    do k = 0, 100  !!! 100は反復回数の上限(kmax)
        y = matmul(a,x)
        gamma1 = 0
        gamma2 = dot_product(x,y)
        c = 0
        do i=1,n
            c = c + y(i)**2
        enddo
        c = sqrt(c)
        y = y / c
        if (abs(gamma2-gamma1) < epsilon) then
            exit
        endif
        x = y
        gamma1 = gamma2

    enddo

    write(*,*) 'eigenvalue :',gamma2   !!! 絶対値の一番大きい固有値
    write(*,*) 'eigenvector :',x   !!! 固有ベクトル

    deallocate(x,y,a_raw,a)
endprogram answer