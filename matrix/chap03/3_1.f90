program answer
    !!! この解答では、別のファイル(dat形式など)に行列やベクトルの情報を与えて
    !!! それを読み込むことで答えを返す、という想定をしている。
    implicit none
    real(8),allocatable :: a_raw(:,:), a(:,:), a1(:,:)
    real(8),allocatable :: b(:), b1(:), x(:)
    real(8),allocatable :: y(:)   !!! jacobiで使う配列
    integer :: n,i,j,k
    
    real(8),parameter :: epsilon=0.001  !!! 打ち切り誤差
    real(8) :: c,d1,d2   !!! 反復計算の際に用いる変数
    integer :: e_jacobi, e_gauss   !!! 繰り返し回数を表示するための変数

    read(*,*) n
    allocate(a_raw(n,n))
    read(*,*) a_raw
    allocate(b(n))
    read(*,*) b
    a =transpose(a_raw)
    a1 = a
    b1 = b   !!! a1,a2 はLU分解にて行列やベクトルの情報を使うときに呼び出す配列。
    allocate(x(n),y(n))
    
!!! Jacobi反復法の実装

    x = 0   !!! 初期値の設定

    do k=1,1000   !!! 1000は反復回数の上限(kmax)
        d1 = 0   !!! 収束判定に用いるノルムの初期値
        do i=1,n
            d1 = d1 + (x(i))**2   !!! ノルムの絶対値の2乗和
        enddo
        
        !!! xを漸化式に従って更新する
        do i=1,n
            c = 0
            do j=1,n
                c = c + a(i,j)*x(j)
            enddo
            y(i) = c
        enddo
        do i=1,n
            x(i) = x(i) + (b(i) - y(i)) / a(i,i)
        enddo
        !!! 更新終了

        d2 = 0
        do i=1,n
            d2 = d2 + (x(i))**2   !!! ノルムの絶対値の2乗和
        enddo

        !!! 収束判定
        if(abs(d2-d1) < epsilon**2) then
            e_jacobi = k
            exit
        endif
    enddo

    write(*,*) 'Jacobi反復法による繰り返し回数は、', e_jacobi
    write(*,*) x
    write(*,*)   !!! 改行

!!! Gauss-Seidel反復法の実装 

    x = 0   !!! 初期値の設定

    do k=1,1000   !!! 1000は反復回数の上限(kmax)
        d1 = 0   !!! 収束判定に用いるノルムの初期値
        do i=1,n
            d1 = d1 + (x(i))**2   !!! ノルムの絶対値の2乗和
        enddo
        
        !!! xを漸化式に従って更新する
        do i =1,n
            c=0
            do j=1,n
                c = c + a(i,j)*x(j)
            enddo
            x(i) = x(i) + (b(i) - c) / a(i,i)
        enddo

        d2 = 0
        do i=1,n
            d2 = d2 + (x(i))**2   !!! ノルムの絶対値の2乗和
        enddo

        !!! 収束判定
        if(abs(d2-d1) < epsilon**2) then
            e_gauss = k
            exit
        endif
    enddo

    write(*,*) 'Gauss-Seidel反復法による繰り返し回数は、', e_gauss
    write(*,*) x
    write(*,*)

!!! ガウスの消去法による解

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

    write(*,*) 'ガウスの消去法による解は、'
    write(*,*) x
    write(*,*)


!!! LU分解による解
    a = a1
    b = b1   !!! aとbに元の行列(a1)とベクトル(b1)を代入する
    
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

    write(*,*) 'LU分解法による解は、'
    write(*,*) x

endprogram answer