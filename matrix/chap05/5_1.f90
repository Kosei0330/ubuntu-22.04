program answer
    !!! この解答では、別のファイル(txt形式など)に行列の情報を与えて
    !!! それを読み込むことで答えを返す、という想定をしている。
        
    implicit none
    real(8),allocatable :: A_raw(:,:), A(:,:), A_k(:,:)
    real(8),allocatable :: a_i(:), u_i(:), Q(:,:),R(:,:)
    integer :: n,i,j,k,l
    real(8),parameter :: epsilon=0.001  !!! 打ち切り誤差
    real(8) :: c,d   !!! 収束判定などの変数
    
    read(*,*) n
        
    allocate(A_raw(n,n))
    read(*,*) A_raw
    A = transpose(A_raw)

    allocate(a_i(n))
    allocate(u_i(n))
    allocate(Q(n,n))
    allocate(R(n,n))
    Q = 0   !!!   適当な値を代入しておく
    R = 0

    do k=0,100   !!! 100はkmax（反復試行の最大回数）
        A_k = A
        do i=1,n   !!! Gram-Schmidt 直交化
            a_i = A(:,i)   !!! Aのi列を取り出す。
            
            u_i = a_i   !!! u_i は更新されていく変数
            do j=1,i-1
                u_i = u_i - dot_product(a_i, Q(:,j))*Q(:,j)
            enddo

            c = 0   !!! 規格化する
            do l=1,n
            c = c + u_i(l)**2
            enddo
            Q(:,i) = u_i / sqrt(c)   !!! 正規直交基底の行列Qのi列に値を代入する。
            R(i,i) = sqrt(c)   !!! ここがポイント、u_i が更新される前にR(i,i)に代入してしまう。
        enddo

        do i=1,n   !!! 上三角行列の計算
            do j=i+1,n
                R(i,j) = dot_product(A(:,j),Q(:,i))
            enddo
        enddo

        A = matmul(R,Q)

        d = 0   !!! 誤差の計算
        do i=1,n
            d = d + (A(i,i) - A_k(i,i))**2
        enddo
        if (d < epsilon**2) then
            exit
        endif

    enddo

    do i=1,n   !!! 固有値の出力
        write(*,*) A(i,i)
    enddo

endprogram answer