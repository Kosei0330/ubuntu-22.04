program answer
    !!! この解答では、別のファイル(txt形式など)に行列の情報を与えて
    !!! それを読み込むことで答えを返す、という想定をしている。
        
    implicit none
    real(8),allocatable :: A_raw(:,:), A(:,:), A_k(:,:)
    real(8),allocatable :: a_i(:), u_i(:), Q(:,:),R(:,:)
    integer :: n,i,j,k,l
    real(8),parameter :: epsilon=0.001  !!! 打ち切り誤差
    real(8) :: c,d   !!! 収束判定などの変数
    
    !!! プログラム後半で使う配列
    real(8), allocatable :: x(:),y(:),z(:),x_k(:)
    real(8), allocatable :: A_lambda(:)   !!! 固有値を入れておく配列
    real(8), allocatable :: A_original(:,:)   !!! もとのAを記憶する

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

    allocate(x(n),y(n),z(n),x_k(n), A_lambda(n))
    A_original = A

!!! QR分解の実行

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
            d = d + abs(A(i,i) - A_k(i,i))
        enddo
        if (d < epsilon) then
            exit
        endif

    enddo

    !!! 固有値の出力。
    !!! この時点でAの対角成分はAの固有値になっていることを利用する。
    do i=1,n   
        write(*,*) A(i,i)
    enddo
    
    !!! A_lambdaに固有値を記憶させる。
    do i=1,n
        A_lambda(i) = A(i,i)
    enddo

!!! LU分解、固有ベクトルを求める
    
    do l=1,n
        !!! A-μIの設定、以下ではAをA-μIとして使う。
        A = A_original
        do k=1,n
            A(k,k) = A(k,k) - A_lambda(l)
        enddo
        !!! この時点でAはA-μIになっている

        !!! LU分解の実行、aはAと同じ（fortranでは大文字小文字の区別がない）
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
        !!! この時点でAは、A-μIをLU分解した際の、LとUが同居した行列。

        !!! 初期ベクトルの設定
        x = 0
        x(1) = 1

        do k=0,1000   !!! 1000はkmax
            x_k = x !!! xを記憶しておく（収束判定の際に使う）
            
            !!! z(n) への値の代入（LU分解による解法の、中間生成物）
            z(1) = x(1)
            do i=2,n
                c=0   !!! 和を整理するための変数
                do j=1,i-1
                    c = c + a(i,j)*z(j)
                enddo
                z(i) = x(i)-c
            enddo
            
            !!! y(n)への値の代入
            y(n) = z(n) / a(n,n)
            do i= n-1,1,-1
                c=0
                do j=i+1,n
                    c = c + a(i,j)*y(j)
                enddo
                y(i) = (z(i)-c) / a(i,i)
            enddo
            
            !!! y(n)の規格化をする
            c=0
            do i=1,n
                c = c + y(i)**2
            enddo
            y = y /sqrt(c)
            !!! 規格化完了

            d = 0 !!! 誤差の計算
            do i =1,n
                d = d + (y(i)-x_k(i))**2
            enddo
            if (d < epsilon**2) then
                exit
            endif
            

            x = y   !!! 更新
        enddo
        
        !!! 固有ベクトルの出力
        write(*,*) y

    enddo
    

endprogram answer