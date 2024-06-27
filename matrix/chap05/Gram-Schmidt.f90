program answer
    !!! この解答では、別のファイル(txt形式など)に行列やベクトルの情報を与えて
    !!! それを読み込むことで答えを返す、という想定をしている。
        
    implicit none
    real(8),allocatable :: A_raw(:,:), A(:,:)
    real(8),allocatable :: a_i(:), u_i(:), Q(:,:)
    integer :: n,i,j,k,l
    real(8) :: c
    
    read(*,*) n
        
    allocate(A_raw(n,n))
    read(*,*) A_raw
    A = transpose(A_raw)

    allocate(a_i(n))
    allocate(u_i(n))
    allocate(Q(n,n))
    Q = 0   !!!   適当な値を代入しておく

        do i=1,n   !!! Gram-Schmidt 直交化
            a_i = A(:,i)   !!! Aのi列を取り出す。
            
            u_i = a_i
            do j=1,i-1
                u_i = u_i - dot_product(a_i, Q(:,j))*Q(:,j)
            enddo

            c = 0
            do l=1,n
            c = c + u_i(l)**2
            enddo
            c = sqrt(c)
            Q(:,i) = u_i / c
            
        enddo

    do i=1,n
        write(*,*) Q(:,i)
    enddo

endprogram answer