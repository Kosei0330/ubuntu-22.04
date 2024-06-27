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
        
        implicit none
        real(8),allocatable :: a_raw(:,:), a(:,:)
        real(8),allocatable :: b(:), x(:)
        integer :: n,i,j,k
    
    !!! a(j,i)の最大値を探査するための変数    
        integer :: m
        integer :: max_place_a_ji
        real(8) :: max_a_ji   
        real(8) :: p,q,r,s

        read(*,*) n
        
        allocate(a_raw(n,n))
        read(*,*) a_raw
        allocate(b(n))
        read(*,*) b
        a =transpose(a_raw)
    
        allocate(x(n))
    
        do i=1,n-1
            !!! 1_1との違いはここから !!!
            max_a_ji = 0
            max_place_a_ji = i
            do j=i+1,n
                if(abs(a(j,i)) >= max_a_ji) then
                    max_a_ji = abs(a(j,i))
                    max_place_a_ji = j
                endif
            enddo

            do m=1,n !!! 行の入れ替え
                p = a(i,m)
                q = a(max_place_a_ji,m)
                a(i,m) = q
                a(max_place_a_ji,m) = p
                r = b(i)   !!! 三浦先生からfeedbackあり。この行含め下4行はdoループの外でok
                s = b(max_place_a_ji)
                b(i) = s
                b(max_place_a_ji) = r
            enddo        
            !!! ここまで !!!
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