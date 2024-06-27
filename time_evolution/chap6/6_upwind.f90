program answer
    implicit none
    integer :: i,j
    integer,parameter :: N=300, mg=2 !!! 格子点の数, 計算領域外に用意する格子点の数
    real(8),parameter :: dx = 1.0_8 / real(N)   !!! 格子点幅
    real(8),parameter :: dt = 1.0e-3_8 , t_end = 1.0_8   !!! 時間刻みと計算終了時刻
    real(8) :: t
    real(8) :: x(-mg:N-1) , u(-mg:N-1)   !!! 位置座標、解くべき変数
    real(8) :: uo(-mg:N-1) , du(-mg:N-1)   !!! 1ステップ前の値を保存する配列、uの増分の配列
    integer :: n_out   !!! t=0のデータは"data0000.dat"に、t=dtのデータは"data0001.dat"に格納されるようにするための変数。（n_out = n_out +1 的な感じで使う）
    character(128) :: filename   !!! ファイル名を表す変数

    !!! 初期条件
    do i=0,ceiling(((float(N)-1)/2))
        u(i) = 1
    enddo
    do i = ceiling(((float(N)-1)/2))+1 ,N-1
        u(i) = -1
    enddo
    u(-mg) = u(N-mg)
    u(-mg+1) = u(N-mg+1)
    
    do i=-mg,N-1
        x(i) = i*dx
    enddo
    n_out = 0
    t = 0.0_8

    !!! 初期状態のデータを書き込み
    !!! dataというフォルダはあらかじめ作っておくこと。
    write(filename, '("./data_upwind/data", i4.4, ".dat")') n_out
    
    open(10, file=filename)
    do i=0,N-1
        write(10,*) x(i), u(i)
    enddo
    close(10)
    n_out = n_out + 1

    do while(t<t_end)  !!! 大枠はtでループを回す
        uo = u   !!! 記憶させる
        
        do j=0,N-1   !!! duは2次精度風上差分法で以下のように書ける。
            du(j) = -1.0_8 * (u(j-2) - 4.0_8*u(j-1) + 3.0_8*u(j)) * dt / (2.0_8*dx) 
        enddo
        do j=0,N-1   !!! dt/2 だけ更新
            u(j) = uo(j) + 0.5_8*du(j)
        enddo

        !!! 境界条件の更新忘れずに。
        u(-mg) = u(N-mg)
        u(-mg+1) = u(N-mg+1)

        do j=0,N-1   !!! 上で更新されたuを使ってduを定義。
            du(j) = -1.0_8 * (u(j-2) - 4.0_8*u(j-1) + 3.0_8*u(j)) * dt / (2.0_8*dx) 
        enddo
        do j=0,N-1   !!! dtだけ更新
            u(j) = uo(j) + du(j)
        enddo

        !!! 境界条件の更新忘れずに。
        u(-mg) = u(N-mg)
        u(-mg+1) = u(N-mg+1)

        t = t + dt

        write(filename, '("./data_upwind/data", i4.4, ".dat")') n_out

        open(20, file=filename)
        do i=0,N-1
            write(20,*) x(i), u(i)
        enddo
        close(20)

        n_out = n_out + 1
    enddo

endprogram answer