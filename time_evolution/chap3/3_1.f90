program answer
    !!! ルンゲクッタ法による時間発展
    
    implicit none
    real(8) :: t,x,y,vx,vy
    real(8) :: dx,dy,dvx,dvy
    real(8) :: to, xo, yo, vxo,vyo  !!! nステップ目の値を記憶させる変数
    real(8),parameter :: dt = 2.0e-1_8, t_end=1.0e3_8   !!! 時間刻み0.2, 終了時刻1000
    
    !!! 初期条件設定
    t = 0.0_8
    x = 1.0_8
    y = 0.0_8
    vx = 0.0_8
    vy = 1.0_8

    !!! テキストファイル形式でデータを出力
    open(10, file="3_1.txt", status='new')   !!! 'new'はファイルの新規作成
    write(10,*) t,x,y  !!! 初期条件出力しておく
    close(10)

    do while(t<t_end)
        !!! まずnステップ目の値を記憶させる
        xo = x
        yo = y
        vxo = vx
        vyo = vy
        to = t

        !!! 1次精度オイラー法でdt/2 だけ更新。
        dx = vxo*dt
        dy = vyo*dt
        dvx = -xo / (xo**2+yo**2)**(1.5_8)*dt
        dvy = -yo / (xo**2+yo**2)**(1.5_8)*dt
        x = xo + 0.5_8*dx
        y = yo + 0.5_8*dy
        vx = vxo + 0.5_8*dvx
        vy = vyo + 0.5_8*dvy
        t = to + 0.5_8*dt
        
        !!! dtだけ更新。
        dx = vx*dt
        dy = vy*dt
        dvx = -x / (x**2+y**2)**(1.5_8)*dt
        dvy = -y / (x**2+y**2)**(1.5_8)*dt
        x = xo + dx
        y = yo + dy
        vx = vxo + dvx
        vy = vyo + dvy
        t = to + dt

        open(10, file="3_1.txt", status='old', position='append')
        write(10,*) t,x,y
        close(10)

    enddo

endprogram answer