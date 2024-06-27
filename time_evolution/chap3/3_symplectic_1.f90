program answer
    !!! 1次精度シンプレクティック法による時間発展
    
    implicit none
    real(8) :: t,x,y,vx,vy
    real(8) :: dx,dy,dvx,dvy
    real(8),parameter :: dt = 2.0e-1_8, t_end=1.0e3_8   !!! 時間刻み0.2, 終了時刻1000
    
    !!! 初期条件設定
    t = 0.0_8
    x = 1.0_8
    y = 0.0_8
    vx = 0.0_8
    vy = 1.0_8

    !!! テキストファイル形式でデータを出力
    open(10, file="3_symplectic_1.txt", status='new')   !!! 'new'はファイルの新規作成
    write(10,*) x,y  !!! 初期条件出力しておく
    close(10)

    do while(t<t_end)   !!! dx,dyなど使わなくてもいける
        !!! 更新
        x = x + vx*dt
        y = y + vy*dt
        vx = vx - (x / (x**2+y**2)**(1.5_8))*dt
        vy = vy - (y / (x**2+y**2)**(1.5_8))*dt
        t = t + dt

        open(10, file="3_symplectic_1.txt", status='old', position='append')
        write(10,*) x,y
        close(10)

    enddo

endprogram answer