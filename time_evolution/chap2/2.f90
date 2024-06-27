program answer
    implicit none
    real(8) :: t,x,v
    real(8) :: dx,dv
    real(8) :: to, xo, vo   !!! nステップ目の値を記憶させる変数
    real(8),parameter :: pi=acos(-1.d0), omega=2.d0*pi   !!! dは倍精度をあらわす
    real(8),parameter :: dt = 1.0e-2_8, t_end=1.0e1_8   !!! 時間刻み0.01, 終了時刻10   !!! 倍精度の表し方はこちらでもok

    !!! 初期条件設定
    t = 0.0_8
    x = 0.0_8
    v = 2.0_8*pi

    !!! テキストファイル形式でデータを出力
    !!! 予め1.txtを作っておく必要はない。コード実行時に1.txtが作成される。
    open(10, file="2.txt")
    write(10,*) t,x,v   !!! 初期条件出力しておく
    close(10)
    
    do while(t<t_end)
        !!! まずnステップ目の値を記憶させる
        xo = x
        vo = v
        to = t
        
        !!! 1次精度オイラー法でdt/2 だけ更新。
        dx = vo*dt
        dv = -(omega**2)*xo*dt
        x = xo + 0.5_8*dx
        v = vo + 0.5_8*dv
        t = to + 0.5_8*dt

        !!! dtだけ更新。
        dx = v*dt
        dv = -(omega**2)*x*dt
        x = xo + dx
        v = vo + dv
        t = to + dt

        open(10, file="2.txt", position='append')
        write(10,*) t,x,v
        close(10)
    enddo
    
endprogram answer

!!! 時間刻み0.025, 計算終了40 などとすると、すこしずつ誤差が増えていくのが見える