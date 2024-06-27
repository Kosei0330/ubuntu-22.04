program answer
    implicit none
    real(8) :: t,x,v
    real(8) :: dx,dv
    real(8),parameter :: pi=acos(-1.d0), omega=2.d0*pi   !!! dは倍精度をあらわす
    real(8),parameter :: dt = 1.0e-2_8, t_end=1.0e1_8   !!! 時間刻み0.01, 終了時刻10   !!! 倍精度の表し方はこちらでもok

    !!! 初期条件設定
    t = 0.0_8
    x = 0.0_8
    v = 2.0_8*pi

    !!! テキストファイル形式でデータを出力
    !!! 予め1.txtを作っておく必要はない。コード実行時に1.txtが作成される。
    open(10, file="1.txt")
    write(10,*) t,x,v   !!! 初期条件出力しておく
    close(10)
    
    do while(t<t_end)
        dx = v*dt
        dv = -(omega**2)*x*dt
        x = x + dx
        v = v + dv
        t = t + dt

        !!! 上のコードの代わりに以下の3行のコードをかくと、
        !!! 1次精度シンプレクティック数値積分法（p.26に記載あり）になる
        ! x = x + v*dt
        ! v = v - (omega**2)*x*dt
        ! t = t + dt

        open(10, file="1.txt", position='append')
        write(10,*) t,x,v
        close(10)
    enddo
    
endprogram answer

!!! 位置も速度も発散しているのが分かる！