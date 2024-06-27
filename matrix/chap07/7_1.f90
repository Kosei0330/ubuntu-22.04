program answer
    !!! 参照 https://sxauroratsubasa.sakura.ne.jp/documents/sdk/SDK_NLC/UsersGuide/man/dstev.html
    !!! 6_1_2.txt の10*10行列を想定する。
    !!! この課題のみ、計算機演習室のPCにてコンパイル。ライブラリをダウンロードした環境と、ファイルをコンパイルする環境は揃える。
    implicit none
    character :: JOBZ
    integer,parameter :: N=10
    real(8) :: D(N)   !!! 対角要素
    real(8) :: E(N-1)   !!! 非対角要素??
    real(8) :: Z(N,N)
    integer :: LDZ
    real(8) :: WORK(2*N-2)
    integer :: INFO

    JOBZ = 'V'
    D = -2.0_8
    E = 1.0_8
    LDZ = 10

    call DSTEV(JOBZ, N, D, E, Z, LDZ, WORK, INFO)

    write(*,*) INFO

endprogram answer