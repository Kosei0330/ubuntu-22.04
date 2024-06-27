program poisson_solver
    !!! 中心差分法を用いて電位の値を繰り返し更新していく
    !!! 「繰り返し」というのが時間発展の課題と異なるところ。
    
    implicit none
    integer, parameter :: nx = 100, ny = 100
    real(8), parameter :: dx = 0.01, dy = 0.01
    integer, parameter :: kmax = 1000   !!! 計算打ち切り回数
    real(8), parameter :: tol = 1.0e-6   !!! 誤差を測るための変数
    real(8) :: error   !!! 誤差を測るための変数
    real(8) :: phi(0:nx, 0:ny), rho(0:nx, 0:ny)   !!! 電位と電荷の2次元配列(101*101行列)
    integer :: i, j, k

    !!! 境界条件と電荷の設定
    call initialize(phi, rho, nx, ny)

    !!! ポアソン方程式を実際に解く
    do k=1,kmax
        error = 0.0_8
        call poisson_step(phi, rho, nx, ny, dx, dy, error)
        if (error < tol) then
            exit
        endif
    enddo

    !!! 値の出力 
    !!! 配列の特性上、doループの回し方がj→iの順番
    open(10, file="report_1a.txt",status="replace")
    do j=0,ny
        do i=0,nx
            write(10,'(f6.3)', advance="no") phi(i,j)
            write(10,'(a)', advance='no') ' '
        enddo
        write(10,*)
    enddo
    close(10)

    stop
contains

subroutine initialize(phi, rho, nx, ny)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(out) :: phi(0:nx, 0:ny), rho(0:nx, 0:ny)
    integer :: i, j

    phi = 0.0_8   !!! いったん全部0にする
    phi(:,0) = -1.0_8
    phi(:,ny) = 1.0_8
    do j = 0,nx
        phi(0,j) = -1.0_8 + j/50.0_8
        phi(nx,j) = -1.0_8 + j/50.0_8
    enddo

    rho = 0.0_8   !!! いったん全部0にする
    do i = 0,nx
        do j = 0,ny
            if ((i/100.0_8 - 0.25)**2 + (j/100.0_8 - 0.25)**2 <= 1.0e-2_8) then   !!! 半径0.1の円盤上の電荷
                rho(i,j) = 1.0e2_8
            endif
        enddo
    enddo
endsubroutine initialize

subroutine poisson_step(phi, rho, nx, ny, dx, dy, error)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: dx, dy
    real(8), intent(inout) :: phi(0:nx, 0:ny)
    real(8), intent(in) :: rho(0:nx, 0:ny)
    real(8), intent(out) :: error
    integer :: i, j
    real(8) :: phi_new(0:nx, 0:ny)

    error = 0.0_8

    !!! 電位の更新＆誤差の計算
    do i = 1,nx-1
        do j = 1,ny-1
            phi_new(i, j) = 0.25d0 * (phi(i+1, j) + phi(i-1, j) + phi(i, j+1) + phi(i, j-1) + (dx*dy) * rho(i, j))
            error = max(error, abs(phi_new(i,j)-phi(i,j)))
        enddo
    enddo

    !!! 電位を更新
    do i = 1,nx-1
        do j = 1,ny-1
            phi(i, j) = phi_new(i, j)
        enddo
    enddo
endsubroutine poisson_step

endprogram poisson_solver
