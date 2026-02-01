program kaiseki
    implicit none !暗黙の型宣言禁止

    integer, parameter :: DivisionNumber = 720 !ローター分割数 回転方向の分割
    integer, parameter :: MaxDepth = 10 !奥行き分割量　軸方向への分割
    real(8), parameter :: RotorWeight = 29.15 !ローター質量 [kg]
    real(8), parameter :: SilicaWeight = 0.8 * RotorWeight !シリカゲルのみの重さ ローター全体の重さのうち8割がシリカゲルの質量だと考えている。
    real(8), parameter :: RotorSurfaceArea = 268.083 !ローター全体の表面積 [m^2]　ローターの幾何学的な形状から算出

    !物性値等の定数部　
    real(8), parameter :: SilicaWeightDivided = SilicaWeight / dble(DivisionNumber) !回転方向分割領域のシリカゲルの質量 [kg]
    real(8), parameter :: RotorWeightDivided = RotorWeight /dble(DivisionNumber) !回転方向分割領域のローターの質量 [kg]
    real(8), parameter :: A_ra = RotorSurfaceArea / dble(DivisionNumber) !回転方向分割領域における空気とローターが接する表面積 [m^2]
    real(8) :: SilicaWeight_cell !軸方向、回転方向に分割した領域のシリカゲルの質量
    real(8) :: RotorWeight_cell !軸方向、回転方向に分割した領域のローターの質量
    real(8) :: A_cell, Weight_cell 
    real(8), parameter :: Cp_r = 921.0 !ローターの比熱 [J/kg·K]
    real(8), parameter :: Cr_w = 4200.0 !水の比熱 [J/kg·K]
    real(8), parameter :: Cp_air = 1006.0 !空気比熱 [J/kg·K]
    real(8), parameter :: h_sa = 60.13 !熱伝達率 [W/m^2·K]
    real(8), parameter :: k = 0.06375 !物質移動係数 [kg/m^2·s(kg/kgDA)]
    real(8), parameter :: L = 2257000.0 !蒸発潜熱 [J/kg]
    real(8), parameter :: kp = 0.0056 !Polanyi DR 定数 [(cm^3/mol·K)^2]
    real(8), parameter :: Vm = 18.0 !水分子のモル容量 [cm^3/mol]

    !運転条件
    real(8) :: omega = 50.0 !ローターの1時間当たりの回転数 [rph]
    real(8) :: qa_process = 1030.0 !処理側風量 [m^3/h]
    real(8) :: qa_regenerate = 1000.0 !再生側風量 [m^3/h]
    real(8) :: qa_purge = 400 !パージ風量 [m^3/h]
    real(8) :: Tai_process = 28.3   !処理側流入温度 [℃]
    real(8) :: Tai_regenerate = 68.6 !再生側流入温度 [℃]
    real(8) :: Tai_purge = 18.8
    real(8) :: X_inProcess = 0.0192 !処理側流入絶対湿度 [kg/kgDA]
    real(8) :: X_inRegenerate = 0.0153 !再生側流入絶対湿度 [Kg/kgDA]
    real(8) :: X_inPurge = 0.0092
    real(8) :: operation_time = 3600.0 !稼働時間 [s]
    real(8) :: rho, P_atm, Cp_mix  !空気密度 (空気中の湿度と温度を考慮) [kg/m^3]

    !初期値
    real(8) :: X_in
    real(8) :: X_ao = 0.010998655 !出口絶対湿度 [kg/kgDA]
    real(8) :: X_s = 0.007422933 !シリカゲル表面絶対湿度 [kg/kgDA]
    real(8) :: T_r = 30.0 !ローター温度 [℃]
    real(8) :: T_ao = 30.0 !出口温度 [℃]
    real(8) :: W = 0.5
    real(8) :: W0 = 1.0
    real(8) :: RHs = 37.72545 !シリカゲル表面相対湿度 [%]
    real(8) :: es, e0, e
    real(8) :: Tai 
    real(8) :: time_n !n分割した際の1ステップの時間 [s]
    real(8) :: qa
    real(8) :: W_all
    real(8) :: W_all_old
    real(8) :: m_dry_process
    real(8) :: m_dry_regenerate
    real(8) :: m_dry_purge
    real(8) :: m_dry_flow
    real(8) :: m_water_proces
    real(8) :: m_water_regenerate
    real(8) :: m_water_purge
    real(8) :: error_m_water
    real(8) :: error_W_all
    real(8) :: Q_mix
    real(8) :: T_mix_ai
    real(8) :: T_mix_ao
    real(8) :: X_mix_ai
    real(8) :: X_mix_ao
    real(8) :: m_dry_mix_ai
    real(8) :: m_water_mix

    !インデックス, 配列　等
    integer :: i, remainder, idx
    integer :: cal_count 
    integer :: depth
    integer :: StartDepth
    integer :: EndDepth
    integer :: Direction

    real(8) ::  Tao_process_Ave
    real(8) :: Tao_regenerate_Ave
    real(8) :: Tao_purge_Ave
    real(8) :: Xao_process_Ave
    real(8) :: Xao_regenerate_Ave
    real(8) :: Xao_purge_Ave

    real(8) :: Tr_list(DivisionNumber, MaxDepth) !ローター温度の配列
    real(8) :: W_list(DivisionNumber, MaxDepth) !含湿量の配列
    real(8) :: Tao_list(DivisionNumber,MaxDepth)
    real(8) :: Xao_list(DivisionNumber,MaxDepth)
    real(8) :: T_r_list_temporary(DivisionNumber,MaxDepth)
    real(8) :: W_list_temporary(DivisionNumber,MaxDepth)

    !関数
    real(8) density

    !ファイルの読み込み
    open(20, file='error.dat', status='replace')
    !open(21, file='W_position.csv', status='replace')

    time_n = 3600.0d0 / (dble(DivisionNumber) * omega) !１ステップあたりの時間[s]
    cal_count = int(operation_time / time_n) !稼働時間分のステップ数
    P_atm = 101325.0d0 !標準大気圧 [Pa]

    !軸方向の分割について
    A_cell = A_ra / dble(MaxDepth)
    RotorWeight_cell = RotorWeightDivided / dble(MaxDepth)
    SilicaWeight_cell = SilicaWeightDivided / dble(MaxDepth)

    !ローターの初期値設定

    Tr_list = 30.0d0
    W_list = 0.05d0

    !タイムステップの計算 開始
    do i = 1, cal_count

        !n分割分それぞれの計算 開始
        do idx = 1, DivisionNumber
                
            !******************************************************************************************************
            !ローター入口温湿度の決定
            remainder = mod(idx - 1, DivisionNumber)
            !余りが0から0.1までは、パージ域、0.1～0.5まではqa = qa_process(処理空気)、 それ以外は qa = qa_regenerate(再生空気)
            if (remainder < int(DivisionNumber * 0.1)) then 
                qa = (qa_purge / 3600.0) / (dble(DivisionNumber) * 0.1)
                X_in = X_inPurge
                Tai = Tai_purge
                StartDepth = 1
                EndDepth = MaxDepth
                Direction = 1
            elseif (remainder < int(DivisionNumber * 0.5)) then
                qa = (qa_process / 3600.0d0) / ( dble(DivisionNumber) * 0.4 )
                X_in = X_inProcess
                Tai = Tai_process
                StartDepth = 1
                EndDepth = MaxDepth
                Direction = 1
            else
                qa = (qa_regenerate / 3600.0d0) / ( dble(DivisionNumber) * 0.5 )
                X_in = X_inRegenerate
                Tai = Tai_regenerate
                StartDepth = MaxDepth
                EndDepth = 1
                Direction = -1
            endif
            !******************************************************************************************************

            !定数
            !湿り空気密度の計算
            !e = P * X / (0.622 + X)
            e0 = P_atm * X_in / (0.622d0 + X_in) !水蒸気分圧 [Pa]
            rho = (28.966d0 * P_atm - (28.966d0 - 18.0d0) * e0) / (8.314d0 * (Tai + 273.15d0)) / 1000.0d0 !流入空気の密度 [kg/m^3]
            m_dry_flow = (rho * qa) / (1.0d0 + X_in) !乾燥空気ベースで計算を行う。kg(DA)/s

            !奥行き分割領域　計算開始
            do depth = StartDepth, EndDepth, Direction
                !パージ風量0の時、パージ域は何の変化もしないので、回転させるだけ。
                if(remainder < int(DivisionNumber * 0.1) .and. qa_purge <= tiny(1.0d0)) then
                    Tao_list(idx,depth) = 0
                    Xao_list(idx,depth) = 0
                    T_r_list_temporary(idx,depth) = Tr_list(idx,depth)
                    W_list_temporary(idx,depth) = W_list(idx,depth)
                    cycle
                endif

                !初期値を持ってくる
                T_r = Tr_list(idx,depth)!ローター温度の初期値
                W = W_list(idx,depth) !含湿量の初期値

                !Polanyi DR式を相対湿度について解いた式を用いて、含湿量からシリカゲルの表面絶対湿度Xsを計算する。
                Es = 6.1078 * 10.0**(7.5 * T_r / (T_r + 237.3)) * 100.0 !飽和水蒸気圧 [Pa]
                RHs = 100*(W / W0)**(Vm**2 / (2*kp*((273.16+T_r)**2)) ) !RHsについて

                if (RHs > 100) then
                    RHs = 99.999
                    write(20, *) "RHs is higher than 100 %. Time step:", i, "Position:", &
                    idx, "Axis position:", depth
                elseif (RHs <= 0) then
                    RHs = 0.0001
                    write(20, *) "RHs is less than 0 %. Time step:", i, "Position:", idx, "Axis position:", depth 
                endif

                e = RHs * Es / 100 !空気の水蒸気分圧 [Pa]
                X_s = 18*e / (29*(101325.0 - e)) !シリカゲル表面絶対湿度

                !出口湿度も計算する
                X_ao = (m_dry_flow * X_in + k * A_cell * X_s) / (m_dry_flow + k * A_cell)
                !まずはひとつ前のローターの温度を取ってきて、それを用いて、出口温度を計算する
                T_ao = ( (rho * Cp_air * qa * Tai) + (h_sa * A_cell * T_r) ) /(rho * Cp_air * qa + h_sa * A_cell) !出口温度の計算

                !ここで、次の時刻のTr(i), W(i)を計算する
                Weight_cell = RotorWeight_cell + SilicaWeight_cell * W
                Cp_mix = (RotorWeight_cell * Cp_r + SilicaWeight_cell * W * Cr_w) / Weight_cell !シリカゲル比熱と吸着水分の比熱を足し合わせた
                T_r = ( Weight_cell * Cp_mix / time_n * T_r + h_sa * A_cell * T_ao + k * A_cell * L * (X_ao - X_s)) / &
                ( Weight_cell * Cp_mix / time_n + h_sa * A_cell) 

                W = time_n / SilicaWeight_cell * k * A_cell * (X_ao - X_s) + W
                
                if (W > W0) then
                    W = W0 * 0.999d0
                    write(20, *) "water content is higher than water content capacity in silicagel."
                elseif (W <= 0 ) then
                    W = 1.0e-8
                    write(20, *) "the amount of water content is minus. Time:", i, "Position:", idx, "Axis position:", depth 
                endif

                !そのまま現時点の位置にとどめていいもの
                Tao_list(idx,depth) = T_ao
                Xao_list(idx,depth) = X_ao

                !回転させる必要があるもの
                T_r_list_temporary(idx,depth) = T_r
                W_list_temporary(idx,depth) = W

                X_in = X_ao !ひとつ前の分割領域の出口温度を持ってくる。
                Tai = T_ao !同上


            enddo
            !軸方向の計算　終了
            
            !確認用
            !if (idx == int(DivisionNumber * 0.3)) then
                !write(*, *) i*time_n/60.0, ",", Vm**2 / (2*kp*(273.16+T_r)**2), ",", RHs, ",", X_s, ",", X_ao - X_s
            !endif
        enddo
        !回転方向分割の計算　終了

        !あるタイムステップでの出口の平均温度、平均絶対湿度
        Tao_process_Ave = sum(Tao_list(int(DivisionNumber*0.1)+1:DivisionNumber/2,MaxDepth)) / int(DivisionNumber*0.4)
        Tao_regenerate_Ave = sum(Tao_list(DivisionNumber/2+1:DivisionNumber,1)) / (DivisionNumber / 2)
        Tao_purge_Ave =sum(Tao_list(1:int(DivisionNumber*0.1),MaxDepth)) / int(DivisionNumber*0.1)
        Xao_process_Ave = sum(Xao_list(int(DivisionNumber*0.1)+1:DivisionNumber/2,MaxDepth)) / int(DivisionNumber*0.4)
        Xao_regenerate_Ave = sum(Xao_list(DivisionNumber/2+1:DivisionNumber,1)) / (DivisionNumber / 2)
        Xao_purge_Ave = sum(Xao_list(1:int(DivisionNumber*0.1),MaxDepth)) / int(DivisionNumber*0.1)


        !湿気の質量保存の確認 ここではW_listがひとつ前のステップの含湿量となる。
        W_all = (sum(W_list_temporary * SilicaWeight_cell) - sum(W_list * SilicaWeight_cell) ) / time_n 
        
        !処理、再生、パージそれぞれを単体で見る。
        m_dry_process =  density(X_inProcess, Tai_process) * (qa_process/3600.0) / (1 + X_inProcess) !乾燥空気は吸脱着の過程でも一定のはずなのでこれを指標とする
        m_dry_regenerate = density(X_inRegenerate, Tai_regenerate) * (qa_regenerate/3600.0) / (1 + X_inRegenerate)
        m_dry_purge = density(X_inPurge, Tai_purge) * (qa_purge/3600.0) / (1 + X_inPurge)
        m_water_proces = m_dry_process * (X_inProcess - Xao_process_Ave) !ローターに湿気が入ってくる方向を正とする
        m_water_regenerate = m_dry_regenerate * (X_inRegenerate -Xao_regenerate_Ave) !ローターに湿気が入ってくる方向を正とする
        m_water_purge = m_dry_purge * (X_inPurge - Xao_purge_Ave) !上と同じ

        error_m_water = m_water_proces + m_water_regenerate + m_water_purge
        error_W_all = W_all - W_all_old 

        W_all_old = W_all !今のやつを保存

        !パージ域と処理側域の混合
        !Q_mix = qa_process + qa_purge
        !T_mix_ai = (qa_process * Tai_process + qa_purge * Tai_purge) / Q_mix !ローター通過前
        !X_mix_ai = (qa_process * X_inProcess + qa_purge * X_inPurge) / Q_mix
        !T_mix_ao = (qa_process * Tao_process_Ave + qa_purge * Tao_purge_Ave) / Q_mix !ローター通過後
        !X_mix_ao = (qa_process * Xao_process_Ave + qa_purge * Xao_purge_Ave) / Q_mix
        !m_dry_mix_ai = (density(X_mix_ai, T_mix_ai) * qa_process / 3600.0) / (1.0d0 + X_mix_ai)
        !m_water_mix = m_dry_mix_ai * (X_mix_ai - X_mix_ao) !処理側の除湿量

        !次のステップで使うものを回す
        Tr_list = cshift(T_r_list_temporary, shift=-1, dim=1)
        W_list = cshift(W_list_temporary, shift=-1, dim=1)

        !結果出力　(時間（分）)、(処理側出口空気平均温度)、(再生側出口空気平均温度)、(処理側出口平均絶対湿度)、(再生側出口平均絶対湿度)
        !(吸脱着量の差), (ステップ間の全含湿量の差), (除湿量 kg/h), (加湿量 kg/h)
        if (i*time_n/30 == real(int(i*time_n/30))) then
            write(*, *) i*time_n/60.0, ",", Tao_process_Ave, ",", Tao_regenerate_Ave,&
            ",", Xao_process_Ave, ",", Xao_regenerate_Ave, ",", error_m_water, ",", error_W_all, &
            ",", (m_water_proces + m_water_purge) * 3600.0, ",", m_water_regenerate * 3600.0
        endif

        !if (i == 17900) then
            !do idx = 1, 360
                !write(*, *) idx, ",", Tr_list(idx), ",", W_list(idx)
            !enddo
        !endif

        !if (remainder == 90) then
            !write(*, *)  X_in, X_ao, X_in - X_ao, k * A_cell * L *(X_in - X_ao)
        !elseif (remainder == 270 ) then
            !write(*, *)  X_in, X_ao, X_in - X_ao, rho*Cp_mix*qa*(T_ao - Tai)
        !endif

    enddo
    !最後の含湿量の状態、出口絶対湿度の値から質量保存の式が成り立っているか確認したい。
    !結果出力 最後のタイムステップ　(密度)、(一時間当たりの除湿量)、(一時間当たりの加湿量)
    !write(*, *) rho, ",", qa_process * rho * (X_inProcess - Xao_process_Ave), ",", &
    !qa_regenerate * rho * (Xao_regenerate_Ave - X_inRegenerate)

    !ファイル関連
    close(20)
    !close(21)

    stop
end program kaiseki  

real(8) function density(x, t)
    real(8) x, t
    real(8) :: P_atm = 101325.0
    real(8) c
    c = P_atm * x / (0.622d0 + x) !水蒸気分圧 [Pa]
    density = (28.966d0 * P_atm - (28.966d0 - 18.0d0) * c) / (8.314d0 * (t + 273.15d0)) / 1000.0d0 !流入空気の密度 [kg/m^3]
end function