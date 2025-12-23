program kaiseki
    implicit none !暗黙の型宣言禁止

    integer, parameter :: DivisionNumber = 360 !ローター分割数 
    real(8), parameter :: m = 29.15 !ローター質量 [kg]
    real(8), parameter :: A = 268.083 !ローター全体の表面積 [m^2]

    !物性値等の定数部　
    real(8), parameter :: m_s = 0.8 * m / dble(DivisionNumber) !シリカゲルの質量 [kg]
    real(8), parameter :: A_sa = A / dble(DivisionNumber) !空気とローターが接する表面積 [m^2]
    real(8), parameter :: Cp_r = 921.0 !ローターの比熱 [J/kg·K]
    real(8), parameter :: Cp_a = 1006.0 !空気比熱 [J/kg·K]
    real(8), parameter :: h_sa = 120 !熱伝達率 [W/m^2·K]
    real(8), parameter :: k = 0.120 !物質移動係数 [kg/m^2·s(kg/kgDA)]
    real(8), parameter :: L = 2253293.2 !蒸発潜熱 [J/kg]
    real(8), parameter :: rho = 1.206 !空気密度 [kg/m^3]
    real(8), parameter :: Cad_0 = 1.0 !使っていない
    real(8), parameter :: kp = 0.0056 !Polanyi DR 定数 [(cm^3/mol·K)^2]
    real(8), parameter :: Vm = 18.0 !水分子のモル容量 [cm^3/mol]

    !運転条件
    real(8) :: omega = 20.0 !ローターの1時間当たりの回転数 [rph]
    real(8) :: qa_process = 1500.0 !処理側風量 [m^3/h]
    real(8) :: qa_regenerate = 1500.0 !再生側風量 [m^3/h]
    real(8) :: Tai_process = 30.0 !処理側流入温度 [℃]
    real(8) :: Tai_regenerate = 70 !再生側流入温度 [℃]
    real(8) :: X_inProcess = 0.0204 !処理側流入絶対湿度 [kg/kgDA]
    real(8) :: X_inRegenerate = 0.0160 !再生側流入絶対湿度 [Kg/kgDA]

    real(8) :: operation_time = 3600 !稼働時間 [s]
    

    !初期値
    real(8) :: X_in
    real(8) :: X_ao = 0.010998655 !出口絶対湿度 [kg/kgDA]
    real(8) :: X_s = 0.007422933 !シリカゲル表面絶対湿度 [kg/kgDA]
    real(8) :: T_r = 25.0 !ローター温度 [℃]
    real(8) :: T_r_list(DivisionNumber) = 25.0
    real(8) :: X_s_list(DivisionNumber) = 0.007422933
    real(8) :: T_ao = 25.0 !出口温度 [℃]
    real(8) :: P_old(DivisionNumber) = 0.05 !
    real(8) :: P_current, P_step_start 
    real(8) :: P0 = 1.0 
    real(8) :: RHs = 37.72545 !シリカゲル表面相対湿度 [%]
    real(8) :: es  
    real(8) :: Tai 
    real(8) :: time_n !n分割した際の1ステップの時間 [s]
    real(8) :: error = 1.0e-8
    real(8) :: polanyi 
    real(8) :: qa, e, dRHs, dPolanyi, dP,de, dP_content, dXao, P_content
    real(8) :: Tr_estimate,Tr_cal, Xao_cal, Xs_new, Xs_assumed
    real(8) :: X_ao_list(DivisionNumber), T_ao_list(DivisionNumber)
    real(8) :: X_ao_list_temporary(DivisionNumber), T_r_list_temporary(DivisionNumber),&
    X_s_list_temporary(DivisionNumber), T_ao_list_temporary(DivisionNumber), P_old_temporary(DivisionNumber)
    Logical complete !while文終了判定 

    !その他の変数 
    integer :: iter_polanyi, iter !無限ループ対策　カウンター
    integer :: j, i, remainder
    integer :: cal_count 
    real(8) :: Xao_process_sum = 0 
    real(8) :: Xao_regenerate_sum = 0 
    real(8) :: Tao_process_sum = 0
    real(8) :: Tao_regenerate_sum = 0

    time_n = 3600.0 / (dble(DivisionNumber) * omega) !１ステップあたりの時間[s]
    cal_count = int(operation_time / time_n) !稼働時間分のステップ数

    do j = 1, cal_count 
        !確認用print*, "timecount:", j
        !確認用print *, T_r_list
        
        !n分割分それぞれの計算
        do i = 1, DivisionNumber

            Xs_assumed = X_s_list(i)!仮のXs (ひとつ前のローターのシリカゲル表面の絶対湿度)
            !確認用print*, "position:", i
            remainder = mod(i - 1, DivisionNumber)
            !余りが0からDivisionNumberの半分まではqa = qa_process(処理空気)、 それ以外は qa = qa_regenerate(再生空気)
            if (remainder < DivisionNumber / 2) then
                qa = (qa_process / 3600.0) / ( dble(DivisionNumber) / 2.0 )
                X_in = X_inProcess
                Tai = Tai_process
            else
                qa = (qa_regenerate / 3600.0) / ( dble(DivisionNumber) / 2.0 )
                X_in = X_inRegenerate
                Tai = Tai_regenerate
            endif

            complete = .false. !while文の終了判定を初期化

            Tr_estimate = T_r_list(i) !仮のTrを仮定、一タイムステップ前の確定温度を使う
            iter = 0 !ステップが次に進む時に無限ループ対策に使用しているカウンターをリセットする。
            P_step_start = P_old(i)

            do while(.not. complete) 
                !無限ループ対策
                iter = iter + 1
                if (iter > 100) then
                    exit
                    print *, "Convergence failed at time count", j, "position", i
                endif

                Polanyi = 100.0 !Polanyi の初期値を100にしておく
                iter_polanyi = 0

                !仮のXrを決める
                do while(abs(Polanyi) > error)

                    Xao_cal = (rho * qa * X_in + k * A_sa * Xs_assumed)  / (rho * qa + k * A_sa) !仮のXaoを決める
                    !無限ループ対策            
                    iter_polanyi = iter_polanyi + 1
                    if (iter_polanyi > 100) then
                        print *, "Convergence failed in the Xr search. time count", j, "position", i
                        exit 
                    endif

                    e = 28.964 * Xs_assumed * 1013.25 / (Xs_assumed * 28.964 + 18.015) !空気の水蒸気分圧
                    if (e < 0.0d0) then
                        e = 0.0d0
                    endif

                    es = 6.1078 * ( 10.0**(7.5*Tr_estimate / (Tr_estimate+ 237.3) ) ) !飽和水蒸気分圧
                    if (es > 0.0d0) then
                        RHs = 100.0d0 * e / es
                    else           
                        RHs = 0.001d0
                    endif
                    
                    if (RHs > 99.99d0) RHs = 99.99d0
                    if (RHs < 0.01d0) RHs = 0.01d0

                    !c(Xs)を更新したい　→　ニュートン法を使う 解析的に微分をする(誤差をできるだけ小さくしたい)

                    de = (1013.25 * 18.015 / 28.964) / (18.015 / 28.964 + Xs_assumed)**2 !空気の水蒸気分圧をXsで微分
                    dRHs = 100 * de / es!相対湿度をXsで微分
                    P_content = -kp * ( ( (273.16 + Tr_estimate) / Vm )**2.0) * log(100 / RHs)**2.0 !P
                    dP_content = kp*( ( (273.16 + Tr_estimate) / Vm )**2.0) * 2 * log(100 / RHs) / RHs * dRHs
                    dP = P0 * exp(P_content) * dP_content !含湿量をXsで微分
                    dXao = k * A_sa / (rho * qa + k * A_sa)

                    P_current = P0 * exp( -kp * ( ((273.16 + Tr_estimate)/Vm)**2.0 ) * log(100.0/RHs)**2.0 )
                    Polanyi = m_s * (P_current - P_step_start) / time_n - k * A_sa * (Xao_cal - Xs_assumed)
                    
                    dPolanyi = (m_s * dP / time_n ) - k * A_sa * dXao + k * A_sa !Polanyiの微分

                    Xs_new = Xs_assumed - 0.5 * Polanyi / dPolanyi !Xsの更新式 0.5は振れ幅を小さくするための係数
                    !Xsガードレール
                    if (Xs_new < 0.0000001) then
                        Xs_new  = 0.0000001
                    endif
                    Xs_assumed = Xs_new !新しい候補となったXs(Xs_new)をXs_assumedとし、次の計算で使う。
                    
                    
                    !print *, "計算回数",iter_polanyi, "Xs_assumed", Xs_assumed, "polanyi", polanyi,&
                    !"dPolanyi:", dPolanyi, "RHs:", RHs, "e:", e, "es"
                enddo
                X_s = Xs_assumed
                X_ao = (rho * qa * X_in + k * A_sa * Xs_new)  / (rho * qa + k * A_sa)!いったん確定したXsからのXaoを決める

                T_ao = ( (rho * Cp_a * qa * Tai) + (h_sa * A_sa * Tr_estimate) ) / (rho * Cp_a * qa + h_sa * A_sa) !確定のTao
                !次は一番最初の式からTrを算出
                Tr_cal = (m_s * Cp_r / time_n * T_r_list(i) + h_sa * A_sa * T_ao + k * A_sa * L * (X_ao - Xs_new)) / &
                (m_s * Cp_r / time_n + h_sa * A_sa) 

                if (abs(Tr_estimate - Tr_cal) > error) then
                    !確認用write(*, *) "カウント", iter, "誤差",abs(Tr_estimate - Tr_cal), "判定", complete           
                    Tr_estimate = Tr_cal
                else
                    T_r = Tr_cal
                    complete = .true.
                    !確認用write(*, *) "カウント", iter, "誤差",abs(Tr_estimate - Tr_cal), "判定", complete
                endif
            enddo
            P_old_temporary(i) = P_current

            T_ao = ( (rho * Cp_a * qa * Tai) + (h_sa * A_sa * T_r) ) / (rho * Cp_a * qa + h_sa * A_sa) !確定のTao
            !確認用write(*, *)"ステップ:",i, "X_ao:",X_ao, " X_s:", X_s," T_r",T_r," T_ao", T_ao
            X_ao_list_temporary(i) = X_ao
            X_s_list_temporary(i) = X_s
            T_r_list_temporary(i) = T_r
            T_ao_list_temporary(i) = T_ao
        enddo

        !あるタイムステップの出口温度湿度の値 配列のすべての値の平均値をとる
        Xao_process_sum = 0
        Xao_regenerate_sum = 0
        Tao_process_sum = 0
        Tao_regenerate_sum = 0
        do i = 1, DivisionNumber !Xaoの平均
            if(i <= DivisionNumber / 2) then
                Xao_process_sum = Xao_process_sum + X_ao_list_temporary(i)
                Tao_process_sum = Tao_process_sum + T_ao_list_temporary(i)
            else
                Xao_regenerate_sum = Xao_regenerate_sum + X_ao_list_temporary(i)
                Tao_regenerate_sum = Tao_regenerate_sum + T_ao_list_temporary(i)
            endif
        enddo
        Xao_process_sum = Xao_process_sum / (DivisionNumber / 2)
        Tao_process_sum = Tao_process_sum / (DivisionNumber / 2)
        Xao_regenerate_sum = Xao_regenerate_sum / (DivisionNumber / 2)
        Tao_regenerate_sum = Tao_regenerate_sum / (DivisionNumber / 2)

        write(*, *) J * time_n,",", Xao_process_sum, ",", Xao_regenerate_sum, ",", Tao_process_sum, ",", Tao_regenerate_sum

        X_ao_list = cshift(X_ao_list_temporary, shift=-1)
        X_s_list = cshift(X_s_list_temporary, shift=-1)
        T_r_list = cshift(T_r_list_temporary, shift=-1)
        T_ao_list = cshift(T_ao_list_temporary, shift=-1)
        P_old = cshift(P_old_temporary, shift=-1)

    enddo
    stop
end program kaiseki  