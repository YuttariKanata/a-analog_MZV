using Base.MPFR

function azeta(k_index::Vector{Int64},e_index::Vector{Int64},a_abstractfloat::AbstractFloat,iterations::Int,precision::Int64)
    #println(precision*log2(BigFloat(10)))



    #ここでは精度の設定をしている
    #10進数でn桁の精度が必要なら、2進数ではn×log_2(10) bitの精度がいるからね
    setprecision(BigFloat, Int(ceil(precision*log2(BigFloat(10)))) )


    #入力された小数のaをさっき指定した精度でa_valueに入れている
    #これはちょっと遠回りしているけど、これでうまくいくから別にいいよね！
    a_value = parse(BigFloat,string(a_abstractfloat))
    

    #入力された内容の表示をする
    println("\nk: ",k_index," length: ",length(k_index))
    
    println("e: ",e_index," length: ",length(e_index))
    
    println("a: ",a_value)
    #println(a_abstractfloat)
    println("iterations: ",iterations)
    println("precision:  ",precision)



    #もし途中でプログラムを止めたいときは、Ctrl + C を押して止めてね
    #もしjuliaから出たかったら、exit() と入力してね
    println("\nif you want to stop this program, press \"Ctrl + C\".")
    println("if you want to exit from julia,   enter \"exit()\" ")

    println("\n\nstart program!!!\n")




    
    ##########ここからが計算##########

    #rをkの深さ(=eの深さ) pは反復回数
    r = length(k_index)
    p = iterations

    zeta_array = zeros(BigFloat,p)          #zeta関数の値を入れて置く一次元配列
    bino_array = zeros(BigFloat,p+r-1)      #   a類似の値を入れて置く一次元配列

    bino_array[1] = a_value
    a_minus = big(a_value-1.0)

    mul = bino_array[1]
    

    for i in 2:p+r-1
        mul *= (1.0+(a_minus/(i)))
        bino_array[i] = mul
    end

    #Base.print_array(stderr,bino_array);

    

    
    zeta_array[1] = a_value^e_index[1]   #ここで(p,r)=(1,1)を設定。すべてはここから始まる。
    sum1 = zeta_array[1]

    for i in 2:p
        sum1 += ((bino_array[i])^e_index[1]*1/(i^big(k_index[1])))
        zeta_array[i] = sum1
        #println(sum1)
    end

    #Base.print_array(stderr,zeta_array);

    println("setup finish")

    
    co_bino_array = view(bino_array ,1:p)       #bino_arrayの配列の1 ~ pをSubArrayとしておく
    
    
    #Base.print_array(stderr,co_bino_array)
    println()
    

    

    
    for r in 2:length(k_index)

        circshift!(bino_array,-1)   #一つもどるシフト([i] <- [i+1])ただしrotate shiftに注意(circleってそういうこと)
        #[r-1] <- [r] だから [1] <- [r] になっている


        #Base.print_array(stderr,co_bino_array)

    

        #いちいち配列からとってくるより変数に入れといたら早いやろコナミ
        expconst1 = big(e_index[r])
        expconst2 = big(k_index[r])
        p_set = [r:p+r-1;]
        @show expconst2
        @show expconst1


        #これがミソだと思った？残念！これは下準備です！
        zeta_array .*= (co_bino_array.^expconst1)./(p_set.^expconst2)

        println("zeta_array1");
        Base.print_array(stderr,zeta_array);
        println()
        #こっちが本体でした！
        for i in 2:p
            zeta_array[i] += zeta_array[i-1]
        end

        println("zeta_array2");
        Base.print_array(stderr,zeta_array);
        println()
    
        println()
        println("r:",r)
        #Base.print_array(stderr,zeta_array)
        
    end


    println("result:\n",zeta_array[p])

    #println(zeta_array)
    
    
    
    #ここからε算法  頑張るぞい



    epsilon_array1 = zeros(BigFloat,Int(floor(log2(p))))
    epsilon_array2 = zeros(BigFloat,Int(floor(log2(p))))
    @show p
    @show Int(floor(log2(p)))
    for i in 1:Int(floor(log2(p)))
        epsilon_array2[i] = zeta_array[2^i]
        println(epsilon_array2[i])
    end
    
    
    for i in 1:Int(floor(log2(p) / 2.0)-1)
        @show i
        epsilon_subarray = 1 ./(circshift(epsilon_array2,-1) .- epsilon_array2)

        epsilon_array1 = circshift(epsilon_array1,-1) .+ epsilon_subarray

        println("epsilon_array1")
        Base.print_array(stderr,epsilon_array1)
        println()

        pop!(epsilon_array1)
        pop!(epsilon_array2)

        epsilon_subarray = 1 ./(circshift(epsilon_array1,-1) .- epsilon_array1)

        epsilon_array2 = circshift(epsilon_array2,-1) .+ epsilon_subarray

        println("epsilon_array2")
        Base.print_array(stderr,epsilon_array2)
        println()

        pop!(epsilon_array2)
        pop!(epsilon_array1)

    end

    println("epsilon:")
    println(stderr,epsilon_array2[end])



end