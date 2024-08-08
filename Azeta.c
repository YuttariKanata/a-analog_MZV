#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <string.h>
#include <stdlib.h>
//#include <sys/resource.h>
#include <emscripten.h>



/*
input_ptr_Com
[0]:精度(decimal)
[1]:上限
[2]/[3]:a
[4]-[13]:kインデックス
[14]-[23]:eインデックス
[24]-[33]:eインデックスの正負 (負なら1 正なら0)

output_ptr
[0]:普通の計算結果の指数
[1]:ε算法での指数
[2]-[n] :普通の計算結果
[n+1] : 11
[n]-[m] :ε算法での計算結果
*/
long azeta(unsigned long* input_ptr_Com,unsigned long* output_ptr){

    //fprintf(stderr,"start\n");
    
    long exp_ptr=987;//てきとうな数　あとで変わるしいいじゃろ
    

    //精度(2進数)
    unsigned long prec10 = input_ptr_Com[0];
    //fprintf(stderr,"prec:%lu\n",prec);

    //精度(10進数)
    unsigned long prec = (unsigned long)(roundl((long double)(prec10 * 3.32192809489L))+2UL); // 念押しの+2
    //fprintf(stderr,"prec10:%lu\n",prec10);

    
    //次のmallocに入れるためだけの変数 数字の部分+符号1bit+終端記号1bit
    size_t buffer_size = prec10+2;

    //出力用の文字列1
    char *str_value1 = (char *)malloc(buffer_size * sizeof(char));
    //もし帰ってきたのがNULLなら失敗してる
    if (str_value1 == NULL) {
        //fprintf(stderr, "メモリの割り当てに失敗しました。\n");
        return -987;
    }else{
        //fprintf(stderr, "メモリの割り当てに成功しました。\n");
    }
    //出力用の文字列2
    char *str_value2 = (char *)malloc(buffer_size * sizeof(char));
    //もし帰ってきたのがNULLなら失敗してる
    if (str_value2 == NULL) {
        //fprintf(stderr, "メモリの割り当てに失敗しました。\n");
        return -654;
    }else{
        //fprintf(stderr, "メモリの割り当てに成功しました。\n");
    }

    //gmpの小数の精度設定
    mpf_set_default_prec(prec);



    //多目的変数(主にfor文の回すやつ)
    unsigned long i=0;

    //インデックスの深さ
    unsigned long r=0;
    for ( i = 4; i < 14; i++)
    {
        if(input_ptr_Com[i] != 0){
            r++;
        }
    }
    //fprintf(stderr,"r:%lu\n",r);


    //反復回数の上限
    unsigned long p = input_ptr_Com[1];
    fprintf(stderr,"p:%lu\n",p);


    //ε算法ができる回数を求める
    unsigned long length_buf = (unsigned long)floor(log2((double)p));
    fprintf(stderr,"length_buf:%lu\n",length_buf);



    //########## mpz_t,mpf_t型変数の初期化 ##########
    mpf_t a;
    mpf_t a_minus;
    mpf_t bino_param;
    mpf_t buf1;
    mpf_t buf3;
    mpf_t sum1;
    mpf_inits(a,a_minus,bino_param,buf1,buf3,sum1,NULL);

    mpz_t buf2;
    mpz_init(buf2);

    //0~pの配列 MZVの途中計算を入れるが、本来は1~p(julia)だった。
    mpf_t zeta_array[p+1];          //0~pの配列を作る
    for ( i = 0; i < p+1; i++){
        mpf_init(zeta_array[i]);    //初期化
    }
    //a類似の値を入れて置く配列。本来は1~p+r-1(julia)だった。
    mpf_t bino_array[r];            //0~p+r-1の配列を作る
    for ( i = 0; i < p+r; i++)
    {
        mpf_init(bino_array[i]);
    }
    
    mpf_t epsilon_array1[length_buf+1];
    for ( i = 0; i <= length_buf; i++){
        mpf_init(epsilon_array1[i]);    //初期化
    }
    
    mpf_t epsilon_array2[length_buf+1];
    for ( i = 0; i <= length_buf; i++)
    {
        mpf_init(epsilon_array2[i]);    //初期化
    }

    mpf_t epsilon_subarray[length_buf+1];
    for ( i = 0; i <= length_buf; i++)
    {
        mpf_init(epsilon_subarray[i]);  //初期化
    }
    //########## 終わり ##########


    //aの値(どうするか悩んだ末、分数)
    double a_double = (double)input_ptr_Com[2]/(double)input_ptr_Com[3];

    //fprintf(stderr,"a_double:%f\n",a_double);


    mpf_set_d(a,a_double);//mpfとして入れて置く
    //gmp_fprintf(stderr,"a:%.*Ff\n",50,a);



    //#################ここから計算#######################################


    mpf_set(bino_array[1],a);
    //gmp_fprintf(stderr,"bino_array[1]:%.*Ff\n",50,bino_array[1]);


    mpf_set_d(a_minus,a_double-1.0);
    //gmp_fprintf(stderr,"a_minus:%.*Ff\n",50,a_minus);

    
    mpf_set(bino_param,a);
    //gmp_fprintf(stderr,"mul:%.*Ff\n",50,bino_param);

    
    for ( i = 2; i <= r; i++){
        mpf_div_ui(buf1,a_minus,i);
        mpf_add_ui(buf1,buf1,1UL);
        mpf_mul(bino_param,bino_param,buf1);
        mpf_set(bino_array[i],bino_param); //bino_arrayの設定
        ////gmp_fprintf(stderr,"bino_array[%10lu]:%.*Ff\n",i,prec10,bino_param);
    }
    

    mpf_pow_ui(zeta_array[1],a,input_ptr_Com[14]);  //ここで(p,r)=(1,1)を設定。すべてはここから始まる。
    /*  zeta_array[1] = a^input_ptr_Com[14] 
    本来はa^eなので正負を考えないといけない
    そのためのif文が今後もいっぱい出てくる     */
    if(input_ptr_Com[24]){
        mpf_ui_div(zeta_array[1],1UL,zeta_array[1]);
    }
    ////gmp_fprintf(stderr,"zeta_array[%10lu]:%.*Ff\n",1UL,50,zeta_array[1]);

    
    mpf_set(sum1,zeta_array[1]);

    mpf_set(bino_param,a);
    for ( i = 2; i <=p; i++){
        mpf_div_ui(buf1,a_minus,i);
        mpf_add_ui(buf1,buf1,1UL);
        mpf_mul(bino_param,bino_param,buf1);
        mpf_pow_ui(buf1,bino_param,input_ptr_Com[14]);
        if(input_ptr_Com[24]){
            mpf_ui_div(zeta_array[1],1UL,zeta_array[1]);
        }
        mpz_ui_pow_ui(buf2,i,input_ptr_Com[4]);
        mpf_set_z(buf3,buf2);
        mpf_div(buf1,buf1,buf3);
        mpf_add(sum1,sum1,buf1);
        mpf_set(zeta_array[i],sum1);
        ////gmp_fprintf(stderr,"zeta_array[%10lu]:%.*Ff\n",i,50,sum1);
    }

    //co_array = array[1~p]
    
    unsigned long r_const = r;
    //fprintf(stderr,"r_const:%lu\n",r_const);
    unsigned long exp_const1;
    unsigned long exp_const2;
    unsigned long exp_flag;
    for ( r = 2; r <= r_const; r++)
    {
        //co_array = array[r~p+r-1] (p個)
        exp_const1 = input_ptr_Com[r+13];
        exp_const2 = input_ptr_Com[r+3];
        exp_flag   = input_ptr_Com[r+23];

        /*
        ここは計算のメインとなる部分の一つであり、
        極力計算量を落としたいためにこのようなネストを行う。
        ネストするのはわかりやすいようにしただけなんだけどね
        rは小さい(最大10)だが、pが馬鹿でかいので計算量を削ったほうがいい。
        */
        //fprintf(stderr,"exp_const2:%lu exp_const1:%lu exp_flag:%lu\n",exp_const2,exp_const1,exp_flag);

        
        if (exp_const1){//eインデックスが0でない場合
            mpf_set(bino_param,bino_array[r-1]);
            if (exp_flag){//eインデックスが負の時
                for ( i = 1; i <= p; i++)
                {
                    mpf_div_ui(buf1,a_minus,r-1+i);
                    mpf_add_ui(buf1,buf1,1UL);
                    mpf_mul(bino_param,bino_param,buf1);
                    mpf_pow_ui(buf1,bino_param,exp_const1);
                    mpf_ui_div(buf1,1UL,buf1);  //このひっくり返す作業がいる
                    mpz_ui_pow_ui(buf2,r-1+i,exp_const2);
                    mpf_set_z(buf3,buf2);
                    mpf_div(buf1,buf1,buf3);
                    mpf_mul(zeta_array[i],zeta_array[i],buf1);
                }
            }else{//eインデックスが正の時
                for ( i = 1; i <= p; i++)
                {
                    mpf_div_ui(buf1,a_minus,r-1+i);
                    mpf_add_ui(buf1,buf1,1UL);
                    mpf_mul(bino_param,bino_param,buf1);
                    mpf_pow_ui(buf1,bino_param,exp_const1);
                    //mpf_ui_div(buf1,1UL,buf1);  //このひっくり返す作業がいらない
                    mpz_ui_pow_ui(buf2,i,exp_const2);
                    mpf_set_z(buf3,buf2);
                    mpf_div(buf1,buf1,buf3);
                    mpf_mul(zeta_array[i],zeta_array[i],buf1);
                }
            }
        }else{//eインデックスが0の時
            for ( i = 1; i <= p; i++)
            {
                //mpf_pow_ui(buf1,bino_array[r-1+i],exp_const1);    //ここら辺の計算がいらないから楽
                //mpf_ui_div(buf1,1UL,buf1);
                mpz_ui_pow_ui(buf2,i,exp_const2);
                mpf_set_z(buf3,buf2);
                mpf_ui_div(buf1,1UL,buf3);
                mpf_mul(zeta_array[i],zeta_array[i],buf1);
            }
        }

        //fprintf(stderr,"zeta_array1\n");
        //for ( i = 1; i <= p; i++)
        //{
        //    gmp_fprintf(stderr,"    %.*Ff\n",prec10,zeta_array[i]);
        //}
        
        

        for ( i = 2; i <= p; i++)
        {
            mpf_add(zeta_array[i],zeta_array[i],zeta_array[i-1]);
        }
        //fprintf(stderr,"zeta_array2\n");
        //for ( i = 1; i <= p; i++)
        //{
        //    gmp_fprintf(stderr,"    %.*Ff\n",prec10,zeta_array[i]);
        //}
        //fprintf(stderr,"r_now:%lu\n",r);
    }

    gmp_fprintf(stderr,"result1:%.*Ff\n",prec10,zeta_array[p]);

    //gmp_fprintf(stderr,"%.*Ff\n",30,zeta_array[p]);//結果を表示
    mpf_get_str(str_value1, &exp_ptr, 10, prec10, zeta_array[p]);
    if (exp_ptr>0){
        output_ptr[0]=(unsigned long)(2UL*abs(exp_ptr));
    }
    else if (exp_ptr<0){
        output_ptr[0] = (unsigned long)(2UL*abs(exp_ptr)+1);
    }
    else{
        output_ptr[0] = 0UL;
    }
    
    size_t length_sigma = strlen(str_value1);    //prec10+2
    fprintf(stderr,"length_sigma:%ld\n",length_sigma);
    for ( i = 0; i < prec10; ++i) {       //0~prec10+1   i+2は2~prec10+3
        if(str_value1[i] == '.'){
            output_ptr[i+2] = 10;
            //もし小数点だったら10を入れる。
            continue;
        }else if (str_value1[i] == 0)
        {
            output_ptr[i+2] = 0;
        
        }else{
            output_ptr[i+2] = (unsigned long)(str_value1[i]-48);
            //ASCIIで数字nはn+48とあらわされるから48を引くとその数字が手に入る。
        }
        //fprintf(stderr,"i+2:%ld[%ld]\n",i+2,output_ptr[i+2]);
    }

    output_ptr[length_sigma+2] = 11UL;      //prec10+4
    //fprintf(stderr,"length_sigma+2:%ld\n",length_sigma+2);
    //区切りを表すための数字として11 (10は小数点として使う)
    
    //ここからε算法(2冪)




    

    unsigned long count = (unsigned long)(length_buf/2.0-1);
    //fprintf(stderr,"length_buf:%lu\n",length_buf);
    //fprintf(stderr,"count:%lu\n",count);
    unsigned long t=1;
    for ( i = 1; i <= length_buf; i++)
    {
        t *= 2;
        mpf_set(epsilon_array2[i],zeta_array[t]);
        ////gmp_fprintf(stderr,"    epsilon2_%d:%.*Ff\n",i,prec10,epsilon_array2[i]);
    }

    unsigned long array_dec = 0;
    unsigned long j=1;

    for ( i = 1; i <= count; i++)
    {
        //fprintf(stderr,"epsilon1_i:%lu\n",i);
        for ( j = 1; j <= length_buf-1-array_dec; j++)
        {
            mpf_sub(buf1,epsilon_array2[j+1],epsilon_array2[j]);
            mpf_ui_div(buf1,1UL,buf1);
            mpf_add(epsilon_array1[j],epsilon_array1[j+1],buf1);
            ////gmp_fprintf(stderr,"    epsilon1_%d:%.*Ff\n",j,50,epsilon_array1[j]);
        }
        array_dec++;
        //fprintf(stderr,"epsilon2_i:%lu\n",i);
        for ( j = 1; j <= length_buf-1-array_dec; j++)
        {
            mpf_sub(buf1,epsilon_array1[j+1],epsilon_array1[j]);
            mpf_ui_div(buf1,1UL,buf1);
            mpf_add(epsilon_array2[j],epsilon_array2[j+1],buf1);
            ////gmp_fprintf(stderr,"    epsilon2_%d:%.*Ff\n",j,50,epsilon_array2[j]);
        }
        array_dec++;
    }
    
    //gmp_fprintf(stderr,"epsilon2_1:%.*Ff\n",50,epsilon_array2[1]);
    //gmp_fprintf(stderr,"epsilon2_2:%.*Ff\n",50,epsilon_array2[2]);
    gmp_fprintf(stderr,"epsilon2_%d:%.*Ff\n",length_buf-array_dec,prec10,epsilon_array2[length_buf-array_dec]);
    //gmp_fprintf(stderr,"epsilon1_1:%.*Ff\n",50,epsilon_array1[1]);
    //gmp_fprintf(stderr,"epsilon1_2:%.*Ff\n",50,epsilon_array1[2]);

    //########## mpz_t,mpf_t型変数の解放 ##########

    //gmp_fprintf(stderr,"%.*Ff\n",30,epsilon_array2[length_buf-array_dec]);//結果を表示
    mpf_get_str(str_value2, &exp_ptr, 10, prec10, epsilon_array2[length_buf-array_dec]);
    if (exp_ptr>0){
        output_ptr[1]=(unsigned long)(2UL*abs(exp_ptr));
    }
    else if (exp_ptr<0){
        output_ptr[1] = (unsigned long)(2UL*abs(exp_ptr)+1);
    }
    else{
        output_ptr[1] = 0UL;
    }
    size_t length_epsilon = strlen(str_value2);  //prec10+2
    fprintf(stderr,"length_epsilon:%ld\n",length_epsilon);
    for ( i = 0; i < prec10; ++i) {     //0~prec10+1    iはprec10+5~2*prec10+6
        if(str_value2[i] == '.'){
            output_ptr[i+length_sigma+3] = 10UL;
            continue;
        }else if (str_value2[i] == 0)
        {
            output_ptr[i+length_sigma+3] = 0;
        
        }else{
            output_ptr[i+length_sigma+3] = (unsigned long)(str_value2[i]-48);
        }
        //fprintf(stderr,"i+length_sigma+3:%ld[%ld]\n",i+length_sigma+3,output_ptr[i+length_sigma+3]);
    }

    mpf_clears(a,a_minus,bino_param,buf1,buf3,sum1,NULL);

    mpz_clear(buf2);
   
    for ( i = 0; i < p+1; i++){
        mpf_clear(zeta_array[i]);    //解放
    }

    
    for ( i = 0; i < r_const; i++)
    {
        mpf_clear(bino_array[i]);
    }
    
    for ( i = 0; i <= length_buf; i++){
        mpf_clear(epsilon_array1[i]);    //解放
    }

    for ( i = 0; i <= length_buf; i++)
    {
        mpf_clear(epsilon_array2[i]);    //解放
    }

    for ( i = 0; i <= length_buf; i++)
    {
        mpf_clear(epsilon_subarray[i]);  //解放
    }

    free(str_value1);
    free(str_value2);
    //########## 終わり ##########

    return 0L;

}

/*
input_ptr_Com
[0]:精度(bit)
[1]:上限
[2]/[3]:a
[4]-[13]:kインデックス
[14]-[23]:eインデックス
*/

/*
int main(int argc, char const *argv[])
{
    struct rlimit rl;
    
    // 現在のリソース制限を取得
    if (getrlimit(RLIMIT_STACK, &rl) == -1) {
        perror("getrlimit");
        return 1;
    }
    
    printf("Current stack size limit: %lu\n", rl.rlim_cur);
    
    // スタックサイズを 32 MiB に設定
    rl.rlim_cur = 32 * 1024 * 1024;
    rl.rlim_max = 32 * 1024 * 1024;
    
    if (setrlimit(RLIMIT_STACK, &rl) == -1) {
        perror("setrlimit");
        return 1;
    }
    
    printf("Stack size limit successfully updated.\n");


    unsigned long i = 3906;
    unsigned long input_array[34] = {200UL,i,1UL,2UL  ,  2UL,2UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL  ,  0UL,1UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL , 0UL,1UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL};
    unsigned long output_array[200];
    azeta(input_array,output_array);
    fprintf(stderr,"size:%ld\n",sizeof(output_array)/sizeof(output_array[0]));
    unsigned int size_length = sizeof(output_array)/sizeof(output_array[0]);
    for ( i = 0; i < size_length; i++)
    {
        fprintf(stderr,"output_array[%ld]:%lu\n",i,output_array[i]);
    }
    

    return 0;
}
*/

