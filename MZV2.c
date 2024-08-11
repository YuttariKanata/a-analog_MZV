#include <stdio.h>
#include <math.h>
#include <gmp.h>
#include <string.h>
#include <stdlib.h>
#include <emscripten.h>


long Mzeta2(unsigned long* input_ptr_Com,unsigned long* output_ptr){
    long exp_ptr = 987;
    unsigned long prec = input_ptr_Com[12];//精度(bit)
    unsigned long prec10 = (unsigned long)round((double)(prec * 0.30102999))+10UL;
    size_t buffer_size = prec;
    char *str_value = (char *)malloc(buffer_size * sizeof(char));
    if (str_value == NULL) {
        fprintf(stderr, "メモリの割り当てに失敗しました。\n");
        return -987;
    }
    

    mpf_set_default_prec(prec);
    unsigned long p = input_ptr_Com[10];
    long u = 0;
    unsigned long flag = 0;
    for(u = 0; u < 10; u++){
        if (input_ptr_Com[u] != 0){
            flag++;
        }
    }
    unsigned long r = flag;

    unsigned long a = input_ptr_Com[11];
    if(r == 1){
        unsigned long i1=0;
        unsigned long s1 = (unsigned long)(input_ptr_Com[0]);
        
        mpz_t y1;
        mpz_inits(y1,NULL);
        mpf_t yf1;
        mpf_t sum1;
        mpf_inits(yf1,sum1,NULL);
        

        for(i1 = a+1; i1 <= p; i1++){
            mpz_ui_pow_ui(y1,i1,s1);    //i1^s1
            mpf_set_z(yf1,y1);
            mpf_ui_div(yf1,1UL,yf1);   //i1^(-s1)
            mpf_add(sum1,sum1,yf1);
        }

        gmp_fprintf(stderr,"part1:%.*Ff\n",30,sum1);
        mpf_get_str(str_value, &exp_ptr, 10, prec, sum1);
        size_t length = strlen(str_value);
        for (size_t i = 0; i < length; ++i) {
            if(str_value[i] == '.'){
                output_ptr[i] = 10;
                continue;
            }else{
                output_ptr[i] = (unsigned long)(str_value[i]-48);
            }
        }
        free(str_value);

        mpz_clears(y1,NULL);
        mpf_clears(yf1,sum1,NULL);

        //return sum;               ←###ここを考える必要がある###


        return exp_ptr;
    }

    

    mpf_t array1[p-a-1];
    for(u = 0; u < p-a-1; u++){
        mpf_init(array1[u]);
    }

    
    
    mpf_t firstf;
    mpf_t add1f;
    mpf_t add2f;

    mpf_t beff;
    mpf_inits(firstf,add1f,add2f,beff,NULL);

    mpz_t add1z;
    mpz_t add2z;

    mpz_t befz;
    mpz_inits(add1z,add2z,befz,NULL);

    mpz_ui_pow_ui(add1z,p+r-2,input_ptr_Com[r-1]);
    mpz_ui_pow_ui(add2z,p+r-1,input_ptr_Com[r-1]);
    mpf_set_z(add1f,add1z);
    mpf_set_z(add2f,add2z);
    mpf_ui_div(add1f,1UL,add1f);
    mpf_ui_div(add2f,1UL,add2f);
    mpf_add(firstf,add1f,add2f);
    mpf_set(array1[p-a-2],firstf);

    //fprintf(stderr,"p:%lu a:%lu\n",p,a);


    for(u = p-a-3; u >= 0; --u){
        
        //fprintf(stderr,"hey%lu\n",u);
        mpz_ui_pow_ui(add1z,u+a+r,input_ptr_Com[r-1]);
        mpf_set_z(add1f,add1z);
        mpf_ui_div(add1f,1UL,add1f);

        mpf_add(firstf,firstf,add1f);

        mpf_set(array1[u],firstf);
        //gmp_fprintf(stderr,"array1:%.*Ff\n",20,array1[u]);
        //fprintf(stderr,"u:%lu",u);
    }
    //fprintf(stderr,"hey r:%lu\n",r);
    //fprintf(stderr,"hey r:%lu\n",input_ptr_Com[r-1]);

    mpz_ui_pow_ui(befz,p+r-1,input_ptr_Com[r-1]);
    mpf_set_z(beff,befz);
    mpf_ui_div(beff,1UL,beff);

    mpf_t coc_arrf[p-a-1];
    for(u = 0; u < p-a-1; u++){
        mpf_init(coc_arrf[u]);
    }

    


    mpf_t coc_beforef;
    mpf_t df;
    mpf_t sumf;
    mpf_t beforef;
    mpf_t buff;
    mpf_inits(coc_beforef,df,sumf,beforef,buff,NULL);

    mpz_t dz;
    mpz_t bufz;
    mpz_inits(dz,bufz,NULL);

    mpf_t outarrf[p-a-1];
    for(u = 0; u < p-a-1; u++){
        mpf_init(outarrf[u]);
    }

    

    long an = 0;
    unsigned long i = r-2;
    unsigned long msi1 = input_ptr_Com[r-2];
    mpf_set(beforef,beff);


    an = p-a-2;

    mpz_ui_pow_ui(dz,an+a+i+1,msi1);
    mpf_set_z(df,dz);
    mpf_ui_div(df,1UL,df);

    mpf_mul(df,df,array1[an]);

    

    mpz_ui_pow_ui(bufz,p+i,msi1);
    mpf_set_z(buff,bufz);
    mpf_ui_div(buff,1UL,buff);

    mpf_mul(buff,buff,beforef);


    mpf_add(sumf,df,buff);

    mpf_set(outarrf[an],sumf);

    //fprintf(stderr,"hey1\n");
    //fprintf(stderr,"p:%lu a:%lu an:%lu i:%lu msi1:%lu u:%lu r:%lu\n",p,a,an,i,msi1,u,r);
    //fprintf(stderr,"p-a-3:%lu\n",p-a-3);
    

    for(an = p-a-3; an >= 0; an--){

        //fprintf(stderr,"hey an:%lu\n",an);

        mpz_ui_pow_ui(dz,an+a+i+1,msi1);
        //fprintf(stderr,"%lu^%lu",an+a+i+1,msi1);
        mpf_set_z(df,dz);
        mpf_ui_div(df,1UL,df);
        //gmp_fprintf(stderr,"f:%.*Ff\n",20,df);

        mpf_mul(df,df,array1[an]);
        //gmp_fprintf(stderr,"arr:%.*Ff\n",20,array1[an]);
        //gmp_fprintf(stderr,"df:%.*Ff\n",20,df);


        mpf_add(sumf,sumf,df);

        mpf_set(outarrf[an],sumf);

    }

    

    mpf_set(coc_beforef,buff);

    ////////////////////////////////////////////////

    //fprintf(stderr,"hey2\n");
    for(u = r-3; u >= 0; u--){
        if(u == -1){
            break;
        }
        //fprintf(stderr,"hey3\n");
        msi1 = input_ptr_Com[u];
        i = u;
        mpf_set(beforef,coc_beforef);

        an = p-a-2;
        
        //fprintf(stderr,"p:%lu a:%lu an:%lu i:%lu msi1:%lu u:%lu r:%lu\n",p,a,an,i,msi1,u,r);
        mpz_ui_pow_ui(dz,an+a+i+1,msi1);
        mpf_set_z(df,dz);
        mpf_ui_div(df,1UL,df);

        //fprintf(stderr,"hey4\n");
        mpf_mul(df,df,outarrf[an]);


        mpz_ui_pow_ui(bufz,p+i,msi1);
        mpf_set_z(buff,bufz);
        mpf_ui_div(buff,1UL,buff);

        mpf_mul(buff,buff,beforef);

        mpf_add(sumf,df,buff);

        mpf_set(outarrf[an],sumf);

        for(an = p-a-3; an >= 0; an--){

            mpz_ui_pow_ui(dz,an+a+i+1,msi1);
            mpf_set_z(df,dz);
            mpf_ui_div(df,1UL,df);

            mpf_mul(df,df,outarrf[an]);


            mpf_add(sumf,sumf,df);

            mpf_set(outarrf[an],sumf);
        }
        mpf_set(coc_beforef,buff);

    }

    gmp_fprintf(stderr,"%.*Ff\n",30,outarrf[0]);//結果を表示
    mpf_get_str(str_value, &exp_ptr, 10, prec, outarrf[0]);
    size_t length = strlen(str_value);
    for (size_t i = 0; i < length; ++i) {
        if(str_value[i] == '.'){
            output_ptr[i] = 10;
            continue;
        }else{
            output_ptr[i] = (unsigned long)(str_value[i]-48);
        }
    }
    free(str_value);

    
    for(u = 0; u < p-a-1; u++){
        mpf_clear(array1[u]);
    }

    mpf_clears(firstf,add1f,add2f,beff,NULL);
    mpz_clears(add1z,add2z,befz,NULL);

    for(u = 0; u < p-a-1; u++){
        mpf_clear(coc_arrf[u]);
    }

    mpf_clears(coc_beforef,df,sumf,beforef,buff,NULL);

    mpz_clears(dz,bufz,NULL);

    
    for(u = 0; u < p-a-1; u++){
        mpf_clear(outarrf[u]);
    }

    return exp_ptr;

}




