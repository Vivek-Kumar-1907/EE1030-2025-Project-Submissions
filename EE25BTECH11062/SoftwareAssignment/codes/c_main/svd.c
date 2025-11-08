#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MAT struct Matrix
#define epsilon 0.00000001

struct Matrix{
    int rows;
    int cols;
    float *data;
};

void matInit(MAT *a, int m, int n){
    a->rows = m;
    a->cols = n;
    a->data = (float *)calloc(m*n, sizeof(float));
}

void matFree(MAT *a){
    free(a->data);
    a->rows = 0;
    a->cols = 0;
}

void cp(MAT *a, MAT *b){
    matFree(a);
    matInit(a, b->rows, b->cols);
    for(int i = 0; i<a->rows*a->cols; i++)
        a->data[i] = b->data[i];
}

void scalarmul(MAT *a, float b, MAT *c){
    matInit(c, a->rows, a->cols);
    for(int i = 0; i<c->rows*c->cols; i++)
        c->data[i] = a->data[i]*b;
}

void mul(MAT *a, MAT *b, MAT *c){
        matInit(c, a->rows, b->cols);
        int m = a->rows;
        int n = a->cols;
        int p = b->cols;
        for(int i = 0; i<m; i++)
            for(int j = 0; j<p; j++)
                for(int k = 0; k<n; k++)
                    c->data[p*i+j] += a->data[n*i+k] * b->data[p*k+j];
}

void T(MAT *a, MAT *b){
    matInit(b, a->cols, a->rows);
    for(int i = 0; i<a->rows; i++)
        for(int j=0; j<a->cols; j++)
            b->data[a->rows*j+i] = a->data[a->cols*i+j];
}

void add(MAT *a, MAT *b, MAT *c){
    matInit(c, a->rows, a->cols);
    for(int i = 0; i<a->rows*a->cols; i++)
        c->data[i] = a->data[i] + b->data[i];
}

void sub(MAT *a, MAT *b, MAT *c){
    matInit(c, a->rows, a->cols);
    for(int i = 0; i<a->rows*a->cols; i++)
        c->data[i] = a->data[i] - b->data[i];
}

void displayMat(MAT *a){
    for(int i = 0; i<a->rows; i++){
        for(int j = 0; j<a->cols; j++)
            printf("%.5f ", a->data[a->cols*i+j]);
        printf("\n");
    }
}

void randomize(MAT *a){
    for(int i = 0; i<a->rows*a->cols; i++)
            a->data[i] = 50*(rand()/(float)RAND_MAX);
}

float norm(MAT *a){
    float t;
    for(int i = 0; i<a->rows*a->cols; i++){
            t += (a->data[i])*(a->data[i]);
    }
    return sqrt(t);
}

void findEig(MAT *a, MAT *v){
    clock_t t1;
    static int t = 0;
    MAT d;
    MAT *c = &d;
    t1 = clock();
    MAT b;
    T(a, &b);
   // mul(&b, a, c);
   // matFree(&b);
    matInit(c, a->cols, a->cols);
    for(int i = 0; i<a->cols; i++)
        for(int j = 0; j<a->cols; j++)
             for(int k = 0; k<a->rows; k++)
                    c->data[a->cols*i+j] += b.data[a->rows*i+k] * b.data[a->rows*j+k];
    matFree(&b);
    int m = c->rows;
    int n = c->cols;
    MAT random, prev, curr;
    matInit(&random, n, 1);
    matInit(&prev, n, 1);
    randomize(&random); 
    scalarmul(&random, 1.0/norm(&random), &curr);
    int itr = 0;
    while(1){
        itr++;
        cp(&prev, &curr);
        MAT d;
        mul(c, &prev, &d);
        cp(&curr, &d);
        matFree(&d);
        MAT e;
        scalarmul(&curr, 1.0/norm(&curr), &e);
        cp(&curr, &e);
        matFree(&e);
        MAT curr_T, dot;
        T(&curr, &curr_T);
        mul(&curr_T, &prev, &dot);
        if(abs(dot.data[0]) > 1-epsilon || itr >= 2000){
            matFree(&prev);
            matFree(&dot);
            matFree(&curr_T);
            matFree(&random);
            matInit(v, curr.rows, curr.cols);
            cp(v, &curr);
            t++;
            t1 = clock()-t1;
            printf("Found Singular Value No.%d in %d iterations after %lf seconds\n", t, itr, (double)t1/CLOCKS_PER_SEC);
            return;
        }
    }
}

void svd(MAT *a, int k, MAT *res){
    float arr[k];
    MAT P;
    MAT *b = &P;
    matInit(b, a->rows, a->cols);
    cp( b, a);
    for(int i = 0; i<k; i++){
        MAT u, v;
        matInit(&u, a->rows, 1);
        findEig(b, &v);
        mul(a, &v, &u);
        MAT u_n, v_t;
        float sigma = norm(&u);
        scalarmul(&u, 1.0/sigma, &u_n);
        cp(&u, &u_n);
        T(&v, &v_t);
        MAT sing_n;
        MAT sing;
        mul(&u, &v_t, &sing_n);
        scalarmul(&sing_n, sigma, &sing);
        MAT bmsing;
        sub(b, &sing, &bmsing);
        cp(b, &bmsing);
        matFree(&bmsing);
        MAT respsing;
        add(res, &sing, &respsing);
        cp(res, &respsing);
        matFree(&respsing);
        matFree(&u_n);
        matFree(&v_t);
        matFree(&u);
        matFree(&sing_n);
        matFree(&sing); 
    }
}

int checkwhitespace(char character){
    if(character == ' ' || character == '\t' || character == '\n' || character == '\v'|| character == '\f'|| character == '\r'){
        return 1;
    }
    else
        return 0;
}

int power(int a, int b){
    return (b == 0)? 1 : a* power(a, b-1);
}

int convertstring(char *num){
    int len = 0;
    int n = 0;
    while(*(num+len) != '\0')
        len++;
    for(int i = 0; i<len; i++){
        n += (*(num+i) - '0')*power(10, len-i-1);
    }
    return n;

}

int get_ascii(char *curr, FILE *fptr){
    while(checkwhitespace(*curr))
        fread(curr, sizeof(char), 1, fptr);
    char *val = (char *)calloc(5, sizeof(char));
    int idx = 0;

    while(!checkwhitespace(*curr)){
        *(val+idx) = *curr;
        idx++;
        fread(curr, sizeof(char), 1 ,fptr);
    }
    *(val+idx) = '\0';
    int n = convertstring(val);
    free(val);
    return n;
}

void process_metadata(FILE *fptr, int *w, int *h, int *max){
    char *magic = (char *)malloc(2*sizeof(char));
    int idx = 0;
    fread(magic, sizeof(char), 2, fptr);
    idx += 2;
    printf("%s\n", magic);

    char *curr = (char *)malloc(1*sizeof(char));
    fread(curr, sizeof(char), 1, fptr);

    *w = get_ascii(curr, fptr);
    printf("Width of image: %d\n", *w);

    *h = get_ascii(curr, fptr);
    printf("Height of image: %d\n", *h);

    *max = get_ascii(curr, fptr);
    printf("Maximum value of a pixel: %d\n", *max);

    if(*curr == '#'){
        while(*curr != '\n')
            fread(curr, sizeof(char), 1, fptr);
    }

    fread(curr, sizeof(char), 1, fptr);

    free(curr);
    free(magic);
}

void process_data(FILE *fptr, int w, int h, int max, int **data, MAT *out){
    int bitsize = (max<256) ? 1 : 2;
    matInit(out, h, w);
    for(int i = 0; i<h; i++){
        for(int j = 0; j<w; j++){
            fread(*(data+i)+j, bitsize, 1, fptr);
            out->data[w*i+j] = *(*(data+i)+j);
        }
    }
}

void free_data(int **data, int h){
    for(int i = 0; i<h; i++)
        free(*(data+i));
    free(data);
}

int main(){
    char input_filename[100], input_name[100];
    char output_filename[100], output_name[100];
    printf("Enter complete input file name:- \n");
    scanf("%s", input_filename);
    int len_inp_filename = 0, len_out_filename = 0;
    while(input_filename[len_inp_filename] != '\0') len_inp_filename++; 
    for(int i = 0; i<len_inp_filename; i++){
        if(input_filename[i] == '.'){
            input_name[i] = '.';
            input_name[i+1] = 'p';
            input_name[i+2] = 'g';
            input_name[i+3] = 'm';
            input_name[i+4] = '\0';
            break;
        }
        input_name[i] = input_filename[i];
    }
    char conv_cmd1[300];
    sprintf(conv_cmd1, "convert %s %s", input_filename, input_name);
    printf("%s\n", conv_cmd1);
    system(conv_cmd1);
    printf("Enter complete output file name: \n");
    scanf("%s", output_filename);
    while(output_filename[len_out_filename] != '\0') len_out_filename++;
    printf("%d\n", len_out_filename);
    for(int i = 0; i<len_out_filename; i++){
        if(output_filename[i] == '.'){
            output_name[i] = '.';
            output_name[i+1] = 'p';
            output_name[i+2] = 'g';
            output_name[i+3] = 'm';
            output_name[i+4] = '\0';
            break;
        }
        output_name[i] = output_filename[i];
    }
    printf("%s\n", output_name);
    FILE *fptr = fopen(input_name, "rb");
    int mode = 0;
    if(fptr == NULL)
        printf("Image does not exist/not found\n");
    else{
        int w, h, max, k;
        process_metadata(fptr, &w, &h, &max);
        int **data = (int **)calloc(h, sizeof(int *));
        for(int i = 0; i<h; i++)
            *(data+i) = (int *)calloc(w, sizeof(int));
        MAT out, res, restricted_res, diff;
        process_data(fptr, w, h, max, data, &out);
        matInit(&res, out.rows, out.cols);
        printf("Enter a k value between 1 and %d:\n", (w-h>0)?h:w);
        scanf("%d", &k);
        svd(&out, k, &res);
        matInit(&restricted_res, res.rows, res.cols);
        for(int i = 0; i<h; i++){
            for(int j = 0; j<w; j++){
                int a = (int)res.data[w*i+j];
                a = (a>0)?a:0;
                a = (a<max)?a:max;
                restricted_res.data[w*i+j] = a;
            }
        }
        sub(&out, &res, &diff);
        float frobenius_norm = norm(&diff);
        printf("Frobenius norm ||A-A_k|| = %f\n", frobenius_norm);
        fclose(fptr);
        FILE *fptr2 = fopen(output_name, "wb+");
        fprintf(fptr2, "P5 %d %d %d\n", w, h, max);
        int bitsize = (max<256) ? 1 : 2;
        for(int i = 0; i<h; i++){
            for(int j = 0; j<w; j++){
                int a = (int)restricted_res.data[w*i+j];
                fwrite(&a, bitsize, 1, fptr2);
            }
        }
        //char conv_cmd_2[300];
        //sprintf(conv_cmd_2, "convert %s %s", output_name, output_filename);
        MAT d;
        //system(conv_cmd_2);
        matFree(&out);
        matFree(&res);
        matFree(&restricted_res);
        matFree(&diff);
        free_data(data, h);
        fclose(fptr2);
    }
    return 0;
}
