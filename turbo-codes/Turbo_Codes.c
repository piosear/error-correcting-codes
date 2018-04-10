#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define INFO_LENGTH 1000  // information length
#define NU 3              // number of register
#define TOTAL_LENGTH (INFO_LENGTH + NU)
#define CODE_LENGTH (3 * TOTAL_LENGTH + NU)  // codeword length
#define CODE_RATE ((double)INFO_LENGTH / CODE_LENGTH)
#define NUM_OF_STATE (1 << NU)
#define ITERATE 8  // number of iteration
#define INF 1E+20

#define TRIAL 1000
#define LOOP_MAX 100
#define SNR_MAX 1.0
#define SNR_MIN 0.0
#define SNR_WIDTH 0.2

#define BER_END_CONDITION 1E-6
#define ERROR_END_CONDITION 500
#define ERROR_MIN 10

#define ROW 64
#define COLUMN (INFO_LENGTH / ROW)
#define SPREAD ((int)sqrt(INFO_LENGTH / 2.0))  // spread of interleaver


void generatorPolynomial();
void trelisInit();
void turboEncoder(int *infoArray, int *codedArray);
void recursiveSystematicConvolution(int *infoArray, int *output, bool interleave, int *termination);
void llrCalculator(double complex *transmittedArray, double *llr, int n);
void turboDecoder(double *llr, int *estimatedArray);
void logMAP_decoder(double output1[][2], double output2[][2], double *Lambda, double *Lambda_e,
                    bool interleave);
double Max_fc(double delta1, double delta2);

void bpsk(int *infoArray, double complex *transmittedArray, int num);
void noise(double complex *transmittedArray, int num);
int errorCount(int *infoArray, int *estimatedArray, int num);

void srandomInterleaver();
bool sorting1(int start);
bool sorting2(int start);

double var;
double sigma;

int **interleaveIndex;
int *deinterleaveIndex;
int **forward_path;
int **backward_path;
int **output_conv;
int g0[NU + 1];
int g1[NU + 1];



int main(void) {

    int i, try
        , loop, *infoArray, *codedArray, *estimatedArray, error, bitError, frameError;
    double Eb_N0, BER, FER, *llrArray;
    double complex *transmittedArray;
    char file_w[256], file_w2[256];
    FILE *fp, *fp2;


    fprintf(stderr, "\n-----Simulation Parameters-----\n");
    fprintf(stderr, "\nCODE_LENGTH==%d\nINFO_LENGTH==%d\nNU==%d\nITERATE==%d\n\n", CODE_LENGTH,
            INFO_LENGTH, NU, ITERATE);
    fprintf(stderr, "TRIAL== %d × %d\n\n", TRIAL, LOOP_MAX);

    srandom((unsigned)time(NULL));

    sprintf(file_w, "BER_turbo_K%d_N%d_SR_IT%d_NU%d.txt", INFO_LENGTH, CODE_LENGTH, ITERATE, NU);
    sprintf(file_w2, "FER_turbo_K%d_N%d_SR_IT%d_NU%d.txt", INFO_LENGTH, CODE_LENGTH, ITERATE, NU);

    fp  = fopen(file_w, "w");
    fp2 = fopen(file_w2, "w");


    infoArray        = (int *)calloc(INFO_LENGTH, sizeof(int));
    codedArray       = (int *)calloc(CODE_LENGTH, sizeof(int));
    estimatedArray   = (int *)calloc(INFO_LENGTH, sizeof(int));
    transmittedArray = (double complex *)calloc(CODE_LENGTH, sizeof(double complex));
    llrArray         = (double *)calloc(CODE_LENGTH, sizeof(double));


    // generatorPolynomial
    generatorPolynomial();

    //トレリスのpathの設定
    trelisInit();


    //インターリーバの設定
    interleaveIndex = (int **)calloc(TOTAL_LENGTH, sizeof(int *));
    for (i = 0; i < TOTAL_LENGTH; i++)
        interleaveIndex[i] = (int *)calloc(2, sizeof(int));
    deinterleaveIndex = (int *)calloc(TOTAL_LENGTH, sizeof(int));


    // S-random interleaver
    srandomInterleaver();
    for (i = INFO_LENGTH; i < TOTAL_LENGTH; i++) {
        interleaveIndex[i][0] = i;
        interleaveIndex[i][1] = i;
    }
    for (i = 0; i < TOTAL_LENGTH; i++)
        deinterleaveIndex[interleaveIndex[i][1]] = i;


    // Eb_N0
    for (Eb_N0 = SNR_MIN; Eb_N0 <= SNR_MAX; Eb_N0 += SNR_WIDTH) {

        var        = 1.0 / (2.0 * CODE_RATE * pow(10.0, Eb_N0 / 10.0));
        sigma      = sqrt(var);
        bitError   = 0;
        frameError = 0;
        loop       = 0;

        fprintf(stderr, "Eb/N0==%.2lf\n", Eb_N0);

        while (1) {

            for (try = 0; try < TRIAL; try ++) {

                // Generate information sequence
                for (i = 0; i < INFO_LENGTH; i++)
                    infoArray[i] = random() & 1;

                // Turbo coding
                turboEncoder(infoArray, codedArray);

                // BPSK
                bpsk(codedArray, transmittedArray, CODE_LENGTH);

                // Add noise
                noise(transmittedArray, CODE_LENGTH);

                // LLR
                llrCalculator(transmittedArray, llrArray, CODE_LENGTH);

                // iterative decoder
                turboDecoder(llrArray, estimatedArray);

                // error counter
                error = errorCount(infoArray, estimatedArray, INFO_LENGTH);
                bitError += error;
                if (error > 0)
                    frameError++;
                fprintf(stderr, "try==%d  bitError==%d  frameError==%d\r", (loop * TRIAL + try),
                        bitError, frameError);
            }

            loop++;

            BER = (double)bitError / ((double)INFO_LENGTH * TRIAL * loop);
            FER = (double)frameError / ((double)TRIAL * loop);
            fprintf(stderr, "try==%d  bitError==%d  frameError==%d BER==%.2e FER==%.2e\r",
                    loop * TRIAL, bitError, frameError, BER, FER);

            if (bitError > ERROR_END_CONDITION)
                break;
            if (loop == LOOP_MAX)
                break;
        }


        // BER測定
        BER = (double)bitError / ((double)INFO_LENGTH * TRIAL * loop);
        FER = (double)frameError / ((double)TRIAL * loop);

        fprintf(fp, "%.2lf %e\n", Eb_N0, BER);
        fflush(fp);
        fprintf(fp2, "%.2lf %e\n", Eb_N0, FER);
        fflush(fp2);

        fprintf(stderr, "\nBER==%.3e FER==%.3e\n\n", BER, FER);

        if (BER < BER_END_CONDITION)
            break;
    }

    fclose(fp);
    fclose(fp2);

    return 0;
}


// generatorPolynomial
void generatorPolynomial() {

    int i, temp0, temp1;

    if (NU == 2) {
        temp1 = 05;
        temp0 = 07;

    } else if (NU == 3) {
        temp1 = 015;
        temp0 = 013;

    } else if (NU == 4) {
        temp1 = 021;
        temp0 = 037;

    } else {
        fprintf(stderr, "generator polynominal is NOT found.\n");
        exit(1);
    }
    for (i = 0; i < NU + 1; i++) {
        g1[NU - i] = (temp1 >> i) & 1;
        g0[NU - i] = (temp0 >> i) & 1;
    }
    fprintf(stderr, "g1\n");
    for (i = 0; i < NU + 1; i++)
        fprintf(stderr, "%2d", g1[i]);
    fprintf(stderr, "\ng0\n");
    for (i = 0; i < NU + 1; i++)
        fprintf(stderr, "%2d", g0[i]);
    fprintf(stderr, "\n\n");
}

void srandomInterleaver() {

    int i, counter, temp, idx;
    bool success = true;

    for (i = 0; i < INFO_LENGTH; i++) {
        interleaveIndex[i][0] = i;
        interleaveIndex[i][1] = i;
    }
    for (i = INFO_LENGTH - 1; i > 0; i--) {
        idx = random() % (i + 1);

        temp                    = interleaveIndex[i][1];
        interleaveIndex[i][1]   = interleaveIndex[idx][1];
        interleaveIndex[idx][1] = temp;
    }

    while (1) {

        success = true;
        for (counter = 0; counter < INFO_LENGTH - 1; counter++) {
            if (sorting1(counter) == false) {
                success = false;
                break;
            }
        }

        if (success == true)
            break;

        success = true;
        for (counter = INFO_LENGTH - 1; counter > 0; counter--) {
            if (sorting2(counter) == false) {
                success = false;
                break;
            }
        }

        if (success == true)
            break;
    }
}

bool sorting1(int start) {

    int i, j, temp;

    for (i = start + 1; i < INFO_LENGTH; i++) {
        for (j = 0; j < SPREAD; j++) {
            if (abs(interleaveIndex[start - j][1] - interleaveIndex[i][1]) <= SPREAD)
                break;
            if (start - j > 0 && j != SPREAD - 1)
                continue;

            temp                          = interleaveIndex[start + 1][1];
            interleaveIndex[start + 1][1] = interleaveIndex[i][1];
            interleaveIndex[i][1]         = temp;
            return true;
        }
    }
    return false;
}

bool sorting2(int start) {

    int i, j, temp;

    for (i = start - 1; i >= 0; i--) {
        for (j = 0; j < SPREAD; j++) {
            if (abs(interleaveIndex[start + j][1] - interleaveIndex[i][1]) <= SPREAD)
                break;
            if (start + j < INFO_LENGTH - 1 && j != SPREAD - 1)
                continue;

            temp                          = interleaveIndex[start - 1][1];
            interleaveIndex[start - 1][1] = interleaveIndex[i][1];
            interleaveIndex[i][1]         = temp;
            return true;
        }
    }
    return false;
}


//トレリスの設��
void trelisInit() {

    int bit, i, state, next_state, xor;


    forward_path = (int **)calloc(NUM_OF_STATE, sizeof(int *));
    for (i = 0; i < NUM_OF_STATE; i++)
        forward_path[i] = (int *)calloc(2, sizeof(int));
    backward_path = (int **)calloc(NUM_OF_STATE, sizeof(int *));
    for (i = 0; i < NUM_OF_STATE; i++)
        backward_path[i] = (int *)calloc(2, sizeof(int));
    output_conv = (int **)calloc(NUM_OF_STATE, sizeof(int *));
    for (i = 0; i < NUM_OF_STATE; i++)
        output_conv[i] = (int *)calloc(2, sizeof(int));

    for (state = 0; state < NUM_OF_STATE; state++) {

        for (bit = 0; bit < 2; bit++) {
            xor = bit;
            for (i = 0; i < NU; i++)
                xor ^= (state >> i) & g0[NU - i];
            next_state = (state >> 1) + (xor << (NU - 1));

            forward_path[next_state][bit] = state;
            backward_path[state][bit]     = next_state;

            output_conv[state][bit] = xor&g1[0];
            for (i = 0; i < NU; i++)
                output_conv[state][bit] ^= (state >> i) & g1[NU - i];
        }
    }
}

//ターボ符号化
void turboEncoder(int *infoArray, int *codedArray) {

    int i, output_RSC[2][TOTAL_LENGTH] = {}, termination[2][NU] = {};


    recursiveSystematicConvolution(infoArray, output_RSC[0], false, termination[0]);
    recursiveSystematicConvolution(infoArray, output_RSC[1], true, termination[1]);


    for (i = 0; i < INFO_LENGTH; i++) {
        codedArray[3 * i]     = infoArray[i];
        codedArray[3 * i + 1] = output_RSC[0][i];
        codedArray[3 * i + 2] = output_RSC[1][i];
    }
    for (i = 0; i < NU; i++) {
        codedArray[3 * INFO_LENGTH + 2 * i]     = termination[0][i];
        codedArray[3 * INFO_LENGTH + 2 * i + 1] = output_RSC[0][INFO_LENGTH + i];
    }
    for (i = 0; i < NU; i++) {
        codedArray[3 * INFO_LENGTH + 2 * (NU + i)]     = termination[1][i];
        codedArray[3 * INFO_LENGTH + 2 * (NU + i) + 1] = output_RSC[1][INFO_LENGTH + i];
    }
}


//再帰畳み込み
void recursiveSystematicConvolution(int *infoArray, int *output, bool interleave,
                                    int *termination) {

    int i, j, memory[NU + 1] = {};

    for (i = 0; i < TOTAL_LENGTH; i++) {
        if (i < INFO_LENGTH)
            memory[0] = infoArray[interleaveIndex[i][interleave]];
        else {
            memory[0] = 0;
            for (j = 1; j < NU + 1; j++)
                memory[0] ^= memory[j] & g0[j];
            termination[i - INFO_LENGTH] = memory[0];
        }

        for (j = 1; j < NU + 1; j++)
            memory[0] ^= memory[j] & g0[j];
        output[i] = memory[0] & g1[0];
        for (j = 1; j < NU + 1; j++)
            output[i] ^= memory[j] & g1[j];

        for (j = NU; j > 0; j--)
            memory[j] = memory[j - 1];
    }
}

// LLRの計算
void llrCalculator(double complex *transmittedArray, double *llr, int num) {

    int i;

    for (i = 0; i < num; i++)
        llr[i] = 2.0 * creal(transmittedArray[i]) / var;
}


//反復復号
void turboDecoder(double *llr, int estimatedArray[]) {

    int iterate, t;
    double Lambda1[TOTAL_LENGTH] = {}, Lambda1_e[TOTAL_LENGTH] = {}, Lambda2[TOTAL_LENGTH] = {},
           Lambda2_e[TOTAL_LENGTH] = {};
    double output1[TOTAL_LENGTH][2], output2[TOTAL_LENGTH][2], output3[TOTAL_LENGTH][2],
        termination1[NU][2], termination2[NU][2], y1_llr[INFO_LENGTH];


    for (t = 0; t < INFO_LENGTH; t++) {
        output1[t][1] = -log(1.0 + exp(llr[3 * t]));
        output1[t][0] = llr[3 * t] + output1[t][1];
        y1_llr[t]     = llr[3 * t];

        output2[t][1] = -log(1.0 + exp(llr[3 * t + 1]));
        output2[t][0] = llr[3 * t + 1] + output2[t][1];
        output3[t][1] = -log(1.0 + exp(llr[3 * t + 2]));
        output3[t][0] = llr[3 * t + 2] + output3[t][1];
    }
    for (t = 0; t < NU; t++) {
        termination1[t][1] = -log(1.0 + exp(llr[3 * INFO_LENGTH + 2 * t]));
        termination1[t][0] = llr[3 * INFO_LENGTH + 2 * t] + termination1[t][1];
        termination2[t][1] = -log(1.0 + exp(llr[3 * INFO_LENGTH + 2 * (t + NU)]));
        termination2[t][0] = llr[3 * INFO_LENGTH + 2 * (t + NU)] + termination2[t][1];

        output2[INFO_LENGTH + t][1] = -log(1.0 + exp(llr[3 * INFO_LENGTH + 2 * t + 1]));
        output2[INFO_LENGTH + t][0] =
            llr[3 * INFO_LENGTH + 2 * t + 1] + output2[INFO_LENGTH + t][1];
        output3[INFO_LENGTH + t][1] = -log(1.0 + exp(llr[3 * INFO_LENGTH + 2 * (t + NU) + 1]));
        output3[INFO_LENGTH + t][0] =
            llr[3 * INFO_LENGTH + 2 * (t + NU) + 1] + output3[INFO_LENGTH + t][1];
    }


    for (iterate = 1;; iterate++) {

        // MAP decoder 1
        for (t = 0; t < NU; t++) {
            output1[INFO_LENGTH + t][0] = termination1[t][0];
            output1[INFO_LENGTH + t][1] = termination1[t][1];
        }

        logMAP_decoder(output1, output2, Lambda1, Lambda2_e, false);


        for (t = 0; t < INFO_LENGTH; t++)
            Lambda1_e[t] = Lambda1[t] - y1_llr[t] - Lambda2_e[deinterleaveIndex[t]];
        //ターミネーションビットのLLRの初期化
        for (t = INFO_LENGTH; t < TOTAL_LENGTH; t++)
            Lambda1_e[t] = 0.0;


        // MAP decoder 2
        //ターミネーションビットの入れ替え
        for (t = 0; t < NU; t++) {
            output1[INFO_LENGTH + t][0] = termination2[t][0];
            output1[INFO_LENGTH + t][1] = termination2[t][1];
        }

        logMAP_decoder(output1, output3, Lambda2, Lambda1_e, true);

        if (iterate == ITERATE)
            break;


        for (t = 0; t < INFO_LENGTH; t++)
            Lambda2_e[t] =
                Lambda2[t] - y1_llr[interleaveIndex[t][1]] - Lambda1_e[interleaveIndex[t][1]];
        //ターミネーションビットのLLRの初期化
        for (t = INFO_LENGTH; t < TOTAL_LENGTH; t++)
            Lambda2_e[t] = 0.0;
    }


    for (t = 0; t < INFO_LENGTH; t++) {
        if (Lambda2[deinterleaveIndex[t]] >= 0.0)
            estimatedArray[t] = 0;
        else
            estimatedArray[t] = 1;
    }
}

void logMAP_decoder(double output1[][2], double output2[][2], double *Lambda, double *Lambda_e,
                    bool interleave) {

    int t, l, l_dash, bit;
    double alpha[TOTAL_LENGTH + 1][NUM_OF_STATE] = {}, beta[TOTAL_LENGTH + 1][NUM_OF_STATE],
                                gamma[TOTAL_LENGTH + 1][NUM_OF_STATE][2],
                                log_priori[TOTAL_LENGTH][2], zero, one, likelihood1, likelihood2;


    // calculate log priori probability
    if (interleave == true) {
        for (t = 0; t < TOTAL_LENGTH; t++) {
            log_priori[t][1] = -log(1.0 + exp(Lambda_e[interleaveIndex[t][1]]));
            log_priori[t][0] = Lambda_e[interleaveIndex[t][1]] + log_priori[t][1];
        }
    } else if (interleave == false) {
        for (t = 0; t < TOTAL_LENGTH; t++) {
            log_priori[t][1] = -log(1.0 + exp(Lambda_e[deinterleaveIndex[t]]));
            log_priori[t][0] = Lambda_e[deinterleaveIndex[t]] + log_priori[t][1];
        }
    }


    // gamma initialization
    for (t = 1; t <= TOTAL_LENGTH; t++) {
        for (l = 0; l < NUM_OF_STATE; l++) {
            for (bit = 0; bit < 2; bit++) {

                likelihood1 = output1[interleaveIndex[t - 1][interleave]][bit];
                likelihood2 = output2[t - 1][output_conv[forward_path[l][bit]][bit]];

                gamma[t][l][bit] = log_priori[t - 1][bit] + likelihood1 + likelihood2;
            }
        }
    }


    // forward recursion (alpha calculation)
    alpha[0][0] = 0.0;
    for (l = 1; l < NUM_OF_STATE; l++)
        alpha[0][l] = -INF;

    for (t = 1; t < TOTAL_LENGTH + 1; t++) {
        for (l = 0; l < NUM_OF_STATE; l++) {
            alpha[t][l] = Max_fc(alpha[t - 1][forward_path[l][0]] + gamma[t][l][0],
                                 alpha[t - 1][forward_path[l][1]] + gamma[t][l][1]);
        }
    }

    // backward recursion (beta calculation)
    beta[TOTAL_LENGTH][0] = 0.0;
    for (l = 1; l < NUM_OF_STATE; l++)
        beta[TOTAL_LENGTH][l] = -INF;

    for (t = TOTAL_LENGTH - 1; t > 0; t--) {
        for (l = 0; l < NUM_OF_STATE; l++) {
            // backward_path[l][bit] is equivalent to l_dash
            beta[t][l] =
                Max_fc(beta[t + 1][backward_path[l][0]] + gamma[t + 1][backward_path[l][0]][0],
                       beta[t + 1][backward_path[l][1]] + gamma[t + 1][backward_path[l][1]][1]);
        }
    }

    // LLR calculation
    for (t = 1; t < TOTAL_LENGTH + 1; t++) {

        zero = -INF;
        one  = -INF;

        for (l = 0; l < NUM_OF_STATE; l++) {
            l_dash = forward_path[l][0];
            zero   = Max_fc(zero, alpha[t - 1][l_dash] + gamma[t][l][0] + beta[t][l]);
            l_dash = forward_path[l][1];
            one    = Max_fc(one, alpha[t - 1][l_dash] + gamma[t][l][1] + beta[t][l]);
        }

        Lambda[t - 1] = zero - one;
    }
}

// max with correcting function(Jacobian function)
double Max_fc(double delta1, double delta2) {

    return (fmax(delta1, delta2) + log(1.0 + exp(-fabs(delta2 - delta1))));
}


// BPSK modulation
void bpsk(int *infoArray, double complex *transmittedArray, int num) {

    int i;

    for (i = 0; i < num; i++)
        transmittedArray[i] = 1.0 - 2.0 * infoArray[i];
}


// Add noise
void noise(double complex *transmittedArray, int num) {

    int i;
    double u1, u2, real, imag;


    for (i = 0; i < num; i++) {
        u1 = (random() + 1.0) / (RAND_MAX + 1.0);
        u2 = (random() + 1.0) / (RAND_MAX + 1.0);

        real = creal(transmittedArray[i]) + sqrt(-2.0 * log(u1)) * sigma * cos(2.0 * M_PI * u2);
        imag = cimag(transmittedArray[i]) + sqrt(-2.0 * log(u1)) * sigma * sin(2.0 * M_PI * u2);

        transmittedArray[i] = real + I * imag;
    }
}


//エラーカウンタ
int errorCount(int *infoArray, int *estimatedArray, int num) {

    int i, error = 0;

    for (i = 0; i < num; i++) {
        if (infoArray[i] != estimatedArray[i])
            error++;
    }

    return error;
}
