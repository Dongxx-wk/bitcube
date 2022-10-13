//
//  bitcube.hpp
//  bitcube
//
//  Created by Kosuke Shimizu on 2019/04/03.
//  Copyright © 2019 Kosuke Shimzu. All rights reserved.
//

#ifndef bitcube_h
#define bitcube_h

#ifndef stdio_h
#include <stdio.h>
#endif
#ifndef stdlib_h
#include <stdlib.h>
#endif
#ifndef math_h
#include <math.h>
#endif

#define PRIMAL_BITCUBE_SIDE 2
#define QDCTBLOCK8x8_SIDE 8
#define QDCTBLOCK8x8_LINEARSIDE 64
#define COEF_BITLENGTH 10
#define NUM_OF_ONES_TO_SCRAMBLE 1
#define NUM_OF_BIT_IN_BITCUBE 8
#define ENCRYPTION_PERCEMTAGE 100
#define BITCUBE_SIDE 2

#define signbit(a) ((a < 0)?(1):(0))
#define sign(a)    ((a)?(-1):(1))
#define abs(a)     ((a < 0)?(-a):(a))
#define abssub(a,b)     (abs(abs(a)-abs(b)))

typedef unsigned char uchar;
typedef unsigned short ushort;

class QDCT {
private:
    short **qdct[3], *qdct_tmp_linear[3];
    unsigned qdct_height[3], qdct_width[3], len[3], mb_ctr, mb_scale[3][2];
    int v_samp_factor[3], h_samp_factor[3];
    FILE *fp[3];
    uchar key[256];
public:
    QDCT(const char *theY, const char *theU, const char *theV){
        //file reading.
        fp[0] = fopen(theY, "rb"); fp[1] = fopen(theU, "rb"); fp[2] = fopen(theV, "rb");
        if(fp[0] == NULL || fp[1] == NULL || fp[2] == NULL){
            printf("NOT exist %s, %s and/or %s\n",theY, theU, theV);
            exit(2);
        }
        for(uchar t = 0; t < 3; t++){//Reading QDCT coefficient files
            char dummy;
            // 从fp读取数据并赋值
            fscanf(fp[t], "%d%c%d%c%d%c%d%c",\
                   &qdct_height[t], &dummy, &qdct_width[t], &dummy,\
                   &v_samp_factor[t], &dummy, &h_samp_factor[t], &dummy);
            qdct_tmp_linear[t] = new short[qdct_height[t] * qdct_width[t]];
            fread(qdct_tmp_linear[t], sizeof(short), qdct_height[t] * qdct_width[t], fp[t]);
            qdct[t] = new short*[qdct_height[t]];
            for (unsigned h = 0; h < qdct_height[t]; h++)
                qdct[t][h] = new short[qdct_width[t]];
            //Following procedure rearranges the linear QDCT data to 2D planes
            //Height and width of each Macroblock (MB)
            // 每一个小块垂直的长度
            mb_scale[t][0] = 8 * v_samp_factor[t];
            // 每一个小块横向的长度
            mb_scale[t][1] = 8 * h_samp_factor[t];
            //Vertical/horizontal numbers of the MBs
            unsigned num_mb[2] =\
            {qdct_height[t] / mb_scale[t][0], qdct_width[t] / mb_scale[t][1]};
            //One-time reading amount of the linera data
            unsigned mb_factor = mb_scale[t][0] * mb_scale[t][1];
            //Total amount of the linear amount
            len[t] = num_mb[0] * num_mb[1] * mb_factor;
            mb_ctr = 0;
            for (unsigned l = 0; l < len[t]; l += mb_factor){
                unsigned l_tmp = l;
                //Computing initial-coordinate of each MB
                unsigned first_mb[2] =\
                {(mb_ctr / num_mb[1]) * mb_scale[t][0], (mb_ctr % num_mb[1]) * mb_scale[t][1]};
                mb_ctr++;
                for (uchar j_mb = 0; j_mb < v_samp_factor[t]; j_mb++){
                    //Computing initial-coordinate within each MB
                    unsigned first_in_mb[2] = {8 * j_mb + first_mb[0], 0};
                    for (uchar i_mb = 0; i_mb < h_samp_factor[t]; i_mb++){
                        first_in_mb[1] = 8 * i_mb + first_mb[1];
                        for (unsigned j_blk = first_in_mb[0]; j_blk < first_in_mb[0] + 8; j_blk++)
                            for (unsigned i_blk = first_in_mb[1]; i_blk < first_in_mb[1] + 8; i_blk++)
                                qdct[t][j_blk][i_blk] = qdct_tmp_linear[t][l_tmp++];
                    }
                }
            }
            //Finarization of the 1D-to-2D transform
        }
        fclose(fp[0]); fclose(fp[1]); fclose(fp[2]);
    }
    ~QDCT(){
        for (uchar t = 0; t < 3; t++) delete[] qdct_tmp_linear[t];
        for (uchar t = 0; t < 3; t++){
            for (unsigned h = 0; h < qdct_height[t]; h++)
                delete[] qdct[t][h];
            delete[] qdct[t];
        }
    }
    void write(const char *theY, const char *theU, const char *theV){

        for(uchar t = 0; t < 3; t++){
            //2D-to-1D transform
            //Height and width of each Macroblock (MB)
            mb_scale[t][0] = 8 * v_samp_factor[t];
            mb_scale[t][1] = 8 * h_samp_factor[t];
            //Vertical/horizontal numbers of the MBs
            unsigned num_mb[2] =\
            {qdct_height[t] / mb_scale[t][0], qdct_width[t] / mb_scale[t][1]};
            //One-time reading amount of the linera data
            unsigned mb_factor = mb_scale[t][0] * mb_scale[t][1];
            //Total amount of the linear amount
            len[t] = num_mb[0] * num_mb[1] * mb_factor;
            mb_ctr = 0;
            for (unsigned l = 0; l < len[t]; l += mb_factor){
                unsigned l_tmp = l;
                //Computing initial-coordinate of each MB
                unsigned first_mb[2] =\
                {(mb_ctr / num_mb[1]) * mb_scale[t][0], (mb_ctr % num_mb[1]) * mb_scale[t][1]};
                mb_ctr++;
                for (uchar j_mb = 0; j_mb < v_samp_factor[t]; j_mb++){
                    //Computing initial-coordinate within each MB
                    unsigned first_in_mb[2] = {8 * j_mb + first_mb[0], 0};
                    for (uchar i_mb = 0; i_mb < h_samp_factor[t]; i_mb++){
                        first_in_mb[1] = 8 * i_mb + first_mb[1];
                        for (unsigned j_blk = first_in_mb[0]; j_blk < first_in_mb[0] + 8; j_blk++)
                            for (unsigned i_blk = first_in_mb[1]; i_blk < first_in_mb[1] + 8; i_blk++)
                                qdct_tmp_linear[t][l_tmp++] = qdct[t][j_blk][i_blk];
                    }
                }
            }
            //Finarization of the 2D-to-1D transform
        }
        fp[0] = fopen(theY, "wb"); fp[1] = fopen(theU, "wb"); fp[2] = fopen(theV, "wb");
        if(fp[0] == NULL || fp[1] == NULL || fp[2] == NULL){
            printf("NOT created %s, %s and/or %s\n",theY, theU, theV);
            exit(3);
        }
        for(uchar t = 0; t < 3; t++){
            fprintf(fp[t], "%d,%d,%d,%d,", qdct_height[t], qdct_width[t], v_samp_factor[t], h_samp_factor[t]);
            fwrite(qdct_tmp_linear[t], sizeof(short), len[t], fp[t]);
        }
        fclose(fp[0]); fclose(fp[1]); fclose(fp[2]);
    }
    void input_key(unsigned k[8]){
        for(uchar u = 0; u < 8; u++)
            for(uchar t = 0; t < 32; t++)
                key[u * 32 + t] = (k[u] >> (32 - t - 1)) & 0x01;
    }
    void apply_seed(){
        unsigned theseed = 0;
        static unsigned state = 0;

        for(uchar t = 0; t < 32; t++)
            theseed += key[(state + t) % 256] << (31 - t);
        state++;
        
        srand(theseed);
    }
    void bitcube_scrambling(){
        uchar bitcube_8x8[QDCTBLOCK8x8_SIDE][QDCTBLOCK8x8_SIDE][COEF_BITLENGTH];
        short block_8x8[QDCTBLOCK8x8_SIDE][QDCTBLOCK8x8_SIDE];
        unsigned tmp_height, tmp_width;

        for(uchar t = 0; t < 3; t++){
            tmp_height = qdct_height[t] - (qdct_height[t] % QDCTBLOCK8x8_SIDE);
            tmp_width = qdct_width[t] - (qdct_width[t] % QDCTBLOCK8x8_SIDE);
            for(unsigned j = 0; j < tmp_height; j += QDCTBLOCK8x8_SIDE)
                for(unsigned i = 0; i < tmp_width; i += QDCTBLOCK8x8_SIDE){
                    for(uchar j_blk = 0; j_blk < QDCTBLOCK8x8_SIDE; j_blk++)
                        for(uchar i_blk = 0; i_blk < QDCTBLOCK8x8_SIDE; i_blk++)
                            block_8x8[j_blk][i_blk] = qdct[t][j + j_blk][i + i_blk];
                    bitcuboid_decompose_8x8(bitcube_8x8, block_8x8);
                    apply_seed();
                    bitcube_scramble_8x8(bitcube_8x8);
                    bitcuboid_idecompose_8x8(bitcube_8x8, block_8x8);
                    for(uchar j_blk = 0; j_blk < QDCTBLOCK8x8_SIDE; j_blk++)
                        for(uchar i_blk = 0; i_blk < QDCTBLOCK8x8_SIDE; i_blk++)
                            qdct[t][j + j_blk][i + i_blk] = block_8x8[j_blk][i_blk];
                }
        }

    }
    void bitcuboid_decompose_8x8(uchar bitcube[][QDCTBLOCK8x8_SIDE][COEF_BITLENGTH],\
                                 short qdct_coefblock[][QDCTBLOCK8x8_SIDE]){

        for(uchar j = 0; j < QDCTBLOCK8x8_SIDE; j++)
            for(uchar i = 0; i < QDCTBLOCK8x8_SIDE; i++)
                fixedlength_signed_binarization(bitcube[j][i], qdct_coefblock[j][i], COEF_BITLENGTH);
    }
    void bitcuboid_idecompose_8x8(uchar bitcube[][QDCTBLOCK8x8_SIDE][COEF_BITLENGTH],\
                                  short qdct_coefblock[][QDCTBLOCK8x8_SIDE]){

        for(uchar j = 0; j < QDCTBLOCK8x8_SIDE; j++)
            for(uchar i = 0; i < QDCTBLOCK8x8_SIDE; i++)
                fixedlength_signed_redecimalization(bitcube[j][i], qdct_coefblock[j][i], COEF_BITLENGTH);
    }
    void fixedlength_signed_binarization(uchar binstring[],\
                             short thedecimal,\
                             uchar bitlength = COEF_BITLENGTH){
        uchar the_sign;
        ushort tmp_coef;
        the_sign = signbit(thedecimal);
        tmp_coef = (ushort)abs((int)thedecimal);
        binstring[0] = the_sign;
        for(uchar b = 0; b < bitlength - 1; b++)
            binstring[bitlength - b - 1] = (tmp_coef >> b) & 0x01;
    }
    void fixedlength_signed_redecimalization(uchar binstring[],\
                               short &thedecimal,\
                               uchar bitlength = COEF_BITLENGTH){
        uchar the_sign;
        ushort tmp_coef;
        the_sign = binstring[0];
        tmp_coef = 0;
        for(uchar b = 0; b < bitlength - 1; b++)
            tmp_coef += binstring[bitlength - b - 1] * (1 << b);
        if(the_sign == 1 && tmp_coef == 0)
            thedecimal = -(1 << (bitlength - 1));
        else
            thedecimal = sign(the_sign) * tmp_coef;
    }
    void bitcube_scramble_8x8(uchar bitcube[][QDCTBLOCK8x8_SIDE][COEF_BITLENGTH]){
        uchar the_bitcube[NUM_OF_BIT_IN_BITCUBE];
        uchar *indices = new uchar[NUM_OF_BIT_IN_BITCUBE];
        uchar detected_number_of_ones = NUM_OF_ONES_TO_SCRAMBLE;
        uchar tmp_num_of_ones;

        for(uchar k = 1; k < COEF_BITLENGTH - 1; k += 2)
            for(uchar j = 0; j < QDCTBLOCK8x8_SIDE; j += 2)
                for(uchar i = 0; i < QDCTBLOCK8x8_SIDE; i += 2){
                    tmp_num_of_ones  = the_bitcube[0] = bitcube[j    ][i    ][k    ];
                    tmp_num_of_ones += the_bitcube[1] = bitcube[j    ][i + 1][k    ];
                    tmp_num_of_ones += the_bitcube[2] = bitcube[j + 1][i    ][k    ];
                    tmp_num_of_ones += the_bitcube[3] = bitcube[j + 1][i + 1][k    ];
                    tmp_num_of_ones += the_bitcube[4] = bitcube[j    ][i    ][k + 1];
                    tmp_num_of_ones += the_bitcube[5] = bitcube[j    ][i + 1][k + 1];
                    tmp_num_of_ones += the_bitcube[6] = bitcube[j + 1][i    ][k + 1];
                    tmp_num_of_ones += the_bitcube[7] = bitcube[j + 1][i + 1][k + 1];
                    if((tmp_num_of_ones > 0 && tmp_num_of_ones <= detected_number_of_ones) &&\
                       (rand() % 100 + 1) <= ENCRYPTION_PERCEMTAGE){
                        //If number of the counted ones is matched
                        for(uchar t = 0; t < NUM_OF_BIT_IN_BITCUBE; t += 2){
                            indices[t] = t;
                            indices[t + 1] = t + 1;
                        }
                        uchar tmp_idx, tmp;
                        for(uchar t = 1; t < NUM_OF_BIT_IN_BITCUBE; t++){
                            tmp_idx = rand() % t;
                            tmp = indices[tmp_idx];
                            indices[tmp_idx] = indices[0];
                            indices[0] = tmp;
                        }
                        //permuting the bits within the bitcube
                        for(uchar t = 0; t < NUM_OF_BIT_IN_BITCUBE; t += 2)
                        {
                            tmp = the_bitcube[indices[t]];
                            the_bitcube[indices[t]] = the_bitcube[indices[t + 1]];
                            the_bitcube[indices[t + 1]] = tmp;
                        }
                        bitcube[j    ][i    ][k    ] = the_bitcube[0];
                        bitcube[j    ][i + 1][k    ] = the_bitcube[1];
                        bitcube[j + 1][i    ][k    ] = the_bitcube[2];
                        bitcube[j + 1][i + 1][k    ] = the_bitcube[3];
                        bitcube[j    ][i    ][k + 1] = the_bitcube[4];
                        bitcube[j    ][i + 1][k + 1] = the_bitcube[5];
                        bitcube[j + 1][i    ][k + 1] = the_bitcube[6];
                        bitcube[j + 1][i + 1][k + 1] = the_bitcube[7];
                    }
                }
        delete [] indices;
    }
};

#endif /* bitcube_h */
