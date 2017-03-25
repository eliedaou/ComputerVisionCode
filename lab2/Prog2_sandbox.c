/* Program 1 - Elie Daou
 * Due (was) Wed 2 February 2017 9:00am
 * This program reads an image and outputs either:
 * - output1 will be the left half of a 256x256 image
 * - output2 will be the entire 256x256 image scaled down in half with an output of 128x128
 * */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

const unsigned char black = 0;

void four1(data,nn,isign)
float data[];
int nn,isign;{
    int n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    float tempr,tempi;
    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) {
        if (j > i) {
            SWAP(data[j],data[i]);
            SWAP(data[j+1],data[i+1]);
        }
        m=n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax=2;
    while (n > mmax) {
        istep=2*mmax;
        theta=6.28318530717959/(isign*mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) {
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}


int main (int argc, char *argv[]){
    int row, column, n;
    int k = 0;
    //float **data, **spectrum, **output;
    float data[256][256];
    char image[256][256];
    float spectrum[256][512];
    float output[256][512];
    float temp[512];
    unsigned char pixel;
    //FILE *output1;
    
    FILE *input;
    
    //input image file
    input = fopen(argv[1], "rb");
    if (input == NULL){
        return 1;
    }
    
    n=256;

    
    // first, you must read in the data from the original picture. same as program 1
    for (row = 0; row < n; row ++) {
        for (column = 0; column < n; column++) {
            fread(&data[row][column], sizeof(char), 1, input);
            image[row][column] = (char)data[row][column];
        }
    }
    fclose(input);

 
    /* since fourl destroys the input array with the output */
    /* use the output of the FFT as the input value */
//    for (row = 0; row < n; row ++) {
//        for (column = 0; column < n; column++) {
//            spectrum[row][k+=2] = data[row][column];
//        }
//    }
//
//    for (row=0; row<n; row++) {
//        for (column = 0; column < 2*n; column++) {
//            output[row][column] = spectrum[row][column];
//        }
//    }
    

    
    for (row=0; row<n; row++) {
        printf("%d ", row);
        for (column = 0; column < n; column++) {
            if ((row % 10 == 0) && (column % 10 == 0) ) {
                printf ("%c ", image[row][column]);
            }
        }
        printf("\n");
        
    }


    
    return 0;
}

