#include <math.h>
#include <stdlib.h>
#include <stdio.h>


void four2 (char **data){
    int row, column, n, k = 0;
    float spectrum[256][512]
    float temp[512];
    FILE *output1;
    FILE *input;
    n=256;
    
    // same as dynamic 2-d array
    for (row = 0; row < n; row ++) {
        for (column = 0; column < n; column++) {
            spectrum[row][k+=2] = data[row][column];
        }
    }
    for (row = 0; row < n; row ++) {
        for (column = 0; column < 2*n; column++) {
            temp[column] = spectrum[row][column];
        }
        four1(spectrum[row]-1, 2*n, 1);
        //four1(temp-1, 2*n, 1);
        //for (k=0; k<512; k++) {
        //    spectrum[row][k] = temp[k];
        //}
    }
    
    for (column = 0; column < 512; column++) {
        for (row=0; row<256; row++) {
            temp[row] = spectrum[row][column];
        }
        
    }
    
    return 0;
}

