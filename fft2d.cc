    // Distributed two-dimensional Discrete FFT transform
    // Marcellus Pleasant

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Complex.h"
#include "InputImage.h"
#include <thread>
#include <math.h>
constexpr unsigned int NUMTHREADS = 4;

using namespace std;

//undergrad students can assume NUMTHREADS will evenly divide the number of rows in tested images
//graduate students should assume NUMTHREADS will not always evenly divide the number of rows in tested images.
// I will test with a different image than the one given

//--------------------  Transpose, inspired by online examples  ------------------//
void transpose(Complex* oldA, Complex* newA, int w){
    for (int i = 0; i < w; i++ ) {
        for (int j = 0; j < w; j++ ) {
            int index1 = i*w+j;
            int index2 = j*w+i;
            newA[index2] = oldA[index1];
        }
    }

    for (int i=0; i<w*w; i++) {
        oldA[i] = newA[i];
    }
}

void Transform1D(Complex* h, int w, Complex* H, bool second, int Th)
{
    // Implement a simple 1-d DFT using the double summation equation
    // given in the assignment handout.  h is the time-domain input
    // data, w is the width (N), and H is the output array.
    Complex sum;

    if(second == false)
    for(int row=0;row<Th;row++){
        for(int col=0;col<w;col++){ //n = i
            for(int k=0;k<w;k++){
                /*real = h[(i*w)+k].real;
                imag = h[(i*w)+k].imag;
                hold.real = cos((2*M_PI*i*k)/w) * real;
                hold.imag = -1*sin((2*M_PI*i*k)/w) * imag;*/
                Complex hold(cos((2*M_PI*col*k)/w), -sin((2*M_PI*col*k)/w) );
                //cout << H[k].imag << endl;
                sum = sum + (hold * h[(row*w)+k]);
            }
            H[(row*w)+col] = sum;
            //sum.Print();
            //printf("\n");
            sum.real = 0;
            sum.imag = 0;
        }
    }


    //If we are ready to do the second transform, before I used transpose
    /*if(second == true) {
        int bigass = 0;
        for (int col = 0; col < Th; col++) {
            for (int row = 0; row < w; row++) { //n = i
                for (int k = 0; k < w; k++) {
                    Complex hold(cos((2 * M_PI * (row) * k) / w), -sin((2 * M_PI * (row) * k) / w));
                    sum = sum + (hold * h[(w * k) + col]);
                }
                H[bigass] = sum; // insert
                bigass++;
                //sum.Print();
                //printf("\n");
                sum.real = 0;
                sum.imag = 0;
            }
        }

        cout << bigass << endl;
    }*/
}


//Reverse func here

void R_Transform1D(Complex* h, int w, Complex* H, bool doCol, int Th){

        // Implement a simple 1-d DFT using the double summation equation
        // given in the assignment handout.  h is the time-domain input
        // data, w is the width (N), and H is the output array.
        Complex sum; double NV = (double)1/w; Complex N(NV,0);

        //for rows
        if(doCol == false)
            for(int row=0;row<Th;row++){ //cout << "here" << endl;
                for(int col=0;col<w;col++){ //n = i
                    for(int k=0;k<w;k++){
                        /*real = h[(i*w)+k].real;
                        imag = h[(i*w)+k].imag;
                        hold.real = cos((2*M_PI*i*k)/w) * real;
                        hold.i mag = -1*sin((2*M_PI*i*k)/w) * imag;*/
                        Complex hold(cos((2*M_PI*col*k)/w), sin((2*M_PI*col*k)/w) );
                        //cout << H[k].imag << endl;//
                        sum = sum + (hold * h[(row*w)+k]);
                    }
                    H[(row*w)+col] = ( sum*N );
                    //sum.Print();
                    //printf("\n");
                    sum.real = 0;
                    sum.imag = 0;
                }
            }

        // for columns
        /*if(doCol == true) {
            int bigass = 0;
            for (int col = 0; col < Th; col++) {
                for (int row = 0; row < w; row++) { //n = i
                    for (int k = 0; k < w; k++) {
                        Complex hold(cos((2 * M_PI * (row) * k) / w), sin((2 * M_PI * (row) * k) / w));
                        sum = sum + (hold * h[(w * k) + col]);
                    }
                    H[bigass] = (sum * N); // insert
                    bigass++;
                    //sum.Print();
                    //printf("\n");
                    sum.real = 0;
                    sum.imag = 0;
                }
            }
        }*/
}

void Transform2D(const char* inputFN)
    { // Do the 2D transform here.
        // 1) Use the InputImage object to read in the Tower.txt file and
        //    find the width/height of the input image.
        // 2) Create a vector of complex objects of size width * height to hold
        //    values calculated
        // 3) Do the individual 1D transforms on the rows assigned to each thread
        // 4) Force each thread to wait until all threads have completed their row calculations
        //    prior to starting column calculations
        // 5) Perform column calculations
        // 6) Wait for all column calculations to complete
        // 7) Use SaveImageData() to output the final results

        // Step (1) in the comments is the line above.
        // Your code here, steps 2-7

        InputImage image(inputFN);  // Create the helper object for reading the image
        Complex* temp = image.GetImageData();
        double widthHeight = image.GetWidth();
        int Mp = (int)(widthHeight*widthHeight);
        int Th = (int) widthHeight/NUMTHREADS; //How many rows there are
        int Wh = (int) widthHeight * Th; //The actual starting position in the array
        Complex* savetemp;
        savetemp = new Complex[Mp];

        //------------------------------------  Transform rows  ------------------------------//
        vector<thread*> threadL;
        for(int k=0; k<NUMTHREADS;k++){
            threadL.push_back(new thread{Transform1D,&temp[k*Wh],widthHeight ,&savetemp[k*Wh], false, Th } );

        }

        for(auto &t : threadL) {
            t->join();
        }

        threadL.clear();

        //------------------------------------  Transform columns  -------------------------------//
        transpose(savetemp, temp, widthHeight);
        for(int k=0; k<NUMTHREADS;k++){
            threadL.push_back(new thread{Transform1D,&temp[k*Wh],widthHeight ,&savetemp[k*Wh], false, Th } );
        }

        for(auto &t : threadL) {
            t->join();
        }

        threadL.clear();
        transpose(savetemp, temp, widthHeight);

        //-------------------------------------  Print objects to file  -------------------------//
        string fn2("MyAfter2D.txt");
        image.SaveImageData(fn2.c_str(), temp, widthHeight, widthHeight);


        //-----------------------------------  Reverse columns  ------------------------------------//
        transpose(temp, savetemp, widthHeight);
        for(int k=0; k<NUMTHREADS;k++){
            threadL.push_back(new thread{R_Transform1D,&savetemp[k*Wh],widthHeight ,&temp[k*Wh], false, Th } );
        }

        for(auto &t : threadL) {
            t->join();
        }

        threadL.clear();
        transpose(temp, savetemp, widthHeight);

        //-------------------------------------  Reverse rows  -----------------------------------//
        for(int k=0; k<NUMTHREADS;k++){
            threadL.push_back(new thread{R_Transform1D,&savetemp[k*Wh],widthHeight ,&temp[k*Wh], false, Th } );
        }

        for(auto &t : threadL) {
            t->join();
        }

        threadL.clear();

        //-------------------------------------  Print objects to file  -------------------------//
        string fn3("MyAfterInverse.txt");
        image.SaveImageData(fn3.c_str(), temp, widthHeight, widthHeight);

        delete savetemp;

    }


int main(int argc, char** argv)
{
    string fn("Tower.txt"); // default file name
    if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
    Transform2D(fn.c_str()); // Perform the transform.

}



