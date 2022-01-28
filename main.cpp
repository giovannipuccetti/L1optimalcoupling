//
//  main.cpp
//  Computation of the minimal L1 euclidean w distance between
//  two Gaussian distribution
//  computed via the swapping algorithm
//
//  Created by Giovanni Puccetti on 4/02/20.
//  Copyright © 2016 Giovanni Puccetti. All rights reserved.
//

#include <random>
#include <iostream>
#include <cstdio>
#include <ctime>
int main()
{
    // utility variables
    float buf1,buf2,fin_cost,WW;
    // INPUT VARIABLES
    // number of discretization points
    int d=100000;
    // number of simulations
    int MM=50;
    // accuracy
    float acc=0.00001;
    // computation times
    double T[MM];
    // estimates
    float est[MM];
    // iteration times
    float I[MM];
    //initialize variables
    for( int t = 0; t < MM; t++ ){
        T[t]=0;
        est[t]=0;
        I[t]=0;
    }
    //
    float X[d][2];
    float Y[d][2];
    //iteration by M
    for(int t=0;t<MM;t++){
        std::clock_t start;
        start = std::clock();
        //construction of multivariate marginal samples using the Cholesky decomposition
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> dd(0,1);
       for( int i = 0; i < d; i++ ){
           //X[i][0]=sqrt(2)*dd(gen);
           X[i][0]=1.341640786499874*dd(gen);
           X[i][1]=-0.298142396999972*X[i][0]+1.308094458023239*dd(gen);
           //X[i][1]=dd(gen);
           Y[i][0]=dd(gen);
            //Y[i][1]=(Y[i][0]+sqrt(3)*dd(gen))/2;
           Y[i][1]=0.400000000000000*Y[i][0]+0.916515138991168*dd(gen);
          // Y[i][1]=1.000000000000000*dd(gen);
           //Y[i][1]=sqrt(2)*dd(gen);
        }
        fin_cost=0;
        for( int i = 0; i < d; i++ ){
            fin_cost=fin_cost+(sqrt((X[i][0]-Y[i][0])*(X[i][0]-Y[i][0])+(X[i][1]-Y[i][1])*(X[i][1]-Y[i][1])))/d;
            //end of for
        }
        WW=100000;
        while(fabs(WW-fin_cost)>acc)
        {
            I[t]=I[t]+1;
            WW=fin_cost;
            //swapping procedure
            for ( int i = 0; i < (d-1); i++ )
                for ( int j = i+1; j < d; j++ )
                {
                    if ( sqrt((X[i][0]-Y[i][0])*(X[i][0]-Y[i][0])+(X[i][1]-Y[i][1])*(X[i][1]-Y[i][1]))+
                        sqrt((X[j][0]-Y[j][0])*(X[j][0]-Y[j][0])+(X[j][1]-Y[j][1])*(X[j][1]-Y[j][1]))>
                        sqrt((X[i][0]-Y[j][0])*(X[i][0]-Y[j][0])+(X[i][1]-Y[j][1])*(X[i][1]-Y[j][1]))+
                        sqrt((X[j][0]-Y[i][0])*(X[j][0]-Y[i][0])+(X[j][1]-Y[i][1])*(X[j][1]-Y[i][1])) )
                    {
                        buf1=Y[i][0];
                        buf2=Y[i][1];
                        Y[i][0]=Y[j][0];
                        Y[i][1]=Y[j][1];
                        Y[j][0]=buf1;
                        Y[j][1]=buf2;
                        //end of if
                    }
                    //end of 2-for
                }
            fin_cost=0;
            for( int i = 0; i < d; i++ ){
                fin_cost=fin_cost+(sqrt((X[i][0]-Y[i][0])*(X[i][0]-Y[i][0])+(X[i][1]-Y[i][1])*(X[i][1]-Y[i][1])))/d;
                //end of for
            }
            //end of while
        }
        T[t]=( std::clock() - start ) / (double) CLOCKS_PER_SEC;
        est[t]=fin_cost;
        //end of iteration by M
    }
    //computation of final mean estimates
    
    float estimate=0;
    float time=0;
    float iteration=0;
    for( int t = 0; t < MM; t++ ){
        estimate=estimate+est[t]/MM;
        time=time+T[t]/MM;
        iteration=iteration+I[t]/MM;
    }
    printf("Sample number is %i", d);
    printf("      ");
    printf("Mean estimate is %4.8f", estimate);
    printf("      ");
    printf("Mean number of while iteration is %4.4f", iteration);
    printf("      ");
    printf("Mean computation time is %4.4f", time);
    //float buf3=sqrt(1.3);
    //printf("Number is %4.8f", buf3);
    return 0;
    //end of procedure
}
