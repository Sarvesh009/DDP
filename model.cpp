#include <iostream>
#include <vector>
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <math.h>
#include <complex>
#include <cstdlib>

using namespace itpp;
using namespace std;

float factor  = 0.6 ; //beta fator(temporary)
int No_of_bits = 10000 ;
int bst = 7, usr = 4, cel = 7 ; //no. of base stations, users and cells (assuming 1bst per cell)
int usrno = 0 ; //simulation for interference to be calculate for this user number_temporary 
vec EbN0dB, EbN0, N0, noise_variance, bit_error_rate; //vec is a vector containing double
bvec transmitted_bits, received_bits, rxbits; //bvec is a vector containing bits
cvec transmitted_symbols(5000), buff(5000), cnoise(5000), cbuff, noise, noise1, noise2(5000);  //cvec is a vector containing double_complex
std::complex<double> received_symbols[10][5][10][5000], rxnoise_symbols[10][5][10][5000];
std::complex<double> interference[5000];
float beta[7][4][1];
AWGN_Channel awgn_channel;     //The AWGN channel class

		for(int z = 0 ; z < cel ; y++)
		{
				for(int d = 0 ; d < 10000 ; d++)
				{
					pil_out[x][d] = sqrt((beta[x][y][z])*power)*Lsfade[x][y][z][d]*si[x][y]; //to be edited
				}
		}
			
	}
	

	for(int i=0;i<5000;i++)
	{
		noise2[i] = pil_out[0][i];
	}

	noise1 = awgn_channel(noise2);	

	for(int i = 0 ; i < 50000; i++)
	{
		g_hat[i] = noise2[i]/K ;
		norm += pow(g_hat[i],2) ;
    }

    norm = sqrt(norm) ;


	//beamforming
	// w = g^/||g||
    for(int i = 0 ; i < 50000; i++)
    	w[i] = g_hat[i]*si_inv/norm  ;
	

//plot uplink

//plot downlink


/*
int buf_d = 0;
alpha[k][l] = pow*beta[l][j[k] + 1/K ;
for(int l=0; l<cel ; l++)
{
sinr_d = P[k][l]*pow(beta[l][k][i],2)/pow(alpha,2);
buf_d+= sinr_d ;
}
i = 0;
sinr_d = P[k][l]*pow(beta[l][k][i],2)/pow(alpha,2);
sinr_d = sinr_d/buf_d ;
*/
	return;

}


int main(void)
{	
	void Lsfadev(void) ;
	void Betav(void) ;
	//void uplink(void)

	Cell c(0.0, 0.0, 1);
	Mobile m(1.0, 0.3);
	c.addUser(m); 
	bst=7 ; cel = 7 ; usr =1 ;   //usr is users per cell    
	Betav();
	Lsfadev() ;
        //uplink() ;
	return 0;
}
