#include <iostream>
#include <vector>
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <math.h>
#include <complex>
#include <cstdlib>

using namespace itpp;
using namespace std;

float factor  = 1.0 ; //beta fator(temporary)
int No_of_bits = 10000 ;
int bst = 7, usr = 7, cel = 7 ; //no. of base stations, users and cells (assuming 1bst per cell)
int usrno = 0 ; //simulation for interference to be calculate for this user number_temporary 
vec EbN0dB, EbN0, N0, noise_variance, bit_error_rate; //vec is a vector containing double
bvec transmitted_bits, received_bits, rxbits; //bvec is a vector containing bits
cvec transmitted_symbols(5000), buff(5000), cnoise(5000), cbuff, noise, noise1, noise2(10000);  //cvec is a vector containing double_complex
std::complex<double> received_symbols[10][5][10][5000], rxnoise_symbols[10][5][10][5000], rx[7][7][5000], w[7][10000];
std::complex<double> interference[5000];
float beta[7][1][7];
AWGN_Channel awgn_channel;     //The AWGN channel class

vector<vector<int> > si;
int si_inv =  1;

class Mobile
{
public:
	Mobile(float x, float y) { this->x = x; this->y = y;}
private:
	float x, y;
};

class Cell
{
public:
	Cell(float x, float y, float radius) : x(x), y(y), radius(radius) {};
	void addUser(Mobile &mobile) { mobiles.push_back(mobile);}
private:
	std::vector<Mobile> mobiles;
	float x, y, radius;
};

class BTS
{
public:
	BTS(float x, float y) : x(x), y(y) {} ;	
private:
	float x, y;
};

class lsfade
{
public:
	lsfade(int i, int k, int l) : i(i), k(k), l(l) {} ;   
private:
	int i, k, l ;
};

void Betav(void)
{
	float base_loc[7][2] = {{0,0},{1,1.73},{0,2},{1, -1.73},{-1,-1.73},{0,-2},{-1,1.73}};
	float user_loc[7][2] = {{0,0.5}, {1.25,2.165}, {0,2.5}, {1.25,-2.165}, {-1.25, -2.165}, {0,-2.5}, {-1.25,2.165}};	
	int D = 0; int cl =0;

	for(int i = 0; i < bst; i++)
	{
		for(int j = 0 ; j <1; j++)
		{	
			for(int cl = 0 ; cl <7; cl++)
			{	
				D = sqrt((pow((base_loc[i][0]-user_loc[cl][0]),2)) + (pow((base_loc[i][1]-user_loc[cl][1]),2)));
				beta[i][j][cl] = pow((4*3.14*D*3),2) ; //friis					 
			}
		}
	}
        //out<<beta[0][0][0];	
	/* float n=0 ;
   	n = 4*3.14*D*3 ; //friis
	n = pow(n, 2);	
	vector < vector < vector<int> > > Beta; */
}

void Lsfadev(void)
{
        void interfer(void) ;
	void zf(void) ;
	void uplink(complex<double> Lsfade[7][1][7][10000]);
	//void uplink(std::complex<double> Lsfade);
	//Declarations of classes:
	
	BPSK bpsk;
	QPSK qpsk;                     //The QPSK modulator class
	it_file ff;                    //For saving the results to file
	BERC berc;                     //Used to count the bit errors
	Real_Timer tt;                 //The timer used to measure the execution time
	//Reset and start the timer:
	
	tt.tic();	
	//Init
	double Ec = 1.0;                      //The transmitted energy per QPSK symbol is 1.
	double Eb = Ec/2.0;                   //The transmitted energy per bit is 0.5.
	vec EbN0dB = linspace(0.0, 30.0, 18); //Simulate for 10 Eb/N0 values from 0 to 30 dB.
	vec EbN0 = inv_dB(EbN0dB);            //Calculate Eb/N0 in a linear scale instead of dB.
	vec N0 = Eb*pow(EbN0, -1.0);          //N0 is the variance of the (complex valued) noise.
	bit_error_rate.set_size(EbN0dB.length(), false);
	//Randomize the random number generators in it++:
	RNG_randomize();
	//RNG_randomize();

	for (int m = 0; m < EbN0dB.length(); m++)
	{		
		//Generate a vector of random bits to transmit:
		transmitted_bits = randb(No_of_bits);
		//Modulate the bits to QPSK symbols:
		transmitted_symbols = qpsk.modulate_bits(transmitted_bits);
		//Set the noise variance of the AWGN channel:
		awgn_channel.set_noise(N0(m));		
	
		std::complex<double> Lsfade[7][1][7][10000];
		for(int p = 0; p < bst; p++)
		{
			for(int q = 0; q < usr; q++)
			{
				for(int r = 0; r < cel; r++)
				{   
					// vec c1buff = randn(No_of_bits) ;
					cbuff = randn_c(No_of_bits);
					for(int b=0; b < No_of_bits; b++)
					{	//std::complex<double> mycomplex (c1buff[b],c1buff[b]);
						//mycomplex = mycomplex*1/sqrt(2);
						Lsfade[p][q][r][b] = cbuff[b];					 
					}
				}
			}
		}
	
		uplink(Lsfade) ;

//------------------------------------------------------------------------------downlink--------------------------------------------------------------------------- //

		//for (int p = 0; p < bst; p++) //from every base station
	//	{
			int p = 0, power1 =1  ;
			for ( int r = 0; r < cel ; r++) //for every cell
			{				
				for (int q = 0; q < usr; q++) //for every user
				{				
					buff.clear() ;

					for (int b = 0; b < transmitted_symbols.size(); b++)
					{
						received_symbols[r][q][p][b] += sqrt(power1*beta[r][0][0])*Lsfade[r][0][0][b]*w[r][b]*transmitted_symbols[b];
						//cout<<"Lsfade[" <<p<< "][" << q << "][" << r << "][" << b << "] = " << rxnoise_symbols[6][q][r][b] << endl; 
					}						
				//cout << "Lsfade[" << p << "][" << q << "][" << r << "][" << b << "] = " << Lsfade[p][q][r][b] << endl; 					
				}
			}							

				for(int r = 0 ; r < 7 ; r++)
				{
					for(int q = 0 ; q < 1 ; q++)
					{
						for (int b = 0; b < transmitted_symbols.size(); b++)
						{
							rx[0][0][b] += received_symbols[r][q][p][b] ;
						}
					}
				}	
		
				for (int b = 0; b < transmitted_symbols.size(); b++)
					buff[b] = rx[0][0][b];

				//Run the transmited symbols through the channel using the () operator for adding white noise:
				noise = awgn_channel(buff);	
				
				for (int b = 0; b < transmitted_symbols.size(); b++)
				{
					//rxnoise_symbols[p][q][r][b] = noise[b]/Lsfade[p][q][r][b]; //zero forcing
					rx[0][0][b] = noise[b]; 
					cnoise[b] = noise[b];						
					//cout << "Lsfade[" << p << "][" << q << "][" << r << "][" << b << "] = " << rxnoise_symbols[p][q][r][b] << endl;
				}
				//cout<<received_bits<<endl ;

				noise.clear() ;
		
//------------------------------------------------------------------------------------------------------------------------------------------------------------------//

		//interfer() ;
		//zf() ;
		//mmse() ;
		//Demodulate the received QPSK symbols into received bits:
		received_bits = qpsk.demodulate_bits(cnoise);
		//cout<< cnoise <<endl;
		//cout<<transmitted_symbols[4888]<<endl ;
		//cout<<Lsfade[6][0][6][4888]<<endl ;
		//cout<<received_bits.size()<<endl;
		//Calculate the bit error rate:

		berc.clear();                               //Clear the bit error rate counter
		berc.count(transmitted_bits, received_bits); //Count the bit errors
		bit_error_rate(m) = berc.get_errorrate();   //Save the estimated BER in the res
	}

	tt.toc();
	//Print the results:
	cout << endl;
	//cout<<interference[466] - cnoise[466] <<endl;
        //  cout << received_symbols[6][0][0][916];
	cout << "EbN0dB = " << EbN0dB << " [dB]" << endl;
	cout << "BER = " << bit_error_rate << endl;
	cout << "Saving results to ./result.it" << endl;	
	//Save the results to file
	ff.open("result.it");     
	ff << Name("EbN0dB") << EbN0dB;
	ff << Name("ber") << bit_error_rate;
	ff.close();
	//Exit_program
}

void interfer(void)
{    
	for (int b = 0; b < transmitted_symbols.size(); b++)
	{
		for (int p = 0; p < bst; p++) //for every base station
		{
			for ( int r = 0; r < cel ; r++) //for every cell
			{
				interference[b]+= rxnoise_symbols[p][usrno][r][b];				
			}		
		}
	}

	for (int b = 0; b < transmitted_symbols.size(); b++)
	{
		cnoise[b] = interference[b];
		interference[b] = 0;                  
	}

}


void zf(void)
{
	for (int b = 0; b < transmitted_symbols.size(); b++)
	{
		cnoise[b] =  rxnoise_symbols[0][0][0][b];
	}
}

/*/
void mmse(void)
{
float h_mmse[bst][cel];
int tau = 0 , pwr = 0;
for(int j = 0 ; j < cel ; j++)
{
	for(int l = 0 ; l < bst ; l++)
	{
		h_mmse[j][l] = sqrt(beta[j][0][l])*sqrt(pwr*tau)/(1+ beta[j][0][l]);
}
}
return;
}
*/

void uplink(complex<double> Lsfade[7][1][7][10000])
{
	
	int power = 1, K = 1;
	float norm1 = 0, alpha[7];
	std::complex<double> pil_out[7][10000];
	std::complex<double> norm_f = 1, pilot_bst[10000];
	std::complex<double> g_hat[7][10000];
	for(int i = 0 ; i< usr; i++)
	
	{
		vector <int> v;
		si.push_back(v);
		for(int j = 0 ; j< K; j++)
		{
			if(i==j)
				si[i].push_back(1);
			else
				si[i].push_back(0);
		}
	}

	for(int x = 0 ; x < 7 ; x++)
	{
		for(int d = 0 ; d < 10000 ; d++)
		{
			for(int z = 0 ; z < 7 ; z++)
			{
				for(int y = 0 ; y < 1 ; y++)
				{
					pil_out[x][d] += sqrt((beta[x][y][z])*power)*Lsfade[x][y][z][d]*si[x][y]; //to be edited			
					}			
				}			
			}
	
	

	for(int i=0;i<10000;i++)
	{
		noise2[i] = pil_out[x][i];
	}
	noise1 = awgn_channel(noise2);	

	for(int i = 0 ; i < 10000; i++)
	{
		g_hat[x][i] = noise1[i]/K ;
		norm1 += pow(abs(g_hat[x][i]),2);
        }
	norm1 = sqrt(norm1) ;	

	alpha[x] = norm1/sqrt(10000); //scalar alpha - normalization factor


//beamforming: w = g/||g||

   	for(int i = 0 ; i < 10000; i++)
		w[x][i] = g_hat[x][i]/norm1;

	}
	

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
	bst = 7 ; 
	cel = 7 ; 
	usr = 1 ;   //usr is users per cell    
	Betav();
	Lsfadev() ;
        //uplink() ;
	return 0;
}
