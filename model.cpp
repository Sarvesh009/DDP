#include <iostream>
#include <vector>
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <math.h>
#include <complex>
#include <cstdlib>

using namespace itpp;
using namespace std;

float factor  = 0 ; //beta fator(temporary)
int No_of_bits = 10000 ;
int bst = 1, usr = 1, cel = 1 ; //no. of base stations, users and cells (assuming 1bst per cell)
int usrno = 0 ; //simulation for interference to be calculate for this user number_temporary 
vec EbN0dB, EbN0, N0, noise_variance, bit_error_rate; //vec is a vector containing double
bvec transmitted_bits, received_bits, rxbits; //bvec is a vector containing bits
cvec transmitted_symbols(5000), buff(5000), cnoise(5000), cbuff, noise;  //cvec is a vector containing double_complex
std::complex<double> received_symbols[10][5][10][5000], rxnoise_symbols[10][5][10][5000];
std::complex<double> interference[5000];
float beta[7][1][1] ;
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
	float base_loc[][2] = {{0,0},{0,2},{1.73,1},{1.73,-1},{0,-2},{-1.73,-1},{-1.73,1}};
	float user_loc[][2] =  {{0.5,0.5}};	
	for(int i = 0; i < bst; i++)
	{
		for(int j = 0 ; j <usr; j++)
			{
				for(int k= 0 ; k<cel; k++)
					{
					beta[i][j][k]= sqrt((pow((base_loc[i][0]-user_loc[0][0]),2)) + (pow((base_loc[i][1]-user_loc[0][1]),2)));
					beta[i][j][k]= pow((4*3.14*beta[i][j][k]*3),2) ; //friis					 
					}
			}
	}
        //8cout<<beta[0][0][0] ;	
	int D = 2;
	float n=0 ;
        n= 4*3.14*D*3 ; //friis
	n = pow(n, 2);
	
	
	vector < vector < vector<int> > > Beta;
	
}

void Lsfadev(void)
{
        void interfer(void) ;
	void zf(void) ;
	//Declarations of classes:
	BPSK bpsk;
	QPSK qpsk;                     //The QPSK modulator class
	AWGN_Channel awgn_channel;     //The AWGN channel class
	it_file ff;                    //For saving the results to file
	BERC berc;                     //Used to count the bit errors
	Real_Timer tt;                 //The timer used to measure the execution time

	//Reset and start the timer:
	tt.tic();	

	//Init:
	double Ec = 1.0;                      //The transmitted energy per QPSK symbol is 1.
	double Eb = Ec / 2.0;                 //The transmitted energy per bit is 0.5.
	vec EbN0dB = linspace(0.0, 30.0, 18); //Simulate for 10 Eb/N0 values from 0 to 30 dB.
	vec EbN0 = inv_dB(EbN0dB);         //Calculate Eb/N0 in a linear scale instead of dB.
	vec N0 = Eb * pow(EbN0, -1.0);     //N0 is the variance of the (complex valued) noise.

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
		std::complex<double> Lsfade[bst][usr][cel][No_of_bits];
		for(int p = 0; p < bst; p++)
		{
			for(int q = 0; q < usr; q++)
			{
				for(int r = 0; r < cel; r++)
				{   
					// vec c1buff = randn(No_of_bits) ;
					 cbuff = randn_c(No_of_bits);
					for(int b=0; b < No_of_bits; b++)
					{
						//std::complex<double> mycomplex (c1buff[b],c1buff[b]);
						//mycomplex= mycomplex*1/sqrt(2);
						Lsfade[p][q][r][b] = cbuff[b];
					 
					}
				}
			}
		}

		//for h*x cvec
		for (int p = 0; p < bst; p++) //for every base station
		{
			for (int q = 0; q <  usr; q++) //for every user
			{
				for ( int r = 0; r < cel ; r++) //for every cell
				{
					buff.clear() ;
					for (int b = 0; b < transmitted_symbols.size(); b++)
					{
						received_symbols[p][q][r][b] = Lsfade[p][q][r][b]*transmitted_symbols[b] ;
						if(p!=0)
						{
							received_symbols[p][q][r][b] = factor*received_symbols[p][q][r][b] ;
							//cout << "Lsfade[" << p << "][" << q << "][" << r << "][" << b << "] = " << rxnoise_symbols[6][q][r][b] << endl; 
								
						}
						
						//cout << "Lsfade[" << p << "][" << q << "][" << r << "][" << b << "] = " << Lsfade[p][q][r][b] << endl; 
						buff[b] = received_symbols[p][q][r][b];
					}
					
					//Run the transmited symbols through the channel using the () operator for adding white noise:
					noise = awgn_channel(buff);

					for (int b = 0; b < transmitted_symbols.size(); b++)
					{
						rxnoise_symbols[p][q][r][b] = noise[b]/Lsfade[p][q][r][b]; //zero forcing
						//              rxnoise_symbols[p][q][r][b] = noise[b]; 
						
						//cout << "Lsfade[" << p << "][" << q << "][" << r << "][" << b << "] = " << rxnoise_symbols[p][q][r][b] << endl;
					}
					//cout<<received_bits<<endl ;
					noise.clear() ;
				}
			}
		}

		interfer() ;

		//zf() ;

		mmse() ;

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
	cout << endl;
	//Save the results to file:
	ff.open("result.it");
	ff << Name("EbN0dB") << EbN0dB;
	ff << Name("ber") << bit_error_rate;
	ff.close();
	//Exit program:


}



void interfer(void)
{
        
	for (int b = 0; b < transmitted_symbols.size(); b++)
		{
			for (int p = 0; p < bst; p++) //for every base station
			{
				for ( int r = 0; r < cel ; r++) //for every cell
				{
				  interference[b]+= rxnoise_symbols[p][usrno][r][b] ;
					
				}		
			}
		}



	for (int b = 0; b < transmitted_symbols.size(); b++)
	{
		cnoise[b] =  interference[b];
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


void mmse(void)
{

float h_mmse[bst][cel];

for(int j = 0 ; j < cel ; j++)
{
	for(int l = 0 ; l < bst ; l++)
	{
		h_mmse[j][l] = sqrt(beta[j][0][l])*sqrt(pwr*tau)/(1+ beta[j][0][l]);
	}
}


}



int main(void)
{	
	void Lsfadev(void) ;
	void Betav(void) ;

	Cell c(0.0, 0.0, 1);
	Mobile m(1.0, 0.3);
	c.addUser(m); 
	bst=7 ; cel = 1 ; usr =1 ;       
	Betav();
	Lsfadev() ;
        
	return 0;
}
