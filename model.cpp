#include <iostream>
#include <vector>
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <math.h>
#include <complex>
#include <cstdlib>


using namespace itpp;
using namespace std;
float factor  = 0.2	 ;
int No_of_bits = 10000 ;
double Ec, Eb;
vec EbN0dB, EbN0, N0, noise_variance, bit_error_rate; //vec is a vector containing double
bvec transmitted_bits, received_bits, rxbits;                 //bvec is a vector containing bits
cvec transmitted_symbols(5000),  cbuff, noise;           //cvec is a vector containing double_complex
cvec buff(5000), cnoise(5000) ;
std::complex<double> received_symbols[10][5][10][5000];
std::complex<double> rxnoise_symbols[10][5][10][5000];


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



void Lsfadev(void)
{
int i =2, k=1 ,l=2;
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
;
 //RNG_randomize();

for (int m = 0; m < EbN0dB.length(); m++) {
    //Show how the simulation progresses:
    //cout << "Now simulating Eb/N0 value number " << m + 1 << " of " << EbN0dB.length() << endl;
    //Generate	 a vector of random bits to transmit:
    transmitted_bits = randb(No_of_bits);


    //Modulate the bits to QPSK symbols:
    transmitted_symbols = qpsk.modulate_bits(transmitted_bits);

    //Set the noise variance of the AWGN channel:
       awgn_channel.set_noise(N0(m));
    //Run the transmited symbols through the channel using the () operator:
               //cvec received_symbols = awgn_channel(transmitted_symbols);
    //Demodulate the received QPSK symbols into received bits:
  //  bvec 
    //Calculate the bit error rate:
   

vector< vector < vector < vector< complex<double> > > > > Lsfade;
  for(int p = 0; p <i; p++)
  {
    vector <vector< vector < complex<double> > > > w;
    Lsfade.push_back( w );
    for(int q = 0; q < k; q++)
    {
      vector< vector < complex<double> > > v;
      Lsfade[p].push_back( v );
      for(int r = 0; r < l; r++)
     {   
	vector< complex<double> >  y;
      	Lsfade[p][q].push_back( y );

        // vec c1buff = randn(No_of_bits) ;
	 cbuff = randn_c(No_of_bits);

	for(int b=0; b < No_of_bits; b++)
	{
	   //std::complex<double> mycomplex (c1buff[b],c1buff[b]);
	   //mycomplex= mycomplex*1/sqrt(2);
           Lsfade[p][q][r].push_back(cbuff[b]);

          
	}
      }
    }
  }



//for h*x cvec
  for (int p = 0; p < Lsfade.size(); p++)
{
    for (int q = 0; q < Lsfade[p].size(); q++)
{
      for ( int r = 0; r < Lsfade[p][q].size(); r++)
{
	buff.clear() ;

	for (int b = 0; b < transmitted_symbols.size(); b++)
{
		received_symbols[p][q][r][b] = factor*Lsfade[p][q][r][b]*transmitted_symbols[b] ;
 		//cout << "Lsfade[" << p << "][" << q << "][" << r << "][" << b << "] = " << Lsfade[p][q][r][b] << endl; 
		buff[b] = received_symbols[p][q][r][b] ;
}

		noise = awgn_channel(buff);
     
	for (int b = 0; b < transmitted_symbols.size(); b++)
{
		rxnoise_symbols[p][q][r][b] = noise[b]/Lsfade[p][q][r][b]	 ;
		//cout << "Lsfade[" << p << "][" << q << "][" << r << "][" << b << "] = " << rxnoise_symbols[p][q][r][b] << endl; 

}

//cout<<received_bits<<endl ;
noise.clear() ;
}
}
}

for (int b = 0; b < transmitted_symbols.size(); b++)
{
cnoise[b] =  rxnoise_symbols[0][0][0][b];
//cout<< noise <<endl ;
//cout<<"rxnoise"<<rxnoise_symbols[1][0][1][4999]<<endl;
}



received_bits = qpsk.demodulate_bits(cnoise);

//cout<< cnoise <<endl;
//cout<<transmitted_symbols[4888]<<endl ;
//cout<<Lsfade[6][0][6][4888]<<endl ;

//cout<<received_bits.size()<<endl;
  berc.clear();                               //Clear the bit error rate counter
    berc.count(transmitted_bits, received_bits); //Count the bit errors
    bit_error_rate(m) = berc.get_errorrate();   //Save the estimated BER in the res

}


tt.toc();
  //Print the results:
  cout << endl;
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

class Beta
{
public:
	Beta(int i, int k, int l) : i(i), k(k), l(l) {} ;   

void Betav()
{

  vector < vector < vector<int> > > Beta;
  for(int a = 0; a < i; a++)
  {
    vector < vector < int > > w;
    Beta.push_back( w );
    for(int b = 0; b < k; b++)
    {
      vector <int> v;
      Beta[i].push_back( v );
      for(int c = 0; c < l; c++)
      {
	//k = 1/sqrt(2)*
        Beta[a][b].push_back(0);
      }
    }
  }
}

private:
	int i, k, l ;
	
};

int main(void)
{	
void Lsfadev(void) ;
	Cell c(0.0, 0.0, 1);
	Mobile m(1.0, 0.3);
	c.addUser(m);
        
	Beta(7,1,7);
       Lsfadev() ;

	return 0;
}
