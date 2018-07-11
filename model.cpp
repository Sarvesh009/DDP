#include <iostream>
#include <vector>
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
#include <math.h>
#include <complex>
#include <cstdlib>


using namespace itpp;
using namespace std;

extern int No_of_bits = 10000 ;
extern double Ec, Eb;
extern vec EbN0dB, EbN0, N0, noise_variance, bit_error_rate; //vec is a vector containing double
extern bvec transmitted_bits, received_bits;                 //bvec is a vector containing bits
extern 	cvec transmitted_symbols, received_symbols, buff;           //cvec is a vector containing double_complex

 

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

class beta
{

public:
	beta(int i, int k, int l) : i(i), k(k), l(l) {} ;   
void Betav(void)
{

//Declarations of classes:
  BPSK bpsk;
  QPSK qpsk;                     //The QPSK modulator class
  AWGN_Channel awgn_channel;     //The AWGN channel class
  it_file ff;                    //For saving the results to file
  BERC berc;                     //Used to count the bit errors
  Real_Timer tt;                 //The timer used to measure the execution time
  

 double Ec = 1.0;                      //The transmitted energy per QPSK symbol is 1.
  double Eb = Ec / 2.0;                 //The transmitted energy per bit is 0.5.
  vec EbN0dB = linspace(0.0, 30.0, 10); //Simulate for 10 Eb/N0 values from 0 to 30 dB.
  vec EbN0 = inv_dB(EbN0dB);         //Calculate Eb/N0 in a linear scale instead of dB.
  vec N0 = Eb * pow(EbN0, -1.0);     //N0 is the variance of the (complex valued) noise.
  vec bit_error_rate;
  bit_error_rate.set_size(EbN0dB.length(), false);

 RNG_randomize();

for (int m = 0; m < EbN0dB.length(); m++) {
    //Show how the simulation progresses:
    cout << "Now simulating Eb/N0 value number " << m + 1 << " of " << EbN0dB.length() << endl;
    //Generate	 a vector of random bits to transmit:
    bvec transmitted_bits = randb(No_of_bits);

    //Modulate the bits to QPSK symbols:
    cvec transmitted_symbols = qpsk.modulate_bits(transmitted_bits);
    //Set the noise variance of the AWGN channel:
               // awgn_channel.set_noise(N0(i));
    //Run the transmited symbols through the channel using the () operator:
               cvec received_symbols = awgn_channel(transmitted_symbols);
    //Demodulate the received QPSK symbols into received bits:
    bvec received_bits = qpsk.demodulate_bits(received_symbols);
    //Calculate the bit error rate:
    berc.clear();                               //Clear the bit error rate counter
    berc.count(transmitted_bits, received_bits); //Count the bit errors
    bit_error_rate(i) = berc.get_errorrate();   //Save the estimated BER in the result vector

}

vector < vector < vector< complex<double> > > > Beta;
  for(int p = 0; p <i; p++)
  {
    vector < vector < complex<double> > > w;
    Beta.push_back( w );
    for(int q = 0; q < k; q++)
    {
      vector< complex<double> > v;
      Beta[p].push_back( v );
      for(int r = 0; r < l; r++)
     {   

   vec cbuff = randn(No_of_bits) ;

           std::complex<double> mycomplex (1,2);
           mycomplex= mycomplex*1/sqrt(2);
        Beta[p][q].push_back(mycomplex);
      }
    }
  }

  for (size_t p = 0; p < Beta.size(); p++)
    for (size_t q = 0; q < Beta[p].size(); q++)
      for (size_t r = 0; r < Beta[p][q].size(); r++)
        cout << "Beta[" << p<< "][" << q << "][" << r << "] = " << Beta[p][q][r] << endl;  
cout<<No_of_bits<<endl ;
}
private:
	int i, k, l ;
	
};

class lsfade
{
public:
	lsfade(int i, int k, int l) : i(i), k(k), l(l) {} ;   

void lsfadev()
{

  vector < vector < vector<int> > > Lsfade;
  for(int a = 0; a < i; a++)
  {
    vector < vector < int > > w;
    Lsfade.push_back( w );
    for(int b = 0; b < k; b++)
    {
      vector <int> v;
      Lsfade[i].push_back( v );
      for(int c = 0; c < l; c++)
      {
	//k = 1/sqrt(2)*
        Lsfade[a][b].push_back(0);
      }
    }
  }
}

private:
	int i, k, l ;
	
};

int main(void)
{	
	Cell c(0.0, 0.0, 1);
	Mobile m(1.0, 0.3);
	c.addUser(m);
        beta mybeta(7,1,7);
	lsfade(7,1,7);
       	mybeta.Betav() ;

	return 0;
}
