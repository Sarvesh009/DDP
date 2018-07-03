#include <iostream>
#include <vector>
#include <itpp/itbase.h>
#include <cstdlib>

using namespace itpp;
using namespace std;

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
        Beta[a][b].push_back(0);
      }
    }
  }
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
        Lsfade[a][b].push_back(0);
      }
    }
  }
}

private:
	int i, k, l ;
	
};

int
main(void)
{	Cell c(0.0, 0.0, 1);
	Mobile m(1.0, 0.3);
	c.addUser(m);
        beta(1,2,6);
	lsfade(1,2,6);
	return 0;
}
