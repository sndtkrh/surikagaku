#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include "matrix.hpp"
#include "homology.hpp"
#include "Smith_normal_form.hpp"

// 単体複体のホモロジーを計算する

typedef int Point;
typedef std::vector<Point> Simplex;
typedef std::vector<std::vector<Simplex> > Complex; // 単体を次元ごとに分けて持っておく

// 足りない単体を追加して単体複体を作る
Complex get_complex(Complex c){
  int dim = c.size() - 1;
  std::vector< std::set<Simplex> > S(dim+1);
  // 今ある単体を全部入れておく
  for(int d = 0; d <= dim; d++){
    for(int i = 0; i < c[d].size(); i++){
      std::sort( c[d][i].begin(), c[d][i].end() );
      S[d].insert( c[d][i] );
    }
  }
  for(int d = dim-1; d >= 0; d--){
    // add all primary face of each d+1-simplex
    for( Simplex s : S[d+1] ){
      Simplex t(d+1);
      for(int i = 1; i < s.size(); i++){
	t[i-1] = s[i];
      }
      S[d].insert( t );
      for(int i = 0; i < t.size(); i++){
	t[i] = s[i];
	S[d].insert( t );
      }
    }
  }

  Complex ret(dim+1);
  for(int d = 0; d <= dim; d++)
    for(Simplex s : S[d])
      ret[d].push_back( s );
  return ret;
}

std::ostream& operator<< (std::ostream & os, const Complex & c){
  int dim = c.size() - 1;
  for(int d = 0; d <= dim; d++){
    os << "dim=" << d << ", #=" << c[d].size() << "\n";
    for(int i = 0; i < c[d].size(); i++){
      os << "|";
      for(int j = 0; j < c[d][i].size(); j++){
	os << c[d][i][j];
	if( j == d )
	  os << "|, ";
	else
	  os << ",";
      }
    }
    os << "\n";
  }
  return os;
}

class SimplicialHomology : public Homology{
  Complex cpx;
public:
  std::vector< std::map< Simplex, int> > index; // index of d-simplex
  SimplicialHomology( Complex c ) : cpx(c) {
    dim = cpx.size() - 1;
    index.resize( dim + 1 );
    for(int d = 0; d <= dim; d++){
      C.push_back( cpx[d].size() );
      for(int i = 0; i < cpx[d].size(); i++){
	index[d][ cpx[d][i] ] = i;
      }
    }
    boundary_map.resize( dim+2 );
    boundary_map[0] = {}; // boundary_map[0] : C_0 -> 0 is a zero map
    boundary_map[dim+1] = {}; // boundary_map[dim+1] : 0 -> C_{dim} is a zero map
    boundary_map_snf.resize( dim+2 );
    boundary_map_snf[0] = {};
    boundary_map_snf[dim+1] = {};
    for(int d = 1; d <= dim; d++){
      // compute boundary_map_snf[d] : C_d -> C_{d-1}
      Matrix M( C[d-1], std::vector<long long>( C[d], 0 ) ); // 境界準同型の表現行列
      for(Simplex s : cpx[d]){
	int j = index[d][s];
	Simplex t;
	for(int k = 1; k < s.size(); k++){
	  t.push_back( s[k] );
	}
	int i = index[d-1][t];
	M[i][j] = 1;
	int sgn = -1;
	for(int k = 0; k < t.size(); k++){
	  t[k] = s[k];
	  i = index[d-1][t];
	  M[i][j] = sgn;
	  sgn *= -1;
	}
      }
      boundary_map[d] = M;
      boundary_map_snf[d] = compute_snf(M);
    }
    compute_homology();
  }
};


int main(){
  Complex T2 = // 2-dim Torus
    {
      {}, {},
      {
	{0,1,4}, {0,1,6}, {0,2,3}, {0,2,8}, {0,3,4}, {0,6,8},
	{1,2,5}, {1,2,7}, {1,4,5}, {1,6,7},
	{2,3,5}, {2,7,8},
	{3,4,7}, {3,5,6}, {3,6,7},
	{4,7,8}, {4,5,8},
	{5,6,8}
      }
    };
  Complex K2 = // 2-dim Klein bottle
    {
      {}, {},
      {
	{0,1,4}, {0,2,6}, {0,2,3}, {0,1,8}, {0,3,4}, {0,6,8},
	{1,2,5}, {1,2,7}, {1,4,5}, {2,6,7},
	{2,3,5}, {1,7,8},
	{3,4,7}, {3,5,6}, {3,6,7},
	{4,7,8}, {4,5,8},
	{5,6,8}
      }
    };
  Complex P2 = // 2-dim Projective space
    {
      {}, {},
      {
	{0,1,4}, {0,1,8}, {0,3,4}, {0,3,8},
	{1,2,5}, {1,2,7}, {1,4,5}, {1,7,8},
	{2,5,6}, {2,6,7},
	{3,4,7}, {3,5,6}, {3,6,7}, {3,5,8},
	{4,5,8}, {4,7,8}
      }
    };
  Complex S2 = { {}, {}, { {0,1,2}, {0,1,3}, {0,2,3}, {1,2,3} } }; // 2-dim Sphere
  
  T2 = get_complex( T2 );
  std::cout << T2 << std::endl;
  SimplicialHomology H( T2 );
  std::cout << H << std::endl;
}
