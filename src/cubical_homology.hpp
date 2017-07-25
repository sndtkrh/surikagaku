#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <set>
#include "Smith_normal_form.hpp"
#include "homology.hpp"

// interval [l,l+1] or [l,l]
struct Interval {
  int coord;
  bool degenerate;
  Interval(int c, bool d) : coord(c), degenerate(d) {
    return;
  }
  bool is_in( const Interval & i ){
    if( i.degenerate ){
      if( degenerate ){
	return coord == i.coord;
      }else{
	return (coord == i.coord || coord + 1 == i.coord);
      }
    }
  }
  friend std::ostream& operator<<(std::ostream & os, const Interval & i){
    if( i.degenerate ){
      os << "[" << i.coord << "]";
    }else{
      os << "[" << i.coord << "," << i.coord + 1 << "]";
    }
    return os;
  }
};
bool operator==(const Interval & i1, const Interval & i2){
  return (i1.coord == i2.coord) && (i1.degenerate == i2.degenerate);
}
bool operator<(const Interval & i1, const Interval & i2){
  return (i1.coord < i2.coord) || ( (i1.coord == i2.coord) && i1.degenerate && ! i2.degenerate);
}
bool operator>(const Interval & i1, const Interval & i2){
  return !(i1 == i2) && !(i1 < i2);
}

// cube
template<int DIM>
class Cube{
public:
  std::vector<Interval> intervals;
  Cube(std::vector<Interval> is) : intervals(is) {
    if( intervals.size() != DIM ){
      throw "dimension size error";
    }
    return;
  }
  int dim(){
    int d = 0;
    for(Interval i : intervals){
      if( ! i.degenerate ){
	d++;
      }
    }
    return d;
  }
  int emb(){
    return intervals.size();
  }
  friend std::ostream& operator<<(std::ostream & os, const Cube<DIM> & c){
    int d = 1;
    for(Interval i : c.intervals){
      os << i;
      if( d < DIM ){
	os << "x";
      }
      d++;
    }
    return os;
  }
};
template <int DIM>
bool operator==(const Cube<DIM> & c1, const Cube<DIM> & c2){
  return c1.intervals == c2.intervals;
}
template <int DIM>
bool operator<(const Cube<DIM> & c1, const Cube<DIM> & c2){
  for(int i = 0; i < DIM; i++){
    if( c1.intervals[i] < c2.intervals[i] ){
      return true;
    }else if( c1.intervals[i] > c2.intervals[i] ){
      return false;
    }
  }
  return false;
}


template <int DIM>
class CubicalHomology : public Homology {
  typedef std::vector< std::vector<Cube<DIM> > > CubicalSet;
  CubicalSet cubicalset;
  std::vector< std::map<Cube<DIM>, int> > index;
public:
  CubicalHomology(CubicalSet cs) : cubicalset(cs) {
    dim = cubicalset.size() - 1;
    index.resize( dim + 1 );
    for(int d = 0; d <= dim; d++){
      C.push_back( cubicalset[d].size() );
      for(int i = 0; i < cubicalset[d].size(); i++){
	index[d][ cubicalset[d][i] ] = i;
      }
    }
    boundary_map.resize( dim + 2 );
    boundary_map[0] = {};
    boundary_map[dim+1] = {};
    boundary_map_snf.resize( dim + 2 );
    boundary_map_snf[0] = {};
    boundary_map_snf[dim+1] = {};
    for(int d = 1; d <= dim; d++){
      Matrix M( C[d-1], std::vector<long long>( C[d], 0 ) );
      for( Cube<DIM> & c : cubicalset[d] ){
	int idx = index[d][c];
	long long sign = 1;
	for( Interval & i : c.intervals ){
	  if( i.degenerate ) continue;
	  int boundary;
	  i.degenerate = true;
	  i.coord++;
	  boundary = index[d-1][c];
	  M[boundary][idx] = sign;
	  i.coord--;
	  boundary = index[d-1][c];
	  M[boundary][idx] = -sign;
	  i.degenerate = false;
	  sign *= -1;
	}
      }
      boundary_map[d] = M;
      boundary_map_snf[d] = compute_snf(M);
    }
    compute_homology();
  }
};

// 足りない基本方体を追加する
template <int DIM>
std::vector<std::vector<Cube<DIM> > > get_cubicalset( std::vector<std::vector<Cube<DIM> > > cs ){
  int dim = cs.size() - 1;
  std::vector<std::set<Cube<DIM> > > S(dim+1);
  // 今ある基本方体を全部入れておく
  for(int d = 0; d <= dim; d++){
    for(int i = 0; i < cs[d].size(); i++){
      S[d].insert( cs[d][i] );
    }
  }
  for(int d = dim-1; d >= 0; d--){
    // add all primary face of each d+1-cube
    for( Cube<DIM> c : S[d+1] ){
      for(int i = 0; i < DIM; i++){
	if( ! c.intervals[i].degenerate ){
	  c.intervals[i].degenerate = true;
	  S[d].insert( c );
	  c.intervals[i].coord++;
	  S[d].insert( c );
	  c.intervals[i].coord--;
	  c.intervals[i].degenerate = false;
	}
      }
    }
  }
  std::vector<std::vector<Cube<DIM>>> ret(dim+1);
  for(int d = 0; d <= dim; d++)
    for(Cube<DIM> c : S[d])
      ret[d].push_back( c );
  return ret;
}

template <int DIM>
std::ostream& operator<<(std::ostream & os, std::vector<std::vector<Cube<DIM> > > & cubicalset){
  int dim = 0;
  for( auto & cs : cubicalset ){
    os << "dim=" << dim << ", #=" << cs.size() << "\n";
    dim++;
    for( Cube<DIM> & c : cs ){
      os << c << ", ";
    }
    os << "\n";
  }
  return os;
}
