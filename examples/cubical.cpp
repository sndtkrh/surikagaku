#include <vector>
#include <algorithm>
#include <iostream>
#include "../src/cubical_homology.hpp"

int main(){
  // square
  Cube<2> e1( {Interval(0,true), Interval(0,false)} );
  Cube<2> e2( {Interval(1,true), Interval(0,false)} );
  Cube<2> e3( {Interval(0,false), Interval(0,true)} );
  Cube<2> e4( {Interval(0,false), Interval(1,true)} );
  std::vector< std::vector< Cube<2> > > square = { {}, {e1,e2,e3,e4} };

  // regular hexahedron
  Cube<3> f1( {Interval(0,false), Interval(0,false), Interval(0,true)} );
  Cube<3> f2( {Interval(0,false), Interval(0,true), Interval(0,false)} );
  Cube<3> f3( {Interval(0,true), Interval(0,false), Interval(0,false)} );
  Cube<3> f4( {Interval(0,false), Interval(0,false), Interval(1,true)} );
  Cube<3> f5( {Interval(0,false), Interval(1,true), Interval(0,false)} );
  Cube<3> f6( {Interval(1,true), Interval(0,false), Interval(0,false)} );
  std::vector< std::vector< Cube<3> > > hexahedron = { {}, {}, {f1,f2,f3,f4,f5,f6} };

  square = get_cubicalset( square );
  std::cout << square << std::endl;
  CubicalHomology<2> cubical_homology(square);
  std::cout << cubical_homology << std::endl;
}

