#include <vector>
#include <iostream>
#include "../src/simplicial_homology.hpp"

int main(){
  Complex Point = { { {0} } }; // one point

  Complex S1 = // 1-dim Sphere
    { {}, { {0,1}, {1,2}, {2,0} } };
  S1 = get_complex( S1 );
  Complex S2 = // 2-dim Sphere
    { {}, {}, { {0,1,2}, {0,1,3}, {0,2,3}, {1,2,3} } };
  S2 = get_complex( S2 );
  Complex S3 = // 3-dim Sphere
    { {}, {}, {}, { {0,1,2,3}, {0,1,2,4}, {0,1,3,4}, {0,2,3,4}, {1,2,3,4} } };
  S3 = get_complex( S3 );

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
  T2 = get_complex( T2 );
  Complex T2T2 = // T^2 # T^2
    {
      {}, {},
      {
	{0,1,4}, {0,1,6}, {0,2,3}, {0,2,8}, {0,3,4}, {0,6,8},
	{1,2,5}, {1,2,7}, {1,4,5}, {1,6,7},
	{2,3,5}, {2,7,8},
	{3,4,7}, {3,5,6}, {3,6,7},
	{5,6,8},
	{9,10,4}, {9,10,13}, {9,11,12}, {9,11,8}, {9,12,4}, {9,13,8},
	{10,11,5}, {10,11,7}, {10,4,5}, {10,13,7},
	{11,12,5}, {11,7,8},
	{12,4,7}, {12,5,13}, {12,13,7},
	{5,13,8}
      }
    };
  T2T2 = get_complex( T2T2 );

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
  K2 = get_complex( K2 );

  Complex K2K2 = // K^2 # K^2
    {
      {}, {},
      {
	{9,1,4}, {9,2,6}, {9,2,3}, {9,1,8}, {9,3,4}, {9,6,8},
	{1,2,5}, {1,2,7}, {1,4,5}, {2,6,7},
	{2,3,5}, {1,7,8},
	{3,4,7}, {3,5,6}, {3,6,7},
	{5,6,8},
	{-9,-1,-4}, {-9,-2,-6}, {-9,-2,-3}, {-9,-1,-8}, {-9,-3,-4}, {-9,-6,-8},
	{-1,-2,-5}, {-1,-2,-7}, {-1,-4,-5}, {-2,-6,-7},
	{-2,-3,-5}, {-1,-7,-8},
	{-3,-4,-7}, {-3,-5,-6}, {-3,-6,-7},
	{-5,-6,-8},
	{-4,4,5}, {-4,-5,5}, {-7,4,7}, {-4,-7,4},
	{-7,-8,7}, {-8,8,7}, {-5,-8,8}, {-5,5,8}
      }
    };
  K2K2 = get_complex( K2K2 );

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
  P2 = get_complex( P2 );

  std::cout << T2T2 << std::endl;
  SimplicialHomology H( T2T2 );
  std::cout << H << std::endl;
}
