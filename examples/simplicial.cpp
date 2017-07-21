#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include "../src/matrix.hpp"
#include "../src/homology.hpp"
#include "../src/Smith_normal_form.hpp"
#include "../src/simplicial_homology.hpp"

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
