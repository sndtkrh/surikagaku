#ifndef HOMOLOGY
#define HOMOLOGY
#include <tuple>
#include <vector>
#include <iostream>
#include <iomanip>

#include "matrix.hpp"

class Homology{
protected:
  std::vector<int> C; // chain complex
  int dim;
  std::vector< Matrix > boundary_map;
  std::vector< std::vector<long long> > boundary_map_snf; // Smith normal form
  std::vector< int > Hrank; // rank of H^d
  std::vector< std::vector<int> > Htor; // torsion of H^d

  void compute_homology(){
    Hrank.resize(dim+1);
    Htor.resize(dim+1);
    for(int d = 0; d <= dim; d++){
      Hrank[d] = C[d] - boundary_map_snf[d].size() - boundary_map_snf[d+1].size();
      for(int a : boundary_map_snf[d+1]){
	if( a != 1 ){
	  Htor[d].push_back( a );
	}
      }
    }
  }
public:
  int get_dim(){
    return dim;
  }
  std::vector<int> get_chain_complex(){
    return C;
  }
  std::tuple<std::vector<int>, std::vector<std::vector<int>>> get_homology(){
    return std::make_tuple( Hrank, Htor );
  }

  void print_boundary_map(int d){
    std::cout << "boundary map : C^" << d << " -> C^"<< d-1 << std::endl;
    if( boundary_map[d].size() == 0 ){
      std::cout << "size = ( 0 * 0 )" << std::endl;
    }else{
      std::cout << "size = (" << boundary_map[d].size() << " * " << boundary_map[d][0].size() << ")" << std::endl;
      std::cout << boundary_map[d] << std::endl;
    }
  }
  void print_boundary_map_snf(int d){
    std::cout << "boundary map : C^" << d << " -> C^"<< d-1 << std::endl;
    std::cout << "Smith normal form (size=" <<  boundary_map_snf[d].size() << ")" << std::endl;
    std::cout << boundary_map_snf[d] << std::endl;
  }
  void print_boundary_map(){
    for(int d = 0; d <= dim; d++){
      print_boundary_map(d);
    }
  }
  void print_boundary_map_snf(){
    for(int d = 0; d <= dim; d++){
      print_boundary_map_snf(d);
    }
  }
  void print_homology(int d){
    printf("H_%d = Z^{%d} ", d, Hrank[d]);
    for(int tor : Htor[d] ){
      printf("(+) Z/%dZ ", tor );
    }
    printf("\n");
  }
  void print_homology(){
    for(int d = 0; d <= dim; d++){
      print_homology(d);
    }
  }

  friend std::ostream& operator<< (std::ostream & os, const Homology & h){
    for(int d = 1; d <= h.dim; d++){
      os << "boundary map : C^" << d << " -> C^"<< d-1 << "\n";
      // matrix
      os << "size = (" << h.boundary_map[d].size() << " * " << h.boundary_map[d][0].size() << ")\n";
      os << h.boundary_map[d] << "\n";
      // Smith normal form
      os << "Smith normal form (size=" <<  h.boundary_map_snf[d].size() << ")\n";
      os << h.boundary_map_snf[d] << "\n";
    }
    // homology group
    os << "Homology group\n";
    for(int d = 0; d <= h.dim; d++){
      os << "H_{" << d << "} = Z^{" << h.Hrank[d] << "}";
      for(int tor : h.Htor[d] ){
	os << " (+) Z/"<< tor << "Z";
	}
      os << "\n";
    }
    return os;
  }
};


#endif
