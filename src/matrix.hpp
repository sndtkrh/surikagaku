#ifndef MATRIX
#define MATRIX
#include <vector>
#include <iomanip>
#include <iostream>

typedef std::vector< std::vector<long long> > Matrix;
std::ostream& operator<< (std::ostream & os, const Matrix & m){
  for(int i = 0; i < m.size(); i++){
    for(int j = 0; j < m[i].size(); j++){
      os << std::setw(4) << std::right;
      os << m[i][j];
      if( j == m[i].size() - 1 ){
	os << "\n";
      }else{
	os << " ";
      }
    }
  }
  return os;
}
std::ostream& operator<<(std::ostream & os, const std::vector<long long> & v){
  for(int i = 0; i < v.size(); i++){
    os << v[i];
    if( i == v.size() - 1 ){
      os << "\n";
    }else{
      os << " ";
    }
  }
  return os;
}


#endif
