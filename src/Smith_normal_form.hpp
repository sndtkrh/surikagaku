#ifndef SMITHNORMALFORM
#define SMITHNORMALFORM
#include <vector>
#include "Matrix.hpp"
/*
 * 整数係数行列のSmith標準形を求める
 * 単体複体のホモロジーを求めるときなどに使用する
 */

// 拡張Euclid互除法
// ax + by = gcd(a,b);
long long extgcd(long long a, long long b, long long &x, long long &y){
  long long g = a; x = 1; y = 0;
  if(b != 0){
    g = extgcd( b, a % b, y, x );
    y -= (a / b) * x;
  }
  return g;
}

void swap_row(Matrix &M, int s, int t){ // 列を交換
  if( s == t ) return;
  for(int j = 0; j < M[s].size(); j++) std::swap( M[s][j], M[t][j] );
}

std::vector<long long> compute_snf(Matrix M){
  int n = M.size();
  int m = M[0].size();
  std::vector<long long> diag;
  int t = 0;
  for(int j = 0; j < m; j++){
    // find pivot
    bool pivot = false;
    for(int i = t; i < n; i++){
      if( M[i][j] != 0 ){
	swap_row(M, t, i );
	pivot = true;
	break;
      }
    }
    if( ! pivot ) continue; // no pivot
    
    // (t,j) is pivot
    // 第j列のなかで(t,j)が他の全てを割り切るようにする
    for(int k = t + 1; k < n; k++){
      if( M[k][j] % M[t][j] != 0 ){
	long long p = 0, q = 0;
	long long g = extgcd( M[t][j], M[k][j], p, q );
	long long u = M[t][j] / g, v = M[k][j] / g;
	for(int l = j; l < m; l++){
	  long long a = M[t][l];
	  long long b = M[k][l];
	  M[t][l] = a * p + b * q;
	  M[k][l] = -a * v + b * u;
	}
      }
    }
    // 第j列を(t,j)を残して残りは0にする
    for(int k = t + 1; k < n; k++){
      if( M[k][j] != 0 ){
	int d = M[k][j] / M[t][j];
	for(int l = j; l < m; l++)
	  M[k][l] -= d * M[t][l];
      }
    }
    // 第t行のなかで(t,j)が他の全てを割り切るようにする
    for(int k = j + 1; k < m; k++){
      if( M[t][k] % M[t][j] != 0 ){
	long long p = 0, q = 0;
	long long g = extgcd( M[t][j], M[t][k], p, q );
	long long u = M[t][j] / g, v = M[t][k] / g;
	for(int l = t; l < n; l++){
	  long long a = M[l][t];
	  long long b = M[l][k];
	  M[l][t] = a * p + b * q;
	  M[k][k] = -a * v + b * u;
	}
      }
      // 第t行を(t,j)を残して残りは0にする
      for(int k = j + 1; k < m; k++){
	if( M[t][k] != 0 ){
	  int d = M[t][k] / M[t][j];
	  for(int l = t; l < n; l++)
	    M[l][k] -= d * M[l][j];
	}
      }
    }
    diag.push_back( M[t][j] );
    t++;
  }
  // 対角成分に関して diag[i] | diag[i+1] for any i が成立するように調整
  for(int i = 0; i < diag.size() - 1; i++){
    if( diag[i+1] % diag[i] != 0 ){
      int g = std::__gcd( diag[i], diag[i+1] ); // gcd__(x,y) はgcc拡張機能
      int p = diag[i];
      diag[i] = g;
      diag[i+1] *= -p / g;
    }
    diag[i] = abs( diag[i] ); // positive
  }
  diag[t-1] = abs( diag[t-1] );
  return diag;
}

#endif
