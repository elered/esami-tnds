#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>  // se se vogliono usare algoritmi STL
#include <numeric>    // per inner_product

using namespace std;

// ===============================================================================               
// somma di due vettori : somma componente per componente
// ===============================================================================                                                                                               
template <typename T> vector<T> operator+(const vector<T> &a, const vector<T> &b) {
  
  if ( a.size() != b.size() ) { 
    cout << "Trying to sum vectors with different size, exiting" << endl ; 
    exit(-1); 
  } ; 
  
  vector<T> result(a.size());

  transform(a.begin(), a.end(), b.begin(), result.begin() , plus<T>()); 
  
  return result;
  
}

// ===============================================================================          
// differenza di due vettori componente per componente
// ===============================================================================                                                                                                
template <typename T> vector<T> operator-(const vector<T> &a, const vector<T> &b) {
                                                                                             
  if ( a.size() != b.size() ) {
    cout << "Trying to sum vectors with different size, exiting" << endl ;
    exit(-1);
  } ;

  vector<T> result(a.size());
  
  transform(a.begin(), a.end(), b.begin(), result.begin() , minus<T>()); 
  
  return result;
  
}

// ===============================================================================  
// prodotto scalare tra due vettori 
// ===============================================================================

template <typename T> T operator*(const vector<T> &a, const vector<T> &b) {
                                                                                       
  if ( a.size() != b.size() ) {
    cout << "Trying to sum vectors with different size, exiting" << endl ;
    exit(-1);
  } ;

  T sum = 0 ;

  sum = inner_product(begin(a), end(a), begin(b), 0.0);
  
  return sum;
  
}

// ===============================================================================
// prodotto tra uno scalare e un vettore
// ===============================================================================

template <typename T> vector<T> operator*( T c , const vector<T> &a) {
  
  vector<T> result(a.size());
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = c * a[i];  
  
  return result;
  
}

// ===============================================================================
// prodotto tra un vettore e uno scalare 
// ===============================================================================

template <typename T> vector<T> operator*( const vector<T> &a , T c) {
  
  vector<T> result(a.size());
  
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = c * a[i];  
  
  return result;
  
}

// ===============================================================================
// divisione tra un vettore e uno scalare 
// ===============================================================================

template <typename T> vector<T> operator/( const vector<T> &a , T c) {
  
  vector<T> result(a.size());
  for (int i = 0; i < static_cast<int>(a.size()); i++) result[i] = a[i] / c ;  
  
  return result;
  
}

// ===============================================================================
// somma ad a un vettore b e il risultato viene messo in a
// ===============================================================================

template <typename T> vector<T>& operator+=(vector<T>& a, const vector<T>& b) {
                                                                                     
  if ( a.size() != b.size() ) {
    cout << "Trying to sum vectors with different size, exiting" << endl ;
    exit(-1);
  } ;

  transform(a.begin(), a.end(), b.begin(), a.begin() , plus<T>()); 
  
  return a;
}

// ===============================================================================
// sottrai ad a un vettore b e il risultato viene messo in a
// ===============================================================================

template <typename T> vector<T>& operator-=(vector<T>& a, const vector<T>& b) {
                                                                                        
  if ( a.size() != b.size() ) {
    cout << "Trying to sum vectors with different size, exiting" << endl ;
    exit(-1);
  } ;
  
  transform(a.begin(), a.end(), b.begin(), a.begin() , minus<T>()); 
  
  return a;
}

// ===============================================================================
// Possiamo usare funzioni generiche che agiscono sui vettori
// ===============================================================================
  
template <typename T> void Print( const vector<T> &v ) {
  
  cout << "Printing vector" << endl;
  for ( auto it = v.begin() ; it != v.end() ; it++ ) cout << *it << " " ;
  cout << endl;
  cout << "End of printing vector" << endl;
  
};

template <typename T> double Distanza(vector <T> &v,vector <T> &w){
  T sum = 0;
 
  for(int i = 0; i < v.size(); i++){
    sum += (v[i]-w[i])*(v[i]-w[i]);
  }

  return sqrt(sum);
}