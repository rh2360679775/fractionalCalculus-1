// c++ 20
#include <iostream>
#include <iomanip> // for setprecision
#include <array>
#include <ranges>
#include <numeric>


// data types for templates
#include "concepts.hpp"

// list of templates 

// print all elements of a list (array)
template <FloatTypesAllowed T, std::size_t N>
inline void printList(const std::array<T, N>& elements, const int prec)
{
  std::cout << std::fixed << std::setprecision( prec );
  for (const auto& el : elements) {
    std::cout << el << std::endl;
  }
}

template <FloatTypesAllowed T, std::size_t N>
inline  T polySum( const std::array<T, N>& as, const T x ){
	T result = T(0);
  T xn = T(1);
	
	for (const auto& a : as ) {
			result += a*xn;
      xn*= x ;
	}
	return result;
}

template <FloatTypesAllowed T, std::size_t N>
inline  T polySumUp( const std::array<T, N>& as, const T x ){
  
	T result = as[0];   
	T xn = T(1);

	for(int n = 1; n < N; n++){
		xn *= x;
		result += as[n]*xn;
	}
	return result;
}

template <FloatTypesAllowed T, std::size_t N>
inline  T horner( const std::array<T, N>& as, const T x ){
	T result = T(0);
	for (const auto& a : as | std::views::reverse) {
			result = result * x + a;
	}
	return result;
}

// special: works for as[20] only 
template <FloatTypesAllowed T, std::size_t N>
inline  T horner45( const std::array<T, N>& as, const T x ){
	T a[4];
	a[0] = as[ 0] + as[ 1]*x + x*x*(as[ 2] + x*(as[ 3] + as[ 4]*x));
	a[1] = as[ 5] + as[ 6]*x + x*x*(as[ 7] + x*(as[ 8] + as[ 9]*x));
	a[2] = as[10] + as[11]*x + x*x*(as[12] + x*(as[13] + as[14]*x));
	a[3] = as[15] + as[16]*x + x*x*(as[17] + x*(as[18] + as[19]*x));
  T x5 = x*x*x*x*x;
	return (a[0] + a[1]*x5 + x5*x5*(a[2] + a[3]*x5));
}

// special: works for as[20] only 
template <FloatTypesAllowed T, std::size_t N>
inline  T horner120( const std::array<T, N>& as, const T x ){
	return       as[ 0] + x* as[ 1] + x*x*(as[2]  + x*(as[ 3] 
           + x*(as[ 4] + x*(as[ 5] +  x*(as[ 6] + x*(as[ 7] 
           + x*(as[ 8] + x*(as[ 9] + 	x*(as[10] + x*(as[11] 
           + x*(as[12] + x*(as[13] +  x*(as[14] + x*(as[15] 
           + x*(as[16] + x*(as[17] +  x*(as[18] + x*as[19])))
           ))))
           ))))
           ))))
           ));
}

// a==b for double / float / complex
template < typename T >
inline bool almostEqual(T a, T b) {
  return (a == 0.0 || b == 0.0) ?
    (abs(a - b) <= std::numeric_limits<T>::epsilon()) :
    ((abs(a - b) <= std::numeric_limits<T>::epsilon() * abs(a)) || (abs(a - b) <= std::numeric_limits<T>::epsilon() * abs(b)));
}


template <FloatTypesAllowed T>
inline  T trapez(T ( *funcptr)(T), const  T a, const  T b, const  int intervals ){
	// funcptr    function pointer to integrand
  // a          lower bound  of integral
  // b          upper bound  of integral 
  // intervals  number of intervals      where intervals > 0
  
  T h = (b - a) / intervals; // interval size
    
  // start with endpoints , so we can treat all values simultaneously
  T total = T(0.5) * (-funcptr(a) + funcptr(b)); // init total = summing up
  
  for (int i = 0; i < intervals; i++) {
    T x = a + h * i;     
    total += funcptr(x); // get function values f(x), a <= f(x) <= b
  }

  total *= h;            // sum times interval size

  return total;
}

template <FloatTypesAllowed T>
inline  T trapez(T ( *funcptr)(T), const  T a, const  T b, const  int intervals , const  double power){
	// integrates \int_a^b dx x^power f(x)
  //
  // power      x**power
  // funcptr    function pointer to integrand
  // a          lower bound  of integral
  // b          upper bound  of integral 
  // intervals  number of intervals  
  
  T h = (b - a) / intervals; // interval size
    
  // start with endpoints , so we can treat all values simultaneously
  T total = T(0.5) * (-pow(a,power)*funcptr(a) + pow(b,power)*funcptr(b)); // init total = summing up
  
  for (int i = 0; i < intervals; i++) {
    T x = a + h * i;     
    total += pow(x,power)*funcptr(x); // get function values f(x), a <= f(x) <= b
  }

  total *= h;            // sum times interval size

  return total;
}


template <FloatTypesAllowed T>  // normalize
inline void shiftAndFactor( T* x,  T* factor) {
  T one  = T(1);
	T zero = T(0);
	while(*x > one ){ *factor *= *x; *x -= one; } 
	while(*x < zero){ *factor /= one + *x; *x += one; } 
}






