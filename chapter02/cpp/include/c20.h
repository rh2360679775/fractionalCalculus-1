// c++20 version
#include <array>
#include <string>

// intended data type 
using FloatType = double;

// precision achieved e-16
const int prec = 32;

// id 
const std::string name = "seriesDouble"; 

// order k of the polynomial => array-size = k+1 = 20
const std::array<FloatType, 20> as = {
    1.0E0, 
    -5.7721566490112839913612516377983E-1, 
    9.8905599526004546488561968743339E-1, 
    -9.0747907205154988315825417112244E-1, 
    9.8172796446769040219142487649730E-1, 
    -9.8199282456119311153012614684595E-1, 
    9.9312179432071013302064800481304E-1, 
    -9.9576722393326028337686469083995E-1, 
    9.9662498968300581987064380393876E-1, 
    -9.9193816665161287149311120508649E-1, 
    9.7319991296155473144438879723251E-1, 
    -9.2254446683703155299124831330242E-1, 
    8.1800576403405562896572614793714E-1, 
    -6.5063455752880833774388636414585E-1, 
    4.4314628622691644741140621352038E-1, 
    -2.4602120185468968949442836220557E-1, 
    1.0527746780281570888539551947701E-1, 
    -3.2239556919065630956147031443838E-2, 
    6.2446860911293176575928532278774E-3, 
    -5.72125609583894452654455305E-4
};	


