#ifndef STRING_H
#define STRING_H

#include <string>
#include <sstream>

namespace patch
{
 template < typename T > std::string to_string( const T& n )
 {
  std::ostringstream stm ;
  stm << n ;
  return stm.str() ;
 }
}

#endif
