#ifndef _arr1d_h
#define _arr1d_h

#include <array>

namespace mc
{
template<typename T, unsigned int N>
class one_dim_array
{
private:
  std::array<T,N> _data;
public:
  // get element i
  T& operator[](const unsigned int& i)
  {
    return _data[i];
  };
  // get element i
  const T& operator[](const unsigned int& i) const
  {
    return _data[i];
  };
  // size of the array
  const unsigned int size() const
  {
    return N;
  }
  // find the norm
  T norm() const
  {
    T norm2 = 0;
    for (const auto& d: _data)
    {
      norm2 += d*d;
    }
    return std::sqrt(norm2);
  };
  // calculate norm^2 of an array
  T norm2() const
  {
    T norm2 = 0;
    for (const auto& d: _data)
    {
      norm2 += d*d;
    }
    return norm2;
  };
};

} // end of namespace mc


#endif //_arr1d_h
