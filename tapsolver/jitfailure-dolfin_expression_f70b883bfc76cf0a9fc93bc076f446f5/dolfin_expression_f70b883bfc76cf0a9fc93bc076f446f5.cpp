
// Based on https://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
    #ifdef __GNUC__
        #define DLL_EXPORT __attribute__ ((dllexport))
    #else
        #define DLL_EXPORT __declspec(dllexport)
    #endif
#else
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
#endif

#include <dolfin/function/Expression.h>
#include <dolfin/math/basic.h>
#include <Eigen/Dense>


// cmath functions
using std::cos;
using std::sin;
using std::tan;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cosh;
using std::sinh;
using std::tanh;
using std::exp;
using std::frexp;
using std::ldexp;
using std::log;
using std::log10;
using std::modf;
using std::pow;
using std::sqrt;
using std::ceil;
using std::fabs;
using std::floor;
using std::fmod;
using std::max;
using std::min;

const double pi = DOLFIN_PI;


namespace dolfin
{
  class dolfin_expression_f70b883bfc76cf0a9fc93bc076f446f5 : public Expression
  {
     public:
       

       dolfin_expression_f70b883bfc76cf0a9fc93bc076f446f5()
       {
            
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          values[0] = 0 <= x[0] && x[0] < 0.495049504950495 ? TAPobject_data.reactor.voids[0] : 0 ^ 0.495049504950495 <= x[0] && x[0] < 0.5049504950495048 ? TAPobject_data.reactor.voids[1] : 0 ^ 0.5049504950495048 <= x[0] && x[0] < 0.9999999999999998 ? TAPobject_data.reactor.voids[2] : 0;

       }

       void set_property(std::string name, double _value) override
       {

       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {

       throw std::runtime_error("No such property");
       return 0.0;
       }

       void set_generic_function(std::string name, std::shared_ptr<dolfin::GenericFunction> _value) override
       {

       throw std::runtime_error("No such property");
       }

       std::shared_ptr<dolfin::GenericFunction> get_generic_function(std::string name) const override
       {

       throw std::runtime_error("No such property");
       }

  };
}

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_f70b883bfc76cf0a9fc93bc076f446f5()
{
  return new dolfin::dolfin_expression_f70b883bfc76cf0a9fc93bc076f446f5;
}
