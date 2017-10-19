#ifndef scatter_h
#define scatter_h

#include <iostream>
#include <array>

#include "utility.h"

namespace mc
{

class particle;

class scatter
{
private:

protected:

public:

	virtual mc::t_float ff_time() const = 0; // get random free flight time
	virtual void update_state(mc::particle* p) = 0; // update the final state of the particle
}; // end class scatter

} // end namespace mc

#endif // scatter_h
