#include <stdexcept>
#include <string>
#include "AutoDiff.h"

// Error handling of the automatic differentiation classes.
// These are kept out of the headers so that the (inlined) arithmetic only
// carries a function call for the rare case that something is wrong.

namespace AutoDiff {

  namespace Impl {

    void throw_size_mismatch(const int sz1, const int sz2) {
      std::string error_message = "AutoDiff::Fad::DFad - combining two variables with different numbers of derivative components ("
        + std::to_string(sz1) + " and " + std::to_string(sz2)
        + "); all variables must be either constants or derived from the same set of independent variables.\n";
      throw std::logic_error( error_message );
    }

    void throw_bad_index(const int i, const int sz) {
      std::string error_message = "AutoDiff::Fad::DFad - independent variable index " + std::to_string(i)
        + " is out of range; it must be between 0 and " + std::to_string(sz-1)
        + " for " + std::to_string(sz) + " independent variables.\n";
      throw std::logic_error( error_message );
    }

  } // namespace Impl

} // namespace AutoDiff
