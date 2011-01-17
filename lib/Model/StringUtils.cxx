#include <StringUtils.h>


std::string SqueezeString(const std::string& str) {
  std::string squeezed(str);
  size_t pos;
  while ((pos = squeezed.find_first_of(' ')) != std::string::npos) {
    squeezed = squeezed.erase(pos, 1);
  }

  return squeezed;
}
