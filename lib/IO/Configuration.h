#ifndef __CONFIGURATION_H_
#define __CONFIGURATION_H_


#include <string>
#include <map>

/** 
 * Loads and parses a configuration file and provides methods for accessing
 * settings values. Also writes configuration files. The format of a 
 * configuration file follows the INI de facto standard
 * (http://en.wikipedia.org/wiki/Ini_file). Keys must contain no spaces,
 * but values can.
 */

class Configuration {
 public:
  /** Useful typedefs. */
  typedef std::map<std::string, std::string> SectionType;
  typedef SectionType::iterator              SectionIteratorType;

  typedef std::map<std::string, SectionType> MapType;
  typedef MapType::iterator                  MapIteratorType;


  Configuration();
  ~Configuration();

  /** Parse a file with the given name. */
  void Parse(std::string fileName);

  /** Write out the contents of the configuration to an output stream. */
  void Write(std::ostream& os);

  void SetValue(std::string section, std::string key, std::string value);
  std::string GetValue(std::string section, std::string key);

  void SetValueFromBool(std::string section, std::string key, bool value);
  bool GetValueAsBool(std::string section, std::string key, bool defaultValue=false);

  void SetValueFromInt(std::string section, std::string key, int value);
  int GetValueAsInt(std::string section, std::string key, int defaultValue=0);

  void SetValueFromFloat(std::string section, std::string key, float value);
  float GetValueAsFloat(std::string section, std::string key, float defaultValue=0.0f);

  void SetValueFromDouble(std::string section, std::string key, double value);
  double GetValueAsDouble(std::string section, std::string key, double defaultValue=0.0);

  void SetValueFromDoubleArray(std::string section, std::string key, double values[2], unsigned int length);
  void GetValueAsDoubleArray(std::string section, std::string key, double values[2], unsigned int length);


 protected:


 private:
  MapType m_SectionMap;

};

// _CONFIGURATION_
#endif
