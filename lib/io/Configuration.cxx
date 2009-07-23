#include "Configuration.h"

#include <sstream>
#include <cstdio>


Configuration
::Configuration() {

}


Configuration
::~Configuration() {
  
}


void
Configuration
::Parse(std::string fileName) {
  FILE* fp = fopen(fileName.c_str(), "r");
  if (fp == NULL)
    return;

  std::string currentSection;
  char line[1024];
  while (fgets(line, 1023, fp)) {

    if (strlen(line) == 1) {
      // Empty line
      continue;
    } else if (line[0] == ';') {
      // Comment
      continue;
    } else if (line[0] == '[') {
      // Start new section
      char* rightBracket = strstr(line, "]");
      *rightBracket = '\0';

      SectionType newStruct;
      currentSection = std::string(line+1);
      m_SectionMap[currentSection] = newStruct;
    } else {
      // New entry

      // Zap newline from the string
      line[strlen(line)-1] = '\0';

      char* equalSign = strstr(line, "=");
      *equalSign = '\0';

      std::string key = std::string(line);
      std::string val = std::string(equalSign+1);

      m_SectionMap[currentSection][key] = val;
    }
  }

  fclose(fp);
}


void
Configuration
::Write(std::ostream& os) {

  // Iterate over sections
  for (MapIteratorType section = m_SectionMap.begin();
       section != m_SectionMap.end(); section++) {

    // Print section name
    os << "[" << section->first << "]" << std::endl;

    // Iterate over key/value pairs.
    for (SectionIteratorType keyValue = section->second.begin();
	 keyValue != section->second.end(); keyValue++) {
      os << keyValue->first << "=" << keyValue->second << std::endl;
    }
    os << std::endl;
  }
  
}


void
Configuration
::SetValue(std::string section, std::string key, std::string value) {
  m_SectionMap[section][key] = value;
}


std::string
Configuration
::GetValue(std::string section, std::string key) {
  return m_SectionMap[section][key];
}


void
Configuration
::SetValueFromBool(std::string section, std::string key, bool value) {
  std::stringstream stream;
  stream << value;
  m_SectionMap[section][key] = stream.str();
}


bool
Configuration
::GetValueAsBool(std::string section, std::string key, bool defaultValue) {
  std::string strValue = this->GetValue(section, key);
  bool value = defaultValue;
  if (strValue == "true")
    value = true;
  else if (strValue == "false")
    value = false;

  return value;
}


void
Configuration
::SetValueFromInt(std::string section, std::string key, int value) {
  std::stringstream stream;
  stream << value;
  m_SectionMap[section][key] = stream.str();
}


int
Configuration
::GetValueAsInt(std::string section, std::string key, int defaultValue) {
  std::string strValue = this->GetValue(section, key);
  int value = defaultValue;
  sscanf(strValue.c_str(), "%d", &value);

  return value;
}


void
Configuration
::SetValueFromDouble(std::string section, std::string key, double value) {
  std::stringstream stream;
  stream << value;
  m_SectionMap[section][key] = stream.str();
}


double
Configuration
::GetValueAsDouble(std::string section, std::string key, double defaultValue) {
  std::string strValue = this->GetValue(section, key);
  double value = defaultValue;
  sscanf(strValue.c_str(), "%lf", &value);

  return value;
}


void
Configuration
::SetValueFromFloat(std::string section, std::string key, float value) {
  std::stringstream stream;
  stream << value;
  m_SectionMap[section][key] = stream.str();
}


float
Configuration
::GetValueAsFloat(std::string section, std::string key, float defaultValue) {
  std::string strValue = this->GetValue(section, key);
  float value = defaultValue;
  sscanf(strValue.c_str(), "%f", &value);

  return value;
}


void
Configuration
::SetValueFromDoubleArray(std::string section, std::string key, double* values, unsigned int length) {
  std::stringstream stream;
  for (unsigned int i = 0; i < length; i++) {
    stream << values[i] << " ";
  }

  m_SectionMap[section][key] = stream.str();
}



void
Configuration
  ::GetValueAsDoubleArray(std::string section, std::string key, double* values, unsigned int length) {
  std::string value = GetValue(section, key);
  std::stringstream stream(value);
  for (unsigned int i = 0; i < length; i++) {
    stream >> values[i];
  }

}
