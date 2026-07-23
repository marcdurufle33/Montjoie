#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Montjoie.hxx"

using namespace Montjoie;

// #include <map>
// #include <string>
// #include <fstream>
// #include <iostream>
// #include <stdlib.h>
 
// using namespace std;

class Configuration
{
public:
  // clear all values
  void Clear();
 
  // load a configuration file
  bool Load(const string& File);
 
  // check if value associated with given key exists
  bool Contains(const string& key) const;
 
  // get value associated with given key
  bool Get(const string& key, string& value, bool go_on=false) const;
  bool Get(const string& key, int&    value, bool go_on=false) const;
  bool Get(const string& key, long&   value, bool go_on=false) const;
  bool Get(const string& key, double& value, bool go_on=false) const;
  bool Get(const string& key, bool&   value, bool go_on=false) const;

  void Display_keys() {for(map<string,string>::iterator it=data.begin();it!=data.end();++it) cout << it->first << endl;};

  // get parameters deduced from config file
  double GethOmega();
  int GetNt();
  int GetLimit();
  VectReal_wp GetOmega();
  VectReal_wp GetTime();
  VectReal_wp GetPs();
 
private:
  // the container
  map<string,string> data;
 
  // remove leading and trailing tabs and spaces
  static string Trim(const string& str);
};
 
 
// ----------------------------------
// method implementations
// ----------------------------------
 
void Configuration::Clear()
{
  data.clear();
}
 
bool Configuration::Load(const string& file)
{
  ifstream inFile;
  inFile.open(file.c_str());
 
  if (!inFile.is_open())
    {
      cout << "Cannot read configuration file " << file << endl;
      return false;
    }

  string line;
  while (getline(inFile,line))
    {
 
      // filter out comments
      if (!line.empty())
        {
          line = line.substr(0,line.find_first_of('#'));
        }
 
      // split line into key and value
      if (!line.empty())
        {
	  int pos = line.find('=');
 
	  if (pos != string::npos)
            {
	      string key     = Trim(line.substr(0, pos));
	      string value   = Trim(line.substr(pos + 1));

	      if (!key.empty() && !value.empty())
                {
		  data[key] = value;
                }
            }
        }
    }
 
  return true;
}
 
bool Configuration::Contains(const string& key) const
{
  return data.find(key) != data.end();
}
 
bool Configuration::Get(const string& key, string& value, bool go_on) const
{
  map<string,string>::const_iterator iter = data.find(key);
 
  if (iter != data.end())
    {
      value = iter->second;
      return true;
    }
  else
    {
      if (!go_on) 
        {
          cerr << "Missing entry " << key << ". Will stop now." << endl;
          abort();
        }
      return false;
    }
}
 
bool Configuration::Get(const string& key, int& value, bool go_on) const
{
  string str;
 
  if (Get(key,str,go_on))
    {
      value = atoi(str.c_str());
      return true;
    }
  else
    {
      return false;
    }
}
 
bool Configuration::Get(const string& key, long& value, bool go_on) const
{
  string str;
 
  if (Get(key,str,go_on))
    {
      value = atol(str.c_str());
      return true;
    }
  else
    {
      return false;
    }
}
 
bool Configuration::Get(const string& key, double& value, bool go_on) const
{
  string str;
 
  if (Get(key,str,go_on))
    {
      value = atof(str.c_str());
      return true;
    }
  else
    {
      return false;
    }
}
 
bool Configuration::Get(const string& key, bool& value, bool go_on) const
{
  string str;
 
  if (Get(key,str,go_on))
    {
      value = ((str == "true") || (str == "1"));
      return true;
    }
  else
    {
      return false;
    }
}
 
string Configuration::Trim(const string& str)
{
  int first = str.find_first_not_of(" \t");
 
  if (first != string::npos)
    {
      int last = str.find_last_not_of(" \t");
 
      return str.substr(first, last - first + 1);
    }
  else
    {
      return "";
    }
}


double Configuration::GethOmega()
{
  int T;
  Get("T", T);
  return 2 * pi_wp / T;
}

int Configuration::GetNt()
{
  int T, dt;
  Get("T", T);
  Get("dt", dt);
  return T / dt;
}

int Configuration::GetLimit()
{
  int Nt = GetNt();
  bool useCausality;
  Get("useCausality", useCausality);
  double limit;
 if (useCausality)
    limit = (Nt + 1) / 2;
  else
    limit = Nt;

 return limit;
}


VectReal_wp Configuration::GetOmega()
{
  int Nt = GetNt();
  int limit = GetLimit();
  double homega = GethOmega();
  VectReal_wp omega(Nt);
  for (int i = 0; i < limit; i++)
    omega(i) = i * homega; 
  for (int i = limit; i < Nt; i++)
    omega(i) = - (i - limit + 1) * homega; 
  return omega;
}

VectReal_wp Configuration::GetPs()
{
  int Nt = GetNt(); 
  VectReal_wp omega = GetOmega();
  double w0, sig, gamma;
  Get("w0", w0);
  Get("sig", sig);
  Get("gamma", gamma);

  VectReal_wp Ps(Nt);
  for (int i = 0; i < Nt; i++)
    Ps(i) = exp(-pow(abs(omega(i))-w0, 2) / (2 * pow(sig, 2) )) * abs(omega(i) / w0);
  return Ps;
 }

VectReal_wp Configuration::GetTime()
{
  int Nt = GetNt();
  int limit = GetLimit();
  double dt;
  Get("dt", dt);
  VectReal_wp time(Nt);
  for (int i = 0; i < Nt; i++)
    time(i) = (i - limit + 1) * dt; 
  return time;
}
