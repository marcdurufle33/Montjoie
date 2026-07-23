namespace Montjoie
{


  Vector<string> split(const string& s, char delimiter)
  {
    Vector<string> tokens;
    string token;
    istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
      {
	tokens.PushBack(token);
      }
    return tokens;
  }

  string removeSpaceBegin(string s)
  {
    vector<char> vec;
    bool beg = true;
    for (string::iterator it = s.begin(); it != s.end(); ++it)
      {
	string crt = string(1,*it);
	if (not (crt == " ") or not beg)
	  {
	    vec.push_back(*it);
	    beg = false;
	  }
      }
    string str(vec.begin(),vec.end());
    return str;
  }

  string removeSpaceBeginEnd(string s)
  {
    s = removeSpaceBegin(s);
    reverse(s.begin(), s.end());
    s = removeSpaceBegin(s);
    reverse(s.begin(), s.end());
    return s;
  }

  Vector<string> RealLine(const string line)
  {
    // Read keyword and value for each line of the config file
    Vector<string> ans(2); // contains keyword and value
    char delim = '=';
    size_t current = line.find(delim);
    size_t previous = 0;
    string keyword = line.substr(0,current);
    //keyword.erase(remove_if(keyword.begin(), keyword.end(), ::isspace), keyword.end()); // remove spaces if necessary
    keyword = removeSpaceBeginEnd(keyword);
    previous = current + 1;
    current = line.find(delim, previous);
    string value = line.substr(previous, current - previous);
    //value.erase(remove_if(value.begin(), value.end(), ::isspace), value.end());
    value = removeSpaceBeginEnd(value);
    ans(0) = keyword;
    ans(1) = value;
    return ans;
  }

  Vector<Real_wp> ReadRealVect(const string filename, const int nbVal)
  {
    /*
      Reads nbVal values of a real vector from a file. 
      !!! It only works if the data in the file are in scientific format e.g. 1.000000000000000000e+02
     */

    ifstream f(filename.c_str());
    if (!f.good())
      {
	cout << "Unable to open file " << filename << endl;
	abort();
      }
    string tmp;
    Vector<Real_wp> ans(nbVal);
    
    for (int i=0; i < nbVal; i++)
      {
	f >> tmp;
	Real_wp valR = to_num<Real_wp>(tmp);
	ans(i) = valR;
      }
    f.close();
    return ans;
  }

  Vector<Complex_wp> ReadComplexVect(const string filename, const int nbVal)
  {
    /*
      Reads nbVal values of a complex vector from a file. 
      !!! It only works if the data in the file are in scientific format e.g. (1.000000000000000000e+02+5.000000000000000000e-01j)
     */

    ifstream f(filename.c_str());
    if (!f.good())
      {
	cout << "Unable to open file " << filename << endl;
	abort();
      }
    string tmp;
    vector<char> vec;
    int indexExp;
    Vector<Complex_wp> ans(nbVal);
  
    for (int i=0; i < nbVal; i++)
      {
	f >> tmp;
	vec.clear();
	indexExp = 0;
	for (string::iterator it = tmp.begin(); it != tmp.end(); ++it)
	  {
	    // Put the line in a vector and locates the position of the first exponent to split later
	    string crt = string(1,*it);
	    if (not (crt == "(" or crt == ")"))
	      {
		vec.push_back(*it);
		if (crt == "e" and indexExp == 0)
		  indexExp = vec.size();
	      }
	  }
	// Split into 2 parts containing real and imaginary part
	string str(vec.begin(),vec.end());
	string str2 = str;
	string valS = str2.erase(indexExp+3,vec.size());
	Real_wp valR = to_num<Real_wp>(valS);
	valS = str.erase(0,indexExp+4);
	Real_wp valI = to_num<Real_wp>(valS);
	ans(i) = Complex_wp(valR,valI);
      }
    f.close();
    return ans;
  }

  Matrix<Real_wp> reshapeVector(Vector<Real_wp> vect, int size1, int size2)
  {
    Matrix<Real_wp> ans(size1,size2);
    int cpt = 0;
    for (int i=0; i < size1; i++)
      {
	for (int j=0; j < size2; j++)
	  {
	    ans(i,j) = vect(cpt);
	    cpt++;
	  }
      }
    return ans;
  }

  Matrix<Complex_wp> reshapeVector(Vector<Complex_wp> vect, int size1, int size2)
  {
    Matrix<Complex_wp> ans(size1,size2);
    int cpt = 0;
    for (int i=0; i < size1; i++)
      {
	for (int j=0; j < size2; j++)
	  {
	    ans(i,j) = vect(cpt);
	    cpt++;
	  }
      }
    return ans;
  }

  Matrix<Complex_wp> conjugate(Matrix<Complex_wp> A)
  {
    Matrix<Complex_wp> B(A.GetM(), A.GetN());
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j=0; j < A.GetN(); j++)
	  B(i,j) = conj(A(i,j));
      }
    return B;
  }

  Vector<Matrix<Complex_wp> > conjugate(Vector<Matrix<Complex_wp> > A)
  {
    Vector<Matrix<Complex_wp> >  B(A.GetM());
    for (int l = 0; l < A.GetM(); l++)
      {
	B(l) = Matrix<Complex_wp>(A(l).GetM(), A(l).GetN());
	for (int i = 0; i < A(l).GetM(); i++)
	  {
	    for (int j=0; j < A(l).GetN(); j++)
	      B(l)(i,j) = conj(A(l)(i,j));
	  }
      }
    return B;
  }

  Vector<Vector<Matrix<Complex_wp> > > conjugate(Vector<Vector<Matrix<Complex_wp> > >A)
  {
    Vector<Vector<Matrix<Complex_wp> > > B(A.GetM());
    for (int k = 0; k < A.GetM(); k++)
      {
	B(k).Reallocate(A.GetM());
	for (int l = 0; l < A(k).GetM(); l++)
	  {
	    B(k)(l) = Matrix<Complex_wp>(A(k)(l).GetM(), A(k)(l).GetN());
	    for (int i = 0; i < A(k)(l).GetM(); i++)
	      {
		for (int j=0; j < A(k)(l).GetN(); j++)
		  B(k)(l)(i,j) = conj(A(k)(l)(i,j));
	      }
	  }
      }
    return B;
  }

}
