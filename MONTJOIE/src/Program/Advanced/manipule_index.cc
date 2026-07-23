#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_TWO_DIM

#include "Output/MontjoieOutput.hxx"

using namespace Montjoie;

template<class T>
void ReadParameters(Vector<T>& param)
{
  string chaine;
  // cin>>flush;
  // getline(cin, chaine ,'\n');
  getline(cin, chaine ,'\n');
  istringstream stream_data(chaine);
  
  param.ReadText(stream_data);
}

void ReadParameters(string& param)
{
  getline(cin, param ,'\n');
}

int main(int argc, char** argv)
{
  if (argc < 2)
    {
      cout << "Entrez le nom du ficier .elb a modifier" << endl;
      abort();
    }
  
  string data_file(argv[1]);
  string output_file("toto.elb");
  
  Vector<Vector<Real_wp> > nu;
  ReadMeshData(nu, data_file);
  
  bool test_loop = true;
  IVect num;
  while (test_loop)
    {
      cout << "1- Normaliser le champ " << endl;
      cout << "2- Sauver " << endl;
      cout << "3- Quitter " << endl;

      ReadParameters(num); int choix = num(0);
      
      switch (choix)
	{
	case 1 :
	  {
	    Real_wp vmax = 0; bool negative_data = false;
	    for (int i = 0; i < nu.GetM(); i++)
	      for (int j = 0; j < nu(i).GetM(); j++)
		{
		  if (nu(i)(j) < -vmax)
		    {
		      vmax = abs(nu(i)(j));
		      negative_data = true;
		    }
		  else if (nu(i)(j) > vmax)
		    {
		      vmax = nu(i)(j);
		      negative_data = false;
		    }
		}
	    
	    Real_wp coef = 1.0/vmax;
	    if (negative_data)
	      coef = -coef;
	    
	    nu *= coef;
	  }
	  break;
	case 2 :
	  WriteMeshData(nu, output_file);
	  break;
	case 3 :
	  test_loop = false;
	  break;
	}
    }
  
  return 0;
}
