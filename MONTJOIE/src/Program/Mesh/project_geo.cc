#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{

  if (argc < 3)
    {
      cout << "Usage ./project_geo source.geo target.geo" << endl;
      abort();
    }

  string ligne;
  
  ifstream file_geo(argv[1]);
  if (!file_geo.is_open())
    {
      cout << "File " << argv[1] << "does not exist ?" << endl;
      abort();
    }
  
  ofstream file_geo_out(argv[2]);
  file_geo_out.precision(15);
  
  Real_wp a, b, c, d, epsilon;
  cout << "Entrez les coefficients du plan a, b, c, d ou les points doivent etre projetes" << endl;
  cin >> a >> b >> c >> d;
  cout << "Entrez le seuil epsilon" << endl;
  cin >> epsilon;

  R3 normale(a, b, c);
  if (Norm2(normale) == 0)
    {
      cout << "Normale must be non-null " << endl;
      abort();
    }
  
  Real_wp coef = 1.0/Norm2(normale);
  Mlt(coef, normale);
  d *= coef;
  
  while (!file_geo.eof())
    {
      getline(file_geo, ligne);
      
      if (file_geo.good())
        {
          if (ligne.substr(0, 5) == string("Point"))
            {
              int pos = ligne.find(')');
              int n = to_num<int>(ligne.substr(6, pos-6));
              Real_wp x, y, z;
              int posA = ligne.find('{');
              pos = ligne.find(',', posA);
              x = to_num<Real_wp>(ligne.substr(posA+1, pos-posA-1));
              posA = pos; pos = ligne.find(',', posA+1);
              y = to_num<Real_wp>(ligne.substr(posA+1, pos-posA-1));
              posA = pos; pos = ligne.find(',', posA+1);
              z = to_num<Real_wp>(ligne.substr(posA+1, pos-posA-1));
              
              R3 point(x, y, z);
              Real_wp scal = DotProd(point, normale) + d;
              if (abs(scal) <= epsilon)
                { 
                  Add(-scal, normale, point);
                  DISP(scal); DISP(DotProd(point, normale)+d);
                }

              file_geo_out << "Point(" << n << ") = {" << point(0) << ", " << point(1) << ", "
                           << point(2) << ligne.substr(pos, ligne.size()-pos) << '\n';              
            }
          else
            file_geo_out << ligne << '\n';
        }
    }

  file_geo.close();
  file_geo_out.close();
  
  return 0;
}
