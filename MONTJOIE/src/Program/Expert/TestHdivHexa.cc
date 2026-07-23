#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_NODAL_HDIV
#define MONTJOIE_WITH_NODAL_H1
#define MONTJOIE_WITH_NODAL_DG
#define MONTJOIE_WITH_HP_HDIV

#include "FiniteElement/MontjoieFiniteElement.hxx"

using namespace Montjoie;

// fonction de test d'inclusion des fonctions definis dans test_elt par les fonctions de first_elt
template<class Matrix1>
void CheckInclusion(const Matrix<R3>& val_test, const Matrix<R3>& val_first, Matrix1& mass_first,
                    const VectR3& points3d, VectReal_wp& weights3d, int & nb_not_inside,
                    IVect& num_psi_not_inside, VectReal_wp& dist_max)
{
  // seuil utilise pour savoir si psi = \sum alpha_j phi_j
  Real_wp threshold = 1e-9;
  
  int nb_dof_test = val_test.GetM();
  int nb_dof_first = val_first.GetM();
  int nb_points = points3d.GetM();
  VectReal_wp B(nb_dof_first), C(nb_dof_first);
  C.Fill(0);
  nb_not_inside = 0;
  num_psi_not_inside.Reallocate(nb_dof_test);
  dist_max.Reallocate(nb_dof_test);
  for (int n = 0; n < nb_dof_test; n++)
    {
      // calcul de B_i = (psi, phi_i) pour tout phi_i de la premiere famille
      B.Fill(0);
      for (int i = 0; i < nb_dof_first; i++)
        for (int k = 0; k < nb_points; k++)
          B(i) += weights3d(k)*DotProd(val_test(n, k), val_first(i, k));
      
      // on inverse par la masse de la premiere famille, C est donc le projete de psi sur la premiere famille
      Mlt(mass_first, B, C);
      
      // et on regarde si psi et C ont les memes valeurs sur les points de quadrature
      bool psi_inside = true;
      Real_wp dist = 0;
      for (int k = 0; k < nb_points; k++)
        {
          R3 val_psi = val_test(n, k);
          R3 val_C;
          for (int j = 0; j < nb_dof_first; j++)
            val_C += val_first(j, k)*C(j);
          
          if (val_C.Distance(val_psi) > threshold)
            {
              //DISP(val_psi); DISP(val_C);
              dist = max(dist, val_C.Distance(val_psi));
              psi_inside = false;
            }
        }
      
      if (!psi_inside)
        {
          num_psi_not_inside(nb_not_inside) = n;
          dist_max(nb_not_inside) = dist;
          nb_not_inside++;
        }
    }
  
  if (nb_not_inside > 0)
    {
      num_psi_not_inside.Resize(nb_not_inside);
      dist_max.Resize(nb_not_inside);
    }
  else
    {
      num_psi_not_inside.Clear();
      dist_max.Clear();
    }
}


int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc < 3)
    {
      cout << "Veuillez preciser l'ordre d'approximation des fonctions test" 
           << " et de l'espace de la premiere famille : " << endl;
      
      cout << "./test.x 2 1" << endl;
      abort();
    }
  
  bool test_contains = true;
  int r = atoi(argv[1]);
  int rf = atoi(argv[2]);
  if (argc > 3)
    test_contains = false;
  
  /*******************************************************************
   * Construction des elements finis test, premiere famille, optimal *
   *******************************************************************/
  
  //typedef Matrix<Real_wp, Symmetric, RowSymPacked> MatrixMass;
  //bool sym_mass = true;
  typedef Matrix<Real_wp, General, RowMajor> MatrixMass;
  bool sym_mass = false;
  
  // on construit d'abord WedgeHcurlOther
  HexahedronHdivOptimalHpFirstFamily test_elt;
  //HexahedronHcurlHpFirstFamily test_elt;
  //HexahedronHcurlFirstFamily test_elt;
  //test_elt.type_basis = test_elt.NIGAM_PHILLIPS_2;
  //test_elt.type_basis = test_elt.NIGAM_PHILLIPS_1;
  test_elt.type_basis = test_elt.SUPER_OPTIMAL;
  test_elt.ConstructFiniteElement(r);

  // on recupere une formule de quadrature d'ordre max(r, rf) + 1
  VectR3 points3d;
  VectReal_wp weights3d;
	
  points3d = test_elt.PointsND();
  weights3d = test_elt.WeightsND();
  
  int nq = test_elt.GetNbPointsQuadratureInside();
  points3d.Resize(nq);
  weights3d.Resize(nq);
  
  // nombre de points de quadrature
  int nb_points = points3d.GetM();
  // on calcule la matrice de masse du nouvel element
  MatrixMass mass_test(test_elt.GetNbDof(), test_elt.GetNbDof());
  Matrix<R3> val_test(test_elt.GetNbDof(), nb_points);
  mass_test.Fill(0);
  VectR3 phi;
  for (int k = 0; k < nb_points; k++)
    {
      test_elt.ComputeValuesPhiRef(points3d(k), phi);
      for (int i = 0; i < mass_test.GetM(); i++)
        {
          int j0 = 0;
          if (sym_mass)
            j0 = i;
          
          for (int j = j0; j < mass_test.GetM(); j++)
            mass_test(i, j) += weights3d(k)*DotProd(phi(i), phi(j));
        }
      
      for (int i = 0; i < test_elt.GetNbDof(); i++)
        val_test(i, k) = phi(i);
    }
  
  // on essaie de l'inverser
  cout << endl;
  if (test_contains)
    {
      cout << "On tente d'inverser la matrice de masse des fonctions tests " << endl;
      GetInverse(mass_test);
      cout << "Inversion reussie " << endl;
    }
  else
    {
      // on regarde les valeurs propres et le rang de la matrice de masse
      /*
      mass_test.WriteText(string("Mh") + to_str(r) + ".dat");
      VectReal_wp lambda; Real_wp threshold = 1e-12;
      cout << "Entrez le seuil pour decider si une valeur propre est nulle " << endl;
      cin >> threshold;
      //GetEigenvalues(mass_test, lambda);
      Sort(lambda);
      cout << "Valeurs propres : " << endl << lambda << endl;
      int rank = 0;
      for (int i = 0; i < mass_test.GetM(); i++)
        if (abs(lambda(i)) > threshold)
          rank++;
      
	  cout << " Rang de la matrice de masse : "  << rank << endl; */
    }
  
  // ensuite on construit une base quasi-ortho de la premiere famille
  HexahedronHdivFirstFamily first_elt;
  //TetrahedronHcurlFirstFamily first_elt;
  first_elt.ConstructOrthogonalBasis(rf);
  //HexahedronHcurlHpFirstFamily first_elt;
  //first_elt.ConstructFiniteElement(rf);
  
  // dont on calcule l'inverse de la matrice de masse
  MatrixMass mass_first(first_elt.GetNbDof(), first_elt.GetNbDof());
  Matrix<R3> val_first(first_elt.GetNbDof(), nb_points);
  mass_first.Fill(0);
  for (int k = 0; k < points3d.GetM(); k++)
    {
      first_elt.ComputeValuesPhiOrthoRef(points3d(k), phi);
      //first_elt.ComputeValuesPhiRef(points3d(k), phi);
      for (int i = 0; i < mass_first.GetM(); i++)
        {
          int j0 = 0;
          if (sym_mass)
            j0 = i;
          
          for (int j = j0; j < mass_first.GetM(); j++)
            mass_first(i, j) += weights3d(k)*DotProd(phi(i), phi(j));
        }
      
      for (int i = 0; i < first_elt.GetNbDof(); i++)
        val_first(i, k) = phi(i);
    }
  
  // on essaie de l'inverser
  GetInverse(mass_first);
  
  // ensuite on construit une base quasi-ortho de la premiere famille optimale
  HexahedronHdivOptimalFirstFamily opt_elt;
  opt_elt.ConstructOrthogonalBasis(rf, true);
  
  // dont on calcule l'inverse de la matrice de masse
  MatrixMass mass_opt(opt_elt.GetNbDof(), opt_elt.GetNbDof());
  Matrix<R3> val_opt(opt_elt.GetNbDof(), nb_points);
  mass_opt.Fill(0);
  for (int k = 0; k < points3d.GetM(); k++)
    {
      opt_elt.ComputeValuesPhiOrthoRef(points3d(k), phi);
      for (int i = 0; i < mass_opt.GetM(); i++)
        {
          int j0 = 0;
          if (sym_mass)
            j0 = i;
          
          for (int j = j0; j < mass_opt.GetM(); j++)
            mass_opt(i, j) += weights3d(k)*DotProd(phi(i), phi(j));
        }
      
      for (int i = 0; i < opt_elt.GetNbDof(); i++)
        val_opt(i, k) = phi(i);
    }
  
  //mass_opt.Write("Mh.dat");
  // on essaie de l'inverser
  GetInverse(mass_opt);
  
  /**********************************************************
   * On teste l'inclusion Test inclus dans Premiere Famille *
   **********************************************************/
  
  cout << "L'espace test contient " << test_elt.GetNbDof() << " ddls" << endl;
  cout << "La premiere famille contient " << first_elt.GetNbDof() << " ddls" << endl;
  cout << "La premiere famille optimale contient " << opt_elt.GetNbDof() << " ddls" << endl;
  
  cout << endl << endl;
  cout << "On regarde si les fonctions sont incluses dans la premiere famille d'ordre " << rf << endl << endl;
  
  // ensuite on essaie d'inclure l'espace test dans la premiere famille
  int nb_not_inside;
  IVect num_psi_not_inside;
  VectReal_wp dist_max;
  CheckInclusion(val_test, val_first, mass_first, points3d, weights3d, 
                 nb_not_inside, num_psi_not_inside, dist_max);
  
  if (nb_not_inside == 0)
    {
      if (first_elt.GetNbDof() == test_elt.GetNbDof())
        {
          cout << "Ces fonctions tests generent exactement la premiere famille d'ordre " << rf << endl << endl;
          return FinalizeMontjoie();
        }
      else
        {
          cout << "Ces fonctions tests sont incluses dans la premiere famille d'ordre " << rf << endl;          
        }
    }
  else
    {
      cout << "Ces fonctions tests ne sont pas incluses dans la premiere famille d'ordre " << rf << endl;
      cout << "On compte " << nb_not_inside << " fonctions qui n'y sont pas, dont voici les numeros : " << endl;
      cout << num_psi_not_inside << endl;
      cout << "Ecart maximal : " << endl << dist_max << endl;
    }
  
  cout << endl << endl;

  /*******************************************************************
   * On teste l'inclusion Test inclus dans Premiere Famille Optimale *
   *******************************************************************/
  
  cout << "On regarde si les fonctions sont incluses dans la premiere famille optimale d'ordre " << rf << endl << endl;
  
  // ensuite on essaie d'inclure l'espace test dans la premiere famille optimale
  CheckInclusion(val_test, val_opt, mass_opt, points3d, weights3d, 
                 nb_not_inside, num_psi_not_inside, dist_max);
  
  if (nb_not_inside == 0)
    {
      if (opt_elt.GetNbDof() == test_elt.GetNbDof())
        {
          cout << "Ces fonctions tests generent exactement la premiere famille optimale d'ordre " << rf << endl;
          return FinalizeMontjoie();
        }
      else
        {
          cout << "Ces fonctions tests sont incluses dans la premiere famille optimale d'ordre " << rf << endl;          
        }
    }
  else
    {
      cout << "Ces fonctions tests ne sont pas incluses dans la premiere famille optimale d'ordre " << rf << endl;
      cout << "On compte " << nb_not_inside << " fonctions qui n'y sont pas, dont voici les numeros : " << endl;
      cout << num_psi_not_inside << endl;
      cout << "Ecart maximal : " << endl << dist_max << endl;
    }
  
  if (!test_contains)
    return FinalizeMontjoie();
  
  /**********************************************************
   * On teste l'inclusion Premiere Famille inclus dans Test *
   **********************************************************/
  
  cout << endl << endl;
  cout << "On regarde si les fonctions contiennent la premiere famille d'ordre " << rf << endl << endl;
  
  // ensuite on essaie d'inclure la premiere famille dans l'espace test
  CheckInclusion(val_first, val_test, mass_test, points3d, weights3d, 
                 nb_not_inside, num_psi_not_inside, dist_max);
  
  if (nb_not_inside == 0)
    {
      cout << "Ces fonctions tests contiennent la premiere famille d'ordre " << rf << endl;          
    }
  else
    {
      cout << "Ces fonctions tests ne contiennent pas la premiere famille d'ordre " << rf << endl;
      cout << "On compte " << nb_not_inside << " fonctions de la premiere famille qui n'y sont pas, dont voici les numeros : " << endl;
      cout << num_psi_not_inside << endl;
      cout << "Ecart maximal : " << endl << dist_max << endl;
    }

  /*******************************************************************
   * On teste l'inclusion Test inclus dans Premiere Famille Optimale *
   *******************************************************************/
  
  cout << endl << endl;
  cout << "On regarde si les fonctions contiennent la premiere famille optimale d'ordre " << rf << endl << endl;
  
  // ensuite on essaie d'inclure la premiere famille optimale dans l'espace test
  CheckInclusion(val_opt, val_test, mass_test, points3d, weights3d, 
                 nb_not_inside, num_psi_not_inside, dist_max);
  
  if (nb_not_inside == 0)
    {
      cout << "Ces fonctions tests contiennent la premiere famille optimale d'ordre " << rf << endl;          
    }
  else
    {
      cout << "Ces fonctions tests ne contiennent pas la premiere famille optimale d'ordre " << rf << endl;
      cout << "On compte " << nb_not_inside << " fonctions de la premiere famille optimale qui n'y sont pas, dont voici les numeros : " << endl;
      cout << num_psi_not_inside << endl;
      cout << "Ecart maximal : " << endl << dist_max << endl;
    }
  
  return FinalizeMontjoie();
}
