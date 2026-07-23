#include "Elliptic/Helmholtz/MontjoieHelmholtz.hxx"

using namespace Montjoie;

class InputModal : public InputDataProblem_Base
{
  int nb_modes;
  Real_wp L;
  Real_wp r; R2 pt;
public:
  InputModal() { nb_modes = 2; L = 1.0; r = 1.0; }
  
  inline int GetNbModes() const { return nb_modes; }
  inline Real_wp GetLengthSuperCell() const { return L; }
  inline Real_wp GetRadiusSource() const { return r; }
  inline const R2& GetCenterSource() const { return pt; }
  
  void SetInputData(const string& nom, const Vector<string>& param)
  {
    if (nom == "FourierModes")
      {
        nb_modes = to_num<int>(param(0));
        L = to_num<Real_wp>(param(1));
      }
    else if (nom == "SourceGaussian")
      {
        pt(0) = to_num<Real_wp>(param(0));
        pt(1) = to_num<Real_wp>(param(1));
        r = to_num<Real_wp>(param(2));
      }
  }

  
};

template<class TypeEquation>
void RunAll(EllipticProblem<TypeEquation>& Vars,
            string& file_name_data, string& type_element, string& type_equation)
{
  InputModal var_mode;
  
  Vector<string> lines_data_file;
  ReadLinesFile(file_name_data, lines_data_file);

  Vars.SetTypeEquation(type_equation);
  Vars.InitIndices(PhysicalConstant::nb_max_indices);

  ReadInputFile(lines_data_file, Vars);
  ReadInputFile(lines_data_file, var_mode);
  int nb_modes = var_mode.GetNbModes();
  Real_wp L = var_mode.GetLengthSuperCell();

  All_LinearSolver* glob_solver;
  glob_solver = Vars.GetNewLinearSolver();
  ReadInputFile(lines_data_file, *glob_solver);

  Vars.ComputeMeshAndFiniteElement(type_element, false);

  Vars.PerformOtherInitializations();

  Vars.ComputeMassMatrix();

  Complex_wp m_omega2;
  Vars.GetMomega2(m_omega2);
  int N = Vars.GetNbDof(); DISP(N);
  
  VectComplex_wp mass_rho(N), mass_mu(N);
  mass_rho.Zero(); mass_mu.Zero();
  int nb_elt = Vars.mesh.GetNbElt();
  TinyVector<Complex_wp, 2> v0;
  TinyMatrix<Complex_wp, Symmetric, 2, 2> tmp1, tmp2;
  for (int i = 0; i < nb_elt; i++)
    {
      int i1 = i - Vars.mesh.GetNbElt() + Vars.GetNbEltPML();
      const ElementReference<Dimension2, 1>& Fb = Vars.GetReferenceElementH1(i);
      bool affine = Vars.mesh.IsElementAffine(i);
      int nb_pts_quad = Fb.GetNbPointsQuadratureInside();
      int ref = Vars.mesh.Element(i).GetReference();
      IVect Nodle = Vars.GetDofNumberOnElement(i);
      for (int j = 0; j < nb_pts_quad; j++)
        {
          int k = Nodle(j);
          Real_wp jacob = Vars.GetWeightedJacobian(i, j, affine, Fb.GetGeometricElement());
          Complex_wp rho = Vars.ref_rho(ref).GetConstant();
          Complex_wp mu = Vars.ref_mu(ref).GetConstant()(0, 0);
          if (Vars.InsidePML(i))
            Vars.ModifyPMLCoefficient(rho, mu, tmp1, tmp2, v0, i, i1, j);
          
          mass_rho(k) += m_omega2*rho*jacob;
          mass_mu(k) += mu*jacob;
        }
    }

  DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse> mat_stiff, Ah;
  GlobalGenericMatrix<Complex_wp> nat_mat;
  nat_mat.SetCoefMass(0.0);
  nat_mat.SetCoefStiffness(1.0);
  nat_mat.SetCoefDamping(0.0);

  Vars.AddMatrixWithBC(mat_stiff, nat_mat);

  VectComplex_wp x(N); x.Zero(); VectR2 s;
  for (int m = 0; m < nb_modes; m++)
    {
      Real_wp kz = pi_wp*m/L;
      Real_wp coef_mu = kz*kz;
      
      // calcul de la matrice
      Ah = mat_stiff;
      for (int i = 0; i < N; i++)
        Ah.Get(i, i) += mass_rho(i) + mass_mu(i)*coef_mu;

      // calcul du second membre
      VectComplex_wp rhs(N);
      rhs.Zero();
      SetPoints<Dimension2> Pts;
      Real_wp dz = 2.0*L / nb_modes;
      Real_wp r2 = square(var_mode.GetRadiusSource());
      for (int i = 0; i < nb_elt; i++)
        {
          const ElementReference<Dimension2, 1>& Fb = Vars.GetReferenceElementH1(i);
          bool affine = Vars.mesh.IsElementAffine(i);
          int nb_pts_quad = Fb.GetNbPointsQuadratureInside();
          Vars.mesh.GetVerticesElement(i, s);
          Fb.FjElemQuadrature(s, Pts, Vars.mesh, i);
          IVect Nodle = Vars.GetDofNumberOnElement(i);
          for (int j = 0; j < nb_pts_quad; j++)
            {
              int num_dof = Nodle(j);
              Real_wp jacob = Vars.GetWeightedJacobian(i, j, affine, Fb.GetGeometricElement());
              R2 pt = Pts.GetPointQuadrature(j);
              Real_wp x = pt(0), y = pt(1);
              Real_wp x0 = var_mode.GetCenterSource()(0);
              Real_wp y0 = var_mode.GetCenterSource()(1);
              Real_wp r_p = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));
              if (r_p > 2.0*var_mode.GetRadiusSource())
                continue;
              
              Complex_wp fm(0, 0);
              Real_wp coef = log(1e-6);
              for (int k = 0; k < nb_modes; k++)
                {
                  Real_wp z = k*dz;
                  if (k >= nb_modes/2)
                    z -= 2.0*L;
                  
                  Real_wp f = exp(coef*((x-x0)*(x-x0) + (y-y0)*(y-y0) + z*z)/r2);
                  fm += f*exp(-Iwp*kz*z);
                }
              
              rhs(num_dof) += fm*jacob;
            }
        }

      DISP(Norm2(rhs));

      // calcul de la solution
      All_MatrixLU<Complex_wp> mat_lu;
      mat_lu.Factorize(Ah);

      mat_lu.Solve(rhs);

      // on somme les solutions (z = 0)
      x = x + rhs;
    }

  Vars.WriteDatas(x);
}

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);

      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);

      EllipticProblem<HelmholtzEquation<Dimension2> > Vars;
      RunAll(Vars, file_name_data, type_element, type_equation);
    }

  return FinalizeMontjoie();
}
