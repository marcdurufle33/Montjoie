#define MONTJOIE_WITH_ONE_DIM
#define MONTJOIE_WITH_NODAL_H1

#include "Harmonic/MontjoieHarmonic.hxx"

#include "Elliptic/Helmholtz/Helmholtz1D.hxx"
#include "Elliptic/Helmholtz/Helmholtz1D_Inline.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "Elliptic/Helmholtz/Helmholtz1D.cxx"
#endif

using namespace Montjoie;

class InputVariables : public InputDataProblem_Base
{
public :
  Real_wp Lr_max, Li_max, Li_min, spectrum_radius_min;
  Real_wp pt_relevement;
  int type_relevement;
  enum { STANDARD, BORIS1, BORIS2};

  int type_source;
  enum {GAUSSIAN, SINUS, INCIDENT};

  int type_projection;
  enum {PROJ_INTEGRAL, PROJ_BIORTHO};
  
  Real_wp fmin, fmax; int nb_freq;
  bool dispersive_pml; Real_wp freq_centrale_pml;

  string file_fem_sol, file_qnm_sol;
  
  InputVariables()
  {
    Lr_max = Real_wp(0.6); Li_max = Real_wp(30); Li_min = Real_wp(0);
    pt_relevement = Real_wp(1e300);
    spectrum_radius_min = Real_wp(0); type_relevement = STANDARD;
    fmin = 0.01; fmax = 2.0; nb_freq = 100;
    type_source = INCIDENT;
    type_projection = PROJ_BIORTHO;
    dispersive_pml = true; freq_centrale_pml = Real_wp(0);
    file_fem_sol = "SolFem"; file_qnm_sol = "SolQnm";
  }
  
  bool SelectEigenvalue(const Complex_wp& z)
  {
    if ((realpart(z) <= Lr_max)
	&& (abs(imagpart(z)) <= Li_max)
	&& (abs(imagpart(z)) >= Li_min)
	&& (abs(z) >= spectrum_radius_min))
      return true;

    return false;
  }
  
  void SetInputData(const string& keyword, const Vector<string>& param)
  {
    if (keyword == "SpectrumBox")
      {
	Lr_max = to_num<Real_wp>(param(0));
	Li_max = to_num<Real_wp>(param(1));
        if (param.GetM() > 2)
          Li_min = to_num<Real_wp>(param(2));
      }
    else if (keyword == "SpectrumRadius")
      spectrum_radius_min = to_num<Real_wp>(param(0));
    else if (keyword == "Relevement")
      {
	pt_relevement = to_num<Real_wp>(param(0));
	if (param.GetM() > 1)
	  {
	    if (param(1) == "BORIS1")
	      type_relevement = BORIS1;
	    else if (param(1) == "BORIS2")
	      type_relevement = BORIS2;
	    else
	      type_relevement = STANDARD;
	  }
      }
    else if (keyword == "RangeFrequency")
      {
	fmin = to_num<Real_wp>(param(0));
	fmax = to_num<Real_wp>(param(1));
	nb_freq = to_num<int>(param(2));
      }
    else if (keyword == "SourceF")
      {
	if (param(0) == "Gaussian")
	  type_source = GAUSSIAN;
	else if (param(0) == "Sinus")
	  type_source = SINUS;
	else if (param(0) == "Incident")
	  type_source = INCIDENT;
      }
    else if (keyword == "ProjectionQNM")
      {
	if (param(0) == "Integral")
	  type_projection = PROJ_INTEGRAL;
	else if (param(0) == "Biortho")
	  type_projection = PROJ_BIORTHO;
	else
	  {
	    cout << "Unknown projection " << param(0) << endl;
	    abort();
	  }
      }
    else if (keyword == "DispersivePML")
      {
	if (param(0) == "YES")
	  dispersive_pml = true;
	else if (param(0) == "NO")
	  {
	    dispersive_pml = false;
	    freq_centrale_pml = to_num<Real_wp>(param(1));
	  }
	else
	  {
	    cout << "Invalid PML" << endl;
	    abort();
	  }
      }
    else if (keyword == "FileSolutionQNM")
      file_qnm_sol = param(0);
    else if (keyword == "FileSolutionFEM")
      file_fem_sol = param(0);
  }

};

class QnmSolver
{
protected:
  EllipticProblem<LaplaceEquation1D> vars;
  InputVariables input_var;

  Matrix<Real_wp> Kh, Mh;
  Matrix<Complex_wp> Mh_cplx;
  VectReal_wp F;

  VectComplex_wp AllLambda;
  Vector<VectComplex_wp> AllEigenVec;
  Matrix<Complex_wp> Vc;
  VectComplex_wp norme_mode_qnm;
  IVect num_extract;
  
  GridInterpolationFull<Dimension1> grid;
  int double_prec ; bool ascii; int Ngrid;
  int ref_cavity;
  
public:
  QnmSolver() { ref_cavity = 2; }

  void ComputeSource()
  {
    int N = vars.GetNbDof();
    VectReal_wp Mh_vec;
    F.Reallocate(N);
    F.Zero();
    const VectReal_wp& coor_dofs = vars.GetCoordinateDof();
    vars.GetMassMatrix(Mh_vec);
    //   on remplit la source comme une gaussienne centree en 0.6
    if (input_var.type_source != input_var.INCIDENT)
      for (int i = 0; i < vars.mesh.GetNbDof(); i++)
	if (abs(coor_dofs(i)) <= 1.0)
	  {
	    if (input_var.type_source == input_var.SINUS)
	      F(i) = sin(pi_wp*coor_dofs(i))*Mh_vec(i);
	    else if (input_var.type_source == input_var.GAUSSIAN)
	      F(i) = exp(-100.0*square(coor_dofs(i) - 0.6))*Mh_vec(i);
	  }
    
    //pour Dirichlet, on met a 0
    if (vars.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      F(0) = 0;    
  }
  
  void ConstructMatrix(const string& input_file, const string& name_element)
  {
    ReadInputFile(input_file, vars);
    ReadInputFile(input_file, input_var);
    
    cout << "The data file has been read" << endl;
    
    // initialisation calcul
    // true : formulation d'ordre 1 (p, u)
    // false : formulation d'ordre 2 (juste u)
    vars.SetFirstOrderFormulation(true);
    
    // ici je rajoute mes points de sortie
    // GridInterpolation<Dimension1>& var_section = vars.GetSectionGrid();
    // VectReal_wp points_sortie;
    // points_sortie.ReadText("points.dat");
    // var_section.Init(points_sortie);
    
    // ici on calcule grille et maillage
    vars.ComputeMeshAndFiniteElement(name_element);    
    vars.ComputePhysicalCoefficients();
    vars.ComputeMassMatrix();
    
    // (-i omega Mh + Kh) X = F
    // calcul des matrices elements finis
    GlobalGenericMatrix<Real_wp> nat_mat;
    nat_mat.SetCoefMass(0);
    nat_mat.SetCoefStiffness(1);
    nat_mat.SetCoefDamping(1);

    int N = vars.GetNbDof();
    int Nscal = vars.mesh.GetNbDof();
    cout << "Order = " << vars.mesh.GetOrder() << endl;
    cout << "Number of dofs = " << N << endl;
    cout << "Number of dofs for E = " << Nscal << endl;

    ComputeSource();

    // on calcule Kh
    if (!input_var.dispersive_pml)
      nat_mat.SetCoefDamping(0);
    
    Kh.Reallocate(N, N); Kh.Zero();
    vars.AddMatrixFEM(Kh, nat_mat);
    vars.AddBoundaryTerms(Kh, nat_mat);
    if (vars.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      Kh.ClearRow(0);

    if (vars.GetRightBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
      Kh.ClearRow(Nscal-1);

    // on calcule Mh
    nat_mat.SetCoefMass(1);
    nat_mat.SetCoefStiffness(0);
    nat_mat.SetCoefDamping(0);

    if (!input_var.dispersive_pml)
      {
	GlobalGenericMatrix<Complex_wp> nat_mat_c;	
	nat_mat_c.SetCoefMass(Complex_wp(1, 0));
	Real_wp omega0 = 2.0*pi_wp*input_var.freq_centrale_pml;
	nat_mat_c.SetCoefDamping(Iwp/omega0);
	nat_mat_c.SetCoefStiffness(Complex_wp(0, 0));
	Mh_cplx.Reallocate(N, N); Mh_cplx.Zero();
	vars.AddMatrixFEM(Mh_cplx, nat_mat_c);

	Mh_cplx.Write("Mh_cplx.dat");
      }
    else
      {
	Mh.Reallocate(N, N); Mh.Zero();
	vars.AddMatrixFEM(Mh, nat_mat);
	if (vars.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
	  {
	    Mh.ClearRow(0);
	    Mh(0, 0) = 1.0;
	  }
	
	if (vars.GetRightBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
	  {
	    Mh.ClearRow(Nscal-1);
	    Mh(Nscal-1, Nscal-1) = 1.0;
	  }

	Mh.WriteText("Mh.dat");
      }
    
    Kh.WriteText("Kh.dat");
  }

  void ComputeEigenmodes()
  {
    int N = vars.GetNbDof();
    int nev = N;
    
    if (input_var.dispersive_pml)
      {
	Matrix<Real_wp> Ah;
	
	// on calcule Ah = Mh^{-1} Kh
	Ah.Reallocate(N, N);
	
	for (int i = 0; i < N; i++)
	  for (int j = 0; j < N; j++)
	    Ah(i, j) = Kh(i, j) / Mh(i, i);
	
	Ah.WriteText("Ah.dat");    
	
	// calcul des valeurs propres et vecteurs propres
	VectReal_wp Lr, Li; Matrix<Real_wp> V;
	GetEigenvaluesEigenvectors(Ah, Lr, Li, V);
	
	Lr.WriteText("Lr.dat");   Li.WriteText("Li.dat");
	V.WriteText("V.dat");
	
	// boucle sur les valeurs propres
	// AllLambda : toutes les valeurs propres
	// AllEigenVec : vecteurs propres en complexe	
	AllLambda.Reallocate(nev);
	AllEigenVec.Reallocate(nev);  
	int nb = 0, num = 0;
	while (nb < nev)
	  {
	    if (Li(nb) == Real_wp(0))
	      {
		AllLambda(num) = Complex_wp(Lr(nb), Li(nb));
		AllEigenVec(num).Reallocate(N);
		for (int j = 0; j < N; j++)
		  AllEigenVec(num)(j) = Complex_wp(V(j, nb), 0);
		
		num++; nb++;
	      }
	    else
	      {
		AllLambda(num) = Complex_wp(Lr(nb), Li(nb));
		AllLambda(num+1) = Complex_wp(Lr(nb), -Li(nb));
		AllEigenVec(num).Reallocate(N);
		AllEigenVec(num+1).Reallocate(N);
		for (int j = 0; j < N; j++)
		  {
		    AllEigenVec(num)(j) = Complex_wp(V(j, nb), V(j, nb+1));
		    AllEigenVec(num+1)(j) = Complex_wp(V(j, nb), -V(j, nb+1));
		  }
		
		/* 
		   if (Li(nb) < 0)
		   {
		   AllLambda(num) = Complex_wp(Lr(nb), -Li(nb));
		   Conjugate(AllEigenVec(num));
		   } */
		
		
		// num++;
		
		num += 2;
		nb += 2;
	      }
	  }
	
	AllLambda.Resize(num); AllEigenVec.Resize(num);
      }
    else
      {
	Matrix<Complex_wp> Ah, V;
	Ah.Reallocate(N, N); V.Reallocate(N, N);
	for (int i = 0; i < N; i++)
	  for (int j = 0; j < N; j++)
	    Ah(i, j) = Kh(i, j) / Mh_cplx(i, i);

	AllLambda.Reallocate(nev);
	AllEigenVec.Reallocate(nev);
	GetEigenvaluesEigenvectors(Ah, AllLambda, V);

	VectReal_wp Lr(nev), Li(nev);
	for (int i = 0; i < nev; i++)
	  {
	    Lr(i) = realpart(AllLambda(i));
	    Li(i) = imagpart(AllLambda(i));
	  }
	
	Lr.WriteText("Lr.dat");   Li.WriteText("Li.dat");
		
	for (int i = 0; i < nev; i++)
	  {
	    AllEigenVec(i).Reallocate(N);
	    for (int j = 0; j < N; j++)
	      AllEigenVec(i)(j) = V(j, i);
	  }
      }
    
    // P : vecteurs propres de A
    // Vc : vecteurs propres biorthogonaux P^{-1}*
    // on peut les calculer directement
    Vc.Reallocate(nev, nev);
    /* if (input_var.dispersive_pml)
      {
	int Nscal = vars.mesh.GetNbDof();
	Vector<Complex_wp> vec_u(N);
	for (int i = 0; i < nev; i++)
	  {
	    for (int j = 0; j < Nscal; j++)
	      vec_u(j) = Mh(j, j)*AllEigenVec(i)(j);

	    for (int j = Nscal; j < N; j++)
	      vec_u(j) = -Mh(j, j)*AllEigenVec(i)(j);

	    Complex_wp coef = Real_wp(1) / DotProd(vec_u, AllEigenVec(i));
	    vec_u *= coef;
	    Conjugate(vec_u);
	    SetCol(vec_u, i, Vc);
	  }
      }
    */
    // ancien calcul des biorthogonaux
    Vc.Reallocate(nev, nev);
    for (int i = 0; i < nev; i++)
      for (int j = 0; j < nev; j++)
	Vc(i, j) = AllEigenVec(j)(i);
    
    GetInverse(Vc);
    Transpose(Vc); Conjugate(Vc);
    
    // num_extract : numero des valeurs propres qui seront conservees
    //num_extract.ReadText("number_extract.dat");
    
    // AllLambda : valeurs propres qui sont dans la boite
    // AllEigenVec : vecteurs propres associes
    
    int num = 0;
    for (int i = 0; i < AllLambda.GetM(); i++)
      if (input_var.SelectEigenvalue(AllLambda(i)))
	num++;
    
    num_extract.Reallocate(num);
    num = 0;
    for (int i = 0; i < AllLambda.GetM(); i++)
      if (input_var.SelectEigenvalue(AllLambda(i)))
	num_extract(num++) = i;
    
    DISP(N);
    DISP(num_extract);  
    for (int i = 0; i < num_extract.GetM(); i++)
      {
	AllLambda(i) = AllLambda(num_extract(i));
	AllEigenVec(i) = AllEigenVec(num_extract(i));      
      }
    
    num = num_extract.GetM();
    AllLambda.Resize(num); AllEigenVec.Resize(num);
    nev = num;

    VectReal_wp Lr, Li;
    Lr.Reallocate(nev);
    Li.Reallocate(nev);
    for (int i = 0; i < nev; i++)
      {
	Lr(i) = realpart(AllLambda(i));
	Li(i) = imagpart(AllLambda(i));
      }
    
    Lr.WriteText("eigenval_real.dat");
    Li.WriteText("eigenval_imag.dat");
  }

  void InitInterpolationGrid()
  {
    double_prec = OutputTypeEnum::DOUBLE_PRECISION;;
    ascii = false;
    GridInterpolation<Dimension1>& var_section = vars.GetSectionGrid();
    int type = GridInterpolationFull<Dimension1>::LINE;
    
    Ngrid = var_section.GetNbPointsGrid();
    grid.Init(type, var_section.GetXmin(), var_section.GetXmax(), Ngrid);
  }

  void DisplayEigenmodes()
  {
    int nev = AllLambda.GetM();
    GridInterpolation<Dimension1>& var_section = vars.GetSectionGrid();
    int N = vars.GetNbDof();
    int Nscal = vars.mesh.GetNbDof();
    
    // affichage des modes
    VectComplex_wp output_vector(Ngrid);
    for (int i = 0; i < nev; i++)
      {
	// on sort le mode
	vars.ComputeInterpolationU(AllEigenVec(i), var_section, output_vector);
	string name_file = "ModeE" + to_str(i) + ".dat";
	WriteMatlab(output_vector, grid, name_file, double_prec, ascii);
	
	int Nv = N - Nscal;
	Vector<Complex_wp> V(Nv);
	for (int j = 0; j < Nv; j++)
	  V(j) = AllEigenVec(i)(Nscal + j);
	
	vars.ComputeInterpolationU(V, var_section, output_vector, true, false);
	name_file = "ModeH"+ to_str(i) + ".dat";
	WriteMatlab(output_vector, grid, name_file, double_prec, ascii);      
      }

  }
  

  void CheckOrthogonality()
  {
    int Nscal = vars.mesh.GetNbDof();
    int N = vars.GetNbDof();
    int nev = AllLambda.GetM();
    norme_mode_qnm.Reallocate(nev);
    for (int i = 0; i < nev; i++)
      {
	// verification de l'orthogonalite
	if (false)
	for (int i2 = 0; i2 < nev; i2++)
	  {
	    Complex_wp prod_scal(0, 0);
	    for (int j = 0; j < Nscal; j++)
	      prod_scal += Mh(j, j)*AllEigenVec(i)(j)*AllEigenVec(i2)(j);

	    if (vars.mesh.GetNbLayersPML() > 0)
	      {
		for (int j = Nscal; j < N; j++)
		  prod_scal -= Mh(j, j)*AllEigenVec(i)(j)*AllEigenVec(i2)(j);
	      }
	    else
	      {
		// terme de bord (pour abc)
		prod_scal -= Real_wp(1)/(AllLambda(i) + AllLambda(i2)) *
		  (AllEigenVec(i)(Nscal-1) * AllEigenVec(i2)(Nscal-1)
		   +AllEigenVec(i)(0) * AllEigenVec(i2)(0));
	      }
	    
	    if ((i != i2) && (abs(prod_scal) > 1e-10))
	    {DISP(i); DISP(i2); DISP(prod_scal); }
	    
	    if (i2 == i) { DISP(i); DISP(prod_scal); }
	  }
	
	// norme du mode 
	Complex_wp normMode(0, 0);
	if (input_var.dispersive_pml)
	  for (int j = 0; j < Nscal; j++)
	    normMode += Mh(j, j)*square(AllEigenVec(i)(j));
	else
	  for (int j = 0; j < Nscal; j++)
	    normMode += Mh_cplx(j, j)*square(AllEigenVec(i)(j));

	if (true)
	  {
	    if (input_var.dispersive_pml)
	      for (int j = Nscal; j < N; j++)
		normMode -= Mh(j, j)*AllEigenVec(i)(j)*AllEigenVec(i)(j);
	    else
	      for (int j = Nscal; j < N; j++)
		normMode -= Mh_cplx(j, j)*AllEigenVec(i)(j)*AllEigenVec(i)(j);
	  }
	else
	  {
	    // avec terme de bord
	    normMode -= Real_wp(1) / (2.0*AllLambda(i)) *
	      (square(AllEigenVec(i)(Nscal-1)) + square(AllEigenVec(i)(0)));

	    normMode *= 2.0;
	  }
	
	norme_mode_qnm(i) = normMode;
	
	//DISP(normMode);
	//for (int j = Nscal; j < N; j++)
	//normMode -= Mh(j, j)*square(AllEigenVec(i)(j));

	//DISP(i); DISP(normMode);
	//cout << "waiting" << endl; int test_input; cin >> test_input;      
      }
  }
  
  void ComputeSourceOmega(const Real_wp& omega, VectComplex_wp& Fcplx)
  {
    int N = vars.GetNbDof();
    VectReal_wp s(2);
    if (input_var.type_source != input_var.INCIDENT)
      for (int k = 0; k < N; k++)
	Fcplx(k) = Complex_wp(F(k), 0);
    else
      {
	Fcplx.Zero();
	for (int j = 0; j < vars.mesh.GetNbElt(); j++)
	  if (vars.mesh.Element(j).GetReference() == ref_cavity)
	    {
	      vars.mesh.GetVerticesElement(j, s);
	      Real_wp h = abs(s(1)-s(0)); //DISP(j); DISP(s);
	      Real_wp eps = vars.ref_rho(ref_cavity).GetConstant();
	      const ElementReference<Dimension1, 1>& Fb = vars.GetReferenceElementH1(j);
	      for (int k = 0; k < Fb.GetNbDof(); k++)
		{
		  Real_wp x = Fb.Points(k)*s(1) + (1.0-Fb.Points(k))*s(0);
		  int num_dof = vars.mesh.GetNumberDof(j, k);
		  Real_wp poids = h*Fb.Weights(k);
		  Fcplx(num_dof) += Iwp*poids*omega*(eps-1.0)*exp(Iwp*omega*x);
		}
	    }
	
	if (vars.GetLeftBoundaryCondition() == BoundaryConditionEnum::LINE_DIRICHLET)
	  Fcplx(0) = 0.0;
      }
  }

  void ComputeProjectionQNM(const VectComplex_wp& invM_Fcplx,
			    const Real_wp& omega, VectComplex_wp& decomp)
  {
    int nev = AllLambda.GetM();
    int N = vars.GetNbDof();
    VectComplex_wp rhs(nev), biortho_qnm(N);
    if (input_var.type_projection == input_var.PROJ_BIORTHO)
      {
	// produit scalaire avec le vecteur biorthogonal
	for (int i = 0; i < nev; i++)
	  {
	    GetCol(Vc, num_extract(i), biortho_qnm);
	    rhs(i) = DotProdConj(biortho_qnm, invM_Fcplx);
	  }
	
	// decomposition en utilisant les vecteurs biorthogonaux
	decomp = rhs;
	for (int i = 0; i < nev; i++)
	  decomp(i) /= -Iwp*omega + AllLambda(i);
      }
    else
      {
	// rhs = 1/2 \int_{-a}^a f(x) u_j(x) dx
	rhs.Zero();
	for (int i = 0; i < nev; i++)
	  {
	    // expression analytique
	    /*VectReal_wp s(2);
	    for (int j = 0; j < vars.mesh.GetNbElt(); j++)
	      if (vars.mesh.Element(j).GetReference() == ref_cavity)
		{
		  vars.mesh.GetVerticesElement(j, s);
		  Real_wp h = abs(s(1)-s(0)); //DISP(j); DISP(s);
		  const ElementReference<Dimension1, 1>& Fb = vars.GetReferenceElement(j);
		  for (int k = 0; k < Fb.GetNbDof(); k++)
		    {
		      int num_dof = vars.mesh.GetNumberDof(j, k);
		      Real_wp poids = h*Fb.Weights(k);
		      rhs(i) += poids*invM_Fcplx(num_dof)*AllEigenVec(i)(num_dof);
		    }
		}

	    // on multiplie ici par rho parce que invM_Fcplx nous donne 1/rho f
	    rhs(i) *= 0.5*vars.ref_rho(ref_cavity).GetConstant(); */

	    // expression plus numerique
	    int Nscal = vars.mesh.GetNbDof();
	    if (input_var.dispersive_pml)
	      for (int j = 0; j < Nscal; j++)
		rhs(i) += invM_Fcplx(j)*Mh(j, j)*AllEigenVec(i)(j);
	    else
	      for (int j = 0; j < Nscal; j++)
		rhs(i) += invM_Fcplx(j)*Mh_cplx(j, j)*AllEigenVec(i)(j);
	  }

	// on applique la formule analytique avec norme_mode_qnm  = << u_j, u_j >>
	decomp = rhs;
	for (int i = 0; i < nev; i++)
	  decomp(i) /= norme_mode_qnm(i)*(AllLambda(i) - Iwp*omega);
      }    
  }

  void ComputeLinearCombination_Relevement(const VectComplex_wp& decomp,
					   const Real_wp& omega, VectComplex_wp& X)
  {
    int nev = AllLambda.GetM();
    int N = vars.GetNbDof();
    const VectReal_wp& coor_dofs = vars.GetCoordinateDof();
    VectComplex_wp biortho_qnm(N);
    
    X.Zero();
    if (input_var.type_relevement == input_var.BORIS1)
      {	
	for (int i = 0; i < nev; i++)
	  {
	    // 1er relevement de Boris
	    if (abs(input_var.pt_relevement) < 1e100)
	      {
		Complex_wp omega_q = -Iwp*AllLambda(i);
		Real_wp xmax = input_var.pt_relevement;
		Complex_wp beta = decomp(i)*exp(-Iwp*(omega-omega_q)*xmax);
		Complex_wp betam = decomp(i)*exp(Iwp*(omega-omega_q)*xmax);
		for (int j = 0; j < N; j++)
		  if (coor_dofs(j) > input_var.pt_relevement)
		    {
		      Complex_wp phase = exp(Iwp*(omega-omega_q)*coor_dofs(j));
		      X(j) += beta*AllEigenVec(i)(j)*phase;
		    }
		  else	if (coor_dofs(j) < -input_var.pt_relevement)
		    {
		      Complex_wp phase = exp(-Iwp*(omega-omega_q)*coor_dofs(j));
		      X(j) += beta*AllEigenVec(i)(j)*phase;
		    }
		  else
		    X(j) += decomp(i)*AllEigenVec(i)(j);
	      }
	    else
	      for (int i = 0; i < nev; i++)
		X += decomp(i)*AllEigenVec(i);
	  }
      }
    else if (input_var.type_relevement == input_var.BORIS2)
      {
	// 2ND RELEVEMENT DE BORIS
	if (abs(input_var.pt_relevement) < 1e100)
	  {
	    Matrix<Complex_wp> Ab(nev,nev);
	    VectComplex_wp uq(N);
	    for (int i = 0; i < nev; i++)
	      for (int j = 0; j < nev; j++)
		{
		  for (int k = 0; k < N; k++)
		    {
		      Complex_wp omega_q = -Iwp*AllLambda(j);
		      Complex_wp phase = exp(Iwp*(omega-omega_q)*coor_dofs(k));
		      uq(k)= AllEigenVec(j)(k)*phase;
		    }
		  
		  GetCol(Vc, num_extract(i), biortho_qnm);
		  Ab(i,j) = DotProdConj(biortho_qnm, uq);
		}

	    Vector<int> pivot;
	    GetLU(Ab, pivot);
	    VectComplex_wp beta(decomp);
	    SolveLU(Ab, pivot, beta);
	    
	    for (int i = 0; i < nev; i++)
	      {
		Complex_wp omega_q = -Iwp*AllLambda(i);
		for (int j = 0; j < vars.mesh.GetNbDof(); j++)
		  if (coor_dofs(j) > input_var.pt_relevement)
		    {
		      Complex_wp phase = exp(Iwp*(omega-omega_q)*coor_dofs(j));
		      X(j) += beta(i)*AllEigenVec(i)(j)*phase;
		    }
		  else
		    X(j) += decomp(i)*AllEigenVec(i)(j);
	      }
	  }
      }
    else
      {
	for (int i = 0; i < nev; i++)
	  X += decomp(i)*AllEigenVec(i);
	
	if (abs(input_var.pt_relevement) < 1e100)
	  {
	    int dof0 = -1, dof1 = -1;
	    for (int i = 0; i < vars.mesh.GetNbDof(); i++)
	      {
		if (abs(coor_dofs(i) - input_var.pt_relevement) <= 1e-10)
		  dof0 = i;
		
		if (abs(coor_dofs(i) + input_var.pt_relevement) <= 1e-10)
		  dof1 = i;
	      }
	    
	    if (dof0 == -1)
	      {
		cout << "Point " << input_var.pt_relevement << " non trouve " << endl;
		abort();
	      }
	    else
	      {
		Complex_wp val0 = X(dof0); //DISP(val0);
		Complex_wp val1(0, 0);
		if (dof1 != -1)
		  val1 = X(dof1);
		
		for (int i = 0; i < vars.mesh.GetNbDof(); i++)
		  {
		    if (coor_dofs(i) > input_var.pt_relevement)
		      X(i) = val0 * exp(Iwp*omega*(coor_dofs(i) - input_var.pt_relevement));
		    else if (coor_dofs(i) < -input_var.pt_relevement)
		      X(i) = val1 * exp(-Iwp*omega*(coor_dofs(i) + input_var.pt_relevement));
		  }
	      }
	  }	
      }
  }
  
  void RunAll(const string& input_file, const string& name_element)
  {
    cout << "On calcule les matrices elements finis" << endl;
    // calcul des matrices Ah, Mh et Kh
    this->ConstructMatrix(input_file, name_element);
    
    // calcul des valeurs propres qu'on conserve (modes QNM)
    cout << "On calcule les modes" << endl;
    this->ComputeEigenmodes();
    
    // initialisation de la grille d'affichage
    cout << "Initialisation de la grille" << endl;
    this->InitInterpolationGrid();

    // on affiche les modes
    cout << "On affiche les modes" << endl;
    this->DisplayEigenmodes();

    // on verifie la relation d'orthogonalite
    cout << "On verifie l'orthogonalite" << endl;
    this->CheckOrthogonality();

    int N = vars.GetNbDof();
    int nev = AllLambda.GetM();
    VectComplex_wp Fcplx(N), X(N), Y(N);
    IVect pivot;
    
    // boucle sur les frequences
    Real_wp err = 0, errI = 0;
    VectComplex_wp output_vector2(Ngrid), output_vector(Ngrid), decomp;
    GridInterpolation<Dimension1>& var_section = vars.GetSectionGrid();
    DISP(nev);
    DISP(input_var.type_relevement);
    cout << "Boucle en frequence" << endl;
    for (int m = 0; m < input_var.nb_freq; m++)
      {
	Real_wp f = input_var.fmin + m*(input_var.fmax-input_var.fmin)/(input_var.nb_freq-1);
	Real_wp omega = 2*pi_wp*f;
	
	// on remplit la source Fcplx
	this->ComputeSourceOmega(omega, Fcplx);
	
	// on inverse par la masse
	VectComplex_wp invM_Fcplx(N);
	if (input_var.dispersive_pml)
	  for (int i = 0; i < N; i++)
	    invM_Fcplx(i) = Fcplx(i) / Mh(i, i);
	else
	  for (int i = 0; i < N; i++)
	    invM_Fcplx(i) = Fcplx(i) / Mh_cplx(i, i);

	//DISP(m); DISP(omega); DISP(invM_Fcplx);
	// on projette sur les QNMs (calcul des alpha_i dans decomp)
	ComputeProjectionQNM(invM_Fcplx, omega, decomp);
	
	decomp.Write("alpha"+ to_str(m) + ".dat");	
	
	//DISP(decomp);
	// on calcule la combinaison lineaire \sum \alpha_i qnm_i
	ComputeLinearCombination_Relevement(decomp, omega, X);

	// on normalise decomp pour en faire une sortie (pour afficher le truc en couleur)
	if (norme_mode_qnm.GetM() >= decomp.GetM())
	  {
	    for (int i = 0; i < decomp.GetM(); i++)
	      decomp(i) *= sqrt(norme_mode_qnm(i));
	    
	    decomp.Write("beta"+ to_str(m) + ".dat");
	  }
	
	//DISP(X);
	// possibilite plus couteuse : calcul direct de (-i omega M_h + K_h)^{-1} F
	// on le fait ici en plein, on pourrait le faire en creux
	Matrix<Complex_wp> A(N, N);
	if (input_var.dispersive_pml)
	  {
	    for (int i = 0; i < N; i++)
	      for (int j = 0; j < N; j++)
		A(i, j) = -Iwp*omega*Mh(i, j) + Kh(i, j);
	  }
	else
	  {
	    for (int i = 0; i < N; i++)
	      for (int j = 0; j < N; j++)
		A(i, j) = -Iwp*omega*Mh_cplx(i, j) + Kh(i, j);
	  }

	GetLU(A, pivot);
	Y = Fcplx;
	SolveLU(A, pivot, Y);
	
	// calcul de l'erreur entre X et Y
	err = max(err, Norm2(X-Y)/Norm2(X));
      
	// on ecrit la solution QNM
	vars.ComputeInterpolationU(X, var_section, output_vector);
	
	string name_file = input_var.file_qnm_sol + to_str(m) + ".dat";
	WriteMatlab(output_vector, grid, name_file, double_prec, ascii);
	
	// on ecrit la solution couteuse
	vars.ComputeInterpolationU(Y, var_section, output_vector2);
	
	name_file = input_var.file_fem_sol + to_str(m) + ".dat";
	WriteMatlab(output_vector2, grid, name_file, double_prec, ascii);
	
	errI = max(errI, Norm2(output_vector - output_vector2) / Norm2(output_vector));
      
	//cout << "waiting" << endl; int test_input; cin >> test_input;
      }

    // on affiche les erreurs
    DISP(err);
    DISP(errI);
  }

};


int main(int argc, char **argv) 
{

  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
      // we add the default path to the file name given
      string file_name_data;
      file_name_data = string(argv[1]);
            
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);

      cout<<" Equation "<<type_equation<<" Type Element "<<type_element<<endl;
      
      QnmSolver solver_qnm;
      solver_qnm.RunAll(file_name_data, type_element);
      
      cout<<" we destroy the variables "<<endl;
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"mode1D.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }
  
  return FinalizeMontjoie();
  
}
