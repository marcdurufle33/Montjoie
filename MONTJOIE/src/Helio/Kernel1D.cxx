#ifndef MONTJOIE_FILE_KERNEL1D_CXX

#include "Helio/Kernel1D.hxx"
#include "Helio/Kernel1DInline.cxx"
#include "Helio/sphericalHarmonics.cxx"

namespace Montjoie
{

  Kernel1D::Kernel1D()
  {
  }


  void Kernel1D::ReadInputFile(const string filename, const int indFreqMin)
  {

    Vector<string> fileG(2), filedG(2), fileC(2), filedC(2);
    string filerho, filec, filer, fileW;
    kernelDimension = "2D";
    Rsun = 696.e6;

    Vector<string> lines_data_file;
    ReadLinesFile(filename, lines_data_file);

    for (int i=0; i < lines_data_file.GetM(); i++)
      {
	// Read keyword and value for each line of the config file
	Vector<string> ans = RealLine(lines_data_file(i));
	string keyword = ans(0);
	string value = ans(1);

	// Initialize the parameters of the class
	if (keyword == "kernelDimension")
	  kernelDimension = to_num<int>(value);	  
	else if (keyword == "lMaxGreen")
	  lMaxG = to_num<int>(value);
	else if (keyword == "lMaxKernel")
	  lMaxK = to_num<int>(value);
	else if (keyword == "nR")
	  nr = to_num<int>(value);
	else if (keyword == "nbKernels")
	  nbKernels = to_num<int>(value);
	else if (keyword == "typeOfKernel")
	  typeOfKernel = value;
	else if (keyword == "th1")
	  {
	    Vector<string> th = split(value, ' ');
	    theta1.Reallocate(nbKernels);
	    if (th(0) == "LIST")
	      {
		for (int i=0; i < nbKernels; i++)
		  {
		    if (th.GetM() == 2)
		      theta1(i) = to_num<Real_wp>(th(1))*pi_wp/180.;
		    else
		      theta1(i) = to_num<Real_wp>(th(i+1))*pi_wp/180.;
		  }
	      }
	    else if (th(0) == "RANGE")
	      {
		Real_wp thMin = to_num<Real_wp>(th(1))*pi_wp/180.;
		Real_wp thSample = to_num<Real_wp>(th(2))*pi_wp/180.;
		for (int i=0; i < nbKernels; i++)
		  {
		    theta1(i) = thMin + 1.*i*thSample;
		  }
	      }
	  }
	else if (keyword == "th2")
	  {
	    Vector<string> th = split(value, ' ');
	    theta2.Reallocate(nbKernels);
	    if (th(0) == "LIST")
	      {
		for (int i=0; i < nbKernels; i++)
		  {
		    if (th.GetM() == 2)
		      theta2(i) = to_num<Real_wp>(th(1))*pi_wp/180.;
		    else
		      theta2(i) = to_num<Real_wp>(th(i+1))*pi_wp/180.;
		  }
	      }
	    else if (th(0) == "RANGE")
	      {
		Real_wp thMin = to_num<Real_wp>(th(1))*pi_wp/180.;
		Real_wp thSample = to_num<Real_wp>(th(2))*pi_wp/180.;
		for (int i=0; i < nbKernels; i++)
		  {
		    theta2(i) = thMin + 1.*i*thSample;
		  }
	      }
	  }
	else if (keyword == "phi1")
	  {
	    Vector<string> phi = split(value, ' ');
	    phi1.Reallocate(nbKernels);
	    for (int i=0; i < nbKernels; i++)
	      {
		if (phi.GetM() == 1)
		  phi1(i) = to_num<Real_wp>(phi(0))*pi_wp/180.;
		else
		  phi1(i) = to_num<Real_wp>(phi(i))*pi_wp/180.;
	      }
	  }
	else if (keyword == "phi2")
	  {
	    Vector<string> phi = split(value, ' ');
	    phi2.Reallocate(nbKernels);
	    for (int i=0; i < nbKernels; i++)
	      {
		if (phi.GetM() == 1)
		  phi2(i) = to_num<Real_wp>(phi(0))*pi_wp/180.;
		else
		  phi2(i) = to_num<Real_wp>(phi(i))*pi_wp/180.;
	      }
	  }
	else if (keyword == "Frequencies")
	  {
	    Vector<string> params = split(value, ' ');
	    Real_wp T = to_num<Real_wp>(params(1));
	    Real_wp dt = to_num<Real_wp>(params(2));
	    int Nthalf = (T/dt-1)/2;
	    Real_wp hf = 1. / T;
	    frequencies.Reallocate(Nthalf+1);
	    for (int i=0; i <= Nthalf; i++)
	      frequencies(i) = hf*i;
	  }
	else if (keyword == "typeOfObservable")
	  {
	    Vector<string> params = split(value, ' ');
	    typeOfObservable = params(0);
	    if (typeOfObservable == "XS")
	      indFreq = indFreqMin + to_num<int>(params(1));
	  }
	else if (keyword == "Filter")
	  {
	    Vector<string> params = split(value, ' ');
	    typeOfFilter = params(0);
	    centerFilter = to_num<Real_wp>(params(1));
	    widthFilter = to_num<Real_wp>(params(2));
	  }

	else if (keyword == "dirOut")
	  dirOut = value;
	else if (keyword == "fileG")
	  fileG = split(value, ' ');
	else if (keyword == "filedG")
	  filedG = split(value, ' ');
	else if (keyword == "fileC")
	  fileC = split(value, ' ');
	else if (keyword == "filedC")
	  filedC = split(value, ' ');
	else if (keyword == "filerho")
	  filerho = value;
	else if (keyword == "filec")
	  filec = value;
	else if (keyword == "filer")
	  filer = value;
	else if (keyword == "fileW")
	  fileW = value;

	else if (keyword == "Damping")
	  {
	    Vector<string> params = split(value, ' ');
	    string damping_type = params(0);
	    if (typeOfObservable == "XS")
	      gamma.Reallocate(1);
	    else
	      {
		gamma.Reallocate(frequencies.GetM());
		gamma.Zero();
	      }
	    if (damping_type == "POWER_LAW")
	      {
		Real_wp gamma0 = to_num<Real_wp>(params(1));
		Real_wp omega0 = to_num<Real_wp>(params(2));
		Real_wp beta = to_num<Real_wp>(params(3));
		if (typeOfObservable == "XS")
		  gamma(0) = gamma0 * pow(abs(2.*pi_wp*frequencies(indFreq)/omega0), beta);
		else
		  {
		    for (int i=0; i < frequencies.GetM(); i++)
		      gamma(i) = gamma0 * pow(abs(2.*pi_wp*frequencies(i)/omega0), beta);
		  }
	      }
	    else if (damping_type == "PROPORTIONAL")
	      {
		Real_wp gamma0 = to_num<Real_wp>(params(1));
		if (typeOfObservable == "XS")
		  gamma(0) =  gamma0 * abs(2.*pi_wp*frequencies(indFreq));
		else
		  {
		    for (int i=0; i < frequencies.GetM(); i++)
		      gamma(i) =  gamma0 * abs(2.*pi_wp*frequencies(i));
		  }
	      }
	    else if (damping_type == "CONSTANT")
	      {
		Real_wp gamma0 = to_num<Real_wp>(params(1));
		if (typeOfObservable == "XS")
		  gamma(0) = gamma0;
		else
		  {
		    for (int i=0; i < frequencies.GetM(); i++)
		      gamma(i) =  gamma0;
		  }
	      }
	    else
	      {
		cout << "This type of damping is not implemented yet" << endl;
		abort();
	      }

	  }
      }

    if (!filer.empty())
      r = ReadRealVect(filer, nr);
    if (!filerho.empty())
      rho = ReadRealVect(filerho, nr);
    if (!filec.empty())
      c = ReadRealVect(filec, nr);

    if (typeOfObservable == "XS")
      {
	omega = 2.*pi_wp*frequencies(indFreq);
	Vector<Real_wp> F(lMaxG+1); F.Zero();
	if (!typeOfFilter.empty())
	  {
	    for (int l=1; l < lMaxG+1; l++)
	      {
		if (typeOfFilter == "PhaseSpeed")
		  F(l) = exp(-pow(abs(omega)/(l/Rsun)-centerFilter,2)/(2*widthFilter*widthFilter));
		else if (typeOfFilter == "Lfilter")
		  {
		    F(l) = 0.5*(1.-tanh(0.1*(l-centerFilter)));
		    //if (l > centerFilter - widthFilter/2 and l < centerFilter + widthFilter/2)
		    //  F(l) = 1;
		    }
	      }
	  }
	if (!fileG(0).empty())
	  {
	    string toRead = fileG(0) + to_str(indFreq) + fileG(1);
	    G =  reshapeVector(ReadComplexVect(toRead,(lMaxG+1)*nr),nr, lMaxG+1);
	    if (!typeOfFilter.empty())
	      {
		for (int l=0; l < lMaxG+1; l++)
		  {
		    for (int i=0; i < nr; i++)
		      G(i,l) *= F(l);
		  }
	      }
	  }
	if (!filedG(0).empty())
	  {
	    string toRead = filedG(0) + to_str(indFreq) + filedG(1);
	    dG =  reshapeVector(ReadComplexVect(toRead,(lMaxG+1)*nr),nr, lMaxG+1);
	    if (!typeOfFilter.empty())
	      {
		for (int l=0; l < lMaxG+1; l++)
		  {
		    for (int i=0; i < nr; i++)
		      dG(i,l) *= F(l);
		  }
	      }
	  }
	if (!fileC(0).empty())
	  {
	    string toRead = fileC(0) + to_str(indFreq) + fileC(1);
	    C =  reshapeVector(ReadComplexVect(toRead,(lMaxG+1)*nr),nr, lMaxG+1);
	    if (!typeOfFilter.empty())
	      {
		for (int l=0; l < lMaxG+1; l++)
		  {
		    for (int i=0; i < nr; i++)
		      C(i,l) *= F(l);
		  }
	      }
	  }
	if (!filedC(0).empty())
	  {
	    string toRead = filedC(0) + to_str(indFreq) + filedC(1);
	    dC =  reshapeVector(ReadComplexVect(toRead,(lMaxG+1)*nr),nr, lMaxG+1);
	    if (!typeOfFilter.empty())
	      {
		for (int l=0; l < lMaxG+1; l++)
		  {
		    for (int i=0; i < nr; i++)
		      dC(i,l) *= F(l);
		  }
	      }
	  }
      }
    else
      {
	clock_t begin_time = clock();
	if (!fileW.empty())
	  W = ReadComplexVect(fileW, frequencies.GetM());

	fllpr = Vector<Matrix<Complex_wp> >(nr);
	gllpr = Vector<Matrix<Complex_wp> >(nr);
	Vector<Matrix<Complex_wp> > fllprCrt = Vector<Matrix<Complex_wp> >(nr);
	Vector<Matrix<Complex_wp> > gllprCrt = Vector<Matrix<Complex_wp> >(nr);
	for (int i=0; i < nr; i++)
	  {
	    fllpr(i) = Matrix<Complex_wp>(lMaxG+1,lMaxG+1); fllpr(i).Zero();
	    gllpr(i) = Matrix<Complex_wp>(lMaxG+1,lMaxG+1); gllpr(i).Zero();
	    fllprCrt(i) = Matrix<Complex_wp>(lMaxG+1,lMaxG+1); fllprCrt(i).Zero();
	    gllprCrt(i) = Matrix<Complex_wp>(lMaxG+1,lMaxG+1); gllprCrt(i).Zero();
	  }

	int nb_proc = 1, rank_proc = 0;
#ifdef SELDON_WITH_MPI    
	MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#endif
	
	// Separate depths  on each processor
	Vector<Vector<int> > listIndFreq(nb_proc);
	
	int aveJobs = frequencies.GetM()/nb_proc;
	int extraJobs = frequencies.GetM()-aveJobs*nb_proc;
	int cpt = 0;
	for (int i=0; i < nb_proc; i++)
	  {
	    int nbJobs = aveJobs;
	    if (i < extraJobs)
	      nbJobs++;
	    for (int j=0; j < nbJobs; j++)
	      {
		listIndFreq(i).PushBack(cpt);
		cpt++;
	      }
	  }
	
	Real_wp Ps = 1.;
	Real_wp homega = frequencies(1)-frequencies(0);
	Real_wp omega0 = 2.*pi_wp*3.e-3;
	//Real_wp sigma = 2.*pi_wp*6.e-4;
	Real_wp sigma = 2.*pi_wp*1.2e-3;
	Complex_wp I_wp = Complex_wp(0.,1.);
	for (int i=0; i < listIndFreq(rank_proc).GetM(); i++)
	  {
	    int ind = listIndFreq(rank_proc)(i);
	    if (ind > 0) // remove f=0
	      {
		if (rank_proc == 0)
		  cout << ind << " / " <<  listIndFreq(rank_proc).GetM() << endl;
		omega = 2.*pi_wp*frequencies(ind);
		string toRead = fileG(0) + to_str(ind) + fileG(1);
		G =  reshapeVector(ReadComplexVect(toRead,(lMaxG+1)*nr),nr, lMaxG+1);
		if (typeOfKernel == "UR")
		  {
		    toRead = filedC(0) + to_str(ind) + filedC(1);
		    dC =  reshapeVector(ReadComplexVect(toRead,(lMaxG+1)*nr),nr, lMaxG+1);
		  }
		else if (typeOfKernel == "UTHETA" or typeOfKernel == "UPHI" or typeOfKernel == "C")
		  {
		    toRead = fileC(0) + to_str(ind) + fileC(1);
		    C =  reshapeVector(ReadComplexVect(toRead,(lMaxG+1)*nr),nr, lMaxG+1);
		  }
		
		if (typeOfObservable == "tau")
		  {
		    //Real_wp Ps = exp(-((omega-omega0)*(omega-omega0) / (2*sigma*sigma)));
		    Ps = (sigma/2.)*(sigma/2.) / ((omega-omega0)*(omega-omega0)+(sigma/2.)*(sigma/2.) );
		  }
		for (int l=0; l<=lMaxG; l++)
		  {
		    for (int lp=0; lp<=lMaxG; lp++)
		      {
			for (int indr=0; indr < nr; indr++)
			  {
			    if (typeOfKernel == "UR")
			      {
				fllprCrt(indr)(l,lp) += -2.* imag(omega * conj(W(ind)) * G(indr,l) * dC(indr,lp)) * 2. * rho(indr) * homega * Ps;
				gllprCrt(indr)(l,lp) += 2.* imag(omega * conj(W(ind)) * conj(G(indr,l)) * conj(dC(indr,lp))) * 2. * rho(indr) * homega * Ps;
			      }
			    else if (typeOfKernel == "UTHETA" or typeOfKernel == "UPHI")
			      {
				fllprCrt(indr)(l,lp) += -2.* imag(omega * conj(W(ind)) * G(indr,l) * C(indr,lp)) * 2. * rho(indr) * homega * Ps;
				gllprCrt(indr)(l,lp) += 2.* imag(omega * conj(W(ind)) * conj(G(indr,l)) * conj(C(indr,lp))) * 2. * rho(indr) * homega * Ps;
			      }
			    else if (typeOfKernel == "C")
			      {
				fllprCrt(indr)(l,lp) += -2.* real(omega * (omega + 2.*I_wp * gamma(ind)) * conj(W(ind)) * G(indr,l) * C(indr,lp)) * 2. * rho(indr) / c(indr) * homega * Ps;
				gllprCrt(indr)(l,lp) += -2.* real(omega * (omega - 2.*I_wp * gamma(ind)) * conj(W(ind)) * conj(G(indr,l)) * conj(C(indr,lp))) * 2. * rho(indr) / c(indr) * homega * Ps;
			      }
			  }
		      }
		  }
	      }
	  }


	Vector<int64_t> xtmp;
	for (int indr=0; indr < nr; indr++)
	  {
	    MpiReduce(MPI_COMM_WORLD, fllprCrt(indr).GetData(), xtmp,
		      fllpr(indr).GetData(), fllprCrt(indr).GetDataSize(), MPI_SUM, 0);
	    MpiReduce(MPI_COMM_WORLD, gllprCrt(indr).GetData(), xtmp,
		      gllpr(indr).GetData(), gllprCrt(indr).GetDataSize(), MPI_SUM, 0);
	  }
	if (rank_proc == 0)
	  {
	    Real_wp seconds = 1.*(clock() - begin_time) / CLOCKS_PER_SEC;
	    cout << "Time to load the Green and cross-covariance functions: " << seconds << endl;
	  }

	//FftInterface<Complex_wp> fft;
	//y = x;
	//fft.ApplyForward(y);
	//x2 = y;
	//fft.ApplyInverse(x2);

      }
    
  }


  void Kernel1D::getGreen(string input_file)
{

  MontjoieTimer var_chrono;

  var_chrono.SetMessage("GlobalComputation", "perform the overall simulation");
  var_chrono.Start("GlobalComputation");
  
  var_chrono.SetMessage("ConstructAll", "initialize the computation");
  var_chrono.Start("ConstructAll");

  string name_elt;// = "EDGE_LOBATTO";
  string type_equation;// = "HELMHOLTZ";
  getElement_Equation(input_file, name_elt, type_equation);
  EllipticProblem<HelmholtzEquationRadial> var;
  
  // initialisation of computations
  var.ConstructAll(input_file, name_elt, var_chrono);
  
  var_chrono.Stop("ConstructAll");

  // computation of the right hand sides
  var_chrono.SetMessage("ComputationRhs", "compute the right hand side");
  var_chrono.Start("ComputationRhs");

  Vector<Vector<Vector<Complex_wp> > > rhs;
  var.ComputeRightHandSide(rhs, var_chrono);
  
  var_chrono.Stop("ComputationRhs");

  var_chrono.SetMessage("FactorisationMatrix", "compute and factorize finite element matrix");
  var_chrono.SetMessage("ComputationSolution", "solves linear systems");

  Matrix<Complex_wp, General, BandedCol> mat_lu;
  Vector<int> pivot;

  //Vector<Real_wp> x(nr);
  Vector<int>     eltS(nr); eltS.Zero();
  Vector<Real_wp> locPos(nr); locPos.Zero();

  for (int i =0; i < nr; i++)
    {
      //x(i) = 1.*i/nR;

      while ((eltS(i)<var.mesh.GetNbVertices()) && (var.mesh.Vertex(eltS(i))<r(i)))
	eltS(i)++;

      // Point at a mesh boundary or very close
      if (eltS(i)==0) {
	if (abs(var.mesh.Vertex(0)-r(i))<=10*epsilon_machine)
	  locPos(i) = 0.e0;
	else {
	  cout << "point " << r(i) << " not in the mesh " << endl;
	  abort();
	}}
      else if (eltS(i)==var.mesh.GetNbVertices()) {
	if (abs(var.mesh.Vertex(eltS(i)-1)-r(i)) <= 10*epsilon_machine) {
	  eltS(i)    = var.mesh.GetNbElt()-1;
	  locPos(i) = 1.e0;
	}
	else {
	  cout << "point " << r(i) << " not in the mesh " << endl;
	  abort();
	}}
      // point inside the mesh
      else {
	eltS(i)--;
	locPos(i) = (r(i)-var.mesh.Vertex(eltS(i)))/(var.mesh.Vertex(eltS(i)+1)-var.mesh.Vertex(eltS(i)));
      }
    }

  G.Reallocate(lMaxG, nr);
  dG.Reallocate(lMaxG, nr);


  for (int l=0; l < lMaxG; l++)
    {
      var_chrono.Start("FactorisationMatrix");
      
      // computation and factorisation of finite element matrix
      var.PerformFactorizationStep(mat_lu, pivot, l);
      
      var_chrono.Stop("FactorisationMatrix");
      
      var_chrono.Start("ComputationSolution");

      if (rhs(l)(0).GetM() > 0)
          {
            var.ComputeSolution(mat_lu, pivot, rhs(l)(0));
          }
	var_chrono.Stop("ComputationSolution");

	// 2. Reconstruct the solutions
	VectReal_wp phi;
	VectReal_wp grad_phi;

	for (int i =0; i < nr; i++)
	  {
	    const ElementReference<Dimension1,1>& Fb = var.GetReferenceElement(eltS(i));
	    //int   nLocDof = Fb.GetNbDof();
	    //grad_phi.Reallocate(nLocDof);
	    Fb.ComputeValuesPhiRef(locPos(i),phi);
	    Fb.ComputeGradientPhiRef(locPos(i),grad_phi);
	    for (int IDloc=0; IDloc<phi.GetM(); IDloc++) 
	      {
		G(l,i) += rhs(l)(0)(var.mesh.GetNumberDof(eltS(i),IDloc)) * phi(IDloc);
		dG(l,i) += rhs(l)(0)(var.mesh.GetNumberDof(eltS(i),IDloc)) * grad_phi(IDloc);
	      }
	  }
    }
  G.Write("Ul.dat");
  dG.Write("drUl.dat");
}

  Vector<Matrix<Complex_wp> > Kernel1D::computeKernel(const Vector<int> lKs, const int mK, int rank_proc)
  {
    // Associated Legendre polynomials
    AssociatedLegendrePolynomial<Real_wp> assocLegendre;
    assocLegendre.Init(lMaxG);
    Vector<Vector<VectReal_wp> > P1(theta1.GetM());
    Vector<Vector<VectReal_wp> > P2(theta2.GetM());
    for (int i=0; i < theta1.GetM(); i++)
      {
	assocLegendre.EvaluatePnm(lMaxG,lMaxG,theta1(i),P1(i));
	assocLegendre.EvaluatePnm(lMaxG,lMaxG,theta2(i),P2(i));
      }

    Vector<Matrix<Complex_wp> > K(theta1.GetM());
    for (int i=0; i<theta1.GetM();i++)
      {
	K(i) = Matrix<Complex_wp>(lKs.GetM(), nr);
	K(i).Zero();
      }

    Real_wp plm1, plm2, plmp1, plmp2;
    Complex_wp ylm1, ylm2, ylmp1, ylmp2;
    Complex_wp I_wp = Complex_wp(0.,1.);
    Vector<Vector<Matrix<Real_wp,General,ArrayRowSparse> > > intPlm_m1, intPlm_p1;
      if (typeOfKernel == "UTHETA" or typeOfKernel == "UPHI")
      {
	//Precompute the required integrals
	clock_t begin_time = clock();
	if (rank_proc == 0)
	  {
	    cout << "Precomputing the Plm integrals ..." << endl; 
	  }
	intPlm_m1 = generalGauntNew(lMaxG,lMaxG,lKs,mK, -1);
	if (mK != 0)
	  intPlm_p1 = generalGauntNew(lMaxG,lMaxG,lKs,mK, 1);
	if (rank_proc == 0)
	  {
	    Real_wp seconds = 1.*(clock() - begin_time) / CLOCKS_PER_SEC;
	    cout << "Time to precompute the Plm integrals: " << seconds << endl;
	  }
      }

      Real_wp timeGaunt = 0.;
      clock_t time0;
      for (int l = 0; l<=lMaxG; l++)
	{
	  if (rank_proc == 0)
	    DISP(l);

	  Real_wp alphal = sqrt((4*pi_wp)/(2.*l+1.));
	  for (int m = -l; m <= l; m++)
	    {
	      int mp = mK-m;
	      Vector<Real_wp> g;

	      for (int indlK=0; indlK<lKs.GetM(); indlK++)
		{
		  int lK = lKs(indlK);
		  if (typeOfKernel == "UR" or typeOfKernel == "C")
		    {
		      if (rank_proc == 0)
			time0 = clock();
		      g = GauntVect(lMaxG, lK, l, mp, -mK, m);
		      if (rank_proc == 0)
			timeGaunt += 1.*(clock() - time0) / CLOCKS_PER_SEC;
		    }
		 
		  int  indMin = 0, indMax = lMaxG;
		  if (typeOfKernel == "UR" or typeOfKernel == "C")
		    {
		      indMin =  abs(l-lK);
		      indMax =  min(lMaxG,abs(l+lK));
		    }
		  else if (typeOfKernel == "UPHI")
		    {
		      if ((l+lK)%2==0)
			indMin = 2;
		      else
			indMin = 1;
		    }
		  else if (typeOfKernel == "UTHETA")
		    {
		      // in order to keep only the non-zero terms in the loop (one term over two is zero)
		      if ((l+lK)%2==0)
			indMin = 1;
		      else
			indMin = 0;
		    }

		  for (int lp = indMin; lp <=indMax; lp+=2)
		    {
		      if (abs(mp) <= lp)
			{
			  Real_wp alphalp = sqrt((4*pi_wp)/(2.*lp+1.));

			  Complex_wp Ij = 0.;
			  if (typeOfKernel == "UR" or typeOfKernel == "C")
			    {
			      Ij = pow(-1,mK) * g(lp);
			    }
			  else if (typeOfKernel == "UTHETA")
			    {
			      Real_wp gm1, gp1;
			      if (lp <= l)
				{
				  gm1 = intPlm_m1(lp)(indlK)(l,lMaxG+m);
				  if (mK == 0)
				    gp1 = -intPlm_m1(lp)(indlK)(l,lMaxG-m);
				  else
				    gp1 = intPlm_p1(lp)(indlK)(l,lMaxG+m);
				}
			      else
				{
				  gm1 = intPlm_m1(l)(indlK)(lp,lMaxG+mp-1);
				  if (mK == 0)
				    gp1 = -intPlm_m1(l)(indlK)(lp,lMaxG-mp-1);
				  else
				    gp1 = intPlm_p1(l)(indlK)(lp,lMaxG+mp+1);
				}
			      Ij = (-gm1 * sqrt((lp+mp)*(lp-mp+1)) + gp1 * sqrt((lp+mp+1)*(lp-mp)) ) / (2.*sqrt(2.*pi_wp));
			    }
			  else if (typeOfKernel == "UPHI")
			    {
			      Real_wp gm1, gp1;
			      if (lp-1 <= l)
				{
				  gm1 = intPlm_m1(lp-1)(indlK)(l,lMaxG+m);
				  if (mK == 0)
				    gp1 = -intPlm_m1(lp-1)(indlK)(l,lMaxG-m);
				  else
				    gp1 = intPlm_p1(lp-1)(indlK)(l,lMaxG+m);
				}
			      else
				{
				  gm1 = intPlm_m1(l)(indlK)(lp-1,lMaxG+mp-1);
				  if (mK == 0)
				    gp1 = -intPlm_m1(l)(indlK)(lp-1,lMaxG-mp-1);
				  else
				    gp1 = intPlm_p1(l)(indlK)(lp-1,lMaxG+mp+1);
				}

			      /*
				if (l==5 and gm1 != 0)
				{
				cout << "m1 = " << m << ", gm1 = " << gm1 << ", gm1_single = " << generalGauntSingle(l,lp-1,lK,m,mp-1,mK) << endl; //" ,symmetry = " << intPlm_m1(l)(0)(lp,lMaxG+mp-1) << endl;
				cout << "m1 = " << m << ", gp1 = " << gp1 << ", gp1_single = " << generalGauntSingle(l,lp-1,lK,m,mp+1,mK) << endl; //<< " ,symmetry = " << intPlm_p1(l)(0)(lp,lMaxG+mp+1) << endl;
				}
			      */

			      Ij = I_wp* (gm1 * sqrt( (2.*lp+1.)*(lp+mp-1)*(lp+mp)/ (2.*lp-1) ) + gp1 * sqrt( (2.*lp+1)*(lp-mp-1)*(lp-mp) / (2.*lp-1) ) ) / (2.*sqrt(2.*pi_wp));
			    }

			  for (int i=0; i < theta1.GetM(); i++)
			    {
			      if (m > 0)
				{
				  plm1 = P1(i)(l)(abs(m));
				  plm2 = P2(i)(l)(abs(m));
				}
			      else
				{
				  plm1 = pow(-1.,m) * P1(i)(l)(abs(m));
				  plm2 = pow(-1.,m) * P2(i)(l)(abs(m));
				}
			      ylm1 = plm1 * exp(I_wp*m*phi1(i));
			      ylm2 = plm2 * exp(I_wp*m*phi2(i));
			      
			      if (mp > 0)
				{
				  plmp2 = P2(i)(lp)(abs(mp));
				  plmp1 = P1(i)(lp)(abs(mp));
				}
			      else
				{
				  plmp2 = pow(-1.,mp) * P2(i)(lp)(abs(mp));		
				  plmp1 = pow(-1.,mp) * P1(i)(lp)(abs(mp));
				}
			      ylmp2 = plmp2 * exp(I_wp*mp*phi2(i));
			      ylmp1 = plmp1 * exp(I_wp*mp*phi1(i));



			      
			      for (int j=0; j < nr; j++)
				{
				  Complex_wp partf, partg;
				  if (typeOfObservable == "XS")
				    {
				      if (typeOfKernel == "UR")
					{
					  partf = 2.*I_wp*rho(j) * omega* dC(j,lp) *  G(j,l);
					  //partg = -2*I_wp*rho(j) * omega* conj(dC(j,lp)) *  conj(G(j,l));
					}
				      else if (typeOfKernel == "UTHETA" or typeOfKernel == "UPHI")
					{
					  partf = 2.*I_wp*rho(j) * omega* C(j,lp) *  G(j,l) / (r(j)*Rsun);
					  //partg = -2*I_wp*rho(j) * omega* conj(C(j,lp)) *  conj(G(j,l));
					}
				      else if (typeOfKernel == "C")
					{
					  partf = rho(j) / c(j) * omega* (omega + 2.*I_wp*gamma(0)) * C(j,lp) *  G(j,l);
					  //partg = rho(j) / c(j) * omega* (omega - 2.*I_wp*gamma(0)) * conj(C(j,lp)) *  conj(G(j,l));
					}
				      partg = conj(partf);
				    }
				  else
				    {
				      partf = fllpr(j)(l,lp);
				      partg = gllpr(j)(l,lp);
				    }
				  K(i)(indlK,j) += alphal * alphalp * (Ij * partf * conj(ylm2) * conj(ylmp1) + conj(Ij) * partg* conj(ylm1)*conj(ylmp2));					  
				}
			    }
			}
		    }
		}
	    }
	}
      if (rank_proc == 0)
	cout << "Time to compute the Gaunt integrals : " << timeGaunt << endl;
      return K;  
  }

  void Kernel1D::compute(const int mK)
  {
    // General function to compute kernels. It distributes the jobs over the different cores, call the appropriate function and gather the results.

    int nb_proc = 1, rank_proc = 0;
#ifdef SELDON_WITH_MPI    
    MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#endif

    // Separate values of lK  on each processor
    Vector<Vector<int> > listL(nb_proc);
    
    /*
    int aveJobs = (lMaxK+1)/nb_proc;
    int extraJobs = lMaxK+1-aveJobs*nb_proc;
    int cpt = 0;
    for (int i=0; i < nb_proc; i++)
      {
	int nbJobs = aveJobs;
	if (i < extraJobs)
	  nbJobs++;
	for (int j=0; j < nbJobs; j++)
	  {
	    listL(i).PushBack(cpt);
	    cpt++;
	  }
      }
    */
    
    
    
    for (int l=0; l <= lMaxK; l++)
      {
	int proc_crt = l%nb_proc;
	//	listL(proc_crt).PushBack(l); 
	listL((proc_crt+1)%nb_proc).PushBack(l); // to guarantee that proc0 has the less jobs
      }
    
    Vector<int64_t> xtmp;
    int tag = 23;
    MPI_Status status;

    if (typeOfObservable != "XS")
      {
	if (rank_proc == 0)
	  {
	    clock_t begin_time = clock();
	    for (int i=1; i < nb_proc; i++)
	      {
		for (int indz=0; indz<nr; indz++)
		  {
		    MpiSend(MPI_COMM_WORLD, fllpr(indz).GetData(), xtmp, fllpr(indz).GetDataSize(), i, tag);
		    MpiSend(MPI_COMM_WORLD, gllpr(indz).GetData(), xtmp, gllpr(indz).GetDataSize(), i, tag);
		  }
	      }
	    Real_wp seconds = 1.*(clock() - begin_time) / CLOCKS_PER_SEC;
	    cout << "Time to send data to other processors: " << seconds << endl;
	  }
	else
	  {
	    for (int indz=0; indz<nr; indz++)
	      {
		MpiRecv(MPI_COMM_WORLD, fllpr(indz).GetData(), xtmp, fllpr(indz).GetDataSize(), 0, tag, status);
		MpiRecv(MPI_COMM_WORLD, gllpr(indz).GetData(), xtmp, gllpr(indz).GetDataSize(), 0, tag, status);
	      }
	  }
      }

    clock_t begin_time = clock();
    Vector<Matrix<Complex_wp> > ansCrt = computeKernel(listL(rank_proc), mK, rank_proc);
    if (rank_proc == 0)
      {
	Real_wp seconds = 1.*(clock() - begin_time) / CLOCKS_PER_SEC;
	cout << "Time to compute the kernels on proc 0: " << seconds << endl;
      }

    Vector<Matrix<Complex_wp> > Klm(theta1.GetM());
    if (rank_proc == 0)
      {
	for (int i=0; i < theta1.GetM(); i++)
	  {
	    Klm(i) = Matrix<Complex_wp>(lMaxK+1, nr); 
	    Klm(i).Zero();
	  }
      }
    
    //#ifdef SELDON_WITH_MPI
    //MPI_Barrier(MPI_COMM_WORLD);
    //#endif
    
    begin_time = clock();

    for (int i=0; i < theta1.GetM(); i++)
      {
	MpiGather(MPI_COMM_WORLD, ansCrt(i).GetData(), xtmp, Klm(i).GetData(), ansCrt(i).GetDataSize(), 0);
	if (rank_proc == 0)
	  {
	    // Regroups the values of lK from the different processors and write the answer
	    string filename;
	    
	    Matrix<Complex_wp> Klm2(Klm(i).GetM(), Klm(i).GetN()); 
	    Klm2.Zero();
	    int cpt = 0;
	    for (int k=0; k < listL.GetM(); k++)
	      {
		for (int n=0; n < listL(k).GetM(); n++)
		  {
		    for (int j=0; j <nr; j++)
		      {
			Klm2(listL(k)(n), j) = Klm(i)(cpt,j);
		      }
		    cpt++;
		  }
		//cpt += listL(k).GetM();
	      }
	    
	    //filename = dirOut + "kernel_lm_th1_" + to_str(int(round(theta1(i)*180./pi_wp))) + "_phi1_" + to_str(int(round(phi1(i)*180./pi_wp))) + "_th2_" + to_str(int(round(theta2(i)*180./pi_wp))) + "_phi2_" + to_str(int(round(phi2(i)*180./pi_wp))) + "_" + typeOfKernel;
	    filename = dirOut + "kernel_lm_th1_" + to_str(theta1(i)*180./pi_wp) + "_phi1_" + to_str(int(round(phi1(i)*180./pi_wp))) + "_th2_" + to_str(theta2(i)*180./pi_wp) + "_phi2_" + to_str(int(round(phi2(i)*180./pi_wp))) + "_" + typeOfKernel;
	    if (typeOfObservable == "XS")
	      filename += "_f" + to_str(indFreq);
	    filename += ".dat";
	    Klm2.Write(filename);
	  }
      }
    
    if (rank_proc == 0)
      {
	Real_wp seconds = 1.*(clock() - begin_time) / CLOCKS_PER_SEC;
	cout << "Time to gather the results and write data: " << seconds << endl;
      }
  }


  void  Kernel1D::computeKernel3D_depth(const Vector<int> listIndz)
  {

    int Nth = 1001;
    Vector<Real_wp> theta(Nth);
    Linspace(0.,pi_wp,Nth,theta);
    int Nphi = 2001;
    Vector<Real_wp> phi(Nphi);
    Linspace(0.,2.*pi_wp,Nphi,phi);

    // Associated Legendre polynomials
    AssociatedLegendrePolynomial<Real_wp> assocLegendre;
    assocLegendre.Init(lMaxG);
    Vector<Vector<VectReal_wp> > P1(theta1.GetM());
    Vector<Vector<VectReal_wp> > P2(theta2.GetM());

    Vector<Vector<VectReal_wp> > P(theta.GetM());
    for (int i=0; i < theta1.GetM(); i++)
      {
	assocLegendre.EvaluatePnm(lMaxG,lMaxG,theta1(i),P1(i));
	assocLegendre.EvaluatePnm(lMaxG,lMaxG,theta2(i),P2(i));
      }

    for (int i=0; i < theta.GetM(); i++)
	assocLegendre.EvaluatePnm(lMaxG,lMaxG,theta(i),P(i));

    Vector<Matrix<Complex_wp> >  K(theta1.GetM());
    Vector<Matrix<Complex_wp> >  partG(theta1.GetM());
    Vector<Matrix<Complex_wp> >  partC(theta1.GetM());
    Vector<Matrix<Complex_wp> >  partGc(theta1.GetM());
    Vector<Matrix<Complex_wp> >  partCc(theta1.GetM());

    for (int i=0; i<theta1.GetM();i++)
      {
	K(i) = Matrix<Complex_wp>(Nth,Nphi);  
	partG(i) = Matrix<Complex_wp>(Nth,Nphi); 
	partC(i) = Matrix<Complex_wp>(Nth,Nphi);  
	partGc(i) = Matrix<Complex_wp>(Nth,Nphi); 
	partCc(i) = Matrix<Complex_wp>(Nth,Nphi); 	  
      }

    Real_wp plm1, plm2, plm;
    Complex_wp ylm1, ylm2, ylm;
    Complex_wp I_wp = Complex_wp(0.,1.);


    for (int ind=0; ind < listIndz.GetM(); ind++)
      {
	int indz = listIndz(ind);
	for (int i=0; i<theta1.GetM();i++)
	  {
	    K(i).Zero();
	    partG(i).Zero();
	    partC(i).Zero();
	    partGc(i).Zero();
	    partCc(i).Zero();	  
	  }
 
	for (int l = 0; l<=lMaxG; l++)
	  {
	    cout << "ind depth = " << indz << " , l = " << l << endl;
	    for (int m = 0; m <= l; m++)
	      {
		
		for (int k=0; k< theta1.GetM(); k++)
		  {
		    plm1 = P1(k)(l)(abs(m));
		    plm2 = P2(k)(l)(abs(m));
		    
		    ylm1 = plm1 * exp(I_wp*m*phi1(k));
		    ylm2 = plm2 * exp(I_wp*m*phi2(k));
		    
		    for (int i=0; i < Nth; i++)
		      {
			plm = P(i)(l)(abs(m));
			
			for (int j=0; j < Nphi; j++)
			  {
			    ylm = plm * exp(I_wp*m*phi(j));
			    
			    
			    int factor;
			    if (m==0)
			      factor = 1;
			    else
			      factor = 2.;
			    
			    partG(k)(i,j) += sqrt((4*pi_wp)/(2.*l+1.)) * factor*real(conj(ylm2)* ylm) *  G(indz,l);
			    partC(k)(i,j) += sqrt((4*pi_wp)/(2.*l+1.)) * factor*real(conj(ylm1)* ylm) * dC(indz,l);
			    partGc(k)(i,j) += sqrt((4*pi_wp)/(2.*l+1.)) * factor*real(ylm1* conj(ylm)) *  conj(G(indz,l));
			    partCc(k)(i,j) += sqrt((4*pi_wp)/(2.*l+1.)) * factor*real(ylm2* conj(ylm)) * conj(dC(indz,l));				    
			    
			  }
		      }		 
		  }
	      }
	  }
	
  
	for (int k=0; k< theta1.GetM(); k++)
	  {
	    
	    for (int i=0; i < Nth; i++)
	      {
		for (int j=0; j < Nphi; j++)
		  {		      
		    K(k)(i,j) = partG(k)(i,j)*partC(k)(i,j)-partGc(k)(i,j)*partCc(k)(i,j);
		    K(k)(i,j) *= 2*I_wp*rho(indz)*omega;
		  }
	      }
	    string filename = dirOut + "kernel_" + to_str(k) + "_z" + to_str(indz) + "_" + typeOfKernel + ".dat";		
	    K(k).Write(filename);
	
	  }
      }

  }


  void Kernel1D::computeKernel3D()
  {
    // General function to compute kernels. It distributes the jobs over the different cores, call the appropriate function and gather the results.

    int nb_proc = 1, rank_proc = 0;
#ifdef SELDON_WITH_MPI    
    MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_proc);
#endif

    // Separate depths  on each processor
    Vector<Vector<int> > listIndz(nb_proc);
    
    int aveJobs = nr/nb_proc;
    int extraJobs = nr-aveJobs*nb_proc;
    int cpt = 0;
    for (int i=0; i < nb_proc; i++)
      {
	int nbJobs = aveJobs;
	if (i < extraJobs)
	  nbJobs++;
	for (int j=0; j < nbJobs; j++)
	  {
	    listIndz(i).PushBack(cpt);
	    cpt++;
	  }
      }
    
    
    computeKernel3D_depth(listIndz(rank_proc));

  }

}

#define MONTJOIE_FILE_KERNEL1D_CXX
#endif
