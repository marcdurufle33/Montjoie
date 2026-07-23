#ifndef MONTJOIE_FILE_SPACE_STRING_SCHEME_CXX

namespace Montjoie
{
  
  //! default constructor
  template<int nb_base, int nb_quad, class TypeEquation>
  SpaceStringScheme<nb_base, nb_quad, TypeEquation>::SpaceStringScheme()
  {
    Nx			= 0;
    Nx_interp			= 0;
    Deltax			= 0.0;
    invDeltax			= 0.0;
    point_observation         = 0.0;
    elt_obs			= 0;
    nb_elt_obs		= 0;
    use_lobatto_quadrature = false;

    file_output_sismo_string = "Sismo_String";
    file_output_last_point = "Sismo_Last_Point";
    file_output_crush = "Sismo_Crush";
    file_output_interp = "Val";
  }
  
  
  //! initialisation des sorties de type "sismogramme"
  template<int nb_base, int nb_quad, class TypeEquation>
  void SpaceStringScheme<nb_base, nb_quad, TypeEquation>
  ::InitBuffer(const string& DOSSIER, int size, int numero_corde, bool remove_file,
               bool double_prec, bool binary)
  {
    output_string.SetDoublePrecision(double_prec);
    output_last_point.SetDoublePrecision(double_prec);
    output_crush.SetDoublePrecision(double_prec);
    
    output_string.SetBinary(binary);
    output_last_point.SetBinary(binary);
    output_crush.SetBinary(binary);
    
    size_buffer = size;
    output_string.Init(DOSSIER+file_output_sismo_string + to_str(numero_corde) +".txt",
                       size_buffer, remove_file);
    
    output_last_point.Init(DOSSIER + file_output_last_point + to_str(numero_corde) +".txt",
                           size_buffer, remove_file);
    
    output_crush.Init(DOSSIER+ file_output_crush + to_str(numero_corde) +".txt",
                      size_buffer, remove_file);
  }
  
  
  //! construction des elements finis et du maillage
  template<int nb_base, int nb_quad, class TypeEquation>
  void SpaceStringScheme<nb_base, nb_quad, TypeEquation>::ConstructMesh(const Real_wp& L)
  {
    // construit les fonctions de base, la quadrature
    base.ConstructQuadrature(nb_base-1, Globatto<Real_wp>::QUADRATURE_LOBATTO); 
    if(use_lobatto_quadrature)
      quad.ConstructQuadrature(nb_quad-1, Globatto<Real_wp>::QUADRATURE_LOBATTO);
    else
      quad.ConstructQuadrature(nb_quad-1, Globatto<Real_wp>::QUADRATURE_GAUSS);
    
    // construction du maillage
    mesh.CreateRegularMesh(0,L,Nx,1);
    mesh.SetOrder(nb_base-1);
    mesh.NumberMesh(false);
    
    // pas de maillage et son inverse
    Deltax = L/(Nx-1);
    invDeltax = 1.0/Deltax;
    
    // calcul de phi_i(\xi_j), d/dx(phi_i)(xi_j), etc
    // boucle sur les fct de base
    for (int i = 0; i < nb_base; i++)
      {
	//boucle sur les points d'integration
	for (int j = 0; j < nb_quad; j++)
	  {
	    valPhi(i,j) = base.EvaluatePhi(i, quad.Points(j));
	    valPhiWeight(i,j) = valPhi(i,j)*quad.Weights(j);
	    Weights (j) = quad.Weights(j);
	  }
	
	valPhiLastPoint(i)=base.EvaluatePhi(i, 1.0);
      }  
    
    // boucle sur les fct de base
    for (int i = 0; i < nb_base; i++)
      { 
	//boucle sur les points d'integration
        for (int j = 0; j < nb_quad; j++)
	  {
	    valdPhi(i,j) = base.EvaluatePhiGrad(i, quad.Points(j));
	    valdPhiWeight(i,j) = valdPhi(i,j)*quad.Weights(j);
	  }
	
        valdPhiLastPoint(i) = base.EvaluatePhiGrad(i, Real_wp(1.0))*invDeltax;
      } 
    
    for (int i = 0; i < nb_base; i++)
      for (int j = 0; j < nb_quad; j++)
          for (int i2 = 0; i2 < nb_base; i2++)
	    {
	      valDoubledPhiWeight(i,j,i2) = valdPhiWeight(i,j)*valdPhi(i2,j)*invDeltax;
	      valDoublePhiWeight(i,j,i2) = valPhiWeight(i,j)*valPhi(i2,j)*Deltax;
	      valMixteWeight(i,j,i2) = valPhiWeight(i,j)*valdPhi(i2,j);
	    }
    
    MassMatrixBase.Fill(0);
      
    //Petite mass matrix
    for (int i = 0; i < nb_base; i++)
      for (int j = i; j < nb_base; j++)
	      for (int k = 0; k < nb_quad; k++)
	        MassMatrixBase(i,j) += quad.Weights(k)*valPhi(i,k)*valPhi(j,k);
  }
  
  
  //! calcule \int phi_j \delta(x-x_0) ou delta est la fonction de distribution du marteau
  /*!
    \param[in] hammer marteau considere (on recupere le point d'impact et l'epaisseur du marteau)
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  void SpaceStringScheme<nb_base, nb_quad, TypeEquation>
  ::ConstruitImpact(Hammer & hammer)
  {
    ComputeRepartitionVector(repartition_marteau, repartition_marteau_sparse,
                             hammer.epaisseur_marteau, hammer.point_impact_marteau,
                             hammer);
  }
  
  
  //! calcule \int phi_j \delta(x-x_0) ou delta est la fonction de distribution du marteau
  /*!
    \param[out] repartition vecteur plein contenant \int phi_j \delta(x-x_0)
    \param[out] repartition_sparse vecteur creux contenant \int phi_j \delta(x-x_0)
    \param[in] delta largeur de la fonction de distribution delta
    \param[in] point_impact centre x0 de la fonction de distrubtion
    \param[in] hammer marteau considere
   */
  template<int nb_base, int nb_quad, class TypeEquation> template<class GenericPb>
  void SpaceStringScheme<nb_base, nb_quad, TypeEquation>
  ::ComputeRepartitionVector(VectReal_wp& repartition,
                             Vector<Real_wp, VectSparse>& repartition_sparse,
                             const Real_wp& delta, const Real_wp& point_impact, GenericPb& hammer)
  { 
    repartition.Reallocate(mesh.GetNbDof());
    repartition.Fill(Real_wp(0));
 
    if (delta != 0)
    {
      // calcul de \int delta \varphi ou delta est la fonction de repartition
      for(int ne = 0; ne < mesh.GetNbElt(); ne++)
      {
        for(int p = 0; p < nb_base; p++)
        {
          int pp = mesh.GetNumberDof(ne,p);
          Real_wp res = 0;
          for (int g = 0; g < nb_quad; g++)
          {
            Real_wp x = (ne + quad.Points(g))*Deltax;
            res += quad.Weights(g)*hammer.Fct_Delta(x - point_impact)*valPhi(p,g);
          }
              
          if(abs(res) > 1e-12)
            repartition(pp) += Deltax*res;
        }
      }
    }
    else
    {
      // cas d'un Dirac
      // elt_x0 : element ou x0 appartient
      // r : coordonnee locale dans l'element
      int elt_x0 = toInteger(point_impact/Deltax);
      Real_wp r = (point_impact - elt_x0*Deltax)/Deltax;
        
      // on remplit repartition_marteau
      for(int j = 0; j < nb_base; j++)
      {
        int jj = mesh.GetNumberDof(elt_x0, j);
        repartition(jj) = base.EvaluatePhi(j, r);
      }
    }
    
    // conversion vers un vecteur creux
    Copy(repartition, repartition_sparse);
  }
  
  
  //! calcul de la grille d'interpolation afin de calculer l'interpolation de u sur la grille
  /*!
    \param[in] L longueur de la grille d'interpolation
    Les points de la grille sont localises sur le maillage, et on calcule
    MatInterp(i, j) = phi_j(xi_i)
    MatInterp est une matrice creuse
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  void SpaceStringScheme<nb_base, nb_quad, TypeEquation>::ConstruitInterp(const Real_wp & L)
  {
    if(Nx_interp !=0 )
      {
	
	grid_interp.Init(0, L, Nx_interp);
	// en utilisant le grid_interp de montjoie
	grid_interp.LocalizePoints(mesh);
	// necessite qu'il y ait eu un grid_interp.Init(xmin,xmax,N);
	
	// MatInterp(j,l) = \phi_j (x_interp_l)
	MatInterp.Reallocate(grid_interp.GetNbPointsGrid(), mesh.GetNbDof());
	
	for (int l = 0; l < grid_interp.GetNbPointsGrid(); l++) 
	  {
	    int ne = grid_interp.GetElementNumber(l);
	    Real_wp x = grid_interp.GetLocalCoordinate(l);
	    for (int jj = 0; jj < nb_base; jj++)
	      {
		int j = mesh.GetNumberDof(ne, jj);
		MatInterp.Get(l, j) = base.EvaluatePhi(jj, x);
	      }
	    
	  }
      }
  }
  
  
  //! Localise sur le maillage le point ou le sismogramme est demande
  /*!
    Le point point_observation est localise sur le maillage
    et on calcule observ = phi_i(point_observation)
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  void SpaceStringScheme<nb_base, nb_quad, TypeEquation>::ConstruitObserv()
  {
    // Construit nb_elt_x0, elt_x0, impact
    nb_elt_obs = 1;
    elt_obs = toInteger(point_observation/Deltax);
    // coordonnee locale dans l'element, r appartient a [0, 1]
    Real_wp r = (point_observation-elt_obs*Deltax)/Deltax;
    
    if(r == 0)
      if(elt_obs != 0)
	nb_elt_obs = 2;
    
    // remplir observ
    observ.Reallocate(nb_base);
    observ.Fill(Real_wp(0));
    if (nb_elt_obs==1)
      {
	for(int j = 0; j < nb_base; j++)
	  observ(j) += base.EvaluatePhi(j, r);
      }
    else
      observ(0) = 1;
    
  } 
  
  
  //! calcule l'interpolation de u sur la grille, et imprime la solution sur le disque
  /*!
    \param[in] n numero de l'instantane a ecrire sur le disque
    \param[in] var objet contenant la solution
    \param[in] numero_corde numero de la corde dont on veut connaitre le deplacement
   */
  template<int nb_base, int nb_quad, class TypeEquation>  template<class GenericPb>
  void SpaceStringScheme<nb_base,nb_quad,TypeEquation>
  ::WriteInterp(int n, GenericPb & var, int numero_corde)
  
  {
    int dim = TypeEquation::dimension;
    Matrix<Real_wp> tmp;
    tmp.Reallocate(Nx_interp, dim);
    tmp.Fill(Real_wp(0));
    // on calcule l'interpolation de u sur les points de la grille
    for (int k = 0; k < dim; k++)
      for (int l = 0; l < Nx_interp; l++)
	for (int p = 0; p < MatInterp.GetRowSize(l); p++)
	  {
	    int indice_ddl = MatInterp.Index(l,p);	    
	    tmp(l,k) += var.GetValue(numero_corde, k, indice_ddl, 0, var.val)*MatInterp.Value(l, p);
	  }
    
    // nom du fichier de sortie
    string name_file = var.DOSSIER + file_output_interp + to_str(n) + "string" + to_str(numero_corde) + ".txt";
    remove(name_file.data());
    ofstream objetfichier;
    objetfichier.open(name_file.data(), ios::out); // on ouvre le fichier en ecriture
    if (objetfichier.bad()) // permet de tester si le fichier s'est ouvert sans probleme
      cout << "pb dans le fichier val "<<endl;
    else
      {
	objetfichier.precision(15);
	objetfichier << tmp << endl;//		
      }
    
    objetfichier.close(); //on ferme le fichier
    tmp.Clear();
  }
  

  //! calcule l'interpolation de u sur le sismo, et imprime la solution
  /*!
    \param[in] temps temps actuel
    \param[in] var objet contenant la solution
    \param[in] numero_corde numero de la corde dont on veut connaitre le deplacement
   */  
  template<int nb_base, int nb_quad, class TypeEquation> template<class GenericPb>
  void SpaceStringScheme<nb_base, nb_quad, TypeEquation>
  ::WriteSismo(Real_wp temps, GenericPb & var, int numero_corde)
  {
    int dim = TypeEquation::dimension;

    // on calcule l'interpolation de u sur le point d'observation
    TinyVector<Real_wp, TypeEquation::dimension+1> temp;
    temp(0) = temps;
    for (int k = 0; k < dim; k++)
      {
	Real_wp ux0p = 0.0;
	for(int p = 0; p < nb_base; p++)
	  {
	    int pp = mesh.GetNumberDof(elt_obs, p);
	    ux0p += var.GetValue(numero_corde, k, pp, 0, var.val)*observ(p);
	  }
	
	temp(k+1) = ux0p;
      }
    
    // on ecrit le sismo sur le disque
    // (l'ecriture peut etre eventuellement differe a cause du buffer)
    output_string.AddTinyVect(temp);
    
    // on fait la meme chose pour le dernier point
    temp.Fill(Real_wp(0));
    temp(0) = temps;
    for (int k = 0; k < dim; k++)
      temp(k+1) = var.GetValue(numero_corde, k, mesh.GetNbDof()-1, 0, var.val);
    
    output_last_point.AddTinyVect(temp);
  }
  
  
  //! calcule res_j = res_j + \int f phi_j
  /*!
    \param[in] valeurs f(xi_j) (valeur de la fonction f aux points de quadrature)
    \param[in] res composantes sur les fonctions de base
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  void SpaceStringScheme<nb_base, nb_quad, TypeEquation>
  ::CalculeIntPhi(const TinyVector<Real_wp, nb_quad>& valeurs,
                  TinyVector<Real_wp, nb_base> & res)
  {
    MltAdd(Deltax, valPhiWeight, valeurs, res);
  }
  
  
  //! calcule res_j = res_j + \int f d/dx(varphi_j)
  /*!
    \param[in] valeurs f(xi_j) (valeur de la fonction f aux points de quadrature)
    \param[in] res composantes sur les fonctions de base
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  void SpaceStringScheme<nb_base, nb_quad, TypeEquation>
  ::CalculeIntGrad(const TinyVector<Real_wp, nb_quad>& valeurs, TinyVector<Real_wp, nb_base> & res)
  {
    MltAdd(1.0, valdPhiWeight, valeurs, res);
  }
  
  
  //! calcule res_j = res_j + \int f phi_j
  /*!
    \param[in] valeurs f(xi_j) (valeur de la fonction f aux points de quadrature)
    \param[in] res composantes sur les fonctions de base
    On fait cette operation sur chaque ligne de res (qui est en fait une matrice)
   */
  template<int nb_base, int nb_quad, class TypeEquation> template<int dimension>
  void SpaceStringScheme<nb_base, nb_quad, TypeEquation>
  ::CalculeIntPhi(const TinyMatrix<Real_wp, General, dimension, nb_quad>& valeurs,
		  TinyMatrix<Real_wp, General, dimension, nb_base> & res)
  {
    for (int k = 0; k < dimension; k++)
      {
	TinyVector<Real_wp, nb_base> res_ligne;
	TinyVector<Real_wp, nb_quad> ligne;
	
	GetRow(valeurs, k, ligne);
	MltAdd(Deltax, valPhiWeight, ligne, res_ligne);
	SetRow(res_ligne, k, res);
      }
  }
  
  
  //! calcule res_j = res_j + \int f d/dx(varphi_j)
  /*!
    \param[in] valeurs f(xi_j) (valeur de la fonction f aux points de quadrature)
    \param[in] res composantes sur les fonctions de base
    On fait cette operation sur chaque ligne de res (qui est en fait une matrice)
   */  
  template<int nb_base, int nb_quad, class TypeEquation> template<int dimension>
  void SpaceStringScheme<nb_base, nb_quad, TypeEquation>
  ::CalculeIntGrad(const TinyMatrix<Real_wp, General, dimension, nb_quad>& valeurs,
		   TinyMatrix<Real_wp, General, dimension, nb_base> & res)
  {
    for (int k = 0; k < dimension; k++)
      {
	TinyVector<Real_wp, nb_base> res_ligne;
	TinyVector<Real_wp, nb_quad> ligne;
	
	GetRow(valeurs, k, ligne);
	MltAdd(1.0, valdPhiWeight, ligne,res_ligne);
	SetRow(res_ligne, k, res);
      }
  }
  
  
  //! computes the values of u on quadrature points of an element
  /*!
    \param[in] numero_corde numero de la corde consideree
    \param[in] ne numero de l'element ou on veut connaitre la valeur de u
                     sur les points de quadrature
    \param[in] instant la solution est dans val(:, instant)
    \param[in] val solution
    \param[in] var objet contenant les cordes
    \param[out] TM valeurs de u sur les points de quadrature
   */
  template<int nb_base, int nb_quad, class TypeEquation> template<int dimension, class GenericPb>
  void SpaceStringScheme<nb_base, nb_quad, TypeEquation>
  ::ComputeTM(int numero_corde, int ne, int instant, Matrix<Real_wp> & val,
	      GenericPb & var, TinyMatrix<Real_wp, General, dimension, nb_quad> & TM )
  {
    TinyVector<Real_wp, nb_quad> quad_values_ligne;
    // boucles sur les composantes de la solution
    for (int k = 0; k < dimension; k++)
      {
	quad_values_ligne.Zero();
	// we store the nodal values of this component
	TinyVector<Real_wp, nb_base> nodal_values;
	for (int j = 0; j < nb_base; j++)
	  {
	    int jj = mesh.GetNumberDof(ne,j);
	    nodal_values(j) = var.GetValue(numero_corde, k, jj, instant, val);
	  }
	
	//we compute the values at the quad points
	MltTrans(valPhi, nodal_values, quad_values_ligne);
	
	// we store it in TM
	SetRow(quad_values_ligne, k, TM);
      }
  }
  

  //! computes the derivative of u on quadrature points of an element
  /*!
    \param[in] numero_corde numero de la corde consideree
    \param[in] ne numero de l'element ou on veut connaitre
                  la derivee de u sur les points de quadrature
    \param[in] instant la solution est dans val(:, instant)
    \param[in] val solution u
    \param[in] var objet contenant les cordes
    \param[out] dTM valeurs de la derivee de u sur les points de quadrature
   */   
  template<int nb_base, int nb_quad, class TypeEquation> template<int dimension, class GenericPb>
  void SpaceStringScheme<nb_base, nb_quad, TypeEquation>
  ::ComputedTM(int numero_corde, int ne, int instant, Matrix<Real_wp>& val,
	       GenericPb & var, TinyMatrix<Real_wp, General, dimension, nb_quad> & dTM )
  {    
    TinyVector<Real_wp,nb_quad> quad_values_ligne;
	
    Real_wp invDeltax = var.piano_strings(numero_corde).invDeltax;
    
    for (int k = 0; k < dimension; k++)
      {
	quad_values_ligne.Zero();
	// on extrait les valeurs nodales sur l'element
	TinyVector<Real_wp, nb_base> nodal_values;
	for (int j = 0; j < nb_base; j++)
	  {
	    int jj = mesh.GetNumberDof(ne, j);
	    nodal_values(j) = var.GetValue(numero_corde, k, jj, instant, val);
	  }
	
	// we compute the derivative at the quad points
	MltTrans(valdPhi,nodal_values,quad_values_ligne);
	
	// we store it in dTM
	SetRow(invDeltax*quad_values_ligne, k, dTM);
      }
  }
  
}

#define MONTJOIE_FILE_SPACE_STRING_SCHEME_CXX
#endif
