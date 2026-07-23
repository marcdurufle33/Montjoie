#ifndef MONTJOIE_FILE_TRANSMISSION_MODEL_HELMHOLTZ_CXX

#include "TransmissionModelHelmholtzInline.cxx"

namespace Montjoie
{

  //! Default constructor
  template<class Complexe, class Dimension>
  void VarTransmission_Helm<Complexe, Dimension>::InitDefaultValues()
  {
    ModeleCentre = false ;
    CasDirichlet = false;
    ModeleTest = false;
    //ModeleTest = true;

    AalphaDelta = 0;
    B1alphaDelta = 0;
    B2alphaDelta = 0;
    A = 0;
    B = 0;
    C = 0;
    mu0 = var_helm.mu0;
    rho0 = var_helm.rho0;
  }
  
  
  //! we read data file
  template<class Complexe, class Dimension>
  void VarTransmission_Helm<Complexe, Dimension>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarTransmission_Base<Dimension>::SetInputData(description_field, parameters);
    
    if (!description_field.compare("ModeleCentre"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarTransmission_Helm" << endl;
	    cout << "ModeleCentre needs more parameters, for instance :" << endl;
	    cout << "ModeleCentre = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (!parameters(0).compare("YES"))
          ModeleCentre = true ;
        else
          ModeleCentre = false ;
      }
    else if (!description_field.compare("CasDirichlet"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of VarTransmission_Helm" << endl;
	    cout << "CasDirichlet needs more parameters, for instance :" << endl;
	    cout << "CasDirichlet = YES" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

        if (!parameters(0).compare("YES"))
          CasDirichlet = true ;
        else
          CasDirichlet = false ;
      }
    else if (!description_field.compare("ReferenceTransmission"))
      {
        int pos = -1;
        int ref_max = 0;
        for (int i = 0; i < parameters.GetM(); i++)
          if ( !  parameters(i).compare("PARAM")) 
            pos = i;
        
        if (pos <= 0)
          {
	    cout << "Needs references before keyword PARAM in ReferenceTransmission" << endl;
            abort();
          }
        
        IVect ref_(pos);
        for (int i = 0; i < pos; i++)
          {
            ref_(i) = to_num<int>(parameters(i));
            ref_max = max (ref_max,ref_(i));
          }
        
        if (ref_max > var_problem.mesh.GetNbReferences())
          var_problem.mesh.ResizeNbReferences(ref_max);
        
	Vector<Vector<Complexe> >& Param_condition = var_boundary.GetParamCondition();
        if (ref_max > Param_condition.GetM())
          Param_condition.Resize(ref_max+1);
        
        if (parameters.GetM() > pos)
          {
            // we store these values in parameters_condition
            Vector<Complexe> parameters_condition(parameters.GetM()- pos-1);
            for (int i = pos+1; i < parameters.GetM(); i++)
              parameters_condition(i - pos-1)
                = to_num<Complexe>(parameters(i));
	    
	    if (parameters_condition.GetM() <= 4)
	      {
		cout << "In SetInputData of VarTransmission_Helm" << endl;
		cout << "ReferenceTransmission needs more parameters, for instance :" << endl;
		cout << "ReferenceTransmission = ref PARAM A B C Decalage delta" << endl;
		cout << "Current parameters are : " << endl << parameters_condition << endl;
		abort();
	      }

            // dans l'ordre, A, B, C, Decalage, delta
            this->A = parameters_condition(0) ;
            this->B = parameters_condition(1) ;
            this->C = parameters_condition(2) ;
            this->Decalage = realpart(parameters_condition(3)) ;
            this->delta = realpart(parameters_condition(4));
          }
        else
          {
            abort();
          }
        
        this->transmission_references = ref_;
      }
  }
  
  
  //! initialization before adding transmission terms
  template<class Complexe, class Dimension>
  void VarTransmission_Helm<Complexe, Dimension>::InitTransmission()
  {
    if (this->transmission_references.GetM() <= 0)
      return;

    var_boundary.ProcRowNumbers_Impedance.Reallocate(var_problem.GetNbDof());
    var_boundary.ProcColumnNumbers_Impedance.Reallocate(var_problem.GetNbDof());
    var_boundary.ProcRowNumbers_Impedance.Fill(0);
    var_boundary.ProcColumnNumbers_Impedance.Fill(0);
    
    var_boundary.NewRowNumbers_Impedance.Reallocate(var_problem.GetNbDof());
    var_boundary.NewColumnNumbers_Impedance.Reallocate(var_problem.GetNbDof());
    var_boundary.NewRowNumbers_Impedance.Fill(-1);
    var_boundary.NewColumnNumbers_Impedance.Fill(-1);
    for (int i = 0; i < this->DdlVolPlus.GetM(); i++)
      {
        var_boundary.NewRowNumbers_Impedance(this->DdlVolPlus(i)) = this->DdlVolMinus(i);
        var_boundary.NewColumnNumbers_Impedance(this->DdlVolPlus(i)) = this->DdlVolMinus(i);
	
        var_boundary.NewRowNumbers_Impedance(this->DdlVolMinus(i)) = this->DdlVolPlus(i);
        var_boundary.NewColumnNumbers_Impedance(this->DdlVolMinus(i)) = this->DdlVolPlus(i);
      }
    
    int nb_dof_surf = this->DdlVolPlus.GetM();
    if ((this->OrdreTransmission == 2) && !ModeleCentre && !CasDirichlet) 
      var_boundary.ResizeNbDof(var_problem.GetNbDof() + 3 * nb_dof_surf);
    
    // delta = 2 pi / N
    this->delta = this->R0 * 2.0 *pi_wp / this->nb_sections;
  }

  
  //! adding terms due to transmission condition
  template<class Complexe, class Dimension> 
  void VarTransmission_Helm<Complexe, Dimension>::
  AddTransmissionTerms(const Complex_wp& alpha, const GlobalGenericMatrix<Complex_wp>& nat_mat,
                       VirtualMatrix<Complex_wp>& mat_sp, int offset_row, int offset_col)
  {
    if (this->transmission_references.GetM() <= 0)
      return;
    
    this->SetModifiedColNumbers(0);
    this->SetModifiedRowNumbers(0);
    
    IVect ref_plus(var_problem.mesh.GetNbReferences()+1);
    IVect ref_minus(var_problem.mesh.GetNbReferences()+1);
    ref_plus.Fill(0); ref_minus.Fill(0);
    ref_plus(this->RefGammaPlus) = 1; ref_minus(this->RefGammaMoins) = 1;

    if (ModeleTest)
      {
	Complexe m_omega2; var_problem.GetMomega2(m_omega2);
	ImpedanceFunction_Base<Complex_wp, Dimension> fimpedance(var_problem);
        
	Complex_wp coef_stiff = nat_mat.GetCoefStiffness();
	Complex_wp coef_mass = nat_mat.GetCoefMass();
	Complex_wp coefA = 1.0/(this->delta*this->A)*coef_stiff;
	Complex_wp coefB0 = -0.25*this->delta*this->B*coef_mass*m_omega2;
	Complex_wp coefB2 = -0.25*this->delta*this->C*coef_stiff;
	
	fimpedance.SetCoefficient(coefA + coefB0, coefB2);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, false, false, var_problem);
	
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, true, var_problem);
	
	fimpedance.SetCoefficient(-coefA + coefB0, coefB2);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, false, true, var_problem);
	
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, false, var_problem);
	
        return;
      }
    
    Complex_wp zero(0, 0);
    if (CasDirichlet)
      {
         this->AalphaDelta = this->delta*(this->A + 2 *this->Decalage/this->R0) ;
         this->B1alphaDelta = this->delta*(this->B + this->Decalage/ (2* this->R0)) ;
         DISP(this->AalphaDelta) ;
         DISP(this->B1alphaDelta) ;
         DISP(this->delta) ;
         DISP(this->R0) ;
         
	 cout << " termes sur le bord Gamma plus" << endl ;
	 ImpedanceFunction_Base<Complex_wp, Dimension> fimpedance(var_problem);
             
	 // (1/(delta R+) (1/A0 + 1/4 B0)) int_gamma plus u+ v+ 
	 fimpedance.SetCoefficient(mu0 * ( 1.0/this->AalphaDelta + 0.25/this->B1alphaDelta)
				   / this->R0, zero);

	 cout << " terme 1 " << endl ;
	 //DISP(fimpedance.coef_phi) ;
	 var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						 offset_row, offset_col, fimpedance, false, false, var_problem);
	 
	 // (1/(delta R+) (-1/A0 + 1/4B0)) int_gamma plus u- v+ 
	 fimpedance.SetCoefficient(mu0 * ( -1.0/this->AalphaDelta + 0.25/this->B1alphaDelta)
				   / this->R0, zero);
	 
	 cout << " terme 2 " << endl ;
	 //DISP(fimpedance.coef_phi) ;
	 var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						 offset_row, offset_col, fimpedance, true, false, var_problem);       
      
	 cout << " termes sur le bord Gamma moins" << endl ;
	 // (1/(delta R-) (1/A0 + 1/4B0)) int_gamma plus u- v- 
	 fimpedance.SetCoefficient(mu0 * ( 1.0/this->AalphaDelta + 0.25/this->B1alphaDelta)
				   / this->R0, zero);
	 
	 cout << " terme 1 bis " << endl ;
	 //DISP(fimpedance.coef_phi) ;
	 var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						 offset_row, offset_col, fimpedance, false, false, var_problem);
	 //DISP(fimpedance.coef_phi) ;
	 
	 fimpedance.SetCoefficient(mu0 * ( -1.0/this->AalphaDelta + 0.25/this->B1alphaDelta)
				   / this->R0 , zero);
	 
	 cout << " terme 2 bis" << endl ;
	 //DISP(fimpedance.coef_phi) ;
	 var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						 offset_row, offset_col, fimpedance, true, false, var_problem);
             
         return;
      }
    
    Complex_wp B(0, 0), C(0, 0), D(0, 0), E(0, 0), F(0, 0), Alpha(0, 0);
    this->AalphaDelta = this->delta*(this->A + 2 *this->Decalage/this->R0) ;
    this->B1alphaDelta = this->delta*(this->B - (2 *this->Decalage * this->R0
                                                 *var_problem.GetSquareOmega()* rho0 ) /mu0);
    
    this->B2alphaDelta = this->delta*(this->C - 2 *this->Decalage/this->R0)  ;
    if (this->OrdreTransmission == 2)
      {
        F = this->AalphaDelta / this->delta ;
        //B = var.Param_condition(ref)(0)/(  this->delta * this->R0)  - 1./ this->R0 
        // * (var.Param_condition(ref)(2)/ this->delta + 2. * Alpha/ this->R0)  ;
        
        //C = (var.GetSquareOmega() * rho0 * this->R0 *var.Param_condition(ref)(0)
        //   / ( this->delta *  mu0))
        //  - (var.Param_condition(ref)(1)/ this->delta + 2. * Alpha * this->R0 
        //  *var.GetSquareOmega() * rho0/ mu0)/this->R0 ;
        
        D = this->B1alphaDelta/ this->delta; //D = - Boalpha * omega2  (notations these)
        E = - this->B2alphaDelta/ this->delta; // E = B2alpha (notationThese)
        B = E *F*0.5 ;
        C= - F * D *0.5; 
        DISP(F);
        DISP(B);
        DISP(C);
        DISP(D);
        DISP(E);
      }
    
    if ((this->OrdreTransmission == 1) && (!ModeleCentre) )
      {
	cout << "on traite un bord de type  gamma-" << endl;
	ImpedanceFunction_Base<Complex_wp, Dimension> fimpedance(var_problem);
	
	// on calcul :(( B1 + 1/(A)) * \mu_inf  )
	// (int_gamma^-  u- v- )- (B2*\mu_inf) (int_gamma^- du\dtheta dv\dtheta)   
	cout << " terme 1 sur gamma moins " << endl ;
	fimpedance.SetCoefficient(mu0 * ( 1.0/this->AalphaDelta + this->B1alphaDelta)
				  / this->R0, -mu0 * this->B2alphaDelta * this->R0);
	cout << "" << endl ;
            
	//DISP(fimpedance.coef_phi) ; DISP(fimpedance.coef_grad) ;
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, false, false, var_problem);
	
	//  on calcul : -(mu_inf /A)  (int_gamma^- (u+ v-))
	cout << " terme 2 sur gamma moins :changement de colonne " << endl ;
	fimpedance.SetCoefficient(-mu0/ ( this->R0  * this->AalphaDelta), zero);
	//DISP(fimpedance.coef_phi) ; DISP(fimpedance.coef_grad) ;
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, false, var_problem);
	cout << "" << endl ;		  
      }
    else if ( ModeleCentre )
      {
	cout << " terme 1 sur gamma moins " << endl ;
	ImpedanceFunction_Base<Complex_wp, Dimension> fimpedance(var_problem);
	
	fimpedance.SetCoefficient(mu0 * ( 1.0/this->AalphaDelta
                                              + this->B1alphaDelta * 0.25) / this->R0,
				  -mu0 * this->B2alphaDelta * 0.25 * this->R0 );
	
	cout << "" << endl ;
	//DISP(fimpedance.coef_phi) ; DISP(fimpedance.coef_grad) ;
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, false, false, var_problem);
            
	fimpedance.SetCoefficient(mu0 * (  this->B1alphaDelta * 0.25) / this->R0,
				  -mu0 * this->B2alphaDelta * 0.25 * this->R0 );
	
	cout << "" << endl ;
	//DISP(fimpedance.coef_phi) ; DISP(fimpedance.coef_grad) ;
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, false, var_problem);
	
	//  on calcul : -(mu_inf /A)  (int_gamma^- (u+ v-))
	cout << " terme 2 sur gamma moins :changement de colonne " << endl ;
	fimpedance.SetCoefficient(-mu0/ ( this->R0  * this->AalphaDelta), zero);
	//DISP(fimpedance.coef_phi) ; DISP(fimpedance.coef_grad) ;
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, false, var_problem);
	cout << "" << endl ;
      }
    else if ((this->OrdreTransmission == 2) & (!ModeleCentre))
      {
	
	// verification realisee pour les lignes et colonnes 
	// (en utilisant boundary condition harmonic) 
	cout << "ordre 2"<< endl; 
	//bord de type Gamma moins
	// termes 1 et 12
	cout <<  "termes 1 et 12" << endl ;
	ImpedanceFunction_Base<Complex_wp, Dimension> fimpedanceMoins1(var_problem);

	fimpedanceMoins1.SetCoefficient(mu0* this->delta * D / this->R0,
					this->delta * mu0 * E * this->R0);
	
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedanceMoins1, false, false, var_problem);
	//DISP(fimpedanceMoins1.coef_transmission) ;
	//DISP(fimpedanceMoins1.coef_gradtransmission) ;
        
	// terme 5
	cout <<  "terme 5" << endl ;
	fimpedanceMoins1.SetCoefficient(-  mu0 / (this->R0 * F *  this->delta ), zero);
	this->SetModifiedColNumbers(2);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedanceMoins1, true, false, var_problem);
	// DISP(fimpedanceMoins1.coef_transmission) ;
	// DISP(fimpedanceMoins1.coef_gradtransmission) ;
	
	// terme 3
	cout <<  "terme 3" << endl ;
	fimpedanceMoins1.SetCoefficient(mu0  / (this->R0 * F *  this->delta ), zero);
	this->SetModifiedColNumbers(3);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedanceMoins1, true, false, var_problem);
	//DISP(fimpedanceMoins1.coef_transmission) ;
	// DISP(fimpedanceMoins1.coef_gradtransmission) ;
	
	//termes 8 et 10
	cout <<  "terme 8" << endl ;
	fimpedanceMoins1.SetCoefficient(mu0 /(this->R0*( 1. + C * this->delta * this->delta)),
					mu0 *  this->R0 *  (this->delta *  this->delta * B )
					/ ( 1. + Alpha * this->delta * this->delta * C) );
	
	this->SetModifiedColNumbers(1);
	this->SetModifiedRowNumbers(1);
	
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedanceMoins1, true, true, var_problem); // terme 8
	cout <<  "terme 10" << endl ;
	
	this->SetModifiedColNumbers(3);
	this->SetModifiedRowNumbers(3);
	
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedanceMoins1, true, true, var_problem);  //terme 10
         
	//termes 9  et 11
	cout <<  "terme 9" << endl ;
	fimpedanceMoins1.SetCoefficient(-mu0/ this->R0, zero);
	
	this->SetModifiedColNumbers(1);
	this->SetModifiedRowNumbers(1);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedanceMoins1, false, true, var_problem); //terme 9
	
	cout <<  "terme 11" << endl ;
	this->SetModifiedRowNumbers(3);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_minus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedanceMoins1, true, true, var_problem); //terme 11
	// DISP(fimpedanceMoins1.coef_transmission) ;
	// DISP(fimpedanceMoins1.coef_gradtransmission) ;
	cout << "________________________________________________________" << endl;
	cout << ""<< endl;
      }
    
    if ((this->OrdreTransmission == 1) & (!ModeleCentre) )
      {
	cout << "on traite un  bord de type  gamma+" << endl;
        
	//on calcul : (mu_inf /A)  (int_gamma^+ (u+ v+))
	cout << " terme 3 sur gamma plus " << endl ;
	ImpedanceFunction_Base<Complex_wp, Dimension> fimpedance(var_problem);
            
	fimpedance.SetCoefficient(mu0/ ( this->R0  * this->AalphaDelta), zero);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, false, false, var_problem);
	//DISP(fimpedance.coef_phi) ;
	cout << "  " << endl ;
	
	//  on calcul : -(mu_inf /A)  (int_gamma^+ (u- v+))
	cout << " terme 4 sur gamma plus " << endl ;
	fimpedance.SetCoefficient(-mu0/ ( this->R0  * this->AalphaDelta), zero);
	//DISP(fimpedance.coef_phi) ; 
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, false, var_problem);
	cout << "  " << endl ;
      }
    else if (ModeleCentre)
      {
	ImpedanceFunction_Base<Complex_wp, Dimension> fimpedance(var_problem);
	fimpedance.SetCoefficient(mu0 * (  this->B1alphaDelta * 0.25) / this->R0,
				  -mu0 * this->B2alphaDelta * 0.25 * this->R0 );
	cout << "" << endl ;
	//DISP(fimpedance.coef_phi) ; DISP(fimpedance.coef_grad) ;
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, false, false, var_problem);
	
	fimpedance.SetCoefficient(mu0 * (  this->B1alphaDelta * 0.25) / this->R0,
				  -mu0 * this->B2alphaDelta * 0.25 * this->R0 );
	cout << "" << endl ;
	//DISP(fimpedance.coef_phi) ; DISP(fimpedance.coef_grad) ;
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, false, var_problem);

	
	//on calcul : (mu_inf /A)  (int_gamma^+ (u+ v+))
	cout << " terme 3 sur gamma plus " << endl ;
	
	fimpedance.SetCoefficient(mu0/ ( this->R0  * this->AalphaDelta), zero);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, false, false, var_problem);
	//DISP(fimpedance.coef_phi) ;
	cout << "  " << endl ;
	//  on calcul : -(mu_inf /A)  (int_gamma^+ (u- v+))
	cout << " terme 4 sur gamma plus " << endl ;
	fimpedance.SetCoefficient(-mu0/ ( this->R0  * this->AalphaDelta), zero);
	//DISP(fimpedance.coef_phi) ; 
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedance, true, false, var_problem);
	cout << "  " << endl ;
      }
    else if ( (this->OrdreTransmission == 2)  & (!ModeleCentre))
      {
	cout << "ordre 2"<< endl ;
	// verification realisee
	ImpedanceFunction_Base<Complex_wp, Dimension> fimpedancePlus1(var_problem);
	// terme 4 
	cout <<  "terme 4" << endl ;
	fimpedancePlus1.SetCoefficient(-mu0/( F * this->R0 * this->delta), zero);
	this->SetModifiedColNumbers(1);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedancePlus1, true, false, var_problem);
	// DISP(fimpedancePlus1.coef_transmission) ;
	// DISP(fimpedancePlus1.coef_gradtransmission) ;
	
	//terme 2 
	cout <<  "terme 2"  << endl ;
	fimpedancePlus1.SetCoefficient(mu0/( F * this->R0 * this->delta), zero);
	
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedancePlus1, false, false, var_problem);
	//DISP(fimpedancePlus1.coef_transmission) ;
	//DISP(fimpedancePlus1.coef_gradtransmission) ;
	
	//terme 6
	cout <<  "terme 6" << endl ;
	fimpedancePlus1.SetCoefficient(mu0/( this->R0*( 1. + C* this->delta * this->delta)),
				       mu0* this->R0 *  this->delta *  this->delta * B
				       / ( 1. + this->delta * this->delta  * C));
	//fimpedancePlus1.coef_transmission =var.->mu0/ this->TabRmoinsRplus(indice)(1);
	//fimpedancePlus1.coef_gradtransmission = 0;
	this->SetModifiedColNumbers(2);
	this->SetModifiedRowNumbers(2);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedancePlus1, true, true, var_problem);
	// DISP(fimpedancePlus1.coef_transmission) ;
	// DISP(fimpedancePlus1.coef_gradtransmission);
	
	//terme 7
	cout <<  "terme 7" << endl ;
	fimpedancePlus1.SetCoefficient(-mu0/ this->R0, zero);
	var_boundary.AddMatrixImpedanceBoundary(alpha, ref_plus, 1, nat_mat, mat_sp,
						offset_row, offset_col, fimpedancePlus1, false, true, var_problem);
	//DISP(fimpedancePlus1.coef_transmission) ;
	// DISP(fimpedancePlus1.coef_gradtransmission) ;
	
	
	cout << "________________________________________________________" << endl;
	cout << ""<< endl;
      }
  }
  
}// end namespace

#define MONTJOIE_FILE_TRANSMISSION_MODEL_HELMHOLTZ_CXX
#endif

  
