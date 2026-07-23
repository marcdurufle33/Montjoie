#ifndef MONTJOIE_FILE_EQUATION_STRING_HXX

namespace Montjoie
{
  
  //! Liste d'equations disponibles pour modeliser le deplacement d'une corde
  /*!
    L'equation modele est decrite dans la these de Juliette Chabassier :
    M d^2 q/dt^2 + d/dt( R q - d/dx(H dq/dx))
      - d/dx( A dq/dx + B q + grad U (dq/dx)) + C q + B^T dq/dx = S(x)
    
    q regroupe les differentes inconnues (u, v, phi pour le modele le plus complet)
    M : matrice de masse, ApplyM est la methode implementant le produit M x
    R : matrice d'amortissement constant, ApplyR est la methode implementant le produit R x
    H : matrice d'amortissement quadratique,
        ApplyAdamping est la methode implementant le produit H x
    A : matrice de rigidite, Le produit A x est divise en deux A x = A^grad x + A^stiff x
        A^grad x est implemente dans ApplyAGrad, A^stiff x dans ApplyAStiffGrad
	des theta differents peuvent ensuite etre appliques sur chaque operateur
    B : matrice de "gradient", ApplyB implemente le produit B x, ApplyD le produit B^T x
    U : non-linear term, ComputeU computes U(x), ComputeDeltakU computes grad U(x)
        in a conservative manner
        ComputeDerivativeDeltakU computes the derivative of ComputeDeltakU with respect to x
    C : matrice, ApplyC est la methode implementant le produit C x
    
    Sur toutes les equations on a des termes d'amortissement similaires :
    2 rho A b1 \frac{\partial u}{\partial t} + 2 T0 b3 \frac{\partial^3 u}{\partial t \partial x^2}
    b1 est l'amortissement constant (note amo), et b3 l'amortissement quadratique (note amoB2)
    
    Par consequent pour chaque documentation de classe,
    on ne mettra que la version sans amortissement
   */
  class GenericEquationString
  {
  public :
    enum {WAVE_LIN, WAVE_LIN_TIMO, WAVE_NL, WAVE_NL_TIMO, WAVE_PIPO, WAVE_NL_2T,WAVE_NL_2T_TIMO,WAVE_LIN_2T};    

    //! produit scalaire pour separer energie u et energie v
    template<class T, int m, int n>
    static void DotProdColUV(const TinyMatrix<T, General, m, n>& A,
			     const TinyMatrix<T, General, m, n>& B,
			     TinyVector<T, n>& res_u, TinyVector<T, n>& res_v)
    {
      // cas general, on remplit juste res_u
      DotProdCol(A, B, res_u);
    }
    
    
    //! produit scalaire pour separer energie u et energie v
    template<class T, int m>
    static void DotProdColUV(const TinyMatrix<T, General, 1, m>& A,
			     const TinyMatrix<T, General, 1, m>& B,
			     TinyVector<T, m>& res_u, TinyVector<T, m>& res_v)
    {
      for (int p = 0; p < m; p++)
	res_u(p) += A(0, p)*B(0, p);
    }
    
    
    //! produit scalaire pour separer energie u et energie v
    template<class T, int m>
    static void DotProdColUV(const TinyMatrix<T, General, 2, m>& A,
			     const TinyMatrix<T, General, 2, m>& B,
			     TinyVector<T, m>& res_u, TinyVector<T, m>& res_v)
    {
      for (int p = 0; p < m; p++)
	{
	  res_u(p) += A(0, p)*B(0, p);
	  res_v(p) += A(1, p)*B(1, p);
	}
    }
    
    
    //! produit scalaire pour separer energie u et energie v
    template<class T, int m>
    static void DotProdColUV(const TinyMatrix<T, General, 3, m>& A,
			     const TinyMatrix<T, General, 3, m>& B,
			     TinyVector<T, m>& res_u, TinyVector<T, m>& res_v)
    {
      for (int p = 0; p < m; p++)
	{
	  res_u(p) += A(0, p)*B(0, p) + A(2, p)*B(2, p);
	  res_v(p) += A(1, p)*B(1, p);
	}
    }
    
  };
  

  // ***********************************
  // Equation de corde linéaire scalaire
  // ***********************************

  
  //! Equation d'onde lineaire (equation d'Alembert)
  /*!
    Le deplacement transverse u est solution de l'equation des ondes suivantes : 
    rho A \frac{\partial^2 u}{\partial t^2} - T0 \frac{\partial^2 u}{\partial x^2} = 0
   */
  class WaveEquationString : public GenericEquationString
  {
  public :
    
    enum {dimension = 1}; // connu a la compil, super pour templater
    enum {type_equation = WAVE_LIN};
    static const bool linear = true;
    
    TinyVector<bool, dimension> Dirichlet_Agraffe;    
    TinyVector<bool, dimension> Dirichlet_Bridge;    
    bool only_one_LM;
    
    //! default constructor
    WaveEquationString()
    {
      Dirichlet_Agraffe(0) = true;
      Dirichlet_Bridge(0) = true;
      only_one_LM = true;
    }

    //! application de la matrice de masse, vect_out = M vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyM(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out=piano_string.rho*piano_string.A*vect_in;
    } 
    
    //! application de la matrice d'amortissement, vect_out = R vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyR(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out = 2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*vect_in;
    } 
    
    //! application de l'amortissement quadratique, vect_out = H vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAdamping(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
		  TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
		  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out = 2.0*piano_string.T0*piano_string.amoB2(0)*vect_in;
    }
    
    //! application de la rigidite non-raide, vect_out = A^grad vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	       TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	       PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out=piano_string.T0*vect_in;
    } 
    
    //! application de la rigidite raide, vect_out = A^stiff vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAStiffGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
		    TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
		    PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
    } 
    
    //! application de la matrice B, vect_out = B vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyB(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
    } 
    
    //! application de la matrice C, vect_out = C vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyC(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
    }
    
    //! application de la transposee de la matrice B, vect_out = B^T vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyD(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
    } 
    
    //! Calcul de vect_out = U(Q)
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeU( TinyMatrix<Real_wp, General, dimension, nb_quad> & Q,
	      TinyVector<Real_wp, nb_quad> & vect_out,
              PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();      
    }
        
    //! calcul de vect_out = grad(U) (q) a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus, TinyVector<Real_wp, nb_quad> & Qminus,
		    TinyMatrix<Real_wp, General, dimension-1, nb_quad>& Qsigma,
		    TinyVector<Real_wp, nb_quad>& vect_out,
                    PianoString<nb_base, nb_quad, TypeEquation> & piano_string, int k)
    {
      vect_out.Fill(0.0);
    }
    
    //! calcul de vect_out = d/dq( grad(U)(q))  a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDerivativeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,
                              TinyVector<Real_wp, nb_quad> & Qminus,
			      TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
			      TinyVector<Real_wp, nb_quad> & vect_out,
			      PianoString<nb_base, nb_quad, TypeEquation> & piano_string,
                              int k1, int k2)
    {
      vect_out.Zero();      
    }

    //! application of the scheme    
    template<int nb_base> static inline void
    ApplyMatrixEquation(TinyMatrix<TinyMatrix<Real_wp, General, nb_base, nb_base>,
                        General, dimension, dimension>& mat_Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& res_Un)
    {
      Mlt(mat_Un(0, 0), Un(0), res_Un(0));
    }

    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {
    }
    
  };
  
  
  // *******************************
  // Equation de corde linéaire raide
  // *******************************

  
  //! Equation de Timoshenko (inconnues u, phi)
  /*!
    rho A d^2 u/dt^2 - (A G K + T0) d^2 u/dx^2 + A G K dphi/dx = S(x, t)
    rho I d^2 phi/dt^2 - E I d^2 phi/dx^2 - A G K du/dx + A G K phi = 0
   */
  class WaveEquationStiffString : public GenericEquationString
  {
  public :
    enum{dimension = 2}; // connu a la compil, super pour templater
    enum { type_equation = WAVE_LIN_TIMO};
    static const bool linear = true;
    
    TinyVector<bool,dimension> Dirichlet_Agraffe;
    TinyVector<bool,dimension> Dirichlet_Bridge;
    
    // true si on a une equation en (u, phi) seulement ou u seulement
    bool only_one_LM;
    
    //! default constructor
    WaveEquationStiffString()
    {
      Dirichlet_Agraffe(0) = true;
      Dirichlet_Agraffe(1) = false;
      Dirichlet_Bridge(0) = true;
      Dirichlet_Bridge(1) = false;
      
      only_one_LM = true;
    }
        
    //! application de la matrice de masse, vect_out = M vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyM(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.rho*piano_string.I*ligne_temp, 1, vect_out);
    } 
    
    //! application de la matrice d'amortissement, vect_out = R vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyR(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.I*piano_string.amo(2)*ligne_temp, 1, vect_out);
    }
    
    //! application de l'amortissement quadratique, vect_out = H vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyAdamping(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                  TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                  const GenericPb& piano_string)
    {
      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.T0*piano_string.amoB2(0)*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.I*piano_string.amoB2(2)*ligne_temp, 1, vect_out);
    }
    
    //! application de la rigidite non-raide, vect_out = A^grad vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyAGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
               TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
               const GenericPb& piano_string)
    {      
      vect_out.Zero();
      
      TinyVector<Real_wp,nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.T0)*ligne_temp, 0, vect_out);
    } 
    
    //! application de la rigidite raide, vect_out = A^stiff vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyAStiffGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                    TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                    const GenericPb& piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp,nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.A*piano_string.G*piano_string.k_prime)*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow((piano_string.E*piano_string.I)*ligne_temp, 1, vect_out);
    } 
    
    //! application de la matrice B, vect_out = B vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyB(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      vect_out.Zero();
      TinyVector<Real_wp,nb_quad> ligne_temp;
      GetRow(vect_in, 1, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 0, vect_out);
      
      //  GetRow(vect_in, 0, ligne_temp);
      //  SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 1, vect_out);
    } 
    
    //! application de la matrice C, vect_out = C vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyC(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      vect_out.Fill(0.0);
      TinyVector<Real_wp,nb_quad> ligne_temp;
      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 1, vect_out);
      
    } 
        
    //! application de la transposee de la matrice B, vect_out = B^T vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyD(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      vect_out.Zero();
      TinyVector<Real_wp,nb_quad> ligne_temp;
      // GetRow(vect_in, 1, ligne_temp);
      //  SetRow(piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 1, vect_out);
    } 
    
    //! Calcul de vect_out = U(Q)
    template<int nb_quad, class GenericPb> static inline void
    ComputeU( TinyMatrix<Real_wp,General,dimension, nb_quad> & Q,
              TinyVector<Real_wp,nb_quad> & vect_out, GenericPb & piano_string)
    {
      vect_out.Zero();
      
    }
    
    //! calcul de vect_out = grad(U) (q) a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_quad, class GenericPb> static inline void
    ComputeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus, TinyVector<Real_wp, nb_quad> & Qminus,
                    TinyMatrix<Real_wp,General,dimension-1,nb_quad> & Qsigma,
                    TinyVector<Real_wp,nb_quad> & vect_out, GenericPb & piano_string, int k)
    {
      vect_out.Fill(0.0);
    }

    //! calcul de vect_out = d/dq( grad(U)(q))  a partir de Q^n-1, Q^n+1 et Q^n    
    template<int nb_quad, class GenericPb> static inline void
    ComputeDerivativeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,
                              TinyVector<Real_wp, nb_quad> & Qminus,
                              TinyMatrix<Real_wp,General,dimension-1,nb_quad> & Qsigma,
                              TinyVector<Real_wp, nb_quad> & vect_out,
                              GenericPb & piano_string, int k1, int k2)
    {
      vect_out.Zero();      
    }

    //! application of the scheme    
    template<int nb_base> static inline void
    ApplyMatrixEquation(TinyMatrix<TinyMatrix<Real_wp, General, nb_base, nb_base>,
                        General, dimension, dimension>& mat_Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& res_Un)
    {
      Mlt(mat_Un(0, 0), Un(0), res_Un(0));
      MltAdd(mat_Un(0, 1), Un(1), res_Un(0));        
      Mlt(mat_Un(1, 0), Un(0), res_Un(1));
      MltAdd(mat_Un(1, 1), Un(1), res_Un(1));
    }

    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {
    }
    
  };
  
  //! Classes pour deroulement de boucles
  template <int g> 
  class WaveEquationStiffNLStringLoop
  {
  public:
    template<int nb_base, int nb_quad, class TypeEquation, int dimensionM1>
    static inline void ComputeDeltakUk0( TinyVector<Real_wp, nb_quad> & Qplus,
				  TinyVector<Real_wp, nb_quad> & Qminus,
				  TinyMatrix<Real_wp, General, dimensionM1, nb_quad> & Qsigma,
				  TinyVector<Real_wp, nb_quad> & vect_out,
				  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      
      Real_wp s1 = sqrt(Qplus(g) *Qplus(g)  + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
      Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
      vect_out(g) = (piano_string.E*piano_string.A - piano_string.T0)
        * ( (0.5 - 1.0/(s1+s2))* (Qplus(g)+Qminus(g)) );
      /*
	  Real_wp s1 = sqrt(Qplus(g) *Qplus(g)  + (Qsigma(0,g))*(Qsigma(0,g)));
	  Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (Qsigma(0,g))*(Qsigma(0,g)));
	  vect_out(g) = (piano_string.E*piano_string.A - piano_string.T0)
		* ( (0.5 - 1.0/(s1+s2))* (Qplus(g)+Qminus(g)) );
	  */	
      WaveEquationStiffNLStringLoop<g-1>::ComputeDeltakUk0(Qplus, Qminus, Qsigma,
                                                           vect_out, piano_string);
    }
    
    template<int nb_base, int nb_quad, class TypeEquation, int dimensionM1>
    static inline void ComputeDeltakUk1( TinyVector<Real_wp, nb_quad> & Qplus,
                                         TinyVector<Real_wp, nb_quad> & Qminus,
				  TinyMatrix<Real_wp, General, dimensionM1, nb_quad> & Qsigma,
				  TinyVector<Real_wp, nb_quad> & vect_out,
				  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      
      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qplus (g))*(1.0+Qplus (g)) );
      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qminus(g))*(1.0+Qminus(g)) );
      vect_out(g) = (piano_string.E*piano_string.A - piano_string.T0)
        * (1.0 - (2.0+Qplus(g)+Qminus(g))/(s1+s2));
      /*
      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qplus (g))*(Qplus (g)) );
	  Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qminus(g))*(Qminus(g)) );
	  vect_out(g) = (piano_string.E*piano_string.A - piano_string.T0)
		* ( (0.5 - 1.0/(s1+s2))* (Qplus(g)+Qminus(g)) );
	  */
      WaveEquationStiffNLStringLoop<g-1>::ComputeDeltakUk1(Qplus, Qminus, Qsigma,
                                                           vect_out, piano_string);
    }
	
    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {
      Real_wp norm = 1.0 / sqrt(dTm(0, g)*dTm(0, g) + square(1.0 + dTm(1, g)));
      nablaU(0, g) = coefNL*dTm(0, g)*(1.0 - norm);
      nablaU(1, g) = coefNL*(1.0 - (1.0+dTm(1, g))*norm);
      WaveEquationStiffNLStringLoop<g-1>::ApplyNonLinearTerm(coefNL, dTm, nablaU);
    }
    
  };
  
  
  //! specialisation pour terminer la boucle
  template <> 
  class WaveEquationStiffNLStringLoop<-1>
  {
  public:
    template<int nb_base, int nb_quad, class TypeEquation, int dimensionM1>
    static inline void ComputeDeltakUk0(TinyVector<Real_wp, nb_quad> & Qplus,
                                        TinyVector<Real_wp, nb_quad> & Qminus,
                                        TinyMatrix<Real_wp, General, dimensionM1, nb_quad>& Qsigma,
                                        TinyVector<Real_wp, nb_quad> & vect_out,
                                        PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    { }
    
    template<int nb_base, int nb_quad, class TypeEquation, int dimensionM1>
    static inline void ComputeDeltakUk1( TinyVector<Real_wp, nb_quad> & Qplus,
                                         TinyVector<Real_wp, nb_quad> & Qminus,
                                         TinyMatrix<Real_wp, General,
                                         dimensionM1, nb_quad> & Qsigma,
                                         TinyVector<Real_wp,nb_quad> & vect_out,
                                         PianoString<nb_base, nb_quad, TypeEquation>& piano_string)
    { }	

    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {}

  };


  template <int g> 
  class WaveEquationStiffNL2TStringLoop
  {
  public:
    template<int nb_base, int nb_quad, class TypeEquation, int dimensionM1>
    static inline void ComputeDeltakUk0( TinyVector<Real_wp, nb_quad> & Qplus,
				  TinyVector<Real_wp, nb_quad> & Qminus,
				  TinyMatrix<Real_wp, General, dimensionM1, nb_quad> & Qsigma,
				  TinyVector<Real_wp, nb_quad> & vect_out,
				  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {

	      Real_wp s1 = sqrt(Qplus (g)*Qplus (g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
	      Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
	      vect_out(g) = (piano_string.E*piano_string.A - piano_string.T0)
                * ( (0.5 - 1.0/(s1+s2))* (Qplus(g)+Qminus(g)) );
	  	
      WaveEquationStiffNL2TStringLoop<g-1>::ComputeDeltakUk0(Qplus, Qminus, Qsigma,
                                                           vect_out, piano_string);
    }
    
    template<int nb_base, int nb_quad, class TypeEquation, int dimensionM1>
    static inline void ComputeDeltakUk1( TinyVector<Real_wp, nb_quad> & Qplus,
                                         TinyVector<Real_wp, nb_quad> & Qminus,
				  TinyMatrix<Real_wp, General, dimensionM1, nb_quad> & Qsigma,
				  TinyVector<Real_wp, nb_quad> & vect_out,
				  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qplus (g))*(Qplus (g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qminus(g))*(Qminus(g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
	      vect_out(g) = (piano_string.E*piano_string.A - piano_string.T0)
                * ( (0.5 - 1.0/(s1+s2))* (Qplus(g)+Qminus(g)) );
	  
      WaveEquationStiffNL2TStringLoop<g-1>::ComputeDeltakUk1(Qplus, Qminus, Qsigma,
                                                           vect_out, piano_string);
    }
    
	template<int nb_base, int nb_quad, class TypeEquation, int dimensionM1>
    static inline void ComputeDeltakUk2( TinyVector<Real_wp, nb_quad> & Qplus,
                                         TinyVector<Real_wp, nb_quad> & Qminus,
				  TinyMatrix<Real_wp, General, dimensionM1, nb_quad> & Qsigma,
				  TinyVector<Real_wp, nb_quad> & vect_out,
				  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qplus (g))*(1+Qplus (g)));
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qminus(g))*(1+Qminus(g)));
	      vect_out(g) = (piano_string.E*piano_string.A - piano_string.T0)
                * (1.0 - (2.0+Qplus(g)+Qminus(g))/(s1+s2));
                
      WaveEquationStiffNL2TStringLoop<g-1>::ComputeDeltakUk2(Qplus, Qminus, Qsigma,
                                                           vect_out, piano_string);
    }
	
    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {

      
	  Real_wp norm = 1.0 / sqrt(dTm(0, g)*dTm(0, g) + dTm(1, g)*dTm(1, g) + square(1.0 + dTm(2, g)));
      nablaU(0, g) = coefNL*dTm(0, g)*(1.0 - norm);
	  nablaU(1, g) = coefNL*dTm(1, g)*(1.0 - norm);
      nablaU(2, g) = coefNL*(1.0 - (1.0+dTm(2, g))*norm);
      WaveEquationStiffNL2TStringLoop<g-1>::ApplyNonLinearTerm(coefNL, dTm, nablaU);
    }
    
  };
  
  
  //! specialisation pour terminer la boucle
  template <> 
  class WaveEquationStiffNL2TStringLoop<-1>
  {
  public:
    template<int nb_base, int nb_quad, class TypeEquation, int dimensionM1>
    static inline void ComputeDeltakUk0(TinyVector<Real_wp, nb_quad> & Qplus,
                                        TinyVector<Real_wp, nb_quad> & Qminus,
                                        TinyMatrix<Real_wp, General, dimensionM1, nb_quad>& Qsigma,
                                        TinyVector<Real_wp, nb_quad> & vect_out,
                                        PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    { }
    
    template<int nb_base, int nb_quad, class TypeEquation, int dimensionM1>
    static inline void ComputeDeltakUk1( TinyVector<Real_wp, nb_quad> & Qplus,
                                         TinyVector<Real_wp, nb_quad> & Qminus,
                                         TinyMatrix<Real_wp, General,
                                         dimensionM1, nb_quad> & Qsigma,
                                         TinyVector<Real_wp,nb_quad> & vect_out,
                                         PianoString<nb_base, nb_quad, TypeEquation>& piano_string)
    { }	
    template<int nb_base, int nb_quad, class TypeEquation, int dimensionM1>
    static inline void ComputeDeltakUk2( TinyVector<Real_wp, nb_quad> & Qplus,
                                         TinyVector<Real_wp, nb_quad> & Qminus,
                                         TinyMatrix<Real_wp, General,
                                         dimensionM1, nb_quad> & Qsigma,
                                         TinyVector<Real_wp,nb_quad> & vect_out,
                                         PianoString<nb_base, nb_quad, TypeEquation>& piano_string)
    { }	
    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {}

  };


  
  // *******************************
  // Equation de corde non linéaire // mais pas raide
  // *******************************
  
  
  //! Equation d'onde non-lineaire (inconnues u, v)
  /*!
    Le deplacement transverse u et le deplacement longitudinal v sont solutions de
    rho A d^2 u/dt^2 - T0 d^2 u/dx^2
         - d/dx( (E A - T0) du/dx (1 - 1 / sqrt( (du/dx)^2 + (1 + dv/dv)^2 ) ) ) = 0
    rho A d^2 v/dt^2 - E A d^2 v/dx^2
         - d/dx( (E A - T0) (1 - 1 / sqrt( (du/dx)^2 + (1 + dv/dv)^2 ) ) ) = 0
   */
  class WaveEquationNLString : public GenericEquationString 
  {
  public :
    enum{dimension = 2}; // connu a la compil, super pour templater
    enum{ type_equation = WAVE_NL};
    static const bool linear = false;
    TinyVector<bool, dimension> Dirichlet_Agraffe;
    TinyVector<bool, dimension> Dirichlet_Bridge;
	
    bool only_one_LM;
    
    WaveEquationNLString()
    {
      Dirichlet_Agraffe(0) = true;
      Dirichlet_Agraffe(1) = true;
      Dirichlet_Bridge(0) = true;
      Dirichlet_Bridge(1) = true;
      
      only_one_LM = false;
    }
    
    //! application de la matrice de masse, vect_out = M vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyM(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 1, vect_out);
    } 
    
    //! application de la matrice d'amortissement, vect_out = R vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyR(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*ligne_temp, 0, vect_out);
    
      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(1)*ligne_temp, 1, vect_out);
    } 
    
    //! application de l'amortissement quadratique, vect_out = H vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAdamping(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                  TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.T0*piano_string.amoB2(0)*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.A*piano_string.amoB2(1)*ligne_temp, 1, vect_out);
      
    }
    
    //! application de la rigidite non-raide, vect_out = A^grad vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	       TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	       PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp,nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.T0)*ligne_temp, 0, vect_out);
    } 

    //! application de la rigidite raide, vect_out = A^stiff vect_in	
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAStiffGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
		    TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
		    PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
	  
      GetRow(vect_in, 1, ligne_temp);
      SetRow((piano_string.E*piano_string.A)*ligne_temp, 1, vect_out);
    } 
    
    //! application de la matrice B, vect_out = B vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyB(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
    } 
    
    //! application de la matrice C, vect_out = C vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyC(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
    } 
	
    //! application de la transposee de la matrice B, vect_out = B^T vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyD(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
    } 
    
    //! Calcul de vect_out = U(Q)
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeU( TinyMatrix<Real_wp, General, dimension, nb_quad> & Q,
	      TinyVector<Real_wp, nb_quad> & vect_out,
              PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    { 
      Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
      for (int g = 0; g < nb_quad; g++)
	{
	  Real_wp racine = sqrt(Q(0,g) *Q(0,g)  + (1.0+Q(1,g))*(1.0+Q(1,g)));	
	  vect_out(g) = coef*(Q(0,g) *Q(0,g) *0.5 + (1.0 + Q(1,g)) - racine);
	  
	}  
    }
    
    //! calcul de vect_out = grad(U) (q) a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,TinyVector<Real_wp, nb_quad> & Qminus,
		    TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
		    TinyVector<Real_wp, nb_quad> & vect_out,
                    PianoString<nb_base, nb_quad, TypeEquation> & piano_string, int k)
    {
      vect_out.Fill(0.0);
      if (k==0)
	{ 
	  // ligne u
	  for (int g = 0; g < nb_quad; g++)
	    {
	      // A DEROULER
	      Real_wp s1 = sqrt(Qplus(g) *Qplus(g)  + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
	      Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
	      vect_out(g) = (piano_string.E*piano_string.A - piano_string.T0)
                * ( (0.5 - 1.0/(s1+s2))* (Qplus(g)+Qminus(g)) );
	    }
	}
      
      if (k==1)
	{
	  // ligne v
	  for (int g = 0; g < nb_quad; g++)
	    {
	    // A DEROULER
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qplus (g))*(1.0+Qplus (g)) );
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qminus(g))*(1.0+Qminus(g)) );
	      vect_out(g) = (piano_string.E*piano_string.A - piano_string.T0)
                * (1.0 - (2.0+Qplus(g)+Qminus(g))/(s1+s2));
	    }
	}
    }
    
    //! calcul de vect_out = d/dq( grad(U)(q))  a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDerivativeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,
                              TinyVector<Real_wp, nb_quad> & Qminus,
			      TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
			      TinyVector<Real_wp, nb_quad> & vect_out,
			      PianoString<nb_base, nb_quad, TypeEquation> & piano_string,
                              int k1, int k2)
    {
      vect_out.Zero();
      if (k1 == 0)
	{
	  if(k2 == 0)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
		  Real_wp s1 = sqrt(Qplus(g) *Qplus(g)  + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
		  Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( 0.5 - 1.0/(s1+s2)
                                        +  (Qplus(g)+Qminus(g))*Qplus(g)/(s1*(s1+s2)*(s1+s2)) );
		}
	    }
	  else if (k2 == 1)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
		  Real_wp s1 = sqrt(Qplus(g) *Qplus(g)  + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
		  Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (Qplus(g)+Qminus(g))*(1.0+Qsigma(0,g))
                                        *(1.0/s1+1.0/s2)/( (s1+s2)*(s1+s2) ) );
		}
	    }
	}      
      else if (k1 == 1)
	{
	  if (k2 == 0)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
		  Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qplus (g))*(1.0+Qplus (g)) );
		  Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qminus(g))*(1.0+Qminus(g)) );
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (2.0+Qplus(g)+Qminus(g))*Qsigma(0,g)
                                        *(1.0/s1 + 1.0/s2) )/( (s1+s2)*(s1+s2) );
		}
	    }
	  else if (k2 ==1)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
		  Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qplus (g))*(1.0+Qplus (g)) );
		  Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qminus(g))*(1.0+Qminus(g)) );
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef * (-1.0/(s1+s2) +  (2.0+Qplus(g)+Qminus(g))
                                        *(1.0+Qplus(g))/(s1*(s1+s2)*(s1+s2)) );
		}
	    }
	}      
    }	

    //! application of the scheme    
    template<int nb_base> static inline void
    ApplyMatrixEquation(TinyMatrix<TinyMatrix<Real_wp, General, nb_base, nb_base>,
                        General, dimension, dimension>& mat_Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& res_Un)
    {
      Mlt(mat_Un(0, 0), Un(0), res_Un(0));
      Mlt(mat_Un(1, 1), Un(1), res_Un(1));
    }
    
    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {
      WaveEquationStiffNLStringLoop<nb_quad-1>::ApplyNonLinearTerm(coefNL, dTm, nablaU);
    }

  };
  
  
  // *******************************
  // Equation de corde NL PIPO // du lineaire dans le U
  // *******************************
  
  
  // classe pipo
  class WavePipo : public GenericEquationString
  {
  public :
    enum{dimension = 2}; // connu a la compil, super pour templater
    enum{ type_equation = WAVE_PIPO};
    static const bool linear = false;
    TinyVector<bool, dimension> Dirichlet_Agraffe;
    TinyVector<bool, dimension> Dirichlet_Bridge;
    
    bool only_one_LM;
    
    WavePipo()
    {
      Dirichlet_Agraffe(0) = false;
      Dirichlet_Agraffe(1) = false;
      Dirichlet_Bridge(0) = false;
      Dirichlet_Bridge(1) = false;
      
      only_one_LM = true;
    }
    
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyM(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
    } 
    
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyR(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
    } 
    
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAdamping(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                  TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
    }
	
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	       TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	       PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
    } 
    
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAStiffGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
		    TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
		    PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
    } 
	
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyB(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
    } 
    
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyC(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
    } 
    
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeU( TinyMatrix<Real_wp, General, dimension, nb_quad> & Q,
	      TinyVector<Real_wp, nb_quad> & vect_out,
              PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    { 
      //impossible sauf si G = A
      vect_out . Zero();
    }
    
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus, TinyVector<Real_wp, nb_quad> & Qminus,
		    TinyMatrix<Real_wp, General, dimension-1, nb_quad> Qsigma,
                    TinyVector<Real_wp, nb_quad> & vect_out,
		    PianoString<nb_base, nb_quad, TypeEquation> & piano_string, int k)
    {
      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(Qsigma, 0, ligne_temp);
      
      if (k == 0)
	vect_out = 0.5*(piano_string.E)*(Qplus +Qminus) + piano_string.G*ligne_temp;
      else if (k == 1)
	vect_out = (piano_string.A)*ligne_temp + 0.5*piano_string.I*(Qplus+Qminus);      
    }
	
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDerivativeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,
                              TinyVector<Real_wp, nb_quad> & Qminus,
			      TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
			      TinyVector<Real_wp, nb_quad> & vect_out,
                              PianoString<nb_base, nb_quad, TypeEquation>& piano_string,
			      int k1, int k2)
    {
      if (k1 == 0)
	{
	  if(k2 == 0)
	    vect_out.Fill(piano_string.E*0.5);
	  else if (k2 == 1)
	    vect_out.Fill(piano_string.G);
	}
      else if(k1 == 1)
	{
	  if(k2 == 0)
	    vect_out.Fill(piano_string.A);
	  else if(k2 == 1)
	    vect_out.Fill(piano_string.I*0.5);
	}
    }

    //! application of the scheme    
    template<int nb_base> static inline void
    ApplyMatrixEquation(TinyMatrix<TinyMatrix<Real_wp, General, nb_base, nb_base>,
                        General, dimension, dimension>& mat_Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& res_Un)
    {
      Mlt(mat_Un(0, 0), Un(0), res_Un(0));
      MltAdd(mat_Un(0, 1), Un(1), res_Un(0));        
      Mlt(mat_Un(1, 0), Un(0), res_Un(1));
      MltAdd(mat_Un(1, 1), Un(1), res_Un(1));
    }

    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {
      WaveEquationStiffNLStringLoop<nb_quad-1>::ApplyNonLinearTerm(coefNL, dTm, nablaU);
    }
    
  };
  
  
  // *******************************
  // Equation de corde non linéaire tout NL
  // *******************************
  
  
  // classe pipo
  class WaveEquationNLPipo : public GenericEquationString
  {
  public :
    enum{dimension = 2}; // connu a la compil, super pour templater
    enum { type_equation = WAVE_NL};
    static const bool linear = false;
    TinyVector<bool, dimension> Dirichlet_Agraffe;
    TinyVector<bool, dimension> Dirichlet_Bridge;
    
    bool only_one_LM;
    
    WaveEquationNLPipo()
    {
      Dirichlet_Agraffe(0) = true;
      Dirichlet_Agraffe(1) = true;
      Dirichlet_Bridge(0) = true;
      Dirichlet_Bridge(1) = true;
      
      only_one_LM = false;
    }
    
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ApplyM(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation>& piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 1, vect_out);
    } 
    
    template<int nb_base, int nb_quad, class TypeEquation> inline void
    ApplyR(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      for(int i = 0; i < dimension; i++)
	{
	  TinyVector<Real_wp,nb_quad> ligne_temp;
	  GetRow(vect_in, i, ligne_temp);
	  SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(i)*ligne_temp, i, vect_out);
	}
    } 
	
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ApplyAdamping(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                  TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
    }
	
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ApplyAGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	       TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	       PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
    } 
    
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ApplyAStiffGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
		    TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
		    PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
    } 
    
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ApplyB(TinyMatrix<Real_wp, General, dimension, nb_quad>& vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad>& vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation>& piano_string)
    {
      vect_out.Zero();
    } 
    
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ApplyC(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
    } 
    
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ComputeU( TinyMatrix<Real_wp, General, dimension, nb_quad> & Q,
	      TinyVector<Real_wp, nb_quad> & vect_out,
              PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
      for (int g = 0; g < nb_quad; g++)
	{
	  Real_wp racine = sqrt(Q(0,g) *Q(0,g)  + (1.0+Q(1,g))*(1.0+Q(1,g)));	
	  vect_out(g) = piano_string.T0*Q(0,g) *Q(0,g) *0.5
            + piano_string.E*piano_string.A*Q(1,g) *Q(1,g) *0.5
	    + coef*(Q(0,g) *Q(0,g) *0.5 + (1.0 + Q(1,g)) - racine);
	}      
    }
    
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ComputeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus, TinyVector<Real_wp, nb_quad> & Qminus,
		    TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
		    TinyVector<Real_wp, nb_quad> & vect_out,
                    PianoString<nb_base, nb_quad, TypeEquation> & piano_string, int k)
    {
      vect_out.Fill(0.0);
      if (k == 0)
	{
	  // ligne u
	  for (int g = 0; g < nb_quad; g++)
	    {
	      // A DEROULER
	      Real_wp s1 = sqrt(Qplus(g) *Qplus(g)  + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
	      Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
	      vect_out(g) = piano_string.E*piano_string.A*(Qplus(g)+Qminus(g))*0.5
		- (piano_string.E*piano_string.A - piano_string.T0) 
                * ( ( 1.0/(s1+s2))* (Qplus(g)+Qminus(g)) );
	    }
	}
      if (k == 1)
	{
	  // ligne v
	  for (int g = 0; g < nb_quad; g++)
	    {
	      // A DEROULER
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qplus (g))*(1.0+Qplus (g)) );
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qminus(g))*(1.0+Qminus(g)) );
	      vect_out(g) =  piano_string.E*piano_string.A*(Qplus(g)+Qminus(g))*0.5 
		+ (piano_string.E*piano_string.A - piano_string.T0) 
                * (1.0 - (2.0+Qplus(g)+Qminus(g))/(s1+s2));
	    }
	}
    }
    
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ComputeDerivativeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,
                              TinyVector<Real_wp, nb_quad> & Qminus,
			      TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
			      TinyVector<Real_wp, nb_quad> & vect_out,
			      PianoString<nb_base, nb_quad, TypeEquation> & piano_string,
                              int k1, int k2)
    {
      vect_out.Zero();
      if (k1 == 0)
	{
	  if (k2 == 0)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
		  Real_wp s1 = sqrt(Qplus(g) *Qplus(g)  + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
		  Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = piano_string.E*piano_string.A*0.5 -coef
		    * (  1.0/(s1+s2) -  (Qplus(g)+Qminus(g))*Qplus(g)/(s1*(s1+s2)*(s1+s2)) );
		}
	    }
	  else if (k2 == 1)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
		  Real_wp s1 = sqrt(Qplus(g) *Qplus(g)  + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
		  Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (Qplus(g)+Qminus(g))*(1.0+Qsigma(0,g))
                                        *(1.0/s1+1.0/s2)/( (s1+s2)*(s1+s2) ) );
		}
	    }
	}      
      else if (k1 == 1)
	{
	  if (k2 == 0)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
		  Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qplus (g))*(1.0+Qplus (g)) );
		  Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qminus(g))*(1.0+Qminus(g)) );
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (2.0+Qplus(g)+Qminus(g))*Qsigma(0,g)
                                        *(1.0/s1 + 1.0/s2) )/( (s1+s2)*(s1+s2) );
		}
	    }
	  else if (k2 == 1)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
		  Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qplus (g))*(1.0+Qplus (g)) );
		  Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qminus(g))*(1.0+Qminus(g)) );
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = piano_string.E*piano_string.A*0.5  + coef
		    * (-1.0/(s1+s2) +  (2.0+Qplus(g)+Qminus(g))
                       *(1.0+Qplus(g))/(s1*(s1+s2)*(s1+s2)) );
		}
	    }
	}      
    }
    
    //! application of the scheme    
    template<int nb_base> static inline void
    ApplyMatrixEquation(TinyMatrix<TinyMatrix<Real_wp, General, nb_base, nb_base>,
                        General, dimension, dimension>& mat_Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& res_Un)
    {
      Mlt(mat_Un(0, 0), Un(0), res_Un(0));
      MltAdd(mat_Un(0, 1), Un(1), res_Un(0));        
      Mlt(mat_Un(1, 0), Un(0), res_Un(1));
      MltAdd(mat_Un(1, 1), Un(1), res_Un(1));
    }

    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {
      WaveEquationStiffNLStringLoop<nb_quad-1>::ApplyNonLinearTerm(coefNL, dTm, nablaU);
    }

  };
  
    
  // *******************************
  // Equation de corde non linéaire Raide
  // *******************************
  
  
  //! Modele complet de la these de Juliette Chabassier (inconnues : u, v phi)
  /*!
    Le deplacement transverse u, le deplacement longitudinal v et l'angle phi sont solutions de
    rho A d^2 u/dt^2 - (A G K + T0) d^2 u/dx^2 + A G K dphi/dx 
             - d/dx( (E A - T0) du/dx (1 - 1 / sqrt( (du/dx)^2 + (1 + dv/dv)^2 ) ) ) = S(x)
    rho A d^2 v/dt^2 - E A d^2 v/dx^2
             - d/dx( (E A - T0) (1 - 1 / sqrt( (du/dx)^2 + (1 + dv/dv)^2 ) ) ) = 0
    rho I d^2 phi/dt^2 - E I d^2 phi/dx^2 - A G K du/dx + A G K phi = 0
   */
  class WaveEquationStiffNLString : public GenericEquationString
  {
  public :
    enum{dimension = 3}; // connu a la compil, super pour templater
    enum{ type_equation = WAVE_NL_TIMO};
    static const bool linear = false;
    TinyVector<bool, dimension> Dirichlet_Agraffe;
    TinyVector<bool, dimension> Dirichlet_Bridge;
    
    bool only_one_LM;
    
    //! default constructor
    WaveEquationStiffNLString()
    {
      Dirichlet_Agraffe(0) = true;
      Dirichlet_Agraffe(1) = true;
      Dirichlet_Agraffe(2) = false;
      Dirichlet_Bridge(0) = true;
      Dirichlet_Bridge(1) = true;
      Dirichlet_Bridge(2) = false;
      
      only_one_LM = false;
    }
    
    //! application de la matrice de masse, vect_out = M vect_in
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ApplyM(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 1, vect_out);
      
      GetRow(vect_in, 2, ligne_temp);
      SetRow(piano_string.rho*piano_string.I*ligne_temp, 2, vect_out);
    } 
    
    //! application de la matrice d'amortissement, vect_out = R vect_in
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ApplyR(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(1)*ligne_temp, 1, vect_out);

      GetRow(vect_in, 2, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.I*piano_string.amo(2)*ligne_temp, 2, vect_out);
    } 
    
    //! application de l'amortissement quadratique, vect_out = H vect_in
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ApplyAdamping(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                  TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.T0*piano_string.amoB2(0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.A*piano_string.amoB2(1)*ligne_temp, 1, vect_out);

      GetRow(vect_in, 2, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.I*piano_string.amoB2(2)*ligne_temp, 2, vect_out);
    }
    
    //! application de la rigidite non-raide, vect_out = A^grad vect_in
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ApplyAGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	       TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	       PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.T0)*ligne_temp, 0, vect_out);
    } 
    
    //! application de la rigidite raide, vect_out = A^stiff vect_in
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ApplyAStiffGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
		    TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
		    PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.A*piano_string.G*piano_string.k_prime)*ligne_temp, 0, vect_out);
      
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow((piano_string.E*piano_string.A)*ligne_temp, 1, vect_out);
      
      GetRow(vect_in, 2, ligne_temp);
      SetRow((piano_string.E*piano_string.I)*ligne_temp, 2, vect_out);
    } 
    
    //! application de la matrice B, vect_out = B vect_in
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ApplyB(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 2, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 0, vect_out);
      
      //  GetRow(vect_in, 0, ligne_temp);
      //  SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 2, vect_out);
    } 
        
    //! application de la matrice C, vect_out = C vect_in
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ApplyC(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 2, ligne_temp);
      SetRow(piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 2, vect_out);
    } 
	
    //! application de la transposee de la matrice B, vect_out = B^T vect_in
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ApplyD(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      //	  GetRow(vect_in, 2, ligne_temp);
      //	  SetRow(- piano_string.A*piano_string.G*piano_string
      //               .k_prime*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 2, vect_out);
    }

    //! Calcul de vect_out = U(Q)	
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ComputeU( TinyMatrix<Real_wp, General, dimension, nb_quad> & Q,
              TinyVector<Real_wp, nb_quad> & vect_out,
	      PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
      for (int g = 0; g < nb_quad; g++)
	{
	  Real_wp racine = sqrt(Q(0,g) *Q(0,g)  + (1.0+Q(1,g))*(1.0+Q(1,g)));	
	  vect_out(g) = coef*(Q(0,g) *Q(0,g) *0.5 + (1.0 + Q(1,g)) - racine);
	}
    }
    
    //! calcul de vect_out = grad(U) (q) a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ComputeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus, TinyVector<Real_wp, nb_quad> & Qminus,
		    TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
		    TinyVector<Real_wp, nb_quad> & vect_out,
		    PianoString<nb_base, nb_quad, TypeEquation> & piano_string, int k)
    {
      vect_out.Fill(0.0);
      if (k == 0)
	WaveEquationStiffNLStringLoop<nb_quad-1>::ComputeDeltakUk0(Qplus, Qminus, Qsigma,
                                                                   vect_out, piano_string);
      
      if(k == 1)
	WaveEquationStiffNLStringLoop<nb_quad-1>::ComputeDeltakUk1(Qplus, Qminus, Qsigma,
                                                                   vect_out, piano_string);
    }
    
    //! calcul de vect_out = d/dq( grad(U)(q))  a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_base, int nb_quad, class TypeEquation> inline static void
    ComputeDerivativeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,
                              TinyVector<Real_wp, nb_quad> & Qminus,
			      TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
			      TinyVector<Real_wp, nb_quad> & vect_out,
			      PianoString<nb_base, nb_quad, TypeEquation> & piano_string,
                              int k1, int k2)
    {
      vect_out.Zero();
      if (k1 == 0)
	{
	  if (k2 == 0)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
		  Real_wp s1 = sqrt(Qplus(g) *Qplus(g)  + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
		  Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( 0.5 - 1.0/(s1+s2) 
                                        +  (Qplus(g)+Qminus(g))*Qplus(g)/(s1*(s1+s2)*(s1+s2)) );
		}
	    }
	  else if (k2 == 1)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
		  Real_wp s1 = sqrt(Qplus(g) *Qplus(g)  + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
		  Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (1.0+Qsigma(0,g))*(1.0+Qsigma(0,g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (Qplus(g)+Qminus(g))
                                        *(1.0+Qsigma(0,g))*(1.0/s1+1.0/s2)/( (s1+s2)*(s1+s2) ) );
		}
	    }
	}      
      else if (k1 == 1)
	{
	  if (k2 == 0)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
		  Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qplus (g))*(1.0+Qplus (g)) );
		  Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qminus(g))*(1.0+Qminus(g)) );
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (2.0+Qplus(g)+Qminus(g))*Qsigma(0,g)
                                        *(1.0/s1 + 1.0/s2) )/( (s1+s2)*(s1+s2) );
		}
	    }
	  else if (k2 == 1)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
		  Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qplus (g))*(1.0+Qplus (g)) );
		  Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (1.0+Qminus(g))*(1.0+Qminus(g)) );
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef * (-1.0/(s1+s2) +  (2.0+Qplus(g)+Qminus(g))
                                        *(1.0+Qplus(g))/(s1*(s1+s2)*(s1+s2)) );
		}
	    }
	}
    }

    //! application of the scheme    
    template<int nb_base> static inline void
    ApplyMatrixEquation(TinyMatrix<TinyMatrix<Real_wp, General, nb_base, nb_base>,
                        General, dimension, dimension>& mat_Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& res_Un)
    {
      Mlt(mat_Un(0, 0), Un(0), res_Un(0));
      MltAdd(mat_Un(0, 2), Un(2), res_Un(0));        
      Mlt(mat_Un(1, 1), Un(1), res_Un(1));
      Mlt(mat_Un(2, 0), Un(0), res_Un(2));
      MltAdd(mat_Un(2, 2), Un(2), res_Un(2));
    }
    
    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {
      WaveEquationStiffNLStringLoop<nb_quad-1>::ApplyNonLinearTerm(coefNL, dTm, nablaU);
    }
    
  };


  // *******************************
  // Equation de corde linéaire Raide mais a trois inconnues pour feinter le only_one_LM qui foire
  // *******************************
  
  
  //! Modele complet de corde de Juliette Chabassier sans non-linearite (3 inconnues : u, v, phi)
  /*!
    Le deplacement transverse u, le deplacement longitudinal v et l'angle phi sont solutions de
    rho A d^2 u/dt^2 - (A G K + T0) d^2 u/dx^2 + A G K dphi/dx = S(x)
    rho A d^2 v/dt^2 - E A d^2 v/dx^2 = 0
    rho I d^2 phi/dt^2 - E I d^2 phi/dx^2 - A G K du/dx + A G K phi = 0
   */
  class WaveEquationStiffLINString : public GenericEquationString
  {
  public :
    enum{dimension = 3}; // connu a la compil, super pour templater
    enum{ type_equation = WAVE_NL_TIMO};
    static const bool linear = true;
    TinyVector<bool, dimension> Dirichlet_Agraffe;
    TinyVector<bool, dimension> Dirichlet_Bridge;
    
    bool only_one_LM;
    
    //! default constructor
    WaveEquationStiffLINString()
    {
      Dirichlet_Agraffe(0) = true;
      Dirichlet_Agraffe(1) = true;
      Dirichlet_Agraffe(2) = false;
      Dirichlet_Bridge(0) = true;
      Dirichlet_Bridge(1) = true;
      Dirichlet_Bridge(2) = false;
      
      only_one_LM = false;
    }

    
    //! application de la matrice de masse, vect_out = M vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyM(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 1, vect_out);
      
      GetRow(vect_in, 2, ligne_temp);
      SetRow(piano_string.rho*piano_string.I*ligne_temp, 2, vect_out);
    } 
    
    //! application de la matrice d'amortissement, vect_out = R vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyR(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {      
      TinyVector<Real_wp,nb_quad> ligne_temp;

      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(1)*ligne_temp, 1, vect_out);

      GetRow(vect_in, 2, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.I*piano_string.amo(2)*ligne_temp, 2, vect_out);
    } 
    
    //! application de l'amortissement quadratique, vect_out = H vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAdamping(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                  TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.T0*piano_string.amoB2(0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.A*piano_string.amoB2(1)*ligne_temp, 1, vect_out);

      GetRow(vect_in, 2, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.I*piano_string.amoB2(2)*ligne_temp, 2, vect_out); 
    }
    
    //! application de la rigidite non-raide, vect_out = A^grad vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	       TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	       PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.T0)*ligne_temp, 0, vect_out);
    } 
    
    //! application de la rigidite raide, vect_out = A^stiff vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAStiffGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
		    TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
		    PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.A*piano_string.G*piano_string.k_prime)*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow((piano_string.E*piano_string.A)*ligne_temp, 1, vect_out);
      
      GetRow(vect_in, 2, ligne_temp);
      SetRow((piano_string.E*piano_string.I)*ligne_temp, 2, vect_out);
    } 
    
    //! application de la matrice B, vect_out = B vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyB(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 2, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 0, vect_out);
      
      //  GetRow(vect_in, 0, ligne_temp);
      //  SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 2, vect_out);
    } 
    
    //! application de la matrice C, vect_out = C vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyC(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 2, ligne_temp);
      SetRow(piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 2, vect_out);      
    } 
    
    //! application de la transposee de la matrice B, vect_out = B^T vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyD(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      //	  GetRow(vect_in, 2, ligne_temp);
      //	  SetRow(- piano_string.A*piano_string.G*piano_string
      //    .k_prime*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 2, vect_out);
    }
    
    //! Calcul de vect_out = U(Q)
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeU( TinyMatrix<Real_wp, General, dimension, nb_quad> & Q,
              TinyVector<Real_wp, nb_quad> & vect_out,
	      PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();	
    }
    
    //! calcul de vect_out = grad(U) (q) a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus, TinyVector<Real_wp, nb_quad> & Qminus,
		    TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
		    TinyVector<Real_wp, nb_quad> & vect_out,
		    PianoString<nb_base, nb_quad, TypeEquation> & piano_string, int k)
    {
      vect_out.Fill(0.0);
    }
    
    //! calcul de vect_out = d/dq( grad(U)(q))  a partir de Q^n-1, Q^n+1 et Q^n	
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDerivativeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,
                              TinyVector<Real_wp, nb_quad> & Qminus,
			      TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
			      TinyVector<Real_wp, nb_quad> & vect_out,
			      PianoString<nb_base, nb_quad, TypeEquation> & piano_string,
                              int k1, int k2)
    {
      vect_out.Zero();
    }

    //! application of the scheme    
    template<int nb_base> static inline void
    ApplyMatrixEquation(TinyMatrix<TinyMatrix<Real_wp, General, nb_base, nb_base>,
                        General, dimension, dimension>& mat_Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& res_Un)
    {
      Mlt(mat_Un(0, 0), Un(0), res_Un(0));
      MltAdd(mat_Un(0, 2), Un(2), res_Un(0));        
      Mlt(mat_Un(1, 1), Un(1), res_Un(1));
      Mlt(mat_Un(2, 0), Un(0), res_Un(2));
      MltAdd(mat_Un(2, 2), Un(2), res_Un(2));
    }
    
    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {
    }

  };
  
  
  // *******************************
  // Equation de corde non linéaire Raide
  // *******************************
  
  
  //! comme WaveEquationStiffLINString !
  class WaveEquationStiffLinString : public GenericEquationString 
  {
  public :
    enum{dimension = 3}; // connu a la compil, super pour templater
    enum{ type_equation = WAVE_NL_TIMO};
    static const bool linear = true;
    TinyVector<bool, dimension> Dirichlet_Agraffe;
    TinyVector<bool, dimension> Dirichlet_Bridge;
    
    bool only_one_LM;
    
    //! default constructor
    WaveEquationStiffLinString()
    {
      Dirichlet_Agraffe(0) = true;
      Dirichlet_Agraffe(1) = true;
      Dirichlet_Agraffe(2) = false;
      Dirichlet_Bridge(0) = true;
      Dirichlet_Bridge(1) = true;
      Dirichlet_Bridge(2) = false;
      
      only_one_LM = false;
    }
    
    
    //! application de la matrice de masse, vect_out = M vect_in	
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyM(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 0, vect_out);
	  
      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 1, vect_out);
      
      GetRow(vect_in, 2, ligne_temp);
      SetRow(piano_string.rho*piano_string.I*ligne_temp, 2, vect_out);
    } 
    
    //! application de la matrice d'amortissement, vect_out = R vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyR(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(1)*ligne_temp, 1, vect_out);
      
      GetRow(vect_in, 2, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.I*piano_string.amo(2)*ligne_temp, 2, vect_out);
    } 
    
    //! application de l'amortissement quadratique, vect_out = H vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAdamping(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                  TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.T0*piano_string.amoB2(0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.A*piano_string.amoB2(1)*ligne_temp, 1, vect_out);

      GetRow(vect_in, 2, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.I*piano_string.amoB2(2)*ligne_temp, 2, vect_out);
    }
	
    //! application de la rigidite non-raide, vect_out = A^grad vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	       TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	       PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();   
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.T0)*ligne_temp, 0, vect_out);
    } 
    
    //! application de la rigidite raide, vect_out = A^stiff vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAStiffGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
		    TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
		    PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.A*piano_string.G*piano_string.k_prime)*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow((piano_string.E*piano_string.A)*ligne_temp, 1, vect_out);
      
      GetRow(vect_in, 2, ligne_temp);
      SetRow((piano_string.E*piano_string.I)*ligne_temp, 2, vect_out);
    } 
    
    //! application de la matrice B, vect_out = B vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyB(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp,nb_quad> ligne_temp;
      GetRow(vect_in, 2, ligne_temp);
      SetRow(- piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 0, vect_out);
      
      // GetRow(vect_in, 0, ligne_temp);
      // SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 2, vect_out);      
    } 
    
    //! application de la matrice C, vect_out = C vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyC(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
      
      TinyVector<Real_wp,nb_quad> ligne_temp;
      GetRow(vect_in, 2, ligne_temp);
      SetRow(piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 2, vect_out);
    } 
    
    //! application de la transposee de la matrice B, vect_out = B^T vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyD(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp,nb_quad> ligne_temp;
      //	  GetRow(vect_in, 2, ligne_temp);
      //	  SetRow(- piano_string.A*piano_string.G*piano_string
      //               .k_prime*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 2, vect_out);
    }
    
    //! Calcul de vect_out = U(Q)
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeU( TinyMatrix<Real_wp,General,dimension, nb_quad> & Q,
              TinyVector<Real_wp, nb_quad> & vect_out,
	      PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();      
    }
    
    //! calcul de vect_out = grad(U) (q) a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus, TinyVector<Real_wp, nb_quad> & Qminus,
		    TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
		    TinyVector<Real_wp, nb_quad> & vect_out,
		    PianoString<nb_base, nb_quad, TypeEquation> & piano_string, int k)
    {
      vect_out.Fill(0.0);
    }
	
    //! calcul de vect_out = d/dq( grad(U)(q))  a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDerivativeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,
                              TinyVector<Real_wp, nb_quad> & Qminus,
			      TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
			      TinyVector<Real_wp, nb_quad> & vect_out,
			      PianoString<nb_base, nb_quad, TypeEquation> & piano_string,
                              int k1, int k2)
    {
      vect_out.Zero();
    }
    
    //! application of the scheme
    template<int nb_base> static inline void
    ApplyMatrixEquation(TinyMatrix<TinyMatrix<Real_wp, General, nb_base, nb_base>,
                        General, dimension, dimension>& mat_Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& res_Un)
    {
      Mlt(mat_Un(0, 0), Un(0), res_Un(0));
      MltAdd(mat_Un(0, 2), Un(2), res_Un(0));        
      Mlt(mat_Un(1, 1), Un(1), res_Un(1));
      Mlt(mat_Un(2, 0), Un(0), res_Un(2));
      MltAdd(mat_Un(2, 2), Un(2), res_Un(2));
    }

    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {
    }

  };


  class WaveEquationStringLINTL : public GenericEquationString
  {
  public :
    enum{dimension = 2}; // connu a la compil, super pour templater
    enum { type_equation = WAVE_NL};
    static const bool linear = true;
    
    TinyVector<bool,dimension> Dirichlet_Agraffe;
    TinyVector<bool,dimension> Dirichlet_Bridge;
    
    bool only_one_LM;
    
    //! default constructor
    WaveEquationStringLINTL()
    {
      Dirichlet_Agraffe(0) = true;
      Dirichlet_Agraffe(1) = true;
      Dirichlet_Bridge(0) = true;
      Dirichlet_Bridge(1) = true;
      
      only_one_LM = false;
    }

    //! application de la matrice de masse, vect_out = M vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyM(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 1, vect_out);
    } 
    
    //! application de la matrice d'amortissement, vect_out = R vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyR(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*ligne_temp, 1, vect_out);
    }
    
    //! application de l'amortissement quadratique, vect_out = H vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyAdamping(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                  TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                  const GenericPb& piano_string)
    {
      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.T0*piano_string.amoB2(0)*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.A*piano_string.amoB2(0)*ligne_temp, 1, vect_out);
    }
    
    //! application de la rigidite non-raide, vect_out = A^grad vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyAGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
               TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
               const GenericPb& piano_string)
    {      
      vect_out.Zero();
      
      TinyVector<Real_wp,nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.T0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow((piano_string.E*piano_string.A)*ligne_temp, 1, vect_out);
    } 
    
    //! application de la rigidite raide, vect_out = A^stiff vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyAStiffGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                    TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                    const GenericPb& piano_string)
    {
      vect_out.Zero();
      
      vect_out.Fill(0.0);
    } 
    
    //! application de la matrice B, vect_out = B vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyB(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      vect_out.Fill(0.0);
    } 
    
    //! application de la matrice C, vect_out = C vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyC(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      vect_out.Fill(0.0);
      
    } 
        
    //! application de la transposee de la matrice B, vect_out = B^T vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyD(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      vect_out.Fill(0.0);
    } 
    
    //! Calcul de vect_out = U(Q)
    template<int nb_quad, class GenericPb> static inline void
    ComputeU( TinyMatrix<Real_wp,General,dimension, nb_quad> & Q,
              TinyVector<Real_wp,nb_quad> & vect_out, GenericPb & piano_string)
    {
      vect_out.Zero();
      
    }
    
    //! calcul de vect_out = grad(U) (q) a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_quad, class GenericPb> static inline void
    ComputeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus, TinyVector<Real_wp, nb_quad> & Qminus,
                    TinyMatrix<Real_wp,General,dimension-1,nb_quad> & Qsigma,
                    TinyVector<Real_wp,nb_quad> & vect_out, GenericPb & piano_string, int k)
    {
      vect_out.Fill(0.0);
    }

    //! calcul de vect_out = d/dq( grad(U)(q))  a partir de Q^n-1, Q^n+1 et Q^n    
    template<int nb_quad, class GenericPb> static inline void
    ComputeDerivativeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,
                              TinyVector<Real_wp, nb_quad> & Qminus,
                              TinyMatrix<Real_wp,General,dimension-1,nb_quad> & Qsigma,
                              TinyVector<Real_wp, nb_quad> & vect_out,
                              GenericPb & piano_string, int k1, int k2)
    {
      vect_out.Zero();      
    }

    //! application of the scheme    
    template<int nb_base> static inline void
    ApplyMatrixEquation(TinyMatrix<TinyMatrix<Real_wp, General, nb_base, nb_base>,
                        General, dimension, dimension>& mat_Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& res_Un)
    {
      Mlt(mat_Un(0, 0), Un(0), res_Un(0));
      MltAdd(mat_Un(0, 1), Un(1), res_Un(0));        
      Mlt(mat_Un(1, 0), Un(0), res_Un(1));
      MltAdd(mat_Un(1, 1), Un(1), res_Un(1));
    }

    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {
    }
    
  };

// two linear wave equations
 class WaveEquationStringLIN2T : public GenericEquationString
  {
  public :
    enum{dimension = 2}; // connu a la compil, super pour templater
    enum { type_equation = WAVE_LIN_2T};
    static const bool linear = true;
    
    TinyVector<bool,dimension> Dirichlet_Agraffe;
    TinyVector<bool,dimension> Dirichlet_Bridge;
    
    bool only_one_LM;
    
    //! default constructor
    WaveEquationStringLIN2T()
    {
      Dirichlet_Agraffe(0) = true;
      Dirichlet_Agraffe(1) = true;
      Dirichlet_Bridge(0) = true;
      Dirichlet_Bridge(1) = true;
      
      only_one_LM = false;
    }

    //! application de la matrice de masse, vect_out = M vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyM(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.rho*piano_string.rho_detuning*piano_string.A*ligne_temp, 1, vect_out);
    } 
    
    //! application de la matrice d'amortissement, vect_out = R vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyR(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.rho_detuning*piano_string.A*piano_string.amo(0)*ligne_temp, 1, vect_out);
    }
    
    //! application de l'amortissement quadratique, vect_out = H vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyAdamping(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                  TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                  const GenericPb& piano_string)
    {
      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.T0*piano_string.amoB2(0)*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.T0*piano_string.amoB2(0)*ligne_temp, 1, vect_out);
    }
    
    //! application de la rigidite non-raide, vect_out = A^grad vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyAGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
               TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
               const GenericPb& piano_string)
    {      
      vect_out.Zero();
      
      TinyVector<Real_wp,nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.T0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow((piano_string.T0)*ligne_temp, 1, vect_out);
    } 
    
    //! application de la rigidite raide, vect_out = A^stiff vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyAStiffGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                    TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                    const GenericPb& piano_string)
    {
      vect_out.Zero();
      
      vect_out.Fill(0.0);
    } 
    
    //! application de la matrice B, vect_out = B vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyB(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      vect_out.Fill(0.0);
    } 
    
    //! application de la matrice C, vect_out = C vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyC(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      vect_out.Fill(0.0);
      
    } 
        
    //! application de la transposee de la matrice B, vect_out = B^T vect_in
    template<int nb_quad, class GenericPb> static inline void
    ApplyD(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
           TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
           const GenericPb& piano_string)
    {
      vect_out.Fill(0.0);
    } 
    
    //! Calcul de vect_out = U(Q)
    template<int nb_quad, class GenericPb> static inline void
    ComputeU( TinyMatrix<Real_wp,General,dimension, nb_quad> & Q,
              TinyVector<Real_wp,nb_quad> & vect_out, GenericPb & piano_string)
    {
      vect_out.Zero();
      
    }
    
    //! calcul de vect_out = grad(U) (q) a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_quad, class GenericPb> static inline void
    ComputeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus, TinyVector<Real_wp, nb_quad> & Qminus,
                    TinyMatrix<Real_wp,General,dimension-1,nb_quad> & Qsigma,
                    TinyVector<Real_wp,nb_quad> & vect_out, GenericPb & piano_string, int k)
    {
      vect_out.Fill(0.0);
    }

    //! calcul de vect_out = d/dq( grad(U)(q))  a partir de Q^n-1, Q^n+1 et Q^n    
    template<int nb_quad, class GenericPb> static inline void
    ComputeDerivativeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,
                              TinyVector<Real_wp, nb_quad> & Qminus,
                              TinyMatrix<Real_wp,General,dimension-1,nb_quad> & Qsigma,
                              TinyVector<Real_wp, nb_quad> & vect_out,
                              GenericPb & piano_string, int k1, int k2)
    {
      vect_out.Zero();      
    }

    //! application of the scheme    
    template<int nb_base> static inline void
    ApplyMatrixEquation(TinyMatrix<TinyMatrix<Real_wp, General, nb_base, nb_base>,
                        General, dimension, dimension>& mat_Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& res_Un)
    {
      Mlt(mat_Un(0, 0), Un(0), res_Un(0));
      MltAdd(mat_Un(0, 1), Un(1), res_Un(0));        
      Mlt(mat_Un(1, 0), Un(0), res_Un(1));
      MltAdd(mat_Un(1, 1), Un(1), res_Un(1));
    }

    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {
    }
    
  };

//non-planar non-linear wave equation//

  class WaveEquationStringNL2T : public GenericEquationString
  {
  public :
    enum{dimension = 3}; // connu a la compil, super pour templater
    enum { type_equation = WAVE_NL_2T};
    static const bool linear = false;
    
    TinyVector<bool,dimension> Dirichlet_Agraffe;
    TinyVector<bool,dimension> Dirichlet_Bridge;
    
    bool only_one_LM;
    
    //! default constructor
    WaveEquationStringNL2T()
    {
      Dirichlet_Agraffe(0) = true;
      Dirichlet_Agraffe(1) = true;
      Dirichlet_Agraffe(2) = true;
      Dirichlet_Bridge(0) = true;
      Dirichlet_Bridge(1) = true;
      Dirichlet_Bridge(2) = true;
      
      only_one_LM = false;
    }

    //! produit scalaire pour separer energie u et energie v
    template<class T, int m>
    static void DotProdColUV(const TinyMatrix<T, General, 3, m>& A,
			     const TinyMatrix<T, General, 3, m>& B,
			     TinyVector<T, m>& res_u, TinyVector<T, m>& res_v)
    {
      for (int p = 0; p < m; p++)
	{
	  res_u(p) += A(0, p)*B(0, p) + A(1, p)*B(1, p);
	  res_v(p) += A(2, p)*B(2, p);
	}
    }
        
    //! application de la matrice de masse, vect_out = M vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyM(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.rho_detuning*piano_string.rho*piano_string.A*ligne_temp, 1, vect_out);
      
	  GetRow(vect_in, 2, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 2, vect_out);
    } 
    
    //! application de la matrice d'amortissement, vect_out = R vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyR(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.rho_detuning*2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*ligne_temp, 1, vect_out);
      
      GetRow(vect_in, 2, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(1)*ligne_temp, 2, vect_out);
    }
    
    //! application de l'amortissement quadratique, vect_out = H vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAdamping(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                  TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.T0*piano_string.amoB2(0)*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.T0*piano_string.amoB2(0)*ligne_temp, 1, vect_out);  
          
      GetRow(vect_in, 2, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.A*piano_string.amoB2(1)*ligne_temp, 2, vect_out);
    }
    
    //! application de la rigidite non-raide, vect_out = A^grad vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	       TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	       PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {      
      vect_out.Zero();
      
      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.T0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow((piano_string.T0)*ligne_temp, 1, vect_out);
      

    } 
    
    //! application de la rigidite raide, vect_out = A^stiff vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAStiffGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
		    TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
		    PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
	  TinyVector<Real_wp, nb_quad> ligne_temp;
	  GetRow(vect_in, 2, ligne_temp);
      SetRow((piano_string.E*piano_string.A)*ligne_temp, 2, vect_out);
      
      //~ vect_out.Fill(0.0);
    } 
    
    //! application de la matrice B, vect_out = B vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyB(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
    } 
    
    //! application de la matrice C, vect_out = C vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyC(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
      
    } 
        
    //! application de la transposee de la matrice B, vect_out = B^T vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyD(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
    } 
    
    //! Calcul de vect_out = U(Q)
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeU( TinyMatrix<Real_wp, General, dimension, nb_quad> & Q,
	      TinyVector<Real_wp, nb_quad> & vect_out,
              PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    { 
      Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
      for (int g = 0; g < nb_quad; g++)
	{
	  Real_wp racine = sqrt(Q(0,g) *Q(0,g)  + (Q(1,g))*(Q(1,g))+(1+Q(2,g))*(1+Q(2,g)));	
	  vect_out(g) = coef*(Q(0,g) *Q(0,g) *0.5 + Q(1,g) *Q(1,g) *0.5 + (1+Q(2,g))- racine);

	}  

    }
    
    //! calcul de vect_out = grad(U) (q) a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,TinyVector<Real_wp, nb_quad> & Qminus,
		    TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
		    TinyVector<Real_wp, nb_quad> & vect_out,
                    PianoString<nb_base, nb_quad, TypeEquation> & piano_string, int k)
    {
      vect_out.Fill(0.0);
      if (k==0)
	{ 
	  // ligne u
	  for (int g = 0; g < nb_quad; g++)
	    {
	      // A DEROULER
	      Real_wp s1 = sqrt(Qplus (g)*Qplus (g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
	      Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
	      vect_out(g) = (piano_string.E*piano_string.A - piano_string.T0)
                * ( (0.5 - 1.0/(s1+s2))* (Qplus(g)+Qminus(g)) );
	    }

	}
      
      if (k==1)
	{
	  // ligne v
	  for (int g = 0; g < nb_quad; g++)
	    {
	    // A DEROULER
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qplus (g))*(Qplus (g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qminus(g))*(Qminus(g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
	      vect_out(g) = (piano_string.E*piano_string.A - piano_string.T0)
                * ( (0.5 - 1.0/(s1+s2))* (Qplus(g)+Qminus(g)) );

	    }
	}
	
	  if (k==2)
	{
	  // ligne w
	  for (int g = 0; g < nb_quad; g++)
	    {
         // A DEROULER
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qplus (g))*(1+Qplus (g)));
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qminus(g))*(1+Qminus(g)));
	      vect_out(g) = (piano_string.E*piano_string.A - piano_string.T0)
                * (1.0 - (2.0+Qplus(g)+Qminus(g))/(s1+s2));

		}
	}
    }
    
    //! calcul de vect_out = d/dq( grad(U)(q))  a partir de Q^n-1, Q^n+1 et Q^n    
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDerivativeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,
                              TinyVector<Real_wp, nb_quad> & Qminus,
			      TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
			      TinyVector<Real_wp, nb_quad> & vect_out,
			      PianoString<nb_base, nb_quad, TypeEquation> & piano_string,
                              int k1, int k2)
    {
      vect_out.Zero();
      if (k1 == 0)
	{
	  if(k2 == 0)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	      Real_wp s1 = sqrt(Qplus (g)*Qplus (g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
	      Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( 0.5 - 1.0/(s1+s2)
                                        +  (Qplus(g)+Qminus(g))*Qplus(g)/(s1*(s1+s2)*(s1+s2)) );

                                        
		}
			    //~ cout << "what is DerivativedeltakU00?  " << vect_out << endl;
	    }

	  else if (k2 == 1)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	      Real_wp s1 = sqrt(Qplus (g)*Qplus (g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
	      Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (Qplus(g)+Qminus(g))*(Qsigma(0,g))*(1.0/s1+1.0/s2)
                                        /( (s1+s2)*(s1+s2) ) );
		}
			    //~ cout << "what is DerivativedeltakU01?  " << vect_out << endl;
	    }

	  else if (k2 == 2)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	      Real_wp s1 = sqrt(Qplus (g)*Qplus (g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
	      Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (Qplus(g)+Qminus(g))*(1.0+Qsigma(1,g))
                                        *(1.0/s1+1.0/s2)/( (s1+s2)*(s1+s2) ) );
		}
					    //~ cout << "what is DerivativedeltakU02?  " << vect_out << endl;
	    }
	}      
      else if (k1 == 1)
	{
	  if (k2 == 0)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qplus (g))*(Qplus (g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qminus(g))*(Qminus(g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (Qplus(g)+Qminus(g))*(Qsigma(0,g))*(1.0/s1+1.0/s2)
                                        /( (s1+s2)*(s1+s2) ) );
		}
			    //~ cout << "what is DerivativedeltakU10?  " << vect_out << endl;
	    }
	  else if (k2 ==1)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qplus (g))*(Qplus (g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qminus(g))*(Qminus(g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( 0.5 - 1.0/(s1+s2)
                                        +  (Qplus(g)+Qminus(g))*Qplus(g)/(s1*(s1+s2)*(s1+s2)) );
		}
			    //~ cout << "what is DerivativedeltakU11?  " << vect_out << endl;
	    }
	  else if (k2 == 2)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qplus (g))*(Qplus (g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qminus(g))*(Qminus(g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (Qplus(g)+Qminus(g))*(1.0+Qsigma(1,g))
                                        *(1.0/s1+1.0/s2)/( (s1+s2)*(s1+s2) ) );
		}
	    //~ cout << "what is DerivativedeltakU12?  " << vect_out << endl;
	    }

	}      
      else if (k1 == 2)
	{
	  if (k2 == 0)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qplus (g))*(1+Qplus (g)));
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qminus(g))*(1+Qminus(g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (2.0+Qplus(g)+Qminus(g))*Qsigma(0,g)*(1.0/s1+1.0/s2)
                                        /( (s1+s2)*(s1+s2) ) );
		}
			    //~ cout << "what is DerivativedeltakU20?  " << vect_out << endl;
	    }
	  else if (k2 ==1)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qplus (g))*(1+Qplus (g)));
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qminus(g))*(1+Qminus(g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (2.0+Qplus(g)+Qminus(g))*Qsigma(1,g)*(1.0/s1+1.0/s2)
                                        /( (s1+s2)*(s1+s2) ) );
		}
			    //~ cout << "what is DerivativedeltakU21?  " << vect_out << endl;
	    }
	  else if (k2 == 2)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qplus (g))*(1+Qplus (g)));
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qminus(g))*(1+Qminus(g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef * (-1.0/(s1+s2) +  (2.0+Qplus(g)+Qminus(g))
                                        *(1.0+Qplus(g))/(s1*(s1+s2)*(s1+s2)) );
		}
			    //~ cout << "what is DerivativedeltakU22?  " << vect_out << endl;
	    }

	}            
    }	


    //! application of the scheme    
    template<int nb_base> static inline void
    ApplyMatrixEquation(TinyMatrix<TinyMatrix<Real_wp, General, nb_base, nb_base>,
                        General, dimension, dimension>& mat_Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& res_Un)
    {
      Mlt(mat_Un(0, 0), Un(0), res_Un(0));     
      Mlt(mat_Un(1, 1), Un(1), res_Un(1));
	  Mlt(mat_Un(2, 2), Un(2), res_Un(2));
    }
    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {
      WaveEquationStiffNL2TStringLoop<nb_quad-1>::ApplyNonLinearTerm(coefNL, dTm, nablaU);
    }
 
  };

class WaveEquationStringStiffNL2T : public GenericEquationString
  {
  public :
    enum{dimension = 5}; // connu a la compil, super pour templater
    enum { type_equation = WAVE_NL_2T_TIMO};
    static const bool linear = false;
    
    TinyVector<bool,dimension> Dirichlet_Agraffe;
    TinyVector<bool,dimension> Dirichlet_Bridge;
    
    bool only_one_LM;
    
    //! default constructor
    WaveEquationStringStiffNL2T()
    {
      Dirichlet_Agraffe(0) = true;
      Dirichlet_Agraffe(1) = true;
      Dirichlet_Agraffe(2) = true;
      Dirichlet_Agraffe(3) = false;
      Dirichlet_Agraffe(4) = false;
      Dirichlet_Bridge(0) = true;
      Dirichlet_Bridge(1) = true;
      Dirichlet_Bridge(2) = true;
      Dirichlet_Bridge(3) = false;
      Dirichlet_Bridge(4) = false;
      
      only_one_LM = false;
    }

    //! produit scalaire pour separer energie u et energie v
    template<class T, int m>
    static void DotProdColUV(const TinyMatrix<T, General, 5, m>& A,
			     const TinyMatrix<T, General, 5, m>& B,
			     TinyVector<T, m>& res_u, TinyVector<T, m>& res_v)
    {
      for (int p = 0; p < m; p++)
	{
	  res_u(p) += A(0, p)*B(0, p) + A(1, p)*B(1, p) + A(3, p)*B(3, p) + A(4, p)*B(4, p);
	  res_v(p) += A(2, p)*B(2, p);
	}
    }
        
    //! application de la matrice de masse, vect_out = M vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyM(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.rho_detuning*piano_string.rho*piano_string.A*ligne_temp, 1, vect_out);
      
	  GetRow(vect_in, 2, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 2, vect_out);
      
      GetRow(vect_in, 3, ligne_temp);
      SetRow(piano_string.rho*piano_string.I*ligne_temp, 3, vect_out);
      
	  GetRow(vect_in, 4, ligne_temp);
      SetRow(piano_string.rho*piano_string.I*ligne_temp, 4, vect_out);   
    } 
    
    //! application de la matrice d'amortissement, vect_out = R vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyR(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.rho_detuning*2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*ligne_temp, 1, vect_out);
      
      GetRow(vect_in, 2, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(1)*ligne_temp, 2, vect_out);
      
      GetRow(vect_in, 3, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.I*piano_string.amo(2)*ligne_temp, 3, vect_out);
      
      GetRow(vect_in, 4, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.I*piano_string.amo(2)*ligne_temp, 4, vect_out);
          }
    
    //! application de l'amortissement quadratique, vect_out = H vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAdamping(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                  TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.T0*piano_string.amoB2(0)*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.T0*piano_string.amoB2(0)*ligne_temp, 1, vect_out);  
          
      GetRow(vect_in, 2, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.A*piano_string.amoB2(1)*ligne_temp, 2, vect_out);

      GetRow(vect_in, 3, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.I*piano_string.amoB2(2)*ligne_temp, 3, vect_out);

      GetRow(vect_in, 4, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.I*piano_string.amoB2(2)*ligne_temp, 4, vect_out);


    }
    
    //! application de la rigidite non-raide, vect_out = A^grad vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	       TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	       PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {      
      vect_out.Zero();
      
      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.T0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow((piano_string.T0)*ligne_temp, 1, vect_out);
      

    } 
    
    //! application de la rigidite raide, vect_out = A^stiff vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAStiffGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
		    TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
		    PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();

      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.A*piano_string.G*piano_string.k_prime)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow((piano_string.A*piano_string.G*piano_string.k_prime)*ligne_temp, 1, vect_out);
      
            
      GetRow(vect_in, 2, ligne_temp);
      SetRow((piano_string.E*piano_string.A)*ligne_temp, 2, vect_out);
      
      GetRow(vect_in, 3, ligne_temp);
      SetRow((piano_string.E*piano_string.I)*ligne_temp, 3, vect_out);
      
      GetRow(vect_in, 4, ligne_temp);
      SetRow((piano_string.E*piano_string.I)*ligne_temp, 4, vect_out);
    } 
    
    //! application de la matrice B, vect_out = B vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyB(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 3, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 0, vect_out);

      GetRow(vect_in, 4, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 1, vect_out);

 
    } 
    
    //! application de la matrice C, vect_out = C vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyC(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 3, ligne_temp);
      SetRow(piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 3, vect_out);

      GetRow(vect_in, 4, ligne_temp);
      SetRow(piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 4, vect_out);

    } 
        
    //! application de la transposee de la matrice B, vect_out = B^T vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyD(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 3, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 4, vect_out);

 
    } 
    
    //! Calcul de vect_out = U(Q)
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeU( TinyMatrix<Real_wp, General, dimension, nb_quad> & Q,
	      TinyVector<Real_wp, nb_quad> & vect_out,
              PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    { 
      Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
      for (int g = 0; g < nb_quad; g++)
	{
	  Real_wp racine = sqrt(Q(0,g) *Q(0,g)  + (Q(1,g))*(Q(1,g))+(1+Q(2,g))*(1+Q(2,g)));	
	  vect_out(g) = coef*(Q(0,g) *Q(0,g) *0.5 + Q(1,g) *Q(1,g) *0.5 + (1+Q(2,g))- racine);

	}  

    }
    
    //! calcul de vect_out = grad(U) (q) a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,TinyVector<Real_wp, nb_quad> & Qminus,
		    TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
		    TinyVector<Real_wp, nb_quad> & vect_out,
                    PianoString<nb_base, nb_quad, TypeEquation> & piano_string, int k)
    {
      vect_out.Fill(0.0);
      if (k==0)
	WaveEquationStiffNL2TStringLoop<nb_quad-1>::ComputeDeltakUk0(Qplus, Qminus, Qsigma,
                                                                   vect_out, piano_string);

      
      if (k==1)
	WaveEquationStiffNL2TStringLoop<nb_quad-1>::ComputeDeltakUk1(Qplus, Qminus, Qsigma,
                                                                   vect_out, piano_string);

	
	  if (k==2)
	WaveEquationStiffNL2TStringLoop<nb_quad-1>::ComputeDeltakUk2(Qplus, Qminus, Qsigma,
                                                                   vect_out, piano_string);		

    }
    
    //! calcul de vect_out = d/dq( grad(U)(q))  a partir de Q^n-1, Q^n+1 et Q^n    
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDerivativeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,
                              TinyVector<Real_wp, nb_quad> & Qminus,
			      TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
			      TinyVector<Real_wp, nb_quad> & vect_out,
			      PianoString<nb_base, nb_quad, TypeEquation> & piano_string,
                              int k1, int k2)
    {
      vect_out.Zero();
      if (k1 == 0)
	{
	  if(k2 == 0)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	      Real_wp s1 = sqrt(Qplus (g)*Qplus (g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
	      Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( 0.5 - 1.0/(s1+s2)
                                        +  (Qplus(g)+Qminus(g))*Qplus(g)/(s1*(s1+s2)*(s1+s2)) );

                                        
		}
			    //~ cout << "what is DerivativedeltakU00?  " << vect_out << endl;
	    }

	  else if (k2 == 1)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	      Real_wp s1 = sqrt(Qplus (g)*Qplus (g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
	      Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (Qplus(g)+Qminus(g))*(Qsigma(0,g))*(1.0/s1+1.0/s2)
                                        /( (s1+s2)*(s1+s2) ) );
		}
			    //~ cout << "what is DerivativedeltakU01?  " << vect_out << endl;
	    }

	  else if (k2 == 2)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	      Real_wp s1 = sqrt(Qplus (g)*Qplus (g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
	      Real_wp s2 = sqrt(Qminus(g)*Qminus(g) + (Qsigma(0,g))*(Qsigma(0,g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (Qplus(g)+Qminus(g))*(1.0+Qsigma(1,g))
                                        *(1.0/s1+1.0/s2)/( (s1+s2)*(s1+s2) ) );
		}
					    //~ cout << "what is DerivativedeltakU02?  " << vect_out << endl;
	    }
	}      
      else if (k1 == 1)
	{
	  if (k2 == 0)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qplus (g))*(Qplus (g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qminus(g))*(Qminus(g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (Qplus(g)+Qminus(g))*(Qsigma(0,g))*(1.0/s1+1.0/s2)
                                        /( (s1+s2)*(s1+s2) ) );
		}
			    //~ cout << "what is DerivativedeltakU10?  " << vect_out << endl;
	    }
	  else if (k2 ==1)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qplus (g))*(Qplus (g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qminus(g))*(Qminus(g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( 0.5 - 1.0/(s1+s2)
                                        +  (Qplus(g)+Qminus(g))*Qplus(g)/(s1*(s1+s2)*(s1+s2)) );
		}
			    //~ cout << "what is DerivativedeltakU11?  " << vect_out << endl;
	    }
	  else if (k2 == 2)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qplus (g))*(Qplus (g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qminus(g))*(Qminus(g)) + (1+Qsigma(1,g))*(1+Qsigma(1,g)) );
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (Qplus(g)+Qminus(g))*(1.0+Qsigma(1,g))
                                        *(1.0/s1+1.0/s2)/( (s1+s2)*(s1+s2) ) );
		}
	    //~ cout << "what is DerivativedeltakU12?  " << vect_out << endl;
	    }

	}      
      else if (k1 == 2)
	{
	  if (k2 == 0)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qplus (g))*(1+Qplus (g)));
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qminus(g))*(1+Qminus(g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (2.0+Qplus(g)+Qminus(g))*Qsigma(0,g)*(1.0/s1+1.0/s2)
                                        /( (s1+s2)*(s1+s2) ) );
		}
			    //~ cout << "what is DerivativedeltakU20?  " << vect_out << endl;
	    }
	  else if (k2 ==1)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qplus (g))*(1+Qplus (g)));
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qminus(g))*(1+Qminus(g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef* ( (2.0+Qplus(g)+Qminus(g))*Qsigma(1,g)*(1.0/s1+1.0/s2)
                                        /( (s1+s2)*(s1+s2) ) );
		}
			    //~ cout << "what is DerivativedeltakU21?  " << vect_out << endl;
	    }
	  else if (k2 == 2)
	    {
	      for (int g = 0; g < nb_quad; g++)
		{
		  // A DEROULER
	      Real_wp s1 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qplus (g))*(1+Qplus (g)));
	      Real_wp s2 = sqrt( Qsigma(0,g)*Qsigma(0,g) + (Qsigma(1,g))*(Qsigma(1,g)) + (1+Qminus(g))*(1+Qminus(g)));
		  Real_wp coef = (piano_string.E*piano_string.A - piano_string.T0);
		  vect_out(g) = coef * (-1.0/(s1+s2) +  (2.0+Qplus(g)+Qminus(g))
                                        *(1.0+Qplus(g))/(s1*(s1+s2)*(s1+s2)) );
		}
			    //~ cout << "what is DerivativedeltakU22?  " << vect_out << endl;
	    }

	}            
    }	


    //! application of the scheme    
    template<int nb_base> static inline void
    ApplyMatrixEquation(TinyMatrix<TinyMatrix<Real_wp, General, nb_base, nb_base>,
                        General, dimension, dimension>& mat_Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& res_Un)
    {
      Mlt(mat_Un(0, 0), Un(0), res_Un(0));  
      MltAdd(mat_Un(0, 3), Un(3), res_Un(0)); 
        
      Mlt(mat_Un(1, 1), Un(1), res_Un(1));
      MltAdd(mat_Un(1, 4), Un(4), res_Un(1));
      
	  Mlt(mat_Un(2, 2), Un(2), res_Un(2));
	  
	  Mlt(mat_Un(3, 0), Un(0), res_Un(3));
	  MltAdd(mat_Un(3, 3), Un(3), res_Un(3));	  
	  
	  Mlt(mat_Un(4, 1), Un(1), res_Un(4));
	  MltAdd(mat_Un(4, 4), Un(4), res_Un(4));	  

    }
    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {
      WaveEquationStiffNL2TStringLoop<nb_quad-1>::ApplyNonLinearTerm(coefNL, dTm, nablaU);
    }
 
  };

class WaveEquationStringStiffLIN2T : public GenericEquationString
  {
  public :
    enum{dimension = 5}; // connu a la compil, super pour templater
    enum { type_equation = WAVE_NL_2T_TIMO};
    static const bool linear = true;
    
    TinyVector<bool,dimension> Dirichlet_Agraffe;
    TinyVector<bool,dimension> Dirichlet_Bridge;
    
    bool only_one_LM;
    
    //! default constructor
    WaveEquationStringStiffLIN2T()
    {
      Dirichlet_Agraffe(0) = true;
      Dirichlet_Agraffe(1) = true;
      Dirichlet_Agraffe(2) = true;
      Dirichlet_Agraffe(3) = false;
      Dirichlet_Agraffe(4) = false;
      Dirichlet_Bridge(0) = true;
      Dirichlet_Bridge(1) = true;
      Dirichlet_Bridge(2) = true;
      Dirichlet_Bridge(3) = false;
      Dirichlet_Bridge(4) = false;
      
      only_one_LM = false;
    }

    //! produit scalaire pour separer energie u et energie v
    template<class T, int m>
    static void DotProdColUV(const TinyMatrix<T, General, 5, m>& A,
			     const TinyMatrix<T, General, 5, m>& B,
			     TinyVector<T, m>& res_u, TinyVector<T, m>& res_v)
    {
      for (int p = 0; p < m; p++)
	{
	  res_u(p) += A(0, p)*B(0, p) + A(1, p)*B(1, p) + A(3, p)*B(3, p) + A(4, p)*B(4, p);
	  res_v(p) += A(2, p)*B(2, p);
	}
    }
        
    //! application de la matrice de masse, vect_out = M vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyM(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.rho_detuning*piano_string.rho*piano_string.A*ligne_temp, 1, vect_out);
      
	  GetRow(vect_in, 2, ligne_temp);
      SetRow(piano_string.rho*piano_string.A*ligne_temp, 2, vect_out);
      
      GetRow(vect_in, 3, ligne_temp);
      SetRow(piano_string.rho*piano_string.I*ligne_temp, 3, vect_out);
      
	  GetRow(vect_in, 4, ligne_temp);
      SetRow(piano_string.rho*piano_string.I*ligne_temp, 4, vect_out);   
    } 
    
    //! application de la matrice d'amortissement, vect_out = R vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyR(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow(piano_string.rho_detuning*2.0*piano_string.rho*piano_string.A*piano_string.amo(0)*ligne_temp, 1, vect_out);
      
      GetRow(vect_in, 2, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.A*piano_string.amo(1)*ligne_temp, 2, vect_out);
      
      GetRow(vect_in, 3, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.I*piano_string.amo(2)*ligne_temp, 3, vect_out);
      
      GetRow(vect_in, 4, ligne_temp);
      SetRow(2.0*piano_string.rho*piano_string.I*piano_string.amo(2)*ligne_temp, 4, vect_out);
          }
    
    //! application de l'amortissement quadratique, vect_out = H vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAdamping(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
                  TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
                  PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow(2.0*piano_string.T0*piano_string.amoB2(0)*ligne_temp, 0, vect_out);
      
      GetRow(vect_in, 1, ligne_temp);
      SetRow(2.0*piano_string.T0*piano_string.amoB2(0)*ligne_temp, 1, vect_out);  
          
      GetRow(vect_in, 2, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.A*piano_string.amoB2(1)*ligne_temp, 2, vect_out);

      GetRow(vect_in, 3, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.I*piano_string.amoB2(2)*ligne_temp, 3, vect_out);

      GetRow(vect_in, 4, ligne_temp);
      SetRow(2.0*piano_string.E*piano_string.I*piano_string.amoB2(2)*ligne_temp, 4, vect_out);


    }
    
    //! application de la rigidite non-raide, vect_out = A^grad vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	       TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	       PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {      
      vect_out.Zero();
      
      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.T0)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow((piano_string.T0)*ligne_temp, 1, vect_out);
      

    } 
    
    //! application de la rigidite raide, vect_out = A^stiff vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyAStiffGrad(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
		    TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
		    PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();

      TinyVector<Real_wp,nb_quad> ligne_temp;
      
      GetRow(vect_in, 0, ligne_temp);
      SetRow((piano_string.A*piano_string.G*piano_string.k_prime)*ligne_temp, 0, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow((piano_string.A*piano_string.G*piano_string.k_prime)*ligne_temp, 1, vect_out);
      
            
      GetRow(vect_in, 2, ligne_temp);
      SetRow((piano_string.E*piano_string.A)*ligne_temp, 2, vect_out);
      
      GetRow(vect_in, 3, ligne_temp);
      SetRow((piano_string.E*piano_string.I)*ligne_temp, 3, vect_out);
      
      GetRow(vect_in, 4, ligne_temp);
      SetRow((piano_string.E*piano_string.I)*ligne_temp, 4, vect_out);
    } 
    
    //! application de la matrice B, vect_out = B vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyB(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 3, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 0, vect_out);

      GetRow(vect_in, 4, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 1, vect_out);

 
    } 
    
    //! application de la matrice C, vect_out = C vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyC(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Fill(0.0);
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 3, ligne_temp);
      SetRow(piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 3, vect_out);

      GetRow(vect_in, 4, ligne_temp);
      SetRow(piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 4, vect_out);

    } 
        
    //! application de la transposee de la matrice B, vect_out = B^T vect_in
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ApplyD(TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_in,
	   TinyMatrix<Real_wp, General, dimension, nb_quad> & vect_out,
	   PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
    {
      vect_out.Zero();
      
      TinyVector<Real_wp, nb_quad> ligne_temp;
      GetRow(vect_in, 0, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 3, vect_out);

      GetRow(vect_in, 1, ligne_temp);
      SetRow(-piano_string.A*piano_string.G*piano_string.k_prime*ligne_temp, 4, vect_out);

 
    } 
    
    //! Calcul de vect_out = U(Q)
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeU( TinyMatrix<Real_wp, General, dimension, nb_quad> & Q,
	      TinyVector<Real_wp, nb_quad> & vect_out,
              PianoString<nb_base, nb_quad, TypeEquation> & piano_string)
     {
      vect_out.Zero();
	 }
    
    //! calcul de vect_out = grad(U) (q) a partir de Q^n-1, Q^n+1 et Q^n
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,TinyVector<Real_wp, nb_quad> & Qminus,
		    TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
		    TinyVector<Real_wp, nb_quad> & vect_out,
                    PianoString<nb_base, nb_quad, TypeEquation> & piano_string, int k)
    {
      vect_out.Fill(0.0);
    }
    //! calcul de vect_out = d/dq( grad(U)(q))  a partir de Q^n-1, Q^n+1 et Q^n    
    template<int nb_base, int nb_quad, class TypeEquation> static inline void
    ComputeDerivativeDeltakU( TinyVector<Real_wp, nb_quad> & Qplus,
                              TinyVector<Real_wp, nb_quad> & Qminus,
			      TinyMatrix<Real_wp, General, dimension-1, nb_quad> & Qsigma,
			      TinyVector<Real_wp, nb_quad> & vect_out,
			      PianoString<nb_base, nb_quad, TypeEquation> & piano_string,
                              int k1, int k2)
     {
      vect_out.Zero();
	 }

    //! application of the scheme    
    template<int nb_base> static inline void
    ApplyMatrixEquation(TinyMatrix<TinyMatrix<Real_wp, General, nb_base, nb_base>,
                        General, dimension, dimension>& mat_Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& Un,
                        TinyVector<TinyVector<Real_wp, nb_base>, dimension>& res_Un)
    {
      Mlt(mat_Un(0, 0), Un(0), res_Un(0));  
      MltAdd(mat_Un(0, 3), Un(3), res_Un(0)); 
        
      Mlt(mat_Un(1, 1), Un(1), res_Un(1));
      MltAdd(mat_Un(1, 4), Un(4), res_Un(1));
      
	  Mlt(mat_Un(2, 2), Un(2), res_Un(2));
	  
	  Mlt(mat_Un(3, 0), Un(0), res_Un(3));
	  MltAdd(mat_Un(3, 3), Un(3), res_Un(3));	  
	  
	  Mlt(mat_Un(4, 1), Un(1), res_Un(4));
	  MltAdd(mat_Un(4, 4), Un(4), res_Un(4));	  

    }
    template<int dimension, int nb_quad> static inline void
    ApplyNonLinearTerm(const Real_wp& coefNL,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& dTm,
                       TinyMatrix<Real_wp, General, dimension, nb_quad>& nablaU)
    {
    }
 
  };



}
#define MONTJOIE_FILE_EQUATION_STRING_HXX
#endif
