#ifndef MONTJOIE_FILE_MULTI_STRING_INLINE_CXX

namespace Montjoie
{

  inline MultiString_Base::MultiString_Base(int d) { dimension_ = d; }
  
  inline int MultiString_Base::GetNbStrings() const { return nb_corde; }
  
  inline int MultiString_Base::GetNbUnknownsPerString() const { return dimension_; }

  //! retourne le numero du ddl j de la corde i et la composante k du deplacement
  inline int
  MultiString_Base::GetStringDofNumber(int i, int k, int j)
  {
    return (nb_dof_before(i) + j)*dimension_ + k; 
  }
  
  
  //! retourne le numero du ddl du multiplicateur de Lagrange i (direction nu)
  inline int
  MultiString_Base::GetLMDofNumber(int i)
  {
#ifdef SELDON_DEBUG_LEVEL_4
    if (i < nb_corde)
      {
#endif
        return offset_lagrange_mult(0) + i;
#ifdef SELDON_DEBUG_LEVEL_4
      }
    else
      {
	cout << "You are trying to access to the lagrange multiplier 0 "
             << "of the string " << i << " but i is larger than the number of strings " << nb_corde << endl;
        
	abort();
      }
#endif

  }

  
  //! retourne le numero du ddl du multiplicateur de Lagrange i (direction tau)
  inline int
  MultiString_Base::GetLMGDofNumber(int i)
  {
#ifdef SELDON_DEBUG_LEVEL_4
    if (i < nb_corde)
      {
#endif
        return offset_lagrange_mult(1) + i;
#ifdef SELDON_DEBUG_LEVEL_4
      }
    else
      {
	cout << "You are trying to access to the lagrange multiplier 1 "
             << "of the string " << i << " but i is larger than the number of strings " << nb_corde << endl;
        
	abort();
      }
#endif
  }


  //! retourne le numero du ddl du multiplicateur de Lagrange i (direction nu)
  inline int
  MultiString_Base::GetLagrangeM_DofNumber(int i, int p)
  {
#ifdef SELDON_DEBUG_LEVEL_4
    if ((i < nb_corde) && (p < nb_lagrange_mult))
      {
#endif
        return offset_lagrange_mult(p) + i;
#ifdef SELDON_DEBUG_LEVEL_4
      }
    else
      {
	cout << "You are trying to access to the lagrange multiplier " << p
             << " of the string " << i << " but i is larger than the number of strings " << nb_corde
             << " or p is larger that the number of multipliers " << nb_lagrange_mult << endl;
        
	abort();
      }
#endif
  }
  
  
  //! retourne le numero du ddl associe au deplacement du marteau
  inline int
  MultiString_Base::GetHammerDofNumber()
  {
    return offset_hammer;
  }


  //! retourne le numero du ddl associe à la variable auxiliaire du marteau associée à la corde i
  inline int
  MultiString_Base::GetHammerAuxDofNumber(int i)
  {
    //return offset_hammer + i;
    return offset_hammer_aux + i;
  }

  
  //! Retourne la valeur de val2(n) pour la corde i, direction k et point xj.
  inline Real_wp 
  MultiString_Base::
  GetValue(int i, int k, int j, int n, Matrix<Real_wp> & val2)
  {
    return val2(GetStringDofNumber(i, k, j), n); 
  }
  
  
  //! retourne la valeur de val2(n) pour le multiplicateur de Lagrange de la corde i (direction nu)
  inline Real_wp 
  MultiString_Base::GetLM(int i, int n, Matrix<Real_wp> & val2) 
  {
    return val2(GetLMDofNumber(i), n);
  }
  
  
  //! retourne la valeur de val2(n) pour le multiplicateur de Lagrange
  //! de la corde i (direction tau)
  inline Real_wp 
  MultiString_Base::GetLMG(int i, int n, Matrix<Real_wp> & val2) 
  {
    return val2(GetLMGDofNumber(i), n);
  }
  

  //! retourne la valeur de val2(n) pour le multiplicateur de Lagrange de la corde i (direction nu)
  inline Real_wp 
  MultiString_Base
  ::GetLagrangeM(int i, int p, int n, Matrix<Real_wp> & val2) 
  {
    return val2(GetLagrangeM_DofNumber(i, p), n);
  }

    
  //! retourne la valeur de val2(n) pour le marteau (inconnue xi)  
  inline Real_wp 
  MultiString_Base::GetHammer(int n, Matrix<Real_wp> & val2) 
  {
    return val2(GetHammerDofNumber(), n);
  }


  //! retourne la valeur de val2(n) pour la variable auxiliaire du marteau
  inline Real_wp 
  MultiString_Base::GetHammerAux(int i, int n, Matrix<Real_wp> & val2) 
  {
    return val2(GetHammerAuxDofNumber(i), n);
  }

  
  //! Actualise la valeur de val2(n) pour la corde i, direction k et point xj.
  inline void 
  MultiString_Base::SetValue(int i, int k, int j, int n, Real_wp valeur,
                                                        Matrix<Real_wp>& val2)  
  {
    val2(GetStringDofNumber(i, k, j), n) = valeur;
  } 
  

  //! Actualise la valeur de val(n) pour les multiplicateurs de Lagrange
  inline void 
  MultiString_Base::SetLM(int i, int n, Real_wp valeur,
                                                     Matrix<Real_wp>& val2)
  {
    val2(GetLMDofNumber(i), n) = valeur;
  }
  
  
  //! Actualise la valeur de val(n) pour le deplacement du marteau
  inline void 
  MultiString_Base
  ::SetHammer(int n, Real_wp valeur, Matrix<Real_wp>& val2)
  { 
    val2(GetHammerDofNumber(), n) = valeur;
  }


  //! Actualise la valeur de val(n) pour la variable auxiliaire du marteau
  inline void 
  MultiString_Base
  ::SetHammerAux(int i, int n, Real_wp valeur, Matrix<Real_wp>& val2)
  { 
    val2(GetHammerAuxDofNumber(i), n) = valeur;
  }
  
 
  //! Actualise la valeur de scheme pour la corde i, direction k et point xj. 
  inline void 
  MultiString_Base
  ::SetSchemeString(int i, int k, int j, Real_wp valeur, VectReal_wp& scheme) 
  {
    scheme(GetStringDofNumber(i, k, j)) = valeur; 
  }
  
  
  //! Actualise la valeur de scheme pour le multiplicateur de Lagrange de la corde i (direction nu)
  inline void 
  MultiString_Base
  ::SetSchemeLM(int i, Real_wp valeur, VectReal_wp& scheme)
  {
    scheme(GetLMDofNumber(i)) = valeur;
  }
  
  
  //! Actualise la valeur de scheme pour le multiplicateur de Lagrange de la corde i
  //! (direction tau)
  inline void 
  MultiString_Base
  ::SetSchemeLMG(int i, Real_wp valeur, VectReal_wp& scheme)
  {
    scheme(GetLMGDofNumber(i)) = valeur;
  }
  

  //! Actualise la valeur de scheme pour le deplacement du marteau
  inline void 
  MultiString_Base::SetSchemeHammer(Real_wp valeur, VectReal_wp& scheme)
  
  {
    scheme(GetHammerDofNumber()) = valeur;
  }


  //! Actualise la valeur de scheme pour la variable auxiliaire du marteau
  inline void 
  MultiString_Base::SetSchemeHammerAux(int i, Real_wp valeur, VectReal_wp& scheme)
  
  {
    scheme(GetHammerAuxDofNumber(i)) = valeur;
  }
  
  
  //! Ajoute a la valeur de scheme pour la corde i, direction k et point xj
  inline void 
  MultiString_Base
  ::AddSchemeString(int i, int k, int j, Real_wp valeur, VectReal_wp& scheme)
  {
    scheme(GetStringDofNumber(i, k, j)) += valeur;
  }
  
  
  //! Ajoute a la valeur de scheme pour le multiplicateur de Lagrange i
  inline void 
  MultiString_Base
  ::AddSchemeLM(int i, Real_wp valeur, VectReal_wp& scheme)
  {
    scheme(GetLMDofNumber(i)) += valeur;
  }
  
  
  //! Ajoute a la valeur de scheme pour le deplacement du marteau
  inline void 
  MultiString_Base::AddSchemeHammer(Real_wp valeur, VectReal_wp& scheme)
  {
    scheme(GetHammerDofNumber()) += valeur;
  }


  //! Ajoute a la valeur de scheme pour la variable auxiliaire du marteau
  inline void 
  MultiString_Base::AddSchemeHammerAux(int i, Real_wp valeur, VectReal_wp& scheme)
  {
    scheme(GetHammerAuxDofNumber(i)) += valeur;
  }
  
  
  //! Actualise la valeur de la differentielle entre la corde i1, position j1, composante k1
  //! et la corde i2, position j2, composante k2
  inline void 
  MultiString_Base::
  SetDiffStringString(int i1, int k1, int j1, Real_wp valeur,
		      int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.SetEntry(GetStringDofNumber(i1, k1, j1), GetStringDofNumber(i2, k2, j2), valeur);
  }
  
  
  //! Actualise la valeur de la differentielle entre le multiplicateur 
  //! de Lagrange de la corde i (direction nu)
  //! et la corde i2, position j2, composante k2
  inline void 
  MultiString_Base::
  SetDiffLMString(int i, Real_wp valeur, int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.SetEntry(GetLMDofNumber(i), GetStringDofNumber(i2, k2, j2), valeur);	
  }
  
  
  //! Actualise la valeur de la differentielle entre le multiplicateur 
  //! de Lagrange de la corde i (direction tau)
  //! et la corde i2, position j2, composante k2  
  inline void 
  MultiString_Base::
  SetDiffLMGString(int i, Real_wp valeur, int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.SetEntry(GetLMGDofNumber(i), GetStringDofNumber(i2, k2, j2), valeur);	
  }
  

  //! Actualise la valeur de la differentielle entre la corde i2, position j2, composante k2  
  //! et le multiplicateur de Lagrange de la corde i (direction nu)
  inline void 
  MultiString_Base::
  SetDiffStringLM(int i2, int k2, int j2, Real_wp valeur, int i, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.SetEntry(GetStringDofNumber(i2, k2, j2), GetLMDofNumber(i), valeur);
  }
  

  //! Actualise la valeur de la differentielle entre la corde i2, position j2, composante k2  
  //! et le multiplicateur de Lagrange de la corde i (direction tau) 
  inline void 
  MultiString_Base
  ::SetDiffStringLMG(int i2, int k2, int j2, Real_wp valeur, int i, VirtualMatrix<Real_wp>& DiffMatrix)
  {    
    DiffMatrix.SetEntry(GetStringDofNumber(i2, k2, j2), GetLMGDofNumber(i), valeur);	    
  }
  
  
  //! Actualise la valeur de la differentielle entre la corde i, position j, composante k
  //! et le deplacement du marteau
  inline void 
  MultiString_Base
  ::SetDiffStringHammer(int i, int k, int j, Real_wp valeur, VirtualMatrix<Real_wp>& DiffMatrix)  
  {
    DiffMatrix.SetEntry(GetStringDofNumber(i, k, j), GetHammerDofNumber(), valeur);	
  }


  //! Actualise la valeur de la differentielle entre la corde i, position j, composante k
  //! et la variable auxiliaire du marteau pour la corde l
  inline void 
  MultiString_Base
  ::SetDiffStringHammerAux(int i, int k, int j, Real_wp valeur, int l, VirtualMatrix<Real_wp>& DiffMatrix)  
  {
    DiffMatrix.SetEntry(GetStringDofNumber(i, k, j), GetHammerAuxDofNumber(l), valeur); 
  }
  

  //! Actualise la valeur de la differentielle entre le deplacement du marteau  
  //! et la corde i, position j, composante k
  inline void 
  MultiString_Base
  ::SetDiffHammerString(Real_wp valeur, int i, int k, int j, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.SetEntry(GetHammerDofNumber(), GetStringDofNumber(i, k, j), valeur);	
  }


  //! Actualise la valeur de la differentielle entre la variable auxiliaire du marteau  
  //! et la corde i, position j, composante k
  inline void 
  MultiString_Base
  ::SetDiffHammerAuxString(int l, Real_wp valeur, int i, int k, int j, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.SetEntry(GetHammerAuxDofNumber(l), GetStringDofNumber(i, k, j), valeur); 
  }
  
  
  //! Actualise la valeur de la differentielle entre le multiplicateur de Lagrange de la corde i1,
  //! et le multiplicateur de la corde i2 (direction nu)
  inline void 
  MultiString_Base
  ::SetDiffLMLM(int i1, Real_wp valeur, int i2, VirtualMatrix<Real_wp>& DiffMatrix) 
  {  
    DiffMatrix.SetEntry(GetLMDofNumber(i1), GetLMDofNumber(i2), valeur);
  }
  

  //! Actualise la valeur de la differentielle entre le multiplicateur de Lagrange de la corde i1,
  //! et le multiplicateur de la corde i2 (direction tau)  
  inline  void 
  MultiString_Base
  ::SetDiffLMGLMG(int i1, Real_wp valeur, int i2, VirtualMatrix<Real_wp>& DiffMatrix) 
  {
    DiffMatrix.SetEntry(GetLMGDofNumber(i1), GetLMGDofNumber(i2), valeur);
  }
  
  
  //! Actualise la valeur de la differentielle entre le marteau et le marteau
  inline void 
  MultiString_Base
  ::SetDiffHammerHammer(Real_wp valeur, VirtualMatrix<Real_wp>& DiffMatrix) 
  {	
    DiffMatrix.SetEntry(GetHammerDofNumber(), GetHammerDofNumber(), valeur);
  }


  //! Actualise la valeur de la differentielle entre le marteau et la variable auxiliaire du marteau
  inline void 
  MultiString_Base
  ::SetDiffHammerHammerAux(Real_wp valeur, int i, VirtualMatrix<Real_wp>& DiffMatrix) 
  { 
    DiffMatrix.SetEntry(GetHammerDofNumber(), GetHammerAuxDofNumber(i), valeur);
  }


  //! Actualise la valeur de la differentielle entre la variable auxiliaire du marteau et le marteau
  inline void 
  MultiString_Base
  ::SetDiffHammerAuxHammer(int i, Real_wp valeur, VirtualMatrix<Real_wp>& DiffMatrix) 
  { 
    DiffMatrix.SetEntry(GetHammerAuxDofNumber(i), GetHammerDofNumber(), valeur);
  }


  //! Actualise la valeur de la differentielle entre la variable auxiliaire du marteau
  //! et la variable auxiliaire du marteau
  inline void 
  MultiString_Base
  ::SetDiffHammerAuxHammerAux(int i1, Real_wp valeur, int i2, VirtualMatrix<Real_wp>& DiffMatrix) 
  { 
    DiffMatrix.SetEntry(GetHammerAuxDofNumber(i1), GetHammerAuxDofNumber(i2), valeur);
  }
  
  
  //! Incremente la valeur de la differentielle entre la corde i1, position j1, composante k1
  //! et la corde i2, position j2, composante k2
  inline void 
  MultiString_Base
  ::AddDiffStringString(int i1, int k1, int j1, Real_wp valeur,
			int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.AddInteraction(GetStringDofNumber(i1, k1, j1),
                              GetStringDofNumber(i2, k2, j2), valeur);
  }
  
  
  //! Incremente la valeur de la differentielle entre le multiplicateur de Lagrange de la corde i
  //! (direction nu) et la corde i2, position j2, composante k2  
  inline void 
  MultiString_Base
  ::AddDiffLMString(int i, Real_wp valeur, int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.AddInteraction(GetLMDofNumber(i), GetStringDofNumber(i2, k2, j2), valeur);
  }
  
  
  //! Incremente la valeur de la differentielle entre le multiplicateur de Lagrange de la corde i
  //!  (direction tau) et la corde i2, position j2, composante k2    
   inline void 
  MultiString_Base
  ::AddDiffLMGString(int i, Real_wp valeur, int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.AddInteraction(GetLMGDofNumber(i), GetStringDofNumber(i2, k2, j2), valeur);
  }
  
  
  //! Incremente la valeur de la differentielle entre la corde i2, position j2, composante k2  
  //! et le multiplicateur de Lagrange de la corde i (direction nu)
   inline void 
  MultiString_Base
  ::AddDiffStringLM(int i2, int k2, int j2, Real_wp valeur, int i, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.AddInteraction(GetStringDofNumber(i2, k2, j2), GetLMDofNumber(i), valeur); 
  }
  
  
  //! Incremente la valeur de la differentielle entre la corde i2, position j2, composante k2  
  //! et le déplacement du marteau 
  inline void 
  MultiString_Base
  ::AddDiffStringHammer(int i, int k, int j, Real_wp valeur, VirtualMatrix<Real_wp>& DiffMatrix) 
  {
    DiffMatrix.AddInteraction(GetStringDofNumber(i, k, j), GetHammerDofNumber(), valeur); 
  }


  //! Incremente la valeur de la differentielle entre la corde i2, position j2, composante k2  
  //! et la variable auxiliaire du marteau 
  inline void 
  MultiString_Base
  ::AddDiffStringHammerAux(int i, int k, int j, Real_wp valeur, int l, VirtualMatrix<Real_wp>& DiffMatrix) 
  {
    DiffMatrix.AddInteraction(GetStringDofNumber(i, k, j), GetHammerAuxDofNumber(l), valeur); 
  }


  //! Incremente la valeur de la differentielle entre le déplacement du marteau
  //! et la corde i, position j, composante k
  inline void 
  MultiString_Base
  ::AddDiffHammerString(Real_wp valeur, int i, int k, int j, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.AddInteraction(GetHammerDofNumber(), GetStringDofNumber(i, k, j), valeur); 
  }


  //! Incremente la valeur de la differentielle entre la variable auxiliaire du marteau 
  //! et la corde i, position j, composante k
  inline void 
  MultiString_Base
  ::AddDiffHammerAuxString(int l, Real_wp valeur, int i, int k, int j, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.AddInteraction(GetHammerAuxDofNumber(l), GetStringDofNumber(i, k, j), valeur); 
  }
  
  
  //! Incremente la valeur de la differentielle entre le marteau et le marteau
  inline void MultiString_Base
  ::AddDiffHammerHammer(Real_wp valeur, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.AddInteraction(GetHammerDofNumber(), GetHammerDofNumber(), valeur); 
  }


  //! Incremente la valeur de la differentielle entre le marteau et la variable auxiliaire du marteau
  inline void MultiString_Base
  ::AddDiffHammerHammerAux(Real_wp valeur, int i, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.AddInteraction(GetHammerDofNumber(), GetHammerAuxDofNumber(i), valeur); 
  }


  //! Incremente la valeur de la differentielle entre la variable auxiliaire du marteau et le marteau
  inline void MultiString_Base
  ::AddDiffHammerAuxHammer(int i, Real_wp valeur, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.AddInteraction(GetHammerAuxDofNumber(i), GetHammerDofNumber(), valeur); 
  }


  //! Incremente la valeur de la differentielle entre la variable auxiliaire du marteau
  //! et la variable auxiliaire du marteau
  inline void MultiString_Base
  ::AddDiffHammerAuxHammerAux(int i1, Real_wp valeur, int i2, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.AddInteraction(GetHammerAuxDofNumber(i1), GetHammerAuxDofNumber(i2), valeur); 
  } 


  template<int nb_base, int nb_quad, class TypeEquation>
  inline bool MultiString<nb_base, nb_quad, TypeEquation>::IsAttachedString(int g)
  {
    return piano_strings(g).is_attached;
  }
  
  //~ template<int nb_base, int nb_quad, class TypeEquation>
  //~ inline Real_wp MultiString<nb_base, nb_quad, TypeEquation>::GetNuString(int g, int k2)
  //~ {
    //~ return piano_strings(g).NU(k2);
  //~ }
  //~ 
  //~ template<int nb_base, int nb_quad, class TypeEquation>
  //~ inline Real_wp MultiString<nb_base, nb_quad, TypeEquation>::GetNuOrthogonalString(int g, int k2)
  //~ {
    //~ return piano_strings(g).NUORTH(k2);
  //~ }
  //~ 
  //~ template<int nb_base, int nb_quad, class TypeEquation>
  //~ inline Real_wp MultiString<nb_base, nb_quad, TypeEquation>::GetNuOrthogonal2String(int g, int k2)
  //~ {
    //~ return piano_strings(g).NUORTH2(k2);
  //~ }
    
  template<int nb_base, int nb_quad, class TypeEquation>
  inline int MultiString<nb_base, nb_quad, TypeEquation>::GetNbStringDofs(int g)
  {
    return piano_strings(g).mesh.GetNbDof();
  }
  
  //get the type of equation
  template<int nb_base, int nb_quad, class TypeEquation>
  inline int MultiString<nb_base, nb_quad, TypeEquation>::GetStringType()
  {
    return piano_strings(0).var_eq.type_equation;
  }

	// extract transmission chevalet of row k2, col k3
	//made GetNuString, GetNuOrthogonalString, GetNuOrthogonal2String obselete
	//string number g, k1={NU,NUORTH,NUORTH2...} at dimension k2
  template<int nb_base, int nb_quad, class TypeEquation>
  inline Real_wp MultiString<nb_base, nb_quad, TypeEquation>::GetTransmissionChevaletString(int g, int k1, int k2)
  {
    return piano_strings(g).transmission_chevalet(k1,k2);
  }  
  
  //~ template<int nb_base, int nb_quad, class TypeEquation>
  //~ inline Real_wp MultiString<nb_base, nb_quad, TypeEquation>::GetNuOrthogonal2String(int g, int k2)
  //~ {
    //~ return piano_strings(g).NUORTH2(k2);
  //~ }
  
  
    template<int nb_base, int nb_quad, class TypeEquation>
  inline Real_wp MultiString<nb_base, nb_quad, TypeEquation>::GetBridgeHeight(int i)
  {
    return piano_strings(0).vect_bridge_height(i);
  }  
    


  //! Same as AddDiffLMString but you can choose the LM
   inline void 
  MultiString_Base
  ::AddDiffLaMuString(int n_LaMu, int i, Real_wp valeur, int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.AddInteraction(GetLaMuDofNumber(n_LaMu, i), GetStringDofNumber(i2, k2, j2), valeur);
  }
  
   //! AddDiff for LaMu of choice
   inline void 
  MultiString_Base
  ::AddDiffStringLaMu(int n_LaMu, int i2, int k2, int j2, Real_wp valeur, int i, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.AddInteraction(GetStringDofNumber(i2, k2, j2), GetLaMuDofNumber(n_LaMu, i), valeur); 
  }
  
    //! Same as AddSchemeLM but now you can define the LM of choice
  inline void 
  MultiString_Base
  ::AddSchemeLaMu(int n_LaMu, int i, Real_wp valeur, VectReal_wp& scheme)
  {
    scheme(GetLaMuDofNumber(n_LaMu, i)) += valeur;    
  }
  
    //! get Lagrange Multiplier of own choice
  inline Real_wp 
  MultiString_Base::GetLaMu(int n_LaMu, int i, int n, Matrix<Real_wp> & val2) 
  {
    return val2(GetLaMuDofNumber(n_LaMu, i), n);
  }
  
    //! get Lagrange multiplier Dof number of own choice
  inline int
  MultiString_Base::GetLaMuDofNumber(int n_LaMu, int i)
  {
#ifdef SELDON_DEBUG_LEVEL_4
    if (i < nb_corde)
      {
#endif
        return offset_lagrange_mult(n_LaMu) + i;
#ifdef SELDON_DEBUG_LEVEL_4
      }
    else
      {
	cout << "You are trying to access to the lagrange multiplier 0 "
             << "of the string " << i << " but i is larger than the number of strings " << nb_corde << endl;
        
	abort();
      }
#endif

  }
  
  //! Set diff matrix for the Langrange Multiplier that we want (in case there's more than one)
   inline void 
  MultiString_Base
  ::SetDiffLaMuLaMu(int n_LaMu, int i1, Real_wp valeur, int i2, VirtualMatrix<Real_wp>& DiffMatrix) 
  {  
    DiffMatrix.SetEntry(GetLaMuDofNumber(n_LaMu, i1), GetLaMuDofNumber(n_LaMu, i2), valeur);
  }
  
  //! Set Diff for LaMu of choice to string
  inline void 
  MultiString_Base::
  SetDiffLaMuString(int n_LaMu, int i, Real_wp valeur, int i2, int k2, int j2, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.SetEntry(GetLaMuDofNumber(n_LaMu, i), GetStringDofNumber(i2, k2, j2), valeur);	
  }
  

  
  //! Set diff value for desired LaMu
  inline void 
  MultiString_Base::
  SetDiffStringLaMu(int n_LaMu, int i2, int k2, int j2, Real_wp valeur, int i, VirtualMatrix<Real_wp>& DiffMatrix)
  {
    DiffMatrix.SetEntry(GetStringDofNumber(i2, k2, j2), GetLaMuDofNumber(n_LaMu, i), valeur);
  }


  
    //! Setscheme for Lagrange Multiplier of choice
  inline void 
  MultiString_Base
  ::SetSchemeLaMu(int n_LaMu, int i, Real_wp valeur, VectReal_wp& scheme)
  {
    scheme(GetLaMuDofNumber(n_LaMu, i)) = valeur;
  }
  
    //! Set LM of choice like SetLM
  inline void 
  MultiString_Base::SetLaMu(int n_LaMu, int i, int n, Real_wp valeur,
                                                     Matrix<Real_wp>& val2)
  {
    val2(GetLaMuDofNumber(n_LaMu, i), n) = valeur;
  }  
  
}

#define MONTJOIE_FILE_MULTI_STRING_INLINE_CXX
#endif
