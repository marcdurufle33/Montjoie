#ifndef MONTJOIE_FILE_COMPUTE_SCHEME_AND_DIFF_LM_CXX

namespace Montjoie
{
  
  /// *************
  // * Schema   *
  // ************
  
  //! Calcul du systeme non-lineaire F(X) = 0
  //! a resoudre pour calculer u^{n+1} en fonction de u^{n-1} et u^n
  /*!
    \param[inout] scheme resultat F(X)
    \param[in] val2 donnee X, en pratique val2(:, 0) = u^{n-1},
                    val2(:, 1) = u^{n-1}, val2(:, 2) = u^{n+1}
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  void MultiString<nb_base, nb_quad, TypeEquation>
  ::ComputeScheme(VectReal_wp & scheme, VectReal_wp& rhs, Matrix<Real_wp> & val2)
  {
    if (print_level >= 6)
      chrono_string.Start(VirtualTimer::SCHEME);
    
    // initializing scheme = rhs
    for (int i = 0; i < scheme.GetM(); i++)
      scheme(i) = rhs(i);
    
    VectReal_wp somme_corde_0(nb_corde);
    VectReal_wp somme_corde_2(nb_corde);
    somme_corde_0.Fill(0);
    somme_corde_2.Fill(0);
    
    // Loop : with openMP
    Real_wp sum = 0.0;
    Real_wp somme0;
    Real_wp somme2;
#ifndef MONTJOIE_WITH_MULTIPLE
#pragma omp parallel for schedule(static) private(somme2,somme0) reduction(+:sum)
#endif
    for (int i = 0; i < nb_corde; i++)
    {
      // Each string does its part and computes its interaction with the hammer
      // sum will contain the sum of the forces applied to the hammer (\sum F_i)
      sum += piano_strings(i).ComputeScheme(scheme, val2, *this);
      piano_strings(i).GetScalarProductAtBridge(val2, *this, somme0, somme2);
      somme_corde_0(i) = somme0;
      somme_corde_2(i) = somme2;  
    }
          
    // on stocke la force qui s'applique sur le marteau
    // pour ensuite la mettre dans un fichier de sortie
    somme_interaction_hammer = sum;
    
    // Contribution of the bridge
    bridge.AddScheme(*this, somme_corde_0, somme_corde_2, val2, scheme);
    
    // Contribution of the hammer with RHS coming from strings
    Real_wp coef_hammer = hammer.masse_marteau*invDeltat;
    if (shank.is_interacting_string)
    {    
      // contribution of the shank
      shank.t_courant = t_courant;
      shank.AddScheme(scheme, val2, this->GetHammerDofNumber(), somme_interaction_hammer,
                       DOSSIER + file_storage_param);
    }
    else
    {
      if (!drop_hammer)
      {
        // equation masse_marteau (xi^(n+1) + xi^(n-1) - 2 xi^n)/dt - dt (\sum F_i) = 0
        SetSchemeHammer((GetHammer(2,val2) + GetHammer(0,val2) - 2.0*GetHammer(1,val2))
                              *coef_hammer, scheme);
                          
        if (hammer.is_interacting)
        {
          AddSchemeHammer(-somme_interaction_hammer*Deltat, scheme);
        }
      }
    }
    if ((hammer.with_auxiliary_variable) && (!drop_hammer))
    {
      for (int i = 0; i < nb_corde; i++)
      {
        Real_wp ux0p  = 0.0;
        Real_wp ux0   = 0.0;
        Real_wp ux0m  = 0.0;
  
        // on calcule le deplacement moyen de la corde au niveau de l'impact du marteau
        for (int j = 0; j < piano_strings(i).repartition_marteau_sparse.GetM(); j++)
        {
          int jj = piano_strings(i).repartition_marteau_sparse.Index(j);
          Real_wp repj = piano_strings(i).repartition_marteau_sparse.Value(j);
            
          if (hammer.strike_at_an_angle)
          {
            ux0p += GetValue(i, 0, jj, 2, val2)*repj*cos(hammer.strike_angle*pi_wp/180.0)
                     + GetValue(i, 1, jj, 2, val2)*repj*sin(hammer.strike_angle*pi_wp/180.0);
            ux0  += GetValue(i, 0, jj, 1, val2)*repj*cos(hammer.strike_angle*pi_wp/180.0)
                     + GetValue(i, 1, jj, 1, val2)*repj*sin(hammer.strike_angle*pi_wp/180.0);
            ux0m += GetValue(i, 0, jj, 0, val2)*repj*cos(hammer.strike_angle*pi_wp/180.0)
                     + GetValue(i, 1, jj, 0, val2)*repj*sin(hammer.strike_angle*pi_wp/180.0); 
          }
          else
          {
            ux0p += GetValue(i, 0, jj, 2, val2)*repj;
            ux0  += GetValue(i, 0, jj, 1, val2)*repj;
            ux0m += GetValue(i, 0, jj, 0, val2)*repj; 
          }
        }
        SetSchemeHammerAux(i, GetHammerAux(i, 2, val2) - GetHammerAux(i, 1, val2)
                                  - 0.5*hammer.DeriveeGAux(ux0-GetHammer(1, val2), i)
                                  *((ux0p-GetHammer(2, val2)) - (ux0m-GetHammer(0, val2))), scheme);
      }
    }

    if (print_level >= 6)
      chrono_string.Stop(VirtualTimer::SCHEME);
  }

  
  //! forcing constraint equations to be true in the new solution x
  template<int nb_base, int nb_quad, class TypeEquation> void 
  MultiString<nb_base, nb_quad, TypeEquation>
  ::ForceConstraint(Matrix<Real_wp>& x)
  {
    if (drop_hammer)
      return;
    
    // on force d'abord la position du marteau
    if (shank.is_interacting_string)
      shank.ForceConstraint(x, this->GetHammerDofNumber());
    
    if (hammer.is_interacting)
    {
      for (int i = 0; i < nb_corde; i++)
      {
        // calcul de la force du marteau
        Real_wp interaction_hammer = this->ComputeHammerInteraction(piano_strings(i), x);
          
        // enforcing force_marteau = interaction_hammer
        x(this->offset_force_hammer+i, 2) = interaction_hammer;
      }
    }
  }
  
  
  //! computation of right hand side for linear scheme
  /*!
    scheme : part of F(X) which depends only on X^n and X^n-1
    val2 : val2(:, 2) is X, val2(:, 1) is X^n, val2(:, 0) is X^n-1
    newton_rhs : true if called by MegaNewtonSolver (during a Newton solve)
               false if called by Advance
   */
  template<int nb_base, int nb_quad, class TypeEquation> void 
  MultiString<nb_base, nb_quad, TypeEquation>
  ::ComputeRightHandSide(VectReal_wp & scheme, Matrix<Real_wp> & val2,
                        bool newton_rhs)
  {
    if (print_level >= 6)
      chrono_string.Start(VirtualTimer::SCHEME);
    
    scheme.Fill(Real_wp(0));
    Real_wp somme0;
    Real_wp somme2;

    VectReal_wp somme_corde_0(nb_corde);
    somme_corde_0.Fill(0);
    
    for (int i = 0; i < nb_corde; i++)
    {
      piano_strings(i).ComputeRightHandSide(scheme, val2, *this, newton_rhs);
      piano_strings(i).GetScalarProductAtBridge(val2, *this, somme0, somme2);
      somme_corde_0(i) = somme0;
    }
    
    // Contribution of the LM
    if (!newton_rhs)
      bridge.AddRightHandSide(*this, somme_corde_0, val2, scheme);
    
    if (print_level >= 6)
      chrono_string.Stop(VirtualTimer::SCHEME);
  }
  
  
  //! calcul de la jacobienne du systeme non-lineaire F(X) = 0 a resoudre
  /*!
    \param[inout] DiffMatrix matrice jacobienne
    \param[in] val2 donnee X, en pratique val2(:, 0) = u^{n-1},
                    val2(:, 1) = u^{n-1}, val2(:, 2) = u^{n+1}
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  void MultiString<nb_base, nb_quad, TypeEquation>
  ::ComputeDiff(VirtualMatrix<Real_wp> & DiffMatrix, Matrix<Real_wp> & val2)
  {
    if (print_level >= 6)
      chrono_string.Start(VirtualTimer::JACOBIAN);
    
    DiffMatrix.Zero();
    
    Real_wp interaction_hammer_deriv = 0.0;
    int num_xsi = this->GetHammerDofNumber();
    if (shank.is_interacting_string)
      num_xsi = shank.Fcoupl_dof_number+1;
    
    // Loop : to do with openMP
    // TOTALEMENT INUTILE DE REMPLIR UNE VARIABLE ICI NON ? ET LE NUM XSI SERT A RIEN !
    for (int i = 0; i < nb_corde; i++)
    {
      // Chaque corde remplit sa partie
      interaction_hammer_deriv += piano_strings(i).ComputeDiff(DiffMatrix, val2, num_xsi, *this);
    }

    bridge.AddDiff(*this, DiffMatrix, val2);
    
    if (shank.is_interacting_string)
    {    
      // contribution of the shank
      shank.AddDiff(DiffMatrix, val2, this->GetHammerDofNumber(), this->offset_force_hammer);
    }
    else
    {
      if (!drop_hammer)
      {
        AddDiffHammerHammer(hammer.masse_marteau*invDeltat, DiffMatrix);
      }
    }
    
    if (print_level >= 6)
      chrono_string.Stop(VirtualTimer::JACOBIAN);    
  }


  //! mise à jour de la jacobienne du systeme non-lineaire F(X) = 0 a resoudre
  /*!
    \param[inout] DiffMatrix matrice jacobienne
    \param[in] val2 donnee X, en pratique val2(:, 0) = u^{n-1},
                    val2(:, 1) = u^{n-1}, val2(:, 2) = u^{n+1}
   */
  template<int nb_base, int nb_quad, class TypeEquation>
  void MultiString<nb_base, nb_quad, TypeEquation>
  ::UpdateDiff(VirtualMatrix<Real_wp> & DiffMatrix, Matrix<Real_wp> & val2)
  {
    for (int i = 0; i < nb_corde; i++)
    {
      // Chaque corde met à jour sa partie
      piano_strings(i).UpdateDiff(DiffMatrix, val2, *this);
    }
  }
  
  
  //! computes linear matrix applied to Un and Unm1
  template<int nb_base, int nb_quad, class TypeEquation>
  void MultiString<nb_base, nb_quad, TypeEquation>::ConstructLinearMatrix()
  {
    //mat_Unm1.Reallocate(size, size);
    //mat_Un.Reallocate(size, size);
    //mat_Unm1.Fill(0);
    //mat_Un.Fill(0);

    for (int i = 0; i < nb_corde; i++)
    {
      piano_strings(i).ComputeLinearMatrix(*this);
    }
  }
  
}

#define MONTJOIE_FILE_COMPUTE_SCHEME_AND_DIFF_LM_CXX
#endif
