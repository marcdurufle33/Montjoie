#ifndef MONTJOIE_FILE_AXISYM_GALBRUN_CXX

#include "AxiSymGalbrun.hxx"

namespace Montjoie
{

  /***********************
   * VarGalbrunIndex_Axi *
   ***********************/


  //! additionnal input datas
  void VarGalbrunIndex_Axi
    ::SetInputData(const string& description_field, const VectString& parameters)
    {
      if (description_field == "DropUnstableTerms")
      {
        if (parameters(0) == "Convective")
        {
          drop_unstable_terms = DROP_CONVECTIVE;
          if (parameters.GetM() > 1)
            coef_convective_term = to_num<Real_wp>(parameters(1));
          else
            coef_convective_term = 0.0;
        }
        else if (parameters(0) == "NonUniform")
          drop_unstable_terms = DROP_NON_UNIFORM;
        else
          drop_unstable_terms = DROP_NONE;
      }
      else if (description_field == "ApplyConvectiveCorrectionSource")
      {
        if (parameters(0) == "YES")
          apply_convective_derivate_source = true;
        else
          apply_convective_derivate_source = false;
      }
      else if (description_field == "ComputeGravity") // AJOUT NATHAN
      {
        if (parameters(0) == "YES")
        {
          compute_gravity = true;
          std::cout << "Gravity will be computed" << std::endl;
        }
        else
          compute_gravity = false;
      }
      else
        VarAxisymProblem::SetInputData(description_field, parameters);

    }


  //! allocation of arrays containing physical indices
  void VarGalbrunIndex_Axi::InitIndices(int n)
  {
    ref_rho0.Reallocate(n+1);
    ref_sigma.Reallocate(n+1);
    ref_c0.Reallocate(n+1);
    ref_p0.Reallocate(n+1);
    ref_m0.Reallocate(n+1);
    //ref_g0.Reallocate(n+1); // AJOUT NATHAN
    for (int i = 0; i <= n; i++)
    {
      ref_rho0(i).SetIdentity();
      ref_sigma(i).Zero();
      ref_c0(i).SetIdentity();
      ref_p0(i).SetIdentity();
      ref_m0(i).Zero();
      //ref_g0(i).Zero(); // AJOUT NATHAN
    }
  }


  //! we list the indices that are varying
  void VarGalbrunIndex_Axi
    ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
        IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
        Vector<bool>& compute_hess)
    {
      int nb = 0;
      for (int i = 0; i < ref_rho0.GetM(); i++)
      {
        nb += ref_rho0(i).GetNbVaryingMedia();
        nb += ref_c0(i).GetNbVaryingMedia();
        nb += ref_p0(i).GetNbVaryingMedia();
        nb += ref_m0(i).GetNbVaryingMedia();
        nb += ref_sigma(i).GetNbVaryingMedia();
        //nb += ref_g0(i).GetNbVaryingMedia(); // AJOUT NATHAN
      }

      rho_real.Reallocate(nb);
      num_ref.Reallocate(nb);
      num_index.Reallocate(nb);
      num_component.Reallocate(nb);
      compute_grad.Reallocate(nb);
      compute_hess.Reallocate(nb);
      compute_grad.Fill(false);
      compute_hess.Fill(false);
      nb = 0;
      for (int i = 0; i < ref_rho0.GetM(); i++)
      {
        int nb0 = nb;
        ref_rho0(i).GetVaryingMedia(nb, rho_real, num_component);
        for (int j = nb0; j < nb; j++)
        {
          if ((type_model == GALBRUN) || (type_model == LEE_MODEL))
            compute_grad(j) = true;
          if(compute_gravity) // AJOUT NATHAN
            compute_grad(j) = true;

          num_index(j) = 0;
          num_ref(j) = i;
        }

        nb0 = nb;
        ref_c0(i).GetVaryingMedia(nb, rho_real, num_component);
        for (int j = nb0; j < nb; j++)
        {
          if (type_model == LEE_MODEL)
            compute_grad(j) = true;

          num_index(j) = 1;
          num_ref(j) = i;
        }

        nb0 = nb;
        ref_p0(i).GetVaryingMedia(nb, rho_real, num_component);
        for (int j = nb0; j < nb; j++)
        {
          num_index(j) = 2;
          num_ref(j) = i;
          if ((type_model == LEE_MODEL) || (type_model == GALBRUN) || (type_model == GALBRUN_DIV))
            compute_grad(j) = true;
          if(compute_gravity) // AJOUT NATHAN
          {
            compute_hess(j) = true;
            compute_grad(j) = true;
          }
        }

        nb0 = nb;
        ref_m0(i).GetVaryingMedia(nb, rho_real, num_component);
        for (int j = nb0; j < nb; j++)
        {
          num_index(j) = 3;
          num_ref(j) = i;
          if ((type_model == LEE_MODEL) || (type_model == SIMPLIFIED_LEE) || (type_model == GALBRUN))
            compute_grad(j) = true;
          if(compute_gravity) //AJOUT NATHAN
          {
            compute_hess(j) = true;
            compute_grad(j) = true;
          }
        }

        nb0 = nb;
        ref_sigma(i).GetVaryingMedia(nb, rho_real, num_component);
        for (int j = nb0; j < nb; j++)
        {
          if (type_model == GALBRUN)
            compute_grad(j) = true;

          num_index(j) = 4;
          num_ref(j) = i;
        }

        // AJOUT NATHAN
//        nb0 = nb;
//        ref_g0(i).GetVaryingMedia(nb, rho_real, num_component);
//        for(int j=nb0; j<nb; j++)
//        {
//          num_index(j) = 5;
//          num_ref(j) = i;
//          //if(type_model == GALBRUN)
//            compute_grad(j) = true;
//        }
      }
    }


  //! computes physical properties on quadrature points
  void VarGalbrunIndex_Axi::ComputePhysicalCoefficients()
  {
    int nb_elt = var_problem.mesh.GetNbElt();
    eval_flow.Reallocate(nb_elt);
    eval_rho.Reallocate(nb_elt);
    eval_c0.Reallocate(nb_elt);
    eval_sigma.Reallocate(nb_elt);
    if ((type_model == LEE_MODEL) || (type_model == SIMPLIFIED_LEE) || (type_model == GALBRUN))
      grad_flow.Reallocate(nb_elt);
    if (compute_gravity) // AJOUT NATHAN
    {
      eval_gravity.Reallocate(nb_elt);
      grad_gravity.Reallocate(nb_elt);
    }

    bool compute_grad_gamma = false;
    if (type_model == LEE_MODEL)
    {
      compute_grad_gamma = true;
      div_flow.Reallocate(nb_elt);
      eval_gamma.Reallocate(nb_elt);
      grad_gamma.Reallocate(nb_elt);
      grad_rho.Reallocate(nb_elt);
    }

    if ((type_model == LEE_MODEL) || (type_model == GALBRUN) || (type_model == GALBRUN_DIV))
      grad_p0.Reallocate(nb_elt);

    if (type_model == GALBRUN)
    {
      grad_sigma.Reallocate(nb_elt);
      grad_rho.Reallocate(nb_elt);
    }

    bool compute_grad_rho = false, compute_grad_c0 = false;
    if (grad_rho.GetM() > 0)
      compute_grad_rho = true;

    Real_wp c, p0, r, rho0;
    R2 nabla_rho0, nabla_c0;
    for (int i = 0; i < nb_elt; i++)
    {
      int ref = var_problem.mesh.Element(i).GetReference();
      int N = var_problem.Glob_PointsQuadrature(i).GetM();
      eval_flow(i).Reallocate(N);
      eval_rho(i).Reallocate(N);
      if (grad_flow.GetM() > 0)
        grad_flow(i).Reallocate(N);

      eval_sigma(i).Reallocate(N);
      eval_c0(i).Reallocate(N);

      if (div_flow.GetM() > 0)
        div_flow(i).Reallocate(N);

      if (eval_gamma.GetM() > 0)
        eval_gamma(i).Reallocate(N);

      if (grad_p0.GetM() > 0)
        grad_p0(i).Reallocate(N);

      if (grad_sigma.GetM() > 0)
        grad_sigma(i).Reallocate(N);

      if (grad_rho.GetM() > 0)
        grad_rho(i).Reallocate(N);

      if (compute_grad_gamma)
        grad_gamma(i).Reallocate(N);

      if (compute_gravity) // AJOUT NATHAN
      {
        eval_gravity(i).Reallocate(N);
        grad_gravity(i).Reallocate(N); 
      }

      for (int j = 0; j < N; j++)
      {
        if (grad_flow.GetM() > 0)
          this->ref_m0(ref).GetCoefGradient(var_problem, i, j, 
              eval_flow(i)(j), grad_flow(i)(j));
        else
          eval_flow(i)(j) = this->ref_m0(ref).GetCoefficient(var_problem, i, j); 

        if (grad_p0.GetM() > 0)
        {
          this->ref_p0(ref).GetCoefGradient(var_problem, i, j, p0, grad_p0(i)(j));
        }
        else
          p0 = this->ref_p0(ref).GetCoefficient(var_problem, i, j);

        r = var_problem.Glob_PointsQuadrature(i)(j)(0);
        if (div_flow.GetM() > 0)
          div_flow(i)(j) = grad_flow(i)(j)(0, 0) + grad_flow(i)(j)(2, 1) + eval_flow(i)(j)(0)/r;

        if (compute_grad_rho || compute_grad_gamma)
          this->ref_rho0(ref).GetCoefGradient(var_problem, i, j, rho0, nabla_rho0);
        else
          rho0 = this->ref_rho0(ref).GetCoefficient(var_problem, i, j);

        eval_rho(i)(j) = rho0;
        if (grad_rho.GetM() > 0)
          grad_rho(i)(j) = nabla_rho0;

        if (grad_sigma.GetM() > 0)
          this->ref_sigma(ref).GetCoefGradient(var_problem, i, j, eval_sigma(i)(j), grad_sigma(i)(j));
        else
          eval_sigma(i)(j) = this->ref_sigma(ref).GetCoefficient(var_problem, i, j);

        if (compute_grad_c0 || compute_grad_gamma)
          this->ref_c0(ref).GetCoefGradient(var_problem, i, j, c, nabla_c0);
        else
          c = this->ref_c0(ref).GetCoefficient(var_problem, i, j);

        eval_c0(i)(j) = c;

        if (eval_gamma.GetM() > 0)
          eval_gamma(i)(j) = c*c*eval_rho(i)(j)/p0;

        if (grad_gamma.GetM() > 0)
          grad_gamma(i)(j) = 2.0*c*rho0/p0*nabla_c0 + c*c/p0*nabla_rho0 - c*c*rho0/(p0*p0)*grad_p0(i)(j);

        // AJOUT NATHAN
        if (compute_gravity)
        {
          R3 g;
          TinyMatrix<Real_wp, General, 3, 2> grad_g;

          Real_wp rho;
          R2 grad_rho;
          this->ref_rho0(ref).GetCoefGradient(var_problem,i,j,rho,grad_rho);

          R3 v;
          TinyMatrix<Real_wp, General, 3, 2> grad_v;
          TinyVector<TinyMatrix<Real_wp,Symmetric,2,2>,3> hess_v;
          this->ref_m0(ref).GetCoefHessian(var_problem,i,j,v,grad_v,hess_v);

          Real_wp p;
          R2 grad_p;
          TinyMatrix<Real_wp,Symmetric,2,2> hess_p;
          this->ref_p0(ref).GetCoefHessian(var_problem,i,j,p,grad_p,hess_p);
          // Add (\nabla p)/rho
          g(0) = grad_p(0)/rho;
          g(1) = 0.0;
          g(2) = grad_p(1)/rho;

          // Evaluate v\cdot\nabla v /!\ material derivative in polar coordinates
          g(0) += v(0)*grad_v(0,0)+v(2)*grad_v(0,1) -v(1)*v(1)/r;
          g(1) += v(0)*grad_v(1,0)+v(2)*grad_v(1,1) +v(1)*v(0)/r;
          g(2) += v(0)*grad_v(2,0)+v(2)*grad_v(2,1);


          // Compute the gradient of the gravity, needed for galbrun
          for(int l=0; l<3; l++)
          {
            // d_r g_l
            grad_g(l,0) = v(0)*hess_v(l)(0,0)+v(2)*hess_v(l)(1,0);
            grad_g(l,0) += grad_v(0,0)*grad_v(l,0)+grad_v(2,0)*grad_v(l,1);

            // d_z g_l
            grad_g(l,1) = v(0)*hess_v(l)(1,0)+v(2)*hess_v(l)(1,1);
            grad_g(l,1) += grad_v(0,1)*grad_v(l,0)+grad_v(2,1)*grad_v(l,1);
          }

          grad_g(0,0) += hess_p(0,0)/rho-grad_rho(0)*grad_p(0)/(rho*rho);
          grad_g(0,1) += hess_p(1,0)/rho-grad_rho(1)*grad_p(0)/(rho*rho);
          grad_g(2,0) += hess_p(0,1)/rho-grad_rho(0)*grad_p(1)/(rho*rho);
          grad_g(2,1) += hess_p(1,1)/rho-grad_rho(1)*grad_p(1)/(rho*rho);


          grad_g(0,0) -= 2.0*v(1)*grad_v(1,0)/r+v(1)*v(1)/(r*r);
          grad_g(0,1) -= 2.0*v(1)*grad_v(1,1)/r;
          grad_g(1,0) += v(0)*grad_v(1,0)/r+grad_v(0,0)*v(1)/r-v(0)*v(1)/(r*r);
          grad_g(1,1) += v(0)*grad_v(1,1)/r+grad_v(0,1)*v(1)/r;

          // save the solution
          eval_gravity(i)(j) = g;
          grad_gravity(i)(j) = grad_g;

//          this->ref_g0(ref).GetCoefGradient(var_problem,i,j,g,grad_g);
//          Real_wp rc = var_problem.Glob_PointsQuadrature(i)(j)(0);
//          Real_wp z = var_problem.Glob_PointsQuadrature(i)(j)(1);
//          Real_wp r = sqrt(rc*rc+z*z);
//          eval_gravity(i)(j)(0) = (rc/r)*g(0)+(z/r)*g(1);
//          eval_gravity(i)(j)(2) = (z/r)*g(0)-(rc/r)*g(1);
//          eval_gravity(i)(j)(1) = g(2);
//          grad_gravity(i)(j)(0,0) = (1/r-(rc*rc)/(r*r*r))*g(0)+(rc/r)*grad_g(0,0)-(rc*z)/(r*r*r)*g(1)+(z/r)*grad_g(1,0);
//          grad_gravity(i)(j)(0,1) = -(rc*z)/(r*r*r)*g(0)+(rc/r)*grad_g(0,1)-(1/z-(z*z)/(r*r*r))*g(1)+(z/r)*grad_g(1,1);
//          grad_gravity(i)(j)(1,0) = grad_g(2,0);
//          grad_gravity(i)(j)(1,1) = grad_g(2,1);
//          grad_gravity(i)(j)(2,0) = -(rc*z)/(r*r*r)*g(0)+(z/r)*grad_g(0,0)-(1/r-(rc*rc)/(r*r*r))*g(1)-(rc/r)*grad_g(1,0);
//          grad_gravity(i)(j)(2,1) = (1/r-(z*z)/(r*r*r))*g(0)+(z/r)*grad_g(0,1)+(rc*z)/(r*r*r)*g(1)-(rc/r)*grad_g(1,1);

        }
        if (rho0 == Real_wp(0))
        {
          cout << "specify a non-null density " << endl;
          abort();
        }

        if (p0 == Real_wp(0))
        {
          cout << "specify a non-null pressure " << endl;
          abort();
        }            

        if (c == Real_wp(0))
        {
          cout << "specify a non-null speed of sound " << endl;
          abort();
        }            
      }
    }

    for (int ref = 0; ref < ref_m0.GetM(); ref++)
    {
      this->ref_m0(ref).Clear();
      this->ref_p0(ref).Clear();
      this->ref_c0(ref).Clear();
      this->ref_rho0(ref).Clear();
      this->ref_sigma(ref).Clear();
      //this->ref_g0(ref).Clear(); // AJOUT NATHAN
    }

    //AJOUT NATHAN
    //Write gravity on file
    if(compute_gravity)
    {
      int N = var_problem.GetMeshNumbering(0).GetNbDof();
      std::cout << "Writing the gravity..." << std::endl;
      VectReal_wp f(3*N);
      f.Zero();
      for( int i=0; i<var_problem.mesh.GetNbElt(); i++)
      {
        const ElementReference<Dimension2,1>& Fb = var_problem.GetReferenceElementH1(i);
        IVect Nodle = var_problem.GetDofNumberOnElement(i);
        for(int j=0; j<Fb.GetNbDof(); j++)
        {
          int numDof = Nodle(j);
          f(numDof) = eval_gravity(i)(j)(0);
          f(numDof+N) = eval_gravity(i)(j)(1);
          f(numDof+2*N) = eval_gravity(i)(j)(2);
        }
      }
      var_problem.GetOutputProblem().WriteOutputFile(f,"grav",3);
      std::cout << "Done" << std::endl;
    }
  }


  //! sets physical indices of media with reference i
  void VarGalbrunIndex_Axi
    ::SetIndices(int i, const VectString& parameters)
    {
      if (i >= ref_rho0.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_rho0.GetM() << endl;
        abort();
      }

      int nb = 1;
      ref_m0(i).SetInputData(nb, parameters, parameters(0));
      ref_sigma(i).SetInputData(nb, parameters, parameters(0));
      ref_rho0(i).SetInputData(nb, parameters, parameters(0));
      ref_c0(i).SetInputData(nb, parameters, parameters(0));
      ref_p0(i).SetInputData(nb, parameters, parameters(0));
      //ref_g0(i).SetInputData(nb, parameters, parameters(0)); // AJOUT NATHAN
    }


  //! reading of a physical index
  /*!
    \param[in] name_media name of the physical index
    \param[in] i physical domain domain
    \param[in] parameters parameters of the matching line of the data file
    the data file contains a line like PhysicalMedia = ...
    */
  void VarGalbrunIndex_Axi
    ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
    {
      int nb = 1;
      if (i >= ref_rho0.GetM())
      {
        cout << "Not enough indices stored, call InitIndices with a higher N" << endl;
        cout << "Current reference i : " << i << " < " << ref_rho0.GetM() << endl;
        abort();
      }

      if (name_media == "M")
        ref_m0(i).SetInputData(nb, parameters, parameters(0));
      else if (name_media == "sigma")
        ref_sigma(i).SetInputData(nb, parameters, parameters(0));
      else if (name_media == "rho0")
        ref_rho0(i).SetInputData(nb, parameters, parameters(0));
      else if (name_media == "c0")
        ref_c0(i).SetInputData(nb, parameters, parameters(0));
      else if (name_media == "p0")
        ref_p0(i).SetInputData(nb, parameters, parameters(0));
      //else if (name_media == "g0") // AJOUT NATHAN
        //ref_g0(i).SetInputData(nb, parameters, parameters(0));
      else
      {
        cout << "Unknown media : " << name_media << endl;
        abort();
      }
    }


  //! returns the name associated with the physical index m
  string VarGalbrunIndex_Axi::GetPhysicalIndexName(int m) const
  {
    switch (m)
    {
      case 3: return string("M");
      case 4: return string("sigma");
      case 0: return string("rho0");
      case 1: return string("c0");
      case 2: return string("p0");
      case 5: return string("g0"); // AJOUT NATHAN
    }

    return string();
  }


  //! returns the number of physical indices
  int VarGalbrunIndex_Axi::GetNbPhysicalIndices() const
  {
    return ref_rho0.GetM();
  }


  //! returns true if the domain i contains a varying media
  bool VarGalbrunIndex_Axi::IsVaryingMedia(int i) const
  {
    if (ref_rho0(i).IsVarying())
      return true;

    if (ref_c0(i).IsVarying())
      return true;

    if (ref_p0(i).IsVarying())
      return true;

    if (ref_m0(i).IsVarying())
      return true;

    if (ref_sigma(i).IsVarying())
      return true;

    //if (ref_g0(i).IsVarying()) // AJOUT NATHAN
     // return true;

    return false;
  }


  //! returns velocity of media ref
  Real_wp VarGalbrunIndex_Axi::GetVelocityOfMedia(int ref) const
  {
    return ref_c0(ref).GetConstant();
  }


  //! returns velocity at infinity
  Real_wp VarGalbrunIndex_Axi::GetVelocityOfInfinity() const
  {
    return 1.0;
  }


  void VarGalbrunIndex_Axi
  ::ConstructAll(const string& name_file, const string& name_element,
                 Vector<string>& lines_data_file)
  {
    var_problem.SetTypeEquation("none");
    this->InitIndices(PhysicalConstant::nb_max_indices);
    
#ifdef SELDON_WITH_MPI
    int rank_proc; MPI_Comm_rank(var_problem.comm_group_mode, &rank_proc);
    int nb_proc; MPI_Comm_size(var_problem.comm_group_mode, &nb_proc);
    ReadLinesFile(name_file, lines_data_file, var_problem.comm_group_mode);
#else
    int rank_proc(0), nb_proc(1);
    ReadLinesFile(name_file, lines_data_file);
#endif
    
    ReadInputFile(lines_data_file, var_output);
    
    // construction du numero de dossier
    // si le dossier n'est pas connu
    if (rank_proc == 0)
      if (var_output.DOSSIER_output.size() == 0)
          {
            int num = -1;
            var_output.DOSSIER_output = "[STIFFOUT]/"; 
            EcritDossier(var_output.DOSSIER_output, name_file, num);
          }
    
#ifdef SELDON_WITH_MPI
    if (nb_proc > 1)
      MPI_Bcast_string(var_output.DOSSIER_output, 0, var_problem.comm_group_mode);
#endif    
    
    // we read mesh and construct reference element (for finite element method)
    bool split_mesh = true;
    if (nb_proc == 1)
      split_mesh = false;

    var_problem.ComputeMeshAndFiniteElement(name_element, split_mesh);
    
    if (nb_proc == 1)
      var_problem.mesh.Write("test.mesh");
    
    var_problem.PerformOtherInitializations();
    
  }


  /******************
   * VarGalbrun_Axi *
   ******************/


  //! allocation of arrays needed for the computation of elementary matrices
  template<class Complexe>
    void VarGalbrun_Axi<Complexe>::AllocateMassMatrices()
    {
      Glob_rtilde.Reallocate(var_problem.mesh.GetNbElt());
      Glob_radius.Reallocate(var_problem.mesh.GetNbElt());
    }


  //! computation of quantities needed for elementary matrices
  template<class Complexe>
    void VarGalbrun_Axi<Complexe>::ComputeLocalMassMatrix(int i)
    {
      int N = var_problem.Glob_PointsQuadrature(i).GetM();
      Glob_rtilde(i).Reallocate(N);
      Glob_radius(i).Reallocate(N);
      if (!var_problem.InsidePML(i))
      {
        for (int j = 0; j < N; j++)
        {
          Glob_rtilde(i)(j) = var_problem.Glob_PointsQuadrature(i)(j)(0);    
          Glob_radius(i)(j) = var_problem.Glob_PointsQuadrature(i)(j)(0);    
        }
      }
      else
      {
        Complexe rtilde;
        int i1 = i - var_problem.mesh.GetNbElt() + var_boundary.GetNbEltPML();
        for (int j = 0; j < N; j++)
        {
          rtilde = var_boundary.GetPrimitiveTauPML(i1, j, 0);
          Real_wp radius = var_problem.Glob_PointsQuadrature(i)(j)(0);

          Glob_radius(i)(j) = radius;
          Glob_rtilde(i)(j) = rtilde;
        }
      }
    }


  /*******************************
   * VolumetricSource_AxiGalbrun *
   *******************************/


  //! volumetric source
  bool VolumetricSource_AxiGalbrun
    ::IsNonNullVolumetricSource(const VectR2& s)
    {
      return true;
    }


  //! Evaluation of volumetric source f (term \int f \varphi)
  void VolumetricSource_AxiGalbrun
    ::EvaluateVolumetricSource(int i, int j, const R2& x, Vector<Complex_wp>& fvec)
    {
      R3 pt3D;
      Real_wp f;
      R3 grad_f;

      pt3D.Init(x(0), 0, x(1));
      fsrc.GetGradAmplitude(pt3D, f, grad_f);

      Complex_wp feval(f, 0);
      R2_Complex_wp grad_feval(grad_f(0), grad_f(2));
      if (var_galbrun.apply_convective_derivate_source)
      {
        Complex_wp m_iomega;
        var_problem.GetMiomega(m_iomega);
        Real_wp sigma = var_galbrun.eval_sigma(i)(j);
        Real_wp m_ = var_boundary.GetCurrentModeNumber();
        Real_wp coef = var_galbrun.coef_convective_term;
        feval *= m_iomega + sigma - Iwp*m_*coef*var_galbrun.eval_flow(i)(j)(1)/x(0);
        grad_feval *= coef;

        feval += var_galbrun.eval_flow(i)(j)(0)*grad_feval(0)
          + var_galbrun.eval_flow(i)(j)(2)*grad_feval(1);
      }    

      feval *= coef_vol*x(0);
      for (int p = 0; p < polar.GetM(); p++)
        fvec(p) = feval*polar(p);

      if (var_boundary.GetCurrentModeNumber() == 0)
      {
        fvec(0) = 0.0;
        fvec(1) = 0.0;
        fvec(3) = 0.0;
        fvec(4) = 0.0;
      }
      else
      {
        Complex_wp fr(0, 0), ftheta(0, 0), gr(0, 0), gtheta(0, 0);
        if (var_boundary.GetCurrentModeNumber() == -1)
        {
          fr = 0.5*(fvec(3) - Iwp*fvec(4));
          ftheta = 0.5*(Iwp*fvec(3) + fvec(4));            
          gr = 0.5*(fvec(0) - Iwp*fvec(1));
          gtheta = 0.5*(Iwp*fvec(0) + fvec(1));            
        }
        else if (var_boundary.GetCurrentModeNumber() == 1)
        {
          fr = 0.5*(fvec(3) + Iwp*fvec(4));
          ftheta = 0.5*(-Iwp*fvec(3) + fvec(4));
          gr = 0.5*(fvec(0) + Iwp*fvec(1));
          gtheta = 0.5*(-Iwp*fvec(0) + fvec(1));
        }

        fvec(0) = gr;
        fvec(1) = gtheta;
        fvec(2) = 0;
        fvec(3) = fr;
        fvec(4) = ftheta;
        fvec(5) = 0;
        fvec(6) = 0;
      }    
    }


  /***************************
   * GalbrunAxiEquation_Base *
   ***************************/


  // for compatbility purpose
  template<class T> template<class TypeEquation>
  void GalbrunAxiEquation_Base<T>
  ::ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                      int num_elem, const ElementReference_Dim<Dimension>& Fb)
  {
  }



  template<class T> template<class GenericPb, class T0, class Vector1>
    void GalbrunAxiEquation_Base<T>
    ::GetNeededDerivative(const GenericPb& vars,
        const GlobalGenericMatrix<T0>& nat_mat,
        Vector1& unknown_to_derive, Vector1& fct_test_to_derive)
    {
      fct_test_to_derive.Fill(true);
      unknown_to_derive.Fill(true);
    }


  template<class T> template<class TypeEquation, class T0, class Vector1>
    void GalbrunAxiEquation_Base<T>
    ::ApplyTensorStiffness(const EllipticProblem<TypeEquation>& var,
        int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
        int ref, Vector1& dU, Vector1& dV)
    {
      abort();
    }


  template<class T> template<class TypeEquation, class T0, class MatStiff>
    void GalbrunAxiEquation_Base<T>
    ::GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
        int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat,
        int ref, MatStiff& Cgrad_grad)
    {
      FillZero(Cgrad_grad);
    }


  template<class T> template<class TypeEquation, class T0, class MatStiff>
    void GalbrunAxiEquation_Base<T>
    ::GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
        int num_elem, int jloc, const GlobalGenericMatrix<T0>& nat_mat, int ref,
        MatStiff& Dgrad_phi, MatStiff& Ephi_grad)
    {
      FillZero(Dgrad_phi); FillZero(Ephi_grad);

      T0 coef = nat_mat.GetCoefStiffness();
      T radius = vars.Glob_rtilde(num_elem)(jloc);
      Real_wp c02 = square(vars.eval_c0(num_elem)(jloc));
      Real_wp rho = vars.eval_rho(num_elem)(jloc);
      R2 dp0 = vars.grad_p0(num_elem)(jloc);

      if (vars.type_model == vars.GALBRUN)
      {
        // flow term :  m_r dq/dr (for all equations)
        Ephi_grad(0, 0)(0) = radius*rho*vars.eval_flow(num_elem)(jloc)(0)*coef;
        Ephi_grad(1, 1)(0) = Ephi_grad(0, 0)(0);
        Ephi_grad(2, 2)(0) = Ephi_grad(0, 0)(0);
        Ephi_grad(3, 3)(0) = vars.coef_convective_term*Ephi_grad(0, 0)(0);
        Ephi_grad(4, 4)(0) = vars.coef_convective_term*Ephi_grad(0, 0)(0);
        Ephi_grad(5, 5)(0) = vars.coef_convective_term*Ephi_grad(0, 0)(0);
        Ephi_grad(6, 6)(0) = Ephi_grad(0, 0)(0);

        // flow term m_z dq/dz (for all equations)
        Ephi_grad(0, 0)(1) = radius*rho*vars.eval_flow(num_elem)(jloc)(2)*coef;
        Ephi_grad(1, 1)(1) = Ephi_grad(0, 0)(1);
        Ephi_grad(2, 2)(1) = Ephi_grad(0, 0)(1);
        Ephi_grad(3, 3)(1) = vars.coef_convective_term*Ephi_grad(0, 0)(1);
        Ephi_grad(4, 4)(1) = vars.coef_convective_term*Ephi_grad(0, 0)(1);
        Ephi_grad(5, 5)(1) = vars.coef_convective_term*Ephi_grad(0, 0)(1);
        Ephi_grad(6, 6)(1) = Ephi_grad(0, 0)(1);

        // term -grad( r p)/dr in equation of u
        Dgrad_phi(0, 6)(0) = radius*coef;
        Dgrad_phi(2, 6)(1) = radius*coef;

        if (vars.drop_unstable_terms != vars.DROP_NON_UNIFORM)
        {
          // m_dot_drho = (M \cdot \nabla rho_0) / rho_0
          Real_wp m_dot_drho = vars.eval_flow(num_elem)(jloc)(0)*vars.grad_rho(num_elem)(jloc)(0)
            + vars.eval_flow(num_elem)(jloc)(2)*vars.grad_rho(num_elem)(jloc)(1);

          m_dot_drho /= rho;

          // term - (M \cdot \nabla rho_0)/rho_0 \grad(p) in equation of q
          Ephi_grad(3, 6)(0) = -radius*coef*m_dot_drho;
          Ephi_grad(5, 6)(1) = -radius*coef*m_dot_drho;

          // term - dmr/dr dp/dr - dmz/dr dp/dz in equation of q_r
          Ephi_grad(3, 6)(0) -= radius*coef*vars.grad_flow(num_elem)(jloc)(0, 0);
          Ephi_grad(3, 6)(1) -= radius*coef*vars.grad_flow(num_elem)(jloc)(2, 0);

          // term duz/dz dp0/dr - duz/dr dp0/dz in equation of q_r
          Ephi_grad(3, 2)(0) = -radius*coef*vars.grad_p0(num_elem)(jloc)(1);
          Ephi_grad(3, 2)(1) = radius*coef*vars.grad_p0(num_elem)(jloc)(0);

          // term - dmr/dz dp/dr - dmz/dz dp/dz in equation of q_z
          Ephi_grad(5, 6)(0) -= radius*coef*vars.grad_flow(num_elem)(jloc)(0, 1);
          Ephi_grad(5, 6)(1) -= radius*coef*vars.grad_flow(num_elem)(jloc)(2, 1);

          // term dur/dr dp0/dz - dur/dz dp0/dr
          Ephi_grad(5, 0)(0) = radius*coef*vars.grad_p0(num_elem)(jloc)(1);
          Ephi_grad(5, 0)(1) = -radius*coef*vars.grad_p0(num_elem)(jloc)(0);                        

          // AJOUT NATHAN
          // term mt/r dp/dr in equation of q_t
          Ephi_grad(4, 6)(0) = coef*vars.eval_flow(num_elem)(jloc)(1);
        }

        // term -(rho0 c0)^2 (dur/dr + duz/dz) in equation of p
        Ephi_grad(6, 0)(0) = -coef*rho*rho*c02*radius;
        Ephi_grad(6, 2)(1) = -coef*rho*rho*c02*radius;

        return;
      }

      // term -1/r d/dr( r u) - d/dz( u_z) in equation of p
      Dgrad_phi(6, 0)(0) = radius*coef;
      Dgrad_phi(6, 2)(1) = Dgrad_phi(6, 0)(0);

      // term - \grad (rho_0 c_0^2 p) in equation of v
      Dgrad_phi(3, 6)(0) = radius*rho*c02*coef;
      Dgrad_phi(5, 6)(1) = Dgrad_phi(3, 6)(0);

      // flow term :  m_r dq/dr
      Ephi_grad(0, 0)(0) = radius*rho*vars.eval_flow(num_elem)(jloc)(0)*coef;
      Ephi_grad(1, 1)(0) = Ephi_grad(0, 0)(0);
      Ephi_grad(2, 2)(0) = Ephi_grad(0, 0)(0);
      Ephi_grad(3, 3)(0) = Ephi_grad(0, 0)(0);
      Ephi_grad(4, 4)(0) = Ephi_grad(0, 0)(0);
      Ephi_grad(5, 5)(0) = Ephi_grad(0, 0)(0);

      // flow term m_z dq/dz
      Ephi_grad(0, 0)(1) = radius*rho*vars.eval_flow(num_elem)(jloc)(2)*coef;
      Ephi_grad(1, 1)(1) = Ephi_grad(0, 0)(1);
      Ephi_grad(2, 2)(1) = Ephi_grad(0, 0)(1);
      Ephi_grad(3, 3)(1) = Ephi_grad(0, 0)(1);
      Ephi_grad(4, 4)(1) = Ephi_grad(0, 0)(1);
      Ephi_grad(5, 5)(1) = Ephi_grad(0, 0)(1);

      // term du_z/dz dp0/dr  in equation of v_r
      Ephi_grad(3, 2)(1) = radius*coef*dp0(0);

      // term -du_z/dr dp0/dz / (-i omega + sigma) in equation of v_r
      Ephi_grad(3, 2)(0) = -radius*coef*dp0(1);        

      // term du_r/dr dp0/dz / (-i omega + sigma) in equation of v_z
      Ephi_grad(5, 0)(0) = radius*coef*dp0(1);

      // term -du_r/dz dp0/dr / (-i omega + sigma) in equation of v_z
      Ephi_grad(5, 0)(1) = -radius*coef*dp0(0);        
    }


  template<class T>
    template<class TypeEquation, class T0, class Vector1, class Vector2>
    void GalbrunAxiEquation_Base<T>
    ::ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
        int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
        int ref, Vector1& dU, Vector2& V)
    {
      abort();
    }


  template<class T> template<class TypeEquation, class T0, class MatMass>
    void GalbrunAxiEquation_Base<T> 
    ::GetTensorMass(const EllipticProblem<TypeEquation>& vars,
        int num_elem, int jloc,
        const GlobalGenericMatrix<T0>& nat_mat, int ref, MatMass& mass)
    {
      mass.Fill(0);

      Complex_wp im = Iwp*Real_wp(vars.GetCurrentModeNumber());
      T radius = vars.Glob_rtilde(num_elem)(jloc);
      Real_wp c02 = square(vars.eval_c0(num_elem)(jloc));
      Real_wp rho = vars.eval_rho(num_elem)(jloc);
      Real_wp sigma = vars.eval_sigma(num_elem)(jloc);
      T0 coef = nat_mat.GetCoefStiffness(), s = coef;

      // term rho (-i omega + sigma) in u and v equations
      T m_iomega; vars.GetMiomega(m_iomega);
      mass(0, 0) = radius*(m_iomega*nat_mat.GetCoefMass() + sigma*nat_mat.GetCoefDamping())*rho;
      mass(1, 1) = mass(0, 0);
      mass(2, 2) = mass(0, 0);
      mass(3, 3) = mass(0, 0);
      mass(4, 4) = mass(0, 0);
      mass(5, 5) = mass(0, 0);

      // term -i m rho m_phi p/r in equations of u and v
      coef *= Complex_wp(0, -vars.GetCurrentModeNumber())*vars.eval_flow(num_elem)(jloc)(1)*rho;    
      mass(0, 0) += coef;
      mass(1, 1) += coef;
      mass(2, 2) += coef;
      coef *= vars.coef_convective_term;
      mass(3, 3) += coef;
      mass(4, 4) += coef;
      mass(5, 5) += coef;

      // AJOUT NATHAN
      // terms from material derivarive in polar coordinates
      mass(0,1) -= s*vars.eval_flow(num_elem)(jloc)(2);
      mass(1,0) += s*vars.eval_flow(num_elem)(jloc)(2);
      mass(3,4) -= s*vars.eval_flow(num_elem)(jloc)(2);
      mass(4,3) += s*vars.eval_flow(num_elem)(jloc)(2);

      // AJOUT NATHAN
      // term from (\nabla u)^T\nabla p_0 in polar coordinates in equation of q_thetal or v_theta
      mass(4, 2) += s*vars.grad_p0(num_elem)(jloc)(0)*radius;


      if(vars.compute_gravity) // AJOUT NATHAN
      {
        TinyMatrix<Real_wp, General, 3, 2> grad_g = vars.grad_gravity(num_elem)(jloc);
        //term u\cdot\nabla g in equation of q or v
        mass(3,0) += s*rho*grad_g(0,0)*radius;
        mass(3,2) += s*rho*grad_g(0,1)*radius;
        mass(3,1) -= s*rho*vars.eval_gravity(num_elem)(jloc)(1);

        mass(4,0) += s*rho*grad_g(1,0)*radius;
        mass(4,2) += s*rho*grad_g(1,1)*radius;
        mass(4,1) += s*rho*vars.eval_gravity(num_elem)(jloc)(0);

        mass(5,0) += s*rho*grad_g(2,0)*radius;
        mass(5,2) += s*rho*grad_g(2,1)*radius;

      }

      if (vars.type_model == vars.GALBRUN)
      {
        // evolution equation in p as well
        mass(6, 6) += mass(2, 2);

        // term p/r in equation of u_r 
        // this term appears because -dp/dr = -d/dr( r p ) + p
        mass(0, 6) += s;

        // term -rho0 q in equation of u
        mass(0, 3) += -rho*radius*s;
        mass(1, 4) += mass(0, 3);
        mass(2, 5) += mass(0, 3);

        // term + im p / r in equation of u_theta
        mass(1, 6) += im*s;

        // term - grad(sigma) p in equation of q
        mass(3, 6) = -vars.grad_sigma(num_elem)(jloc)(0)*radius*s;
        mass(5, 6) = -vars.grad_sigma(num_elem)(jloc)(1)*radius*s;

        if (vars.drop_unstable_terms != vars.DROP_NON_UNIFORM)
        {        
          // term + im/r d mtheta/dr p in equation of q_r
          mass(3, 6) += im*vars.grad_flow(num_elem)(jloc)(1, 0)*s;

          // term (u_r - i m u_theta)/r dp0/dr in equation of q_r
          mass(3, 0) += vars.grad_p0(num_elem)(jloc)(0)*s;
          mass(3, 1) -= im*vars.grad_p0(num_elem)(jloc)(0)*s;

          // term im/r (M \cdot \nabla rho0) / rho0 in equation of q_theta
          Real_wp m_dot_drho = vars.eval_flow(num_elem)(jloc)(0)*vars.grad_rho(num_elem)(jloc)(0)
            + vars.eval_flow(num_elem)(jloc)(2)*vars.grad_rho(num_elem)(jloc)(1);

          mass(4, 6) = im*m_dot_drho/rho*s;

          // term im/r (u_r dp0/dr + u_z dp0/dz) in equation of q_theta
          mass(4, 0) = im*vars.grad_p0(num_elem)(jloc)(0)*s;
          mass(4, 2) = im*vars.grad_p0(num_elem)(jloc)(1)*s;            

          // term + im/r d mtheta/dz p in equation of q_z
          mass(5, 6) += im*vars.grad_flow(num_elem)(jloc)(1, 1)*s;

          // term (u_r - i m u_theta)/r dp0/dz in equation of q_z
          mass(5, 0) += vars.grad_p0(num_elem)(jloc)(1)*s;
          mass(5, 1) -= im*vars.grad_p0(num_elem)(jloc)(1)*s;            

          // AJOUT NATHAN
          // terms from material derivative in equation of q_r
          // (-mt/(r^2)+mt/r)*dp/dth
          mass(3, 0) += s*vars.eval_flow(num_elem)(jloc)(1)*im*radius;
          mass(3, 0) -= s*vars.eval_flow(num_elem)(jloc)(1)*im;

        }


        // term -(rho0 c0)^2 (u_r - im u_theta)/r in equation of p
        mass(6, 0) += -rho*rho*c02*s;
        mass(6, 1) += rho*rho*c02*im*s;



        return;
      }

      // term -rho_0 v in equation of u
      mass(0, 3) = -rho*radius*s;
      mass(1, 4) = mass(0, 3);
      mass(2, 5) = mass(0, 3);

      // term + i m rho_0 c_0^2/r p in equation of v_theta
      mass(4, 6) = im*rho*c02*s;

      // term + i m/r u_theta in equation of p
      mass(6, 1) = im*s;

      // term rho_0 c_0^2 p (due to integration by parts of grad(rho_0 c_0^2 p)
      mass(3, 6) = rho*c02*s;

      R2 dp0 = vars.grad_p0(num_elem)(jloc);
      // term (u_r - im u_\theta)/ r dp0/dr in equation of v_r
      mass(3, 0) += s*dp0(0);
      mass(3, 1) -= im*s*dp0(0);

      // term (u_r - im u_\theta)/ r dp0/dz in equation of v_z
      mass(5, 0) += s*dp0(1);
      mass(5, 1) -= im*s*dp0(1);          

      // term im/r (u_r dp0/dr + u_z dp0/dz) in equation of v_theta
      mass(4, 0) = im*s*dp0(0);
      mass(4, 2) = im*s*dp0(1);

      // term p in equation of p
      mass(6, 6) = radius*s;
    }


  template<class T> template<class TypeEquation, class T0, class Vector1>
    void GalbrunAxiEquation_Base<T>
    ::ApplyTensorMass(const EllipticProblem<TypeEquation>& var, int i, int j,
        const GlobalGenericMatrix<T0>& nat_mat, int ref, Vector1& U, Vector1& V)
    {
      abort();
    }


  //! fills matrix |D| from matrix D
  template<class T> template<class T0>
    void GalbrunAxiEquation_Base<T>::
    GetAbsoluteD(TinyMatrix<T0, General, 7, 7>& Dtest, const R2& normale,
        const Real_wp& c0, const T0& rtilde, bool non_unif, const Real_wp& alpha1_,
        const Real_wp& alpha, const Real_wp& gamma, const Real_wp& br, const Real_wp& bz)
    {
      // axisymmetric case
      Real_wp nr = normale(0), nz = normale(1);
      Real_wp a = abs(alpha), d = 0.5*(abs(alpha+c0) - abs(alpha-c0)), s = 0.5*(abs(alpha+c0) + abs(alpha-c0));

      if (alpha > c0)
      {
        // nothing to do, Dtest is already positive
        // Dtest = Dtest;
      }
      else if (alpha < -c0)
        Dtest *= -1.0;
      else
      {
        a = abs(alpha); d = 0.5*(abs(alpha+c0) - abs(alpha-c0)); s = 0.5*(abs(alpha+c0) + abs(alpha-c0));
        Dtest(0, 0) = a + nr*nr*(s-a); Dtest(0, 1) = 0; Dtest(0, 2) = nr*nz*(s-a);
        Dtest(0, 3) = 0; Dtest(0, 4) = 0; Dtest(0, 5) = 0; Dtest(0, 6) = -nr*d/c0;

        Dtest(1, 0) = 0; Dtest(1, 1) = a; Dtest(1, 2) = 0; Dtest(1, 3) = 0;
        Dtest(1, 4) = 0; Dtest(1, 5) = 0; Dtest(1, 6) = 0;

        Dtest(2, 0) = nr*nz*(s-a); Dtest(2, 1) = 0; Dtest(2, 2) = a + nz*nz*(s-a);
        Dtest(2, 3) = 0; Dtest(2, 4) = 0; Dtest(2, 5) = 0; Dtest(2, 6) = -nz*d/c0;

        if (non_unif)
        {            
          if (alpha1_ != Real_wp(0))
          {
            Real_wp alpha1 = alpha1_*alpha;
            Real_wp a11 = (-gamma*nz + br*c0) / (alpha + c0 - alpha1);
            Real_wp a12 = (gamma*nz + br*c0) / (alpha - c0 - alpha1);
            Real_wp a21 = (gamma*nr + bz*c0) / (alpha + c0 - alpha1);
            Real_wp a22 = (-gamma*nr + bz*c0) / (alpha - c0 - alpha1);

            Real_wp sa = sign(alpha), sAc = 0.5*(abs(alpha1) - abs(alpha+c0)), dAc = 0.5*(abs(alpha-c0)-abs(alpha1));
            Dtest(3, 0) = -gamma*nr*nz*sa + nr*(a11*sAc + a12*dAc); Dtest(3, 1) = 0;
            Dtest(3, 2) = gamma*nr*nr*sa + nz*(a11*sAc + a12*dAc); Dtest(3, 3) = abs(alpha1);
            Dtest(3, 4) = 0; Dtest(3, 5) = 0; Dtest(3, 6) = (-a11*sAc + a12*dAc)/c0;

            Dtest(4, 0) = 0; Dtest(4, 1) = 0; Dtest(4, 2) = 0; Dtest(4, 3) = 0;
            Dtest(4, 4) = abs(alpha1); Dtest(4, 5) = 0; Dtest(4, 6) = 0;

            Dtest(5, 0) = -gamma*nz*nz*sa + nr*(a21*sAc + a22*dAc); Dtest(5, 1) = 0;
            Dtest(5, 2) = gamma*nr*nz*sa + nz*(a21*sAc + a22*dAc); Dtest(5, 3) = 0;
            Dtest(5, 4) = 0; Dtest(5, 5) = abs(alpha1); Dtest(5, 6) = (-a21*sAc + a22*dAc)/c0;
          }
          else
          {
            Real_wp sa = sign(alpha), sAc = 0.5*(sign(alpha+c0)+sign(alpha-c0)), dAc = 0.5*(sign(alpha+c0)-sign(alpha-c0));
            Dtest(3, 0) = nr*nz*gamma*(-sa+sAc)-nr*br*c0*dAc; Dtest(3, 1) = 0;
            Dtest(3, 2) = gamma*nr*nr*sa+nz*nz*gamma*sAc-br*nz*c0*dAc; Dtest(3, 3) = 0;
            Dtest(3, 4) = 0; Dtest(3, 5) = 0; Dtest(3, 6) = -gamma*nz/c0*dAc+br*sAc;

            Dtest(4, 0) = 0; Dtest(4, 1) = 0; Dtest(4, 2) = 0; Dtest(4, 3) = 0;
            Dtest(4, 4) = 0; Dtest(4, 5) = 0; Dtest(4, 6) = 0;

            Dtest(5, 0) = -nz*nz*gamma*sa-nr*nr*gamma*sAc-bz*nr*c0*dAc; Dtest(5, 1) = 0;
            Dtest(5, 2) = nr*nz*gamma*(sa-sAc)-bz*nz*c0*dAc;
            Dtest(5, 3) = 0; Dtest(5, 4) = 0; Dtest(5, 5) = 0; Dtest(5, 6) = gamma*nr/c0*dAc+bz*sAc;
          }
        }
        else
        {
          Dtest(3, 0) = 0; Dtest(3, 1) = 0; Dtest(3, 2) = 0; Dtest(3, 3) = a;
          Dtest(3, 4) = 0; Dtest(3, 5) = 0; Dtest(3, 6) = 0;

          Dtest(4, 0) = 0; Dtest(4, 1) = 0; Dtest(4, 2) = 0; Dtest(4, 3) = 0;
          Dtest(4, 4) = a; Dtest(4, 5) = 0; Dtest(4, 6) = 0;

          Dtest(5, 0) = 0; Dtest(5, 1) = 0; Dtest(5, 2) = 0; Dtest(5, 3) = 0;
          Dtest(5, 4) = 0; Dtest(5, 5) = a; Dtest(5, 6) = 0;
        }

        Dtest(6, 0) = -nr*d*c0; Dtest(6, 1) = 0; Dtest(6, 2) = -nz*d*c0;
        Dtest(6, 3) = 0; Dtest(6, 4) = 0; Dtest(6, 5) = 0; Dtest(6, 6) = s;

        Dtest *= rtilde;
      }
    }


  template<class T> template<class Matrix1, class GenericPb, class T0>
    void GalbrunAxiEquation_Base<T>
    ::GetAbsoluteD(Matrix1& Dtest, const R2& normale, int iquad, int k,
        const GlobalGenericMatrix<T0>& nat_mat, const GenericPb& vars)
    {
      bool non_unif = false;
      if ((vars.drop_unstable_terms == vars.DROP_CONVECTIVE)
          || (vars.drop_unstable_terms == vars.DROP_NONE))
        non_unif = true;

      // we put non-unif = false because of problems of convergence
      non_unif = false;
      Real_wp rho = vars.eval_rho(iquad)(k);
      Real_wp c0 = vars.eval_c0(iquad)(k);

      Real_wp nr = normale(0), nz = normale(1);
      Real_wp alpha = vars.eval_flow(iquad)(k)(0)*nr + vars.eval_flow(iquad)(k)(2)*nz;
      alpha *= rho;

      Real_wp m_dot_drho = vars.eval_flow(iquad)(k)(0)*vars.grad_rho(iquad)(k)(0)
        + vars.eval_flow(iquad)(k)(2)*vars.grad_rho(iquad)(k)(1);

      m_dot_drho /= rho;
      T r = vars.Glob_rtilde(iquad)(k);
      Real_wp gamma = vars.grad_p0(iquad)(k)(0)*nz - vars.grad_p0(iquad)(k)(1)*nr;
      Real_wp br = -(vars.grad_flow(iquad)(k)(0,0)*nr + vars.grad_flow(iquad)(k)(2,0)*nz + m_dot_drho*nr);
      Real_wp bz = -(vars.grad_flow(iquad)(k)(0,1)*nr + vars.grad_flow(iquad)(k)(2,1)*nz + m_dot_drho*nz);

      GetAbsoluteD(Dtest, normale, rho*c0, r, non_unif, vars.coef_convective_term,
          alpha, gamma, br, bz);
    }


  template<class T> template<class Matrix1, class GenericPb, class T0>
    void GalbrunAxiEquation_Base<T>
    ::GetNabc(Matrix1& Nabc, const R2& normale,
        int ref, int iquad, int k, const GlobalGenericMatrix<T0>& nat_mat, int ref_d, 
        const GenericPb& vars, const ElementReference<Dimension, 1>& Fb)
    {
      T0 s = nat_mat.GetCoefStiffness();
      Real_wp c0 = vars.eval_c0(iquad)(k);
      Real_wp rho = vars.eval_rho(iquad)(k);
      T radius = vars.Glob_rtilde(iquad)(k);
      switch (vars.mesh.GetBoundaryCondition(ref))
      {
        case BoundaryConditionEnum::LINE_DIRICHLET :
          {
            if (vars.type_model == vars.GALBRUN)
            {
              for (int i = 0; i < 6; i++)
                Nabc(i, 6) *= -1.0;

              Nabc *= -1.0;
              return;
            }

            Nabc.Fill(0);
            Nabc(3, 6) = -radius*rho*c0*c0*normale(0)*s;
            Nabc(5, 6) = -radius*rho*c0*c0*normale(1)*s;

            Nabc(6, 0) = radius*normale(0)*s;
            Nabc(6, 2) = radius*normale(1)*s;
          }
          break;
        case BoundaryConditionEnum::LINE_NEUMANN :
          {
            if (vars.type_model == vars.GALBRUN)
            {
              for (int i = 0; i < 6; i++)
                Nabc(i, 6) *= -1.0;

              return;
            }

            Nabc.Fill(0);
            Nabc(3, 6) = radius*rho*c0*c0*normale(0)*s;
            Nabc(5, 6) = radius*rho*c0*c0*normale(1)*s;

            Nabc(6, 0) = -radius*normale(0)*s;
            Nabc(6, 2) = -radius*normale(1)*s;
          }
          break;
        case BoundaryConditionEnum::LINE_ABSORBING :
          {
            if (vars.type_model == vars.GALBRUN)
            {
              GetAbsoluteD(Nabc, normale, iquad, k, nat_mat, vars);
              Nabc *= nat_mat.GetCoefStiffness();

              /* Nabc.Fill(0);
                 typename NatureMatrix::value_type coef = rho*c0*radius*nat_mat.GetCoefStiffness();
                 Nabc(0, 0) = coef*normale(0)*normale(0);
                 Nabc(0, 2) = coef*normale(0)*normale(1);
                 Nabc(2, 0) = coef*normale(0)*normale(1);
                 Nabc(2, 2) = coef*normale(1)*normale(1);

                 Nabc(6, 6) = coef; */

              return;
            }

            Nabc.Fill(0);
            T0 coef;
            T m_iomega; vars.GetMiomega(m_iomega);
            coef = m_iomega*rho*c0*radius*s;
            Nabc(3, 0) = coef*normale(0)*normale(0);
            Nabc(3, 2) = coef*normale(0)*normale(1);
            Nabc(5, 0) = coef*normale(1)*normale(0);
            Nabc(5, 2) = coef*normale(1)*normale(1);

            coef = c0/m_iomega*radius*s;
            Nabc(6, 6) = coef;
          }
          break;
      }
    }


  template<class T> template<class Vector1, class TypeEquation, class T0>
    void GalbrunAxiEquation_Base<T>
    ::MltNabc(const R2& normale, int ref,
        const Vector1& Vn, Vector1& Un, int num_elem1,
        int num_point, const GlobalGenericMatrix<T0>& nat_mat, int ref_d, 
        const EllipticProblem<TypeEquation>& vars,
        const ElementReference<Dimension, 1>& )
    {
      abort();
    }


  template<class T> template<class Matrix1, class TypeEquation, class T0>
    void GalbrunAxiEquation_Base<T>
    ::GetPenalDG(Matrix1& Nabc, const R2& normale, int iquad, int k,
        int nf, const GlobalGenericMatrix<T0>& nat_mat, int ref, int ref2,
        const EllipticProblem<TypeEquation>& vars,
        const ElementReference<Dimension, 1>& Fb)
    {
      Real_wp c0 = vars.eval_c0(iquad)(k);
      Real_wp rho = vars.eval_rho(iquad)(k);
      T radius = vars.Glob_radius(iquad)(k);
      T0 s = nat_mat.GetCoefStiffness();

      if (vars.type_model == vars.GALBRUN)
      {
        if (vars.upwind_fluxes)
        {
          GetAbsoluteD(Nabc, normale, iquad, k, nat_mat, vars);
          Nabc *= -nat_mat.GetCoefStiffness();
          return;
        }

        Nabc.Fill(0);
        T0 coef;
        coef = rho*c0*radius*s;

        Nabc(6, 6) = coef*vars.alpha_penalization;
        coef *= vars.delta_penalization;
        Nabc(0, 0) = coef*normale(0)*normale(0);
        Nabc(0, 2) = coef*normale(0)*normale(1);
        Nabc(2, 0) = coef*normale(0)*normale(1);
        Nabc(2, 2) = coef*normale(1)*normale(1);

        return;
      }

      Nabc.Fill(0);

      s *= vars.alpha_penalization;

      T0 coef;
      T m_iomega; vars.GetMiomega(m_iomega);
      coef = m_iomega*rho*c0*radius*s;
      Nabc(3, 0) = coef*normale(0)*normale(0);
      Nabc(3, 2) = coef*normale(0)*normale(1);
      Nabc(5, 0) = coef*normale(0)*normale(0);
      Nabc(5, 2) = coef*normale(1)*normale(1);

      coef = c0/m_iomega*radius*s;
      Nabc(6, 6) = coef;
    }


  template<class T> template<class Vector1, class Vector2, class GenericPb,
    class T0>
      void GalbrunAxiEquation_Base<T>
      ::MltPenalDG(const R2& normale, const Vector1& Vn, Vector2& Un,
          int i, int k, int nf, const GlobalGenericMatrix<T0>& nat_mat,
          int ref, int ref2, const GenericPb& vars,
          const ElementReference<Dimension, 1>& Fb)
      {
        abort();
      }


  /********************
   * VarGalbrunAxi_Eq *
   ********************/


  template<class TypeEquation>  
    VarGalbrunAxi_Eq<TypeEquation>::VarGalbrunAxi_Eq()
    : VarHarmonic<TypeEquation>(), VarGalbrun_Axi<Complexe>(this->GetLeafClass())
    {
    }


  template<class TypeEquation>  
    const R3& VarGalbrunAxi_Eq<TypeEquation>::GetPhaseOrigin() const
    {
      return VarGalbrun_Axi<Complexe>::GetPhaseOrigin();
    }


  //! returns true if numerical integration must be used
  template<class TypeEquation>  
    bool VarGalbrunAxi_Eq<TypeEquation>::UseNumericalIntegration(int i) const
    {
      return true;
    }


  //! adds Dirac for Galbrun's equation
  template<class TypeEquation>
    void VarGalbrunAxi_Eq<TypeEquation>
    ::AddDiracSource(const Real_wp& alpha, Vector<Vector<Real_wp> > & b_source,
        Vector<VirtualSourceFEM<Real_wp, Dimension2>* >& f) const
    {
      cout << "Not possible" << endl;
      abort();
    }


  //! adds Dirac for Galbrun's equation
  template<class TypeEquation>
    void VarGalbrunAxi_Eq<TypeEquation>
    ::AddDiracSource(const Complex_wp& alpha, Vector<Vector<Complex_wp> > & b_source,
        Vector<VirtualSourceFEM<Complex_wp, Dimension2>* >& f) const
    {
      VectR3 pt_source3D(1);
      pt_source3D(0) = this->origine_phase3D;

      Vector<bool> vec_unknown(7);
      vec_unknown.Fill(true);
      vec_unknown(2) = false; vec_unknown(5) = false; vec_unknown(6) = false;
      this->AddDiracSourceAxisym(alpha, vec_unknown, this->number_mode, b_source,
          pt_source3D, f);
    }


  template<class TypeEquation>
  void VarGalbrunAxi_Eq<TypeEquation>
  ::ModifyOutputUnknown(VectReal_wp& val_u, VectReal_wp& grad_u,
                        int i, const GridInterpolation<Dimension2>& var_interp,
                        int iquad, bool compute_grad) const
  {
  }
  
  
  template<class TypeEquation>
  void VarGalbrunAxi_Eq<TypeEquation>
  ::ModifyOutputUnknown(VectComplex_wp& val_u, VectComplex_wp& grad_u,
                        int i, const GridInterpolation<Dimension2>& var_interp,
                        int iquad, bool compute_grad) const
  {
      Real_wp teta = var_interp.GetTheta(i);
      Complex_wp coef = exp(-Iwp*Complex_wp(this->number_mode)*teta);
      Real_wp cos_teta = cos(teta), sin_teta = sin(teta);

      val_u(0) *= coef; val_u(1) *= coef; val_u(2) *= coef;
      val_u(3) *= coef; val_u(4) *= coef; val_u(5) *= coef; val_u(6) *= coef;
      if (compute_grad)
        for (int m = 0; m < grad_u.GetM(); m++)
          grad_u(m) *= coef;

      Complex_wp Er = val_u(0), Eteta = val_u(1), Ez = val_u(2);

      // expression of E in cartesian coordinates
      val_u(0) = Er*cos_teta - Eteta*sin_teta;
      val_u(1) = Er*sin_teta + Eteta*cos_teta;
      val_u(2) = Ez;
    }


  template<class TypeEquation>
  void VarGalbrunAxi_Eq<TypeEquation>
  ::ModifyOutputUnknown(Vector<VectReal_wp>&, Vector<VectReal_wp>&,
                        int, bool, bool) const
  {
  }
  
  
  template<class TypeEquation>
  void VarGalbrunAxi_Eq<TypeEquation>
  ::ModifyOutputUnknown(Vector<VectComplex_wp>&, Vector<VectComplex_wp>&,
                        int, bool, bool) const
  {
  }
  
  
  template<class TypeEquation>
  void VarGalbrunAxi_Eq<TypeEquation>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VarHarmonic<TypeEquation>::SetInputData(description_field, parameters);
    VarGalbrun_Axi<Complexe>::SetInputData(description_field, parameters);
    
    if (!description_field.compare("EnergyConservingAeroacousticModel"))
      {
        if (parameters(0) == "Galbrun")
          this->type_model = this->GALBRUN;
      }
  }
  
  
  template<class TypeEquation>
  void VarGalbrunAxi_Eq<TypeEquation>::InitIndices(int n)
  {
    VarGalbrun_Axi<Complexe>::InitIndices(n);
  }
  
  
  template<class TypeEquation>
  void VarGalbrunAxi_Eq<TypeEquation>
  ::GetVaryingIndices(Vector<PhysicalVaryingMedia<Dimension2, Complex_wp>* >& rho_complex,
                      Vector<PhysicalVaryingMedia<Dimension2, Real_wp>* >& rho_real, IVect& num_ref,
                      IVect& num_index, IVect& num_component, Vector<bool>& compute_grad,
                      Vector<bool>& compute_hess)
  {
    VarGalbrun_Axi<Complexe>::GetVaryingIndices(rho_real, num_ref, num_index, num_component,
                                                compute_grad, compute_hess);
  }
  
  
  template<class TypeEquation>
  void VarGalbrunAxi_Eq<TypeEquation>::ComputePhysicalCoefficients()
  {
    VarHarmonic<TypeEquation>::ComputePhysicalCoefficients();
    VarGalbrun_Axi<Complexe>::ComputePhysicalCoefficients();
  }
  
  
  template<class TypeEquation>
  void VarGalbrunAxi_Eq<TypeEquation>::AllocateMassMatrices()
  {
    VarGalbrun_Axi<Complexe>::AllocateMassMatrices();
  }
  
  
  template<class TypeEquation>
  void VarGalbrunAxi_Eq<TypeEquation>
  ::ComputeLocalMassMatrix(int i, int N, bool linear_sparse,
                           SetPoints<Dimension2>& PointsElem,
                           SetMatrices<Dimension2>& MatricesElem,
                           IVect& OrderFace, const ElementGeomReference<Dimension2>& Fb)
  {
    VarHarmonic<TypeEquation>::ComputeLocalMassMatrix(i, N, linear_sparse, PointsElem,
                                                      MatricesElem, OrderFace, Fb);
    
    VarGalbrun_Axi<Complexe>::ComputeLocalMassMatrix(i);
  }
  
  template<class TypeEquation>
  void VarGalbrunAxi_Eq<TypeEquation>::SetIndices(int i, const VectString& parameters)
  {
    VarGalbrun_Axi<Complexe>::SetIndices(i, parameters);
  }
  
  
  template<class TypeEquation>
  void VarGalbrunAxi_Eq<TypeEquation>
  ::SetPhysicalIndex(const string& name_media, int i, const VectString& parameters)
  {
    VarGalbrun_Axi<Complexe>::SetPhysicalIndex(name_media, i, parameters);
  }
  
  
  template<class TypeEquation>
  string VarGalbrunAxi_Eq<TypeEquation>::GetPhysicalIndexName(int m) const
  {
    return VarGalbrun_Axi<Complexe>::GetPhysicalIndexName(m);
  }
  
  
  template<class TypeEquation>
  int VarGalbrunAxi_Eq<TypeEquation>::GetNbPhysicalIndices() const
  {
    return VarGalbrun_Axi<Complexe>::GetNbPhysicalIndices();
  }
  
  
  template<class TypeEquation>
  bool VarGalbrunAxi_Eq<TypeEquation>::IsVaryingMedia(int i) const
  {
    return VarGalbrun_Axi<Complexe>::IsVaryingMedia(i);
  }
  
  
  template<class TypeEquation>
  bool VarGalbrunAxi_Eq<TypeEquation>::IsVaryingMedia(int m, int i) const
  {
    return VarGalbrun_Axi<Complexe>::IsVaryingMedia(i);
  }
  
  
  template<class TypeEquation>
  Real_wp VarGalbrunAxi_Eq<TypeEquation>::GetVelocityOfMedia(int ref) const
  {
    return VarGalbrun_Axi<Complexe>::GetVelocityOfMedia(ref);
  }
  
  
  template<class TypeEquation>
  Real_wp VarGalbrunAxi_Eq<TypeEquation>::GetVelocityOfInfinity() const
  {
    return VarGalbrun_Axi<Complexe>::GetVelocityOfInfinity();
  }

  
  //! retrieve and treat referenced edges
  template<class TypeEquation>
  void VarGalbrunAxi_Eq<TypeEquation>::PerformOtherInitializations()
  {        
    this->ComputeDofOnAxe(*this);
  }
  
  
  template<class TypeEquation>
  void VarGalbrunAxi_Eq<TypeEquation>::CheckInputMesh()
  {
    VarGalbrun_Axi<Complexe>::CheckSectionMeshAxi();
  }
  
  
  void EllipticProblem<HarmonicGalbrunEquationAxi>
  ::AddElementaryFluxesDG(VirtualMatrix<Real_wp>& mat_sp,
                          const GlobalGenericMatrix<Real_wp>& nat_mat,
                          int offset_row, int offset_col)
  {
    cout << "Not possible" << endl;
    abort();
  }
  
  
  void EllipticProblem<HarmonicGalbrunEquationAxi>
  ::AddElementaryFluxesDG(VirtualMatrix<Complex_wp>& mat_sp,
                          const GlobalGenericMatrix<Complex_wp>& nat_mat,
                          int offset_row, int offset_col)
  {
    Montjoie::AddElementaryFluxesDG(mat_sp, nat_mat, *this, offset_row, offset_col);
  }
  
  
  void EllipticProblem<HarmonicGalbrunEquationAxi>
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Real_wp>& A,
                            CondensationBlockSolver_Base<Real_wp>&,
                            const GlobalGenericMatrix<Real_wp>& nat_mat)
  {
    cout << "Not possible" << endl;
    abort();
  }
  
  
  void EllipticProblem<HarmonicGalbrunEquationAxi>
  ::ComputeElementaryMatrix(int iquad, IVect& num_dof, VirtualMatrix<Complex_wp>& A,
                            CondensationBlockSolver_Base<Complex_wp>&,
                            const GlobalGenericMatrix<Complex_wp>& nat_mat)
  {
    Montjoie::ComputeElementaryMatrix(iquad, num_dof, A, nat_mat, *this,
                                      this->GetReferenceElementH1(iquad));
  }
  
}

#define MONTJOIE_FILE_AXISYM_GALBRUN_CXX
#endif
