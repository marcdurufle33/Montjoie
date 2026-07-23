

// version simplifiee de ConstructFiniteElement (dans VarProblem.cxx)
void ModeEz_Solver::ConstructFiniteElement(const string& name_elt)
{
  // previous finite elements are cleared if present
  this->ClearFiniteElement();

  // orders present in the mesh
  TinyVector<IVect, 4> order;
  this->mesh_num.GetOrder(order);

  Vector<bool> change_elt(this->mesh.GetNbElt());
  change_elt.Fill(true);

  this->AddFiniteElement(name_elt, order, change_elt, this->mesh, this->mesh_num,
                         0, false);

  this->mesh_num.GetOrderQuadrature(order);
  RemoveDuplicate(order(0));

  int type_integration_edge, type_integration_tri, type_integration_quad;
  Vector<MeshNumbering<Dimension2>* > other_mesh_num;
  this->GetTypeIntegrationBoundary(this->mesh, this->mesh_num, other_mesh_num,
                                   type_integration_edge,
           type_integration_tri, type_integration_quad);

  this->mesh_num.number_map.ConstructQuadrature2D(order, type_integration_edge);
}

// traitement de Dirichlet, quasi-periodique et ddl conductors
// À COMMENTER POUR BIEN COMPRENDRE
void ModeEz_Solver::RenumberDofs()
{
  // Sert à dire si on applique dirichlet sur l'interface correspondant à la ref (surfacique)
  Vector<bool> is_dirichlet(this->mesh.GetNbReferences()+1);
  is_dirichlet.Fill(false);
  for (int i = 0; i < this->ref_bord.GetM(); i++)
    is_dirichlet(ref_bord(i)) = true;

  // on met le ddl 0 en Dirichlet (dans le cas quasi-periodique avec m <> 0)
  Vector<int> DofDirichlet(this->mesh_num.GetNbDof());
  DofDirichlet.Fill(-1);
  // Sur les dof : -1 = pas de condition
  //                0 = Dirichlet homogène
  //                1 = Quasi-périodique

  int nb_dof_periodic = this->mesh_num.GetNbPeriodicDof();
  int nb_dof_dirichlet = 0;

  // en cas de cable tiers, on fixe le point central à 0 (condition nécessaire la quasi-périodicité)
  if ((nb_dof_periodic > 0) && (this->num_mode_periodic != 0))
    {
      DofDirichlet(0) = 0;
      nb_dof_dirichlet++;
    }

  for (int i = 0; i < nb_dof_periodic; i++)
    {
      int n = this->mesh_num.GetPeriodicDof(i);
      if (DofDirichlet(n) == -1)
        DofDirichlet(n) = 1;
    }
  cout << "Number of periodic dofs = " << nb_dof_periodic << endl;

  // là on a un array sur les éléments, contenant un array de booléens sur leurs dofs
  // on initialise tous ces points à faux (ils ne sont a priori pas des points périodiques)
  this->is_dof_with_phase.Reallocate(this->mesh.GetNbElt());
  for (int i = 0; i < this->mesh.GetNbElt(); i++)
    {
      this->is_dof_with_phase(i).Reallocate(this->mesh_num.Element(i).GetNbDof());
      this->is_dof_with_phase(i).Fill(false);
    }

  // on renumerote
  if ((nb_dof_dirichlet > 0) || (nb_dof_periodic > 0))
    {
      Vector<int> IndexDof(this->mesh_num.GetNbDof());
      IndexDof.Fill(-1);
      int nodl = 0;
      for (int i = 0; i < this->mesh_num.GetNbDof(); i++)
        {
          if (DofDirichlet(i) == -1)
            IndexDof(i) = nodl++;
        }

      for (int i = 0; i < nb_dof_periodic; i++)
        {
          int ndof = this->mesh_num.GetPeriodicDof(i);
          int n0 = this->mesh_num.GetOriginalPeriodicDof(i);
          if (DofDirichlet(ndof) == 1)
            {
              //DISP(i); DISP(n0+1); DISP(ndof+1);
              IndexDof(ndof) = -IndexDof(n0) - 2;
            }
        }

      //DISP(IndexDof);

      for (int i = 0; i < this->mesh.GetNbElt(); i++)
        {
          //DISP(i); DISP(this->mesh_num.Element(i).GetNodle());
          for (int j = 0; j < this->mesh_num.Element(i).GetNbDof(); j++)
            {
              int num_dof = this->mesh_num.Element(i).GetNumberDof(j);
              if (IndexDof(num_dof) <= -2)
                {
                  this->mesh_num.Element(i).SetNumberDof(j, -IndexDof(num_dof)-2);
                  this->is_dof_with_phase(i)(j) = true;
                }
              else
                this->mesh_num.Element(i).SetNumberDof(j, IndexDof(num_dof));
            }

          // DISP(this->is_dof_with_phase(i));
          //DISP(this->mesh_num.Element(i).GetNodle());
        }

      this->mesh_num.SetNbDof(nodl);
    }

  // boucle sur les aretes de bord pour lister les ddls Dirichlets
  nb_dof_dirichlet = 0;
  DofDirichlet.Fill(-1);
  for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
    {
      int ref = this->mesh.BoundaryRef(i).GetReference();
      if (is_dirichlet(ref))
        {
          int num_elem = this->mesh.BoundaryRef(i).numElement(0);
          int num_loc = this->mesh.Element(num_elem).GetPositionBoundary(i);
          const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(num_elem);
          const IVect& Nodle = this->mesh_num.Element(num_elem).GetNodle();
          int nb_dof = Fb.GetNbDofBoundary(num_loc);
          for (int j = 0; j < nb_dof; j++)
            {
              int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
              int num_dof = Nodle(num_dof_loc);
              if ((num_dof >= 0) && (DofDirichlet(num_dof) == -1))
                {
                  DofDirichlet(num_dof) = 0;
                  nb_dof_dirichlet++;
                }
            }
        }
    }

  cout << "Number of border dofs = " << nb_dof_dirichlet << endl;

  int nb_dof_border = 0;
  IndexNonDirichlet.Reallocate(this->mesh_num.GetNbDof());
  IndexNonDirichlet.Fill(-1);
  IndexLocalBorder.Reallocate(this->mesh_num.GetNbDof());
  IndexLocalBorder.Fill(-1);
  rhoInfini = Complex_wp(0.0, 0.0);
  DofKeptNonDir.Reallocate(this->mesh_num.GetNbDof()-nb_dof_dirichlet);
  nb_dof_dirichlet = 0;
  for (int i = 0; i < this->mesh_num.GetNbDof(); i++)
    if (DofDirichlet(i) == -1)
      {
        IndexNonDirichlet(i) = nb_dof_dirichlet;
        DofKeptNonDir(nb_dof_dirichlet) = i;
        nb_dof_dirichlet++;
      }

  // on recupere les ddls des conducteurs
  DofDirichlet.Fill(-1);
  ddl_conductor.Reallocate(ref_conductor.GetM());
  Vector<int> index_conductor(this->mesh.GetNbReferences()+1);
  index_conductor.Fill(-1);
  for (int i = 0; i < ref_conductor.GetM(); i++)
    {
      index_conductor(ref_conductor(i)) = i;
      ddl_conductor(i).Clear();
    }

  for (int i = 0; i < this->mesh.GetNbBoundaryRef(); i++)
    {
      int ref = this->mesh.BoundaryRef(i).GetReference();
      if (index_conductor(ref) >= 0)
        {
          int ic = index_conductor(ref);
          int num_elem = this->mesh.BoundaryRef(i).numElement(0);
          int num_loc = this->mesh.Element(num_elem).GetPositionBoundary(i);
          const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(num_elem);
          const IVect& Nodle = this->mesh_num.Element(num_elem).GetNodle();
          int nb_dof = Fb.GetNbDofBoundary(num_loc);
          for (int j = 0; j < nb_dof; j++)
            {
              int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
              int num_dof = Nodle(num_dof_loc);
              if ((num_dof >= 0) && (DofDirichlet(num_dof) == -1))
                {
                  DofDirichlet(num_dof) = 0;
                  ddl_conductor(ic).PushBack(num_dof);
                }
            }
        }

      if (ref == ref_bord(0))
        {
          int num_elem = this->mesh.BoundaryRef(i).numElement(0);
          int num_loc = this->mesh.Element(num_elem).GetPositionBoundary(i);
          const ElementReference<Dimension2, 1>& Fb = this->GetReferenceElementH1(num_elem);
          const IVect& Nodle = this->mesh_num.Element(num_elem).GetNodle();
          if (realpart(rhoInfini) == 0.0)
            {
              ref_water = this->mesh.Element(num_elem).GetReference();
              rhoInfini = this->ref_epsilon(ref_water)*this->ref_mu(ref_water)
                          + Iwp * this->ref_sigma(ref_water) * this->ref_mu(ref_water) / this->omega;
            }
          int nb_dof = Fb.GetNbDofBoundary(num_loc);
          for (int j = 0; j < nb_dof; j++)
            {
              int num_dof_loc = Fb.GetLocalNumber(num_loc, j);
              int num_dof = Nodle(num_dof_loc);
              // cout << "Arête " << i << ", ddl local " << j << endl;
              // DISP(num_dof_loc);
              // DISP(num_dof);
              if (IndexLocalBorder(num_dof) == -1)
                {
                  IndexLocalBorder(num_dof) = nb_dof_border;
                  NumGlobBorder.PushBack(num_dof);
                  nb_dof_border++;
                }
            }
        }
    }

  for (int i = 0; i < ref_conductor.GetM(); i++)
    {
      DISP(ref_conductor(i));
      DISP(ref_vol_conductor(i));
      //DISP(ddl_conductor(i));
    }
}
