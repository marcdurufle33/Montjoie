#define MONTJOIE_WITH_ONE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

Real_wp threshold;

template<class T>
bool EqualVector(const Vector<T>& x, const Vector<T>& y)
{
  if (x.GetM() != y.GetM())
    {
      return false;
    }

  for (int i = 0; i < x.GetM(); i++)
    if (isnan(abs(x(i)-y(i))) || (abs(x(i)-y(i)) > threshold))
      {
	DISP(i); DISP(x(i)); DISP(y(i));
	return false;
      }

  return true;
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  threshold = 1e4*epsilon_machine;
  
  Mesh<Dimension1> mesh;

  int ref = 2; int N = 150;
  Real_wp xmin(-3), xmax(2);
  mesh.CreateRegularMesh(xmin, xmax, N, ref);

  if (mesh.Vertex().GetM() != N)
    {
      cout << "CreateRegularMesh incorrect" << endl;
      abort();
    }

  mesh.Write("points.dat");

  VectReal_wp x_step;
  Linspace(xmin, xmax, N, x_step);
  if (!EqualVector(x_step, mesh.Vertex()))
    {
      cout << "Vertex or CreateRegularMesh incorrect" << endl;
      abort();
    }

  if (mesh.GetNbVertices() != N)
    {
      cout << "GetNbVertices incorrect"  << endl;
      abort();
    }
  
  for (int i = 0; i < mesh.GetNbVertices(); i++)
    if (isnan(mesh.Vertex(i)) || (abs(mesh.Vertex(i) - x_step(i)) > threshold))
      {
	cout << "Vertex incorrect"  << endl;
	abort();
      }

  if (mesh.GetNbElt() != N-1)
    {
      cout << "GetNbElt incorrect"  << endl;
      abort();      
    }

  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      if ((mesh.Element(i).numVertex(0) != i)
	  || (mesh.Element(i).numVertex(1) != i+1)
	  || (mesh.Element(i).GetReference() != ref))
	{
	  cout << "Element or CreateRegularMesh incorrect" << endl;
	  abort();
	}
    }
  
  if (isnan(mesh.GetXmin()) || (abs(mesh.GetXmin()-xmin) > threshold))
    {
      cout << "GetXmin incorrect"  << endl;
      abort();
    }

  if (isnan(mesh.GetXmax()) || (abs(mesh.GetXmax()-xmax) > threshold))
    {
      cout << "GetXmax incorrect"  << endl;
      abort();
    }

  int r = 8;
  mesh.SetOrder(r);

  if (mesh.GetOrder() != r)
    {
      cout << "GetOrder/SetOrder incorrect" << endl;
      abort();
    }

  mesh.NumberMesh();
  if (mesh.GetNbDof() != (N-1)*r + 1)
    {
      cout << "GetNbDof/NumberMesh incorrect" << endl;
      abort();
    }
  
  int offset = 0;
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      if (mesh.GetOrderElement(i) != r)
	{
	  cout << "GetOrderElement incorrect" << endl;
	  abort();
	}

      for (int j = 0; j <= r; j++)
	if (mesh.GetNumberDof(i, j) != offset+j)
	  {
	    cout << "GetNumberDof incorrect" << endl;
	    abort();
	  }
      
      offset += r;
    }

  mesh.SetPeriodicExtremity();
  r += 2;
  mesh.SetOrder(r);
  mesh.NumberMesh(false);

  if (mesh.GetNbDof() != (N-1)*r)
    {
      cout << "GetNbDof/NumberMesh incorrect" << endl;
      abort();
    }
  
  offset = 0;
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      if (mesh.GetOrderElement(i) != r)
	{
	  cout << "GetOrderElement incorrect" << endl;
	  abort();
	}

      int jmax = r-1;
      if (i != mesh.GetNbElt()-1)
	jmax = r;
      
      for (int j = 0; j <= jmax; j++)
	if (mesh.GetNumberDof(i, j) != offset+j)
	  {
	    cout << "GetNumberDof incorrect" << endl;
	    abort();
	  }
      
      offset += r;
    }

  if (mesh.GetNumberDof(0, 0) != mesh.GetNumberDof(N-2, r))
    {
      cout << "SetPeriodicExtremity incorrect" << endl;
      abort();
    }

  mesh.NumberMesh(true);

  if (mesh.GetNbDof() != (N-1)*(r+1))
    {
      cout << "GetNbDof/NumberMesh incorrect" << endl;
      abort();
    }
  
  offset = 0;
  for (int i = 0; i < mesh.GetNbElt(); i++)
    {
      if (mesh.GetOrderElement(i) != r)
	{
	  cout << "GetOrderElement incorrect" << endl;
	  abort();
	}

      for (int j = 0; j <= r; j++)
	if (mesh.GetNumberDof(i, j) != offset+j)
	  {
	    cout << "GetNumberDof incorrect" << endl;
	    abort();
	  }
      
      offset += r+1;
    }

  if (abs(mesh.GetMeshSize() - (x_step(1)-x_step(0))) > threshold)
    {
      cout << "GetMeshSize() incorrect" << endl;
      abort();
    }
  
  VectReal_wp pos;
  pos.ReadText("points.dat");

  if (!EqualVector(pos, mesh.Vertex()))
    {
      cout << "Write incorrect" << endl;
      abort();
    }

  pos.Reallocate(N+10);
  pos.FillRand(); pos *= 1e-9;
  Sort(pos);
  pos.WriteText("points.dat");

  mesh.Read("points.dat");
  if (!EqualVector(pos, mesh.Vertex()))
    {
      cout << "Read incorrect" << endl;
      abort();
    }

  cout << "Taille prise par le maillage = " << GetHumanReadableMemory(mesh.GetMemorySize()) << endl;

  Mesh<Dimension1> mesh2; IVect ref_domain(pos.GetM());
  ref_domain.FillRand();
  mesh2.CreateNonRegularMesh(pos, ref_domain);

  if ((mesh2.GetNbElt() != pos.GetM()-1) || (mesh2.GetNbVertices() != pos.GetM())
      || (!EqualVector(mesh2.Vertex(), pos)))
    {
      cout << "CreateNonRegularMesh incorrect " << endl;
      abort();
    }

  for (int i = 0; i < mesh2.GetNbElt(); i++)
    if (mesh2.Element(i).GetReference() != ref_domain(i))
      {
	cout << "CreateNonRegularMesh incorrect " << endl;
	abort();
      }

  int nb_layer = 12; IVect num_cell_layer(nb_layer);
  ref_domain.Reallocate(nb_layer);  
  pos.Reallocate(nb_layer+1);
  pos.FillRand(); pos *= 1e-9; Sort(pos);
  VectReal_wp pos_ref;
  for (int i = 0; i < nb_layer; i++)
    {
      ref_domain(i) = i;
      num_cell_layer(i) = rand()%20+1;
      Linspace(pos(i), pos(i+1), num_cell_layer(i)+1, x_step);
      if (i < nb_layer-1)
	x_step.Resize(num_cell_layer(i));
      
      pos_ref.PushBack(x_step);
    }
  
  mesh.CreateLayeredMesh(pos, num_cell_layer, ref_domain);
  if ((mesh.GetNbVertices() != pos_ref.GetM()) || (mesh.GetNbElt() != pos_ref.GetM()-1)
      || (!EqualVector(mesh.Vertex(), pos_ref)))
    {
      cout << "CreateLayeredMesh incorrect" << endl;
      abort();
    }

  offset = 0;
  for (int i = 0; i < nb_layer; i++)
    {
      for (int j = 0; j < num_cell_layer(i); j++)
	if (mesh.Element(offset+j).GetReference() != ref_domain(i))
	  {
	    cout << "CreateLayeredMesh incorrect" << endl;
	    abort();	    
	  }
      
      offset += num_cell_layer(i);
    }

  Mesh<Dimension1> fine_mesh;
  mesh2.CreateNonRegularMesh(pos, ref_domain);
  mesh2.SubdivideMesh(fine_mesh, num_cell_layer);

  if (!EqualVector(fine_mesh.Vertex(), mesh.Vertex()))
    {
      cout << "SubdivideMesh incorrect" << endl;
      abort();
    }

  Mesh<Dimension1> sub_mesh;
  int n0 = rand()%(mesh.GetNbVertices()-2);
  int n1 = rand()%(mesh.GetNbVertices()-1-n0) + n0;
  mesh.CreateSubmesh(sub_mesh, n0, n1);

  if (sub_mesh.GetNbVertices() != n1 - n0 + 1)
    {
      DISP(n0); DISP(n1); DISP(n1-n0); DISP(mesh.GetNbVertices());
      cout << "CreateSubmesh incorrect" << endl;
      abort();
    }

  for (int j = n0; j <= n1; j++)
    if (sub_mesh.Vertex(j-n0) != mesh.Vertex(j))
      {
	cout << "CreateSubmesh incorrect" << endl;
	abort();
      }

  sub_mesh.Clear();
  if (sub_mesh.GetNbVertices() != 0)
    {
      cout << "Clear incorrect" << endl;
      abort();
    }

  Real_wp delta(0.2);
  mesh.SetThicknessPML(delta);
  if (mesh.GetThicknessPML() != delta)
    {
      cout << "GetThicknessPML/SetThicknessPML incorrect" << endl;
      abort();
    }

  nb_layer = 5;
  mesh.SetAdditionPML(mesh.PML_BOTH_SIDES, nb_layer);
  if (mesh.GetTypeAdditionPML(0) != mesh.PML_BOTH_SIDES)
    {
      cout << "GetTypeAdditionPML/SetAdditionPML incorrect" << endl;
      abort();
    }

  if (mesh.GetNbLayersPML() != nb_layer)
    {
      cout << "GetNbLayersPML incorrect" << endl;
      abort();
    }

  VectReal_wp new_vertex, add_vertex;
  Linspace(mesh.GetXmin()-delta, mesh.GetXmin(), nb_layer+1, new_vertex);
  new_vertex.Resize(nb_layer);
  new_vertex.PushBack(mesh.Vertex());
  new_vertex.Resize(new_vertex.GetM()-1);
  
  Linspace(mesh.GetXmax(), mesh.GetXmax()+delta, nb_layer+1, add_vertex);
  new_vertex.PushBack(add_vertex);
  
  mesh.AddPMLElements(nb_layer);
  if (!EqualVector(mesh.Vertex(), new_vertex))
    {
      cout << "AddPMLElements incorrect" << endl;
      abort();
    }

  VectString parameters;
  parameters.PushBack("LAYERED");
  parameters.PushBack("0.1");
  parameters.PushBack("0.33");
  parameters.PushBack("0.45");
  parameters.PushBack("AUTO");
  parameters.PushBack("0.05");
  parameters.PushBack("REFINED");
  parameters.PushBack("2");

  mesh.level_refinement_vertex.PushBack(4);
  mesh.ratio_refinement_vertex.PushBack(4.0);
  mesh.position_refinement_vertex.PushBack(0.33);

  mesh.nb_layers_to_add.PushBack(7);
  mesh.position_layers_to_add.PushBack(0.6);
  mesh.reference_layers_to_add.PushBack(13);
  
  mesh.ConstructMesh(parameters, xmin, xmax);
  ref_domain.Reallocate(mesh.GetNbElt());
  for (int i = 0; i < mesh.GetNbElt(); i++)
    ref_domain(i) = mesh.Element(i).GetReference();

  IVect list_ref;
  list_ref.ReadText("src/Program/Unit/Mesh/reference_mesh1D.dat");
  pos.ReadText("src/Program/Unit/Mesh/points_mesh1D.dat");
  if (!EqualVector(pos, mesh.Vertex()))
    {
      cout << "ConstructMesh incorrect" << endl;
      abort();
    }

  if (!EqualVector(list_ref, ref_domain))
    {
      cout << "ConstructMesh incorrect" << endl;
      abort();
    }
  
  std::remove("points.dat");
  
  cout << "All tests passed successfully" << endl;
  
  return FinalizeMontjoie();
}
