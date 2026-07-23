#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Mesh/MontjoieMesh.hxx"

using namespace Montjoie;

class InputVariables : public InputDataProblem_Base
{
public :
  // types de mailleur 3-D
  enum {MESHER_GHS3D, MESHER_GMSH, MESHER_GEOMPACK};
  // types de mailleur surfacique 3-D
  enum {SURFER_MODULEF, SURFER_MIT, SURFER_GMSH};
  // type d'adapteur surfacique 3-D
  enum {ADAPTER_MIT, ADAPTER_YAMS};
  
  // mailleur 3-D utilise pour le volume
  int type_mailleur;
  // mailleur 3-D utilise pour la peau
  int type_pavage;
  // outil utilise pour raffiner ou deraffiner les maillages de peau
  int type_adapter;
  // ordre d'approximation de la geometrie
  int order;
  
  // commandes pour ghs3d, gmsh, etc
  string gmsh_command, convert_command, apnoxx_command, tranc_command;
  string ghs3d_command, noboiteb2meshb_command, yams_command;
  
  string file_input, file_output;
  
  // pas de maillage pour toutes les primitives creees
  double step_mesh;
  // epaisseur des PMLs si on veut en rajouter
  double thickness_PML;
  // nombre de mailles dans les PMLs (dans la direction des PML)
  int nb_div;
  // surface devant etre raffinee ou deraffinee
  int type_refinement;
  // splitting of tetrahedra into hexahedra
  bool cut_tetra_asked;
  
  // nombre de surfaces importees
  int nb_meshes;
  // noms des fichiers de maillages
  VectString file_names;
  // oui si la surface i est une primitive
  VectBool primitive_base;
  // parametres des surfaces formees a partir d'une primitive
  Vector<VectString> param_primitive;
  // pas de maillage sur chaque surface
  Vector<double> local_space_step;
  
  // liste des volumes
  IVect ref_volumes;
  // pour chaque volume on liste les references composant la frontiere exterieure
  Vector<IVect> numeros_surfaces_exterieures;
  // pour chaque volume on liste les references composant chaque trou
  Vector<Vector<IVect> > numeros_surfaces_internes;

  // changement des references
  // physical_ref(i) -> numero de la surface finale
  // list_mesh_reference(i) : liste des references pour constituer la surface finale
  IVect physical_ref;
  Vector<IVect> list_mesh_reference;
  
  // constructeur par defaut
  InputVariables()
  {
    type_mailleur = MESHER_GMSH;
    type_pavage = SURFER_GMSH;
    type_adapter = ADAPTER_YAMS;
    type_refinement = 0;
    
    gmsh_command = "gmsh";
    convert_command = "convert_mesh";
    apnoxx_command = "apnoxx";
    tranc_command = "tranc";
    ghs3d_command = "ghs3d -m 200";
    noboiteb2meshb_command = "noboiteb2meshb";
    yams_command = "yams -m 200";
    
    step_mesh = 0.5;
    thickness_PML = 1.0;
    nb_div = 1; 
    nb_meshes = 0;
    cut_tetra_asked = false;
    
    local_space_step.Reallocate(100);
    local_space_step.Zero();
  }
  
  
  // lecture d'une ligne fichier donnees : InputMesh = 1 jet.mesh
  // description_field : mot-cle (dans ex : InputMesh )
  // parameters : liste des parametres (dans ex : 1 et jet.mesh)
  //              parameters(0) : 1
  //              parameters(1) : jet.mesh
  void SetInputData(const string& description_field, const VectString& parameters)
  {
    int nb_param = parameters.GetM();
    
    if (!description_field.compare("InputMesh"))
      {
        // on donne un maillage de peau en entree
	int N = nb_meshes;
	nb_meshes++;
	file_names.Resize(nb_meshes);
	primitive_base.Resize(nb_meshes);
	if (parameters(0).compare("AUTO"))
	  local_space_step(N) = to_num<double>(parameters(0));
	
	file_names(N) = parameters(1);
	// on regarde si le nom contient un . (extension)
	string::size_type loc = file_names(N).find('.',0);
	if (loc != string::npos)
	  {
	    // soit c'est un fichier de maillage
	    primitive_base(N) = false;
	  }
	else
	  {
            // soit c'est une primitive (SPHERE, PAVE, etc)
            primitive_base(N) = true;
          }
        
	// on stocke la liste des parametres associes 
	param_primitive.Resize(nb_meshes);
	param_primitive(N) = parameters;
      }
    else if (!description_field.compare("Order"))
      {
        // ordre d'approximation (pour Gmsh)
        order = to_num<int>(parameters(0));
      }
    else if (!description_field.compare("Volume"))
      {
        // definition d'un volume
	int N = ref_volumes.GetM();
	ref_volumes.Resize(N+1);
        numeros_surfaces_exterieures.Resize(N+1);
	numeros_surfaces_internes.Resize(N+1);
	
        // reference du volume (Physical Volume)
        ref_volumes(N) = to_num<int>(parameters(0));
	int k = 1;
        
        // on recupere les references pour le bord externe
	numeros_surfaces_exterieures(N).Clear();
	while ((k < nb_param)&&(parameters(k).compare("HOLE")))
	  {
	    numeros_surfaces_exterieures(N).PushBack(to_num<int>(parameters(k)));
	    k++;
	  }
	
        numeros_surfaces_internes(N).Clear();
        // boucle sur les trous
        while (k < nb_param)
          {
            // pour chaque trou, on recupere les references
            k++;
            IVect num_tmp;
            while ((k < nb_param)&&(parameters(k).compare("HOLE")))
              {
                num_tmp.PushBack(to_num<int>(parameters(k)));
                k++;
              }
            
            if (num_tmp.GetM() > 0)
              numeros_surfaces_internes(N).PushBack(num_tmp);
          }	
      }
    else if (!description_field.compare("PhysicalSurface"))
      {
        // pour regrouper des references en une Physical Surface
	int N = physical_ref.GetM();
	physical_ref.Resize(N+1);
	list_mesh_reference.Resize(N+1);
	physical_ref(N) = to_num<int>(parameters(0));
	list_mesh_reference(N).Reallocate(parameters.GetM()-1);
	for (int j = 1; j < parameters.GetM(); j++)
	  list_mesh_reference(N)(j-1) = to_num<int>(parameters(j));
      }
    else if (!description_field.compare("OutputMesh"))
      {
        // nom du fichier de sortie
	file_output = parameters(0);
      }
    else if (!description_field.compare("SpaceStep"))
      {
        // pas de maillage
	step_mesh = atof(parameters(0).data());
	if (!parameters(1).compare("FINE"))
	  type_refinement = 1;
	else if (!parameters(1).compare("COARSE"))
	  type_refinement = -1;
	else
	  type_refinement = 0;
      }
    else if (!description_field.compare("ThicknessPML"))
      {
	thickness_PML = atof(parameters(0).data());
	if (nb_param > 1)
	  nb_div = atoi(parameters(1).data());
      }
    else if (!description_field.compare("MailleurSurfacique"))
      {
	if (!parameters(0).compare("MODULEF"))
	  type_pavage = SURFER_MODULEF;
	else if (!parameters(0).compare("MIT"))
	  type_pavage = SURFER_MIT;
	else if (!parameters(0).compare("GMSH"))
	  type_pavage = SURFER_GMSH;
	// DISP(type_pavage);
      }
    else if (!description_field.compare("MailleurVolumique"))
      {
	if (!parameters(0).compare("GHS3D"))
	  type_mailleur = MESHER_GHS3D;
	else if (!parameters(0).compare("GMSH"))
	  type_mailleur = MESHER_GMSH;
	else if (!parameters(0).compare("GEOMPACK"))
	  type_mailleur = MESHER_GEOMPACK;
      }
    else if (!description_field.compare("AdapteurSurfacique"))
      {
	if (!parameters(0).compare("YAMS"))
	  type_adapter = ADAPTER_YAMS;
	else if (!parameters(0).compare("MIT"))
	  type_adapter = ADAPTER_MIT;
      }
    else if (!description_field.compare("TetrahedralMesh"))
      {
	if (!parameters(0).compare("YES"))
	  cut_tetra_asked = false;
	else
	  cut_tetra_asked = true;
      }
    else if (!description_field.compare("GmshCommand"))
      {
	gmsh_command = parameters(0);
      }
    else if (!description_field.compare("ConvertCommand"))
      {
	convert_command = parameters(0);
      }
    else if (!description_field.compare("ApnoxxCommand"))
      {
	apnoxx_command = parameters(0);
      }
    else if (!description_field.compare("TrancCommand"))
      {
	tranc_command = parameters(0);
      }
    else if (!description_field.compare("GhsCommand"))
      {
	ghs3d_command = parameters(0);
      }
    else if (!description_field.compare("GhsBoiteCommand"))
      {
	noboiteb2meshb_command = parameters(0);
      }
  }
  
};

// a partir d'un maillage triangular (tri_mesh), on recombine les triangles pour avoir 
// un maillage hybride dans quad_mesh
void RecombineSurface(const Mesh<Dimension3>& tri_mesh, Mesh<Dimension3>& quad_mesh)
{
  int ref_max = 0;
  for (int i = 0; i < tri_mesh.GetNbBoundaryRef(); i++)
    ref_max = max(ref_max, tri_mesh.BoundaryRef(i).GetReference());
  
  Mesh<Dimension3> surf_mesh; quad_mesh.Clear();
  for (int i = 1; i <= ref_max; i++)
    {
      // we write mhs file
      IVect IndexVert(tri_mesh.GetNbVertices()); IndexVert.Fill(-1);
      int nb_vertices = 0, nb_faces = 0; IVect ListeVertex(tri_mesh.GetNbVertices());
      for (int nf = 0; nf < tri_mesh.GetNbBoundaryRef(); nf++)
	if (tri_mesh.BoundaryRef(nf).GetReference() == i)
	  {
	    for (int j = 0; j < tri_mesh.BoundaryRef(nf).GetNbVertices(); j++)
	      {
		int nv = tri_mesh.BoundaryRef(nf).numVertex(j);
		if (IndexVert(nv) == -1)
		  {
		    IndexVert(nv) = nb_vertices;
		    ListeVertex(nb_vertices) = nv;
		    nb_vertices++;
		  }
	      }
	    nb_faces++;
	  }
      
      if (nb_faces > 0)
	{
	  ofstream surface_file("old_surface.mhs");
	  surface_file.precision(15);
	  surface_file << nb_vertices << '\n';
	  for (int nv = 0; nv < nb_vertices; nv++)
	    surface_file << tri_mesh.Vertex(ListeVertex(nv)) << " 1 \n";
  
	  surface_file << "\n 0 \n \n";
	  
	  surface_file << "103 "<<nb_faces<<'\n'; 
	  for (int nf = 0; nf < tri_mesh.GetNbBoundaryRef(); nf++)
	    if (tri_mesh.BoundaryRef(nf).GetReference() == i)
	      {
		for (int j = 0; j < tri_mesh.BoundaryRef(nf).GetNbVertices(); j++)
		  {
		    int nv = tri_mesh.BoundaryRef(nf).numVertex(j);
		    surface_file << IndexVert(nv) + 1 <<" ";
		  }
	      }
	  
	  for (int nf = 0; nf < tri_mesh.GetNbBoundaryRef(); nf++)
	    if (tri_mesh.BoundaryRef(nf).GetReference() == i)
	      {
		surface_file << " 5 ";
		for (int j = 0; j < tri_mesh.BoundaryRef(nf).GetNbVertices(); j++)
		  {
		    // int ne = mesh.BoundaryRef(nf).numEdge(j);
		    surface_file <<  "5 ";
		  }
		surface_file << '\n';
	      }
	  
	  surface_file.close();
	  
	  ofstream curve_file("curve_input");
	  curve_file<<"0 0"<<endl;
	  curve_file.close();
	  
	  ofstream main_file("main_input");
	  int cs_in_type(1), mesh_in_type(1), mesh_out_type(1);
	  
	  main_file<<"251 0 "<<cs_in_type<<" "<<mesh_in_type<<" 0 0 "<<mesh_out_type<<'\n';
	  main_file<<"curve_input"<<'\n';
	  main_file<<"old_surface.mhs"<<'\n';
	  main_file<<"new_surface.mhs"<<'\n';
	  main_file<<"1e-10 0 0"<<'\n';
	  int implmidn(0), nodelem(-4), mtype(3), nimpiter(0), nsmpas(0);
	  Real_wp quadmu(0.01), optmu(0.3);
	  main_file<<implmidn<<" "<<nodelem<<" "<<mtype<<" "<<quadmu<<" "<<nimpiter<<" "<<nsmpas<<" "<<optmu<<'\n';
	  main_file.close();
	  
	  // we generate a mesh file with geompack
	  string command_line = string("zgp0803_linux main_input test.out");
	  int ierr = system(command_line.data());
	  if (ierr != 0)
	    {
	      cout<<"Error during the generation of mh3 file"<<endl;
	      abort();
	    }
	  
	  surf_mesh.Read("new_surface.mhs");
	  for (int j = 0; j < surf_mesh.GetNbBoundaryRef(); j++)
	    surf_mesh.BoundaryRef(j).SetReference(i);
	  
	  quad_mesh.AppendMesh(surf_mesh);
          
          // removing temporary files
          remove("old_surface.mhs"); remove("curve_input");
          remove("main_input"); remove("new_surface.mhs");
	}
    }
  
  // we find edges from referenced faces
  quad_mesh.FindConnectivity();
  
  int nb_bounds = quad_mesh.GetNbBoundaryRef();
  quad_mesh.ResizeBoundariesRef(2*nb_bounds); int nb = nb_bounds;
  for (int i = 0; i < nb_bounds; i++)
    {
      int nb_vert = quad_mesh.BoundaryRef(i).GetNbVertices();
      int ref = quad_mesh.BoundaryRef(i).GetReference();
      if (nb_vert == 4)
	{
	  IVect num;
	  for (int j  = 0; j < nb_vert; j++)
	    {
	      int ne = quad_mesh.BoundaryRef(i).numEdge(j);
	      if (quad_mesh.GetEdge(ne).GetNbFaces() > 1)
		{
		  int num_face = quad_mesh.GetEdge(ne).numFace(0);
		  if (num_face == i)
		    num_face = quad_mesh.GetEdge(ne).numFace(1);
		  
		  num.PushBack(num_face);
		}
	    }
	  
	  if (num.GetM() > 1)
	    {
	      Sort(num);
	      for (int j = 0; j < num.GetM()-1; j++)
		if ((num(j) == num(j+1))&&(i < num(j)))
		  {
		    TinyVector<int, 4> nv;
		    for (int k = 0; k < 4; k++)
		      nv(k) = quad_mesh.BoundaryRef(i).numVertex(k);

		    int n_ext = -1;
		    for (int k = 0; k < 4; k++)
		      {
			int ns = quad_mesh.BoundaryRef(i).numVertex(k);
			bool different_vertex = true;
			for (int m = 0; m < 4; m++)
			  if (ns == quad_mesh.BoundaryRef(num(j)).numVertex(m))
			    different_vertex = false;
			
			if (different_vertex)
			  n_ext = k;
		      }
		    
		    if (n_ext%2 == 0)
		      {
			quad_mesh.BoundaryRef(i).InitTriangular(nv(0), nv(1), nv(2), ref);
			quad_mesh.BoundaryRef(nb).InitTriangular(nv(0), nv(2), nv(3), ref);
		      }
		    else
		      {
			quad_mesh.BoundaryRef(i).InitTriangular(nv(0), nv(1), nv(3), ref);
			quad_mesh.BoundaryRef(nb).InitTriangular(nv(1), nv(2), nv(3), ref);
		      }
		    nb++;
		  }
	    }
	      
	}
    }
  
  quad_mesh.ResizeBoundariesRef(nb);
  
}

// raffinement ou deraffinement d'un maillage de peau
void Simplify_or_Refine_Mesh(Mesh<Dimension3>& mesh, InputVariables& vars)
{
  string ligne;
  ofstream input_mit; string command_line;
  
  if (vars.type_adapter==vars.ADAPTER_YAMS)
    {
      // on ecrit le fichier au format mesh
      mesh.Write("totoYams.mesh");
      
      string tmps = "totoYams.bb";
      input_mit.open(tmps.data());
      input_mit<<"3 1 1 "<<mesh.GetNbVertices()<<" 2"<<endl;
      for (int i = 0; i < mesh.GetNbVertices(); i++)
	input_mit<<vars.step_mesh<<endl;
      
      input_mit.close();
      
      command_line = vars.yams_command + " -O " + to_str(vars.type_refinement) + " totoYams outputYams";
      int success_mesh = system(command_line.data());
      if (success_mesh != 0)
	{
	  cout<<"Echec lors de l'adaptation de maillage avec YAMS"<<endl;
	  exit(0);
	}
      
      mesh.Read("outputYams.meshb ");
      
      remove("totoYams.mesh");
      remove("totoYams.bb");
      remove("outputYams.meshb");
    }
}

// generation de primitives de base : pave, sphere, cylindre ...
void GeneratePrimitiveMesh(Mesh<Dimension3>& mesh, VectString& param, InputVariables& vars,
			   int num_surface, int& ref_surf)
{
  string ligne;
  ofstream input_mit; string command_line;
  input_mit.precision(15);
  Real_wp step_mesh = vars.step_mesh;
  if (vars.local_space_step(num_surface) >= 1e-7)
    step_mesh = vars.local_space_step(num_surface);
  
  Mesh<Dimension2> mesh_surf;
  IVect Index_FaceSurf_to_EdgeRef, Index_VertexSurf_to_Vertex; VectReal_wp AngleVertex;
  if (!param(1).compare("PAVE"))
    {
      Real_wp delta_boite(0);
      // DISP(param(3)); DISP(vars.xmin); DISP(vars.xmax);
      Real_wp xmin, xmax, ymin, ymax, zmin, zmax;
      // l'utilisateur rentre les bornes dans le fichier d'entree
      xmin = to_num<Real_wp>(param(2)); xmax = to_num<Real_wp>(param(3));
      ymin = to_num<Real_wp>(param(4)); ymax = to_num<Real_wp>(param(5));
      zmin = to_num<Real_wp>(param(6)); zmax = to_num<Real_wp>(param(7));
      
      int nbx = max(int(ceil((xmax-xmin + 2.0*delta_boite)/step_mesh)), 1);
      int nby = max(int(ceil((ymax-ymin + 2.0*delta_boite)/step_mesh)), 1);
      int nbz = max(int(ceil((zmax-zmin + 2.0*delta_boite)/step_mesh)), 1);
      
      // DISP(vars.type_mailleur); DISP(vars.MESHER_GEOMPACK);
      if (vars.type_mailleur == vars.MESHER_GEOMPACK)
	{
	  int nbx = toInteger(ceil((xmax - xmin + 2.0*delta_boite)/step_mesh)) + 1;
	  TinyVector<int, 6> ref_boundary; ref_boundary.Fill(ref_surf);
	  mesh.CreateRegularMesh(R3(xmin-delta_boite, ymin-delta_boite, zmin-delta_boite),
				 R3(xmax+delta_boite, ymax+delta_boite, zmax+delta_boite),
				 TinyVector<int, 3>(nbx, nbx, nbx), 1, ref_boundary, mesh.HEXAHEDRAL_MESH);
	  
          mesh.ClearElements();
	}
      else
	{
	  if (vars.type_pavage == vars.SURFER_GMSH)
	    {
	      ofstream file_geo;
              file_geo.precision(15);
	      file_geo.open("pave_tmp.geo");
	      file_geo<<"lc = "<<step_mesh<<';'<<endl;
	      file_geo<<"x1 = "<<(xmin-delta_boite)<<';'<<endl;
	      file_geo<<"x2 = "<<(xmax+delta_boite)<<';'<<endl;
	      file_geo<<"y1 = "<<(ymin-delta_boite)<<';'<<endl;
	      file_geo<<"y2 = "<<(ymax+delta_boite)<<';'<<endl;
	      file_geo<<"z1 = "<<(zmin-delta_boite)<<';'<<endl;
	      file_geo<<"z2 = "<<(zmax+delta_boite)<<';'<<endl<<endl;
	      
	      file_geo<<"Point(1) = {x1, y1, z1, lc};"<<endl;
	      file_geo<<"Point(2) = {x1, y1, z2, lc};"<<endl;
	      file_geo<<"Point(3) = {x1, y2, z1, lc};"<<endl;
	      file_geo<<"Point(4) = {x1, y2, z2, lc};"<<endl;
	      file_geo<<"Point(5) = {x2, y1, z1, lc};"<<endl;
	      file_geo<<"Point(6) = {x2, y1, z2, lc};"<<endl;
	      file_geo<<"Point(7) = {x2, y2, z1, lc};"<<endl;
	      file_geo<<"Point(8) = {x2, y2, z2, lc};"<<endl<<endl;
	      
	      file_geo<<"Line(1) = {1,2} ;"<<endl;
	      file_geo<<"Line(2) = {1,3} ;"<<endl;
	      file_geo<<"Line(3) = {2,4} ;"<<endl;
	      file_geo<<"Line(4) = {3,4} ;"<<endl;
	      file_geo<<"Line(5) = {5,6} ;"<<endl;
	      file_geo<<"Line(6) = {5,7} ;"<<endl;
	      file_geo<<"Line(7) = {6,8} ;"<<endl;
	      file_geo<<"Line(8) = {7,8} ;"<<endl;
	      file_geo<<"Line(9) = {1,5} ;"<<endl;
	      file_geo<<"Line(10) = {2,6} ;"<<endl;
	      file_geo<<"Line(11) = {3,7} ;"<<endl;
	      file_geo<<"Line(12) = {4,8} ;"<<endl<<endl;
	      
	      file_geo<<"Line Loop(1) = {1,3,-4,-2};"<<endl;
	      file_geo<<"Line Loop(2) = {4,12,-8,-11} ;"<<endl;
	      file_geo<<"Line Loop(3) = {2,11,-6,-9} ;"<<endl;
	      file_geo<<"Line Loop(4) = {9,5,-10,-1} ;"<<endl;
	      file_geo<<"Line Loop(5) = {10,7,-12,-3} ;"<<endl;
	      file_geo<<"Line Loop(6) = {6,8,-7,-5} ;"<<endl<<endl;
	      
	      file_geo<<"Plane Surface(1) = {1} ;"<<endl;
	      file_geo<<"Plane Surface(2) = {2} ;"<<endl;
	      file_geo<<"Plane Surface(3) = {3} ;"<<endl;
	      file_geo<<"Plane Surface(4) = {4} ;"<<endl;
	      file_geo<<"Plane Surface(5) = {5} ;"<<endl;
	      file_geo<<"Plane Surface(6) = {6} ;"<<endl;
	      
	      file_geo<<"Physical Surface(3) = {1,2,3,4,5,6};"<<endl;
	      
	      file_geo.close();
	      command_line = vars.gmsh_command + " -optimize -3 pave_tmp.geo -o pave_tmp.msh";
	      int success_mesh = system(command_line.data());
	      if (success_mesh != 0)
		{
		  cout<<"Echec lors de la generation du maillage de boite par gmsh"<<endl;
		  abort();
		}
	      
	      mesh.Read("pave_tmp.msh");
              remove("pave_tmp.geo");
              remove("pave_tmp.msh");
	    }
	  else if (vars.type_pavage == vars.SURFER_MODULEF)
	    {
	      Real_wp x1(0), x2(1), y1(0), y2(1); int n1(2), n2(2);
	      Mesh<Dimension2> mesh_carre;
	      // maillage des six faces avec modulef
	      int ref_m(3), ref_p(3); IVect ref(6); ref.Fill(1);
	      for (int i = 0; i < 3; i++)
		{
		  input_mit.open("tmp_nopo");
		  input_mit<<"'Carre '"<<endl;
		  if (i == 0)
		    {
		      // face  x = xmin or face x = xmax
		      x1 = ymin-delta_boite; x2 = ymax+delta_boite;
		      y1 = zmin-delta_boite; y2 = zmax+delta_boite;
		      n1 = nby+1; n2 = nbz+1;
		      ref_m = ref(0); ref_p = ref(5);
		    }
		  else if (i == 1)
		    {
		      // face y = ymin or face y = ymax
		      x1 = xmin-delta_boite; x2 = xmax+delta_boite;
		      y1 = zmin-delta_boite; y2 = zmax+delta_boite;
		      n1 = nbx+1; n2 = nbz+1;
		      ref_m = ref(1); ref_p = ref(4);
		    }
		  else if (i == 2)
		    {
		      // face z = ymin or face z = zmax
		      x1 = xmin-delta_boite; x2 = xmax+delta_boite;
		      y1 = ymin-delta_boite; y2 = ymax+delta_boite;
		      n1 = nbx+1; n2 = nby+1;
		      ref_m = ref(2); ref_p = ref(3);
		    }
		  
		  input_mit<<"!X1 = "<<x1<<endl;
		  input_mit<<"!X2 = "<<x2<<endl;
		  input_mit<<"!Y1 = "<<y1<<endl;
		  input_mit<<"!Y2 = "<<y2<<endl;
		  input_mit<<"!N1 = "<<n1<<endl;
		  input_mit<<"!N2 = "<<n2<<endl;
		  input_mit<<"'POIN  '"<<endl;
		  input_mit<<"1 4    $IMPRE NPOINT $"<<endl;
		  input_mit<<"$ NOP   NOREF(NOP)  X(NOP).  Y(NOP).  $"<<endl;
		  input_mit<<"  1  1  X1   Y1  "<<endl;
		  input_mit<<"  2  1  X2   Y1  "<<endl;
		  input_mit<<"  3  1  X2   Y2  "<<endl;
		  input_mit<<"  4  1  X1   Y2  "<<endl;
		  input_mit<<"'LIGN   '"<<endl;
		  input_mit<<"  1  4  $ IMPRE NDLM $"<<endl;
		  input_mit<<"$   NOLIG NOELIG NEXTR1 NEXTR2 NOREFL NFFRON       RAISON $"<<endl;
		  input_mit<<"       1     N1      1      2      0     0     0.100000E+01 "<<endl;
		  input_mit<<"       2     N2      2      3      0     0     0.100000E+01 "<<endl;
		  input_mit<<"       3     N1      3      4      0     0     0.100000E+01 "<<endl;
		  input_mit<<"       4     N2      4      1      0     0     0.100000E+01 "<<endl;
		  input_mit<<"'TRIA                                                                '"<<endl;
		  input_mit<<"     1     1     1     4     1          $ IMPRE NIVEAU NUDSD NBRELI NPROPA"<<endl;
		  input_mit<<"     1     2     3     4   "<<endl;
		  input_mit<<"     1     0    1                      $ NCOMP NBRINT IOPT $"<<endl;
		  input_mit<<"     4                                 $ COMPOSANTE $"<<endl;
		  input_mit<<"'SAUV             '"<<endl;
		  input_mit<<"     1  1   0   $ IMPRE NINOPO NTNOPO "<<endl;
		  input_mit<<"tmp_nopo.nopo"<<endl;
		  input_mit<<"'FIN    '"<<endl;
		  input_mit.close();
		  
		  input_mit.open("tmp_input");
		  input_mit<<"e"<<endl<<"tmp_nopo"<<endl<<"f"<<endl;
		  input_mit.close();
		  
		  // on lance apnoxx
		  command_line = vars.apnoxx_command + " < tmp_nopo > sortTmp";
		  int success_mesh = system(command_line.data());
		  if (success_mesh != 0)
		    {
		      cout<<"Echec lors du pavage avec apnoxx"<<endl;
		      exit(0);
		    }
		  
		  input_mit.open("tmp_input");
		  input_mit<<"tmp_nopo.t"<<endl<<"tmp_nopo.nopo"<<endl;
		  input_mit.close();
		  
		  // on lance tranc
		  command_line = vars.tranc_command + " < tmp_input > sortTmp";
		  success_mesh = system(command_line.data());
		  if (success_mesh != 0)
		    {
		      cout<<"Echec lors de la traduction avec tranc"<<endl;
		      exit(0);
		    }
		  
		  // on lit le maillage
		  mesh_carre.Read("tmp_nopo.t");
		  
                  // removing temporary files
                  remove("tmp_input"); remove("tmp_nopo.t"); remove("tmp_nopo");
                  remove("tmp_nopo.nopo"); remove("sortTmp");
                  
		  // on rajoute le maillage a mesh_cube
		  // sommets
		  int nb_vertices = mesh.GetNbVertices();
		  mesh.ResizeVertices(mesh.GetNbVertices() + 2*mesh_carre.GetNbVertices());
		  if (i == 0)
		    {
		      // on ajoute sommets des deux faces x = xmin and x = xmax
		      for (int j = 0; j < mesh_carre.GetNbVertices(); j++)
			{
			  mesh.Vertex(2*j + nb_vertices).Init(xmin - delta_boite, mesh_carre.Vertex(j)(0), mesh_carre.Vertex(j)(1));
			  mesh.Vertex(2*j + 1 + nb_vertices).Init(xmax + delta_boite, mesh_carre.Vertex(j)(0), mesh_carre.Vertex(j)(1));
			}
		    }
		  else if (i == 1)
		    {
		      // on ajoute sommets des deux faces y = ymin and y = ymax
		      for (int j = 0; j < mesh_carre.GetNbVertices(); j++)
			{
			  mesh.Vertex(2*j + nb_vertices).Init(mesh_carre.Vertex(j)(0), ymin - delta_boite, mesh_carre.Vertex(j)(1));
			  mesh.Vertex(2*j + 1 + nb_vertices).Init(mesh_carre.Vertex(j)(0), ymax + delta_boite, mesh_carre.Vertex(j)(1));
			}
		    }
		  else if (i == 2)
		    {
		      // on ajoute sommets des deux faces z = zmin and z = zmax
		      for (int j = 0; j < mesh_carre.GetNbVertices(); j++)
			{
			  mesh.Vertex(2*j + nb_vertices).Init(mesh_carre.Vertex(j)(0), mesh_carre.Vertex(j)(1), zmin - delta_boite);
			  mesh.Vertex(2*j + 1 + nb_vertices).Init(mesh_carre.Vertex(j)(0), mesh_carre.Vertex(j)(1), zmax + delta_boite);
			}
		    }
		  
		  // faces de references
		  int nb_faces_ref = mesh.GetNbBoundaryRef();
		  mesh.ResizeBoundariesRef(nb_faces_ref + 2*mesh_carre.GetNbElt()); IVect num(4);
		  for (int j = 0; j < mesh_carre.GetNbElt(); j++)
		    {
		      for (int k = 0; k < mesh_carre.Element(j).GetNbVertices(); k++)
			num(k) = mesh_carre.Element(j).numVertex(k);
		      
		      if (mesh_carre.Element(j).GetNbVertices() == 3)
			{
			  mesh.BoundaryRef(nb_faces_ref + 2*j).InitTriangular(nb_vertices+2*num(0), nb_vertices+2*num(1), nb_vertices+2*num(2), ref_m);
			  mesh.BoundaryRef(nb_faces_ref + 2*j + 1).
			    InitTriangular(nb_vertices+2*num(0)+1, nb_vertices+2*num(1)+1, nb_vertices+2*num(2)+1, ref_p);
			}
		      else
			{
			  mesh.BoundaryRef(nb_faces_ref + 2*j).InitQuadrangular(nb_vertices+2*num(0), nb_vertices+2*num(1),
										nb_vertices+2*num(2), nb_vertices+2*num(3), ref_m);
			  mesh.BoundaryRef(nb_faces_ref + 2*j + 1).InitQuadrangular(nb_vertices+2*num(0)+1, nb_vertices+2*num(1)+1,
										    nb_vertices+2*num(2)+1, nb_vertices+2*num(3)+1, ref_p);
			}
		    }
		  
		}
	    }
	}
      
      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	mesh.BoundaryRef(i).SetReference(ref_surf);
    }
  else if (!param(1).compare("PLATE"))
    {
      // l'utilisateur rentre les points du polygone
      Vector<R3> Points((param.GetM()-2)/3);
      for (int i = 0; i < Points.GetM(); i++)
	Points(i).Init(to_num<double>(param(2+3*i)), to_num<double>(param(3+3*i)), to_num<double>(param(4+3*i)));
      
      // origin -> first point, two tangent vectors, and normale
      R3 center, normale, vec_u, vec_v, point;
      center = Points(0); vec_u = Points(1) - Points(0);
      vec_v = Points(Points.GetM()-1) - Points(0);
      TimesProd(vec_u, vec_v, normale);
      if (Norm2(normale) < 1e-10)
	{
	  cout<<"the two vectors "<<vec_u<<" and "<<vec_v<<" are colinear "<<endl;
	  cout<<"Specify a polygon so that the first, second and last point are not on the same line"<<endl;
	  abort();
	}
      
      GetVectorPlane(normale, vec_u, vec_v);
      
      Vector<R2> Points2D(Points.GetM());
      for (int i = 0; i < Points.GetM(); i++)
	{
	  point = Points(i) - Points(0);
	  Points2D(i).Init(DotProd(point, vec_u), DotProd(point, vec_v));
	}
      
      if (vars.type_pavage == vars.SURFER_GMSH)
	{
	  ofstream file_geo;
          file_geo.precision(15);
	  file_geo.open("tmp_geo.geo");
	  file_geo<<"lc = "<<step_mesh<<';'<<endl;
	  
	  for (int j = 0; j < Points.GetM(); j++)
	    file_geo<<"Point("<<j+1<<") = {"<<Points(j)(0)<<", "<<Points(j)(1)<<", "
		    <<Points(j)(2)<<", lc};"<<endl;
	  
	  for (int j = 0; j < Points.GetM(); j++)
	    file_geo<<"Line("<<j+1<<") = {"<<j+1<<", "<<(j+1)%Points.GetM()+1<<"} ;"<<endl;
	  
	  file_geo<<"Line Loop(1) = {1";
	  for (int j = 1; j < Points.GetM(); j++)
	    file_geo<<", "<<j+1;
	  
	  file_geo<<"};"<<endl;
	  
	  file_geo<<"Plane Surface(1) = {1} ;"<<endl;
	  
	  file_geo<<"Physical Surface("<<ref_surf<<") = {1};"<<endl;
	  
	  file_geo.close();
	  command_line = vars.gmsh_command + " -optimize -3 tmp_geo.geo -o tmp_geo.msh";
	  int success_mesh = system(command_line.data());
	  if (success_mesh != 0)
	    {
	      cout<<"Echec lors de la generation du maillage de boite par gmsh"<<endl;
	      abort();
	    }
	  
	  mesh.Read("tmp_geo.msh");
          remove("tmp_geo.geo");
          remove("tmp_geo.msh");
	}
      else if (vars.type_pavage == vars.SURFER_MODULEF)
	{
	  Mesh<Dimension2> mesh_carre;
	  input_mit.open("tmp_nopo");
	  input_mit<<"'Carre '"<<endl;
	  input_mit<<"'POIN  '"<<endl;
	  input_mit<<"1 "<<Points.GetM()<<"    $IMPRE NPOINT $"<<endl;
	  input_mit<<"$ NOP   NOREF(NOP)  X(NOP).  Y(NOP).  $"<<endl;
	  for (int i = 0; i < Points.GetM(); i++)
	    input_mit<<"  "<<i+1<<" 1  "<<Points2D(i)(0)<<"  "<<Points2D(i)(1)<<endl;
	  
	  input_mit<<"'LIGN   '"<<endl;
	  input_mit<<"  1  "<<Points.GetM()<<" $ IMPRE NDLM $"<<endl;
	  input_mit<<"$   NOLIG NOELIG NEXTR1 NEXTR2 NOREFL NFFRON       RAISON $"<<endl;
	  for (int i = 0; i < Points.GetM(); i++)
	    {
	      int iafter = (i+1)%Points.GetM();
	      int n = toInteger(ceil(Points(i).Distance(Points(iafter))/step_mesh)) + 1;
	      input_mit<<"  "<<  i+1<<"  "<<n<<"  "<<i+1<<"  "<<iafter+1<<"  0     0     0.100000E+01 "<<endl;
	    }
	  
	  input_mit<<"'TRIA                                                                '"<<endl;
	  input_mit<<"     1     1     1    "<<Points.GetM()<<"     1          $ IMPRE NIVEAU NUDSD NBRELI NPROPA"<<endl;
	  for (int i = 0; i < Points.GetM(); i++)
	    input_mit<<"  "<<i+1<<"  ";
	  
	  input_mit<<endl;
	  input_mit<<"     1     0    1                      $ NCOMP NBRINT IOPT $"<<endl;
	  input_mit<<"     "<<Points.GetM()<<"                                 $ COMPOSANTE $"<<endl;
	  input_mit<<"'SAUV             '"<<endl;
	  input_mit<<"     1  1   0   $ IMPRE NINOPO NTNOPO "<<endl;
	  input_mit<<"tmp_nopo.nopo"<<endl;
	  input_mit<<"'FIN    '"<<endl;
	  input_mit.close();
	      
	  input_mit.open("tmp_input");
	  input_mit<<"e"<<endl<<"tmp_nopo"<<endl<<"f"<<endl;
	  input_mit.close();
	  
	  // on lance apnoxx
	  command_line = vars.apnoxx_command + " < tmp_input > tmp_output";
	  int success_mesh = system(command_line.data());
	  if (success_mesh != 0)
	    {
	      cout<<"Echec lors du pavage avec apnoxx"<<endl;
	      abort();
	    }
	      
	  input_mit.open("tmp_input");
	  input_mit<<"tmp_nopo.t"<<endl<<"tmp_nopo.nopo"<<endl;
	  input_mit.close();
	  
	  // on lance tranc
	  command_line = vars.tranc_command + " < tmp_input > tmp_output";
	  success_mesh = system(command_line.data());
	  if (success_mesh != 0)
	    {
	      cout<<"Echec lors de la traduction avec tranc"<<endl;
	      abort();
	    }
	  
	  // on lit le maillage
	  mesh_carre.Read("tmp_nopo.t");
          
          // removing temporary files
          remove("tmp_input"); remove("tmp_nopo.t"); remove("tmp_nopo");
          remove("tmp_nopo.nopo"); remove("tmp_output");
                            
	  // on convertit en maillage 3-D
	  mesh.ReallocateVertices(mesh_carre.GetNbVertices());
	  for (int i = 0; i < mesh_carre.GetNbVertices(); i++)
	    {
	      point = Points(0);
	      Add(mesh_carre.Vertex(i)(0), vec_u, point);
	      Add(mesh_carre.Vertex(i)(1), vec_v, point);
	      mesh.Vertex(i) = point;
	    }
	  
	  mesh.ReallocateBoundariesRef(mesh_carre.GetNbElt());
	  for (int i = 0; i < mesh_carre.GetNbElt(); i++)
	    mesh.BoundaryRef(i).
	      InitTriangular(mesh_carre.Element(i).numVertex(0), mesh_carre.Element(i).numVertex(1),
			     mesh_carre.Element(i).numVertex(2), ref_surf);
	}
      
      for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
	mesh.BoundaryRef(i).SetReference(ref_surf);
    }
  else if (!param(1).compare("SPHERE"))
    {
      Real_wp radius; R3 center;
      
      // l'utilisateur rentre les bornes dans le fichier d'entree
      center(0) = to_num<Real_wp>(param(2)); center(1) = to_num<Real_wp>(param(3));
      center(2) = to_num<Real_wp>(param(4)); radius = to_num<Real_wp>(param(5));
      
      // on genere le maillage d'un demi-disque
      int nb_div_teta = int(ceil(pi_wp*radius/step_mesh));
      Real_wp step_teta = pi_wp/nb_div_teta, teta;
      mesh_surf.ReallocateVertices(nb_div_teta+1);
      for (int i = 0; i <= nb_div_teta; i++)
	{
	  teta = -0.5*pi_wp + i*step_teta;
	  mesh_surf.Vertex(i)(0) = radius*cos(teta);
	  mesh_surf.Vertex(i)(1) = radius*sin(teta);
	}
      
      mesh_surf.ReallocateBoundariesRef(nb_div_teta);
      for (int i = 0; i < nb_div_teta; i++)
	mesh_surf.BoundaryRef(i).Init(i, i+1, ref_surf);
      
      // et on fait tourner, pour avoir la sphere complete
      IVect ref_cond(mesh_surf.GetNbReferences()+1); ref_cond.Fill();
      mesh_surf.GenerateSurfaceOfRevolution(ref_surf, mesh, Index_FaceSurf_to_EdgeRef,
					    Index_VertexSurf_to_Vertex, AngleVertex, ref_cond);
      
      // on decale la sphere
      for (int i = 0; i < mesh.GetNbVertices(); i++)
	Add(1.0, center, mesh.Vertex(i));
      
    }
  else if (!param(1).compare("ELLIPSOIDE"))
    {
      Real_wp radius_x, radius_y, radius_z; R3 center;
      
      // l'utilisateur rentre les bornes dans le fichier d'entree
      center(0) = to_num<Real_wp>(param(2)); center(1) = to_num<Real_wp>(param(3));
      center(2) = to_num<Real_wp>(param(4)); radius_x = to_num<Real_wp>(param(5));
      radius_y = to_num<Real_wp>(param(6)); radius_z = to_num<Real_wp>(param(7));
      
      // on genere le maillage d'une demi-ellipse
      int nb_div_teta = int(ceil(pi_wp*radius_x/step_mesh));
      Real_wp step_teta = pi_wp/nb_div_teta, teta;
      mesh_surf.ReallocateVertices(nb_div_teta+1);
      for (int i = 0; i <= nb_div_teta; i++)
	{
	  teta = -0.5*pi_wp + i*step_teta;
	  mesh_surf.Vertex(i)(0) = radius_x*cos(teta);
	  mesh_surf.Vertex(i)(1) = radius_z*sin(teta);
	}
      
      mesh_surf.ReallocateBoundariesRef(nb_div_teta);
      for (int i = 0; i < nb_div_teta; i++)
	mesh_surf.BoundaryRef(i).Init(i, i+1, ref_surf);
      
      // et on fait tourner, pour avoir l'ellipsoide complet
      IVect ref_cond(mesh_surf.GetNbReferences()+1); ref_cond.Fill();
      mesh_surf.GenerateSurfaceOfRevolution(ref_surf, mesh, Index_FaceSurf_to_EdgeRef,
					    Index_VertexSurf_to_Vertex, AngleVertex, ref_cond);
      
      // on decale l'ellipsoide, et on fait une homothetie
      for (int i = 0; i < mesh.GetNbVertices(); i++)
	{
	  mesh.Vertex(i)(1) *= radius_y/radius_x;
	  Add(1.0, center, mesh.Vertex(i));
	}
      // DISP(vars.zmin); DISP(vars.zmax);
      // mesh.Write("ellipsoide.mesh"); exit(0);
    }
  else if (!param(1).compare("CYLINDRE"))
    {
      Real_wp radius, height; R3 center;
      
      // l'utilisateur rentre les bornes dans le fichier d'entree
      center(0) = to_num<Real_wp>(param(2)); center(1) = to_num<Real_wp>(param(3));
      center(2) = to_num<Real_wp>(param(4)); radius = to_num<Real_wp>(param(5));
      height = to_num<Real_wp>(param(6));
      
      // on genere le maillage de deux segments
      int nb_div_radius = int(ceil(radius/step_mesh)); //DISP(nb_div_radius);
      int nb_div_axis = int(ceil(height/step_mesh)); //DISP(nb_div_axis);
      Real_wp step_radius = radius/nb_div_radius;
      Real_wp step_axis = height/nb_div_axis;
      mesh_surf.ReallocateVertices(2*nb_div_radius + nb_div_axis + 1);
      int nb = 0;
      for (int i = 0; i < nb_div_radius; i++)
	{
	  mesh_surf.Vertex(nb)(0) = step_radius*i;
	  mesh_surf.Vertex(nb)(1) = -height*0.5;
	  nb++;
	}
      
      for (int i = 0; i < nb_div_axis; i++)
	{
	  mesh_surf.Vertex(nb)(0) = radius;
	  mesh_surf.Vertex(nb)(1) = step_axis*i - height*0.5;
	  nb++;
	}
      
      for (int i = 0; i <= nb_div_radius; i++)
	{
	  mesh_surf.Vertex(nb)(0) = radius - step_radius*i;
	  mesh_surf.Vertex(nb)(1) = height*0.5;
	  nb++;
	}
      
      mesh_surf.ReallocateBoundariesRef(mesh_surf.GetNbVertices()-1);
      for (int i = 0; i < mesh_surf.GetNbBoundaryRef(); i++)
	{
	  if ((i < nb_div_radius)||(i >= (nb_div_radius + nb_div_axis)))
	    mesh_surf.BoundaryRef(i).Init(i, i+1, 2);
	  else
	    mesh_surf.BoundaryRef(i).Init(i, i+1, 1);
	}
      
      mesh_surf.SetBoundaryCondition(1, ref_surf);
      mesh_surf.SetBoundaryCondition(2, ref_surf);
      mesh_surf.FindConnectivity();
      
      // et on fait tourner, pour avoir le cylindre
      IVect ref_cond(mesh_surf.GetNbReferences()+1);
      ref_cond.Fill(ref_surf);
      mesh_surf.GenerateSurfaceOfRevolution(ref_surf, mesh, Index_FaceSurf_to_EdgeRef,
					    Index_VertexSurf_to_Vertex, AngleVertex, ref_cond);
      
      // on decale le cylindre
      for (int i = 0; i < mesh.GetNbVertices(); i++)
	Add(1.0, center, mesh.Vertex(i));
    }
  
  ref_surf++;
}

// generation d'un maillage volumique (en tetraedres plutot)
void GenerateVolumeMesh(Mesh<Dimension3>& boundary_mesh, Vector<Mesh<Dimension3> >& hole_mesh, int ref_inside,
                        const VectBool& ref_outside, InputVariables& vars, Mesh<Dimension3>& result_mesh)
{
  string ligne, ligne_sans_espace;
  ofstream input_mit; string command_line;
  input_mit.precision(15);
  //double step_mesh = vars.step_mesh;
  
  int nb_faces_ref = 0; Vector<Face<Dimension3> > faces_ref;
  if (vars.type_mailleur == vars.MESHER_GHS3D)
    {
      // on teste si les fichiers noboiteb et bb existent. Si c'est le cas, on les supprime
      ifstream file_in("tmp_ghs3d.noboiteb");
      if (file_in.is_open())
	{
	  file_in.close();
	  remove("tmp_ghs3d.noboiteb");
	}
      
      file_in.open("tmp_ghs3d.bb");
      if (file_in.is_open())
	{
	  file_in.close();
	  remove("tmp_ghs3d.bb");
	}
      
      // on concatene la frontiere exterieure avec les trous, pour avoir un seul maillage
      Mesh<Dimension3> mesh = boundary_mesh;
      for (int i = 0; i < hole_mesh.GetM(); i++)
        mesh.AppendMesh(hole_mesh(i));
      
      // on stocke les frontieres de references
      nb_faces_ref = mesh.GetNbBoundaryRef();
      faces_ref.Reallocate(nb_faces_ref);
      for (int m = 0; m < nb_faces_ref; m++)
	faces_ref(m) = mesh.BoundaryRef(m);
            
      // on ecrit le maillage
      mesh.Write("tmp_ghs3d.mesh");
      
      input_mit.open("tmp_input");
      input_mit<<"10"<<endl<<"tmp_ghs3d"<<endl<<"0"<<endl;
      input_mit.close();
      
      command_line = vars.ghs3d_command + " < tmp_input";
      int success_mesh = system(command_line.data());
      if (success_mesh != 0)
	{
	  cout<<"Echec lors du maillage ghs3d "<<endl;
	  exit(0);
	}
      
      // conversion noboiteb vers meshb
      command_line = vars.noboiteb2meshb_command + " tmp_ghs3d.noboiteb tmp_ghs3d.meshb ";
      success_mesh = system(command_line.data());
      if (success_mesh != 0)
	{
	  cout<<"Echec lors de la conversion de maillage"<<endl;
	  exit(0);
	}
      
      // avec ghs3d, on lit le maillage resultat
      result_mesh.Read("tmp_ghs3d.meshb");      
      
      // removing temporary files
      remove("tmp_ghs3d.noboiteb"); remove("tmp_ghs3d.mesh"); remove("tmp_ghs3d.meshb");
      remove("tmp_input"); remove("tmp_ghs3d.bb");
    }
  else if ((vars.type_mailleur == vars.MESHER_GMSH)||(vars.type_mailleur == vars.MESHER_GEOMPACK))
    {      
      boundary_mesh.Write("tmp_outer.mesh");
      for (int i = 0; i < hole_mesh.GetM(); i++)
        hole_mesh(i).Write(string("tmp_hole")+to_str(i)+".mesh");
      
      ofstream file_geo("tmp_geo.geo");
      file_geo << "Merge " << '"' << "tmp_outer.mesh" << '"' << ";" << endl;
      for (int i = 0; i < hole_mesh.GetM(); i++)
        file_geo << "Merge " << '"' << "tmp_hole" << i << ".mesh" << '"' << ";" << endl;
      
      file_geo << "Physical Surface(1) = {1};" << endl;
      for (int i = 0; i < hole_mesh.GetM(); i++)
        file_geo << "Physical Surface(" << i+2 << ") = {" << i+2 << "};" << endl;

      file_geo << "Surface Loop(1) = {1};" << endl;
      for (int i = 0; i < hole_mesh.GetM(); i++)
        file_geo << "Surface Loop(" << i+2 << ") = {" << i+2 << "};" << endl;

      file_geo << "Volume(1) = {1";
      for (int i = 0; i < hole_mesh.GetM(); i++)
        file_geo << ", " << i+2;
      
      file_geo << "};" << endl;
      
      file_geo << "Physical Volume(1) = {1};" << endl;
      
      file_geo.close();
      
      command_line = vars.gmsh_command + " -optimize -3 tmp_geo.geo -o tmp_geo.msh";
      int success_mesh = system(command_line.data());
      if (success_mesh != 0)
	{
	  cout<<"Echec lors de la generation du maillage de boite par gmsh"<<endl;
	  exit(0);
	}
      
      result_mesh.Read("tmp_geo.msh");
      exit(0);
      //remove("tmp_geo.geo");
      //remove("tmp_geo.msh");
      //remove("tmp_outer.mesh");
      for (int i = 0; i < hole_mesh.GetM(); i++)
        {
          boundary_mesh.AppendMesh(hole_mesh(i));
          string nom(string("tmp_hole")+to_str(i)+".mesh");
          remove(nom.data());
        }

      // sorting vertex to have direct matching
      VectR3 VertRes = result_mesh.Vertex();
      VectR3 VertBound = boundary_mesh.Vertex();
      IVect permut_res(VertRes.GetM());
      permut_res.Fill();
      IVect permut_bound(VertBound.GetM());
      permut_bound.Fill();
      
      Sort(VertRes, permut_res);
      Sort(VertBound, permut_bound);
      
      IVect original_to_new(boundary_mesh.GetNbVertices());
      original_to_new.Fill(-1);
      int k = 0;
      for (int i = 0; i < VertBound.GetM(); i++)
        {
          while ( (k < VertRes.GetM()) && (VertRes(k) < VertBound(i)))
            k++;
          
          if (k < VertRes.GetM())
            if (VertRes(k) == VertBound(i))
              original_to_new(permut_bound(i)) = permut_res(k);
        }
      
      // replacing referenced boundaries
      result_mesh.ReallocateBoundariesRef(boundary_mesh.GetNbBoundaryRef());
      for (int i = 0; i < boundary_mesh.GetNbBoundaryRef(); i++)
        {
          int nb_vert = boundary_mesh.BoundaryRef(i).GetNbVertices();
          int ref = boundary_mesh.BoundaryRef(i).GetReference(); //DISP(ref);
          IVect num(nb_vert);
          for (int j = 0; j < nb_vert; j++)
            num(j) = original_to_new(boundary_mesh.BoundaryRef(i).numVertex(j));
          
          result_mesh.BoundaryRef(i).Init(num, ref);
        }
    }
  
  if (vars.type_mailleur == vars.MESHER_GHS3D)
    {
      // on met le bon nombre de faces de references
      result_mesh.ReallocateBoundariesRef(nb_faces_ref); // DISP(nb_faces_ref);
      for (int m = 0; m < nb_faces_ref; m++)
	result_mesh.BoundaryRef(m) = faces_ref(m);
      
      // elements interieurs, tous avec une reference de 1 (ref_inside)
      for (int i = 0; i < result_mesh.GetNbElt(); i++)
	result_mesh.Element(i).SetReference(ref_inside);
    }
  else if (vars.type_mailleur == vars.MESHER_GEOMPACK)
    {
      /* 
      for (int i = 0; i < boundary_mesh.GetNbBoundaryRef(); i++)
	{
	  // we compute normale of the given boundary
	  int nf = boundary_mesh.BoundaryRef(i).GetNbVertices()-1;
	  int n0 = boundary_mesh.BoundaryRef(i).numVertex(0);
	  int n1 = boundary_mesh.BoundaryRef(i).numVertex(1);
	  int n2 = boundary_mesh.BoundaryRef(i).numVertex(nf);
	  int ref = boundary_mesh.BoundaryRef(i).GetReference();
	  R3 ptO = boundary_mesh.Vertex(n0), t1, t2, normale, vec_u, center;
	  t1 = boundary_mesh.Vertex(n1) - ptO;
	  t2 = boundary_mesh.Vertex(n2) - ptO;
	  TimesProd(t1, t2, normale);
	  
	  // we search a face in mesh that is contained in the considered face of boundary_mesh
	  int nmin = min(n0, n1); nmin = min(nmin, boundary_mesh.BoundaryRef(i).numVertex(2)); nmin = min(nmin, n2);
	  IVect num_face_near; mesh.GetNeighboringBoundariesAroundVertex(nmin, num_face_near);
	  int num_face = -1; Real_wp dist_min = 1e300;
	  for (int k = 0; k < num_face_near.GetM(); k++)
	    {
	      center.Zero();
	      for (int p = 0; p < mesh.Boundary(num_face_near(k)).GetNbVertices(); p++)
		center += mesh.Vertex(mesh.Boundary(num_face_near(k)).numVertex(p));
	      
	      Mlt(1.0/mesh.Boundary(num_face_near(k)).GetNbVertices(), center);
	      center -= ptO;
	      Real_wp dist_cur = abs(DotProd(center, normale));
	      if ( dist_cur < dist_min)
		{
		  num_face = num_face_near(k);
		  dist_min = dist_cur;
		}
	    }
	  
	  int num_elem = mesh.Boundary(num_face).numElement(0);
	  int ns = -1;
	  for (int k = 0; k < 4; k++)
	    {
	      int n_vert = mesh.Element(num_elem).numVertex(k);
	      vec_u = ptO - mesh.Vertex(n_vert);
	      if (abs(DotProd(vec_u, normale)) > 1e-6)
		ns = n_vert;
	    }
	  
	  vec_u = ptO - mesh.Vertex(ns);
	  if (DotProd(normale, vec_u) < 0)
	    {
	      if (nf == 3)
		{
		  ns = boundary_mesh.BoundaryRef(i).numVertex(2);
		  boundary_mesh.BoundaryRef(i).InitQuadrangular(n0, n2, ns, n1, ref);
		}
	      else
		{
		  boundary_mesh.BoundaryRef(i).InitTriangular(n0, n2, n1, ref);
		}
	    }

	}
      
      // boundary_mesh.Write("surface_quad.mesh");
      boundary_mesh.Write("new_surface.mhs");
      
      // generation of geompack files
      ofstream main_file("main_input");
      ofstream curve_file("curve_input");
      curve_file<<"0 0"<<endl;
      curve_file.close();
      
      main_file<<"351 0 1 1 0 0 1"<<endl;
      main_file<<"curve_input"<<endl;
      main_file<<"new_surface.mhs"<<endl;
      main_file<<"mesh_output.mh3"<<endl;
      main_file<<"1e-10 0 0"<<endl;
      
      int implmidn(0), nodelem(-8), mtype(3), nimpiter(0), nsmpas(0), regcod(0);
      Real_wp quadmu(0.1), optmu(0.3), mergmu(0.05);
      main_file << "3 " << implmidn << " " << regcod << " " << nodelem << " " << mtype
		<< " " << quadmu << " " << nimpiter << " " << nsmpas << " " << optmu << " " << mergmu << endl;
      
      main_file.close();
            
      // we generate a mesh file with geompack
      command_line = string("zgp0803_linux main_input test.out");
      int ierr = system(command_line.data());
      if (ierr != 0)
	{
	  cout<<"Error during the generation of mh3 file"<<endl;
	  abort();
	}
      
      mesh.Read(string("mesh_output.mh3"));
      
      for (int i = 0; i < mesh.GetNbElt(); i++)
	mesh.Element(i).SetReference(ref_inside);
        
      */
    }
  
}

int main(int argc, char** argv)
{

  InitMontjoie(argc, argv);
  
  if (argc != 2)
    {
      cout<<"Cette commande demande un argument"<<endl;
      cout<<"tetmesh config.ini"<<endl;
      cout<<"Elle permet de construire un maillage tetraedrique autour d'un objet a partir "<<endl;
      cout<<"du maillage surfacicque de l'objet "<<endl;
    }
  
  InputVariables vars;
  ReadInputFile(string(argv[1]), vars);
  
  Vector<Mesh<Dimension3> > mesh_surface, mesh_quad;
  Vector<Mesh<Dimension3> > mesh_volume;
  
  Mesh<Dimension3> mesh; mesh.ResizeNbReferences(100);
  //mesh.SetThicknessPML(vars.thickness_PML);
  R3::threshold = 1e-5;
  
  // premiere etape : on lit tous les maillages de surface
  // on les decoupe si demande
  mesh_surface.Reallocate(vars.nb_meshes);
  mesh_quad.Reallocate(vars.nb_meshes);
  int ref_surface = 1;
  for (int i = 0; i < vars.nb_meshes; i++)
    {
      if (vars.primitive_base(i))
	{
	  // primitives de base
	  GeneratePrimitiveMesh(mesh_surface(i), vars.param_primitive(i), vars, i, ref_surface);
	}
      else
	{
	  // lecture du maillage
	  mesh_surface(i).Read(vars.file_names(i));
          
          VectString& param = vars.param_primitive(i);	  
	  if (param.GetM() > 1)
	    {
	      if (!param(2).compare("TRANSLATE"))
		{
		  R3 vec_u;
		  vec_u(0) = to_num<Real_wp>(param(3));
		  vec_u(1) = to_num<Real_wp>(param(4));
		  vec_u(2) = to_num<Real_wp>(param(5));
		  
		  if (param.GetM() > 6)
		    {
		      IVect ref_cond(mesh_surface(i).GetNbReferences()+1);
		      ref_cond.Fill(0);
		      for (int j = 6; j < param.GetM(); j++)
			{
			  int ref = to_num<int>(param(j));
			  if (ref < ref_cond.GetM())
			    ref_cond(ref) = 1;
			}
		      
		      for (int j = 0; j < mesh_surface(i).GetNbBoundaryRef(); j++)
			{
			  int ref = mesh_surface(i).BoundaryRef(j).GetReference();
			  if (ref > 0)
			    if (ref_cond(ref) == 0)
			      mesh_surface(i).BoundaryRef(j).SetReference(0);
			}
		      
		      mesh_surface(i).RemoveReference(0);
		    }
		  
		  TranslateMesh(mesh_surface(i), vec_u);
		}
	      else if (!param(2).compare("ROTATE"))
		{
		  R3 center, axis; Real_wp teta;
		  center(0) = to_num<Real_wp>(param(3));
		  center(1) = to_num<Real_wp>(param(4));
		  center(2) = to_num<Real_wp>(param(5));
		  
		  axis(0) = to_num<Real_wp>(param(6));
		  axis(1) = to_num<Real_wp>(param(7));
		  axis(2) = to_num<Real_wp>(param(8));
		  
		  teta = pi_wp*to_num<Real_wp>(param(9))/180.0;
		  
		  if (param.GetM() > 10)
		    {
		      IVect ref_cond(mesh_surface(i).GetNbReferences()+1);
		      ref_cond.Fill(0);
		      for (int j = 10; j < param.GetM(); j++)
			{
			  int ref = to_num<int>(param(j));
			  if (ref < ref_cond.GetM())
			    ref_cond(ref) = 1;
			}
		      
		      for (int j = 0; j < mesh_surface(i).GetNbBoundaryRef(); j++)
			{
			  int ref = mesh_surface(i).BoundaryRef(j).GetReference();
			  if (ref > 0)
			    if (ref_cond(ref) == 0)
			      mesh_surface(i).BoundaryRef(j).SetReference(0);
			}
		      
		      mesh_surface(i).RemoveReference(0);
		    }
		  
		  RotateMesh(mesh_surface(i), center, axis, teta);
		}
	      else
		{
		  if (param.GetM() > 2)
		    {
		      IVect ref_cond(mesh_surface(i).GetNbReferences()+1);
		      ref_cond.Fill(0);
		      for (int j = 2; j < param.GetM(); j++)
			{
			  int ref = to_num<int>(param(j));
			  if (ref < ref_cond.GetM())
			    ref_cond(ref) = 1;
			}
		      
		      for (int j = 0; j < mesh_surface(i).GetNbBoundaryRef(); j++)
			{
			  int ref = mesh_surface(i).BoundaryRef(j).GetReference();
			  if (ref > 0)
			    if (ref_cond(ref) == 0)
			      mesh_surface(i).BoundaryRef(j).SetReference(0);
			}
		      
		      mesh_surface(i).RemoveReference(0);
		    }
		}
	    }

          // on change les references
	  VectBool ref_used(10); ref_used.Fill(false);
	  for (int j = 0; j < mesh_surface(i).GetNbBoundaryRef(); j++)
	    {
	      int ref = mesh_surface(i).BoundaryRef(j).GetReference();
	      if (ref >= ref_used.GetM())
		ref_used.Resize(ref+1);
	      
	      ref_used(ref) = true;
	    }
	  
	  IVect new_ref_number(ref_used.GetM()); new_ref_number.Fill(-1);
	  for (int j = 0; j < ref_used.GetM(); j++)
	    if (ref_used(j))
	      new_ref_number(j) = ref_surface++;
	  
	  for (int j = 0; j < mesh_surface(i).GetNbBoundaryRef(); j++)
	    {
	      int ref = mesh_surface(i).BoundaryRef(j).GetReference();
	      mesh_surface(i).BoundaryRef(j).SetReference(new_ref_number(ref));
	    }
	  
	  if (vars.type_refinement != 0)
	    Simplify_or_Refine_Mesh(mesh_surface(i), vars);
	}
      
      cout<<"Maillage de surface "<<i<<" traite avec succes "<<endl;
      
      if (vars.type_mailleur == vars.MESHER_GEOMPACK)
	{
	  if (mesh_surface(i).GetNbTrianglesRef() > 0)
	    {
	      RecombineSurface(mesh_surface(i), mesh_quad(i));
	      // mesh_quad(i) = mesh_surface(i);
	      mesh_quad(i).SplitIntoHexahedra();
	      mesh_surface(i) = mesh_quad(i);
	      mesh_surface(i).SplitIntoTetrahedra();
	    }
	  else
	    {
	      mesh_quad(i) = mesh_surface(i);
	      mesh_surface(i).SplitIntoTetrahedra();
	    }
	}
      
      mesh_surface(i).Write(string("surface")+to_str(i)+".mesh");
    }
  
  IVect all_ref(ref_surface); all_ref.Fill();
  VectBool ref_outside(ref_surface);
  
  // deuxieme etape : on genere tous les maillages volumiques
  mesh_volume.Reallocate(vars.ref_volumes.GetM());
  for (int i = 0; i < vars.ref_volumes.GetM(); i++)
    {
      Mesh<Dimension3> mesh_boundary;
      ref_outside.Fill(false);
      // on rajoute les contours exterieurs
      for (int j = 0; j < vars.numeros_surfaces_exterieures(i).GetM(); j++)
	for (int k = 0; k < vars.nb_meshes; k++)
	  {
	    int ref = vars.numeros_surfaces_exterieures(i)(j);
	    ref_outside(ref) = true;
	    bool presence_ref = false;
	    for (int p = 0; p < mesh_surface(k).GetNbBoundaryRef(); p++)
	      if (mesh_surface(k).BoundaryRef(p).GetReference() == ref)
		presence_ref = true;
	    
	    if (presence_ref)
	      {
		SurfacicMesh<Dimension3> mesh_surf;
		mesh_surface(k).GetBoundaryMesh(ref, mesh_surf, all_ref);
                mesh_boundary.AppendMesh(mesh_surf, true);
                /* if (vars.type_mailleur == vars.MESHER_GEOMPACK)
		  {
		    mesh_quad(k).GetBoundaryMesh(ref, mesh_surf, all_ref);
		    mesh_boundary.AppendMesh(mesh_surf, true);
                    }*/
	      }
	  }
      
      // et les trous eventuels
      Vector<Mesh<Dimension3> > mesh_hole(vars.numeros_surfaces_internes(i).GetM());
      for (int nhole = 0; nhole < vars.numeros_surfaces_internes(i).GetM(); nhole++)
        for (int j = 0; j < vars.numeros_surfaces_internes(i)(nhole).GetM(); j++)
          for (int k = 0; k < vars.nb_meshes; k++)
            {
              int ref = vars.numeros_surfaces_internes(i)(nhole)(j);
              bool presence_ref = false;
              for (int p = 0; p < mesh_surface(k).GetNbBoundaryRef(); p++)
                if (mesh_surface(k).BoundaryRef(p).GetReference() == ref)
                  presence_ref = true;
              
              if (presence_ref)
                {
                  SurfacicMesh<Dimension3> mesh_surf;
                  mesh_surface(k).GetBoundaryMesh(ref, mesh_surf, all_ref);
                  mesh_hole(nhole).AppendMesh(mesh_surf, true);
                  /* if (vars.type_mailleur == vars.MESHER_GEOMPACK)
                    {
                      mesh_quad(k).GetBoundaryMesh(ref, mesh_surf, all_ref);
                      quad_boundary.AppendMesh(mesh_surf, true);
                      }*/
                }
            }
      
      // on genere le maillage 3-D
      GenerateVolumeMesh(mesh_boundary, mesh_hole, vars.ref_volumes(i), ref_outside, vars, mesh_volume(i));
    }
  
  // troisieme etape : on concatene tout le monde en un seul maillage
  for (int i = 0; i < vars.ref_volumes.GetM(); i++)
    mesh.AppendMesh(mesh_volume(i), true);
  
  // references a changer
  // si on a une physical surface a 0, on enleve des surfaces
  int nb_surf = 0;
  IVect new_ref(ref_surface+1); new_ref.Fill(-1);
  for (int i = 0; i < vars.physical_ref.GetM(); i++)
    for (int j = 0; j < vars.list_mesh_reference(i).GetM(); j++)
      new_ref(vars.list_mesh_reference(i)(j)) = vars.physical_ref(i);
  
  for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
    {
      int ref = new_ref(mesh.BoundaryRef(i).GetReference());
      if (ref > 0)
	{
	  if (nb_surf != i)
	    mesh.BoundaryRef(nb_surf) = mesh.BoundaryRef(i);
	  
	  mesh.BoundaryRef(nb_surf).SetReference(ref);
	  nb_surf++;
	}
    }
  
  if (nb_surf != mesh.GetNbBoundaryRef())
    mesh.ResizeBoundariesRef(nb_surf);
  
  // quatrieme etape : on decoupe les tetras en hexas et on rajoute la couche de PML
  if (vars.cut_tetra_asked)
    {
      mesh.SplitIntoHexahedra();
      //mesh.AddPMLElements(vars.nb_div);
    }  
  
  // derniere etape : on ecrit le maillage au format demande
  mesh.Write(vars.file_output);
    
  
  return FinalizeMontjoie();
}
