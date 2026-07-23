#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#include "Montjoie.hxx"
#include "Configuration.cxx"
#include "Helio_Data.hxx"

using namespace Montjoie;
typedef Complex_wp Complexe;


void Compute_Index(const R2& coord,const int i, Vector<Complexe>& res,Helio_Data& H,const Configuration& cfg) {

  // Fills the Vector res with all the values of the indices at point coord
  // Depending on the equation
  // All given solar data (from ModelS, or whatever) are stored in the Helio_Data structure.

  string equation_type;
  cfg.Get("TypeEquation",equation_type);

  if (equation_type == "HELMHOLTZ") {

    // Indices for scalar Helmholtz equation :
    // res(0) : 1/(\rho c^2)         (:=rho   in Montjoie)
    // res(1) : 2\gamma / (\rho c^2) (:=sigma in Montjoie)
    // res(2) : 1/\rho               (:=mu    in Montjoie)
    // res(3) : flow(radial component, cylindric coordinates)
    // res(4) : flow(orthoradial component)
    // res(5) : flow(vertical component)

    Real_wp r,theta;
    CartesianToPolar(coord(0),coord(1),r,theta);

    Real_wp rho   = H.Eval("rho"  ,r,i);
    Real_wp c     = H.Eval("c"    ,r,i);
    Real_wp gamma = H.Eval("gamma",r,i);

    Real_wp Mr     = H.Eval("Mr"    ,r,i);
    Real_wp Mtheta = H.Eval("Mtheta",r,i);
    Real_wp Mz     = H.Eval("Mz"    ,r,i);

    res.Reallocate(6);

    res(0) = Complexe(1.e0/(rho*c*c)                        ,0.e0);
    res(1) = Complexe(2.e0*gamma*2.e0*M_PI*H.R_sun/(rho*c*c),0.e0);
    res(2) = Complexe(1.e0/rho                              ,0.e0);
    res(3) = Complexe(Mr                                    ,0.e0);
    res(4) = Complexe(Mtheta                                ,0.e0);
    res(5) = Complexe(Mz                                    ,0.e0);

  }

  else if (  equation_type == "HELIO_HELMHOLTZ"
          || equation_type == "HELIO_HELMHOLTZ_V2") {

    // Indices for scalar Helmholtz equation modified for helioseismology :
    // res(0) : 1.                   (:=rho   in Montjoie)
    // res(1) : 2*\gamma             (:=sigma in Montjoie)
    // res(2) : 1/\rho               (:=mu    in Montjoie)
    // res(3) : \rho c               (:=alpha in Montjoie)
    // res(4) : c                    (:=beta  in Montjoie)
    // res(5) : flow(radial component, cylindric coordinates)
    // res(6) : flow(orthoradial component)
    // res(7) : flow(vertical component)

    Real_wp r,theta;
    CartesianToPolar(coord(0),coord(1),r,theta);

    Real_wp rho   = H.Eval("rho"  ,r,i);
    Real_wp c     = H.Eval("c"    ,r,i);
    Real_wp gamma = H.Eval("gamma",r,i);

    Real_wp Mr     = H.Eval("Mr"    ,r,i);
    Real_wp Mtheta = H.Eval("Mtheta",r,i);
    Real_wp Mz     = H.Eval("Mz"    ,r,i);

    res.Reallocate(8);

    res(0) = Complexe(1.e0/c                        ,0.e0);
    res(1) = Complexe(2.e0*gamma*2.e0*M_PI*H.R_sun/c,0.e0);
    res(2) = Complexe(1.e0/rho                    ,0.e0);
    res(3) = Complexe(rho*c                       ,0.e0);
    res(4) = Complexe(1.e0                           ,0.e0);
    res(5) = Complexe(Mr                          ,0.e0);
    res(6) = Complexe(Mtheta                      ,0.e0);
    res(7) = Complexe(Mz                          ,0.e0);

  }

}

void Compute_values_at_infinity (Vector<Complexe>& res, Helio_Data& H, const Configuration& cfg) {

  string equation_type;
  cfg.Get("TypeEquation",equation_type);

  if (equation_type == "HELMHOLTZ") {

    Real_wp rho   = H.infty_vals["rho"  ];
    Real_wp c     = H.infty_vals["c"    ];
    Real_wp gamma = H.infty_vals["gamma"];
    res.Reallocate(6);

    res(0) = Complexe(1.e0/(rho*c*c)                        ,0.e0);
    res(1) = Complexe(2.e0*gamma*H.R_sun*2.e0*M_PI/(rho*c*c),0.e0);
    res(2) = Complexe(1.e0/rho                              ,0.e0);
    res(3) = Complexe(0.e0,0.e0);
    res(4) = Complexe(0.e0,0.e0);
    res(5) = Complexe(0.e0,0.e0);

  }

  else if ( equation_type == "HELIO_HELMHOLTZ"
         || equation_type == "HELIO_HELMHOLTZ_V2") {

    Real_wp rho   = H.infty_vals["rho"  ];
    Real_wp c     = H.infty_vals["c"    ];
    Real_wp gamma = H.infty_vals["gamma"];

    res.Reallocate(8);

    res(0) = Complexe(1.e0/c                        ,0.e0);
    res(1) = Complexe(2.e0*gamma*H.R_sun*2.e0*M_PI/c,0.e0);
    res(2) = Complexe(1.e0/rho                    ,0.e0);
    res(3) = Complexe(rho*c                       ,0.e0);
    res(4) = Complexe(1.e0                           ,0.e0);
    res(5) = Complexe(0.e0,0.e0);
    res(6) = Complexe(0.e0,0.e0);
    res(7) = Complexe(0.e0,0.e0);

  }

}

void Init_coefficient(const string key, Helio_Data& H, const Configuration& cfg, const Real_wp def_val, bool& is_radial) {

  is_radial = true;

  string options;
  bool key_found = cfg.Get(key,options,true);
  if (!key_found) {
    cout << "\t\033[1m\033[34mWarning\033[0m : key \"" << key << "\" was not found in configuration file. Will use default constant " << key << "=" << def_val << endl;
    H.Set_constant(key,def_val);
  }
  else {
    istringstream iss(options);
    string type,filename;
    iss >> type;
    if (type == "CONSTANT") {
      Real_wp truc = def_val;
      if (iss) { iss >> truc;}
      H.Set_constant(key,truc); 
    }
    else if (type == "BACKGROUND") {
      int column = 2;
      if (iss) {iss >> column;}
      H.Read_modelS(cfg,key,column); 
    }
    else if (type == "RADIAL") {
      iss >> filename;
      int column = 2;
      if (iss) {iss >> column;}
      H.Read_radial_data(filename,key,column);
    }
    else if (type == "NODAL") {
      iss >> filename;
      int column = 1;
      if (iss) {iss >> column;}
      H.Read_nodal_data(cfg,filename,key,column);
      is_radial = false;
    }
    else {
      H.Set_constant(key,def_val); 
    }
    // Check if a perturbation field is added
    if (iss) {
      iss >> type;
      if (type=="ADD") {
        if (iss) {
          iss >> type >> filename;
          if (type == "CONSTANT") {
            H.Set_constant(key,atof(filename.c_str()),true); 
          }
          else if (type == "RADIAL") {
            int column = 2;
            if (iss) iss >> column;
            H.Read_radial_data(filename,key,column,true);
          }
          else if (type == "NODAL") {
            int column = 1;
            if (iss) iss >> column;
            H.Read_nodal_data(cfg,filename,key,column,true);
            is_radial = false;
          }
          else {
            H.Set_constant(key,def_val,true); 
          }
    }}}

  }

}


int main (int argc, char** argv) {

  InitMontjoie(argc, argv);

  if (argc < 3)
    {
      cout << "\t\033[1m\033[34mUsage\033[0m : ./write_index_helio config_file output_data_file" << endl;
      abort();
    }

  Configuration cfg;
  cfg.Load(string(argv[1]));

  Helio_Data H;

  int print_level = 0;
  cfg.Get("PrintLevel",print_level,true);
  bool VERBOSE = (print_level >=5);

  // ============================================================================
  // Read all data
  
  if (VERBOSE) cout << "Reading data..." << endl;

  Vector<bool>   IIR(6); // (Is Input Radial ?)
  Init_coefficient("rho"   ,H,cfg,1.e0,IIR(0));  
  Init_coefficient("c"     ,H,cfg,1.e0,IIR(1));
  Init_coefficient("gamma" ,H,cfg,1.e0,IIR(2));
  Init_coefficient("Mr"    ,H,cfg,0.e0,IIR(3));
  Init_coefficient("Mtheta",H,cfg,0.e0,IIR(4));
  Init_coefficient("Mz"    ,H,cfg,0.e0,IIR(5));

  // ============================================================================
  // Equation properties - Initialize output 

  if (VERBOSE) cout << "Initialize equation..." << endl;

  string          equation_type;
  int             Nb_indices(0);
  Vector<string>  output_files;
  Vector<string>  output_keys;
  Vector<bool>    IOR; // (Is Output Radial ?)

  

  cfg.Get("TypeEquation",equation_type);

  if (equation_type == "HELMHOLTZ") {
    Nb_indices = 6;
    output_files.Resize(Nb_indices);
    output_keys .Resize(Nb_indices);
    IOR         .Resize(Nb_indices);
    output_keys(0) = "rho";    IOR(0) = IIR(0) && IIR(1)          ; // res(0) : 1/(\rho c^2)         (:=rho   in Montjoie)
    output_keys(1) = "sigma";  IOR(1) = IIR(0) && IIR(1) && IIR(2); // res(1) : 2\gamma / (\rho c^2) (:=sigma in Montjoie)
    output_keys(2) = "mu";     IOR(2) = IIR(0)                    ; // res(2) : 1/\rho               (:=mu    in Montjoie)
    output_keys(3) = "Mr";     IOR(3) = IIR(3); // flow
    output_keys(4) = "Mtheta"; IOR(4) = IIR(4);
    output_keys(5) = "Mz";     IOR(5) = IIR(5);
  }

  else if ( equation_type=="HELIO_HELMHOLTZ"
         || equation_type=="HELIO_HELMHOLTZ_V2") {
    Nb_indices = 8; 
    output_files.Resize(Nb_indices);
    output_keys .Resize(Nb_indices);
    IOR         .Resize(Nb_indices);
    output_keys(0) = "rho";    IOR(0) = true            ; // res(0) : 1.e0    (:=rho   in Montjoie)
    output_keys(1) = "sigma";  IOR(1) = IIR(2)          ; // res(1) : 2\gamma (:=sigma in Montjoie)
    output_keys(2) = "mu";     IOR(2) = IIR(0)          ; // res(2) : 1/\rho  (:=mu    in Montjoie)
    output_keys(3) = "alpha";  IOR(3) = IIR(0) && IIR(1); // res(3) : rho*c   (:=alpha in Montjoie)
    output_keys(4) = "beta";   IOR(4) = IIR(1)          ; // res(4) : c       (:=beta  in Montjoie)
    output_keys(5) = "Mr";     IOR(5) = IIR(3); // flow
    output_keys(6) = "Mtheta"; IOR(6) = IIR(4);
    output_keys(7) = "Mz";     IOR(7) = IIR(5);
  }

  else {
    cerr << "Equation type not implemented."   << endl
         << "Implemented equation types are :" << endl
         << "HELMHOLTZ, HELIO_HELMHOLTZ"       << endl
         << "HELIO_HELMHOLTZ_V2"               << endl;
    abort();
  }

  bool all_radial = true;
  for (int index = 0; index<Nb_indices; index++) {
    cfg.Get("Input_Montjoie_"+output_keys(index),output_files(index));
    all_radial = all_radial && IOR(index);
  }

  // ============================================================================ 
  // Generate output file containing values at infinity

  string output_data_file = argv[2];
  ofstream output_file(output_data_file.c_str());
  if (!output_file.is_open()) {
    cerr << "Could not open file " << output_data_file << endl;
    abort();
  }

  Vector<Complexe> vals;
  Compute_values_at_infinity(vals,H,cfg); 
  for (int index=0;index<Nb_indices;index++)
    output_file << real(vals(index)) << " ";
  output_file << endl;
  output_file.close();

  // ============================================================================ 
  // Compute indices if radial data
  if (all_radial) {

    // Build file names
    for (int index=0; index<Nb_indices; index++)
      output_files(index) += ".don";
  
    // Compute coefs
    int NR = H.Radius.GetM();
    Vector<Vector<Complexe> > Coefs(NR);
    for (int i=0; i<NR; i++) {

      Coefs(i).Reallocate(Nb_indices);
      R2 coord(H.Radius(i),0.e0);     

      Compute_Index(coord,i,Coefs(i),H,cfg);
    }

    // Write .don files
    for (int index=0; index<Nb_indices; index++) {
      ofstream of;
      of.open(output_files(index).c_str());
      of.precision(16);
      if (!of.is_open()) {
        cerr << "Could not open file " << output_files(index) << " to write coefficient." << endl;
        abort();
      }
      for (int i=0; i<NR; i++)
        of << H.Radius(i) << " \t" << real(Coefs(i)(index)) << endl;
   
      of.close();
    }

  }

  // ============================================================================ 
  // Compute indices on nodal points
  
  else {

    // Build file names
    for (int index=0; index<Nb_indices; index++)
      output_files(index) += ".elb";

    if (VERBOSE) cout << "Computing coefficients on points..." << endl;
  
    string nodal_points_file;
    cfg.Get("Nodal_points_file",nodal_points_file);
  
    ifstream NPF(nodal_points_file.c_str());
    if (!NPF.is_open()) {
      cerr << "Could not open nodal points file " << nodal_points_file << endl;
      abort();
    }
  
    int Nb_elemts,Nb_nodal_points;
    NPF >> Nb_elemts >> Nb_nodal_points;
  
    Vector<Vector<Complexe> > Coefs(Nb_nodal_points);
    // Value of index #i at point #j : Coefs(j)(i)
  
    for (int ID_pt = 0; ID_pt < Nb_nodal_points; ID_pt++) {
      R2 coord;
      NPF >> coord(0) >> coord(1);
      Compute_Index(coord,ID_pt,Coefs(ID_pt),H,cfg);
    }
  
  
  // ============================================================================ 
  // Build vector that contains the index values element by element
  
    if (VERBOSE) cout << "Build elements vectors..." << endl;
  
    Vector<Vector<Vector<Complexe> > > Coefs_elt(Nb_indices);
    // Value of index #i on the j-th nodal point of element k : Coefs_elt(i)(k)(j)
  
    for (int index = 0; index < Nb_indices; index++)
      Coefs_elt(index).Reallocate(Nb_elemts);
  
  
    for (int ID_elt=0; ID_elt<Nb_elemts; ID_elt++) {
  
      int Nb_nodes_in_elt;
      NPF >> Nb_nodes_in_elt;
      for (int index=0; index<Nb_indices; index++)
        Coefs_elt(index)(ID_elt).Reallocate(Nb_nodes_in_elt);
  
      int ID_pt;
      for (int ID_pt_loc=0; ID_pt_loc<Nb_nodes_in_elt; ID_pt_loc++) {
        NPF >> ID_pt;
        for (int index=0; index<Nb_indices; index++)
          Coefs_elt(index)(ID_elt)(ID_pt_loc) = Coefs(ID_pt)(index);
      }
  
    }
  
    NPF.close();
  
    
  // ============================================================================ 
  // Write binary output files                                     
    
    if (VERBOSE) cout << "Writing index files..." << endl;
  
    for (int index=0; index<Nb_indices; index++)
      WriteMeshData(Coefs_elt(index),output_files(index));

  }

  // ============================================================================ 
  // Ciao                                                          

  return FinalizeMontjoie();

}

