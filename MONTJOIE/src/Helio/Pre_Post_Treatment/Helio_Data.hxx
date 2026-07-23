#ifndef HELIO_DATA
#define HELIO_DATA

using namespace Montjoie;

class Helio_Data {

  public:

    Helio_Data() {
      R_sun   = 6.968e8;
    };
   ~Helio_Data() {};

    // Return functions
    // ==========================================================================

    Real_wp Eval(const string key,const Real_wp r, const int i) {
      double res = 0.e0;
      switch(storage[key]) {
        case 0  : res = cvals[key]  ; break;
        case 1  : res = Eval (key,r); break;
        case 2  : res = Eval (key,i); break;
      }
      if (storage_add[key]) {
        switch(storage_add_type[key]) {
          case 0  : res += cvals_add[key]  ; break;
          case 1  : res += Eval(key,r,true); break;
          case 2  : res += Eval(key,i,true); break;
        }
      }
      return res*adim_coefs[key];
    };

    Real_wp Eval(const string key, const int i, const bool add=false) {
      if (!add) {
        if (storage[key]!=2) {
          cerr << "Data given as radial must be evaluated with a given radius!" << endl;
          abort();
        }
        return nodal_data[key](i);
      }
      else {
        if (storage_add_type[key]!=2) {
          cerr << "Data given as radial must be evaluated with a given radius!" << endl;
          abort();
        }
        return nodal_data_add[key](i);
      }
    };



    Real_wp Eval(const string key, const Real_wp r, const bool add=false) {
      if (!add) {
        switch(storage[key]) {
          case 0  : return cvals[key];            break;
          case 1  : return spls[key].Evaluate(r); break;
          default : cerr << "Nodal points data cannot take a radius as parameter to be evaluated" << endl;
                    abort();
        }
      }
      else {
        switch(storage_add_type[key]) {
          case 0  : return cvals_add[key];            break;
          case 1  : return spls_add[key].Evaluate(r); break;
          default : cerr << "Nodal points data cannot take a radius as parameter to be evaluated" << endl;
                    abort();
        }
      }
     
  
    };

    // Constant media Setters
    // ========================================================================== 
    void Set_constant(const string key, const Real_wp val, const bool add = false) {
      if (!add) { 
        cvals      [key] = val;
        storage    [key] = 0;
        storage_add[key] = false;
        infty_vals [key] = val;
        adim_coefs [key] = 1.e0;
      }
      else { 
        cvals_add       [key]  = val;
        storage_add     [key]  = true;
        storage_add_type[key]  = 0;
        infty_vals      [key] += val;
      }
    }

    // ==========================================================================
    // ModelS file reading methods
    //
    // ModelS file should contain SI values of r/R, c, rho, G and T, in this order
    // Comment lines should start by the '%' or '#' characters

    void Read_modelS(const Configuration& cfg,const string key,const int column) {
      if (key!="rho" && key!="c") {
        cerr << "Value corresponding to key \"" << key << "\" is not present in modelS!" << endl
             << "Accepted keys : rho c " << endl;
        abort(); 
      }
      string filename;
      cfg.Get("Background_File",filename);
      Read_radial_data(filename,key,column);
    };

    // ==========================================================================
    // Read data from radial input file

    void Read_radial_data (const string filename,const string key,const int column=2, const bool add=false) {

      ifstream f;
      // Get number of modelS points
      f.open(filename.c_str());
      if (!f.is_open()) {
        cerr << "Could not open file " << filename << "to read radial data." << endl;
        abort();
      }
      string junk;
      int    Npts = 0;
      while (getline(f,junk)) {
        if (junk.at(0)!='%' && junk.at(0)!='#')
          Npts++;
      }
      f.close();

      // Read data
      f.open(filename.c_str());
      Vector<Real_wp> data_val(Npts);
      Radius.Reallocate(Npts);
      int line = 0;

      while (getline(f,junk)) { if (junk.at(0)!='%' && junk.at(0)!='#') {
        istringstream s(junk);
        s >> Radius(line);
        for (int i=2;i<=column;i++)
          s >> data_val(line);
        line++;
      }} 
      f.close();


      // Init spline
      if (!add) {
        spls[key].Init(Radius,data_val);
        storage    [key] = 1;
        storage_add[key] = false;
        infty_vals [key] = data_val(Npts-1);
        adim_coefs [key] = 1.e0;
      }
      else {
        spls_add[key].Init(Radius,data_val);
        storage_add     [key] = true;
        storage_add_type[key] = 1;
      }
    } 

    // ==========================================================================
    // Read data on nodal points from input file

    void Read_nodal_data(const Configuration& cfg, const string filename,const string key, const int column=1, const bool add=false) {

      // Check if the number of points is the same in the nodal points file and the given data
      string nodal_pts_file;
      cfg.Get("Nodal_points_file",nodal_pts_file);
      ifstream nf;
      nf.open(nodal_pts_file.c_str());
      if (!nf.is_open()) {
        cerr << "Could not open nodal points file " << nodal_pts_file << ". Will end now" << endl;
        abort();
      }
      int Npoints_nodal;
      // The first line of the file contains the number of elements and then the number of nodal points
      nf >> Npoints_nodal >> Npoints_nodal;
      int ind_rmax = -1;
      Real_wp rmax = -1337.e42;
      R2      c;
      for (int i=0;i<Npoints_nodal;i++) {
        nf >> c(0) >> c(1);
        Real_wp r = Norm2(c);
        if (r>rmax) {
           ind_rmax = i;
           rmax     = r;
        }
      }

      nf.close();     
 
      ifstream f;
      f.open(filename.c_str());
      if (!f.is_open()) {
        cerr << "Could not data file on nodal points " << filename << ". Will end now" << endl;
        abort();
      }
      string junk;
      int    Npts = 0;
      while (getline(f,junk)) {
        if (junk.at(0)!='%' && junk.at(0)!='#')
          Npts++;
      }
      f.close();
      
      if (Npts!=Npoints_nodal) {
        cerr << "Number of points in " << nodal_pts_file << "(" << Npoints_nodal << ") "
             << "does not match number of given values in " << filename << "(" << Npts << ")." << endl;
        abort();
      }

      // Read data
      f.open(filename.c_str());
      Vector<Real_wp> data_val(Npts);

      int line = 0;

      while (getline(f,junk)) { if (junk.at(0)!='%' && junk.at(0)!='#') {
        istringstream s(junk);
        for (int i=1;i<=column;i++)
          s >> data_val(line);
        line++;
      }} 
      f.close();

      // Init vectors
      if (!add) {
        nodal_data [key] = data_val;
        storage    [key] = 2;
        storage_add[key] = false;
        infty_vals [key] = data_val(ind_rmax);
        adim_coefs [key] = 1.e0;
      }
      else {
        nodal_data_add  [key] = data_val;
        storage_add     [key] = true;
        storage_add_type[key] = 2;
      }

    };

    void Set_type(const string key,const int i) {
      storage[key] = i;
    };

    void Set_type_add(const string key,const int i) {
      storage_add[key] = true;
      storage_add_type[key] = i;
    };

    // Multiply by a constant value
    void Multiply(const string key, const Real_wp coef) {
      adim_coefs[key] = coef;
    }

    // Values at infinity
    map<string,Real_wp> infty_vals;

    Real_wp R_sun;
    Vector<Real_wp> Radius;

  private:

    // Select the type of storage (0=constant, 1=spline, 2=vector)    
    map<string,int> storage;

    // constants
    map<string,Real_wp> cvals;

    // Splines for radial data (mainly modelS), but other properties could be added that way
    map<string,SplineInterpolation<Real_wp> > spls;

    // Vectors to store data given through nodal points input files
    map<string,Vector<Real_wp> > nodal_data;

    // The same for added fields (perturbations) 
    map<string,bool> storage_add;
    map<string,int> storage_add_type;
    map<string,Real_wp> cvals_add;
    map<string,SplineInterpolation<Real_wp> > spls_add;
    map<string,Vector<Real_wp> > nodal_data_add;

    // Constant coef multiplier
    map<string,Real_wp> adim_coefs;

};

#endif
