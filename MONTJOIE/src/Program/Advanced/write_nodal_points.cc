#define MONTJOIE_WITH_TWO_DIM

#include "Harmonic/MontjoieHarmonic.hxx"

using namespace Montjoie;

template<class Dimension>
class InputVar : public InputDataProblem_Base
{
public :
  typedef typename Dimension::R_N R_N;
  Vector<string> mesh_data;
  Mesh<Dimension> mesh;
  MeshNumbering<Dimension> mesh_num;
  
  InputVar() : mesh_num(mesh) {}
  
  void SetInputData(const string& description_field,
                    const VectString& parameters)
  {
    mesh.SetInputData(description_field, parameters);
    mesh_num.SetInputData(description_field, parameters);

    if (!description_field.compare("FileMesh")) 
      mesh_data = parameters;
  }
    
  void WriteNodalPointsMesh(bool discont)
  {
    Vector<R_N> PosNodes; Vector<IVect> Nodle;
    this->mesh.GetNodesCurvedMesh(PosNodes, Nodle, this->mesh.GetGeometryOrder(), true);
    
    ofstream file_out("nodal_points.dat");
    file_out.precision(16);
    if (discont)
      {
        int nb_nodes_total = 0;
        for (int i = 0; i < this->mesh.GetNbElt(); i++)
          nb_nodes_total += Nodle(i).GetM();
        
        file_out << this->mesh.GetNbElt() << " " << nb_nodes_total << endl;
        
        // on ecrit les points pour chaque element
        for (int i = 0; i < this->mesh.GetNbElt(); i++)
          {
            file_out << Nodle(i).GetM() << '\n';
            for (int j = 0; j < Nodle(i).GetM()-1; j++)
              {
                PrintNoBrackets(file_out, PosNodes(Nodle(i)(j)));
                file_out << " ";
              }

            int jmax = Nodle(i).GetM()-1;
            PrintNoBrackets(file_out, PosNodes(Nodle(i)(jmax)));
            file_out << '\n';
          }

        file_out.close();
      }
    else
      {
        file_out << this->mesh.GetNbElt() << " " << PosNodes.GetM() << '\n';
    
        // on ecrit les points
        for (int i = 0; i < PosNodes.GetM(); i++)
          {
            PrintNoBrackets(file_out, PosNodes(i));
            file_out << '\n';
          }
        
        // on ecrit la numerotation pour chaque element
        for (int i = 0; i < this->mesh.GetNbElt(); i++)
          {
            file_out << Nodle(i).GetM() << '\n';
            int n = -1;
            for (int j = 0; j < Nodle(i).GetM()-1; j++)
              {
                n = Nodle(i)(j);
                file_out << n << " ";
              }
            
            n = Nodle(i)(Nodle(i).GetM()-1);
            file_out << n << '\n';
          }
        
        file_out.close();
      }
  }  
  
};

int main(int argc, char **argv) 
{
  InitMontjoie(argc, argv);
  
  if (argc > 1)
    {
      // data file
      string file_name_data;
      file_name_data = string(argv[1]);

      int order = -1;
      if (argc > 2)
        order = atoi(argv[2]);

      bool discont = false;
      if (argc > 3)
        discont = true;
      
      // we get the type of element selected by the user, and type of equation
      string type_element, type_equation;
      getElement_Equation(file_name_data, type_element, type_equation);

      int dim_N = 3;
      if ((type_element.find("TRIANGLE") == 0) || (type_element.find("QUADRANGLE") == 0))
	dim_N = 2;

      if (dim_N == 2)
        {
          cout << "Data file for a 2-D case" << endl;
          
          // data file is read
          InputVar<Dimension2> var;
          ReadInputFile(file_name_data, var);
          
          if (order != -1)
            var.mesh.SetGeometryOrder(order);
          
          var.mesh.ConstructMesh(0, var.mesh_data);
          
          var.WriteNodalPointsMesh(discont);
          var.mesh.Write("test.mesh");
        }
      else
        {
          cout << "Data file for a 3-D case" << endl;
          
          // data file is read
          InputVar<Dimension3> var;
          ReadInputFile(file_name_data, var);
          
          if (order != -1)
            var.mesh.SetGeometryOrder(order);
          
          var.mesh.ConstructMesh(0, var.mesh_data);
          
          var.WriteNodalPointsMesh(discont);
          var.mesh.Write("test.mesh");
        }
    }
  else
    {
      cout<<"This code needs a data file in argument"<<endl;
      cout<<"main.x nom_fichier"<<endl;
      cout<<"is a good syntax"<<endl;
    }

  cout<<"End of the program"<<endl; 

  return FinalizeMontjoie();
}

