#include "Share/MontjoieCommon.hxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  if (argc != 3)
    {
      cout << "Usage ./test.x n1 n2" << endl;
      abort();
    }

  int nb_iter1 = atoi(argv[1]);
  int nb_iter2 = atoi(argv[2]);

  glob_chrono.Start(VirtualTimer::ALL);
  
  for (int k = 0; k < nb_iter1; k++)
    {
      // Dense vectors
      //Vector<double, VectFull, NewAlloc<double> > num;
      Vector<double> num;
      //vector<double> num;
      
      // testing Resize
      for (int p = 0; p < nb_iter2; p++)
        {
          int taille = rand()%100;
	  //num.resize(taille);
	  //for (unsigned j = 0; j < num.size(); j++)
	  //num[j] = j;
	  
          num.Resize(taille);
          //num.Clear();
          //num.Reallocate(taille);
          num.Fill();
	}

      // testing PushBack
      //num.Reallocate(rand()%30);
      /* num.resize(rand()%30);
      for (int p = 0; p < nb_iter2; p++)
      num.push_back(rand());*/
	//num.PushBack(rand());

      // Sparse vectors
      //Vector<double, VectSparse> val;
      
      // testing AddInteraction
      /*for (int p = 0; p < nb_iter2; p++)
        {
          int n = rand()%1000;
          double x = 1e-9*rand();
          val.AddInteraction(n, x);
          } */

      // testing AddInteractionRow
      /* Vector<int> num_interac(10);
      Vector<double> val_interac(10);
      for (int p = 0; p < nb_iter2; p++)
        {
          for (int q = 0; q < num_interac.GetM(); q++)
            {
              num_interac(q) = rand()%1000;
              val_interac(q) = 1e-9*rand();
            }
              
          val.AddInteractionRow(num_interac.GetM(), num_interac, val_interac);
	  }*/
    }
  
  glob_chrono.Stop(VirtualTimer::ALL);
  cout << "Elapsed time = " << glob_chrono.GetSeconds(VirtualTimer::ALL) << endl;


  return FinalizeMontjoie();
}
