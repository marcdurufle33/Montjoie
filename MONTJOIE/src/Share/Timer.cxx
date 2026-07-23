#ifndef MONTJOIE_FILE_TIMER_CXX

namespace Montjoie
{

  /****************
   * VirtualTimer *
   ****************/

  
  //! default constructor
  VirtualTimer::VirtualTimer()
    : elapsed_time(20), jeton_libre(20)
  {
    jeton_libre.Fill(true);
    elapsed_time.Fill(0.0);
    
    jeton_libre(ALL) = false;
    jeton_libre(COMM) = false;
    jeton_libre(PROD) = false;
    jeton_libre(STIFFNESS) = false;
    jeton_libre(FLUX) = false;
    jeton_libre(MASS) = false;
    jeton_libre(EXTRAPOL) = false;
    jeton_libre(SCHEME) = false;
    jeton_libre(PML) = false;
    jeton_libre(OUTPUT) = false;
    jeton_libre(JACOBIAN) = false;
    jeton_libre(SOLVE) = false;
    jeton_libre(FACTO) = false;
  }
  

  //! Destructor
  VirtualTimer::~VirtualTimer()
  {
  }
  
  
  //! returns the number of timers that are stored
  int VirtualTimer::GetNumOfTimers() const
  {
    return jeton_libre.GetM();
  }

  
  //! returns an integer identifying the new timer
  int VirtualTimer::GetNumber()
  {
    for (int i = 0; i < jeton_libre.GetM(); i++)
      if (jeton_libre(i))
	{
	  jeton_libre(i) = false;
	  return i;
	}
    
    int taille = jeton_libre.GetM()+1;
    jeton_libre.PushBack(false);
    elapsed_time.PushBack(0.0);
    ResizeClock(taille);
    return taille-1;
  }
  
  
  //! we don't want to use anymore the timer identified "i"
  void VirtualTimer::ReleaseNumber(int i)
  {
    if (i < jeton_libre.GetM())
      {
	jeton_libre(i) = true;
	elapsed_time(i) = 0.0;
      }
  }
  
  
  //! The timer i is already used
  void VirtualTimer::ReserveNumber(int i)
  {
    if (i >= jeton_libre.GetM())
      {
	int taille = jeton_libre.GetM();
	jeton_libre.Resize(i+1);
	elapsed_time.Resize(i+1);
	ResizeClock(i+1);
	for (int j = taille; j < i; j++)
	  {
	    jeton_libre(j) = true;
	    elapsed_time(j) = 0.0;
	  }
      }
    
    jeton_libre(i) = false;
    elapsed_time(i) = 0.0;
  }
  
  
  //! we want to reset the timer i (elapsed time is set to 0)
  void VirtualTimer::Reset(int i)
  {
    if (i < jeton_libre.GetM())
      elapsed_time(i) = 0.0;
  }


  //! returns elapsed time of the timer i
  /*!
    The result is expressed in seconds
  */
  double VirtualTimer::GetSeconds(int i)
  {
    double temps(0);
    if (i < jeton_libre.GetM())
      temps = elapsed_time(i);
    
    return temps;
  }


  /**************
   * BasicTimer *
   **************/

  
  //! default constructor
  BasicTimer::BasicTimer() : VirtualTimer(), value_clock(20)
  {
  }
  
  
  //! modifies size of value_clock
  void BasicTimer::ResizeClock(int n)
  {
    value_clock.Resize(n);
  }
  
  
  //! start of timing for the timer i
  void BasicTimer::Start(int i)
  {
    if (i < this->jeton_libre.GetM())
      this->value_clock(i) = clock();
  }
  
  
  //! end of timing of the timer i
  void BasicTimer::Stop(int i)
  {
    if (i < this->jeton_libre.GetM())
      {
	clock_t temps  = clock();
	if ((temps < 0) && (this->value_clock(i) > 0))
	  {
	    // case with an overflow of clock_t type
	    // it appears generally on a 32-bits machine
	    clock_t diff = 1;
	    diff += numeric_limits<clock_t>::max() - this->value_clock(i);
	    diff += (temps - numeric_limits<clock_t>::min());
	    this->elapsed_time(i) += double(diff)/CLOCKS_PER_SEC;
	  }
	else
	  this->elapsed_time(i) += double(temps-this->value_clock(i))/CLOCKS_PER_SEC;
      }
  }
  

  /*****************
   * AccurateTimer *
   *****************/
  

#ifdef MONTJOIE_WITH_REAL_TIMING
  //! default constructor
  AccurateTimer::AccurateTimer() : VirtualTimer(), value_clock(20)
  {
  }
  

  //! modifies size of value_clock
  void AccurateTimer::ResizeClock(int n)
  {
    value_clock.Resize(n);
  }

  
  //! start of timing of the timer i
  void AccurateTimer::Start(int i)
  {
    if (i < this->jeton_libre.GetM())
      {
	getrusage(RUSAGE_SELF, &ru);        
        tim = ru.ru_utime;        
        this->value_clock(i) = double(tim.tv_sec) + double(tim.tv_usec)/ 1000000.0;
      }
  }
  
  
  //! end of timing of the timer i
  void AccurateTimer::Stop(int i)
  {
    if (i < this->jeton_libre.GetM())
      {
	getrusage(RUSAGE_SELF, &ru);        
        tim = ru.ru_utime;        
        double temps = double(tim.tv_sec) + double(tim.tv_usec)/ 1000000.0;
	this->elapsed_time(i) += temps - this->value_clock(i);
      }
  }  


  /*************
   * RealTimer *
   *************/  

  
  //! default constructor
  RealTimer::RealTimer() : VirtualTimer(), value_clock(20)
  {
  }
  
  
  //! modifies size of value_clock
  void RealTimer::ResizeClock(int n)
  {
    value_clock.Resize(n);
  }


  //! start of timing for the timer i
  void RealTimer::Start(int i)
  {
    if (i < jeton_libre.GetM())
      clock_gettime(CLOCK_MONOTONIC, &value_clock(i));
  }
  
  
  //! end of timing for the timer i
  void RealTimer::Stop(int i)
  {
    if (i < jeton_libre.GetM())
      {
        time_t sec = value_clock(i).tv_sec;
        long nsec = value_clock(i).tv_nsec;
        clock_gettime(CLOCK_MONOTONIC, &value_clock(i));
	elapsed_time(i) += value_clock(i).tv_sec - sec;
        elapsed_time(i) += (value_clock(i).tv_nsec - nsec)/1.0e9;
      }
  }
#endif
  

  /*****************
   * MontjoieTimer *
   *****************/
  
#ifdef MONTJOIE_WITH_REAL_TIMING
  int MontjoieTimer::default_timer(MontjoieTimer::REAL_TIMER);
#else
  int MontjoieTimer::default_timer(MontjoieTimer::BASIC_TIMER);
#endif
  
  //! default constructor
  MontjoieTimer::MontjoieTimer() 
  {
    timer = NULL;
    SetTimer(default_timer);
    
    int nb_timers = timer->GetNumOfTimers();
    liste_message.Reallocate(nb_timers);
    
#ifdef SELDON_WITH_MPI
    comm_ = MPI_COMM_SELF;
#endif
  }
  
  
  //! Destructor
  MontjoieTimer::~MontjoieTimer()
  {
    if (timer != NULL)
      delete timer;
  }


  //! changes the timing system to use
  void MontjoieTimer::SetTimer(int type)
  {
    if (timer != NULL)
      delete timer;
    
    if (type == BASIC_TIMER)
      timer = new BasicTimer();

#ifdef MONTJOIE_WITH_REAL_TIMING
    if (type == ACCURATE_TIMER)
      timer = new AccurateTimer();

    if (type == REAL_TIMER)
      timer = new RealTimer();
#endif
    
    if (timer == NULL)
      timer = new BasicTimer();
  }


#ifdef SELDON_WITH_MPI
  //! specifies the processors that share the timer
  void MontjoieTimer::SetCommunicator(const MPI_Comm& comm)
  {
    comm_ = comm;
  }
#endif
  

  //! returns an integer identifying the new timer
  int MontjoieTimer::GetNumber()
  {
    int n = timer->GetNumber();
    
    if (timer->GetNumOfTimers() > liste_message.GetM())
      liste_message.Resize(timer->GetNumOfTimers());
    
    return n;
  }

  
  //! releases the timer i
  void MontjoieTimer::ReleaseNumber(int i)
  {
    timer->ReleaseNumber(i);
  }

  
  //! reserves the timer i
  void MontjoieTimer::ReserveNumber(int i)
  {
    timer->ReserveNumber(i);
  }
  
  
  //! Timer i is started or restarted
  void MontjoieTimer::Start(int i)
  {
    timer->Start(i);
  }
  
  
  //! Timer i is stopped
  void MontjoieTimer::Stop(int i)
  {
    timer->Stop(i);
  }
  
  
  //! Timer i is resetted 
  void MontjoieTimer::Reset(int i)
  {
    timer->Reset(i);
  }
    
  
  //! Timer name is started
  void MontjoieTimer::Start(const string& name)
  {
    // adding a timer if name does not exist
    map<string, int>::iterator it = liste_chrono.find(name);
    int n = -1;
    if (it == liste_chrono.end())
      {
        n = GetNumber();
        liste_chrono[name] = n;
      }
    else
      n = it->second;
    
    Start(n);
  }
  
  
  //! Timer name is stopped
  void MontjoieTimer::Stop(const string& name)
  {
    map<string, int>::iterator it = liste_chrono.find(name);
    int n = -1;
    if (it == liste_chrono.end())
      {
        cout << "The timer " << name << " does not exist " << endl;
        cout << "You can't stop it " << endl;
        abort();
      }
    else
      n = it->second;
    
    Stop(n);
  }
  
  
  //! Timer name is resetted
  void MontjoieTimer::Reset(const string& name)
  {
    // adding a timer if name does not exist
    map<string, int>::iterator it = liste_chrono.find(name);
    int n = -1;
    if (it == liste_chrono.end())
      {
        n = GetNumber();
        liste_chrono[name] = n;
      }
    else
      n = it->second;
    
    Reset(n);
  }

  
  //! returns elapsed time in seconds for timer i
  double MontjoieTimer::GetSeconds(int i)
  {
    return timer->GetSeconds(i);
  }


  //! returns elapsed time in seconds for timer name
  double MontjoieTimer::GetSeconds(const string& name)
  {
    map<string, int>::iterator it = liste_chrono.find(name);
    int n = -1;
    if (it == liste_chrono.end())
      {
        cout << "The timer " << name << " does not exist " << endl;
        abort();
      }
    else
      n = it->second;
    
    return timer->GetSeconds(n);
  }
  
  
  //! sets a message for the timer name
  /*!
    The message will be displayed if DisplayTime or DisplayAll is called
  */
  void MontjoieTimer::SetMessage(const string& name, const string& message)
  {
    // adding a timer if name does not exist
    map<string, int>::iterator it = liste_chrono.find(name);
    int n = -1;
    if (it == liste_chrono.end())
      {
        n = GetNumber();
        liste_chrono[name] = n;
      }
    else
      n = it->second;
    
    liste_message(n) = message;
  }
  
  
  //! adds a timer whose name is name
  void MontjoieTimer::ReserveName(const string& name)
  {
    map<string, int>::iterator it = liste_chrono.find(name);
    if (it == liste_chrono.end())
      {
        int n = GetNumber();
        liste_chrono[name] = n;
      }
    else 
      {
        cout << "The timer " << name << " is already used " << endl;
        abort();
      }    
  }
  
  
  //! clears the timer whose name is name
  void MontjoieTimer::ReleaseName(const string& name)
  {
    map<string, int>::iterator it = liste_chrono.find(name);
    if (it == liste_chrono.end())
      {
        cout << "The timer " << name << " does not exist" << endl;
        abort();
      }
    else
      liste_chrono.erase(it);
  }


  //! returns true if the chrono name exists
  bool MontjoieTimer::NameExists(const string& name) const
  {
    map<string, int>::const_iterator it = liste_chrono.find(name);
    if (it == liste_chrono.end())
      return false;
    
    return true;
  }
  
  
  //! fills local time, sum of of time, minimal and maximal time among processors
  void MontjoieTimer::GetGlobalSeconds(int i, double& dt_loc, double& dt_sum,
				       double& dt_min, double& dt_max)
  {
    dt_loc = timer->GetSeconds(i);
    
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    if (nb_proc > 1)
      {  
        dt_sum = 0;
	dt_max = 0;
        MPI_Reduce(&dt_loc, &dt_sum, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&dt_loc, &dt_min, 1, MPI_DOUBLE, MPI_MIN, 0, comm);
        MPI_Reduce(&dt_loc, &dt_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
      }
    else
      {
        dt_sum = dt_loc;
        dt_min = dt_loc;
        dt_max = dt_loc;
      }
#endif
  }
    

  //! Displays times associated with timer i
  void MontjoieTimer::DisplayTime(int i, const string& message)
  { 
    double t_loc, t_sum, t_min, t_max;
    GetGlobalSeconds(i, t_loc, t_sum, t_min, t_max);
    
#ifdef SELDON_WITH_MPI
    MPI_Comm& comm = comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    if (nb_proc > 1)
      {
	if (rank_proc == 0)
	  {
	    cout << "Global time to " << message << " : " << t_sum << endl;
	    cout << "Maximal time to " << message << " : " << t_max << endl;
	    
	     cout << "Minimal time to " << message << " : " << t_min << endl;
            cout << endl;
	  }
      }
    else
#endif
      cout << " Time to " << message << " : " << t_loc << endl;    
  }


  //! fills local time, sum of of time, minimal and maximal time among processors
  void MontjoieTimer::GetGlobalSeconds(const string& name, double& dt_loc, double& dt_sum,
				       double& dt_min, double& dt_max)
  {
    map<string, int>::iterator it = liste_chrono.find(name);
    int n = -1;
    if (it == liste_chrono.end())
      {
        cout << "The timer " << name << " does not exist " << endl;
        abort();
      }
    else
      n = it->second;

    GetGlobalSeconds(n, dt_loc, dt_sum, dt_min, dt_max);
  }
    

  //! Displays times associated with timer name
  void MontjoieTimer::DisplayTime(const string& name, const string& message)
  {
    map<string, int>::iterator it = liste_chrono.find(name);
    int n = -1;
    if (it == liste_chrono.end())
      {
        cout << "The timer " << name << " does not exist " << endl;
        abort();
      }
    else
      n = it->second;
    
    DisplayTime(n, message);
  }


  //! Displays times associated with timer i
  void MontjoieTimer::DisplayTime(int i)
  {
    DisplayTime(i, liste_message(i));
  }
  
  
  //! Displays times associated with timer name
  void MontjoieTimer::DisplayTime(const string& name)
  {
    map<string, int>::iterator it = liste_chrono.find(name);
    int n = -1;
    if (it == liste_chrono.end())
      {
        cout << "The timer " << name << " does not exist " << endl;
        abort();
      }
    else
      n = it->second;
 
    DisplayTime(n, liste_message(n));
  }
  
  
  //! Displays all timers whose name is known (timers associated with integers are not displayed)
  void MontjoieTimer::DisplayAll()
  {
    map<string, int>::iterator it;
    for (it = liste_chrono.begin(); it != liste_chrono.end(); ++it)
      DisplayTime(it->second, liste_message(it->second));    
  }


  // global variable
  MontjoieTimer glob_chrono; 
    
}

#define MONTJOIE_FILE_TIMER_CXX
#endif
