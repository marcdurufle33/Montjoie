#ifndef MONTJOIE_FILE_TIMER_HXX

#ifdef MONTJOIE_WITH_REAL_TIMING
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <map>

namespace Montjoie
{

  //! base class to do timings
  class VirtualTimer
  {
  protected :
    Vector<double> elapsed_time; //!< elapsed time for each timer
    Vector<bool> jeton_libre; //!< free tokens
    
    virtual void ResizeClock(int n) = 0;
    
  public:
    enum {ALL, COMM, PROD, STIFFNESS, FLUX, MASS, EXTRAPOL, SCHEME, PML,
          OUTPUT, JACOBIAN, SOLVE, FACTO};
    
    VirtualTimer();
    virtual ~VirtualTimer();
    
    int GetNumOfTimers() const;
    int GetNumber();
    void ReleaseNumber(int i);
    void ReserveNumber(int i);
    void Reset(int i);
    
    double GetSeconds(int i);
    
    virtual void Start(int i) = 0;
    virtual void Stop(int i) = 0;
    
  };

  
  //! Timer using clock function
  /*!
    In order to launch a timing, you need to call first GetNumber
    Then, to start the time counting, you call Start, and to stop it, you call Stop
    If you recall Start, the time is incremented (there is no reset)
    To get the elapsed time, you call GetSeconds
  */
  class BasicTimer : public VirtualTimer
  {
  protected:
    Vector<clock_t> value_clock; //!< clock time for each timer that started

    virtual void ResizeClock(int n);
    
  public :
    
    BasicTimer();
    
    void Start(int i);
    void Stop(int i);
    
  };


#ifdef MONTJOIE_WITH_REAL_TIMING
  //! Timer using getrusage
  /*!
    In order to launch a timing, you need to call first GetNumber
    Then, to start the time counting, you call Start, and to stop it, you call Stop
    If you recall Start, the time is incremented (there is no reset)
    To get the elapsed time, you call GetSeconds
  */
  class AccurateTimer : public VirtualTimer
  {
  protected :
    struct timeval tim;        
    struct rusage ru;
    Vector<double> value_clock; //!< clock time for each timer that started

    virtual void ResizeClock(int n);
    
  public :
    
    AccurateTimer();

    void Start(int i);
    void Stop(int i);
    
  };  


  //! Timer using clock_gettime
  /*!
    In order to launch a timing, you need to call first GetNumber
    Then, to start the time counting, you call Start, and to stop it, you call Stop
    If you recall Start, the time is incremented (there is no reset)
    To get the elapsed time, you call GetSeconds
  */
  class RealTimer : public VirtualTimer
  {
  protected:
    Vector<struct timespec> value_clock; //!< clock time for each timer that started
    
    virtual void ResizeClock(int n);
    
  public :
    
    RealTimer();
    
    void Start(int i);
    void Stop(int i);
    
  };
#endif

  
  //! Timer used in Montjoie
  class MontjoieTimer
  {
  protected :
#ifdef SELDON_WITH_MPI
    //! MPI communicator
    MPI_Comm comm_;
#endif
    
    //! pointer to a virtual timer
    VirtualTimer* timer;
    
    //! list of timers defined with strings instead of integers as other timers
    map<string, int> liste_chrono;
    
    //! list of messages associated with timers
    Vector<string> liste_message;
    
  public :
    // available timers
    enum {BASIC_TIMER, ACCURATE_TIMER, REAL_TIMER};
    
    //! default timer when creating an object MontjoieTimer
    static int default_timer;

    MontjoieTimer();
    ~MontjoieTimer();
    
    void SetTimer(int type);
    
#ifdef SELDON_WITH_MPI
    void SetCommunicator(const MPI_Comm& comm);
#endif

    int GetNumber();
    void ReleaseNumber(int i);
    void ReserveNumber(int i);

    void Start(int);
    void Stop(int);
    void Reset(int);
    
    void Start(const string&);
    void Stop(const string&);
    void Reset(const string&);

    double GetSeconds(int i);
    double GetSeconds(const string& name);

    void SetMessage(const string& name, const string& message);
    void ReserveName(const string& name);
    void ReleaseName(const string& name);
    bool NameExists(const string& name) const;

    void GetGlobalSeconds(int i, double& dt_loc, double& dt_sum,
                          double& dt_min, double& dt_max);
    
    void DisplayTime(int i, const string& message);

    void GetGlobalSeconds(const string&, double& dt_loc, double& dt_sum,
                          double& dt_min, double& dt_max);
    
    void DisplayTime(const string&, const string& message);
    
    void DisplayTime(int i);
    void DisplayTime(const string& name);
    void DisplayAll();
    
  };

  
  //! global timer
  extern MontjoieTimer glob_chrono; 
    
}

#define MONTJOIE_FILE_TIMER_HXX
#endif
