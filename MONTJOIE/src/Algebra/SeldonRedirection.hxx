#ifndef SELDON_FILE_SELDON_REDIRECTION_HXX

namespace Seldon
{

  template<class Storage>
  class SeldonDefaultAllocator<Storage, short int>
  {
  public:
    typedef MallocAlloc<short int> allocator;
  };  

  template<class Storage>
  class SeldonDefaultAllocator<Storage, short unsigned>
  {
  public:
    typedef MallocAlloc<short unsigned> allocator;
  };  

#ifdef SELDON_WITH_FLOAT128  
  template<class Storage>
  class SeldonDefaultAllocator<Storage, __float128>
  {
  public:
    typedef MallocAlloc<__float128> allocator;
  };

  template<class Storage>
  class SeldonDefaultAllocator<Storage, complex<__float128> >
  {
  public:
    typedef MallocAlloc<complex<__float128> > allocator;    
  };  
#endif
  
};

#define SELDON_FILE_SELDON_REDIRECTION_HXX
#endif

