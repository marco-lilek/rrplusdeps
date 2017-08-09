#ifndef _VECTORS_H
#define _VECTORS_H

#include <stdlib.h>
#include <assert.h>
#include <iostream>

using namespace std;

// templated vectors class, partially based on Budd,
//                         Classic Data Structures in C++
// written 11/5/93, modified 3/23/94
// changed on 3/9/95 to make all methods inline (defined in class decl)
//
// for a vectors of Items use vectors<Item>, e.g., vectors <int> intvector;
//                           
//                           note: Item must have a default constructor
//
// constructors:
//   vectors()                   -- default, vectors of size 0 (no entries)
//   vectors(int size)           -- vectors with size entries
//   vectors(const vectors & vec) -- copy constructor
//    
//   int Length()                -- returns size of vectors (capacity)
//   void SetSize(int newSize)   -- resizes the vectors to newSize elements
//                                  (can result in losing elements if
//                                   new size < old size)
//   void resize(int newSize)     -- synonym for SetSize
//
//   operator =                   -- assignment operator works properly
//   operator []                  -- indexes both const and non-const vectorss
//    
//

class vectors
{
  public:
    vectors();                               // default constructor 0 elts
    vectors(int size);                       // specify size of vectors
    vectors(int size, double fillValue);     // specify size and fill value
    vectors(const vectors & vec);            // copy constructor
    ~vectors ();                             // free new'd storage
    vectors & operator = (const vectors & vec); // overload assignment
    int Length() const;                      // capacity of vectors
    int length() const;
    void Fill(double fillValue);
    void SetSize(int newSize);               // change size dynamically
    void resize(int newSize);
    double & operator [] (int index);
    const double & operator [] (int index) const; // const index 

  private:
    double * myList;   // the array of items
    int myLength;   // # things in vectors (array), 0,1,...,(myLength-1)
};

#endif                          // _vectors_H not defined


