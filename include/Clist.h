/*****************************************************************************
 * Clist.h:  templated operations for circularly-linked lists.               *
 *****************************************************************************/

/* $Id: Clist.h,v 8.1 2015/04/20 11:14:14 hmb Exp $ */

#ifndef ClistH
#define ClistH

template<class T>
class CircularListIterator;


template<class T>
class CircularList {
public:
  CircularList() : head(0), num_elmts(0) {}
  ~CircularList();


  void add    (T x);
  void remove (T x);

  int  size   () const { return num_elmts; }

  friend class CircularListIterator<T>;

private:
  class Node {
  public:
    Node(T x) : next(0), datum(x) {}
    Node* next;
    T     datum;
  };

  Node* head;
  Node* tail;

  int   num_elmts;

  CircularList(const CircularList<T>&);                     // Prohibit.
  CircularList<T>& operator=(const CircularList<T>&);       // Prohibit.
};






template<class T>
class CircularListIterator {
public:
  CircularListIterator(const CircularList<T>& list) : cur(list.head) {}
  
  int     more    () const { return cur != 0;     }
  T       current () const { return cur -> datum; }
  void    advance ()       { cur =  cur -> next;  }

private:
  CircularList<T>::Node* cur;
};





template<class T>
CircularList<T>::~CircularList() {
  while (num_elmts > 0) {
    Node* p = head -> next;
    delete head;
    head = p;
    num_elmts--;
  }
}


template<class T>
void CircularList<T>::add(T x) {
  if (head == 0) {
    head = tail  = new Node (x);
    head -> next = head;
  } else {
    tail = tail -> next = new Node(x);
    tail -> next = head;
  }
  num_elmts++;
}


template<class T>
void CircularList<T>::remove(T x) {
  Node* prev = 0;
  Node* cur  = head;
  int   n    = num_elmts;


  while (n > 0) {
    if (cur -> datum == x) {
      if (prev == 0) {
	if (head -> next == head) {
	  delete head;
	  head = tail = 0;
	  num_elmts = 0;
	  break;
	} else {
	  head = cur -> next;
	  tail -> next = head;
	  delete cur;
	  cur = head;
	}
      } else {
        prev->next = cur->next;
        delete cur;
        cur = prev->next;
      }
      num_elmts--;
    } else {
      prev = cur;
      cur = cur->next;
    }
    n--;
  }
}


#endif
