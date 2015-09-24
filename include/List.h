#ifndef LIST_H
#define LIST_H
///////////////////////////////////////////////////////////////////////////////
// List.h:  template list operations.
//
// Reference: Barton & Nackman, "Scientific & Engineering C++".
//
// ListNode was made an unnested class to avoid problems with cfront.
//
// $Id: List.h,v 8.1 2015/04/20 11:14:14 hmb Exp $
///////////////////////////////////////////////////////////////////////////////

template<class T> class ListIterator;
template<class T> class ListNode {
  public:
    ListNode(T x) : link(0), datum(x) {}
    ListNode<T>* link;
    T            datum;
};

template<class T> class List {
friend class ListIterator<T>;
private:
  ListNode<T>* head;
  ListNode<T>* tail;
  int          nel;

  List(const List<T>&);                // -- Prohibit, since not implemented. 
  List<T>& operator=(const List<T>&);  // -- Prohibit, since not implemented. 
public:
  List() : head(0), tail(0), nel(0) {}

  ~List() {
    while (head != 0) {
      ListNode<T>* p = head -> link;
      delete head;
      head    = p;
    }
    nel = 0;
  }

  void add (T x) {		// -- Unconditional insertion.
    if (head == 0) {
      head = new ListNode<T> (x);
      tail = head;
    } else
      tail = tail -> link = new ListNode<T> (x);
    nel++;
  }

  int xadd (T x) {		// -- Insert if non-replicative.
    register int   found = 0;
    register ListNode<T>* ptr;

    for (ptr = head; !found && ptr; ptr = ptr -> link) found = x == ptr->datum;
    if   (found)  {          return 0; }
    else          { add (x); return 1; }
  }

  T remove (T x) {		// -- Return datum of removed node.
    T            datum = 0;
    ListNode<T>* prev  = 0;
    ListNode<T>* curr  = head;
    while (curr != 0) {
      if (curr -> datum == x) {
	if (prev == 0) {
	  head  = curr -> link;
	  datum = curr -> datum;
	  delete curr;
	  curr = head;
	} else {
	  if (curr == tail) tail = prev;
	  prev -> link = curr -> link;
	  datum = curr -> datum;
	  delete curr;
	  curr = prev -> link;
	}
	nel--;
      } else {
	prev = curr;
	curr = curr -> link;
      }
    }
    return datum;
  }

  void clear  () { head = tail = 0; nel = 0; }

  int  length () const { return nel;                       }
  T    first  () const { return (nel) ? head -> datum : 0; }

};


template<class T>
class ListIterator {
public: 
  ListIterator (List<T>& list) : cur(list.head), top(list.head) { }
  
  int   more    () const  { return cur != 0;           }
  T     current () const  { return cur -> datum;       }
  void  next    ()        {        cur =  cur -> link; }
  void  reset   ()        {        cur =  top;         }

private:
  ListNode<T>* cur;
  ListNode<T>* top;
};

#endif
