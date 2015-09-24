#ifndef STACK_H
#define STACK_H
///////////////////////////////////////////////////////////////////////////////
// Stack.h: templated operations for LIFO stack.
//
// Summary: creator, destructor, push, pop, depth.
//
// $Id: Stack.h,v 8.1 2015/04/20 11:14:14 hmb Exp $
///////////////////////////////////////////////////////////////////////////////

template<class T> class StackNode {
public:
  StackNode (T x) : next (0), datum (x) {}
  StackNode<T>* next;
  T     datum;
};

template<class T>
class Stack {

private:
  StackNode<T>* top;
  int   stack_depth;

  Stack(const Stack<T>&);                // Prohibit, since not implemented.
  Stack<T>& operator=(const Stack<T>&);  // Prohibit, since not implemented.

public:
  Stack () : top (0), stack_depth (0) {}
  ~Stack() {
    while (stack_depth--) {
      StackNode<T>* p = top -> next;
      delete top;
      top = p;
    }
  }

  void push(T x) {
    if (stack_depth) {
      StackNode<T>* p = new StackNode<T> (x);
      p -> next = top;
      top       = p;
    } else {
      top = new StackNode<T> (x);
    }
    stack_depth++;
  }
  
  T pop() {
    if (stack_depth) {
      StackNode<T>* p = top;
      T     value = top -> datum;
      top         = top -> next;
      delete p;
      stack_depth--;
      return value;
    } else {
      return 0;
    }
  }

  int depth () const { return stack_depth; }
};

#endif
