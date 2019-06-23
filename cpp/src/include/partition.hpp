
#include <iostream>

#include <vector>

#ifndef __PARTITION_HPP
#define __PARTITION_HPP

/**********************************************
 * partition
 **********************************************/
/// Sparse set representation

#define NOTIN 0xfffffff

class partition {
private:
  /// values' indices
  std::vector<size_t> index_;
  std::vector<size_t> bag_;
  //@}

public:
  /*!@name Parameters*/
  //@{
  /// list of values
  std::vector<std::vector<int>> bag;

  /*!@name Constructors*/
  //@{
  explicit partition();

  void clear();
  void resize(const size_t n, const size_t m);

  int size();

  const std::vector<int> &operator[](const int i) const;

  /*!@name List Manipulation*/
  //@{
  void move(const int elt, const int to);
  void add_elt(const int elt, const int to);
  void swap_bag(const int a, const int b);
  void remove_bag(const int a);
  int bag_of(const int elt);
  bool contain(const int elt, const int a);
  //@}

  /*!@name Miscellaneous*/
  //@{
  std::ostream &display(std::ostream &os) const;
};

std::ostream &operator<<(std::ostream &os, const partition &x);

#endif // __PARTITION_HPP
