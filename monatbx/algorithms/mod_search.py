class my_search():
  def binary_search(self, data, item):
    """
    From an ordered list, find the item by looking at the left half and right half.
    """
    if len(data) == 0:
      #base case
      return False
    else:
      mid = len(data)//2
      if data[mid]==item:
        return True
      else:
        if item < data[mid]:
          return self.binary_search(data[:mid], item)
        else:
          return self.binary_search(data[mid+1:], item)

if __name__ == "__main__":
  ms = my_search()
  test_list = [0, 1, 2, 8, 13, 17, 19, 32, 42]
  item = 22
  print ms.binary_search(test_list, item)
