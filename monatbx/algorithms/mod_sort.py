class my_sort():
  def __init__(self):
    dummy = 0

  def insertion_sort(self, data):
    """
    O(N^2)
    takes a 1-d array from the second position,
    find the prior position
    where the value is smaller and insert the key
    """
    n = len(data)
    for i in range(1, n):
      key = data[i]
      j = i
      print
      print i, j, key, data
      while j >0 and data[j-1] > key:
        data[j] = data[j-1]
        print i, j, key, data
        j = j-1
      data[j] = key
      print 'sorted data=', data
    return data

  def bubble_sort(self, data):
    """
    O(N^2)
    Do n * (n-1) pass, and bubble down a smaller item between
    an item pair.
    """
    for i in range(len(data)-1, 0, -1):
      for j in range(i):
        if data[j] > data[j+1]:
          tmp = data[j]
          data[j] = data[j+1]
          data[j+1] = tmp
      print "pass", i, data

  def selection_sort(self, data):
    """
    O(N^2)
    Improves from bubble sort, by selecting the largest item in each pass
    and swap the value with pass no position and its position.
    """
    for i in range(len(data)-1, 0, -1):
      i_max = 0
      for j in range(1, i+1):
        if data[j] > data[i_max]:
          i_max = j
      tmp = data[i]
      data[i] = data[i_max]
      data[i_max] = tmp
      print "pass", i, data


  def merge_sort(self, data):
    """
    O(N log N)
    recursively splitting array in half until
    basecase, that is the array on the left and right has a length of 1
    when that happens, compare the left and the right array last element,
    copy the small value to the array
    """
    print "Splitting", data
    if len(data) > 1:
      mid = len(data)// 2
      left_half = data[:mid]
      right_half = data[mid:]
      self.merge_sort(left_half)
      self.merge_sort(right_half)
      i,j,k = (0,0,0)
      while i< len(left_half) and j<len(right_half):
        if left_half[i] < right_half[j]:
          data[k] = left_half[i]
          i += 1
        else:
          data[k] = right_half[j]
          j += 1
        k += 1
      #take care of the unbalance
      while i < len(left_half):
        data[k] = left_half[i]
        i += 1
        k += 1
      while j < len(right_half):
        data[k] = right_half[j]
        j += 1
        k += 1
    print "merging", data

if __name__ == "__main__":
  ms = my_sort()
  data = [8,6,4,10,9,5,3,2]
  print 'Insertion Sort'
  i_sorted = ms.insertion_sort(data[:])
  print 'Bubble Sort'
  ms.bubble_sort(data[:])
  print 'Selection Sort'
  ms.selection_sort(data[:])
  print 'Merge Sort'
  ms.merge_sort(data[:])
