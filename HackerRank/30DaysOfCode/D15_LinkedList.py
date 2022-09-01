class Node:
    def __init__(self, data):
        self.data = data
        self.next = None


class Solution:
    def display(self, head):
        current = head
        while current:
            print(current.data, end=' ')
            current = current.next

    def insert(self, head, data):
        # Initiate list
        if head == None:
            return Node(data)
        itr = head
        while itr.next:
            itr = itr.next
        itr.next = Node(data)
        return head


mylist = Solution()
T = 4
head = None
for i in range(T):
    data = int(i)
    head = mylist.insert(head, data)
    # print(head)
mylist.display(head)
