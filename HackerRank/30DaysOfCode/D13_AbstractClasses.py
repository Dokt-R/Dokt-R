# HackerRank 30 Days of Code
#  Abstract Classes
# Make a Class that
# Inherits from Book
# Has 3 parameters: title, author, price
# Uses the display method to print the information

from abc import ABCMeta, abstractmethod


class Book(object, metaclass=ABCMeta):
    def __init__(self, title, author):
        self.title = title
        self.author = author

    @abstractmethod
    def display(): pass

# Write MyBook class


class MyBook(Book):
    def __init__(self, title, author, price):
        super().__init__(title, author)
        self.price = price

    @abstractmethod
    def display(self):
        print('Title:', self.title, '\nAuthor:',
              self.author, '\nPrice:', self.price)


title = input()
author = input()
price = int(input())
new_novel = MyBook(title, author, price)
new_novel.display()
