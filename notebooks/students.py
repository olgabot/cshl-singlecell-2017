import random

s = '''Hanrui Zhang
Kenny Yu
Jessica Tsai
Kristof Torkenczy
Rhiana Simon
Aarushi Sharma
Bjoern Papke
Tae Hyun Kim
Aditya Kashyap
Yiran Hou
Stefanie Grunwald
Giovanni Diaz
Andre Deslauriers
Naiyan Chen
Ashley Brandebura
Timothy Blosser
Hongtian “Stanley” Yang
Wendy Aquino Nunez'''
student_names = s.split('\n')


def random_student():
    return random.choice(student_names)