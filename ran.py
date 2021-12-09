from numpy.random import choice
list_of_candidates = ['A','T', 'G', 'C']
probability_distribution = [0.4, 0.0, 0.6, 0.0]

for i in range(10):
    draw = choice(list_of_candidates, 1, p = probability_distribution)
    print(draw)


class Student:

    def __init__(self, name='NONE', student_id='0', courses=[]):
        self.name = name
        self.student_id = student_id
        self.courses = courses
    
    def __str__(self):
        return "Name: "+self.name+"\nID: "+self.student_id+"\nCourses: "+str(self.courses)
    
    def courseCount(self):
        return len(self.courses)

x = Student("Joe", "280785947", ['Comp202', 'Comp273', 'Math133'])
x.name = "Sam"
#print(x)
#print(x.courseCount())

ob = "hello" > 'bye'
print(ob)