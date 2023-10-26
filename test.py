
class One:

    def __init__(self):
        self.thing = 7

class Two:

    def __init__(self, one):
        self.thing = one.thing
        self.one = one

one = One()

two = Two(one)

print(two.thing)
one.thing = 5
print(two.thing)
print(two.one.thing)