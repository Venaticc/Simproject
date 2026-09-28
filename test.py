
class Function:
    def __init__(self, func, name):
        self.func = func
        self.name = name

    def use(self, *args):
        return self.func(*args)
    def getname(self):
        return self.name
    def getfunc(self):
        return self.func

def func1(x):
    return x**2

Function(func1, 'func1')
