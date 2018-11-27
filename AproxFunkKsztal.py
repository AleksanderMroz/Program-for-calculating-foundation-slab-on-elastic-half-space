def n0(a, b, x, y):
    Temp = (1 - x) * (1 - y) * (1 + x + y - 2 * x ** 2 - 2 * y ** 2)
    return Temp

def n1(a, b, x, y):
    Temp = b * (1 - x) * (1 - y) ** 2 * y
    return Temp

def n2(a, b, x, y):
    Temp = -a * (1 - y) * (1 - x) ** 2 * x
    return Temp

def n3(a, b, x, y):
    Temp = (1 - y) * (3 * x + y - 2 * x ** 2 - 2 * y ** 2) * x
    return Temp

def n4(a, b, x, y):
    Temp = b * (1 - y) ** 2 * x * y
    return Temp

def n5(a, b, x, y):
    Temp = a * (1 - x) * (1 - y) * x ** 2
    return Temp

def n6(a, b, x, y):
    Temp = (-1 + 3 * x + 3 * y - 2 * x ** 2 - 2 * y ** 2) * x * y
    return Temp

def n7(a, b, x, y):
    Temp = -b * (1 - y) * x * y ** 2
    return Temp

def n8(a, b, x, y):
    Temp = a * (1 - x) * x ** 2 * y
    return Temp

def n9(a, b, x, y):
    Temp = (1 - x) * (x + 3 * y - 2 * x ** 2 - 2 * y ** 2) * y
    return Temp

def n10(a, b, x, y):
    Temp = -b * (1 - x) * (1 - y) * y ** 2
    return Temp

def n11(a, b, x, y):
    Temp = -a * (1 - x) ** 2 * x * y
    return Temp
