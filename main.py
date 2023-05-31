from sympy import *
from random import *

class Paillier:
    def getprime(self,k):
        p = randprime(2**(k-1), 2**k)
        return p


    def genkeys(self,k):
        p = self.getprime(k)
        q = self.getprime(k)
        while (p == q):
            q = self.getprime(k)

        N = int(p * q)
        Phi = (p-1)*(q-1)
        return [N, mod_inverse(N, Phi)]


    def encrypt(self,m, pk):
        r = randint(1, pk-1)
        N2 = pk*pk
        c = ((1+m*pk) * pow(r, pk, N2)) % N2
        return int(c)


    def decrypt(self,c, pk, sk):
        N2 = pk*pk
        r = pow(c, sk, pk)
        s = mod_inverse(r, pk)
        m = ((c * pow(s, pk, N2)) % N2 - 1)//pk

        return int(m)

class Alice:
    def __init__(self, paillier,xa=0,ya=0):
        self.paillier = paillier
        self.pk, self.sk = self.paillier.genkeys(10)
        print("pk = ", self.pk)
        print("sk = ", self.sk)
        print()
        self.xa = xa
        self.ya = ya

    def encrypt(self, x):
        return self.paillier.encrypt(x, self.pk)

    def decrypt(self, c):
        return self.paillier.decrypt(c, self.pk, self.sk)

    def distance(self, res):
        resu = (self.decrypt(res) + self.xa**2 + self.ya**2)
        return resu

class Bob:
    def __init__(self, paillier, pk, xb=0, yb=0):
        self.pk = pk
        self.xb = xb
        self.yb = yb
        self.paillier = paillier

    def distance(self, xa, ya):
        return self.encrypt(self.xb**2 + self.yb**2 - 2*(xa * self.xb + ya * self.yb))

    def encrypt(self, x):
        return self.paillier.encrypt(x, self.pk)

if __name__ == '__main__':
    paillier = Paillier()
    alice = Alice(paillier,1,1)
    bob = Bob(paillier,alice.pk,2,2)
    print(alice.distance(bob.distance(alice.encrypt(alice.xa), alice.encrypt(alice.ya))))

