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

    def oplus(self,X, Y, pk):
        Z = (X * Y) % (pk * pk)
        return Z

    def produitParConstante(self,X, y, pk):
        Z = pow(X, y, pk * pk)
        return Z

    def oppose(self,X, pk):
        Z = mod_inverse(X, pk * pk)
        return Z

class Alice:
    def __init__(self, xa=0,ya=0):
        self.pk, self.sk = paillier.genkeys(10)
        print("pk = ", self.pk)
        print("sk = ", self.sk)
        print()
        self.xa = xa
        self.ya = ya

    def encrypt(self, x):
        return paillier.encrypt(x, self.pk)

    def decrypt(self, c):
        return paillier.decrypt(c, self.pk, self.sk)

    def distance(self, res):
        resu = (self.decrypt(res) + pow(self.xa,2) + pow(self.ya,2))
        return resu

class Bob:
    def __init__(self, pk, xb=0, yb=0):
        self.pk = pk
        self.xb = xb
        self.yb = yb

    def distance(self, xa, ya):
        part1 = paillier.oplus(self.encrypt(pow(self.xb,2)),self.encrypt(pow(self.yb,2)),self.pk)
        somme1 = paillier.oplus(paillier.produitParConstante(xa,self.xb,self.pk)
                    , paillier.produitParConstante(ya,self.yb,self.pk),self.pk)
        part2 = paillier.oppose(paillier.produitParConstante(somme1,2,self.pk),self.pk)
        return paillier.oplus(part1,part2,self.pk)

    def encrypt(self, x):
        return paillier.encrypt(x, self.pk)

if __name__ == '__main__':
    paillier = Paillier()
    alice = Alice(1,1)
    bob = Bob(alice.pk,2,2)
    print(alice.distance(bob.distance(alice.encrypt(alice.xa), alice.encrypt(alice.ya))))

