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
        Phi = int(N - p - q + 1)
        return [N, mod_inverse(N, Phi)]

    def encrypt(self,m, pk):
        if (m < 0):
            m = m + pk
        r = randint(1, pk)
        N2 = pk * pk
        c = ((1 + m * pk) * pow(r, pk, N2)) % N2
        return int(c)

    def decrypt(self,c, pk, sk):
        N2 = pk * pk
        r = pow(c, sk, pk)
        s = mod_inverse(r, pk)
        m = ((c * pow(s, pk, N2)) % N2 - 1) // pk
        return int(m)

    def oplus(self,X, Y, pk):
        Z = (X * Y) % (pk * pk)
        return Z

    def produitParConstante(self,X, y, pk):
        Z = pow(X, y, (pk * pk))
        return Z

    def oppose(self,X, pk):
        Z = pow(X, pk - 1, (pk * pk))
        return Z

class Alice:
    def __init__(self, xa=0,ya=0,taille = 100):
        self.pk, self.sk = paillier.genkeys(taille)
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
        return pow(resu,0.5)

    def checkvecteur(self, vecteur):
        for i in range(len(vecteur)):
            #print(self.decrypt(vecteur[i]))
            if self.decrypt(vecteur[i]) == 0:
               return True
        return False

    def checkVecteurPaire(self,vecteur):
        res = []
        for i in range(len(vecteur)):
            decryptX = self.decrypt(vecteur[i][0])
            decryptY = self.decrypt(vecteur[i][1])
            decryptZ = self.decrypt(vecteur[i][2])
            if decryptZ == 0:
                res = [decryptX,decryptY]
        return res

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

    def distance100(self, xa, ya, xa2, ya2):
        vecteur = []
        xb2 = self.encrypt(pow(self.xb, 2))
        yb2 = self.encrypt(pow(self.yb, 2))
        xAxB = paillier.produitParConstante(xa,self.xb, self.pk)
        yAyB = paillier.produitParConstante(ya, self.yb, self.pk)
        compos2 = paillier.oppose(paillier.produitParConstante(paillier.oplus(xAxB, yAyB, self.pk), 2, self.pk), self.pk)
        distance = paillier.oplus(xa2, paillier.oplus(ya2, paillier.oplus(xb2, paillier.oplus(yb2, compos2, self.pk), self.pk), self.pk), self.pk)
        for i in range(10000):
            randa = randint(1, 100)
            calcul = paillier.produitParConstante(randa,paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk),self.pk)
            #calcul = paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk)
            vecteur.append(calcul)
        shuffle(vecteur)
        return vecteur

    def distance100upgrade(self, xa, ya, xa2, ya2):
        vecteur = []
        #vecteurZ = []
        xb2 = self.encrypt(pow(self.xb, 2))
        yb2 = self.encrypt(pow(self.yb, 2))
        xAxB = paillier.produitParConstante(xa,self.xb, self.pk)
        yAyB = paillier.produitParConstante(ya, self.yb, self.pk)
        compos2 = paillier.oppose(paillier.produitParConstante(paillier.oplus(xAxB, yAyB, self.pk), 2, self.pk), self.pk)
        distance = paillier.oplus(xa2, paillier.oplus(ya2, paillier.oplus(xb2, paillier.oplus(yb2, compos2, self.pk), self.pk), self.pk), self.pk)
        for i in range(10000):
            randa = randint(1, 100)
            calcul = paillier.produitParConstante(randa,paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk),self.pk)
            calcul2 = paillier.produitParConstante(randa,paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk),self.pk)
            calcul3 = paillier.produitParConstante(randa,paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk),self.pk)

            #calcul = paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk)
            vecteur.append([ paillier.oplus(calcul,paillier.encrypt(self.xb,self.pk),self.pk),paillier.oplus(calcul2,paillier.encrypt(self.yb,self.pk),self.pk),calcul3])
        #valeurZ = paillier.oplus(self.encrypt(self.xb),self.encrypt(self.yb),self.pk)
        shuffle(vecteur)
        #shuffle(vecteurZ)
        return vecteur

    def encrypt(self, x):
        return paillier.encrypt(x, self.pk)

if __name__ == '__main__':
    paillier = Paillier()
    # Coordonnées X et Y et taille de la clé d'ALice
    alice = Alice(1,1,100)
    # Clef privée d'Alice et coordonnées X et Y de Bob
    bob = Bob(alice.pk,5,1)

    # Protcole DistanceBob
    print("Distance entre Alice et Bob : ")
    print(alice.distance(bob.distance(alice.encrypt(alice.xa), alice.encrypt(alice.ya))))

    print(" ")

    # Protocole Distance100Bob
    print("Le protocole Distance100Bob montre que Bob est : ")
    if(alice.checkvecteur(bob.distance100(alice.encrypt(alice.xa),alice.encrypt(alice.ya),
                                       alice.encrypt(pow(alice.xa,2)),alice.encrypt(pow(alice.ya,2))))):
        print("trop proche d'Alice (<100)")
    else:
        print("trop loin d'Alice (>100)")

    print(" ")

    # Protocle Distance100BobUpgrade
    print("Le protocole Distance100BobUpgrade montre que Bob est : ")
    vecteur = bob.distance100upgrade(alice.encrypt(alice.xa),alice.encrypt(alice.ya)
                                                                 ,alice.encrypt(pow(alice.xa,2)),alice.encrypt(pow(alice.ya,2)))
    respossible = alice.checkVecteurPaire(vecteur)
    if respossible != []:
        print("aux coordonnées suivantes : ",respossible)
    else:
        print("trop loin d'Alice (>100)")
