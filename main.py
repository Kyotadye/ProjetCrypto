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

    def oppose(self, X, pk):
        #Z = pow(X, pk - 1, (pk * pk))
        Z = mod_inverse(X, pk*pk)
        return Z

class PaillierMult:
    def __init__(self,sk=0):
        self.sk = sk

    def mult(self,X, Y, pk):
        return paillier.produitParConstante(X, paillier.decrypt(Y, pk, self.sk), pk)

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
            decryptZ = self.decrypt(vecteur[i][2])
            if decryptZ == 0:
                decryptX = self.decrypt(vecteur[i][0])
                decryptY = self.decrypt(vecteur[i][1])
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
            randa = randint(1, self.pk)
            calcul = paillier.produitParConstante(paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk),randa,self.pk)
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
            randa = randint(1, self.pk)
            calcul = paillier.produitParConstante(paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk),randa,self.pk)
            randa = randint(1, self.pk)
            calcul2 = paillier.produitParConstante(paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk),randa,self.pk)
            randa = randint(1, self.pk)
            calcul3 = paillier.produitParConstante(paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk),randa,self.pk)

            #calcul = paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk)
            vecteur.append([ paillier.oplus(calcul,paillier.encrypt(self.xb,self.pk),self.pk),paillier.oplus(calcul2,paillier.encrypt(self.yb,self.pk),self.pk),calcul3])
        #valeurZ = paillier.oplus(self.encrypt(self.xb),self.encrypt(self.yb),self.pk)
        shuffle(vecteur)
        #shuffle(vecteurZ)
        return vecteur

    def encrypt(self, x):
        return paillier.encrypt(x, self.pk)

class AlicePart2:
    def __init__(self,pk,paillierMult, xa=0,ya=0,xbcrypt=0,ybcrypt=0):
        self.pk = pk
        self.xa = xa
        self.ya = ya
        self.xbcrypt = xbcrypt
        self.ybcrypt = ybcrypt
        self.paillierMult = paillierMult
        self.tabRandom = []

    def encrypt(self, x):
        return paillier.encrypt(x, self.pk)

    def distance100(self):
        vecteur = []
        xa = self.xa
        ya =  self.ya
        xa2 = self.encrypt(self.xa**self.xa)
        ya2 = self.encrypt(self.ya**self.ya)
        xb2 = self.paillierMult.mult(self.xbcrypt, self.xbcrypt, self.pk)
        yb2 = self.paillierMult.mult(self.ybcrypt, self.ybcrypt, self.pk)
        xAxB = paillier.produitParConstante(self.xbcrypt,xa , self.pk)
        yAyB = paillier.produitParConstante(self.ybcrypt, ya, self.pk)
        compos2 = paillier.oppose(paillier.produitParConstante(paillier.oplus(xAxB, yAyB, self.pk), 2, self.pk), self.pk)
        distance = paillier.oplus(xa2, paillier.oplus(ya2, paillier.oplus(xb2, paillier.oplus(yb2, compos2, self.pk), self.pk), self.pk), self.pk)
        for i in range(10000):
            randa = (randint(1, self.pk))
            calcul = paillier.produitParConstante(paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk),randa,self.pk)
            #calcul = paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk)
            vecteur.append(calcul)
        shuffle(vecteur)
        return vecteur

    def distance100ver2(self):
        vecteur = []
        xa = self.xa
        ya =  self.ya
        xa2 = self.encrypt(self.xa**self.xa)
        ya2 = self.encrypt(self.ya**self.ya)
        xb2 = self.paillierMult.mult(self.xbcrypt, self.xbcrypt, self.pk)
        yb2 = self.paillierMult.mult(self.ybcrypt, self.ybcrypt, self.pk)
        xAxB = paillier.produitParConstante(self.xbcrypt,xa , self.pk)
        yAyB = paillier.produitParConstante(self.ybcrypt, ya, self.pk)
        compos2 = paillier.oppose(paillier.produitParConstante(paillier.oplus(xAxB, yAyB, self.pk), 2, self.pk), self.pk)
        distance = paillier.oplus(xa2, paillier.oplus(ya2, paillier.oplus(xb2, paillier.oplus(yb2, compos2, self.pk), self.pk), self.pk), self.pk)
        for i in range(10000):
            randa = (randint(1, self.pk))
            randa2 = (randint(1, self.pk))
            calcul = paillier.oplus(paillier.produitParConstante(paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk),randa,self.pk),self.encrypt(randa2),self.pk)
            self.tabRandom.append(randa2)
            #calcul = paillier.oplus(distance, paillier.oppose(paillier.encrypt(i, self.pk), self.pk), self.pk)
            vecteur.append(calcul)
        shuffle(vecteur)
        return vecteur

    def receiveVecteur(self,xbc,ybc):
        if xbc == 0 and ybc == 0:
            return True
        else:
            return False

    def receiveVecteur2(self,vecteur):
        for i in range(len(vecteur)):
            for j in range(len(self.tabRandom)):
                if vecteur[i] - self.tabRandom[j]  == 0:
                    return self.encrypt(0)
        return self.encrypt(1)

    def receiveVecteur3(self,vecteur):
        for i in range(len(vecteur)):
            for j in range(len(self.tabRandom)):
                if vecteur[i] - self.tabRandom[j]  == 0:
                    return True
        return False



class BobPart2:
    def __init__(self, xb=0, yb=0,taille = 100):
        self.pk, self.sk = paillier.genkeys(taille)
        print("pk = ", self.pk)
        print("sk = ", self.sk)
        print()
        self.xb = xb
        self.yb = yb

    def encrypt(self, x):
        return paillier.encrypt(x, self.pk)

    def decrypt(self, c):
        return paillier.decrypt(c, self.pk, self.sk)

    def checkvecteur(self, vecteur):
        for i in range(len(vecteur)):
            if self.decrypt(vecteur[i]) == 0:
                return 0,0
        return randint(1,self.pk),randint(1,self.pk)

    def checkvecteur2(self,vecteur):
        vecteurdecrypt = [self.decrypt(vecteur[i]) for i in range(len(vecteur))]
        return vecteurdecrypt

    def receiveVecteur(self, nombre):
        if self.decrypt(nombre) == 0:
            return True
        else:
            return False


if __name__ == '__main__':
    choix = int(input("Voulez vous faire tourner la partie 1 ou la partie 2 (1 ou 2) ?"))
    if choix==1:
        paillier = Paillier()
        # Coordonnées X et Y et taille de la clé d'ALice
        Xalice = int(input("Entrez la coordonnée X d'Alice : "))
        Yalice = int(input("Entrez la coordonnée Y d'Alice : "))
        taille = int(input("Entrez la taille de la clé : "))
        alice = Alice(Xalice,Yalice,taille)
        # Clef privée d'Alice et coordonnées X et Y de Bob
        Xbob = int(input("Entrez la coordonnée X de Bob : "))
        Ybob = int(input("Entrez la coordonnée Y de Bob : "))
        bob = Bob(alice.pk,Xbob,Ybob)

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

    else:
        paillier = Paillier()
        #Coordonnée de Bob et taille de la clé
        Xbob = int(input("Entrez la coordonnée X de Bob : "))
        Ybob = int(input("Entrez la coordonnée Y de Bob : "))
        taille = int(input("Entrez la taille de la clé : "))
        bob2 = BobPart2(Xbob,Ybob,taille)
        paillierMult = PaillierMult(bob2.sk)
        #Clef publique de Bob et coordonnées de Alice
        Xalice = int(input("Entrez la coordonnée X d'Alice : "))
        Yalice = int(input("Entrez la coordonnée Y d'Alice : "))
        alice2 = AlicePart2(bob2.pk,paillierMult,Xalice,Yalice,bob2.encrypt(bob2.xb),bob2.encrypt(bob2.yb))

        print("Le protocole Distance100 montre que : ")
        vecteur = alice2.distance100()
        respossible = bob2.checkvecteur(vecteur)
        if alice2.receiveVecteur(respossible[0],respossible[1]):
            print("Bob est trop proche d'Alice (<100)")
        else:
            print("Bob est trop loin d'Alice (>100)")

        print(" ")

        vecteur = alice2.distance100ver2()
        respossible = bob2.checkvecteur2(vecteur)
        respossiblesuite = alice2.receiveVecteur2(respossible)
        resfinal = bob2.receiveVecteur(respossiblesuite)
        print("Le protocole Distance100ver2 montre que : ")
        if resfinal:
            print("Bob est trop proche d'Alice (<100)")
        else:
            print("Bob est trop loin d'Alice (>100)")

        print(" ")

        vecteur = alice2.distance100ver2()
        respossible = bob2.checkvecteur2(vecteur)
        print("Le protocole Distance100ver3 montre que : ")
        if alice2.receiveVecteur3(respossible):
            print("Bob est trop proche d'Alice (<100) et seul Alice est au courant")
        else:
            print("Bob est trop loin d'Alice (>100) et seul Alice est au courant")