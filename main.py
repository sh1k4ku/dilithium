from Dilithium import Dilithium

def start():
    d = Dilithium()
    M = b'Hello world!'
    sig = d.sign(M)
    print(d.verify(M, sig))

if __name__ == "__main__":
    start()