import sympy as sym
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import WolframSolver as WS
from scipy.io import savemat

s = sym.Symbol('s', real = True)
A = sym.Symbol('A', real = True)
W = sym.Symbol('W', real = True)

class myTF:
    def __init__(self, name):
        self.TF = None
        self.name = name

    def getTF(self):
        return self.TF

    def sim(self, w_array = np.logspace(.1, 6, 100000)):
        # Interface Tfn
        TF = sym.cancel(self.TF)
        n, d = sym.fraction(TF)
        num = np.array(sym.Poly(n, s).all_coeffs()).tolist()
        den = np.array(sym.Poly(d, s).all_coeffs()).tolist()
        num = list(map(float, num))
        den = list(map(float, den))
        tf_ = signal.TransferFunction(num, den)
        self.w, self.mag, self.phase = signal.bode(tf_, w = w_array)

    def plot(self, title=None, fig=None, ax1=None, ax2=None):
        self.sim()
        if(fig is None):
            fig, (ax1, ax2) = plt.subplots(2, 1)
        ax1.set(ylabel='Magnitude')
        ax2.set(xlabel = 'Frequency', ylabel='Phase (º)')
        if(title == None):
            fig.suptitle('myTF response for: '+self.name)
        else:
            plt.suptitle(title)

        ax1.semilogx(self.w, self.mag)  # Bode magnitude plot5
        ax2.semilogx(self.w, self.phase)  # Bode phase plot

    def show(self):
        plt.show()

    def getSim(self):
        return self.w, self.mag, self.phase
    

class BioZ_base(myTF):
    def __init__(self, Rin, Rct, Cdl, Rs, Rgap=1000,  ff = 0.01):
        myTF.__init__(self, "BioImpedance Base")
        global s
        self.TF = - (Rct + Rs*(1+Rct*Cdl*s))/(Rin*(1+Rct*Cdl*s))
    
class BioZ_Cells(myTF):
    def __init__(self, Rin, Rct, Cdl, Rs, Rgap=1000,  ff = 0.01):
        myTF.__init__(self, "BioImpedance Cell Culture")
        global s
        Ce = (1.0 - ff) * Cdl
        Re = Rct / (1.0 - ff)
        Cc = ff * Cdl
        Rc = Rct / ff

        nc1 = Rgap * Rc * Cc
        nc0 = Rc + Rgap
        dc1 = Cc * Rc
        dc0 = 1.0

        ne1 = 0.0
        ne0 = Re
        de1 = Ce * Re
        de0 = 1.0

        n1 = ne0 * nc1
        n0 = ne0 * nc0
        d2 = nc1 * de1
        d1 = nc1 * de0 + nc0 * de1 + ne0 * dc1
        d0 = nc0 * de0 + ne0 * dc0
        div = Rin * d2
        self.TF = (-Rs*d2/div*s**2 -(Rs*d1+n1)/div*s -(Rs*d0+n0)/div)/(Rin*d2/div*s**2 + Rin*d1/div*s + Rin*d0/div)
        

class Comparator(myTF):
    def __init__(self, Vref):
        global A
        myTF.__init__(self, "Comparator")
        self.TF = (4 * Vref) / (np.pi * A)
        
class BPF(myTF):
    def __init__(self, Q, w0, G0):
        myTF.__init__(self, "Band Pass Filter")
        global s
        self.TF = (G0 * (w0 / Q) * s) / (s ** 2 + (w0 / Q) * s + w0 ** 2)
        
class OBT(myTF):
    def __init__(self, BPF, BioZ, Comp):
        myTF.__init__(self, "OBT")
        global s, W, A
        s = sym.Symbol('s', real = True)
        A = sym.Symbol('A', real = True)
        W = sym.Symbol('W', real = True)
        self.TF = 1 + Comp.getTF() * BPF.getTF() * BioZ.getTF()
        self.WF = WS.Wolfram()

    def getPoleCoefficients(self):
        global s
        Eq = sym.cancel(self.TF)
        poles, zeros = sym.fraction(Eq)
        pcoeff = np.array(sym.Poly(poles, s).all_coeffs())
        self.ncoeff = pcoeff / pcoeff[0]

    def getEqusym(self):
        global s, W, A
        self.getPoleCoefficients()
        # Complex Conjugate Pole pair
        eq_order = len(self.ncoeff) - 3
        first = s ** 2 + W ** 2
        second = 0
        i = 0
        eq_symbols = [W, A]
        while (eq_order >= 0):
            sym_str = "x" + str(i)
            t = sym.Symbol(sym_str, real=True)
            if (i == 0):
                second += s ** eq_order
            else:
                second += t * s ** eq_order
                eq_symbols.append(t)
            i += 1
            eq_order -= 1

        temp = (first * second).expand()
        ccpp = np.array(sym.Poly(temp, s).all_coeffs())
        eq_system = []
        i = 1
        while i < len(ccpp):
            temp = ccpp[i] - self.ncoeff[i]
            eq_system.append(temp.cancel())
            i += 1

        return eq_system, eq_symbols

    def solve(self):
        global s, W, A
        self.getPoleCoefficients()
        # Complex Conjugate Pole pair
        eq_order = len(self.ncoeff) - 3
        first = s ** 2 + W ** 2
        second = 0
        i = 0
        eq_symbols = [W, A]
        while (eq_order >= 0):
            sym_str = "x" + str(i)
            t = sym.Symbol(sym_str, real=True)
            if (i == 0):
                second += s ** eq_order
            else:
                second += t * s ** eq_order
                eq_symbols.append(t)
            i += 1
            eq_order -= 1

        temp = (first * second).expand()
        ccpp = np.array(sym.Poly(temp, s).all_coeffs())
        eq_system = []
        i = 1
        while i < len(ccpp):
            temp = ccpp[i] - self.ncoeff[i]
            eq_system.append( temp.cancel() )
            i += 1

        self.WF.createScript(eq_system, eq_symbols)
        self.WF.execute()
        return self.WF.getResults()

