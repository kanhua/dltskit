import numpy as np
from scipy.constants import e, k

import matplotlib.pyplot as plt

e_over_k = e / k;


class trapLevel():
    """A class that generates the capacitance transient.
    It is essentially a exponential transient generator.
    However, it uses physical parameters as input and variable names"""


    def __init__(self, activeE, cap_rate=1,delta_c=1):
        """activeE: a list of activation energies of traps in the units of eV.
           cap : capture rate (implementation of capture is still under development)"""

        self.activeE = np.array(activeE)

        def set_param(param):
            if isinstance(param,int) or isinstance(param,float):
                return np.ones((len(activeE),))*param
            elif isinstance(param,list):
                return np.array(param)


        self.capRate = set_param(cap_rate)
        self.delta_c=set_param(delta_c)



    def emRateT(self, T=300):
        """output an array of emission rates based on the
        set temperture T and activation energy array
        T: temperature in Kelvins"""

        self.emRate = np.exp(-self.activeE * e_over_k / T)
        return self.emRate


    def getTransient(self, T, t="default", plotGraph=False, gridnum=1000):
        """Generate capacitance transients
        T: temperature in Kelvins
        t: an array of time. Set none if using default
        gridnum: number of points of t"""

        emRate = self.emRateT(T)

        # define the contributions of exponential components
        self.expContri = list()

        if t == "default":
            defaultFraction = np.linspace(0, 10, num=gridnum)

            t = defaultFraction / emRate[0]

        transientSum = np.zeros((len(t),))
        for emIndex, em in enumerate(emRate):
            tmpVal = self.capRate[emIndex] * self.delta_c[emIndex]*np.exp(-em * t)
            self.expContri.append(tmpVal)
            transientSum = transientSum + tmpVal

        self.time = t;
        self.trans = transientSum

        if plotGraph == True:
            self.plotTransient(T)

        return t, transientSum

    def plotTransient(self, T):
        plt.plot(self.time, self.trans)
        plt.xlabel('time (s)')
        plt.ylabel('simulated signal')
        plt.title('%s K' % str(int(T)))
        plt.savefig('./output/transient %s K.png' % str(int(T)))
        plt.close()


    def getConvDLTS(self, rateWindow=0.5,
                    shift=None,
                    TRange=np.linspace(100, 300, num=50),
                    ttime=None, plotTransient=False):
        """
        Calculate the transient(t1)-transient(t2).
        t1 and t2 are determined by the parameters shift and rateWindow

        t1=shift*len(ttime)
        t2=(shift+ratewindow)*len(ttime)

        return: transient(t1)-transient(t2)
        """

        if shift == None:
            shift = 0.2

        if ttime == None:
            ttime = np.linspace(0, 100, 1000) / self.emRateT(500)[0]
        dltsSig = np.zeros((len(TRange), 2))

        ttimelen = len(ttime)
        for idx, temper in enumerate(TRange):
            time, trans = self.getTransient(temper, ttime, plotGraph=plotTransient)
            dltsSig[idx, 0] = temper

            # This line needs to be fixed, this rate window implementation is bad
            dltsSig[idx, 1] = trans[int(ttimelen * shift)] - trans[int(ttimelen * (shift + rateWindow))]

        return dltsSig


def plotTransients():
    tp = trapLevel([0.005])

    ttime = np.linspace(0, 10, 1000) / tp.emRateT(300)

    _, ts1 = tp.getTransient(210, ttime)
    _, ts2 = tp.getTransient(300, ttime)
    _, ts3 = tp.getTransient(350, ttime)

    plt.plot(ttime, ts1, ttime, ts2, ttime, ts3)

    np.savetxt("testout.csv",np.vstack((ttime,ts1)).transpose())
    plt.show()


if __name__ == "__main__":

    #gen_dlts_plot_emrate()
    #testFit()
    #gen_dlts_plot()
    #plotTransients()
    pass







	
